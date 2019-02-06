from dbApp.models import (SymportalFramework, DataSet, ReferenceSequence,
                          DataSetSampleSequence, DataSetSample, CladeCollection)
import sys
import os
import shutil
import subprocess
import glob
import pandas as pd
import general
import json
from multiprocessing import Queue, Manager, Process
from django import db
from analysis_classes import InitialMothurWorker, TaxonomicScreeningWorker, TaxonomicScreeningHandler, BlastnAnalysis
from datetime import datetime
from general import write_list_to_destination, read_defined_file_to_list, create_dict_from_fasta, make_new_blast_db
import pickle

class DataLoading:
    # The clades refer to the phylogenetic divisions of the Symbiodiniaceae. Most of them are represented at the genera
    # level. E.g. All members of clade C belong to the genus Cladocopium.
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    def __init__(self, data_set_uid, user_input_path, datasheet_path, screen_sub_evalue=False, debug=False, num_proc=1):
        self.dataset_object = DataSet.objects.get(id=data_set_uid)
        self.user_input_path = user_input_path
        self.output_directory = self.setup_output_directory()
        self.set_temp_working_directory()
        self.create_temp_wkd()
        self.num_proc = num_proc
        # Used for the pre minimum entropy decomposition (MED), quality controlled (QC), sequences dump
        # Within this directory we will have a directory for each sample that will contain clade
        # separated name and fasta pairs
        self.pre_med_sequence_output_directory_path = self.create_pre_med_write_out_directory_path()
        # data can be loaded either as paired fastq or fastq.gz files or as a single compressed file containing
        # the afore mentioned paired files.
        self.is_single_file_or_paired_input = self.determine_if_single_file_or_paired_input()
        self.datasheet_path = datasheet_path
        self.debug = debug
        self.symclade_db_directory_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB'))
        self.symclade_db_full_path = os.path.join(self.symclade_db_directory_path, 'symClade.fa')
        self.sample_meta_info_df = None
        self.list_of_samples_names = None
        self.list_of_fastq_file_names_in_wkd = None
        self.list_of_fastq_files_in_wkd = None
        self.path_to_mothur_batch_file_for_dot_file_creation = None
        self.fastqs_are_gz_compressed = None
        self.path_to_latest_mothur_batch_file = None
        self.fastq_file_to_sample_name_dict = None
        self.num_of_samples = None
        self.sample_fastq_pairs = None
        self.samples_that_caused_errors_in_qc_list = []
        self.post_initial_qc_name_file_name = None
        self.post_initial_qc_fasta_file_name = None
        # args for the taxonomic screening
        self.screen_sub_evalue = screen_sub_evalue
        self.new_seqs_added_in_iteration = 0
        self.new_seqs_added_running_total = 0
        self.discarded_seqs_fasta = None
        self.checked_samples_with_no_additional_symbiodinium_sequences = []
        self.taxonomic_screening_handler = None
        self.sequences_to_screen_fasta_as_list = []
        self.sequences_to_screen_fasta_path = os.path.join(
            self.temp_working_directory, 'taxa_screening_seqs_to_screen.fa'
        )
        # the number of samples that a sub evalue sequences must be
        # found in for us to carry it through for taxonomic screening
        self.required_sample_support_for_sub_evalue_sequencs = 3
        # the number of sequences in the 10 matches that must be annotated as Symbiodinium or Symbiodiniaceae in order
        # for a sequences to be added into the reference symClade database.
        self.required_symbiodinium_matches = 3

    def execute(self):
        self.copy_and_decompress_input_files_to_temp_wkd()

        # the stability file generated here is used as the base of the initial mothur QC
        if self.datasheet_path:
            self.generate_stability_file_and_data_set_sample_objects_with_datasheet()
        else:
            self.generate_stability_file_and_data_set_sample_objects_without_datasheet()


        self.if_symclade_binaries_not_present_remake_db()

        # First execute the initial mothur QC. This should leave us with a set of fastas, name and groupfiles etc in
        # a directory from each sample.
        self.execute_worker_initial_mothur()
        self.taxonomic_screening()

    def taxonomic_screening(self):
        """
        This screening is to identify the Symbiodinium sequences in the sequences that have already been through the
        initial quality control for each of the samples. We have two versions of this taxonomic screening.
        The version that screens the 'sub_e_value' sequences and the version that doesn't. The sub_e_values
        sequences are those that gave a match to
        the symClade.fa database (reference database containing Symbiodinium samples) but that had e-values lower than the
        required threshold. For local instances the sub e value screening is turned off because the screening requires
        running blastn analyses againt the NCBI nt databse. This is likely to big a computational effort for most local
        users who will be running on a personal machine. On the remote system (i.e. when run on the servers that host
        the remote SymPortal the sub e value screening will be run by default. This sub e value screening is done
        iteratively. If sub evalues are found in the initial taxonomic screening (blastn analysis against the symClade.fa
        databse) then these sub evalues are screened againt the nt datbase. If some of the sequnces are found to be
        symbiodinium, these are added to the symClade datbase, and the initial screening is run again. This processs is
        repeated until no further sub evalues are found or if none of the sub evlaues being screened are found to be
        Symbiodinium in origin.
        To make this faster, rather than check every sample in every iteration we will keep track of when a sample
        has been checked and found to contain no non -ymbiodinium sequences.
        This way we can skip these sequences in the iterative rounds of screening.
        """


        fasta_out = None
        fasta_out_path = None
        if self.screen_sub_evalue:

            self.create_symclade_backup_incase_of_accidental_deletion_of_corruption()

            while 1:
                # This method simply identifies whether there are sequences that need screening.
                # Everytime that the execute_worker is run it pickles out the files needed for the next worker
                self.make_fasta_of_sequences_that_need_taxa_screening()

                if self.sequences_to_screen_fasta_as_list:
                    # Now do the screening
                    # The outcome of this will be an updated symClade.fa that we should then make a blastdb from it.
                    self.screen_sub_e_seqs()
                    if self.new_seqs_added_in_iteration == 0:
                        break
                else:
                    break
                # TODO we have got here with the refactoring and class

        else:
            # if not doing the screening we can simply run the execute_worker_taxa_screening once.
            # During its run it will have output all of the files we need to run the following workers.
            # We can also run the generate_and_write_below_evalue_fasta_for_screening function to write out a
            # fasta of significant sequences we can then report to the user using that object
            self.make_fasta_of_sequences_that_need_taxa_screening()


        # Create the queues that will hold the sample information
        input_q = Queue()

        # The list that contains the names of the samples that returned errors during the initial mothur
        worker_manager = Manager()
        error_sample_list_shared = worker_manager.list(error_sample_list)

        # we want to collect all of the discarded sequences so that we can print this to a fasta and
        # output this in the data_set's submission output directory
        # to do this we'll need a list that we can use to collect all of these sequences
        # once we have collected them all we can then set them and use this to make a fasta to output
        # we want to do this on a sample by sample basis, so we now write out directly from
        # the worker as well as doing the 'total' method.
        # we have already created the output dir early on and I will pass it down to here so that we can pass
        # it to the below method
        list_of_discarded_sequences = worker_manager.list()

        # create the directory that will be used for the output for the output of the throw away sequences on
        # a sample by sample basis (one fasta and one name file per sample)
        throw_away_seqs_dir = '{}/throw_awayseqs'.format(output_dir)
        os.makedirs(throw_away_seqs_dir, exist_ok=True)

        # load up the input q
        for contigPair in sample_fastq_pairs:
            input_q.put(contigPair)

        # load in the STOPs
        for n in range(num_proc):
            input_q.put('STOP')

        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        sys.stdout.write('\nPerforming QC\n')

        for n in range(num_proc):
            p = Process(target=worker_taxonomy_write_out, args=(
                input_q, error_sample_list_shared, wkd, data_set_uid, list_of_discarded_sequences,
                throw_away_seqs_dir))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

        # now create a fasta from the list
        discarded_seqs_fasta = []
        discarded_seqs_name_counter = 0
        for seq in set(list(list_of_discarded_sequences)):
            discarded_seqs_fasta.extend(
                ['>discard_seq_{}_data_sub_{}'.format(discarded_seqs_name_counter, data_set_uid), seq])
            discarded_seqs_name_counter += 1

        if self.screen_sub_evalue:
            # If we are screening then we want to be returning the number of seqs added
            return discarded_seqs_fasta
        else:
            # if we are not screening then we want to return the fasta that contains the significant sub_e_value seqs
            return fasta_out, fasta_out_path, discarded_seqs_fasta

    def screen_sub_e_seqs(self):
        """This function screens a fasta file to see if the sequences are Symbiodinium in origin.
        This fasta file contains the below_e_cutoff sequences that need to be screened. These sequences are the sequences
        that were found to have a match in the initial screening against the symClade.fa database but were below
        the required evalue cut off. Here we run these sequences against the entire NCBI 'nt' database to verify if they
        or of Symbiodinium origin of not.
        The fasta we are screening only contains seuqences that were found in at least 3 samples.
        We will call a sequence symbiodinium if it has a match that covers at least 95% of its sequence at a 60% or higher
        identity. It must also have Symbiodinium or Symbiodiniaceae in the name. We will also require that a
        sub_e_value seq has at least the required_sybiodinium_matches (3 at the moment) before we call it Symbiodinium.
        """

        blastn_analysis_object = BlastnAnalysis(
            input_file_path=self.sequences_to_screen_fasta_path,
            output_file_path=os.path.join(self.temp_working_directory, 'blast.out'),
            max_target_seqs=10,
            num_threads=str(self.num_proc)
        )

        blastn_analysis_object.execute()

        blast_output_dict = blastn_analysis_object.return_blast_results_dict()

        query_sequences_verified_as_symbiodinium_list = self.get_list_of_seqs_in_blast_result_that_are_symbiodinium(
            blast_output_dict
        )

        self.new_seqs_added_in_iteration = len(query_sequences_verified_as_symbiodinium_list)
        self.new_seqs_added_running_total += self.new_seqs_added_in_iteration

        if query_sequences_verified_as_symbiodinium_list:
            self.taxa_screening_update_symclade_db_with_new_symbiodinium_seqs(query_sequences_verified_as_symbiodinium_list)


    def taxa_screening_update_symclade_db_with_new_symbiodinium_seqs(self, query_sequences_verified_as_symbiodinium_list):
        new_symclade_fasta_as_list = self.taxa_screening_make_new_fasta_of_screened_seqs_to_be_added_to_symclade_db(
            query_sequences_verified_as_symbiodinium_list
        )
        combined_fasta = self.taxa_screening_combine_new_symclade_seqs_with_current(new_symclade_fasta_as_list)
        self.taxa_screening_make_new_symclade_db(combined_fasta)

    def taxa_screening_make_new_symclade_db(self, combined_fasta):
        write_list_to_destination(self.symclade_db_full_path, combined_fasta)
        make_new_blast_db(input_fasta_to_make_db_from=self.symclade_db_full_path, db_title='symClade')

    def taxa_screening_combine_new_symclade_seqs_with_current(self, new_symclade_fasta_as_list):
        old_symclade_fasta_as_list = read_defined_file_to_list(self.symclade_db_full_path)
        combined_fasta = new_symclade_fasta_as_list + old_symclade_fasta_as_list
        return combined_fasta

    def taxa_screening_make_new_fasta_of_screened_seqs_to_be_added_to_symclade_db(self,
                                                                                  query_sequences_verified_as_symbiodinium_list):
        screened_seqs_fasta_dict = create_dict_from_fasta(
            fasta_path=self.sequences_to_screen_fasta_path
        )
        new_symclade_fasta_as_list = []
        for name_of_symbiodinium_sequence_to_add_to_symclade_db in query_sequences_verified_as_symbiodinium_list:
            new_symclade_fasta_as_list.extend(
                [
                    f'>{name_of_symbiodinium_sequence_to_add_to_symclade_db}',
                    f'{screened_seqs_fasta_dict[name_of_symbiodinium_sequence_to_add_to_symclade_db]}'
                ]
            )
        return new_symclade_fasta_as_list

    def get_list_of_seqs_in_blast_result_that_are_symbiodinium(self, blast_output_dict):
        query_sequences_verified_as_symbiodinium_list = []
        for query_sequence_name, blast_result_list_for_query_sequence in blast_output_dict.items():
            sym_count = 0
            for result_str in blast_result_list_for_query_sequence:
                if 'Symbiodinium' in result_str or 'Symbiodiniaceae' in result_str:
                    percentage_coverage = float(result_str.split('\t')[4])
                    percentage_identity_match = float(result_str.split('\t')[3])
                    if percentage_coverage > 95 and percentage_identity_match > 60:
                        sym_count += 1
                        if sym_count == self.required_symbiodinium_matches:
                            query_sequences_verified_as_symbiodinium_list.append(query_sequence_name)
                            break
        return query_sequences_verified_as_symbiodinium_list

    def make_fasta_of_seqs_found_in_more_than_two_samples_that_need_screening(self):
        """ The below_e_cutoff_dict has nucleotide sequencs as the
        key and the number of samples that sequences was found in as the value.
        """
        sub_evalue_nuclotide_sequence_to_number_of_samples_found_in_dict = dict(
            self.taxonomic_screening_handler.sub_evalue_sequence_to_num_sampes_found_in_mp_dict)
        self.sequences_to_screen_fasta_as_list = []
        sequence_number_counter = 0
        for nucleotide_sequence, num_samples_found_in in sub_evalue_nuclotide_sequence_to_number_of_samples_found_in_dict.items():
            if num_samples_found_in >= self.required_sample_support_for_sub_evalue_sequencs:
                # then this is a sequences that was found in three or more samples
                clade_of_sequence = self.taxonomic_screening_handler.sub_evalue_nucleotide_sequence_to_clade_mp_dict[
                    nucleotide_sequence
                ]
                self.sequences_to_screen_fasta_as_list.extend(
                    [
                        f'>sub_e_seq_count_{sequence_number_counter}_'
                        f'{self.dataset_object.id}_{num_samples_found_in}_'
                        f'{clade_of_sequence}',
                        nucleotide_sequence
                    ]
                )
                sequence_number_counter += 1
        if self.sequences_to_screen_fasta_as_list:
            write_list_to_destination(self.sequences_to_screen_fasta_path, self.sequences_to_screen_fasta_as_list)




    def create_symclade_backup_incase_of_accidental_deletion_of_corruption(self):
        back_up_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB', 'symClade_backup'))
        os.makedirs(back_up_dir, exist_ok=True)
        src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB')) + '/symClade.fa'
        time_stamp = str(datetime.now()).replace(' ', '_').replace(':', '-')
        dst_fasta_path = back_up_dir + '/symClade_{}.fa'.format(time_stamp)
        dst_readme_path = back_up_dir + '/symClade_{}.readme'.format(time_stamp)
        # then write a copy to it.
        shutil.copyfile(src_path, dst_fasta_path)
        # Then write out a very breif readme
        read_me = [
            'This is a symClade.fa backup created during datasubmission of data_set ID: {}'.format(
                self.dataset_object.id)]
        write_list_to_destination(dst_readme_path, read_me)

    def make_fasta_of_sequences_that_need_taxa_screening(self):
        self.init_taxonomic_screening_handler()

        # the self.taxonomic_screening_handler.sub_evalue_sequence_to_num_sampes_found_in_mp_dict is populated here
        self.execute_taxonomy_screening()

        self.taxa_screening_update_checked_samples_list()

        self.make_fasta_of_seqs_found_in_more_than_two_samples_that_need_screening()

    def taxa_screening_update_checked_samples_list(self):
        self.checked_samples_with_no_additional_symbiodinium_sequences = \
            list(self.taxonomic_screening_handler.checked_samples_mp_list)

    def execute_taxonomy_screening(self):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        sys.stdout.write('\nPerforming QC\n')
        for n in range(self.num_proc):
            new_process = Process(
                target=self.worker_taxonomy_screening,
                args=(
                    self.taxonomic_screening_handler.input_queue,
                    self.taxonomic_screening_handler.sub_evalue_sequence_to_num_sampes_found_in_mp_dict,
                    self.taxonomic_screening_handler.error_samples_mp_list,
                    self.taxonomic_screening_handler.checked_samples_mp_list,
                    self.taxonomic_screening_handler.sub_evalue_nucleotide_sequence_to_clade_mp_dict,
                ))

            all_processes.append(new_process)
            new_process.start()
        for process in all_processes:
            process.join()

    def init_taxonomic_screening_handler(self):
        self.taxonomic_screening_handler = TaxonomicScreeningHandler(
            samples_that_caused_errors_in_qc_list=self.samples_that_caused_errors_in_qc_list,
            checked_samples_list=self.checked_samples_with_no_additional_symbiodinium_sequences,
            list_of_samples_names=self.list_of_samples_names, num_proc=self.num_proc
        )

    def worker_taxonomy_screening(
            self, input_q, e_val_collection_mp_dict, err_smpl_list_mp_list,
            checked_samples_mp_list, sub_evalue_nucleotide_sequence_to_clade_mp_dict):
        """
        :param input_q: The multiprocessing queue that holds a list of the sample names
        :param e_val_collection_dict: This is a managed dictionary where key is a nucleotide sequence that has:
        1 - provided a match in the blast analysis
        2 - is of suitable size
        3 - but has an evalue match below the cuttof
        the value of the dict is an int that represents how many samples this nucleotide sequence was found in
        :param err_smpl_list: This is a managed list containing sample names of the samples that had errors during the
        initial mothur qc and therefore don't require taxonomic screening performed on them
        :param checked_samples: This is a list of sample names for samples that were found to contain only
        Symbiodinium sequences or have already had all potential Symbiodinium sequences screened and so don't
        require any further taxonomic screening
        :return: Whilst this doesn't return anything a number of objects are picked out in each of the local
        working directories for use in the workers that follow this one.
        """
        for sample_name in iter(input_q.get, 'STOP'):

            # If the sample gave an error during the inital mothur then we don't consider it here.
            if sample_name in err_smpl_list_mp_list:
                continue

            # A sample will be in this list if we have already performed this worker on it and none of its sequences
            # gave matches to the symClade database at below the evalue threshold
            if sample_name in checked_samples_mp_list:
                continue

            taxonomic_screening_worker = TaxonomicScreeningWorker(
                sample_name=sample_name, wkd=self.temp_working_directory,
                path_to_symclade_db=self.symclade_db_full_path, debug=self.debug,
                checked_samples_mp_list=checked_samples_mp_list, e_val_collection_mp_dict=e_val_collection_mp_dict,
                sub_evalue_nucleotide_sequence_to_clade_mp_dict=sub_evalue_nucleotide_sequence_to_clade_mp_dict)


            taxonomic_screening_worker.execute()

    def execute_worker_initial_mothur(self):

        if not self.sample_fastq_pairs:
            self.exit_and_del_data_set_sample('Sample fastq pairs list empty')

        # Create the queues that will hold the sample information
        input_queue_containing_pairs_of_fastq_file_paths = Queue()

        worker_manager = Manager()
        samples_that_caused_errors_in_qc_list = worker_manager.list()

        for fastq_path_pair in self.sample_fastq_pairs:
            input_queue_containing_pairs_of_fastq_file_paths.put(fastq_path_pair)

        for n in range(self.num_proc):
            input_queue_containing_pairs_of_fastq_file_paths.put('STOP')
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        sys.stdout.write('\nPerforming QC\n')

        for n in range(self.num_proc):
            p = Process(target=self.worker_initial_mothur,
                        args=(
                            input_queue_containing_pairs_of_fastq_file_paths,
                            samples_that_caused_errors_in_qc_list, self.temp_working_directory,
                            self.dataset_object.id, self.debug))

            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

        self.samples_that_caused_errors_in_qc_list = list(samples_that_caused_errors_in_qc_list)

    def worker_initial_mothur(self, input_q, error_sample_list):
        """
        This worker performs the pre-MED processing that is primarily mothur-based.
        This QC includes making contigs, screening for ambigous calls (0 allowed) and homopolymer (maxhomop=5)
        Discarding singletons and doublets, in silico PCR. It also checks whether sequences are rev compliment.
        This is all done through the use of an InitialMothurWorker class which in turn makes use of the MothurAnalysis
        class that does the heavy lifting of running the mothur commands in sequence.
        """

        for contigpair in iter(input_q.get, 'STOP'):

            initial_morthur_worker = InitialMothurWorker(
                contig_pair=contigpair,
                data_set_object=self.dataset_object,
                temp_wkd=self.temp_working_directory,
                debug=self.debug
            )

            try:
                initial_morthur_worker.execute()
            except RuntimeError as e:
                error_sample_list.append(e.sample_name)

        return

    def if_symclade_binaries_not_present_remake_db(self):
        list_of_binaries_that_should_exist = [
            self.dataset_object.reference_fasta_database_used + extension for extension in ['.nhr', '.nin', '.nsq']
        ]

        contents_of_symclade_directory = os.listdir(self.symclade_db_directory_path)

        binary_count = 0
        for item in contents_of_symclade_directory:
            if item in list_of_binaries_that_should_exist:
                binary_count += 1

        if binary_count != 3:
            # then some of the binaries are not present and we need to remake the blast dictionary
            if not self.debug:
                general.make_new_blast_db(
                    input_fasta_to_make_db_from=self.symclade_db_full_path,
                    db_title='symClade', pipe_stdout_sterr=True)
            elif self.debug:
                general.make_new_blast_db(
                    input_fasta_to_make_db_from=self.symclade_db_full_path,
                    db_title='symClade', pipe_stdout_sterr=False)

            # now verify that the binaries have been successfully created
            list_of_dir = os.listdir(self.symclade_db_directory_path)
            binary_count = 0
            for item in list_of_dir:
                if item in list_of_binaries_that_should_exist:
                    binary_count += 1
            if binary_count != 3:
                self.exit_and_del_data_set_sample('Failure in creating blast binaries')

    def generate_stability_file_and_data_set_sample_objects_without_datasheet(self):

        self.list_of_fastq_file_names_in_wkd = [a for a in os.listdir(self.temp_working_directory) if 'fastq' in a]

        self.generate_and_write_mothur_batch_file_for_dotfile_creation()
        general.execute_mothur_batch_file_with_piped_stoud_sterr(self.path_to_latest_mothur_batch_file)

        self.generate_and_write_new_stability_file_without_datasheet()

        self.create_data_set_sample_objects_in_bulk_without_datasheet()

    def create_data_set_sample_objects_in_bulk_without_datasheet(self):
        list_of_sample_objects = []
        sys.stdout.write('\nCreating data_set_sample objects\n')
        for sampleName in self.list_of_samples_names:
            print('\rCreating data_set_sample {}'.format(sampleName))
            # Create the data_set_sample objects in bulk.
            # The cladal_seq_totals property of the data_set_sample object keeps track of the seq totals for the
            # sample divided by clade. This is used in the output to keep track of sequences that are not
            # included in cladeCollections
            clade_zeroes_list = [0 for _ in self.clade_list]
            empty_cladal_seq_totals = json.dumps(clade_zeroes_list)

            dss = DataSetSample(name=sampleName, data_submission_from=self.dataset_object,
                                cladal_seq_totals=empty_cladal_seq_totals)
            list_of_sample_objects.append(dss)
        # http://stackoverflow.com/questions/18383471/django-bulk-create-function-example
        DataSetSample.objects.bulk_create(list_of_sample_objects)

    def generate_and_write_new_stability_file_without_datasheet(self):
        self.read_in_mothur_dot_file_creation_output()
        # We will use the stability file in the rest of the mothur qc but also to make the DataSetSamples (below)
        end_index = self.identify_sample_names_without_datasheet()
        new_stability_file = self.generate_new_stability_file_without_datasheet(end_index)
        # write out the new stability file
        general.write_list_to_destination(f'{self.temp_working_directory}/stability.files', new_stability_file)


    def generate_new_stability_file_without_datasheet(self, end_index):
        new_stability_file = []
        for stability_file_line in self.sample_fastq_pairs:
            pair_components = stability_file_line.split('\t')
            # I am going to use '[dS]' as a place holder for a dash in the sample names
            # Each line of the stability file is a three column format with the first
            # column being the sample name. The second and third are the full paths of the .fastq files
            # the sample name at the moment is garbage, we will extract the sample name from the
            # first fastq path using the end_index that we determined above

            new_stability_file.append(
                '{}\t{}\t{}'.format(
                    pair_components[1].split('/')[-1][:-end_index].replace('-', '[dS]'),
                    pair_components[1],
                    pair_components[2]))
        return new_stability_file

    def get_num_chars_in_common_with_fastq_names(self):
        i = 1
        while 1:
            list_of_endings = []
            for file in self.list_of_fastq_files_in_wkd:
                list_of_endings.append(file[-i:])
            if len(set(list_of_endings)) > 2:
                break
            else:
                i += 1
                # then this is one i too many and our magic i was i-1
        end_index = i - 1
        return end_index

    def get_sample_names_from_fastq_files_using_index(self, end_index):
        list_of_names_non_unique = []
        for file in self.list_of_fastq_files_in_wkd:
            list_of_names_non_unique.append(file[:-end_index])
        list_of_sample_names = list(set(list_of_names_non_unique))
        if len(list_of_sample_names) != len(self.list_of_fastq_files_in_wkd) / 2:
            warning_str = 'Error in sample name extraction'
            self.exit_and_del_data_set_sample(warning_str)
        self.list_of_samples_names = list_of_sample_names

    def identify_sample_names_without_datasheet(self):
        # I think the simplest way to get sample names is to find what parts are common between all samples
        # well actually 50% of the samples so that we also remove the R1 and R2 parts.
        end_index = self.get_num_chars_in_common_with_fastq_names()
        self.get_sample_names_from_fastq_files_using_index(end_index)

        return end_index
    def generate_stability_file_and_data_set_sample_objects_with_datasheet(self):
        # Create a pandas df from the data_sheet if it was provided
        # allow the data_sheet to be in a .csv format or .xlsx format. This is so that we can store a datasheet
        # in the github repo in a non-binary format
        # The sample_meta_df that is created from the data_sheet should be identical irrespective of whether a .csv
        # or a .xlsx is submitted.
        sample_meta_info_df = self.create_sample_meta_info_dataframe_from_datasheet_path()

        # if we are given a data_sheet then use the sample names given as the DataSetSample object names
        self.list_of_samples_names = sample_meta_info_df.index.values.tolist()

        # we should verify that all of the fastq files listed in the sample_meta_df
        # are indeed found in the directory that we've been given
        self.check_all_fastqs_in_datasheet_exist()

        self.generate_and_write_mothur_batch_file_for_dotfile_creation()

        # noinspection PyPep8
        general.execute_mothur_batch_file_with_piped_stoud_sterr(self.path_to_latest_mothur_batch_file)

        # we will also need to know how to relate the sample names to the fastq files
        # for this we will make a dict of fastq file name to sample
        self.create_fastq_file_to_sample_name_dict()
        # We will use the stability file in the rest of the mothur qc but also to make the DataSetSamples (below)
        self.generate_and_write_new_stability_file_with_data_sheet()

        # TODO we got this far in the refactoring. We are working through this create_data_sub code first
        # we still need to transfer the equivalent into the DataLoading class.

        # Create data_set_sample instances
        self.create_data_set_sample_objects_in_bulk_with_datasheet()


    def create_data_set_sample_objects_in_bulk_with_datasheet(self):
        list_of_data_set_sample_objects = []
        sys.stdout.write('\nCreating data_set_sample objects\n')
        for sampleName in self.list_of_samples_names:
            print('\rCreating data_set_sample {}'.format(sampleName))
            # Create the data_set_sample objects in bulk.
            # The cladal_seq_totals property of the data_set_sample object keeps track of the seq totals for the
            # sample divided by clade. This is used in the output to keep track of sequences that are not
            # included in cladeCollections
            empty_cladal_seq_totals = json.dumps([0 for _ in self.clade_list])

            dss = DataSetSample(name=sampleName, data_submission_from=self.dataset_object,
                                cladal_seq_totals=empty_cladal_seq_totals,
                                sample_type=self.sample_meta_info_df.loc[sampleName, 'sample_type'],
                                host_phylum=self.sample_meta_info_df.loc[sampleName, 'host_phylum'],
                                host_class=self.sample_meta_info_df.loc[sampleName, 'host_class'],
                                host_order=self.sample_meta_info_df.loc[sampleName, 'host_order'],
                                host_family=self.sample_meta_info_df.loc[sampleName, 'host_family'],
                                host_genus=self.sample_meta_info_df.loc[sampleName, 'host_genus'],
                                host_species=self.sample_meta_info_df.loc[sampleName, 'host_species'],
                                collection_latitude=self.sample_meta_info_df.loc[sampleName, 'collection_latitude'],
                                collection_longitude=self.sample_meta_info_df.loc[sampleName, 'collection_longitude'],
                                collection_date=self.sample_meta_info_df.loc[sampleName, 'collection_date'],
                                collection_depth=self.sample_meta_info_df.loc[sampleName, 'collection_depth']
                                )
            list_of_data_set_sample_objects.append(dss)
        # http://stackoverflow.com/questions/18383471/django-bulk-create-function-example
        DataSetSample.objects.bulk_create(list_of_data_set_sample_objects)


    def generate_and_write_new_stability_file_with_data_sheet(self):
        # Convert the group names in the stability.files so that the dashes are converted to '[ds]',
        # So for the mothur we have '[ds]'s. But for all else we convert these '[ds]'s to dashes
        self.read_in_mothur_dot_file_creation_output()

        new_stability_file = self.generate_new_stability_file_with_data_sheet()

        general.write_list_to_destination(
            r'{0}/stability.files'.format(self.temp_working_directory),
            new_stability_file
        )


    def generate_new_stability_file_with_data_sheet(self):
        new_stability_file = []
        for stability_file_line in self.sample_fastq_pairs:
            pair_components = stability_file_line.split('\t')
            # I am going to use '[dS]' as a place holder for a dash in the sample names
            # Each line of the stability file is a three column format with the first
            # column being the sample name. The second and third are the full paths of the .fastq files
            # the sample name at the moment is garbage, we will identify the sample name from the
            # first fastq path using the fastq_file_to_sample_name_dict
            new_stability_file.append(
                '{}\t{}\t{}'.format(
                    self.fastq_file_to_sample_name_dict[pair_components[1].split('/')[-1]].replace('-', '[dS]'),
                    pair_components[1],
                    pair_components[2]))
        return new_stability_file

    def read_in_mothur_dot_file_creation_output(self):
        self.sample_fastq_pairs = general.read_defined_file_to_list(f'{self.temp_working_directory}/stability.files')

    def create_fastq_file_to_sample_name_dict(self):
        fastq_file_to_sample_name_dict = {}
        for sample_index in self.sample_meta_info_df.index.values.tolist():
            fastq_file_to_sample_name_dict[self.sample_meta_info_df.loc[sample_index, 'fastq_fwd_file_name']] = sample_index
            fastq_file_to_sample_name_dict[self.sample_meta_info_df.loc[sample_index, 'fastq_rev_file_name']] = sample_index
        self.fastq_file_to_sample_name_dict = fastq_file_to_sample_name_dict

    def execute_mothur_batch_file_with_piped_stoud_sterr(path_to_mothur_batch_file):
        subprocess.run(['mothur', path_to_mothur_batch_file],
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)

    def generate_and_write_mothur_batch_file_for_dotfile_creation(self):
        self.check_if_fastqs_are_gz_compressed()
        mothur_batch_file_as_list = self.generate_mothur_batch_file_for_dotfile_creation_as_list()
        self.path_to_latest_mothur_batch_file = f'{self.temp_working_directory}/mothur_batch_file_makeFile'
        general.write_list_to_destination(self.path_to_latest_mothur_batch_file, mothur_batch_file_as_list)

    def generate_mothur_batch_file_for_dotfile_creation_as_list(self):
        if self.fastqs_are_gz_compressed:
            mothur_batch_file = [
                r'set.dir(input={0})'.format(self.temp_working_directory),
                r'set.dir(output={0})'.format(self.temp_working_directory),
                r'make.file(inputdir={0}, type=gz, numcols=3)'.format(self.temp_working_directory)
            ]
        else:
            mothur_batch_file = [
                r'set.dir(input={0})'.format(self.temp_working_directory),
                r'set.dir(output={0})'.format(self.temp_working_directory),
                r'make.file(inputdir={0}, type=fastq, numcols=3)'.format(self.temp_working_directory)
            ]
        return mothur_batch_file

    def check_if_fastqs_are_gz_compressed(self):
        if self.list_of_fastq_files_in_wkd[0].endswith('fastq.gz'):
            self.fastqs_are_gz_compressed = True
        elif self.list_of_fastq_files_in_wkd[0].endswith('fastq'):
            self.fastqs_are_gz_compressed = False
        else:
            warning_str = f'Unrecognised format of sequecing file: {self.list_of_fastq_files_in_wkd[0]}'
            self.exit_and_del_data_set_sample(warning_str)

    def check_all_fastqs_in_datasheet_exist(self):
        self.list_of_fastq_files_in_wkd = [_ for _ in os.listdir(self.temp_working_directory) if 'fastq' in _]
        list_of_meta_gz_files = self.get_list_of_fastq_file_names_that_should_be_in_directory()
        self.if_fastq_files_missing_sys_exit(list_of_meta_gz_files)

    def if_fastq_files_missing_sys_exit(self, list_of_meta_gz_files):
        for fastq in list_of_meta_gz_files:
            if fastq not in self.list_of_fastq_files_in_wkd:
                warning_str = f'{fastq} not found'
                self.exit_and_del_data_set_sample(warning_str)

    def exit_and_del_data_set_sample(self, warning_str):
        self.dataset_object.delete()
        sys.exit(warning_str)

    def get_list_of_fastq_file_names_that_should_be_in_directory(sample_meta_info_df):
        list_of_meta_gz_files = []
        list_of_meta_gz_files.extend(sample_meta_info_df['fastq_fwd_file_name'].values.tolist())
        list_of_meta_gz_files.extend(sample_meta_info_df['fastq_rev_file_name'].values.tolist())
        return list_of_meta_gz_files

    def create_sample_meta_info_dataframe_from_datasheet_path(self):
        if self.datasheet_path.endswith('.xlsx'):
            self.sample_meta_info_df = pd.read_excel(io=self.datasheet_path, header=0, index_col=0, usecols='A:N', skiprows=[0])
        elif self.datasheet_path.endswith('.csv'):
            self.sample_meta_info_df = pd.read_csv(filepath_or_buffer=self.datasheet_path, header=0, index_col=0, skiprows=[0])
        else:
            sys.exit('Data sheet: {} is in an unrecognised format. '
                     'Please ensure that it is either in .xlsx or .csv format.')


    def copy_and_decompress_input_files_to_temp_wkd(self):
        if not self.is_single_file_or_paired_input:
            self.copy_fastq_files_from_input_dir_to_temp_wkd()
        else:
            self.extract_single_compressed_file_to_temp_wkd()

    def extract_single_compressed_file_to_temp_wkd(self):
        ext_components = self.user_input_path.split('.')
        if ext_components[-1] == 'zip':  # .zip
            subprocess.run(["unzip", self.user_input_path, '-d', self.temp_working_directory], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        elif ext_components[-2] == 'tar' and ext_components[-1] == 'gz':  # .tar.gz
            subprocess.run(["tar", "-xf", self.user_input_path, "-C", self.temp_working_directory],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        elif ext_components[-1] == 'gz':  # .gz
            subprocess.run(["gunzip", "-c", self.user_input_path, ">", self.temp_working_directory],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)

    def copy_fastq_files_from_input_dir_to_temp_wkd(self):
        file = os.listdir(self.user_input_path)[0]
        os.chdir('{}'.format(self.user_input_path))
        # * asterix are only expanded in the shell and so don't work through subprocess
        # need to use the glob library instead
        # https://stackoverflow.com/questions/13875978/python-subprocess-popen-why-does-ls-txt-not-work
        # files could be compressed (fastq.gz) or uncompressed (fastq). Either way, they should contain fastq.
        if 'fastq' in file:
            subprocess.run(['cp'] + glob.glob('*.fastq*') + [self.temp_working_directory], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        elif 'fq' in file:
            subprocess.run(['cp'] + glob.glob('*.fq*') + [self.temp_working_directory], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)

    def determine_if_single_file_or_paired_input(self):
        for file in os.listdir(self.user_input_path):
            if 'fastq' in file or 'fq' in file:
                # Then there is a fastq.gz or fastq file already uncompressed in this folder
                # In this case we will assume that the seq data is not a single file containing the pairs of files
                # rather the pairs of files themselves.
                return False
    def create_pre_med_write_out_directory_path(self):
        pre_med_write_out_directory_path = self.temp_working_directory.replace('tempData', 'pre_MED_seqs')
        os.makedirs(pre_med_write_out_directory_path, exist_ok=True)
        return pre_med_write_out_directory_path

    def setup_output_directory(self):
        output_directory = os.path.join(os.path.dirname(__file__),
                                        'outputs/data_set_submissions/{}'.format(self.dataset_object.id))
        os.makedirs(output_directory, exist_ok=True)
        return output_directory

    def set_temp_working_directory(self):
        # working directory will be housed in a temp folder within the directory in which the sequencing data
        # is currently housed
        if '.' in self.user_input_path.split('/')[-1]:
            # then this path points to a file rather than a directory and we should pass through the path only
            self.temp_working_directory = os.path.abspath(
                '{}/tempData/{}'.format(os.path.dirname(self.user_input_path), self.dataset_object.id))
        else:
            # then we assume that we are pointing to a directory and we can directly use that to make the wkd
            self.temp_working_directory = os.path.abspath('{}/tempData/{}'.format(self.user_input_path, self.dataset_object.id))
        self.dataset_object.working_directory = self.temp_working_directory
        self.dataset_object.save()


    def create_temp_wkd(self):
        # if the directory already exists remove it and start from scratch
        if os.path.exists(self.temp_working_directory):
            shutil.rmtree(self.temp_working_directory)
        os.makedirs(self.temp_working_directory)