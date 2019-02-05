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
from analysis_classes import InitialMothurWorker

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
        self.screen_sub_evalue = screen_sub_evalue
        self.samples_that_caused_errors_in_qc_list = []

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

        :param input_q:
        :param error_sample_list:
        :return:
        """

        for contigpair in iter(input_q.get, 'STOP'):

            initial_morthur_worker = InitialMothurWorker(
                contig_pair=contigpair, data_set_object=self.dataset_object)

            initial_morthur_worker.execute()



            # todo I have written in a make contig method into the Mothur analysis. Incorporate this here.





            # NB We will always crop with the SYMVAR primers as they produce the shortest product
            #todo aim to use the MothurAnalysis class that I made

            # Make the sample by sample directory that we will be working in
            # this will be inside the wkd directory (the temp data directory for the data_set submission)
            # We also need to make the same sample by sample directories for the pre MED sequence dump


            # stability_file = [contigPair]
            # stability_file_name = r'{0}{1}'.format(sample_name, 'stability.files')
            # root_name = r'{0}stability'.format(sample_name)
            # stability_file_path = r'{0}{1}'.format(current_directory, stability_file_name)
            #
            # # write out the stability file. This will be a single pair of contigs with a sample name
            # write_list_to_destination(stability_file_path, stability_file)



            # NB mothur is working very strangely with the python subprocess command. For some
            # reason it is adding in an extra 'mothur' before the filename in the input directory
            # As such we will have to enter all of the paths to files absolutely

            # The mothur batch file that will be run by mothur.
            mothur_batch_file = [
                # r'set.dir(input={0})'.format(current_directory),
                # r'set.dir(output={0})'.format(current_directory),
                # r'make.contigs(file={}{})'.format(current_directory, stability_file_name),
                # r'summary.seqs(fasta={}{}.trim.contigs.fasta)'.format(current_directory, root_name),
                # r'screen.seqs(fasta={0}{1}.trim.contigs.fasta, group={0}{1}.contigs.groups, '
                # r'maxambig=0, maxhomop=5)'.format(current_directory, root_name),
                # r'summary.seqs(fasta={0}{1}.trim.contigs.good.fasta)'.format(current_directory, root_name),
                # r'unique.seqs(fasta={0}{1}.trim.contigs.good.fasta)'.format(current_directory, root_name),
                # r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.fasta, '
                # r'name={0}{1}.trim.contigs.good.names)'.format(current_directory, root_name),
                # r'split.abund(cutoff=2, fasta={0}{1}.trim.contigs.good.unique.fasta, '
                # r'name={0}{1}.trim.contigs.good.names, group={0}{1}.contigs.good.groups)'.format(
                #     current_directory,
                #     root_name),
                # r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.fasta, '
                # r'name={0}{1}.trim.contigs.good.abund.names)'.format(current_directory, root_name),
                # r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.rare.fasta, '
                # r'name={0}{1}.trim.contigs.good.rare.names)'.format(current_directory, root_name),
                r'pcr.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.fasta, group={0}{1}.contigs.good.abund.groups, '
                r'name={0}{1}.trim.contigs.good.abund.names, '
                r'oligos={0}primers.oligos, pdiffs=2, rdiffs=2)'.format(current_directory, root_name)
            ]

            # Write out the batch file
            mothur_batch_file_path = r'{0}{1}{2}'.format(current_directory, 'mothur_batch_file', sample_name)
            write_list_to_destination(mothur_batch_file_path, mothur_batch_file)

            error = False
            # NB the mothur return code doesn't seem to work. We just get None type.
            # apparently they have fixed this in the newest mothur but we have not upgraded to that yet.
            # so for the time being we will check for error by hand in the stdout.
            with subprocess.Popen(['mothur', '{0}'.format(mothur_batch_file_path)], stdout=subprocess.PIPE,
                                  bufsize=1,
                                  universal_newlines=True) as p:
                # Here look for the specific blank fasta name warning (which should be interpreted as an error)
                # and any other error that may be arising
                # if found, log error.
                for line in p.stdout:
                    if debug:
                        print(line)
                    if '[WARNING]: Blank fasta name, ignoring read.' in line:
                        p.terminate()
                        error_reason = 'Blank fasta name'
                        log_qc_error_and_continue(data_set_sample_instance_in_q, sample_name, error_reason)
                        error = True
                        error_sample_list.append(sample_name)
                        break
                    if 'ERROR' in line:
                        p.terminate()
                        error_reason = 'error in inital QC'
                        log_qc_error_and_continue(data_set_sample_instance_in_q, sample_name, error_reason)
                        error = True
                        error_sample_list.append(sample_name)
                        break

            if error:
                continue

            # Here check the outputted files to see if they are reverse complement
            # or not by running the pcr.seqs and checking the results
            # Check to see if there are sequences in the PCR output file
            last_summary = read_defined_file_to_list(
                '{}{}.trim.contigs.good.unique.abund.pcr.fasta'.format(current_directory, root_name))

            # If this file is empty
            #  Then these sequences may well be reverse complement so we need to try to rev first
            if len(last_summary) == 0:

                # RC batch file
                mothur_batch_reverse = [
                    r'reverse.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.fasta)'.format(
                        current_directory, root_name),
                    r'pcr.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.rc.fasta, '
                    r'group={0}{1}.contigs.good.abund.groups, name={0}{1}.trim.contigs.good.abund.names, '
                    r'oligos={0}primers.oligos, pdiffs=2, rdiffs=2)'.format(current_directory, root_name)
                ]
                mothur_batch_file_path = r'{0}{1}{2}'.format(current_directory, 'mothur_batch_file', sample_name)
                # write out RC batch file
                write_list_to_destination(mothur_batch_file_path, mothur_batch_reverse)

                if not debug:
                    subprocess.run(
                        ['mothur', r'{0}'.format(mothur_batch_file_path)],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    # At this point the sequences will be reversed and they will have been renamed so we
                    # can just change the name of the .rc file to the orignal .fasta file that we inputted with
                    # This way we don't need to change the rest of the mothur pipe.
                    subprocess.run(
                        [r'mv', r'{0}{1}.trim.contigs.good.unique.abund.rc.pcr.fasta'.format(current_directory,
                                                                                             root_name),
                         r'{0}{1}.trim.contigs.good.unique.abund.pcr.fasta'.format(current_directory, root_name)],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                elif debug:
                    subprocess.run(
                        ['mothur', r'{0}'.format(mothur_batch_file_path)])
                    subprocess.run(
                        [r'mv', r'{0}{1}.trim.contigs.good.unique.abund.rc.pcr.fasta'.format(current_directory,
                                                                                             root_name),
                         r'{0}{1}.trim.contigs.good.unique.abund.pcr.fasta'.format(current_directory, root_name)])

            # Check again to see if the RC has fixed the problem of having an empty fasta
            # If this file is still empty, then the problem was not solved by reverse complementing
            last_summary = read_defined_file_to_list(
                '{}{}.trim.contigs.good.unique.abund.pcr.fasta'.format(current_directory, root_name))

            if len(last_summary) == 0:
                error_reason = 'error in inital QC'
                log_qc_error_and_continue(data_set_sample_instance_in_q, sample_name, error_reason)
                error_sample_list.append(sample_name)
                continue

            # after having completed the RC checks redo the unique.
            mothur_batch_file_cont = [
                r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.pcr.fasta, '
                r'name={0}{1}.trim.contigs.good.abund.pcr.names)'.format(current_directory, root_name),
                r'unique.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.pcr.fasta, '
                r'name={0}{1}.trim.contigs.good.abund.pcr.names)'.format(current_directory, root_name),
                r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.pcr.unique.fasta, '
                r'name={0}{1}.trim.contigs.good.unique.abund.pcr.names)'.format(current_directory, root_name)
            ]

            mothur_batch_file_path = r'{0}{1}{2}'.format(current_directory, 'mothur_batch_file', sample_name)
            write_list_to_destination(mothur_batch_file_path, mothur_batch_file_cont)

            completed_process = subprocess.run(
                ['mothur', r'{0}'.format(mothur_batch_file_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            if completed_process.returncode == 1 or 'ERROR' in completed_process.stdout.decode('utf-8'):
                if debug:
                    print(completed_process.stdout.decode('utf-8'))
                error_reason = 'error in inital QC'
                log_qc_error_and_continue(data_set_sample_instance_in_q, sample_name, error_reason)
                error_sample_list.append(sample_name)
                continue

            # Check to see if there are sequences in the PCR output file
            try:
                last_summary = read_defined_file_to_list(
                    '{}{}.trim.contigs.good.unique.abund.pcr.unique.fasta'.format(current_directory, root_name))
                if len(last_summary) == 0:  # If this file is empty
                    error_reason = 'error in inital QC'
                    log_qc_error_and_continue(data_set_sample_instance_in_q, sample_name, error_reason)
                    error_sample_list.append(sample_name)
                    continue
            except FileNotFoundError:  # If there is no file then we can assume sample has a problem
                log_qc_error_and_continue(data_set_sample_instance_in_q, sample_name, 'generic_error')
                continue

            # Get unique number of sequences after after sequence QC
            last_summary = read_defined_file_to_list(
                '{}{}.trim.contigs.good.unique.abund.pcr.unique.summary'.format(current_directory, root_name))
            number_of_seqs_contig_unique = len(last_summary) - 1
            data_set_sample_instance_in_q.post_qc_unique_num_seqs = number_of_seqs_contig_unique
            sys.stdout.write(
                '{}: data_set_sample_instance_in_q.post_qc_unique_num_seqs = {}\n'.format(
                    sample_name, number_of_seqs_contig_unique))

            # Get absolute number of sequences after after sequence QC
            last_summary = read_defined_file_to_list(
                '{}{}.trim.contigs.good.unique.abund.pcr.unique.summary'.format(current_directory, root_name))
            absolute_count = 0
            for line in last_summary[1:]:
                absolute_count += int(line.split('\t')[6])
            data_set_sample_instance_in_q.post_qc_absolute_num_seqs = absolute_count
            data_set_sample_instance_in_q.save()
            sys.stdout.write('{}: data_set_sample_instance_in_q.post_qc_absolute_num_seqs = {}\n'.format(
                sample_name, absolute_count))

            sys.stdout.write('{}: Initial mothur complete\n'.format(sample_name))
            # Each sampleDataDir should contain a set of .fasta, .name and .group
            # files that we can use to do local blasts with

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