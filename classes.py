from dbApp.models import (DataSet, ReferenceSequence, DataSetSample)
import sys
import os
import shutil
import subprocess
import glob
import pandas as pd
import general
import json
from analysis_classes import (
    InitialMothurHandler, PotentialSymTaxScreeningHandler, BlastnAnalysis, SymNonSymTaxScreeningHandler,
    PerformMEDHandler, DataSetSampleCreatorHandler, SequenceCountTableCreator, SeqStackedBarPlotter)
from general import write_list_to_destination, read_defined_file_to_list, create_dict_from_fasta, make_new_blast_db
from datetime import datetime


class DataLoading:
    # The clades refer to the phylogenetic divisions of the Symbiodiniaceae. Most of them are represented at the genera
    # level. E.g. All members of clade C belong to the genus Cladocopium.
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    def __init__(
            self, data_set_uid, user_input_path, datasheet_path, screen_sub_evalue=False,
            debug=False, num_proc=1, no_fig=False, no_ord=False, distance_method='braycurtis'):
        self.symportal_root_directory = os.path.abspath(os.path.dirname(__file__))
        self.output_path_list = []
        self.no_fig = no_fig
        self.no_ord = no_ord
        self.distance_method = distance_method
        self.dataset_object = DataSet.objects.get(id=data_set_uid)
        self.user_input_path = user_input_path
        self.output_directory = self._setup_output_directory()
        self.temp_working_directory = self._setup_temp_working_directory()
        self.data_time_string = str(datetime.now()).replace(' ', '_').replace(':', '-')
        # this is the path of the file we will use to deposit a backup copy of the reference sequences
        self.seq_dump_file_path = self._setup_sequence_dump_file_path()
        self.dataset_object.working_directory = self.temp_working_directory
        self.dataset_object.save()
        # This is the directory that sequences that have undergone QC for each sample will be written out as
        # .names and .fasta pairs BEFORE MED decomposition
        self.pre_med_sequence_output_directory_path = self._create_pre_med_write_out_directory_path()
        self.num_proc = num_proc
        # directory that will contain sub directories for each sample. Each sub directory will contain a pair of
        # .names and .fasta files of the non_symbiodinium_sequences that were thrown out for that sample
        self.non_symb_and_size_violation_base_dir_path = os.path.join(
            self.output_directory, 'non_symbiodiniaceae_and_size_violation_sequences'
        )
        os.makedirs(self.non_symb_and_size_violation_base_dir_path, exist_ok=True)
        # data can be loaded either as paired fastq or fastq.gz files or as a single compressed file containing
        # the afore mentioned paired files.
        self.is_single_file_or_paired_input = self._determine_if_single_file_or_paired_input()
        self.datasheet_path = datasheet_path
        self.debug = debug
        self.symclade_db_directory_path = os.path.abspath(os.path.join(self.symportal_root_directory, 'symbiodiniumDB'))
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
        self.initial_mothur_handler = None
        self.post_initial_qc_name_file_name = None
        self.post_initial_qc_fasta_file_name = None
        # args for the taxonomic screening
        self.screen_sub_evalue = screen_sub_evalue
        self.new_seqs_added_in_iteration = 0
        self.new_seqs_added_running_total = 0
        self.checked_samples_with_no_additional_symbiodinium_sequences = []
        self.taxonomic_screening_handler = None
        self.sequences_to_screen_fasta_as_list = []
        self.sequences_to_screen_fasta_path = os.path.join(
            self.temp_working_directory, 'taxa_screening_seqs_to_screen.fa'
        )
        # file names used for writing out
        self.non_sym_fasta_file_name_str = None
        self.non_sym_name_file_name_str = None
        self.sym_fasta_file_name_str = None
        self.sym_name_file_name_str = None
        self.sym_binary_clade_dict_name_str = None
        # the number of samples that a sub evalue sequences must be
        # found in for us to carry it through for taxonomic screening
        self.required_sample_support_for_sub_evalue_sequencs = 3
        # the number of sequences in the 10 matches that must be annotated as Symbiodinium or Symbiodiniaceae in order
        # for a sequences to be added into the reference symClade database.
        self.required_symbiodinium_matches = 3
        # med
        self.list_of_med_output_directories = []
        self.path_to_med_padding_executable = os.path.join(
            self.symportal_root_directory, 'lib/med_decompose/o_pad_with_gaps.py')
        self.path_to_med_decompose_executable = os.path.join(
            self.symportal_root_directory, 'lib/med_decompose/decompose.py')
        self.perform_med_handler_instance = None
        # data set sample creation
        self.data_set_sample_creator_handler_instance = None
        # plotting sequence output
        self.seq_abundance_relative_output_path = None

    def load_data(self):
        self._copy_and_decompress_input_files_to_temp_wkd()

        # the stability file generated here is used as the base of the initial mothur QC
        if self.datasheet_path:
            self._generate_stability_file_and_data_set_sample_objects_with_datasheet()
        else:
            self._generate_stability_file_and_data_set_sample_objects_without_datasheet()

        self._if_symclade_binaries_not_present_remake_db()

        self._do_initial_mothur_qc()

        self._taxonomic_screening()

        self._do_med_decomposition()

        self._create_data_set_sample_sequences_from_med_nodes()

        self._print_sample_successful_or_failed_summary()

        self._perform_sequence_drop()

        self._delete_temp_working_directory_and_log_files()

        self._output_seqs_count_table()

        self._output_seqs_stacked_bar_plot()

        if not self.no_ord:
            print('Calculating between sample pairwise distances')
            pcoa_paths_list = None
            if distance_method == 'unifrac':
                pcoa_paths_list = generate_within_clade_unifrac_distances_samples(
                    data_set_string=data_set_uid,
                    num_processors=num_proc,
                    method='mothur', call_type='submission',
                    date_time_string=date_time_str,
                    output_dir=output_directory)
            elif distance_method == 'braycurtis':
                pcoa_paths_list = generate_within_clade_braycurtis_distances_samples(
                    data_set_string=data_set_uid,
                    call_type='submission',
                    date_time_str=date_time_str,
                    output_dir=output_directory)
            # distance plotting
            output_path_list.extend(pcoa_paths_list)
            if not no_fig:
                if num_samples > 1000:
                    print('Too many samples ({}) to generate plots'.format(num_samples))
                else:
                    for pcoa_path in pcoa_paths_list:
                        if 'PCoA_coords' in pcoa_path:
                            # then this is a full path to one of the .csv files that contains the
                            # coordinates that we can plot we will get the output directory from the passed in pcoa_path
                            sys.stdout.write('\n\nGenerating between sample distance plot clade {}\n'.format(
                                os.path.dirname(pcoa_path).split('/')[-1]))
                            ordination_figure_output_paths_list = plot_between_sample_distance_scatter(
                                csv_path=pcoa_path,
                                date_time_str=date_time_str
                            )
                            output_path_list.extend(ordination_figure_output_paths_list)

    def _output_seqs_stacked_bar_plot(self):
        if not self.no_fig:
            if self.num_of_samples > 1000:
                print(f'Too many samples ({num_samples}) to generate plots')
            else:
                sys.stdout.write('\nGenerating sequence count table figures\n')

                seq_stacked_bar_plotter = SeqStackedBarPlotter(output_directory=self.output_directory,
                                                               seq_relative_abund_count_table_path=self.seq_abundance_relative_output_path)
                seq_stacked_bar_plotter.plot_stacked_bar_seqs()
                self.output_path_list.extend(seq_stacked_bar_plotter.output_path_list)
                # TODO don't for get to add the output path for the non-sym and size violation output
                sys.stdout.write('\nFigure generation complete')
                sys.stdout.write('\nFigures output to:')
                for path in seq_stacked_bar_plotter.output_path_list:
                    sys.stdout.write(f'\n{path}')

    def _output_seqs_count_table(self):
        sys.stdout.write('\nGenerating count tables\n')
        sequence_count_table_creator = SequenceCountTableCreator(
            call_type='submission',
            output_dir=self.output_directory,
            data_set_uids_to_output_as_comma_sep_string=set(self.dataset_object.id),
            num_proc=self.num_proc, time_date_str=self.data_time_string)
        sequence_count_table_creator.execute_output()
        # TODO don't for get to write out where the non-sym and size violation seqs were output
        self.output_directory.extend(sequence_count_table_creator.output_paths_list)
        self._set_seq_abundance_relative_output_path(sequence_count_table_creator)

    def _set_seq_abundance_relative_output_path(self, sequence_count_table_creator):
        for path in sequence_count_table_creator.output_paths_list:
            if 'relative' in path:
                self.seq_abundance_relative_output_path = path

    def _delete_temp_working_directory_and_log_files(self):
        if os.path.exists(self.temp_working_directory):
            shutil.rmtree(self.temp_working_directory)
        # del any log file in the database directory that may have been created by mothur analyses
        log_file_list = [f for f in os.listdir(os.path.dirname(self.symclade_db_full_path)) if f.endswith(".logfile")]
        for f in log_file_list:
            os.remove(os.path.join(os.path.dirname(self.symclade_db_full_path), f))

    def _perform_sequence_drop(self):
        sequence_drop_file = self._generate_sequence_drop_file()
        sys.stdout.write(f'\n\nBackup of named reference_sequences output to {self.seq_dump_file_path}\n')
        write_list_to_destination(self.seq_dump_file_path, sequence_drop_file)

    @staticmethod
    def _generate_sequence_drop_file():
        header_string = ','.join(['seq_name', 'seq_uid', 'seq_clade', 'seq_sequence'])
        sequence_drop_list = [header_string]
        for ref_seq in ReferenceSequence.objects.filter(has_name=True):
            sequence_drop_list.append(','.join([ref_seq.name, str(ref_seq.id), ref_seq.clade, ref_seq.sequence]))
        return sequence_drop_list

    def _print_sample_successful_or_failed_summary(self):
        print(f'\n\n SAMPLE PROCESSING SUMMARY')
        failed_count = 0
        successful_count = 0
        for data_set_sample in DataSetSample.objects.filter(data_submission_from=self.dataset_object):
            if data_set_sample.error_in_processing:
                failed_count += 1
                print(f'{data_set_sample.name}: Error in processing: {data_set_sample.error_reason}')
            else:
                successful_count += 1
                print(f'{data_set_sample.name}: Successful')

        print(f'\n\n{successful_count} out of {successful_count + failed_count} samples successfully passed QC.\n'
              f'{failed_count} samples produced erorrs\n')

    def _create_data_set_sample_sequences_from_med_nodes(self):
        self.data_set_sample_creator_handler_instance = DataSetSampleCreatorHandler()
        self.data_set_sample_creator_handler_instance.execute_data_set_sample_creation(
            data_loading_list_of_med_output_directories=self.list_of_med_output_directories,
            data_loading_debug=self.debug, data_loading_dataset_object=self.dataset_object)
        self.dataset_object.currently_being_processed = False
        self.dataset_object.save()

    def _do_med_decomposition(self):
        self.perform_med_handler_instance = PerformMEDHandler(
            data_loading_temp_working_directory=self.temp_working_directory,
            data_loading_num_proc=self.num_proc,
            data_loading_list_of_samples_names=self.list_of_samples_names)

        self.perform_med_handler_instance.execute_perform_med_worker(
            data_loading_debug=self.debug,
            data_loading_path_to_med_decompoase_executable=self.path_to_med_decompose_executable,
            data_loading_path_to_med_padding_executable=self.path_to_med_padding_executable)

        self.list_of_med_output_directories = self.perform_med_handler_instance.list_of_med_result_dirs

        if self.debug:
            print('MED dirs:')
            for med_output_directory in self.list_of_med_output_directories:
                print(med_output_directory)

    def _do_initial_mothur_qc(self):

        if not self.sample_fastq_pairs:
            self._exit_and_del_data_set_sample('Sample fastq pairs list empty')

        self.initial_mothur_handler = InitialMothurHandler(
            num_proc=self.num_proc, sample_fastq_pairs=self.sample_fastq_pairs)

        self.initial_mothur_handler.execute_worker_initial_mothur(
            data_loading_dataset_object=self.dataset_object,
            data_loading_temp_working_directory=self.temp_working_directory,
            data_loading_debug=self.debug
        )

        self.samples_that_caused_errors_in_qc_list = list(
            self.initial_mothur_handler.samples_that_caused_errors_in_qc_mp_list)

    def _taxonomic_screening(self):
        """
        There are only two things to achieve as part of this taxonomy screening.
        1 - identify the sequences that are Symbiodinium/Symbiodiniaceae in origin
        2 - identify the sequences that are non-Symbiodinium/non-Symbiodiniaceae in origin.
        This may sound very straight forwards: Blast each sequence against a reference database and if we get
        a good enough match, consider this sequence symbiodinium in origin. If not, not.
        BUT, there is a catch. We cannot be sure that every new sequence that we receive is not symbiodinium just
        because we don't get a good match to our reference database. It may be that new diversity is not yet
        represented in our reference database.
        As a results of this we want to do an extra set of taxonomic screening and that is what this first part of code
        concerns. We will run a blast against our reference symClade database. And then, any seuqences that return
        a match to a member of this database, but does not meet the minimum threshold to be directly considered
        Symbiodinium in origin (i.e. not similar enough to a reference sequence in the symClade database) will
        be blasted against the NCBI nt blast database. For a sequence to be blasted against this database, and be
        in contesion for being considered symbiodinium in origin it must also be found in at least three samples.
        We will use an iterative screening process to acheive this symbiodinium identification. The sequences to be
        screened will be referred to as subevalue sequences. Any of these sequences that we deem symbiodinium in
        origin after running against the nt database will be added back into the symClade reference database.
        On the next iteration it will therefore be possible for differnt sequences to be matches given that additional
        sequences may have been added to this reference database. This first part of screening will only be to run
        the symClade blast, and to screen low identity matches. Non-symbiodinium sequence matches will be made in
        the second part of this screening.
        To make this faster, rather than check every sample in every iteration we will keep track of when a sample
        has been found to only contain sequences that are good matches to the symClade database.
        This way we can skip these sequences in the iterative rounds of screening.
        """

        if self.screen_sub_evalue:

            self._create_symclade_backup_incase_of_accidental_deletion_of_corruption()

            while 1:
                self.new_seqs_added_in_iteration = 0
                # This method simply identifies whether there are sequences that need screening.
                # Everytime that the execute_worker is run it pickles out the files needed for the next worker
                self._make_fasta_of_sequences_that_need_taxa_screening()

                if self.sequences_to_screen_fasta_as_list:
                    # Now do the screening
                    # The outcome of this will be an updated symClade.fa that we should then make a blastdb from it.
                    self._screen_sub_e_seqs()
                    if self.new_seqs_added_in_iteration == 0:
                        break
                else:
                    break

        else:
            # if not doing the screening we can simply run the execute_worker_taxa_screening once.
            # During its run it will have output all of the files we need to run the following workers.
            # We can also run the generate_and_write_below_evalue_fasta_for_screening function to write out a
            # fasta of sub_evalue sequences we can then report to the user using that object
            self._make_fasta_of_sequences_that_need_taxa_screening()

        self._do_sym_non_sym_tax_screening()

    def _do_sym_non_sym_tax_screening(self):
        self.sym_non_sym_tax_screening_handler = SymNonSymTaxScreeningHandler(
            data_loading_list_of_samples_names=self.list_of_samples_names,
            data_loading_num_proc=self.num_proc,
            data_loading_samples_that_caused_errors_in_qc_mp_list=self.samples_that_caused_errors_in_qc_list
        )
        self.sym_non_sym_tax_screening_handler.execute_sym_non_sym_tax_screening(
            data_loading_temp_working_directory=self.temp_working_directory,
            data_loading_pre_med_sequence_output_directory_path=self.pre_med_sequence_output_directory_path,
            data_loading_dataset_object=self.dataset_object,
            non_symb_and_size_violation_base_dir_path=self.non_symb_and_size_violation_base_dir_path,
            data_loading_debug=self.debug
        )
        self.samples_that_caused_errors_in_qc_list = list(
            self.sym_non_sym_tax_screening_handler.samples_that_caused_errors_in_qc_mp_list
        )

    def _screen_sub_e_seqs(self):
        """This function screens a fasta file to see if the sequences are Symbiodinium in origin.
        This fasta file contains the below_e_cutoff sequences that need to be screened.
        These sequences are the sequences
        that were found to have a match in the initial screening against the symClade.fa database but were below
        the required evalue cut off. Here we run these sequences against the entire NCBI 'nt' database to verify if they
        or of Symbiodinium origin of not.
        The fasta we are screening only contains seuqences that were found in at least 3 samples.
        We will call a sequence symbiodinium if it has a match that covers at least 95% of its sequence
        at a 60% or higher
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

        query_sequences_verified_as_symbiodinium_list = self._get_list_of_seqs_in_blast_result_that_are_symbiodinium(
            blast_output_dict
        )

        self.new_seqs_added_in_iteration = len(query_sequences_verified_as_symbiodinium_list)
        self.new_seqs_added_running_total += self.new_seqs_added_in_iteration

        if query_sequences_verified_as_symbiodinium_list:
            self._taxa_screening_update_symclade_db_with_new_symbiodinium_seqs(
                query_sequences_verified_as_symbiodinium_list)

    def _taxa_screening_update_symclade_db_with_new_symbiodinium_seqs(
            self, query_sequences_verified_as_symbiodinium_list):
        new_symclade_fasta_as_list = self._taxa_screening_make_new_fasta_of_screened_seqs_to_be_added_to_symclade_db(
            query_sequences_verified_as_symbiodinium_list
        )
        combined_fasta = self._taxa_screening_combine_new_symclade_seqs_with_current(new_symclade_fasta_as_list)
        self._taxa_screening_make_new_symclade_db(combined_fasta)

    def _taxa_screening_make_new_symclade_db(self, combined_fasta):
        write_list_to_destination(self.symclade_db_full_path, combined_fasta)
        make_new_blast_db(input_fasta_to_make_db_from=self.symclade_db_full_path, db_title='symClade')

    def _taxa_screening_combine_new_symclade_seqs_with_current(self, new_symclade_fasta_as_list):
        old_symclade_fasta_as_list = read_defined_file_to_list(self.symclade_db_full_path)
        combined_fasta = new_symclade_fasta_as_list + old_symclade_fasta_as_list
        return combined_fasta

    def _taxa_screening_make_new_fasta_of_screened_seqs_to_be_added_to_symclade_db(
            self, query_sequences_verified_as_symbiodinium_list):
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

    def _get_list_of_seqs_in_blast_result_that_are_symbiodinium(self, blast_output_dict):
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

    def _make_fasta_of_seqs_found_in_more_than_two_samples_that_need_screening(self):
        """ The below_e_cutoff_dict has nucleotide sequencs as the
        key and the number of samples that sequences was found in as the value.
        """
        sub_evalue_nuclotide_sequence_to_number_of_samples_found_in_dict = dict(
            self.taxonomic_screening_handler.sub_evalue_sequence_to_num_sampes_found_in_mp_dict)
        self.sequences_to_screen_fasta_as_list = []
        sequence_number_counter = 0
        for nucleotide_sequence, num_samples_found_in in \
                sub_evalue_nuclotide_sequence_to_number_of_samples_found_in_dict.items():
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

    def _create_symclade_backup_incase_of_accidental_deletion_of_corruption(self):
        back_up_dir = os.path.abspath(os.path.join(self.symportal_root_directory, 'symbiodiniumDB', 'symClade_backup'))
        os.makedirs(back_up_dir, exist_ok=True)
        src_path = os.path.abspath(os.path.join(self.symportal_root_directory, 'symbiodiniumDB')) + '/symClade.fa'
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

    def _make_fasta_of_sequences_that_need_taxa_screening(self):
        self._init_potential_sym_tax_screen_handler()

        # the self.taxonomic_screening_handler.sub_evalue_sequence_to_num_sampes_found_in_mp_dict is populated here
        self.taxonomic_screening_handler.execute_potential_sym_tax_screening(
            data_loading_temp_working_directory=self.temp_working_directory,
            data_loading_path_to_symclade_db=self.symclade_db_full_path,
            data_loading_debug=self.debug
        )

        self._taxa_screening_update_checked_samples_list()

        self._make_fasta_of_seqs_found_in_more_than_two_samples_that_need_screening()

    def _taxa_screening_update_checked_samples_list(self):
        self.checked_samples_with_no_additional_symbiodinium_sequences = \
            list(self.taxonomic_screening_handler.checked_samples_mp_list)

    def _init_potential_sym_tax_screen_handler(self):

        self.taxonomic_screening_handler = PotentialSymTaxScreeningHandler(
            samples_that_caused_errors_in_qc_list=self.samples_that_caused_errors_in_qc_list,
            checked_samples_list=self.checked_samples_with_no_additional_symbiodinium_sequences,
            list_of_samples_names=self.list_of_samples_names, num_proc=self.num_proc
        )

    def _if_symclade_binaries_not_present_remake_db(self):
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
                self._exit_and_del_data_set_sample('Failure in creating blast binaries')

    def _generate_stability_file_and_data_set_sample_objects_without_datasheet(self):

        self.list_of_fastq_file_names_in_wkd = [a for a in os.listdir(self.temp_working_directory) if 'fastq' in a]

        self._generate_and_write_mothur_batch_file_for_dotfile_creation()
        general.execute_mothur_batch_file_with_piped_stoud_sterr(self.path_to_latest_mothur_batch_file)

        self._generate_and_write_new_stability_file_without_datasheet()

        self._create_data_set_sample_objects_in_bulk_without_datasheet()

    def _create_data_set_sample_objects_in_bulk_without_datasheet(self):
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

    def _generate_and_write_new_stability_file_without_datasheet(self):
        self._read_in_mothur_dot_file_creation_output()
        # We will use the stability file in the rest of the mothur qc but also to make the DataSetSamples (below)
        end_index = self._identify_sample_names_without_datasheet()
        new_stability_file = self._generate_new_stability_file_without_datasheet(end_index)
        # write out the new stability file
        general.write_list_to_destination(f'{self.temp_working_directory}/stability.files', new_stability_file)

    def _generate_new_stability_file_without_datasheet(self, end_index):
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

    def _get_num_chars_in_common_with_fastq_names(self):
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

    def _get_sample_names_from_fastq_files_using_index(self, end_index):
        list_of_names_non_unique = []
        for file in self.list_of_fastq_files_in_wkd:
            list_of_names_non_unique.append(file[:-end_index])
        list_of_sample_names = list(set(list_of_names_non_unique))
        if len(list_of_sample_names) != len(self.list_of_fastq_files_in_wkd) / 2:
            warning_str = 'Error in sample name extraction'
            self._exit_and_del_data_set_sample(warning_str)
        self.list_of_samples_names = list_of_sample_names

    def _identify_sample_names_without_datasheet(self):
        # I think the simplest way to get sample names is to find what parts are common between all samples
        # well actually 50% of the samples so that we also remove the R1 and R2 parts.
        end_index = self._get_num_chars_in_common_with_fastq_names()
        self._get_sample_names_from_fastq_files_using_index(end_index)

        return end_index

    def _generate_stability_file_and_data_set_sample_objects_with_datasheet(self):
        # Create a pandas df from the data_sheet if it was provided
        # allow the data_sheet to be in a .csv format or .xlsx format. This is so that we can store a datasheet
        # in the github repo in a non-binary format
        # The sample_meta_df that is created from the data_sheet should be identical irrespective of whether a .csv
        # or a .xlsx is submitted.
        self._create_sample_meta_info_dataframe_from_datasheet_path()

        # if we are given a data_sheet then use the sample names given as the DataSetSample object names
        self.list_of_samples_names = self.sample_meta_info_df.index.values.tolist()

        # we should verify that all of the fastq files listed in the sample_meta_df
        # are indeed found in the directory that we've been given
        self._check_all_fastqs_in_datasheet_exist()

        self._generate_and_write_mothur_batch_file_for_dotfile_creation()

        # noinspection PyPep8
        general.execute_mothur_batch_file_with_piped_stoud_sterr(self.path_to_latest_mothur_batch_file)

        # we will also need to know how to relate the sample names to the fastq files
        # for this we will make a dict of fastq file name to sample
        self._create_fastq_file_to_sample_name_dict()
        # We will use the stability file in the rest of the mothur qc but also to make the DataSetSamples (below)
        self._generate_and_write_new_stability_file_with_data_sheet()

        self._create_data_set_sample_objects_in_bulk_with_datasheet()

    def _create_data_set_sample_objects_in_bulk_with_datasheet(self):
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

    def _generate_and_write_new_stability_file_with_data_sheet(self):
        # Convert the group names in the stability.files so that the dashes are converted to '[ds]',
        # So for the mothur we have '[ds]'s. But for all else we convert these '[ds]'s to dashes
        self._read_in_mothur_dot_file_creation_output()

        new_stability_file = self._generate_new_stability_file_with_data_sheet()

        general.write_list_to_destination(
            r'{0}/stability.files'.format(self.temp_working_directory),
            new_stability_file
        )

    def _generate_new_stability_file_with_data_sheet(self):
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

    def _read_in_mothur_dot_file_creation_output(self):
        self.sample_fastq_pairs = general.read_defined_file_to_list(f'{self.temp_working_directory}/stability.files')

    def _create_fastq_file_to_sample_name_dict(self):
        fastq_file_to_sample_name_dict = {}
        for sample_index in self.sample_meta_info_df.index.values.tolist():
            fastq_file_to_sample_name_dict[
                self.sample_meta_info_df.loc[sample_index, 'fastq_fwd_file_name']] = sample_index
            fastq_file_to_sample_name_dict[
                self.sample_meta_info_df.loc[sample_index, 'fastq_rev_file_name']] = sample_index
        self.fastq_file_to_sample_name_dict = fastq_file_to_sample_name_dict

    def _generate_and_write_mothur_batch_file_for_dotfile_creation(self):
        self._check_if_fastqs_are_gz_compressed()
        mothur_batch_file_as_list = self._generate_mothur_batch_file_for_dotfile_creation_as_list()
        self.path_to_latest_mothur_batch_file = f'{self.temp_working_directory}/mothur_batch_file_makeFile'
        general.write_list_to_destination(self.path_to_latest_mothur_batch_file, mothur_batch_file_as_list)

    def _generate_mothur_batch_file_for_dotfile_creation_as_list(self):
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

    def _check_if_fastqs_are_gz_compressed(self):
        if self.list_of_fastq_files_in_wkd[0].endswith('fastq.gz'):
            self.fastqs_are_gz_compressed = True
        elif self.list_of_fastq_files_in_wkd[0].endswith('fastq'):
            self.fastqs_are_gz_compressed = False
        else:
            warning_str = f'Unrecognised format of sequecing file: {self.list_of_fastq_files_in_wkd[0]}'
            self._exit_and_del_data_set_sample(warning_str)

    def _check_all_fastqs_in_datasheet_exist(self):
        self.list_of_fastq_files_in_wkd = [_ for _ in os.listdir(self.temp_working_directory) if 'fastq' in _]
        list_of_meta_gz_files = self._get_list_of_fastq_file_names_that_should_be_in_directory()
        self._if_fastq_files_missing_sys_exit(list_of_meta_gz_files)

    def _if_fastq_files_missing_sys_exit(self, list_of_meta_gz_files):
        for fastq in list_of_meta_gz_files:
            if fastq not in self.list_of_fastq_files_in_wkd:
                warning_str = f'{fastq} not found'
                self._exit_and_del_data_set_sample(warning_str)

    def _exit_and_del_data_set_sample(self, warning_str):
        self.dataset_object.delete()
        sys.exit(warning_str)

    def _get_list_of_fastq_file_names_that_should_be_in_directory(sample_meta_info_df):
        list_of_meta_gz_files = []
        list_of_meta_gz_files.extend(sample_meta_info_df['fastq_fwd_file_name'].values.tolist())
        list_of_meta_gz_files.extend(sample_meta_info_df['fastq_rev_file_name'].values.tolist())
        return list_of_meta_gz_files

    def _create_sample_meta_info_dataframe_from_datasheet_path(self):
        if self.datasheet_path.endswith('.xlsx'):
            self.sample_meta_info_df = pd.read_excel(
                io=self.datasheet_path, header=0, index_col=0, usecols='A:N', skiprows=[0])
        elif self.datasheet_path.endswith('.csv'):
            self.sample_meta_info_df = pd.read_csv(
                filepath_or_buffer=self.datasheet_path, header=0, index_col=0, skiprows=[0])
        else:
            sys.exit('Data sheet: {} is in an unrecognised format. '
                     'Please ensure that it is either in .xlsx or .csv format.')

    def _copy_and_decompress_input_files_to_temp_wkd(self):
        if not self.is_single_file_or_paired_input:
            self._copy_fastq_files_from_input_dir_to_temp_wkd()
        else:
            self._extract_single_compressed_file_to_temp_wkd()

    def _extract_single_compressed_file_to_temp_wkd(self):
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

    def _copy_fastq_files_from_input_dir_to_temp_wkd(self):
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

    def _determine_if_single_file_or_paired_input(self):
        for file in os.listdir(self.user_input_path):
            if 'fastq' in file or 'fq' in file:
                # Then there is a fastq.gz or fastq file already uncompressed in this folder
                # In this case we will assume that the seq data is not a single file containing the pairs of files
                # rather the pairs of files themselves.
                return False

    def _create_pre_med_write_out_directory_path(self):
        pre_med_write_out_directory_path = os.path.join(self.output_directory, 'pre_med_seqs')
        os.makedirs(pre_med_write_out_directory_path, exist_ok=True)
        return pre_med_write_out_directory_path

    def _setup_output_directory(self):
        output_directory = os.path.join(self.symportal_root_directory,
                                        'outputs', 'data_set_submissions', f'{self.dataset_object.id}')
        os.makedirs(output_directory, exist_ok=True)
        return output_directory

    def _setup_sequence_dump_file_path(self):
        seq_dump_file_path = os.path.join(
            self.symportal_root_directory, f'/dbBackUp/seq_dumps/seq_dump_{self.data_time_string}')
        os.makedirs(os.path.dirname(seq_dump_file_path), exist_ok=True)
        return seq_dump_file_path

    def _setup_temp_working_directory(self):
        # working directory will be housed in a temp folder within the directory in which the sequencing data
        # is currently housed
        if '.' in self.user_input_path.split('/')[-1]:
            # then this path points to a file rather than a directory and we should pass through the path only
            self.temp_working_directory = os.path.abspath(
                f'{os.path.dirname(self.user_input_path)}/tempData/{self.dataset_object.id}')
        else:
            # then we assume that we are pointing to a directory and we can directly use that to make the wkd
            self.temp_working_directory = os.path.abspath(
                f'{self.user_input_path}/tempData/{self.dataset_object.id}')
        self._create_temp_wkd()
        return self.temp_working_directory

    def _create_temp_wkd(self):
        # if the directory already exists remove it and start from scratch
        if os.path.exists(self.temp_working_directory):
            shutil.rmtree(self.temp_working_directory)
        os.makedirs(self.temp_working_directory)
