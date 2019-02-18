from dbApp.models import (DataSet, ReferenceSequence, DataSetSample, DataSetSampleSequence, CladeCollection)
import sys
import os
import shutil
import subprocess
import glob
import pandas as pd
import general
import json
from collections import Counter
from django import db
from multiprocessing import Queue, Manager, Process
from general import write_list_to_destination, read_defined_file_to_list, create_dict_from_fasta, make_new_blast_db, decode_utf8_binary_to_list, return_list_of_file_paths_in_directory, return_list_of_file_names_in_directory
from datetime import datetime
from distance import BrayCurtisDistPCoACreator, UnifracDistPCoACreator
from plotting import DistScatterPlotterSamples, SeqStackedBarPlotter
from symportal_utils import BlastnAnalysis, MothurAnalysis, NucleotideSequence
from output import SequenceCountTableCreator


class DataLoading:
    # The clades refer to the phylogenetic divisions of the Symbiodiniaceae. Most of them are represented at the genera
    # level. E.g. All members of clade C belong to the genus Cladocopium.
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    def __init__(
            self, data_set_object, user_input_path, datasheet_path, screen_sub_evalue, num_proc,no_fig, no_ord,
            distance_method, debug=False):
        self.symportal_root_directory = os.path.abspath(os.path.dirname(__file__))
        self.output_path_list = []
        self.no_fig = no_fig
        self.no_ord = no_ord
        self.distance_method = distance_method
        self.dataset_object = data_set_object
        self.user_input_path = user_input_path
        self.output_directory = self._setup_output_directory()
        self.temp_working_directory = self._setup_temp_working_directory()
        self.date_time_string = str(datetime.now()).replace(' ', '_').replace(':', '-')
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
            self.output_directory, 'non_sym_and_size_violation_sequences'
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

        self._output_data_set_info()

        self._output_seqs_count_table()

        self._output_seqs_stacked_bar_plot()

        if not self.no_ord:
            print('Calculating between sample pairwise distances')
            if self.distance_method == 'unifrac':
                self._do_unifrac_dist_pcoa()
            elif self.distance_method == 'braycurtis':
                self._do_braycurtis_dist_pcoa()

            # distance plotting
            if not self.no_fig:
                if self.num_of_samples > 1000:
                    print('Too many samples ({}) to generate plots'.format(self.num_of_samples))
                else:
                    for output_path in self.output_path_list:
                        if self.this_is_pcoa_path(output_path):
                            clade_of_output = os.path.dirname(output_path).split('/')[-1]
                            sys.stdout.write(f'\n\nGenerating between sample distance plot clade {clade_of_output}\n')
                            dist_scatter_plotter_samples = DistScatterPlotterSamples(csv_path=output_path, date_time_str=self.date_time_string)
                            dist_scatter_plotter_samples.make_sample_dist_scatter_plot()
                            self.output_path_list.extend(dist_scatter_plotter_samples.output_path_list)

    def _output_data_set_info(self):
        print(f'\n\nData loading complete. DataSet UID: {self.dataset_object.id}')

    @staticmethod
    def this_is_pcoa_path(output_path):
        return 'PCoA_coords' in output_path

    def _do_braycurtis_dist_pcoa(self):
        bray_curtis_dist_pcoa_creator = BrayCurtisDistPCoACreator(
            date_time_string=self.date_time_string, symportal_root_directory=self.symportal_root_directory,
            data_set_string=str(self.dataset_object.id), call_type='submission',
            output_dir=self.output_directory)
        bray_curtis_dist_pcoa_creator.compute_braycurtis_dists_and_pcoa_coords()
        self.output_path_list.extend(bray_curtis_dist_pcoa_creator.output_file_paths)

    def _do_unifrac_dist_pcoa(self):
        unifrac_dict_pcoa_creator = UnifracDistPCoACreator(
            call_type='submission', date_time_string=self.date_time_string, output_dir=self.output_directory,
            data_set_string=str(self.dataset_object.id), method='mothur', num_processors=self.num_proc,
            symportal_root_directory=self.symportal_root_directory)
        unifrac_dict_pcoa_creator.compute_unifrac_dists_and_pcoa_coords()
        self.output_path_list.extend(unifrac_dict_pcoa_creator.output_file_paths)

    def _output_seqs_stacked_bar_plot(self):
        if not self.no_fig:
            if self.num_of_samples > 1000:
                print(f'Too many samples ({self.num_samples}) to generate plots')
            else:
                sys.stdout.write('\nGenerating sequence count table figures\n')

                seq_stacked_bar_plotter = SeqStackedBarPlotter(output_directory=self.output_directory,
                                                               seq_relative_abund_count_table_path=self.seq_abundance_relative_output_path, time_date_str=self.date_time_string)
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
            data_set_uids_to_output_as_comma_sep_string=str(self.dataset_object.id),
            num_proc=self.num_proc, time_date_str=self.date_time_string)
        sequence_count_table_creator.make_output_tables()
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
        symclade_current_path = os.path.abspath(os.path.join(self.symportal_root_directory, 'symbiodiniumDB', 'symClade.fa'))

        symclade_backup_path = os.path.join(back_up_dir + f'symClade_{self.date_time_string}.fa')
        symclade_backup_readme_path = os.path.join(back_up_dir, f'symClade_{self.date_time_string}.readme')
        # then write a copy to it.
        shutil.copy(symclade_current_path, symclade_backup_path)
        # Then write out a very breif readme
        read_me = [
            f'This is a symClade.fa backup created during datasubmission of data_set ID: {self.dataset_object.id}']
        write_list_to_destination(symclade_backup_readme_path, read_me)

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
        self.sample_fastq_pairs = new_stability_file
        # write out the new stability file
        write_list_to_destination(os.path.join(self.temp_working_directory, 'stability.files'), self.sample_fastq_pairs)

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

        self.sample_fastq_pairs = new_stability_file
        # write out the new stability file
        write_list_to_destination(os.path.join(self.temp_working_directory, 'stability.files'), self.sample_fastq_pairs)

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
        self.sample_fastq_pairs = read_defined_file_to_list(os.path.join(self.temp_working_directory, 'stability.files'))

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
        self.path_to_latest_mothur_batch_file = os.path.join(self.temp_working_directory, 'mothur_batch_file_makeFile')
        write_list_to_destination(self.path_to_latest_mothur_batch_file, mothur_batch_file_as_list)

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
        if self.list_of_fastq_files_in_wkd[0].endswith('fastq.gz') or self.list_of_fastq_files_in_wkd[0].endswith('fq.gz'):
            self.fastqs_are_gz_compressed = True
        elif self.list_of_fastq_files_in_wkd[0].endswith('fastq') or self.list_of_fastq_files_in_wkd[0].endswith('fq'):
            self.fastqs_are_gz_compressed = False
        else:
            warning_str = f'Unrecognised format of sequecing file: {self.list_of_fastq_files_in_wkd[0]}'
            self._exit_and_del_data_set_sample(warning_str)

    def _check_all_fastqs_in_datasheet_exist(self):
        self.list_of_fastq_files_in_wkd = return_list_of_file_names_in_directory(self.temp_working_directory)
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

    def _get_list_of_fastq_file_names_that_should_be_in_directory(self):
        list_of_meta_gz_files = []
        list_of_meta_gz_files.extend(self.sample_meta_info_df['fastq_fwd_file_name'].values.tolist())
        list_of_meta_gz_files.extend(self.sample_meta_info_df['fastq_rev_file_name'].values.tolist())
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
        """Need to take into account that there may be other non fastq files in the dir. Also need to take into account
        that they could be .fq files rather than fastq. also need to take into account that it could be fastq.gz and
        fq.gz rather than fastq and fq."""
        list_of_files_in_user_input_dir = return_list_of_file_paths_in_directory(self.user_input_path)
        for file_path_outer in list_of_files_in_user_input_dir:
            if 'fastq' in file_path_outer:
                for file_path_inner in [file_path for file_path in list_of_files_in_user_input_dir if
                                        (file_path.endswith('fastq') or file_path.endswith('fastq.gz'))]:
                    # do the copy for each file if fastq is in the file
                    shutil.copy(file_path_inner, self.temp_working_directory)
            elif 'fq' in file_path_outer:
                for file_path_inner in [file_path for file_path in list_of_files_in_user_input_dir if
                                        (file_path.endswith('fq') or file_path.endswith('fq.gz'))]:
                    # do the copy for each file if fastq is in the file
                    shutil.copy(file_path_inner, self.temp_working_directory)

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
            self.symportal_root_directory,'dbBackUp', 'seq_dumps', f'seq_dump_{self.date_time_string}')
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


class InitialMothurHandler:
    def __init__(self, sample_fastq_pairs, num_proc):
        self.input_queue_containing_pairs_of_fastq_file_paths = Queue()
        self.worker_manager = Manager()
        self.samples_that_caused_errors_in_qc_mp_list = self.worker_manager.list()
        self.num_proc = num_proc
        self._populate_input_queue(sample_fastq_pairs)

    def _populate_input_queue(self, sample_fastq_pairs):
        for fastq_path_pair in sample_fastq_pairs:
            self.input_queue_containing_pairs_of_fastq_file_paths.put(fastq_path_pair)

        for n in range(self.num_proc):
            self.input_queue_containing_pairs_of_fastq_file_paths.put('STOP')

    def execute_worker_initial_mothur(
            self, data_loading_dataset_object, data_loading_temp_working_directory, data_loading_debug
    ):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        sys.stdout.write('\nPerforming QC\n')
        for n in range(self.num_proc):
            p = Process(target=self._worker_initial_mothur,
                        args=(data_loading_dataset_object, data_loading_temp_working_directory, data_loading_debug)
                        )

            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

    def _worker_initial_mothur(
            self, data_loading_dataset_object, data_loading_temp_working_directory, data_loading_debug
    ):
        """
        This worker performs the pre-MED processing that is primarily mothur-based.
        This QC includes making contigs, screening for ambigous calls (0 allowed) and homopolymer (maxhomop=5)
        Discarding singletons and doublets, in silico PCR. It also checks whether sequences are rev compliment.
        This is all done through the use of an InitialMothurWorker class which in turn makes use of the MothurAnalysis
        class that does the heavy lifting of running the mothur commands in sequence.
        """

        for contigpair in iter(self.input_queue_containing_pairs_of_fastq_file_paths.get, 'STOP'):

            initial_morthur_worker = InitialMothurWorker(
                contig_pair=contigpair,
                data_set_object=data_loading_dataset_object,
                temp_wkd=data_loading_temp_working_directory,
                debug=data_loading_debug
            )

            try:
                initial_morthur_worker.execute()
            except RuntimeError as e:
                self.samples_that_caused_errors_in_qc_mp_list.append(e.args[0]['sample_name'])
        return


class InitialMothurWorker:
    def __init__(self, contig_pair, data_set_object, temp_wkd, debug):
        self.sample_name = contig_pair.split('\t')[0].replace('[dS]', '-')
        self.data_set_sample = DataSetSample.objects.get(
                name=self.sample_name, data_submission_from=data_set_object
            )
        self.cwd = os.path.join(temp_wkd, self.sample_name)
        os.makedirs(self.cwd, exist_ok=True)
        self.pre_med_seq_dump_dir = self.cwd.replace('tempData', 'pre_MED_seqs')
        os.makedirs(self.pre_med_seq_dump_dir, exist_ok=True)
        self.mothur_analysis_object = MothurAnalysis(
            pcr_analysis_name='symvar',
            input_dir=self.cwd,
            output_dir=self.cwd,
            fastq_gz_fwd_path=contig_pair.split('\t')[1],
            fastq_gz_rev_path=contig_pair.split('\t')[2],
            stdout_and_sterr_to_pipe=(not debug))
        self.debug = debug

    def execute(self):
        sys.stdout.write(f'{self.sample_name}: QC started\n')

        self._do_make_contigs()

        self._set_absolute_num_seqs_after_make_contigs()

        self._do_screen_seqs()

        self._do_unique_seqs()

        self._do_split_abund()

        self._do_unique_seqs()

        self._do_fwd_and_rev_pcr()

        self._do_unique_seqs()

        self._set_unique_num_seqs_after_initial_qc()

        self._set_absolute_num_seqs_after_inital_qc()

        self._save_changes_to_data_set_sample()

        sys.stdout.write(f'{self.sample_name}: Initial mothur complete\n')

        self._write_out_final_name_and_fasta_for_tax_screening()

    def _write_out_final_name_and_fasta_for_tax_screening(self):
        name_file_as_list = read_defined_file_to_list(self.mothur_analysis_object.name_file_path)
        taxonomic_screening_name_file_path = os.path.join(self.cwd, 'name_file_for_tax_screening.names')
        write_list_to_destination(taxonomic_screening_name_file_path, name_file_as_list)
        fasta_file_as_list = read_defined_file_to_list(self.mothur_analysis_object.fasta_path)
        taxonomic_screening_fasta_file_path = os.path.join(self.cwd, 'fasta_file_for_tax_screening.fasta')
        write_list_to_destination(taxonomic_screening_fasta_file_path, fasta_file_as_list)

    def _save_changes_to_data_set_sample(self):
        self.data_set_sample.save()

    def _do_fwd_and_rev_pcr(self):
        self.mothur_analysis_object.execute_pcr(do_reverse_pcr_as_well=True)
        self.check_for_no_seqs_after_pcr_and_raise_runtime_error()

    def _do_split_abund(self):
        self.mothur_analysis_object.execute_split_abund(abund_cutoff=2)
        self.check_for_error_and_raise_runtime_error()

    def _do_unique_seqs(self):
        self.mothur_analysis_object.execute_unique_seqs()
        self.check_for_error_and_raise_runtime_error()

    def _do_screen_seqs(self):
        self.mothur_analysis_object.execute_screen_seqs(argument_dictionary={'maxambig': 0, 'maxhomop': 5})
        self.check_for_error_and_raise_runtime_error()

    def _do_make_contigs(self):
        self.mothur_analysis_object.execute_make_contigs()
        self.check_for_error_and_raise_runtime_error()

    def _set_absolute_num_seqs_after_make_contigs(self):
        number_of_contig_seqs_absolute = len(
            read_defined_file_to_list(self.mothur_analysis_object.latest_summary_path)
        ) - 1
        self.data_set_sample.num_contigs = number_of_contig_seqs_absolute
        sys.stdout.write(
            f'{self.sample_name}: data_set_sample_instance_in_q.num_contigs = {number_of_contig_seqs_absolute}\n')

    def _set_unique_num_seqs_after_initial_qc(self):
        number_of_contig_seqs_unique = len(
            read_defined_file_to_list(self.mothur_analysis_object.latest_summary_path)) - 1
        self.data_set_sample.post_qc_unique_num_seqs = number_of_contig_seqs_unique
        sys.stdout.write(
            f'{self.sample_name}: '
            f'data_set_sample_instance_in_q.post_qc_unique_num_seqs = {number_of_contig_seqs_unique}\n')

    def _set_absolute_num_seqs_after_inital_qc(self):
        last_summary = read_defined_file_to_list(self.mothur_analysis_object.latest_summary_path)
        absolute_count = 0
        for line in last_summary[1:]:
            absolute_count += int(line.split('\t')[6])
        self.data_set_sample.post_qc_absolute_num_seqs = absolute_count
        sys.stdout.write(
            f'{self.sample_name}: data_set_sample_instance_in_q.post_qc_absolute_num_seqs = {absolute_count}\n')

    def check_for_no_seqs_after_pcr_and_raise_runtime_error(self):
        if len(self.mothur_analysis_object.sequence_collection) == 0:
            self.log_qc_error_and_continue(errorreason='No seqs left after PCR')
            raise RuntimeError({'sample_name':self.sample_name})

    def check_for_error_and_raise_runtime_error(self):
        for stdout_line in decode_utf8_binary_to_list(
                self.mothur_analysis_object.latest_completed_process_command.stdout
        ):
            if '[WARNING]: Blank fasta name, ignoring read.' in stdout_line:
                self.log_qc_error_and_continue(errorreason='Blank fasta name')
                raise RuntimeError({'sample_name':self.sample_name})
            if 'ERROR' in stdout_line:
                self.log_qc_error_and_continue(errorreason='error in inital QC')
                raise RuntimeError({'sample_name':self.sample_name})

    def log_qc_error_and_continue(self, errorreason):
        print('Error in processing sample: {}'.format(self.sample_name))
        self.data_set_sample.unique_num_sym_seqs = 0
        self.data_set_sample.absolute_num_sym_seqs = 0
        self.data_set_sample.initial_processing_complete = True
        self.data_set_sample.error_in_processing = True
        self.data_set_sample.error_reason = errorreason
        self._save_changes_to_data_set_sample()


class PotentialSymTaxScreeningHandler:
    """ The purpose of this handler and the executed work is only to get a collection of sequences that will need
    screening against the NCBI database. We also rely on this method to do the blast of our each samples sequences
    against our symclade.fa reference sequences database. We read in the blast.out file in the later functions.
    """
    def __init__(self, samples_that_caused_errors_in_qc_list, checked_samples_list, list_of_samples_names, num_proc):
        self.input_queue = Queue()
        self.manager = Manager()
        self.sub_evalue_sequence_to_num_sampes_found_in_mp_dict = self.manager.dict()
        self.sub_evalue_nucleotide_sequence_to_clade_mp_dict = self.manager.dict()
        self.error_samples_mp_list = self.manager.list(samples_that_caused_errors_in_qc_list)
        self.checked_samples_mp_list = self.manager.list(checked_samples_list)
        self.list_of_sample_names = list_of_samples_names
        self.num_proc = num_proc
        self._load_input_queue()

    def _load_input_queue(self):
        # load up the input q
        for sample_name in self.list_of_sample_names:
            self.input_queue.put(sample_name)
        # load in the STOPs
        for n in range(self.num_proc):
            self.input_queue.put('STOP')

    def execute_potential_sym_tax_screening(
            self, data_loading_temp_working_directory, data_loading_path_to_symclade_db, data_loading_debug):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        sys.stdout.write('\nPerforming QC\n')
        for n in range(self.num_proc):
            new_process = Process(
                target=self._potential_sym_tax_screening_worker,
                args=(data_loading_temp_working_directory, data_loading_path_to_symclade_db, data_loading_debug))

            all_processes.append(new_process)
            new_process.start()
        for process in all_processes:
            process.join()

    def _potential_sym_tax_screening_worker(
            self, data_loading_temp_working_directory, data_loading_path_to_symclade_db, data_loading_debug):
        """
        input_q: The multiprocessing queue that holds a list of the sample names
        e_val_collection_dict: This is a managed dictionary where key is a nucleotide sequence that has:
        1 - provided a match in the blast analysis
        2 - is of suitable size
        3 - but has an evalue match below the cuttof
        the value of the dict is an int that represents how many samples this nucleotide sequence was found in
        err_smpl_list: This is a managed list containing sample names of the samples that had errors during the
        initial mothur qc and therefore don't require taxonomic screening performed on them
        checked_samples: This is a list of sample names for samples that were found to contain only
        Symbiodinium sequences or have already had all potential Symbiodinium sequences screened and so don't
        require any further taxonomic screening
        Whilst this doesn't return anything a number of objects are picked out in each of the local
        working directories for use in the workers that follow this one.
        """
        for sample_name in iter(self.input_queue.get, 'STOP'):

            # If the sample gave an error during the inital mothur then we don't consider it here.
            if sample_name in self.error_samples_mp_list:
                continue

            # A sample will be in this list if we have already performed this worker on it and none of its sequences
            # gave matches to the symClade database at below the evalue threshold
            if sample_name in self.checked_samples_mp_list:
                continue

            taxonomic_screening_worker = PotentialSymTaxScreeningWorker(
                sample_name=sample_name, wkd=data_loading_temp_working_directory,
                path_to_symclade_db=data_loading_path_to_symclade_db, debug=data_loading_debug,
                checked_samples_mp_list=self.checked_samples_mp_list,
                e_val_collection_mp_dict=self.sub_evalue_sequence_to_num_sampes_found_in_mp_dict,
                sub_evalue_nucleotide_sequence_to_clade_mp_dict=self.sub_evalue_nucleotide_sequence_to_clade_mp_dict)

            taxonomic_screening_worker.execute()


class PotentialSymTaxScreeningWorker:
    def __init__(
            self, sample_name, wkd, path_to_symclade_db, debug, e_val_collection_mp_dict,
            checked_samples_mp_list, sub_evalue_nucleotide_sequence_to_clade_mp_dict):
        self.sample_name = sample_name
        self.cwd = os.path.join(wkd, self.sample_name)
        self.fasta_file_path = os.path.join(self.cwd, 'fasta_file_for_tax_screening.fasta')
        self.fasta_dict = create_dict_from_fasta(fasta_path=self.fasta_file_path)
        self.name_file_path = os.path.join(self.cwd, 'name_file_for_tax_screening.names')
        self.name_dict = {a.split('\t')[0]: a for a in read_defined_file_to_list(self.name_file_path)}
        self.path_to_symclade_db = path_to_symclade_db
        self.debug = debug
        # This is a managed dictionary where key is a nucleotide sequence that has:
        # 1 - provided a match in the blast analysis
        # 2 - is of suitable size
        # 3 - but has an evalue match below the cuttof
        # the value of the dict is an int that represents how many samples this nucleotide sequence was found in
        self.e_val_collection_mp_dict = e_val_collection_mp_dict
        # This dictionary will be used outside of the multiprocessing to append the clade of a given sequences
        # that is being added to the symClade reference database
        self.sub_evalue_nucleotide_sequence_to_clade_mp_dict = sub_evalue_nucleotide_sequence_to_clade_mp_dict
        # The potential_non_symbiodinium_sequences_list is used to see if there are any samples, that don't have
        # any potential non symbiodnium sequences. i.e. only definite symbiodinium sequences. These samples are added to
        # the checked list and are not checked in following iterations to speed things up.
        self.potential_non_symbiodinium_sequences_list = []
        self.sequence_name_to_clade_dict = None
        self.blast_output_as_list = None
        self.already_processed_blast_seq_result = []
        # this is a managed list that holds the names of samples from which no sequences were thrown out from
        # it will be used in downstream processes.
        self.checked_samples_mp_list = checked_samples_mp_list

    def execute(self):
        sys.stdout.write(f'{self.sample_name}: verifying seqs are Symbiodinium and determining clade\n')

        blastn_analysis = BlastnAnalysis(
            input_file_path=self.fasta_file_path,
            output_file_path=os.path.join(self.cwd, 'blast.out'), db_path=self.path_to_symclade_db,
            output_format_string="6 qseqid sseqid staxids evalue pident qcovs")

        if self.debug:
            blastn_analysis.execute(pipe_stdout_sterr=False)
        else:
            blastn_analysis.execute(pipe_stdout_sterr=True)

        sys.stdout.write(f'{self.sample_name}: BLAST complete\n')

        self.blast_output_as_list = blastn_analysis.return_blast_output_as_list()

        self._if_debug_warn_if_blast_out_empty_or_low_seqs()

        self.sequence_name_to_clade_dict = {
            blast_out_line.split('\t')[0]: blast_out_line.split('\t')[1][-1] for
            blast_out_line in self.blast_output_as_list}

        self._add_seqs_with_no_blast_match_to_non_sym_list()

        # NB blast results sometimes return several matches for the same seq.
        # as such we will use the already_processed_blast_seq_resulst to make sure that we only
        # process each sequence once.
        self._identify_and_allocate_non_sym_and_sub_e_seqs()

        if not self.potential_non_symbiodinium_sequences_list:
            self.checked_samples_mp_list.append(self.sample_name)

    def _identify_and_allocate_non_sym_and_sub_e_seqs(self):
        for line in self.blast_output_as_list:
            name_of_current_sequence = line.split('\t')[0]
            if name_of_current_sequence in self.already_processed_blast_seq_result:
                continue
            self.already_processed_blast_seq_result.append(name_of_current_sequence)
            identity = float(line.split('\t')[4])
            coverage = float(line.split('\t')[5])

            # noinspection PyPep8,PyBroadException
            # here we are looking for sequences to add to the non_symbiodinium_sequence_list
            # if a sequence fails at any of our if statements it will be added to the non_symbiodinium_sequence_list
            try:
                evalue_power = int(line.split('\t')[3].split('-')[1])
                # With the smallest sequences i.e. 185bp it is impossible to get above the 100 threshold
                # even if there is an exact match. As such, we sould also look at the match identity and coverage
                if evalue_power < 100:  # evalue cut off, collect sequences that don't make the cut
                    self._if_ident_cov_size_good_add_seq_to_non_sym_list_and_eval_dict(
                        coverage, identity, name_of_current_sequence
                    )
            except:
                # here we weren't able to extract the evalue_power for some reason.
                self._if_ident_cov_size_good_add_seq_to_non_sym_list_and_eval_dict(
                    coverage, identity, name_of_current_sequence
                )

    def _if_ident_cov_size_good_add_seq_to_non_sym_list_and_eval_dict(
            self, coverage, identity, name_of_current_sequence):
        """This method will add nucleotide sequences that gave blast matches but that were below the evalue
        or indentity and coverage thresholds to the evalue collection dict and
        record how many samples that sequence was found in.
        It will also take into account the size thresholds that would normally
        happen later in the code during further mothur qc.
        Finally, it will also populate a nucleotide sequence to clade dictionary that will be used outside of
        the MPing to append the clade of a given sequences that is being added to the symClade reference database.
        The potential_non_symbiodinium_sequences_list is used to see if there are any samples, that don't have
        any potential non symbiodnium sequences. i.e. only definite symbiodinium sequences. These samples are added to
        the checked list and are not checked in following iterations to speed things up.
        """

        if identity < 80 or coverage < 95:
            # incorporate the size cutoff here that would normally happen in the further mothur qc later in the code
            self.potential_non_symbiodinium_sequences_list.append(name_of_current_sequence)
            if 184 < len(self.fasta_dict[name_of_current_sequence]) < 310:
                if self.fasta_dict[name_of_current_sequence] in self.e_val_collection_mp_dict.keys():
                    self.e_val_collection_mp_dict[self.fasta_dict[name_of_current_sequence]] += 1
                else:
                    self.e_val_collection_mp_dict[self.fasta_dict[name_of_current_sequence]] = 1
                    self.sub_evalue_nucleotide_sequence_to_clade_mp_dict[
                        self.fasta_dict[name_of_current_sequence]
                    ] = self.sequence_name_to_clade_dict[name_of_current_sequence]

    def _add_seqs_with_no_blast_match_to_non_sym_list(self):
        sequences_with_no_blast_match_as_set = set(self.fasta_dict.keys()) - \
                                               set(self.sequence_name_to_clade_dict.keys())
        self.potential_non_symbiodinium_sequences_list.extend(list(sequences_with_no_blast_match_as_set))
        sys.stdout.write(
            f'{self.sample_name}: {len(sequences_with_no_blast_match_as_set)} sequences thrown out '
            f'initially due to being too divergent from reference sequences\n')

    def _if_debug_warn_if_blast_out_empty_or_low_seqs(self):
        if self.debug:
            if not self.blast_output_as_list:
                print(f'WARNING blast output file is empty for {self.sample_name}')
            else:
                if len(self.blast_output_as_list) < 10:
                    print(
                        f'WARNING blast output file for {self.sample_name} '
                        f'is only {len(self.blast_output_as_list)} lines long')


class SymNonSymTaxScreeningHandler:
    def __init__(
            self, data_loading_samples_that_caused_errors_in_qc_mp_list, data_loading_list_of_samples_names,
            data_loading_num_proc):
        self.sample_name_mp_input_queue = Queue()
        self.sym_non_sym_mp_manager = Manager()
        self.samples_that_caused_errors_in_qc_mp_list = self.sym_non_sym_mp_manager.list(
            data_loading_samples_that_caused_errors_in_qc_mp_list
        )
        self.non_symbiodinium_sequences_list = self.sym_non_sym_mp_manager.list()
        self.num_proc = data_loading_num_proc
        self._populate_input_queue(data_loading_list_of_samples_names)

    def _populate_input_queue(self, data_loading_list_of_samples_names):
        for sample_name in data_loading_list_of_samples_names:
            self.sample_name_mp_input_queue.put(sample_name)
        for n in range(self.num_proc):
            self.sample_name_mp_input_queue.put('STOP')

    def execute_sym_non_sym_tax_screening(
            self, data_loading_temp_working_directory, data_loading_dataset_object,
            non_symb_and_size_violation_base_dir_path, data_loading_pre_med_sequence_output_directory_path,
            data_loading_debug):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        sys.stdout.write('\nPerforming QC\n')

        for n in range(self.num_proc):
            p = Process(target=self._sym_non_sym_tax_screening_worker, args=(
                data_loading_temp_working_directory, data_loading_dataset_object,
                non_symb_and_size_violation_base_dir_path, data_loading_pre_med_sequence_output_directory_path,
                data_loading_debug))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

    def _sym_non_sym_tax_screening_worker(self, data_loading_temp_working_directory, data_loading_dataset_object,
                                          data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path,
                                          data_loading_pre_med_sequence_output_directory_path, data_loading_debug):

        for sample_name in iter(self.sample_name_mp_input_queue.get, 'STOP'):
            if sample_name in self.samples_that_caused_errors_in_qc_mp_list:
                continue

            sym_non_sym_tax_screening_worker_object = SymNonSymTaxScreeningWorker(
                data_loading_temp_working_directory=data_loading_temp_working_directory,
                data_loading_dataset_object=data_loading_dataset_object, sample_name=sample_name,
                data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path=
                data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path,
                data_loading_pre_med_sequence_output_directory_path=data_loading_pre_med_sequence_output_directory_path,
                data_loading_debug=data_loading_debug
            )

            try:
                sym_non_sym_tax_screening_worker_object.identify_sym_non_sym_seqs()
            except RuntimeError as e:
                self.samples_that_caused_errors_in_qc_mp_list.append(e.args[0]['sample_name'])


class SymNonSymTaxScreeningWorker:
    def __init__(
            self, data_loading_temp_working_directory, sample_name, data_loading_dataset_object,
            data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path,
            data_loading_pre_med_sequence_output_directory_path, data_loading_debug
    ):
        self._init_core_class_attributes(
            data_loading_dataset_object, data_loading_temp_working_directory, sample_name, data_loading_debug)
        self._init_fasta_name_blast_and_clade_dict_attributes()
        self._init_non_sym_and_size_violation_output_paths(
            data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path)
        self._init_sym_non_size_violation_output_paths(data_loading_pre_med_sequence_output_directory_path)
        self._init_sets_for_categorizing_sequences()
        self._init_qc_meta_info_counters()

    def _init_qc_meta_info_counters(self):
        self.absolute_number_of_non_sym_sequences = 0
        self.absolute_number_of_sym_no_size_violation_sequences = 0
        self.absolute_number_of_sym_size_violation_sequences = 0

    def _init_sets_for_categorizing_sequences(self):
        self.non_symbiodinium_sequence_name_set_for_sample = set()
        self.sym_size_violation_sequence_name_set_for_sample = set()
        self.sym_no_size_violation_sequence_name_set_for_sample = set()

    def _init_sym_non_size_violation_output_paths(self, data_loading_pre_med_sequence_output_directory_path):
        pre_med_sequence_output_directory_for_sample_path = os.path.join(
            data_loading_pre_med_sequence_output_directory_path, self.sample_name)
        os.makedirs(pre_med_sequence_output_directory_for_sample_path, exist_ok=True)
        self.pre_med_fasta_path = os.path.join(
            pre_med_sequence_output_directory_for_sample_path, f'pre_med_seqs_{self.sample_name}.fasta'
        )
        self.pre_med_names_file_path = os.path.join(
            pre_med_sequence_output_directory_for_sample_path, f'pre_med_seqs_{self.sample_name}.names'
        )

    def _init_non_sym_and_size_violation_output_paths(
            self, data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path):
        non_symbiodiniaceae_and_size_violation_directory_for_sample_path = os.path.join(
            data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path, self.sample_name)
        os.makedirs(non_symbiodiniaceae_and_size_violation_directory_for_sample_path, exist_ok=True)
        self.non_symbiodiniaceae_seqs_fasta_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.sample_name}_non_symbiodiniaceae_sequences.fasta')
        self.non_symbiodiniaceae_seqs_names_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.sample_name}_non_symbiodiniaceae_sequences.names')
        self.symbiodiniaceae_size_violation_seqs_fasta_output_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.sample_name}_symbiodiniaceae_size_violation_sequences.fasta'
        )
        self.symbiodiniaceae_size_violation_seqs_names_output_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.sample_name}_symbiodiniaceae_size_violation_sequences.names')

    def _init_fasta_name_blast_and_clade_dict_attributes(self):
        fasta_file_path = os.path.join(self.cwd, 'fasta_file_for_tax_screening.fasta')
        self.fasta_dict = create_dict_from_fasta(fasta_path=fasta_file_path)
        name_file_path = os.path.join(self.cwd, 'name_file_for_tax_screening.names')
        self.name_dict = {a.split('\t')[0]: a for a in read_defined_file_to_list(name_file_path)}
        blast_output_path = os.path.join(self.cwd, 'blast.out')
        self.blast_dict = {blast_line.split('\t')[0]: blast_line for blast_line in
                           read_defined_file_to_list(blast_output_path)}
        self.sequence_name_to_clade_dict = {
            blast_out_line.split('\t')[0]: blast_out_line.split('\t')[1][-1] for
            blast_out_line in self.blast_dict.values()
        }

    def _init_core_class_attributes(
            self, data_loading_dataset_object, data_loading_temp_working_directory, sample_name, data_loading_debug):
        self.sample_name = sample_name
        self.cwd = os.path.join(data_loading_temp_working_directory, sample_name)
        self.dataset_object = DataSetSample.objects.get(
            name=sample_name, data_submission_from=data_loading_dataset_object
        )
        self.debug = data_loading_debug

    def identify_sym_non_sym_seqs(self):
        """This method completes the pre-med quality control.
        Having gone through the initial mothur qc, and having screened for potential symbiodiniaceae sequences
        (if running on the remote system), this method now identifies
        1 - the non symbiodiniaceae sequences and writes them out to the output dir
        2 - identifies the symbiodiniaceae sequences that violate our size range thresholds and writes them out
        also to the output dir
        3 - and finally the symbiodiniaceae sequences that do not violate our size range thresholds (sequences that
        will be carried through into med decomposition). These are also written out both as clade seperated
        (one redundant fasta for each clade) in the temp working directory (for med processing),
        and as a non clade separated .name
        .file set (pre-med-seqs) in the output directory.
        This method also populates all of the dataset qc metadata attributes accordingly.
        """
        self._identify_and_write_non_sym_seqs_in_sample()

        self._get_size_violation_and_non_size_violations_seq_sets()

        self._write_out_size_violation_seqs()

        if self.sym_no_size_violation_sequence_name_set_for_sample:
            self._write_out_no_size_violation_seqs()
        else:
            self._log_dataset_attr_and_raise_runtime_error()

        self._associate_qc_meta_info_to_data_sheet()

    def _write_out_no_size_violation_seqs(self):
        self._write_out_no_size_violation_seqs_to_pre_med_dirs()
        clades_of_non_violation_seqs = self._get_set_of_clades_represented_by_no_size_violation_seqs()
        self._write_out_no_size_violation_seqs_redundant_fasta_clade_separated(clades_of_non_violation_seqs)

    def _write_out_no_size_violation_seqs_redundant_fasta_clade_separated(self, clades_of_non_violation_seqs):
        for clade_of_sequences_to_write_out in clades_of_non_violation_seqs:
            sequence_names_of_clade = self._get_sequence_names_of_clade_for_no_size_violation_sequences(
                clade_of_sequences_to_write_out
            )

            sample_clade_fasta_path = os.path.join(
                self.cwd,
                clade_of_sequences_to_write_out,
                f'seqs_for_med_{self.sample_name}_clade_{clade_of_sequences_to_write_out}.redundant.fasta'
            )
            os.makedirs(os.path.dirname(sample_clade_fasta_path), exist_ok=True)
            with open(sample_clade_fasta_path, 'w') as f:
                sequence_counter = 0
                for sequence_name in sequence_names_of_clade:
                    # NB that MED will use the last '_' character as the separator for infering what the sample
                    # name is. As such we either need to make sure that '_' are not found after the '_' that
                    # separates the sample name from the rest of the sequence name, or we need to use another character
                    # for doing the sample name inference. This can be provided to med using the -t argument.
                    for i in range(len(self.name_dict[sequence_name].split('\t')[1].split(','))):
                        f.write(f'>{self.sample_name}_{sequence_counter}\n')
                        f.write(f'{self.fasta_dict[sequence_name]}\n')
                        sequence_counter += 1

            self._if_debug_check_length_of_deuniqued_fasta(sample_clade_fasta_path)

    def _if_debug_check_length_of_deuniqued_fasta(self, sample_clade_fasta_path):
        if self.debug:
            deuniqued_fasta = read_defined_file_to_list(sample_clade_fasta_path)
            if deuniqued_fasta:
                if len(deuniqued_fasta) < 100:
                    print(f'{self.sample_name}: WARNING the dequniqed fasta is < {len(deuniqued_fasta)} lines')
            else:
                print(f'{self.sample_name}: ERROR deuniqued fasta is empty')

    def _get_sequence_names_of_clade_for_no_size_violation_sequences(self, clade_of_sequences_to_write_out):
        sequence_names_of_clade = [
            sequence_name for sequence_name in
            self.sym_no_size_violation_sequence_name_set_for_sample if
            self.sequence_name_to_clade_dict[sequence_name] == clade_of_sequences_to_write_out
        ]
        return sequence_names_of_clade

    def _get_set_of_clades_represented_by_no_size_violation_seqs(self):
        clades_of_non_violation_seqs = set(
            [
                clade_value for sequence_name, clade_value in self.sequence_name_to_clade_dict.items()
                if sequence_name in self.sym_no_size_violation_sequence_name_set_for_sample
            ]
        )
        return clades_of_non_violation_seqs

    def _write_out_no_size_violation_seqs_to_pre_med_dirs(self):
        self._write_out_no_size_violation_fasta_to_pre_med_dir()
        self._write_out_no_size_violation_names_file_to_pre_med_dir()

    def _write_out_no_size_violation_names_file_to_pre_med_dir(self):
        with open(self.pre_med_names_file_path, 'w') as f:
            for sequence_name in list(self.sym_no_size_violation_sequence_name_set_for_sample):
                f.write(f'{self.name_dict[sequence_name]}\n')
                self.absolute_number_of_sym_no_size_violation_sequences += len(
                    self.name_dict[sequence_name].split('\t')[1].split(','))

    def _write_out_no_size_violation_fasta_to_pre_med_dir(self):
        with open(self.pre_med_fasta_path, 'w') as f:
            for sequence_name in list(self.sym_no_size_violation_sequence_name_set_for_sample):
                f.write(f'>{sequence_name}\n')
                f.write(f'{self.fasta_dict[sequence_name]}\n')

    def _write_out_size_violation_seqs(self):
        if self.sym_size_violation_sequence_name_set_for_sample:
            self._write_out_size_violation_fasta()

            self._write_out_size_violation_names_file()

    def _write_out_size_violation_names_file(self):
        with open(self.symbiodiniaceae_size_violation_seqs_names_output_path, 'w') as f:
            for sequence_name in list(self.sym_size_violation_sequence_name_set_for_sample):
                f.write(f'{self.name_dict[sequence_name]}\n')
                self.absolute_number_of_sym_size_violation_sequences += len(
                    self.name_dict[sequence_name].split('\t')[1].split(','))

    def _write_out_size_violation_fasta(self):
        with open(self.symbiodiniaceae_size_violation_seqs_fasta_output_path, 'w') as f:
            for sequence_name in list(self.sym_size_violation_sequence_name_set_for_sample):
                f.write(f'>{sequence_name}\n')
                f.write(f'{self.fasta_dict[sequence_name]}\n')

    def _get_size_violation_and_non_size_violations_seq_sets(self):
        for query_sequence_name, blast_line in self.blast_dict.items():
            if query_sequence_name not in self.non_symbiodinium_sequence_name_set_for_sample:
                if 184 < len(self.fasta_dict[query_sequence_name]) < 310:
                    self.sym_no_size_violation_sequence_name_set_for_sample.add(query_sequence_name)
                else:
                    self.sym_size_violation_sequence_name_set_for_sample.add(query_sequence_name)

    def _associate_qc_meta_info_to_data_sheet(self):
        self._associate_non_sym_seq_attributes_to_dataset()
        self._associate_sym_seq_no_size_violation_attributes_to_dataset()
        self._associate_sym_seq_size_violation_attributes_to_dataset()
        self.dataset_object.initial_processing_complete = True
        self.dataset_object.save()
        print(f'{self.sample_name}: pre-med QC complete')

    def _associate_sym_seq_size_violation_attributes_to_dataset(self):
        self.dataset_object.size_violation_absolute = (
                self.dataset_object.post_qc_absolute_num_seqs -
                self.dataset_object.absolute_num_sym_seqs -
                self.dataset_object.non_sym_absolute_num_seqs
        )
        print(f'{self.sample_name}: size_violation_absolute = {self.dataset_object.size_violation_absolute}')
        self.dataset_object.size_violation_unique = (
                self.dataset_object.post_qc_unique_num_seqs -
                self.dataset_object.unique_num_sym_seqs -
                self.dataset_object.non_sym_unique_num_seqs
        )
        print(f'{self.sample_name}: size_violation_unique = {self.dataset_object.size_violation_unique}')

    def _associate_sym_seq_no_size_violation_attributes_to_dataset(self):
        self.dataset_object.unique_num_sym_seqs = len(self.sym_no_size_violation_sequence_name_set_for_sample)
        self.dataset_object.absolute_num_sym_seqs = self.absolute_number_of_sym_no_size_violation_sequences
        print(f'{self.sample_name}: unique_num_sym_seqs = {self.dataset_object.unique_num_sym_seqs}')
        print(f'{self.sample_name}: absolute_num_sym_seqs = {self.absolute_number_of_sym_no_size_violation_sequences}')

    def _associate_non_sym_seq_attributes_to_dataset(self):
        self.dataset_object.non_sym_unique_num_seqs = len(self.non_symbiodinium_sequence_name_set_for_sample)
        self.dataset_object.non_sym_absolute_num_seqs = self.absolute_number_of_non_sym_sequences
        print(f'{self.sample_name}: non_sym_unique_num_seqs = {self.dataset_object.non_sym_unique_num_seqs}')
        print(f'{self.sample_name}: non_sym_absolute_num_seqs = {self.absolute_number_of_non_sym_sequences}')

    def _log_dataset_attr_and_raise_runtime_error(self):
        # if there are no symbiodiniaceae sequenes then log error and associate meta info
        print(f'{self.sample_name}: QC error.\n No symbiodiniaceae sequences left in sample after pre-med QC.')
        self.dataset_object.non_sym_unique_num_seqs = len(self.non_symbiodinium_sequence_name_set_for_sample)
        self.dataset_object.non_sym_absolute_num_seqs = self.absolute_number_of_non_sym_sequences
        self.dataset_object.size_violation_absolute = self.absolute_number_of_sym_size_violation_sequences
        self.dataset_object.size_violation_unique = len(self.sym_size_violation_sequence_name_set_for_sample)
        self.dataset_object.unique_num_sym_seqs = 0
        self.dataset_object.absolute_num_sym_seqs = 0
        self.dataset_object.initial_processing_complete = True
        self.dataset_object.error_in_processing = True
        self.dataset_object.error_reason = 'No symbiodiniaceae sequences left in sample after pre-med QC'
        self.dataset_object.save()
        raise RuntimeError({'sample_name':self.sample_name})

    def _identify_and_write_non_sym_seqs_in_sample(self):
        self._add_seqs_with_no_blast_match_to_non_sym_list()
        self._identify_non_sym_seqs_from_below_match_threshold()
        self._write_out_non_sym_fasta_and_names_files_for_sample()

    def _write_out_non_sym_fasta_and_names_files_for_sample(self):
        if self.non_symbiodinium_sequence_name_set_for_sample:
            self._write_out_non_sym_fasta_for_sample()
            self._write_out_non_sym_names_file_for_sample()

    def _write_out_non_sym_names_file_for_sample(self):
        with open(self.non_symbiodiniaceae_seqs_names_path, 'w') as f:
            for sequence_name in list(self.non_symbiodinium_sequence_name_set_for_sample):
                f.write(f'{self.name_dict[sequence_name]}\n')
                self.absolute_number_of_non_sym_sequences += len(
                    self.name_dict[sequence_name].split('\t')[1].split(','))

    def _write_out_non_sym_fasta_for_sample(self):
        with open(self.non_symbiodiniaceae_seqs_fasta_path, 'w') as f:
            for sequence_name in list(self.non_symbiodinium_sequence_name_set_for_sample):
                f.write(f'>{sequence_name}\n')
                f.write(f'{self.fasta_dict[sequence_name]}\n')

    def _add_seqs_with_no_blast_match_to_non_sym_list(self):
        sequences_with_no_blast_match_as_set = set(self.fasta_dict.keys()) - \
                                               set(self.sequence_name_to_clade_dict.keys())
        self.non_symbiodinium_sequence_name_set_for_sample.update(list(sequences_with_no_blast_match_as_set))
        sys.stdout.write(
            f'{self.sample_name}: {len(sequences_with_no_blast_match_as_set)} sequences thrown out '
            f'initially due to being too divergent from reference sequences\n')

    def _identify_non_sym_seqs_from_below_match_threshold(self):
        """ This method is where the non_symbiodinium sequences for the sample are identified. If they fall below
        the given identity, coverage and evalues for their match to the symClade db then they are conisdered
        non-symbiodinum. This does not take into size at all. We will screen and report size seperately.
        """
        for blast_sequence_name, blast_line in self.blast_dict.items():
            # noinspection PyPep8,PyBroadException
            try:
                evalue_power = int(blast_line.split('\t')[3].split('-')[1])
                # With the smallest sequences i.e. 185bp it is impossible to get above the 100 threshold
                # even if there is an exact match. As such, we sould also look at the match identity and coverage
                if evalue_power < 100:  # evalue cut off, collect sequences that don't make the cut
                    self._if_low_identity_and_coverage_add_to_non_symb_set(
                        name_of_current_sequence=blast_sequence_name, blast_line=blast_line
                    )
            except:
                # here we weren't able to extract the evalue_power for some reason.
                self._if_low_identity_and_coverage_add_to_non_symb_set(
                    name_of_current_sequence=blast_sequence_name, blast_line=blast_line
                )

    def _if_low_identity_and_coverage_add_to_non_symb_set(self, name_of_current_sequence, blast_line):
        percentage_identity = float(blast_line.split('\t')[4])
        percentage_coverage = float(blast_line.split('\t')[5])
        if percentage_identity < 80 or percentage_coverage < 95:
            # incorporate the size cutoff here that would normally happen in the further mothur qc later in the code
            self.non_symbiodinium_sequence_name_set_for_sample.add(name_of_current_sequence)


class PerformMEDHandler:
    def __init__(self, data_loading_list_of_samples_names, data_loading_temp_working_directory, data_loading_num_proc):
        # need to get list of the directories in which to perform the MED
        # we want to get a list of the .
        self.list_of_samples = data_loading_list_of_samples_names
        self.temp_working_directory = data_loading_temp_working_directory
        self.num_proc = data_loading_num_proc
        self.list_of_redundant_fasta_paths = []
        self._populate_list_of_redundant_fasta_paths()
        self.input_queue_of_redundant_fasta_paths = Queue()
        self._populate_input_queue_of_redundant_fasta_paths()
        self.list_of_med_result_dirs = [
            os.path.join(os.path.dirname(path_to_redundant_fasta), 'MEDOUT') for
            path_to_redundant_fasta in self.list_of_redundant_fasta_paths]

    def execute_perform_med_worker(
            self, data_loading_debug, data_loading_path_to_med_padding_executable,
            data_loading_path_to_med_decompoase_executable):
        all_processes = []

        for n in range(self.num_proc):
            p = Process(target=self._perform_med_worker, args=(
                data_loading_debug, data_loading_path_to_med_padding_executable,
                data_loading_path_to_med_decompoase_executable))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

    def _populate_list_of_redundant_fasta_paths(self):
        for dirpath, dirnames, files in os.walk(self.temp_working_directory):
            for file_name in files:
                if file_name.endswith('redundant.fasta'):
                    self.list_of_redundant_fasta_paths.append(os.path.join(dirpath, file_name))

    def _populate_input_queue_of_redundant_fasta_paths(self):
        for redundant_fasta_path in self.list_of_redundant_fasta_paths:
            self.input_queue_of_redundant_fasta_paths.put(redundant_fasta_path)
        for n in range(self.num_proc):
            self.input_queue_of_redundant_fasta_paths.put('STOP')

    def _perform_med_worker(
            self, data_loading_debug, data_loading_path_to_med_padding_executable,
            data_loading_path_to_med_decompose_executable):
        for redundant_fata_path in iter(self.input_queue_of_redundant_fasta_paths.get, 'STOP'):

            perform_med_worker_instance = PerformMEDWorker(
                redundant_fasta_path=redundant_fata_path, data_loading_debug=data_loading_debug, data_loading_path_to_med_padding_executable=data_loading_path_to_med_padding_executable,
                data_loading_path_to_med_decompose_executable=data_loading_path_to_med_decompose_executable)

            perform_med_worker_instance.do_decomposition()


class PerformMEDWorker:
    def __init__(
            self, redundant_fasta_path, data_loading_path_to_med_padding_executable, data_loading_debug,
            data_loading_path_to_med_decompose_executable):
        self.redundant_fasta_path_unpadded = redundant_fasta_path
        self.redundant_fasta_path_padded = self.redundant_fasta_path_unpadded.replace('.fasta', '.padded.fasta')
        self.cwd = os.path.dirname(self.redundant_fasta_path_unpadded)
        self.sample_name = self.cwd.split('/')[-2]
        self.debug = data_loading_debug
        self.path_to_med_padding_executable = data_loading_path_to_med_padding_executable
        self.path_to_med_decompose_executable = data_loading_path_to_med_decompose_executable
        self.med_output_dir = os.path.join(os.path.dirname(self.redundant_fasta_path_unpadded), 'MEDOUT')
        os.makedirs(self.med_output_dir, exist_ok=True)
        self.med_m_value = self._get_med_m_value()

    def do_decomposition(self):
        # TODO check this business about how MED determines the sample name, i.e. whether the underscores matter.
        sys.stdout.write(f'{self.sample_name}: starting MED analysis\n')
        sys.stdout.write(f'{self.sample_name}: padding sequences\n')
        subprocess.run([
            self.path_to_med_padding_executable,
            '-o', self.redundant_fasta_path_padded,
            self.redundant_fasta_path_unpadded], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sys.stdout.write(f'{self.sample_name}: decomposing\n')
        self._do_the_decomposition()
        sys.stdout.write(f'{self.sample_name}: MED analysis complete\n')

    def _do_the_decomposition(self):
        if not self.debug:
            subprocess.run(
                [self.path_to_med_decompose_executable, '-M', str(self.med_m_value), '--skip-gexf-files',
                 '--skip-gen-figures',
                 '--skip-gen-html', '--skip-check-input', '-o',
                 self.med_output_dir, self.redundant_fasta_path_padded], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elif self.debug:
            subprocess.run(
                [self.path_to_med_decompose_executable, '-M', str(self.med_m_value), '--skip-gexf-files',
                 '--skip-gen-figures',
                 '--skip-gen-html',
                 '--skip-check-input', '-o',
                 self.med_output_dir, self.redundant_fasta_path_padded])

    def _get_med_m_value(self):
        # Define MED M value dynamically.
        # The M value is a cutoff that looks at the abundance of the most abundant unique sequence in a node
        # if the abundance is lower than M then the node is discarded
        # we have been working recently with an M that equivaltes to 0.4% of 0.004. This was
        # calculated when working with a modelling project where I was subsampling to 1000 sequences. In this
        # scenario the M was set to 4.
        # We should also take care that M doesn't go below 4, so we should use a max choice for the M
        num_of_seqs_to_decompose = len(read_defined_file_to_list(self.redundant_fasta_path_unpadded)) / 2
        return max(4, int(0.004 * num_of_seqs_to_decompose))


class DataSetSampleSequenceCreatorWorker:
    """This class will be responsible for handling a set of med outputs. Objects will be things like the directory,
    the count table, number of samples, number of nodes, these sorts of things."""
    def __init__(self, med_output_directory,
                 data_set_sample_creator_handler_ref_seq_sequence_to_ref_seq_id_dict,
                 data_set_sample_creator_handler_ref_seq_uid_to_ref_seq_name_dict, data_loading_dataset_obj):
        self.output_directory = med_output_directory
        self.sample_name = self.output_directory.split('/')[-3]
        self.clade = self.output_directory.split('/')[-2]
        self.nodes_list_of_nucleotide_sequences = []
        self._populate_nodes_list_of_nucleotide_sequences()
        self.num_med_nodes = len(self.nodes_list_of_nucleotide_sequences)
        self.node_sequence_name_to_ref_seq_id = {} # TODO pass in a self of the parent to get rid of these long names
        self.ref_seq_sequence_to_ref_seq_id_dict = data_set_sample_creator_handler_ref_seq_sequence_to_ref_seq_id_dict
        self.ref_seq_uid_to_ref_seq_name_dict = data_set_sample_creator_handler_ref_seq_uid_to_ref_seq_name_dict
        self.node_abundance_df = pd.read_csv(
            os.path.join(self.output_directory, 'MATRIX-COUNT.txt'), delimiter='\t', header=0, index_col=0)
        self.total_num_sequences = sum(self.node_abundance_df.iloc[0])
        self.dataset_sample_object = DataSetSample.objects.get(
            data_submission_from=data_loading_dataset_obj, name=self.sample_name)
        self.clade_collection_object = None

    def _populate_nodes_list_of_nucleotide_sequences(self):
        node_file_path = os.path.join(self.output_directory, 'NODE-REPRESENTATIVES.fasta')
        try:
            node_file_as_list = read_defined_file_to_list(node_file_path)
        except FileNotFoundError:
            raise RuntimeError({'med_output_directory':self.output_directory})

        for i in range(0, len(node_file_as_list), 2):
            node_seq_name = node_file_as_list[i].split('|')[0][1:]
            node_seq_abundance = int(node_file_as_list[i].split('|')[1].split(':')[1])
            node_seq_sequence = node_file_as_list[i+1].replace('-', '')
            self.nodes_list_of_nucleotide_sequences.append(
                NucleotideSequence(name=node_seq_name, abundance=node_seq_abundance, sequence=node_seq_sequence))

    def make_data_set_sample_sequences(self):
        for node_nucleotide_sequence_object in self.nodes_list_of_nucleotide_sequences:
            if not self._assign_node_sequence_to_existing_ref_seq(node_nucleotide_sequence_object):
                self._assign_node_sequence_to_new_ref_seq(node_nucleotide_sequence_object)

        if self._two_or_more_nodes_associated_to_the_same_reference_sequence():
            self._make_associations_and_abundances_in_node_abund_df_unique_again()

        self._update_data_set_sample_med_qc_meta_data_and_clade_totals()

        self._if_more_than_200_seqs_create_clade_collection()

        self._create_data_set_sample_sequences()

    def _create_data_set_sample_sequences(self):
        if self._we_made_a_clade_collection():
            self._create_data_set_sample_sequences_with_clade_collection()
        else:
            self._create_data_set_sample_sequences_without_clade_collection()

    def _create_data_set_sample_sequences_without_clade_collection(self):
        data_set_sample_sequence_list = []
        for node_nucleotide_sequence_object in self.nodes_list_of_nucleotide_sequences:
            associated_ref_seq_id = self.node_sequence_name_to_ref_seq_id[node_nucleotide_sequence_object.name]
            associated_ref_seq_object = ReferenceSequence.objects.get(id=associated_ref_seq_id)
            df_index_label = self.node_abundance_df.index.values.tolist()[0]
            dss = DataSetSampleSequence(reference_sequence_of=associated_ref_seq_object,
                                        abundance=self.node_abundance_df.at[
                                            df_index_label, node_nucleotide_sequence_object.name],
                                        data_set_sample_from=self.dataset_sample_object)
            data_set_sample_sequence_list.append(dss)
        DataSetSampleSequence.objects.bulk_create(data_set_sample_sequence_list)

    def _create_data_set_sample_sequences_with_clade_collection(self):
        data_set_sample_sequence_list = []
        associated_ref_seq_uid_as_str_list = []
        for node_nucleotide_sequence_object in self.nodes_list_of_nucleotide_sequences:
            associated_ref_seq_id = self.node_sequence_name_to_ref_seq_id[node_nucleotide_sequence_object.name]
            associated_ref_seq_uid_as_str_list.append(str(associated_ref_seq_id))
            associated_ref_seq_object = ReferenceSequence.objects.get(id=associated_ref_seq_id)
            df_index_label = self.node_abundance_df.index.values.tolist()[0]
            dss = DataSetSampleSequence(
                reference_sequence_of=associated_ref_seq_object,
                clade_collection_found_in=self.clade_collection_object,
                abundance=self.node_abundance_df.at[df_index_label, node_nucleotide_sequence_object.name],
                data_set_sample_from=self.dataset_sample_object)
            data_set_sample_sequence_list.append(dss)
        # Save all of the newly created dss
        DataSetSampleSequence.objects.bulk_create(data_set_sample_sequence_list)
        self.clade_collection_object.footprint = ','.join(associated_ref_seq_uid_as_str_list)
        self.clade_collection_object.save()

    def _we_made_a_clade_collection(self):
        return self.total_num_sequences > 200

    def _if_more_than_200_seqs_create_clade_collection(self):
        if self.total_num_sequences > 200:
            sys.stdout.write(
                f'\n{self.sample_name} clade {self.clade}: '
                f'{self.total_num_sequences} sequences. Creating CladeCollection_object\n')
            new_cc = CladeCollection(clade=self.clade, data_set_sample_from=self.dataset_sample_object)
            new_cc.save()
            self.clade_collection_object = new_cc
        else:
            sys.stdout.write(
                f'\n{self.sample_name} clade {self.clade}: {self.total_num_sequences} sequences. '
                f'Insufficient sequence to create a CladeCollection_object\n')

    def _update_data_set_sample_med_qc_meta_data_and_clade_totals(self):
        self.dataset_sample_object.post_med_absolute += self.total_num_sequences
        self.dataset_sample_object.post_med_unique += len(self.node_abundance_df.iloc[0])
        # Update the cladal_seq_totals
        cladal_seq_abundance_counter = [int(a) for a in json.loads(self.dataset_sample_object.cladal_seq_totals)]
        clade_index = list('ABCDEFGHI').index(self.clade)
        cladal_seq_abundance_counter[clade_index] = self.total_num_sequences
        self.dataset_sample_object.cladal_seq_totals = json.dumps([str(a) for a in cladal_seq_abundance_counter])
        self.dataset_sample_object.save()

    def _two_or_more_nodes_associated_to_the_same_reference_sequence(self):
        """Multiple nodes may be assigned to the same reference sequence. We only want to create one
        data set sample sequence object per reference sequence. As such we need to consolidate abundances
        """
        return len(
            set(self.node_sequence_name_to_ref_seq_id.values())) != len(self.node_sequence_name_to_ref_seq_id.keys())

    def _make_associations_and_abundances_in_node_abund_df_unique_again(self):
        list_of_non_unique_ref_seq_uids = [
            ref_seq_uid for ref_seq_uid, count in
            Counter(self.node_sequence_name_to_ref_seq_id.values()).items() if count > 1]
        for non_unique_ref_seq_uid in list_of_non_unique_ref_seq_uids:

            node_names_to_be_consolidated = self._get_list_of_node_names_that_need_consolidating(non_unique_ref_seq_uid)

            summed_abund_of_nodes_of_ref_seq = self._get_summed_abundances_of_the_nodes_to_be_consolidated(
                node_names_to_be_consolidated)

            self._del_all_but_first_of_non_unique_nodes_from_df_and_node_to_ref_seq_dict(node_names_to_be_consolidated)

            self._update_node_name_abund_in_df_and_nuc_seq_list(node_names_to_be_consolidated, summed_abund_of_nodes_of_ref_seq)

    def _get_list_of_node_names_that_need_consolidating(self, non_unique_ref_seq_uid):
        node_names_to_be_consolidated = [
            node_name for node_name in self.node_sequence_name_to_ref_seq_id.keys() if
            self.node_sequence_name_to_ref_seq_id[node_name] == non_unique_ref_seq_uid]
        return node_names_to_be_consolidated

    def _update_node_name_abund_in_df_and_nuc_seq_list(self, node_names_to_be_consolidated, summed_abund_of_nodes_of_ref_seq):
        df_index_name = self.node_abundance_df.index.values.tolist()[0]
        self.node_abundance_df.at[df_index_name, node_names_to_be_consolidated[0]] = summed_abund_of_nodes_of_ref_seq

        # del consolidated seqs and update abundances from the self.nodes_list_of_nucleotide_sequences
        new_list_of_node_nucleotide_sequences = []
        for i in range(len(self.nodes_list_of_nucleotide_sequences)):
            # update abundance of the representative sequences
            if self.nodes_list_of_nucleotide_sequences[i].name == node_names_to_be_consolidated[0]:
                self.nodes_list_of_nucleotide_sequences[i].abundance = summed_abund_of_nodes_of_ref_seq
            # add nodes not being consolidated to new node nuc seq list
            if self.nodes_list_of_nucleotide_sequences[i].name not in node_names_to_be_consolidated[1:]:
                new_list_of_node_nucleotide_sequences.append(self.nodes_list_of_nucleotide_sequences[i])
        self.nodes_list_of_nucleotide_sequences = new_list_of_node_nucleotide_sequences

    def _del_all_but_first_of_non_unique_nodes_from_df_and_node_to_ref_seq_dict(self, node_names_to_be_consolidated):
        # now delete columns for all but the first of the node_names_to_be_consolidated from df
        self.node_abundance_df = self.node_abundance_df.drop(columns=node_names_to_be_consolidated[1:])
        # and from the node_names_to_ref_seq_id dict
        for node_name in node_names_to_be_consolidated[1:]:
            del self.node_sequence_name_to_ref_seq_id[node_name]

    def _get_summed_abundances_of_the_nodes_to_be_consolidated(self, node_names_to_be_consolidated):
        summed_abund_of_nodes_of_ref_seq = sum(self.node_abundance_df[node_names_to_be_consolidated].iloc[0])
        return summed_abund_of_nodes_of_ref_seq

    def _assign_node_sequence_to_existing_ref_seq(self, node_nucleotide_sequence_object):
        """ We use this to look to see if there is an equivalent ref_seq Sequence for the sequence in question
        This take into account whether the seq_in_q could be a subset or super set of one of the
        ref_seq.sequences.
        Will return false if no ref_seq match is found
        """
        if self._node_sequence_exactly_matches_reference_sequence_sequence(node_nucleotide_sequence_object):
            return self._associate_node_seq_to_ref_seq_by_exact_match_and_return_true(node_nucleotide_sequence_object)
        elif self._node_sequence_matches_reference_sequence_sequence_plus_adenine(node_nucleotide_sequence_object):
            # This was a seq shorter than refseq but we can associate
            return self._associate_node_seq_to_ref_seq_by_adenine_match_and_return_true(node_nucleotide_sequence_object)
        else:
            return self._search_for_super_set_match_and_associate_if_found_else_return_false(
                node_nucleotide_sequence_object)

    def _node_sequence_exactly_matches_reference_sequence_sequence(self, node_nucleotide_sequence_object):
        return node_nucleotide_sequence_object.sequence in self.ref_seq_sequence_to_ref_seq_id_dict

    def _node_sequence_matches_reference_sequence_sequence_plus_adenine(self, node_nucleotide_sequence_object):
        return 'A' + node_nucleotide_sequence_object.sequence in self.ref_seq_sequence_to_ref_seq_id_dict

    def _search_for_super_set_match_and_associate_if_found_else_return_false(self, node_nucleotide_sequence_object):
        # or if the seq in question is bigger than a refseq sequence and is a super set of it
        # In either of these cases we should consider this a match and use the refseq matched to.
        # This might be very coputationally expensive but lets give it a go
        for ref_seq_sequence in self.ref_seq_sequence_to_ref_seq_id_dict.keys():
            if node_nucleotide_sequence_object.sequence in ref_seq_sequence or \
                    ref_seq_sequence in node_nucleotide_sequence_object.sequence:
                # Then this is a match
                self.node_sequence_name_to_ref_seq_id[node_nucleotide_sequence_object.name] = \
                    self.ref_seq_sequence_to_ref_seq_id_dict[ref_seq_sequence]
                name_of_reference_sequence = self.ref_seq_uid_to_ref_seq_name_dict[
                    self.ref_seq_sequence_to_ref_seq_id_dict[ref_seq_sequence]]
                self._print_succesful_association_details_to_stdout(node_nucleotide_sequence_object,
                                                                    name_of_reference_sequence)
                return True
        return False

    def _associate_node_seq_to_ref_seq_by_adenine_match_and_return_true(self, node_nucleotide_sequence_object):
        self.node_sequence_name_to_ref_seq_id[
            node_nucleotide_sequence_object.name] = self.ref_seq_sequence_to_ref_seq_id_dict[
            'A' + node_nucleotide_sequence_object.sequence]
        name_of_reference_sequence = self.ref_seq_uid_to_ref_seq_name_dict[
            self.ref_seq_sequence_to_ref_seq_id_dict['A' + node_nucleotide_sequence_object.sequence]]
        self._print_succesful_association_details_to_stdout(node_nucleotide_sequence_object, name_of_reference_sequence)
        return True

    def _associate_node_seq_to_ref_seq_by_exact_match_and_return_true(self, node_nucleotide_sequence_object):
        self.node_sequence_name_to_ref_seq_id[
            node_nucleotide_sequence_object.name] = self.ref_seq_sequence_to_ref_seq_id_dict[
            node_nucleotide_sequence_object.sequence]
        name_of_reference_sequence = self.ref_seq_uid_to_ref_seq_name_dict[
            self.ref_seq_sequence_to_ref_seq_id_dict[node_nucleotide_sequence_object.sequence]]
        self._print_succesful_association_details_to_stdout(node_nucleotide_sequence_object, name_of_reference_sequence)
        return True

    def _print_succesful_association_details_to_stdout(
            self, node_nucleotide_sequence_object, name_of_reference_sequence):
        sys.stdout.write(f'\r{self.sample_name} clade {self.clade}: '
                         f'Assigning MED node {node_nucleotide_sequence_object.name} '
                         f'to existing reference sequence {name_of_reference_sequence}')

    def _assign_node_sequence_to_new_ref_seq(self, node_nucleotide_sequence_object):
        new_ref_seq = ReferenceSequence(clade=self.clade, sequence=node_nucleotide_sequence_object.sequence)
        new_ref_seq.save()
        new_ref_seq.name = str(new_ref_seq.id)
        new_ref_seq.save()
        self.ref_seq_sequence_to_ref_seq_id_dict[new_ref_seq.sequence] = new_ref_seq.id
        self.node_sequence_name_to_ref_seq_id[node_nucleotide_sequence_object.name] = new_ref_seq.id
        self.ref_seq_uid_to_ref_seq_name_dict[new_ref_seq.id] = new_ref_seq.name

        sys.stdout.write(f'\r{self.sample_name} clade {self.clade}: '
                         f'Assigning MED node {node_nucleotide_sequence_object.name} '
                         f'to new reference sequence {new_ref_seq.name}')


class DataSetSampleCreatorHandler:
    """This class will be where we run the code for creating reference sequences, data set sample sequences and
    clade collections."""
    def __init__(self):
        # dictionaries to save us having to do lots of database look ups
        self.ref_seq_uid_to_ref_seq_name_dict = {
            ref_seq.id: ref_seq.name for ref_seq in ReferenceSequence.objects.all()}
        self.ref_seq_sequence_to_ref_seq_id_dict = {
            ref_seq.sequence: ref_seq.id for ref_seq in ReferenceSequence.objects.all()}

    def execute_data_set_sample_creation(
            self, data_loading_list_of_med_output_directories, data_loading_debug, data_loading_dataset_object):
        for med_output_directory in data_loading_list_of_med_output_directories:
            try:
                data_set_sample_sequence_creator_worker = DataSetSampleSequenceCreatorWorker(
                    med_output_directory=med_output_directory,
                    data_loading_dataset_obj=data_loading_dataset_object,
                    data_set_sample_creator_handler_ref_seq_sequence_to_ref_seq_id_dict=
                    self.ref_seq_sequence_to_ref_seq_id_dict,
                    data_set_sample_creator_handler_ref_seq_uid_to_ref_seq_name_dict=
                    self.ref_seq_uid_to_ref_seq_name_dict)
            except RuntimeError as e:
                non_existant_med_output_dir = e.args[0]['med_output_directory']
                print(f'{non_existant_med_output_dir}: File not found during DataSetSample creation.')
                continue
            if data_loading_debug:
                if data_set_sample_sequence_creator_worker.num_med_nodes < 10:
                    print(
                        f'{med_output_directory}: '
                        f'WARNING node file contains only {data_set_sample_sequence_creator_worker.num_med_nodes} sequences.')
            sys.stdout.write(
                f'\n\nPopulating {data_set_sample_sequence_creator_worker.sample_name} with clade {data_set_sample_sequence_creator_worker.clade} sequences\n')
            data_set_sample_sequence_creator_worker.make_data_set_sample_sequences()