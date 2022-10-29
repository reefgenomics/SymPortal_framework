from dbApp.models import (
    DataSet, ReferenceSequence, DataSetSample, DataSetSampleSequence,
    CladeCollection, DataSetSampleSequencePM)
import sys
import os
import shutil
import subprocess
import pandas as pd
import json
from collections import Counter
from django import db
from multiprocessing import Queue as mp_Queue, Manager, Process, Lock as mp_Lock
from threading import Lock as mt_Lock, Thread, get_ident
from queue import Queue as mt_Queue
from general import ThreadSafeGeneral, file_as_blockiter, hash_bytestr_iter
from datetime import datetime
import distance
from plotting import DistScatterPlotterSamples, SeqStackedBarPlotter
from symportal_utils import BlastnAnalysis, MothurAnalysis, NucleotideSequence
from output import SequenceCountTableCreator
import ntpath
import math
from numpy import NaN
from collections import defaultdict
import itertools
import time
from shutil import which
import sp_config
from django_general import CreateStudyAndAssociateUsers
import logging
import hashlib
from general import check_lat_lon
import re
from calendar import month_abbr, month_name
from psycopg2 import InterfaceError


class DataLoading:
    # The clades refer to the phylogenetic divisions of the Symbiodiniaceae. Most of them are represented at the genera
    # level. E.g. All members of clade C belong to the genus Cladocopium.
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    def __init__(
            self, parent_work_flow_obj, user_input_path, datasheet_path,
            screen_sub_evalue, num_proc,no_fig, no_ord, no_output,
            distance_method, no_pre_med_seqs, multiprocess, start_time, date_time_str, is_cron_loading,
            study_name=None, study_user_string=None,
            debug=False):
        self.parent = parent_work_flow_obj
        self.is_cron_loading = is_cron_loading
        self.thread_safe_general = ThreadSafeGeneral()
        # check and generate the sample_meta_info_df first before creating the DataSet object
        self.sample_meta_info_df = None
        self.user_input_path = user_input_path
        self.datasheet_path = datasheet_path
        self._check_mothur_version()
        if self.datasheet_path:
            self._create_and_check_datasheet()
        self.symportal_root_directory = os.path.abspath(os.path.dirname(__file__))
        self.dataset_object = None
        # the stability file generated here is used as the base of the initial mothur QC
        self.list_of_samples_names = []
        self.list_of_fastq_files_in_wkd = []
        self.sample_fastq_pairs = None
        if self.datasheet_path:
            self._get_sample_names_and_create_new_dataset_object_with_datasheet()
        else:
            end_index = self._get_sample_names_and_create_new_dataset_object_without_datasheet()

        self.num_proc = min(num_proc, len(self.list_of_samples_names))
        self.temp_working_directory = self._setup_temp_working_directory()
        self.date_time_str = date_time_str
        self.output_directory = self._setup_output_directory()
        logging.basicConfig(format='%(levelname)s:%(message)s',
                            filename=os.path.join(self.output_directory, f'{self.date_time_str}_log.log'), filemode='w',
                            level=logging.INFO)
        logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
        # directory for the data_explorer outputs
        self.html_dir = os.path.join(self.output_directory, 'html')
        self.js_file_path = os.path.join(self.html_dir, 'study_data.js')
        os.makedirs(self.html_dir, exist_ok=True)
        # dictionary that will hold the outputfile type to full path of the outputfile
        self.js_output_path_dict = {}
        if self.datasheet_path:
            self._generate_stability_file_and_data_set_sample_objects_with_datasheet()
        else:
            # Dictionary where sample name is value and seq file full paths in fwd, rev order are in list
            self.sample_name_to_seq_files_dict = dict()
            self._generate_stability_file_and_data_set_sample_objects_without_datasheet(end_index)
        self.list_of_dss_objects = DataSetSample.objects.filter(data_submission_from=self.dataset_object)
        if sp_config.system_type == 'remote':
            csaau = CreateStudyAndAssociateUsers(
                date_time_str=self.date_time_str, ds=self.dataset_object,
                list_of_dss_objects=self.list_of_dss_objects, is_cron_loading=self.is_cron_loading,
                study_name=study_name, study_user_string=study_user_string
            )
            csaau.create_study_and_user_objects()
            self.study = csaau.study
        self.output_path_list = []
        self.no_fig = no_fig
        self.no_ord = no_ord
        self.no_output = no_output
        self.distance_method = distance_method
        self.no_pre_med_seqs = no_pre_med_seqs
        # this is the path of the file we will use to deposit a backup copy of the reference sequences
        self.seq_dump_file_path = self._setup_sequence_dump_file_path()
        self.dataset_object.working_directory = self.temp_working_directory
        self.dataset_object.save()
        # This is the directory that sequences that have undergone QC for each sample will be written out as
        # .names and .fasta pairs BEFORE MED decomposition
        # We will delete the directory if it already exists

        self.pre_med_sequence_output_directory_path = self._create_pre_med_write_out_directory_path()
        
        # directory that will contain sub directories for each sample. Each sub directory will contain a pair of
        # .names and .fasta files of the non_symbiodiniaceae_sequences that were thrown out for that sample
        self.non_symb_and_size_violation_base_dir_path = os.path.join(
            self.output_directory, 'non_sym_and_size_violation_sequences'
        )
        os.makedirs(self.non_symb_and_size_violation_base_dir_path, exist_ok=True)
        # data can be loaded either as paired fastq or fastq.gz files or as a single compressed file containing
        # the afore mentioned paired files.
        self.is_single_file_or_paired_input = self._determine_if_single_file_or_paired_input()
        self.debug = debug
        self.symclade_db_directory_path = os.path.abspath(os.path.join(
            self.symportal_root_directory, 'symbiodiniaceaeDB'))
        self.symclade_db_full_path = os.path.join(self.symclade_db_directory_path, 'symClade.fa')

        self.path_to_mothur_batch_file_for_dot_file_creation = None
        self.path_to_latest_mothur_batch_file = None

        self.samples_that_caused_errors_in_qc_list = []
        self.initial_mothur_handler = None
        self.post_initial_qc_name_file_name = None
        self.post_initial_qc_fasta_file_name = None
        # args for the taxonomic screening
        self.screen_sub_evalue = screen_sub_evalue
        self.new_seqs_added_in_iteration = 0
        self.new_seqs_added_running_total = 0
        self.checked_samples_with_no_additional_symbiodiniaceae_sequences = []
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
        self.required_symbiodiniaceae_matches = 3
        # med
        self.list_of_med_output_directories = []
        self.path_to_med_padding_executable = os.path.join(
            self.symportal_root_directory, 'lib/med_decompose/o_pad_with_gaps.py')
        self.path_to_med_decompose_executable = os.path.join(
            self.symportal_root_directory, 'lib/med_decompose/decompose.py')
        self.perform_med_handler_instance = None
        # data set sample creation
        self.data_set_sample_creator_handler_instance = None
        # plotting sequence output from both post-med and pre-med seq outputs
        self.seq_abundance_relative_output_path_post_med = None
        self.seq_abundance_relative_output_path_pre_med = None
        self.seq_abund_relative_df_post_med = None
        self.seq_abund_relative_df_pre_med = None
        # we will use this sequence count table creator when outputting the pre_MED seqs so that the df
        # can be put in the same order
        self.sequence_count_table_creator = None
        # we will use this sequence stacked bar plotter when plotting the pre_MED seqs so that the plotting
        # can be put in the same order
        self.seq_stacked_bar_plotter = None
        # path to executable for multithreading the checking pre-MED sequences against existing ReferenceSequences
        self.path_to_seq_match_executable = os.path.join(
            self.symportal_root_directory, 'seq_match.py')
        # Timers
        # The timers for meausring how long it takes to create the DataSetSampleSequencePM
        self.pre_med_seq_start_time = None
        self.pre_med_seq_stop_time = None
        self.multiprocess = multiprocess
        self.start_time = start_time

    def load_data(self):
        self._copy_and_decompress_input_files_to_temp_wkd()

        self._if_symclade_binaries_not_present_remake_db()

        self._do_initial_mothur_qc()

        self._taxonomic_screening()

        self._do_med_decomposition()

        self._create_data_set_sample_sequences_from_med_nodes()
        if not self.no_pre_med_seqs:
            self._create_data_set_sample_sequence_pre_med_objs()
        else:
            print('\n\nSkipping generation of pre med seq objects at users request\n\n')

        self._print_sample_successful_or_failed_summary()

        self._perform_sequence_drop()

        self._delete_temp_dir__log_files_and_pre_med_dir()

        self._write_data_set_info_to_stdout()

        if not self.no_output:
            self._output_seqs_count_table()

            self._write_sym_non_sym_and_size_violation_dirs_to_stdout()

            self._output_seqs_stacked_bar_plots()

            self._do_sample_ordination()

            # finally write out the dict that holds the output file paths for the DataExplorer
            # covert the full paths to relative paths and then write out dict
            # https://stackoverflow.com/questions/8693024/how-to-remove-a-path-prefix-in-python
            new_dict = {}
            for k, v in self.js_output_path_dict.items():
                new_dict[k] = os.path.relpath(v, self.output_directory)

            self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
                [{'function_name': 'getDataFilePaths', 'python_obj': new_dict}],
                js_outpath=self.js_file_path)

        print('\n\nDATA LOADING COMPLETE')
        print(f'DataSet id: {self.dataset_object.id}')
        print(f'DataSet name: {self.dataset_object.name}')
        self.dataset_object.loading_complete_time_stamp = str(
            datetime.utcnow()).split('.')[0].replace('-','').replace(' ','T').replace(':','')
        self.dataset_object.save()
        print(f'Loading completed in {time.time() - self.start_time}s')
        print(f'DataSet loading_complete_time_stamp: {self.dataset_object.loading_complete_time_stamp}\n\n\n')
        print(f"Log written to {os.path.join(self.output_directory, f'{self.date_time_str}_log.log')}")

    def _check_mothur_version(self):
        mothur_version_cmd = subprocess.run(
            ['mothur', '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        for line in self.thread_safe_general.decode_utf8_binary_to_list(mothur_version_cmd.stdout):
            if "1.43" in line:
                return
        raise RuntimeError('SymPortal currently uses version 1.43 of mothur.\nCheck your version.')

    def _make_new_dataset_object(self):
        self.dataset_object = DataSet(
            name=self.parent.args.name, time_stamp=self.parent.date_time_str,
            reference_fasta_database_used=self.parent.reference_db,
            submitting_user=self.parent.submitting_user,
            submitting_user_email=self.parent.submitting_user_email)
        self.dataset_object.save()
        self.parent.data_set_object = self.dataset_object

    def _write_sym_non_sym_and_size_violation_dirs_to_stdout(self):
        if not self.no_pre_med_seqs:
            print(f'\nPre-MED Symbiodiniaceae sequences written out to:\n'
                  f'{self.pre_med_sequence_output_directory_path}')
        print(f'\nNon-Symbiodiniaceae and size violation sequences written out to:\n'
              f'{self.non_symb_and_size_violation_base_dir_path}')

    def _do_sample_ordination(self):
        if not self.no_ord:
            self._do_sample_dist_and_pcoa()

            self._plot_pcoa()

    def _plot_pcoa(self):
        if not self.no_fig:
            if len(self.list_of_samples_names) > 1000:
                print(f'Too many samples ({len(self.list_of_samples_names)}) to generate plots')
            else:
                for output_path in self.output_path_list:
                    if self.this_is_pcoa_path(output_path):
                        clade_of_output = os.path.dirname(output_path).split('/')[-1]
                        sys.stdout.write(f'\n\nGenerating between sample distance plot clade {clade_of_output}\n')
                        try:
                            dist_scatter_plotter_samples = DistScatterPlotterSamples(csv_path=output_path,
                                                                                     date_time_str=self.date_time_str)
                            dist_scatter_plotter_samples.make_sample_dist_scatter_plot()
                        except RuntimeError:
                            # The error message is printed to stdout at the source
                            continue
                        self.output_path_list.extend(dist_scatter_plotter_samples.output_path_list)

    def _do_sample_dist_and_pcoa(self):
        print('\nCalculating between sample pairwise distances')
        if self.distance_method == 'both':
            if self._check_if_required_packages_found_in_path():
                self._do_unifrac_dist_pcoa()
                self._do_braycurtis_dist_pcoa()
            else:
                print('Changing distance method to braycurtis as one or more of the required '
                      'packages could not be found in your PATH')
                self.distance_method = 'braycurtis'
        if self.distance_method == 'unifrac':
            self._do_unifrac_dist_pcoa()
        elif self.distance_method == 'braycurtis':
            self._do_braycurtis_dist_pcoa()

    def _check_if_required_packages_found_in_path(self):
        """For creating unifrac-derived distances we need
        both iqtree and mafft to be install in the users PATH.
        Here we will check for them. If either of them are not found we will return False"""
        return which('iqtree') and which('mafft')

    def _write_data_set_info_to_stdout(self):
        print(f'\n\nData loading complete. DataSet UID: {self.dataset_object.id}; Name: {self.dataset_object.name}')

    @staticmethod
    def this_is_pcoa_path(output_path):
        return 'PCoA_coords' in output_path

    def _do_braycurtis_dist_pcoa(self):
        braycurtis_dist_pcoa_creator = distance.SampleBrayCurtisDistPCoACreator(
            date_time_str=self.date_time_str,
            data_set_uid_list=[self.dataset_object.id],
            output_dir=self.output_directory, html_dir=self.html_dir,
            js_output_path_dict=self.js_output_path_dict)
        braycurtis_dist_pcoa_creator.compute_braycurtis_dists_and_pcoa_coords()
        self.output_path_list.extend(braycurtis_dist_pcoa_creator.output_path_list)

    def _do_unifrac_dist_pcoa(self):
        unifrac_dict_pcoa_creator = distance.SampleUnifracDistPCoACreator(
            date_time_str=self.date_time_str, output_dir=self.output_directory,
            data_set_uid_list=[self.dataset_object.id], num_processors=self.num_proc,
            html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict)
        unifrac_dict_pcoa_creator.compute_unifrac_dists_and_pcoa_coords()
        self.output_path_list.extend(unifrac_dict_pcoa_creator.output_path_list)

    def _output_seqs_stacked_bar_plots(self):
        """Plot up the post- and pre-MED seqs as both .png and .svg"""
        if not self.no_fig:
            if len(self.list_of_samples_names) > 1000:
                print(f'Too many samples ({len(self.list_of_samples_names)}) to generate plots')
            else:
                sys.stdout.write('\nGenerating sequence count table figures\n')

                self.seq_stacked_bar_plotter = SeqStackedBarPlotter(
                    output_directory=self.output_directory,
                    seq_relative_abund_count_table_path_post_med=self.seq_abundance_relative_output_path_post_med,
                    no_pre_med_seqs=self.no_pre_med_seqs,
                    ordered_seq_list=self.sequence_count_table_creator.clade_abundance_ordered_ref_seq_list,
                    date_time_str=self.date_time_str,
                    seq_relative_abund_df_pre_med=self.seq_abund_relative_df_pre_med
                    )
                self.seq_stacked_bar_plotter.plot_stacked_bar_seqs()
                self.output_path_list.extend(self.seq_stacked_bar_plotter.output_path_list)

    def _output_seqs_count_table(self):
        sys.stdout.write('\nGenerating count tables for post- and pre-MED sequence abundances\n')
        self.sequence_count_table_creator = SequenceCountTableCreator(
            symportal_root_dir=self.symportal_root_directory, call_type='submission',
            no_pre_med_seqs=self.no_pre_med_seqs, ds_uids_output_str=str(self.dataset_object.id),
            num_proc=self.num_proc, date_time_str=self.date_time_str,
            html_dir=self.html_dir,
            js_output_path_dict=self.js_output_path_dict, multiprocess=self.multiprocess)
        self.sequence_count_table_creator.make_seq_output_tables()
        self.seq_abund_relative_df_post_med = self.sequence_count_table_creator.output_df_relative_post_med
        self.output_path_list.extend(self.sequence_count_table_creator.output_paths_list)
        self._set_seq_abundance_relative_output_path(self.sequence_count_table_creator)
        self.seq_abund_relative_df_pre_med = self.sequence_count_table_creator.output_df_relative_pre_med
        self.seq_abundance_relative_output_path_pre_med = self.sequence_count_table_creator.pre_med_relative_df_path

    def _set_seq_abundance_relative_output_path(self, sequence_count_table_creator):
        for path in sequence_count_table_creator.output_paths_list:
            if 'relative.abund_and_meta' in path:
                self.seq_abundance_relative_output_path_post_med = path

    def _delete_temp_dir__log_files_and_pre_med_dir(self):
        if os.path.exists(self.temp_working_directory):
            shutil.rmtree(self.temp_working_directory)
        # Delete any log files that are found anywhere in the SymPortal directory
        subprocess.run(f'find {self.symportal_root_directory} -name "*.logfile" -delete', shell=True, check=True)
        # Delete the pre_med_seq directories holding the .fasta and .names pairs
        # as this information is now stored in the database and output in a count table.
        # This directory will be remade if outputs are being made
        if os.path.exists(self.pre_med_sequence_output_directory_path):
            shutil.rmtree(self.pre_med_sequence_output_directory_path)

    def _perform_sequence_drop(self):
        sequence_drop_file = self._generate_sequence_drop_file()
        sys.stdout.write(f'\n\nBackup of named reference_sequences output to {self.seq_dump_file_path}\n')
        self.thread_safe_general.write_list_to_destination(self.seq_dump_file_path, sequence_drop_file)

    @staticmethod
    def _generate_sequence_drop_file():
        header_string = ','.join(['seq_name', 'seq_uid', 'seq_clade', 'seq_sequence'])
        sequence_drop_list = [header_string]
        for ref_seq in ReferenceSequence.objects.filter(has_name=True):
            sequence_drop_list.append(','.join([ref_seq.name, str(ref_seq.id), ref_seq.clade, ref_seq.sequence]))
        return sequence_drop_list

    def _print_sample_successful_or_failed_summary(self):
        logging.info('SAMPLE PROCESSING SUMMARY')
        failed_count = 0
        successful_count = 0
        for data_set_sample in DataSetSample.objects.filter(data_submission_from=self.dataset_object):
            if data_set_sample.error_in_processing:
                failed_count += 1
                logging.info(f'{data_set_sample.name}: Error in processing: {data_set_sample.error_reason}')
            else:
                successful_count += 1
                logging.info(f'{data_set_sample.name}: Successful')

        logging.info(f'\n\n{successful_count} out of {successful_count + failed_count} '
                     f'samples successfully passed QC.\n'
                     f'{failed_count} samples produced errors\n')

    def _create_data_set_sample_sequences_from_med_nodes(self):
        self.data_set_sample_creator_handler_instance = DataSetSampleCreatorHandler()
        self.data_set_sample_creator_handler_instance.execute_data_set_sample_creation(
            data_loading_list_of_med_output_directories=self.list_of_med_output_directories,
            data_loading_debug=self.debug, data_loading_dataset_object=self.dataset_object)
        self.dataset_object.currently_being_processed = False
        self.dataset_object.save()

    def _create_data_set_sample_sequence_pre_med_objs(self):
        print('\n\nCreating DataSetSampleSequencePM objects')
        self.pre_med_seq_start_time = time.time()
        data_set_sample_pre_med_obj_creator = FastDataSetSampleSequencePMCreator(
            dataset_object=self.dataset_object,
            pre_med_sequence_output_directory_path=self.pre_med_sequence_output_directory_path,
            num_proc=self.num_proc, path_to_seq_match_executable=self.path_to_seq_match_executable,
            temp_working_directory=self.temp_working_directory)
        data_set_sample_pre_med_obj_creator.make_data_set_sample_pm_objects()
        self.pre_med_seq_stop_time = time.time()
        print(f'\n\nCreation of DataSetSampleSequencePM objects took '
              f'{self.pre_med_seq_stop_time-self.pre_med_seq_start_time}s')

    def _do_med_decomposition(self):
        self.perform_med_handler_instance = PerformMEDHandler(
            data_loading_temp_working_directory=self.temp_working_directory,
            data_loading_num_proc=self.num_proc,
            multiprocess=self.multiprocess)

        self.perform_med_handler_instance.execute_perform_med_worker(
            data_loading_debug=self.debug,
            data_loading_path_to_med_decompose_executable=self.path_to_med_decompose_executable,
            data_loading_path_to_med_padding_executable=self.path_to_med_padding_executable)

        self.list_of_med_output_directories = self.perform_med_handler_instance.list_of_med_result_dirs

        if self.debug:
            print('MED dirs:')
            for med_output_directory in self.list_of_med_output_directories:
                print(med_output_directory)

    def _do_initial_mothur_qc(self):

        if not self.sample_fastq_pairs:
            self._exit_and_del_data_set_sample('Sample fastq pairs list empty')

        self.initial_mothur_handler = InitialMothurHandler(data_loading_parent=self)
        self.initial_mothur_handler.execute_worker_initial_mothur()
        self.samples_that_caused_errors_in_qc_list = list(
            self.initial_mothur_handler.samples_that_caused_errors_in_qc_mp_list)

    def _taxonomic_screening(self):
        """
        There are only two things to achieve as part of this taxonomy screening.
        1 - identify the sequences that are Symbiodinium/Symbiodiniaceae in origin
        2 - identify the sequences that are non-Symbiodinium/non-Symbiodiniaceae in origin.
        This may sound very straight forwards: Blast each sequence against a reference database and if we get
        a good enough match, consider this sequence symbiodiniaceae in origin. If not, not.
        BUT, there is a catch. We cannot be sure that every new sequence that we receive is not symbiodiniaceae just
        because we don't get a good match to our reference database. It may be that new diversity is not yet
        represented in our reference database.
        As a results of this we want to do an extra set of taxonomic screening and that is what this first part of code
        concerns. We will run a blast against our reference symClade database. And then, any seuqences that return
        a match to a member of this database, but does not meet the minimum threshold to be directly considered
        Symbiodinium in origin (i.e. not similar enough to a reference sequence in the symClade database) will
        be blasted against the NCBI nt blast database. For a sequence to be blasted against this database, and be
        in contesion for being considered symbiodiniaceae in origin it must also be found in at least three samples.
        We will use an iterative screening process to acheive this symbiodiniaceae identification. The sequences to be
        screened will be referred to as subevalue sequences. Any of these sequences that we deem symbiodiniaceae in
        origin after running against the nt database will be added back into the symClade reference database.
        On the next iteration it will therefore be possible for differnt sequences to be matches given that additional
        sequences may have been added to this reference database. This first part of screening will only be to run
        the symClade blast, and to screen low identity matches. Non-symbiodiniaceae sequence matches will be made in
        the second part of this screening.
        To make this faster, rather than check every sample in every iteration we will keep track of when a sample
        has been found to only contain sequences that are good matches to the symClade database.
        This way we can skip these sequences in the iterative rounds of screening.
        """

        if self.screen_sub_evalue:
            if not len(self.samples_that_caused_errors_in_qc_list) == self.list_of_samples_names:
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
            data_loading_samples_that_caused_errors_in_qc_mp_list=self.samples_that_caused_errors_in_qc_list,
            multiprocess=self.multiprocess, data_loading_dataset_object=self.dataset_object
        )
        self.sym_non_sym_tax_screening_handler.execute_sym_non_sym_tax_screening(
            data_loading_temp_working_directory=self.temp_working_directory,
            data_loading_pre_med_sequence_output_directory_path=self.pre_med_sequence_output_directory_path,

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
        We will call a sequence symbiodiniaceae if it has a match that covers at least 95% of its sequence
        at a 60% or higher
        identity. It must also have Symbiodinium or Symbiodiniaceae in the name. We will also require that a
        sub_e_value seq has at least the required_sybiodinium_matches (3 at the moment) before we call it Symbiodinium.
        """

        blastn_analysis_object = BlastnAnalysis(
            input_file_path=self.sequences_to_screen_fasta_path,
            output_file_path=os.path.join(self.temp_working_directory, 'blast.out'),
            max_target_seqs=10,
            num_threads=str(self.num_proc), pipe_stdout_sterr=(not self.debug)
        )

        blastn_analysis_object.execute_blastn_analysis()

        blast_output_dict = blastn_analysis_object.return_blast_results_dict()

        query_sequences_verified_as_symbiodiniaceae_list = self._get_list_of_seqs_in_blast_result_that_are_symbiodiniaceae(
            blast_output_dict
        )

        self.new_seqs_added_in_iteration = len(query_sequences_verified_as_symbiodiniaceae_list)
        self.new_seqs_added_running_total += self.new_seqs_added_in_iteration

        if query_sequences_verified_as_symbiodiniaceae_list:
            self._taxa_screening_update_symclade_db_with_new_symbiodiniaceae_seqs(
                query_sequences_verified_as_symbiodiniaceae_list)

    def _taxa_screening_update_symclade_db_with_new_symbiodiniaceae_seqs(
            self, query_sequences_verified_as_symbiodiniaceae_list):
        new_symclade_fasta_as_list = self._taxa_screening_make_new_fasta_of_screened_seqs_to_be_added_to_symclade_db(
            query_sequences_verified_as_symbiodiniaceae_list
        )
        combined_fasta = self._taxa_screening_combine_new_symclade_seqs_with_current(new_symclade_fasta_as_list)
        self._taxa_screening_make_new_symclade_db(combined_fasta)

    def _taxa_screening_make_new_symclade_db(self, combined_fasta):
        self.thread_safe_general.write_list_to_destination(self.symclade_db_full_path, combined_fasta)
        self.thread_safe_general.make_new_blast_db(
            input_fasta_to_make_db_from=self.symclade_db_full_path, db_title='symClade')

    def _taxa_screening_combine_new_symclade_seqs_with_current(self, new_symclade_fasta_as_list):
        old_symclade_fasta_as_list = self.thread_safe_general.read_defined_file_to_list(self.symclade_db_full_path)
        combined_fasta = new_symclade_fasta_as_list + old_symclade_fasta_as_list
        return combined_fasta

    def _taxa_screening_make_new_fasta_of_screened_seqs_to_be_added_to_symclade_db(
            self, query_sequences_verified_as_symbiodiniaceae_list):
        screened_seqs_fasta_dict = self.thread_safe_general.create_dict_from_fasta(
            fasta_path=self.sequences_to_screen_fasta_path
        )
        new_symclade_fasta_as_list = []
        for name_of_symbiodiniaceae_sequence_to_add_to_symclade_db in query_sequences_verified_as_symbiodiniaceae_list:
            new_symclade_fasta_as_list.extend(
                [
                    f'>{name_of_symbiodiniaceae_sequence_to_add_to_symclade_db}',
                    f'{screened_seqs_fasta_dict[name_of_symbiodiniaceae_sequence_to_add_to_symclade_db]}'
                ]
            )
        return new_symclade_fasta_as_list

    def _get_list_of_seqs_in_blast_result_that_are_symbiodiniaceae(self, blast_output_dict):
        query_sequences_verified_as_symbiodiniaceae_list = []
        for query_sequence_name, blast_result_list_for_query_sequence in blast_output_dict.items():
            sym_count = 0
            for result_str in blast_result_list_for_query_sequence:
                if 'Symbiodinium' in result_str or 'Symbiodiniaceae' in result_str:
                    percentage_coverage = float(result_str.split('\t')[4])
                    percentage_identity_match = float(result_str.split('\t')[3])
                    if percentage_coverage > 95 and percentage_identity_match > 60:
                        sym_count += 1
                        if sym_count == self.required_symbiodiniaceae_matches:
                            query_sequences_verified_as_symbiodiniaceae_list.append(query_sequence_name)
                            break
        return query_sequences_verified_as_symbiodiniaceae_list

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
            self.thread_safe_general.write_list_to_destination(
                self.sequences_to_screen_fasta_path, self.sequences_to_screen_fasta_as_list)

    def _create_symclade_backup_incase_of_accidental_deletion_of_corruption(self):
        back_up_dir = os.path.abspath(os.path.join(self.symportal_root_directory, 'symbiodiniaceaeDB', 'symClade_backup'))
        os.makedirs(back_up_dir, exist_ok=True)
        symclade_current_path = os.path.abspath(
            os.path.join(self.symportal_root_directory, 'symbiodiniaceaeDB', 'symClade.fa'))

        symclade_backup_path = os.path.join(back_up_dir, f'symClade_{self.date_time_str}.fa')
        symclade_backup_readme_path = os.path.join(back_up_dir, f'symClade_{self.date_time_str}.readme')
        # then write a copy to it.
        shutil.copy(symclade_current_path, symclade_backup_path)
        # Then write out a very breif readme
        read_me = [
            f'This is a symClade.fa backup created during datasubmission of data_set ID: {self.dataset_object.id}']
        self.thread_safe_general.write_list_to_destination(symclade_backup_readme_path, read_me)

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
        self.checked_samples_with_no_additional_symbiodiniaceae_sequences = \
            list(self.taxonomic_screening_handler.checked_samples_mp_list)

    def _init_potential_sym_tax_screen_handler(self):

        self.taxonomic_screening_handler = PotentialSymTaxScreeningHandler(
            samples_that_caused_errors_in_qc_list=self.samples_that_caused_errors_in_qc_list,
            checked_samples_list=self.checked_samples_with_no_additional_symbiodiniaceae_sequences,
            list_of_samples_names=self.list_of_samples_names, num_proc=self.num_proc, multiprocess=self.multiprocess
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
                self.thread_safe_general.make_new_blast_db(
                    input_fasta_to_make_db_from=self.symclade_db_full_path,
                    db_title='symClade', pipe_stdout_sterr=True)
            elif self.debug:
                self.thread_safe_general.make_new_blast_db(
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

    def _get_sample_names_and_create_new_dataset_object_without_datasheet(self):
        for file in os.listdir(self.user_input_path):
            if file.endswith('fastq') or file.endswith('fq') or file.endswith('fastq.gz') or file.endswith('fq.gz'):
                self.list_of_fastq_files_in_wkd.append(file)

        if len(self.list_of_fastq_files_in_wkd) < 3:
            raise RuntimeError(
                f'Cannot auto infer names from {len(self.list_of_fastq_files_in_wkd)} fastq files. '
                f'Please use a datasheet.')

        end_index = self._identify_sample_names_without_datasheet()

        self._make_new_dataset_object()

        return end_index

    def _generate_stability_file_and_data_set_sample_objects_without_datasheet(self, end_index):

        self.make_dot_stability_file_inferred(end_index)

        self._create_data_set_sample_objects_in_bulk_without_datasheet()

    def make_dot_stability_file_inferred(self, end_index):
        """Search for the fastq files that contain the inferred sample names. NB this is not so simple
        as names that are subset of other names will match more than one set of fastq files.
        After identifying the sample direction, write absolute paths of each of the seq files
        as though they are coming from the temp_working_directory. This is where the raw sequencing files
        will be copied over to."""
        print('\nDeducing read direction for the inferred sample names')
        sample_fastq_pairs = []
        for sample_name in self.list_of_samples_names:
            print(f'Sample {sample_name}')
            temp_list = [sample_name.replace('-', '[dS]')]
            fwd_file_path = None
            rev_file_path = None
            for file_path in self.thread_safe_general.return_list_of_file_paths_in_directory(self.user_input_path):
                if sample_name == ntpath.basename(file_path)[:-end_index]:
                    # When here we know which sample the file_path belongs to
                    # but we still need to deduce whether this is the fwd or rev read
                    # If R1 or R2 are in the read, then that is relatively easy
                    # But it may be that there is only a 1 or a 2 in the read.
                    # To test for this we will parse through each of the 1's or 2's
                    # and see if the partner read exists

                    # First search for the simple case of R1 or R2 being present.
                    if 'R1' in file_path or 'R2' in file_path:
                        if 'R1' in file_path:
                            fwd_file_path = file_path
                            self._print_sample_direction_path(direction='fwd', file_path=file_path)
                        if 'R2' in file_path:
                            rev_file_path = file_path
                            self._print_sample_direction_path(direction='rev', file_path=file_path)
                    else:
                        seq_direction_result = self._check_seq_file_path_for_seq_direction(file_path=file_path)
                        if seq_direction_result == 'rev':
                            rev_file_path = file_path
                            self._print_sample_direction_path(direction='rev', file_path=file_path)
                        elif seq_direction_result == 'fwd':
                            fwd_file_path = file_path
                            self._print_sample_direction_path(direction='fwd', file_path=file_path)
                        else:
                            raise RuntimeError(f'Unable to deduce read direction of {file_path}')
                # if we have already found a fwd_file_path and rev_file_path
                # then we can break out of the search
                if fwd_file_path and rev_file_path:
                    break
            # Only add the files to the sample_fastq_pairs if they are above the minium size requirement
            # Check that both files meet the required minimum file size
            if os.path.getsize(fwd_file_path) > 300 and os.path.getsize(rev_file_path) > 300:
                # If so, add them to the dictionary
                self.sample_name_to_seq_files_dict[sample_name] = [fwd_file_path, rev_file_path]
                # Change the current paths so that they don't originate from the temp_working_directory
                # path but rather from the temp working directory
                fwd_file_path = os.path.join(self.temp_working_directory, ntpath.basename(fwd_file_path))
                rev_file_path = os.path.join(self.temp_working_directory, ntpath.basename(rev_file_path))
                temp_list.append(fwd_file_path)
                temp_list.append(rev_file_path)
                if None in temp_list:
                    raise RuntimeError(f'Error in deducing directionality of {sample_name}')
                sample_fastq_pairs.append('\t'.join(temp_list))
            else:
                print(f'WARNING: At least one of the seq files for sample {sample_name} is less than 300 bytes in size')
                print(f'{sample_name} will not be included in the dataloading')
        # Reinit the list_of_samples_names so from the sample_name_to_seq_files_dict so that
        # the samples that had files that were below the 300 byte size threshold are removed form the list
        self.list_of_samples_names = list(self.sample_name_to_seq_files_dict.keys())
        self.thread_safe_general.write_list_to_destination(
            r'{0}/stability.files'.format(self.temp_working_directory), sample_fastq_pairs)
        self.sample_fastq_pairs = sample_fastq_pairs

    @staticmethod
    def _check_seq_file_path_for_seq_direction(file_path):
        """ This will try to deduce whether a given sequencing file is of a given direction"""
        file_name_list = list(ntpath.basename(file_path))

        # for each of the either '1' or '2' s in the file name
        # swap them out for the opposite number and check to see if this
        # file exists. If it does, then we assume that this is the '1' or '2' that we
        # are investigating is the indicator of seq direction
        for char_index, char_element in enumerate(reversed(file_name_list)):
            if char_element == '1':
                search_list = list(reversed(file_name_list))
                search_list[char_index] = '2'
                search_filename = ''.join(reversed(search_list))
                search_path = os.path.join(os.path.dirname(file_path), search_filename)
                if os.path.exists(search_path):
                    return 'fwd'
            elif char_element == '2':
                search_list = list(reversed(file_name_list))
                search_list[char_index] = '1'
                search_filename = ''.join(reversed(search_list))
                search_path = os.path.join(os.path.dirname(file_path), search_filename)
                if os.path.exists(search_path):
                    return 'rev'
        return False

    @staticmethod
    def _print_sample_direction_path(direction, file_path):
        if direction == 'fwd':
            print(f'R1 = {file_path}')
        else:
            print(f'R2 = {file_path}')

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
        for dss_chunk in self.thread_safe_general.chunks(list_of_sample_objects):
            DataSetSample.objects.bulk_create(dss_chunk)

    def _get_num_chars_in_common_with_fastq_names(self):
        i = 1
        while 1:
            list_of_endings = set()
            for file in self.list_of_fastq_files_in_wkd:
                list_of_endings.add(file[-i:])
            if len(list_of_endings) > 2:
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
            warning_str = 'Error in automatic sample name extraction. ' \
                          'Please explicitly supply sample names using a data sheet ' \
                          '(https://github.com/didillysquat/SymPortal_framework/wiki/Running-SymPortal#loading-data)'
            sys.exit(warning_str)
        self.list_of_samples_names = list_of_sample_names

    def _identify_sample_names_without_datasheet(self):
        # I think the simplest way to get sample names is to find what parts are common between all samples
        # well actually 50% of the samples so that we also remove the R1 and R2 parts.
        end_index = self._get_num_chars_in_common_with_fastq_names()
        self._get_sample_names_from_fastq_files_using_index(end_index)
        return end_index

    def _get_sample_names_and_create_new_dataset_object_with_datasheet(self):
        self.list_of_samples_names = self.sample_meta_info_df.index.values.tolist()

        self._make_new_dataset_object()

    def _generate_stability_file_and_data_set_sample_objects_with_datasheet(self):
        # if we are given a data_sheet then use the sample names given as the DataSetSample object names
        self.make_dot_stability_file_datasheet()
        self._create_data_set_sample_objects_in_bulk_with_datasheet()

    def make_dot_stability_file_datasheet(self):
        """Create the .stability file that mothur will use to make contigs.
        This file is the sample name a tab, the fwd full path, a tab, the rev full path.
        The paths should refer to the files that will have been copied over to the temp_working_directory.
        As such the file paths will be the file name of the current file path value in the info_df
        joined with the temp_working_directory."""
        sample_fastq_pairs = []
        for sample_name in self.sample_meta_info_df.index.values.tolist():
            temp_list = [sample_name.replace('-', '[dS]')]
            temp_list.append(
                os.path.join(
                    self.temp_working_directory,
                    ntpath.basename(self.sample_meta_info_df.loc[sample_name, 'fastq_fwd_file_name'])
                )
            )
            temp_list.append(
                os.path.join(
                    self.temp_working_directory,
                    ntpath.basename(self.sample_meta_info_df.loc[sample_name, 'fastq_rev_file_name'])
                )
            )
            sample_fastq_pairs.append('\t'.join(temp_list))

        self.thread_safe_general.write_list_to_destination(
            os.path.join(self.temp_working_directory, 'stability.files'), sample_fastq_pairs)
        self.sample_fastq_pairs = sample_fastq_pairs

    def _create_data_set_sample_objects_in_bulk_with_datasheet(self):
        """
        The proper formatting of the values in the df should already have been taken care of. However, I will
        leave the code in below as a safe guard to make sure that the DataSetSample objects can still be successfuly
        created.
        """
        list_of_data_set_sample_objects = []
        sys.stdout.write('\nCreating data_set_sample objects\n')
        for sampleName in self.list_of_samples_names:
            print('\rCreating data_set_sample {}'.format(sampleName))
            # Create the data_set_sample objects in bulk.
            # The cladal_seq_totals property of the data_set_sample object keeps track of the seq totals for the
            # sample divided by clade. This is used in the output to keep track of sequences that are not
            # included in cladeCollections
            empty_cladal_seq_totals = json.dumps([0 for _ in self.clade_list])

            try:
                sample_type = str(self.sample_meta_info_df.loc[sampleName, 'sample_type'])
                host_phylum = str(self.sample_meta_info_df.loc[sampleName, 'host_phylum'])
                host_class = str(self.sample_meta_info_df.loc[sampleName, 'host_class'])
                host_order = str(self.sample_meta_info_df.loc[sampleName, 'host_order'])
                host_family = str(self.sample_meta_info_df.loc[sampleName, 'host_family'])
                host_genus = str(self.sample_meta_info_df.loc[sampleName, 'host_genus'])
                host_species = str(self.sample_meta_info_df.loc[sampleName, 'host_species'])
                collection_depth = str(self.sample_meta_info_df.loc[sampleName, 'collection_depth'])
                collection_date = str(self.sample_meta_info_df.loc[sampleName, 'collection_date'])
            except:
                sample_type = 'NoData'
                host_phylum = 'NoData'
                host_class = 'NoData'
                host_order = 'NoData'
                host_family = 'NoData'
                host_genus = 'NoData'
                host_species = 'NoData'
                collection_depth = 'NoData'
                collection_date = 'NoData'
            try:
                collection_latitude = float(self.sample_meta_info_df.loc[sampleName, 'collection_latitude'])
                collection_longitude = float(self.sample_meta_info_df.loc[sampleName, 'collection_longitude'])
                if math.isnan(collection_latitude) or math.isnan(collection_longitude):
                    collection_latitude = float(999)
                    print('conversion problem with collection_latitude or collection_longitude, converting both to 999')
                    collection_longitude = float(999)
            except:
                print('conversion problem with collection_latitude or collection_longitude, converting both to 999')
                collection_latitude = float(999)
                collection_longitude = float(999)

            # get the sha256 hash of the fastq.gz files
            fwd_hash = hash_bytestr_iter(file_as_blockiter(open(self.sample_meta_info_df.loc[sampleName, 'fastq_fwd_file_name'], 'rb')), hashlib.sha256(), True)
            rev_hash = hash_bytestr_iter(file_as_blockiter(open(self.sample_meta_info_df.loc[sampleName, 'fastq_rev_file_name'], 'rb')), hashlib.sha256(), True)

            dss = DataSetSample(name=sampleName, data_submission_from=self.dataset_object,
                                cladal_seq_totals=empty_cladal_seq_totals,
                                sample_type=sample_type,
                                host_phylum=host_phylum,
                                host_class=host_class,
                                host_order=host_order,
                                host_family=host_family,
                                host_genus=host_genus,
                                host_species=host_species,
                                collection_latitude=collection_latitude,
                                collection_longitude=collection_longitude,
                                collection_date=collection_date,
                                collection_depth=collection_depth,
                                fastq_fwd_file_name=ntpath.basename(self.sample_meta_info_df.loc[sampleName, 'fastq_fwd_file_name']),
                                fastq_rev_file_name=ntpath.basename(self.sample_meta_info_df.loc[sampleName, 'fastq_rev_file_name']),
                                fastq_fwd_file_hash=fwd_hash,
                                fastq_rev_file_hash=rev_hash
                                )
            list_of_data_set_sample_objects.append(dss)
        # http://stackoverflow.com/questions/18383471/django-bulk-create-function-example
        for dss_chunk in self.thread_safe_general.chunks(list_of_data_set_sample_objects):
            DataSetSample.objects.bulk_create(dss_chunk)

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

    def _create_and_check_datasheet(self):
        # Create a pandas df from the data_sheet if it was provided
        # allow the data_sheet to be in a .csv format or .xlsx format. This is so that we can store a datasheet
        # in the github repo in a non-binary format
        # The sample_meta_df that is created from the data_sheet should be identical irrespective of whether a .csv
        # or a .xlsx is submitted.
        # 1 - Things to check: Check that each of the sample names are unique.
        # 1c - remove any spaces and replace with '_'
        # 1b - Check that there are no white spaces either side of the sample names.
        # 2 - Check that each of the seq files are unique and exist
        # 3 - Check that the lat lon file exists and is in a format that can be used. If not convert to 999.
        # 4 - check that each of the other values can be converted to str

        # Datasheet to pandas df
        self._read_in_datasheet()

        # drop any cells in which the sample name is null
        self.sample_meta_info_df = self.sample_meta_info_df[~pd.isnull(self.sample_meta_info_df['sample_name'])]

        self._check_datasheet_df_vals_unique()

        self.sample_meta_info_df['sample_name'] = self.sample_meta_info_df['sample_name'].astype(str)
        self.sample_meta_info_df['sample_name'] = self.sample_meta_info_df['sample_name'].str.rstrip()\
            .str.lstrip().str.replace(' ', '_').str.replace('/', '_').str.replace('α', 'alpha').str.replace('β', 'beta')

        self.sample_meta_info_df.set_index('sample_name', inplace=True, drop=True)

        self.sample_meta_info_df.index = self.sample_meta_info_df.index.map(str)

        self._check_for_binomial()

        self._replace_null_vals_in_meta_info_df()

        self._check_seq_files_exist()

        self._check_lat_long()

        self._check_vars_can_be_string()

        self._check_date_format()

    def _check_date_format(self):
        """
        Try to coerce some of the common date format errors into YYYYMMDD format
        Common inputs will be DD.MM.YYYY, DD.MM.YY, MM.YYYY or MM.YY
        We shold throw an error if there is some
        other format given (i.e. one that is not just '.' and integers)
        """
        bad_formats = []
        lower_month_abbr = [_.lower() for _ in month_abbr][1:]
        lower_month_name = [_.lower() for _ in month_name][1:]
        for row_name in self.sample_meta_info_df.index:
            current_date_value = self.sample_meta_info_df.at[row_name, 'collection_date']
            if current_date_value == "NoData":
                continue
            if not pd.isnull(current_date_value):
                # sometime a weird float string was coming in e.g. 20220601.0
                try:
                    current_date_value = str(int(float(current_date_value)))
                except ValueError:
                    pass
                if re.findall("[A-Za-z]+", current_date_value):
                    # Then this date is in a bad format
                    # We will try to extract a year and a month from it
                    # We will assume that the year is in YYYY format
                    # and that the month is in either in the common abbreviation form
                    # or written out in form.
                    putative_months = re.findall("[A-Za-z]+", current_date_value)
                    if len(putative_months) == 1:
                        putative_month = putative_months[0].lower()
                        if putative_month in lower_month_abbr:
                            month_ind = lower_month_abbr.index(putative_month) + 1
                        elif putative_month in lower_month_name:
                            month_ind = lower_month_name.index(putative_month) + 1
                        else:
                            # not recognised so log as error
                            bad_formats.append((row_name, current_date_value))
                            continue
                        # If we got here then we have month_id
                        if month_ind < 10:
                            month_ind = f"0{month_ind}"
                        else:
                            month_ind = str(month_ind)

                        # Then we need to pull out the year
                        if len(re.findall("[0-9]{4}", current_date_value)) == 1:
                            year = re.findall("[0-9]{4}", current_date_value)[0]
                        else:
                            bad_formats.append((row_name, current_date_value))
                            continue

                        # Finally check that there is nothing less after we remove the month and the year
                        remaining = current_date_value.lower().replace(putative_month, "").replace(year, "").rstrip()
                        if remaining == "":
                            # Then we can convert
                            new_date_value = f"{year}{month_ind}"
                            print(f'changing {current_date_value} to {new_date_value} for {row_name}')
                            self.sample_meta_info_df.at[row_name, 'collection_date'] = new_date_value
                            continue
                        else:
                            # There is something left and we call it bad
                            bad_formats.append((row_name, current_date_value))
                            continue
                    else:
                        # Add to the bad_formats list so that we can print out and
                        # exit at the end of this
                        bad_formats.append((row_name, current_date_value))
                        continue
                elif "." in current_date_value:
                    if current_date_value.count(".") == 2:
                        # Then this is DD.MM.YYYY or DD.MM.YY
                        new_date_value = current_date_value.replace(".", "")
                        if len(new_date_value) == 6:
                            new_date_value = ''.join(["20", new_date_value[4:], new_date_value[2:4], new_date_value[:2]])
                            print(f'changing {current_date_value} to {new_date_value} for {row_name}')
                            self.sample_meta_info_df.at[row_name, 'collection_date'] = new_date_value
                        elif len(new_date_value) == 8:
                            new_date_value = ''.join([new_date_value[4:], new_date_value[2:4], new_date_value[:2]])
                            print(f'changing {current_date_value} to {new_date_value} for {row_name}')
                            self.sample_meta_info_df.at[row_name, 'collection_date'] = new_date_value
                        else:
                            bad_formats.append((row_name, current_date_value))
                    elif current_date_value.count(".") == 1:
                        # Then this is MM.YY or MM.YYYY
                        new_date_value = current_date_value.replace(".", "")
                        if len(new_date_value) == 4:
                            new_date_value = ''.join(["20", new_date_value[2:], new_date_value[:2]])
                            print(f'changing {current_date_value} to {new_date_value} for {row_name}')
                            self.sample_meta_info_df.at[row_name, 'collection_date'] = new_date_value
                        elif len(new_date_value) == 6:
                            new_date_value = ''.join([new_date_value[2:], new_date_value[:2]])
                            print(f'changing {current_date_value} to {new_date_value} for {row_name}')
                            self.sample_meta_info_df.at[row_name, 'collection_date'] = new_date_value
                        else:
                            bad_formats.append((row_name, current_date_value))
                    else:
                        bad_formats.append((row_name, current_date_value))

                elif len(re.findall("[0-9]{8}", current_date_value)) == 1:
                    # Then this is good: YYYYMMDD
                    continue
                elif len(re.findall("[0-9]{6}", current_date_value)) == 1:
                    # Then this is good: YYYYMM
                    continue
                elif len(re.findall("^[0-9]{4}$", current_date_value)) == 1:
                    # Then this is good: YYYY
                    continue
                else:
                    # Else, something else is going on
                    bad_formats.append((row_name, current_date_value))
            else:
                continue
        if bad_formats:
            print("There are errors in the date_collection formats")
            print("Date format should be YYYYMMDD or YYYYMM")
            for bad_sample, bad_val in bad_formats:
                print(f"{bad_sample}: {bad_val}")
            sys.exit()

    def _read_in_datasheet(self):
        if self.datasheet_path.endswith('.xlsx'):
            self.sample_meta_info_df = pd.read_excel(
                io=self.datasheet_path, header=0, usecols='A:N', skiprows=[0], dtype=str)
        elif self.datasheet_path.endswith('.csv'):
            with open(self.datasheet_path, 'r') as f:
                data_sheet_as_file = [line.rstrip() for line in f]
            if data_sheet_as_file[0].split(',')[0] == 'sample_name':
                self.sample_meta_info_df = pd.read_csv(
                    filepath_or_buffer=self.datasheet_path, dtype=str)
            else:
                self.sample_meta_info_df = pd.read_csv(
                    filepath_or_buffer=self.datasheet_path, skiprows=[0], dtype=str)
        else:
            sys.exit('Data sheet: {} is in an unrecognised format. '
                     'Please ensure that it is either in .xlsx or .csv format.')

    def _replace_null_vals_in_meta_info_df(self):
        self.sample_meta_info_df = self.sample_meta_info_df.replace('N/A', NaN).replace('NA', NaN).replace('na',
                                                                                                           NaN).replace(
            'n/a', NaN)

    def _check_for_binomial(self):
        """People were putting the full binomial in the speices colums. This crops this back to just the
        species component of binomial.
        
        It will only adjust this value if there are two components to the value and if the first component
        matches the genus. I.e. only if the user has put in a genuine binomial. To preven people entering
        things like 'sp. 1' and it being corrected to '1'.
        """
        for row_name in self.sample_meta_info_df.index.values.tolist():
            current_species_val = self.sample_meta_info_df.at[row_name, 'host_species']
            if not pd.isnull(current_species_val):
                components = current_species_val.split(' ')
                if len(components) == 2:
                    if components[0] == self.sample_meta_info_df.at[row_name, 'host_genus']:
                        new_species_val = current_species_val.split(' ')[1]
                        print(f'changing {current_species_val} to {new_species_val} for {row_name}')
                        self.sample_meta_info_df.at[row_name, 'host_species'] = new_species_val

    def _check_vars_can_be_string(self):
        """First convert each of the columns to type string.
        Then make sure that all of the vals are genuine vals or NoData
        """
        for col in ['sample_type', 'host_phylum', 'host_class', 'host_order', 'host_family', 'host_genus',
                    'host_species', 'collection_depth', 'collection_date']:
            self.sample_meta_info_df[col] = self.sample_meta_info_df[col].astype(str)
        for i, sample_name in enumerate(self.sample_meta_info_df.index.values.tolist()):
            for col in ['sample_type', 'host_phylum', 'host_class', 'host_order', 'host_family', 'host_genus',
                        'host_species', 'collection_depth', 'collection_date']:
                try:
                    value = str(self.sample_meta_info_df.at[sample_name, col])
                    if value == 'nan':
                        self.sample_meta_info_df.at[sample_name, col] = 'NoData'
                except:
                    self.sample_meta_info_df.at[sample_name, col] = 'NoData'

    def _check_lat_long(self):
        # check the lat long value for each sample listed
        for i, sample_name in enumerate(self.sample_meta_info_df.index.values.tolist()):
            lat = self.sample_meta_info_df.at[sample_name, 'collection_latitude']
            lon = self.sample_meta_info_df.at[sample_name, 'collection_longitude']

            try:
                lat_float, lon_float = check_lat_lon(lat, lon)
            except Exception:
                print(f'Lat and long are currently nan for {sample_name}. Values will be set to 999')
                self._set_lat_lon_to_999(sample_name)
                continue

            # final check to make sure that the values are in a sensible range
            if (-90 <= lat_float <= 90) and (-180 <= lon_float <= 180):
                self.sample_meta_info_df.at[sample_name, 'collection_latitude'] = lat_float
                self.sample_meta_info_df.at[sample_name, 'collection_longitude'] = lon_float
            else:
                self._set_lat_lon_to_999(sample_name)

        # finally make sure that the lat and long cols are typed as float
        self.sample_meta_info_df['collection_latitude'] = self.sample_meta_info_df['collection_latitude'].astype(float)
        self.sample_meta_info_df['collection_longitude'] = self.sample_meta_info_df['collection_longitude'].astype(
            float)

    def _set_lat_lon_to_999(self, sample_name):
        self.sample_meta_info_df.at[sample_name, 'collection_latitude'] = float(999)
        self.sample_meta_info_df.at[sample_name, 'collection_longitude'] = float(999)

    def _check_seq_files_exist(self):
        """
        Check that all of the sequencing files provided in the datasheet exist.
        If filenames are given in the datasheet then convert these to full paths
        before checking for their existence. This way, all locations of seq files
        are full paths.
        Ensure that we lstrip and rstrip the entries to remove any spaces.
        Also check for size of file and require a 300B minimum. Remove from the sample from the data sheet if
        smaller than this.
        """

        self.sample_meta_info_df['fastq_fwd_file_name'] = self.sample_meta_info_df['fastq_fwd_file_name'].astype(str)
        self.sample_meta_info_df['fastq_fwd_file_name'] = self.sample_meta_info_df['fastq_fwd_file_name'].str.rstrip() \
            .str.lstrip()
        self.sample_meta_info_df['fastq_rev_file_name'] = self.sample_meta_info_df['fastq_rev_file_name'].astype(str)
        self.sample_meta_info_df['fastq_rev_file_name'] = self.sample_meta_info_df['fastq_rev_file_name'].str.rstrip() \
            .str.lstrip()

        file_not_found_list = []
        rows_to_drop = []
        for df_ind in self.sample_meta_info_df.index.values.tolist():
            fwd_file = self.sample_meta_info_df.at[df_ind, 'fastq_fwd_file_name']
            rev_file = self.sample_meta_info_df.at[df_ind, 'fastq_rev_file_name']
            if '/' not in fwd_file:
                fwd_file_path = os.path.join(self.user_input_path, fwd_file)
                self.sample_meta_info_df.at[df_ind, 'fastq_fwd_file_name'] = fwd_file_path
            else:
                fwd_file_path = fwd_file
            if '/' not in rev_file:
                rev_file_path = os.path.join(self.user_input_path, rev_file)
                self.sample_meta_info_df.at[df_ind, 'fastq_rev_file_name'] = rev_file_path
            else:
                rev_file_path = rev_file

            # Check for the fwd read
            if not os.path.exists(fwd_file_path):
                # Take into account that the user might have supplied fastq.gz files but
                # used a .fastq extension. If this is the case. Correct in the df.
                if os.path.exists(fwd_file_path + '.gz'):
                    self.sample_meta_info_df.at[df_ind, 'fastq_fwd_file_name'] = fwd_file_path + '.gz'
                    fwd_file_path = self.sample_meta_info_df.at[df_ind, 'fastq_fwd_file_name']
                # Also take into account that user may have supplied no extension
                elif os.path.exists(fwd_file_path + '.fastq.gz'):
                    self.sample_meta_info_df.at[df_ind, 'fastq_fwd_file_name'] = fwd_file_path + '.fastq.gz'
                    fwd_file_path = self.sample_meta_info_df.at[df_ind, 'fastq_fwd_file_name']
                else:
                    file_not_found_list.append(fwd_file_path)
            # Check for rev read
            if not os.path.exists(rev_file_path):
                # Take into account that the user might have supplied fastq.gz files but
                # used a .fastq extension. If this is the case. Correct in the df.
                if os.path.exists(rev_file_path + '.gz'):
                    self.sample_meta_info_df.at[df_ind, 'fastq_rev_file_name'] = rev_file_path + '.gz'
                    rev_file_path = self.sample_meta_info_df.at[df_ind, 'fastq_rev_file_name']
                # Also take into account that user may have supplied no extension
                elif os.path.exists(rev_file_path + '.fastq.gz'):
                    self.sample_meta_info_df.at[df_ind, 'fastq_rev_file_name'] = rev_file_path + '.fastq.gz'
                    rev_file_path = self.sample_meta_info_df.at[df_ind, 'fastq_rev_file_name']
                else:
                    file_not_found_list.append(rev_file_path)
            # NB if we were unable to find either the fwd or rev read then we will not be able
            # to continue with our checks
            if  rev_file_path in file_not_found_list or fwd_file_path in file_not_found_list:
                continue

            # Check file size
            if os.path.getsize(fwd_file_path) < 300 or os.path.getsize(rev_file_path) < 300:
                print(f'WARNING: At least one of the seq files for sample {df_ind} is less than 300 bytes in size')
                print(f'{df_ind} will be removed from your datasheet and analysis')
                rows_to_drop.append(df_ind)

        # drop the rows that had size violations
        self.sample_meta_info_df.drop(index=rows_to_drop, inplace=True)

        if file_not_found_list:
            print('Some of the sequencing files listed in your datasheet cannot be found:')
            for file_name in file_not_found_list:
                print(f'{file_name}')
            sys.exit()

    def _check_datasheet_df_vals_unique(self):
        # check sample names
        self._check_vals_of_col_unique(column_name='sample_name')
        # check fastq_fwd_file_name
        self._check_vals_of_col_unique(column_name='fastq_fwd_file_name')
        # check fastq_rev_file_name
        self._check_vals_of_col_unique(column_name='fastq_rev_file_name')

    def _check_vals_of_col_unique(self, column_name):
        # check to see that the values held in a column are unique
        sample_name_counter = Counter(self.sample_meta_info_df[column_name].values.tolist())
        non_unique_name_list = []
        for col_val, count in sample_name_counter.items():
            if count != 1:
                non_unique_name_list.append(col_val)
        if non_unique_name_list:
            print(f'There appear to be non unique {column_name}s in your datasheet:')
            for col_val in non_unique_name_list:
                print(f'\t{col_val}')
            print('Please rectify this in your datasheet and reattempt loading')
            sys.exit()

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
        fq.gz rather than fastq and fq.

        Previously we were only allowing files that were all contained in the user_input_path directory.
        However, we are now introducing functionality that allows full paths to be specified in the datasheet
        and these paths should be copied over to the temp_working directory.
        """
        if self.datasheet_path:
            # If working form datasheet this is as easy as copying over the specified file paths
            for sample_name in self.sample_meta_info_df.index.values.tolist():
                # copy over the fwd read
                shutil.copy(self.sample_meta_info_df.loc[sample_name, 'fastq_fwd_file_name'],
                            self.temp_working_directory)
                # copy over the rev read
                shutil.copy(self.sample_meta_info_df.loc[sample_name, 'fastq_rev_file_name'],
                            self.temp_working_directory)
        else:
            # If not working from datasheet then we transfer over all files of the right extension
            for file_path, fwd_rev_list in self.sample_name_to_seq_files_dict.items():
                shutil.copy(fwd_rev_list[0], self.temp_working_directory)
                shutil.copy(fwd_rev_list[1], self.temp_working_directory)

    def _determine_if_single_file_or_paired_input(self):
        for file in os.listdir(self.user_input_path):
            if 'fastq' in file or 'fq' in file:
                # Then there is a fastq.gz or fastq file already uncompressed in this folder
                # In this case we will assume that the seq data is not a single file containing the pairs of files
                # rather the pairs of files themselves.
                return False

    def _create_pre_med_write_out_directory_path(self):
        pre_med_write_out_directory_path = os.path.join(self.output_directory, 'pre_med_seqs')
        if os.path.exists(pre_med_write_out_directory_path):
            shutil.rmtree(pre_med_write_out_directory_path)
        os.makedirs(pre_med_write_out_directory_path)
        return pre_med_write_out_directory_path

    def _setup_output_directory(self):
        output_directory = os.path.join(self.symportal_root_directory,
                                        'outputs', 'loaded_data_sets', f'{self.dataset_object.id}', self.date_time_str)
        os.makedirs(output_directory, exist_ok=True)
        return output_directory

    def _setup_sequence_dump_file_path(self):
        seq_dump_file_path = os.path.join(
            self.symportal_root_directory, 'dbBackUp', 'seq_dumps', f'seq_dump_{self.date_time_str}')
        os.makedirs(os.path.dirname(seq_dump_file_path), exist_ok=True)
        return seq_dump_file_path

    def _setup_temp_working_directory(self):
        """
        Create a working directory that will be housed in a temp folder within the self.user_input_path directory
        We will copy all sequencing files over to this directory before undertaking QC. This temporary directory
        will be deleted upon completion of data loading.
        """
        if '.' in self.user_input_path.split('/')[-1]:
            # then this path points to a file rather than a directory and we should pass through the path only
            self.temp_working_directory = os.path.abspath(
                os.path.join(
                    os.path.dirname(self.user_input_path), 'tempData', str(self.dataset_object.id)
                )
            )
        else:
            # then we assume that we are pointing to a directory and we can directly use that to make the wkd
            self.temp_working_directory = os.path.abspath(os.path.join(self.user_input_path, 'tempData',
                                                                       str(self.dataset_object.id)))
        self._create_temp_wkd()
        return self.temp_working_directory

    def _create_temp_wkd(self):
        # if the directory already exists remove it and start from scratch
        if os.path.exists(self.temp_working_directory):
            shutil.rmtree(self.temp_working_directory)
        os.makedirs(self.temp_working_directory)


class FastDataSetSampleSequencePMCreator:
    def __init__(
            self, pre_med_sequence_output_directory_path, dataset_object, num_proc,
            path_to_seq_match_executable, temp_working_directory):
        # dictionaries to save us having to do lots of database look ups
        self.path_to_seq_match_executable = path_to_seq_match_executable
        self.pre_med_sequence_output_directory_path = pre_med_sequence_output_directory_path
        self.thread_safe_general = ThreadSafeGeneral()
        self.num_proc = num_proc
        self.dataset_object = dataset_object
        clades = list('ABCDEFGHI')
        self.ref_seq_sequence_to_ref_seq_obj_dict = {}
        for clade in clades:
            self.ref_seq_sequence_to_ref_seq_obj_dict[clade] = {
                ref_seq.sequence: ref_seq for ref_seq in ReferenceSequence.objects.filter(clade=clade)}
        self.list_of_pre_med_sample_dirs = self._populate_list_of_pre_med_sample_dirs()
        # This is a dict that will have three levels.
        # The first set of keys will be the clades.
        # The second set of keys will be nucleotide sequences
        # The third set of keys will be DataSetSample objects and the
        # absolute abundnace of the sequence as the value.
        self.consolidated_sequence_to_sample_and_abund_dict = dict()
        self._populated_consolidated_seq_to_sample_and_abund_dict()
        # Now work through the consolidated dictionary matching sequences to reference sequences
        # sequences for which reference sequences have been found will be represented in a new dictionary
        # for which the key will be a reference sequence object that the match was found and the value
        # will be the same dict value that was in the original consolidated dicitonary
        # If a match is not found, then move these sequences into a second dictionary for the non matches
        # This is the dictionary where the matched sequence info will go
        # Each of these two dictionaries will be set within a layer of clade keys
        self.ref_seq_match_obj_to_seq_sample_abundance_dict = defaultdict(dict)
        # This is the dictionary where the non matched sequences will be put
        self.no_match_consolidated_seq_to_sample_and_abund_dict = defaultdict(dict)
        self.temp_working_directory = temp_working_directory

    def _populate_list_of_pre_med_sample_dirs(self):
        return self.thread_safe_general.return_list_of_directory_paths_in_directory(
            self.pre_med_sequence_output_directory_path)

    def _populated_consolidated_seq_to_sample_and_abund_dict(self):
        """Go through the list_of_pre_med_sample_dirs. There will be one per sample.
        Get the list of sequences and their abundances using the fasta and name file pairs.
        Get the sample from the fasta name. Get the current dictionary represented by the sequence key
        and add the sample in question-abunance k, v pairing to it. Then move on to next sample."""
        count = 0
        num_samples_to_process = len(self.list_of_pre_med_sample_dirs)
        print('Populating the consolidated sequence to sample and abundance dictionary'
              'for pre-MED sequence processing')
        for sample_pm_dir in self.list_of_pre_med_sample_dirs:
            count += 1
            print(f'Processing pre-MED seqs for sample {count} of {num_samples_to_process}')
            current_dss_obj = None
            # get list of the fasta files (one per clade) that we will need to process
            # we can deduce the .names file from the fasta file simply by changing the extension
            sample_list_of_fasta_file_paths = [
                f_path for f_path in self.thread_safe_general.return_list_of_file_paths_in_directory(sample_pm_dir) if
                '.fasta' in f_path]
            for f_path in sample_list_of_fasta_file_paths:
                fasta_dict = self.thread_safe_general.create_dict_from_fasta(fasta_path=f_path)
                names_dict = self.thread_safe_general.create_seq_name_to_abundance_dict_from_name_file(
                    name_file_path=f_path.replace('.fasta', '.names'))
                clade = f_path.split('/')[-1].split('_')[3]
                # we only need to get the datasetsample object if we haven't already
                if current_dss_obj is None:
                    sample_name = '_'.join(f_path.split('/')[-1].split('_')[4:]).replace('.fasta', '')
                    current_dss_obj = DataSetSample.objects.get(
                        data_submission_from=self.dataset_object, name=sample_name)
                # Check to see whether a dictionary entry already exists for this clade
                # and if not create a default dict to use
                if clade not in self.consolidated_sequence_to_sample_and_abund_dict.keys():
                    self.consolidated_sequence_to_sample_and_abund_dict[clade] = defaultdict(dict)
                # for each sequence of the fasta log it in the consolidated dictionary
                clade_dict_to_add_to = self.consolidated_sequence_to_sample_and_abund_dict[clade]
                for seq_name, seq_seq in fasta_dict.items():
                    clade_dict_to_add_to[seq_seq][current_dss_obj] = names_dict[seq_name]

    def make_data_set_sample_pm_objects(self):
        print('\nProcessing pre-MED seqs for each clade')
        for clade, seq_dict in self.consolidated_sequence_to_sample_and_abund_dict.items():
            print(f'\nProcessing clade {clade}')
            seq_matcher = self.SeqMatcher(
                clade=clade, seq_dict=seq_dict, rs_dict=self.ref_seq_sequence_to_ref_seq_obj_dict[clade],
                match_dict=self.ref_seq_match_obj_to_seq_sample_abundance_dict[clade],
                non_match_dict=self.no_match_consolidated_seq_to_sample_and_abund_dict[clade],
                num_proc=self.num_proc, path_to_seq_match_executable=self.path_to_seq_match_executable,
                temp_working_directory=self.temp_working_directory
            )
            seq_matcher.match_and_make_ref_seqs()

    class SeqMatcher:
        def __init__(
                self, clade, rs_dict, seq_dict, match_dict, non_match_dict,
                num_proc, path_to_seq_match_executable, temp_working_directory):
            # The current clade we are working with
            self.clade = clade
            # Dict of nucleotide sequence to ref seq obj for all ref seq objs of this clade
            self.rs_dict = rs_dict
            # dict of sequences as keys and dictionaries as value where dict
            # is DataSetSample object as key and the absolute abundance of the sequence as value
            self.seq_dict = seq_dict
            # This dict will be ref seq object to the DataSetSample abundance info from self.seq_dict
            self.match_dict = match_dict
            # This dict will be the same structure as self.seq_dict but for the no matches
            self.non_match_dict = non_match_dict
            # The consolidation path that we will follow to consolidate the non refseq match sequences
            self.consolidation_path_list = []
            self.thread_safe_general = ThreadSafeGeneral()
            self.num_proc = num_proc
            self.path_to_seq_match_executable = path_to_seq_match_executable
            self.temp_working_directory = temp_working_directory

        def match_and_make_ref_seqs(self):
            self._assign_sequence_to_match_or_non_match_dicts_mp()
            # Assess whether this has helped us out of the bottle neck or not
            if self.non_match_dict:
                self._consolidate_non_match_seqs()
                self._make_new_reference_sequences_and_populate_match_dict()
            self._create_data_set_sample_sequence_pm_objects()

        @staticmethod
        def _make_seq_match_dicts(
                pre_med_seq_chunk_path, rs_dict_path, path_to_seq_match_executable, match_dict_out_path):
            print(f'Starting match seq_match with thread {get_ident()}')
            subprocess.run(
                [path_to_seq_match_executable, pre_med_seq_chunk_path, rs_dict_path, match_dict_out_path], check=True)
            print(f'seq_match execution complete for thread {get_ident()}')

        def _assign_sequence_to_match_or_non_match_dicts_mp(self):
            matching_start_time = time.time()
            print('Attempting to match sequences to ReferenceSequence objects')
            # A generator that returns self.seq_dict, split up into chunks
            seqs_to_match = len(self.seq_dict)
            # We add the one so that we always have slightly more sequences per chunk, so that we don't use
            # more than self.num_proc number of sequence
            print(f'Chunking in to {int(1 + (len(self.seq_dict)/self.num_proc))} sequence chunks for {self.num_proc} threads')
            pre_med_seq_chunks = self.thread_safe_general.chunks(self.seq_dict.keys(), n=int(1 + (len(self.seq_dict)/self.num_proc)))
            all_processes = []
            # A digit to make path names unique
            count = 0
            output_dict_paths = []
            print('Multithreading match finding. This make take some time...')
            # Rather than dumping the self.rs_dict that contains ReferenceSequence objects
            # we will work with a list of the ReferenceSequence nucleotide sequences and we will match to these
            # This way we don't have to import Django settings and models for a second time and this will hopefully
            # help us avoid the problems with the multiprocessing and the Django testing framework
            # This will mean that we need to do a further look up to connect the matched sequence back to a
            # ReferenceSequence, but this should still be much faster than the serial way we were doing it before
            # as we can split up the many to many searches across the threads.
            rs_set_to_dump = list(self.rs_dict.keys())
            rs_set_path = os.path.join(self.temp_working_directory, f'rs_set_seqs')
            print(f'Writing out {rs_set_path}')
            with open(rs_set_path, 'w') as f:
                json.dump(rs_set_to_dump, f)
            
            print('Done')
            for pre_med_seq_chunk in pre_med_seq_chunks:
                # Dump out the dictionaries that will be read in by the executable
                # Make the names unique to prevent any need for Lock.
                pre_med_seq_chunk_path = os.path.join(self.temp_working_directory, f'pre_med_seq_chunk_{count}')

                match_dict_out_path = os.path.join(self.temp_working_directory, f'pre_med_match_dict_{count}')
                output_dict_paths.append(match_dict_out_path)
                print(f'Writing out {pre_med_seq_chunk_path}')
                with open(pre_med_seq_chunk_path, 'w') as f:
                    json.dump(pre_med_seq_chunk, f)
                
                print('Done')


                p = Thread(target=self._make_seq_match_dicts,
                            args=(
                                pre_med_seq_chunk_path, rs_set_path,
                                self.path_to_seq_match_executable, match_dict_out_path)
                            )

                all_processes.append(p)
                print(f'Starting thread {count}')
                p.start()
                count += 1

            print('All threads started.\nCollecting threads.')
            for p in all_processes:
                p.join()
                print('Collected a thread')

            print('Processing thread outputs.')
            match_dict = {}
            non_match_list = []
            for output_dict_path in output_dict_paths:
                with open(output_dict_path, 'r') as f:
                    sub_match_dict = json.load(f)
                with open(output_dict_path.replace('match_dict', 'non_match_list'), 'r') as f:
                    sub_non_match_list = json.load(f)

                match_dict.update(sub_match_dict)
                non_match_list.extend(sub_non_match_list)
            assert (seqs_to_match == (len(match_dict) + len(non_match_list)))
            # Here we now know exatly which pre_med seqs had a match and which did not
            # For each of the seqs that had a match, we need to now log the match

            print("Logging matches\n")
            match_count = 0
            tot_matches = len(match_dict)
            for nuc_seq, matching_ref_seq in match_dict.items():
                sys.stdout.write(f'\rprocessing {match_count} out of {tot_matches} matches.')
                self._log_match(nuc_seq, self.rs_dict[matching_ref_seq])
                match_count += 1

            # For those that did not have a match add them to the non_match_dict
            print("\nLogging non_matches\n")
            non_match_count = 0
            tot_non_matches = len(non_match_list)
            for nuc_seq in non_match_list:
                sys.stdout.write(f'\rprocessing {non_match_count} out of {tot_non_matches} matches.')
                self.non_match_dict[nuc_seq] = self.seq_dict[nuc_seq]
                non_match_count += 1

            matching_finish_time = time.time() - matching_start_time
            logging.info(f'pre-MED to ReferenceSequence matching took {matching_finish_time}s to complete '
                         f'using multithreading for clade {self.clade}')

            # Now clean up the files
            for output_dict_path in output_dict_paths:
                os.remove(output_dict_path)
                os.remove(output_dict_path.replace('match_dict', 'non_match_list'))

        def _log_match(self, nuc_seq, rs_obj):
            # Check to see if the rs_obj is already representing in the match
            # dict, and if so combine the value dictionaries
            try:
                # We need to be careful here as it could be that the same DataSetSample
                # object could have had both sequenecs that we are dealing with here.
                # In this case we will need to add together the abundance for that sample
                current_match_dict = self.match_dict[rs_obj]
                seq_dict_to_add = self.seq_dict[nuc_seq]
                new_combined_dict = dict()
                for dss_obj, abundance in seq_dict_to_add.items():
                    try:
                        # If the DataSetSample object is in both dicts, then combine the abundances
                        # and a sincle k, v pair of the DataSetSample object and new abunance
                        # to the match dict
                        new_abund = current_match_dict[dss_obj] + abundance
                        new_combined_dict[dss_obj] = new_abund
                    except KeyError:
                        # If the DataSetSample objects is not in both dicts then simply
                        # add the current k, v pair
                        new_combined_dict[dss_obj] = abundance
                # finally we will need to add the k,v pairs in the current_match_dict
                self.match_dict[rs_obj] = {
                    **new_combined_dict,
                    **{k: v for k, v in current_match_dict.items() if k not in seq_dict_to_add}
                }
            except KeyError:
                # If the rs_obj is not already representing then we can simply
                # add the seq_dict info as the value to the rs_obj key in the match dict
                self.match_dict[rs_obj] = self.seq_dict[nuc_seq]

        def _consolidate_non_match_seqs(self):
            """Here we are going to make what I am calling a consolidation path.
            This path will be a list of tuples in a particular order so that we can consolidate
            the non_match sequences into via super and sub set matches. I.e. consolidating the
            DataSetSample and abundance information of shorter sequences
            with the corresponding information for those sequences that are supersets of the short sequences.
            We will only be doing this for the non_match dictionary as the matched dictionary has already taken
            subsets and supersets into account.

            To do this we will go in order of the shortest sequences first (n) and for each of these
            sequences search through all of the sequences that are larger than them for a superset match.
            We will also search for a match of the sequence + A.
            A short sequence may match multiple longer sequences. We will chose to match the longer sequence
            that is associated with the greaterst number of DataSetSample objects.
            To associate a match we will create a tuple of the sequence and its match.
            When we have created all of these matches we will then be able to follow this path of tuples
            to do the consolidation.

            This should be multiprocessable. (We can give copies of the dictionaries and assign different n
            values to each of the threads. We can then put together the list of tuples for each n in order of asscending
            n values to create the consolidation path.)
            However, I don't think it is necessary to implement this. It already seems to be quite fast."""

            self._make_consolidation_path()

            # At this point we have the consolidation_path_list populated
            # We can now follow this path and create a new dictionary which has the representative sequences
            self._consolidate_non_match_seqs_using_consolidation_path()

        def _make_consolidation_path(self):
            # First get a list of the sequences to work with sorted by order of length
            seq_list = sorted(list(self.non_match_dict.keys()), key=len)
            # len of the shortest element
            start_n = len(seq_list[0])
            # len of the longest element
            finish_n = len(seq_list[-1])
            # Go from smallest n to largest n-1
            print('\nMaking consolidation path for non-ReferenceSequence matching sequences')
            for n in range(start_n, finish_n):
                small_query_seqs = [seq for seq in seq_list if len(seq) == n]
                super_seqs = [seq for seq in seq_list if len(seq) > n]
                count = 0
                tot = len(small_query_seqs)
                for q_seq in small_query_seqs:
                    count += 1
                    sys.stdout.write(f'\rseq {count} out of {tot} for level n={n} of {finish_n}')
                    matches = []
                    for super_seq in super_seqs:
                        if (q_seq in super_seq) or ('A' + q_seq in super_seq):
                            matches.append(super_seq)
                    if matches:
                        if len(matches) > 1:
                            # Here we create a list of tuples where the match sequence and the number of
                            # DataSetSamples that contained that sample.
                            # We then sort it according to the number of DataSetSamples
                            # so that we can get the sequences that was found in the maximum number of DataSetSamples
                            # Finally, we use this sequences as the consolidation representative for the
                            # consolidation path
                            representative_seq = sorted(
                                [(match_seq, len(self.non_match_dict[match_seq].keys())) for match_seq in matches],
                                key=lambda x: x[1], reverse=True)[0][0]
                            self.consolidation_path_list.append((q_seq, representative_seq))
                        else:
                            self.consolidation_path_list.append((q_seq, matches[0]))
                    else:
                        # If there are no matches then there is no entry required in the consolidation path
                        pass

        def _consolidate_non_match_seqs_using_consolidation_path(self):
            for small_seq, super_seq in self.consolidation_path_list:
                small_seq_dict = self.non_match_dict[small_seq]
                super_seq_dict = self.non_match_dict[super_seq]
                new_combined_dict = dict()
                for dss_obj, abundance in small_seq_dict.items():
                    try:
                        # If the DataSetSample object is in both dicts, then combine the abundances
                        # and a single k, v pair of the DataSetSample object and new abunance
                        # to the match dict
                        new_abund = super_seq_dict[dss_obj] + abundance
                        new_combined_dict[dss_obj] = new_abund
                    except KeyError:
                        # If the DataSetSample objects is not in both dicts then simply
                        # add the current k, v pair
                        new_combined_dict[dss_obj] = abundance

                # finally update the k,v dict for the super_seq and delete the entry for the
                self.non_match_dict[super_seq] = {
                    **new_combined_dict,
                    **{k: v for k, v in super_seq_dict.items() if k not in small_seq_dict}
                }
                del self.non_match_dict[small_seq]

        def _make_new_reference_sequences_and_populate_match_dict(self, testing=False):
            """Here we are going to make reference sequences for the non_match representative sequences.
            Essentially we want to end up with the same dictionary that we have for the match sequences,
            i.e. a reference sequence representing the DataSetSample and abundance information so that we can
            eventually get to the business of creating the DataSetSampleSequencePM objects.

            We will make implement a sanity check at this point that can be removed once we have completed testing
            This will check to make sure that none of the sequences of the reference sequences we are about to make
            match or fit into any of the existing reference sequences. We will also check that none of the sequences
            match of fit into any of the other sequences. This will be very slow but worth the price considering
            that we don't want to pollute the referenceSequence pool of the database.
            """

            if testing:
                # first sanity check to see that non of the consolidated sequences fit into any of the other
                # consolidated sequences
                for seq_one, seq_two in itertools.combinations(self.non_match_dict.keys(), 2):
                    if (seq_one in seq_two) or (seq_two in seq_one) or \
                            ('A' + seq_one in seq_two) or ('A' + seq_two in seq_one):
                        raise RuntimeError('Consolidated sequences can be further consolidated')

            # For each consolidated sequences, create a ReferenceSequence after checking to see that the sequence
            # does not already fit into one of the
            print('\ncreating ReferenceSequence objects')
            new_rs_list = []
            for c_seq in self.non_match_dict.keys():
                if testing:
                    if (c_seq in self.rs_dict) or ('A' + c_seq in self.rs_dict):
                        raise RuntimeError(
                            'Consolidated sequence is already found in the ReferenceSequence object collection')
                    for rs_seq in self.rs_dict.keys():
                        if (c_seq in rs_seq) or (rs_seq in c_seq):
                            raise RuntimeError(
                                'Consolidated sequence is already found in the ReferenceSequence object collection')
                # Create the new reference sequence.
                new_rs_list.append(ReferenceSequence(clade=self.clade, sequence=c_seq))

            print(f'\ncreating {len(new_rs_list)} new ReferenceSequence objects in bulk for clade {self.clade}')
            for rs_chunk in self.thread_safe_general.chunks(new_rs_list):
                ReferenceSequence.objects.bulk_create(rs_chunk)

            # Now get the newly create ref seq objects back and create a dict form them
            # with rs sequence as key and the rs object itself as the value
            new_rs_seq_to_obj_dict = {
                rs.sequence: rs for rs in ReferenceSequence.objects.filter(
                    sequence__in=list(self.non_match_dict.keys()))}
            # Now go back through the no match dict and use this dictionary to poulate the match dictionary
            for c_seq in self.non_match_dict.keys():
                self.match_dict[new_rs_seq_to_obj_dict[c_seq]] = self.non_match_dict[c_seq]

        def _create_data_set_sample_sequence_pm_objects(self):
            """Finally now that we have a reference sqeuence object representing
            each of the initial sequences that were found in the DataSetSample objects
            we can create the DataSetSamplePM objects."""
            data_set_sample_sequence_pre_med_list = []
            for rs_rep_obj, dss_abund_dict in self.match_dict.items():
                for dss_obj, abundance in dss_abund_dict.items():
                    dsspm = DataSetSampleSequencePM(reference_sequence_of=rs_rep_obj,
                                                    abundance=abundance,
                                                    data_set_sample_from=dss_obj)
                    data_set_sample_sequence_pre_med_list.append(dsspm)
            print(f'\ncreating {len(data_set_sample_sequence_pre_med_list)} '
                  f'new DataSetSampleSequencePM objects in bulk for clade {self.clade}')
            for dssspm_chunk in self.thread_safe_general.chunks(data_set_sample_sequence_pre_med_list):
                DataSetSampleSequencePM.objects.bulk_create(dssspm_chunk)


class DSSAttributeAssignmentHolder:
    """
    This will hold the values of the attributes of the DataSetSamples that were being set during the
    mothur QC process. The attributes are not able to be set during QC when working within a
    Django TransactionTestCase testing frame work. To get around this problem we will collect the attributes
    and their corresponding values in this class and then set them when we are back outside of the multiprocessing.
    """
    def __init__(self, name, uid):
        self.uid = uid
        self.name = name
        self.num_contigs = 0
        self.post_qc_absolute_num_seqs = 0
        self.post_qc_unique_num_seqs = 0
        self.absolute_num_sym_seqs = 0
        self.unique_num_sym_seqs = 0
        self.non_sym_absolute_num_seqs = 0
        self.non_sym_unique_num_seqs = 0
        self.size_violation_absolute = 0
        self.size_violation_unique = 0
        self.error_in_processing = False
        self.error_reason = 'noError'
        self.cladal_seq_totals = None
        self.initial_processing_complete = False


class InitialMothurHandler:
    def __init__(self, data_loading_parent):
        self.parent = data_loading_parent
        if self.parent.multiprocess:
            self.input_queue_containing_pairs_of_fastq_file_paths = mp_Queue()
            self.worker_manager = Manager()
            self.samples_that_caused_errors_in_qc_mp_list = self.worker_manager.list()
            self.output_queue_for_attribute_data = mp_Queue()
            self._populate_input_queue()
        else:
            self.input_queue_containing_pairs_of_fastq_file_paths = mt_Queue()
            self.samples_that_caused_errors_in_qc_mp_list = []
            self.output_queue_for_attribute_data = mt_Queue()
            self._populate_input_queue()

    def _populate_input_queue(self):
        for fastq_path_pair in self.parent.sample_fastq_pairs:
            sample_name = fastq_path_pair.split('\t')[0].replace('[dS]', '-')
            data_set_sample = DataSetSample.objects.get(
                    name=sample_name, data_submission_from=self.parent.dataset_object
                )
            dss_att_holder = DSSAttributeAssignmentHolder(name=data_set_sample.name, uid=data_set_sample.id)
            self.input_queue_containing_pairs_of_fastq_file_paths.put((fastq_path_pair, dss_att_holder))
            
        for n in range(self.parent.num_proc):
            self.input_queue_containing_pairs_of_fastq_file_paths.put('STOP')

    def execute_worker_initial_mothur(self):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        sys.stdout.write('\nPerforming initial mothur QC\n')
        for n in range(self.parent.num_proc):
            if self.parent.multiprocess:
                p = Process(target=self._worker_initial_mothur,
                            args=(self.input_queue_containing_pairs_of_fastq_file_paths, 
                            self.samples_that_caused_errors_in_qc_mp_list, 
                            self.output_queue_for_attribute_data, 
                            self.parent.temp_working_directory, 
                            self.parent.debug)
                            )
            else:
                p = Thread(target=self._worker_initial_mothur,
                        args=(self.input_queue_containing_pairs_of_fastq_file_paths, 
                        self.samples_that_caused_errors_in_qc_mp_list, 
                        self.output_queue_for_attribute_data, 
                        self.parent.temp_working_directory, 
                        self.parent.debug)
                        )

            all_processes.append(p)
            p.start()

        # Here process the queue as we go because it was getting too large and causing the
        # p.join() to hang. We will put a 'DONE' in queue everytime a process finishes with the
        # list of samples

        self._update_dss_obj_attributes()

        for p in all_processes:
            p.join()

    def _update_dss_obj_attributes(self):
        done_count = 0
        while done_count < self.parent.num_proc:
            dss_proxy = self.output_queue_for_attribute_data.get()
            if dss_proxy == 'DONE':
                done_count += 1
            else:
                # Very large samples are causing an error here when we come to call the .save() method
                # Exception has occurred: InterfaceError
                # connection already closed
                try:
                    dss_obj = DataSetSample.objects.get(id=dss_proxy.uid)
                except InterfaceError:
                    db.connections.close_all()
                    dss_obj = DataSetSample.objects.get(id=dss_proxy.uid)
                except:
                    db.connections.close_all()
                    dss_obj = DataSetSample.objects.get(id=dss_proxy.uid)
                # dss_obj = dss_obj_uid_to_obj_dict[dss_proxy.uid] # TODO delete.
                if dss_proxy.error_in_processing:
                    dss_obj.error_in_processing = dss_proxy.error_in_processing
                    dss_obj.error_reason = dss_proxy.error_reason
                    dss_obj.unique_num_sym_seqs = dss_proxy.unique_num_sym_seqs
                    dss_obj.absolute_num_sym_seqs = dss_proxy.absolute_num_sym_seqs
                    dss_obj.post_qc_absolute_num_seqs = dss_proxy.post_qc_absolute_num_seqs
                    dss_obj.post_qc_unique_num_seqs = dss_proxy.post_qc_unique_num_seqs
                    dss_obj.num_contigs = dss_proxy.num_contigs
                    dss_obj.save()
                else:
                    dss_obj.post_qc_absolute_num_seqs = dss_proxy.post_qc_absolute_num_seqs
                    dss_obj.post_qc_unique_num_seqs = dss_proxy.post_qc_unique_num_seqs
                    dss_obj.num_contigs = dss_proxy.num_contigs
                    dss_obj.save()

    # We will attempt to fix the weakref pickling issue we are having by maing this a static method.
    @staticmethod
    def _worker_initial_mothur(in_q_paths, out_list_error_samples, out_q_attr_data, temp_working_directory, debug):
        """
        This worker performs the pre-MED processing that is primarily mothur-based.
        This QC includes making contigs, screening for ambigous calls (0 allowed), screening for a min 30bp overlap,
        discarding singletons and doublets, in silico PCR. It also checks whether sequences are rev compliment.
        This is all done through the use of an InitialMothurWorker class which in turn makes use of the MothurAnalysis
        class that does the heavy lifting of running the mothur commands in sequence.
        """
        for contigpair, dss_att_holder in iter(in_q_paths.get, 'STOP'):

            initial_morthur_worker = InitialMothurWorker(
                dss_att_holder=dss_att_holder, 
                contig_pair=contigpair, 
                temp_working_directory=temp_working_directory,
                debug=debug, out_q_attr_data=out_q_attr_data
                )

            try:
                initial_morthur_worker.start_initial_mothur_worker()
                out_q_attr_data.put(initial_morthur_worker.dss_att_holder)
            except RuntimeError as e:
                out_list_error_samples.append(e.args[0]['sample_name'])
        out_q_attr_data.put('DONE')
        return


class InitialMothurWorker:
    def __init__(self, dss_att_holder, contig_pair, temp_working_directory, debug, out_q_attr_data):
        self.sample_name = dss_att_holder.name
        self.dss_att_holder = dss_att_holder
        self.cwd = os.path.join(temp_working_directory, self.sample_name)
        self.debug = debug
        os.makedirs(self.cwd, exist_ok=True)
        self.mothur_analysis_object = MothurAnalysis(
            name=self.sample_name,
            output_dir=self.cwd, input_dir=self.cwd,
            fastq_gz_fwd_path=contig_pair.split('\t')[1], fastq_gz_rev_path=contig_pair.split('\t')[2],
            stdout_and_sterr_to_pipe=(not self.debug)
            )
        self.output_queue_for_attribute_data = out_q_attr_data
        self.thread_safe_general = ThreadSafeGeneral()

    def start_initial_mothur_worker(self):
        sys.stdout.write(f'{self.sample_name}: QC started\n')

        self._do_make_contigs()

        # This first screen is for min overlap of 30bp
        self._do_screen_seqs()

        self._do_unique_seqs()

        self._do_fwd_and_rev_pcr()

        self._do_unique_seqs()

        # This second screen is for ambigous sequences (0 allowed)
        self._do_screen_seqs()

        self._do_split_abund()

        self._set_unique_and_abs_num_seqs_after_initial_qc()

        self._write_out_final_name_and_fasta_for_tax_screening()

        self.output_queue_for_attribute_data.put(self.dss_att_holder)

        sys.stdout.write(f'{self.sample_name}: Initial mothur complete\n')

    def _do_make_contigs(self):
        try:
            self.dss_att_holder.num_contigs = self.mothur_analysis_object.execute_make_contigs()
            sys.stdout.write(
                f'{self.sample_name}: data_set_sample_instance_in_q.num_contigs = {self.dss_att_holder.num_contigs}\n')

        except RuntimeError as e:
            if str(e) == 'bad fastq, mothur stuck in loop':
                self.log_qc_error_and_continue(errorreason='Bad fastq, mothur stuck in loop')
                raise RuntimeError({'sample_name': self.sample_name})
            if str(e) == 'error in make.contigs':
                self.log_qc_error_and_continue(errorreason='Error in make.contigs')
                raise RuntimeError({'sample_name': self.sample_name})
            if str(e) == 'Unable to extract fasta path from make.contigs output':
                self.log_qc_error_and_continue(errorreason='Unable to extract fasta path from make.contigs output')
                raise RuntimeError({'sample_name': self.sample_name})
            if str(e) == 'empty fasta':
                self.log_qc_error_and_continue(errorreason='Error in make.contigs')
                raise RuntimeError({'sample_name': self.sample_name})

    def _write_out_final_name_and_fasta_for_tax_screening(self):
        name_file_as_list = self.thread_safe_general.read_defined_file_to_list(
            self.mothur_analysis_object.name_file_path)
        taxonomic_screening_name_file_path = os.path.join(self.cwd, 'name_file_for_tax_screening.names')
        self.thread_safe_general.write_list_to_destination(taxonomic_screening_name_file_path, name_file_as_list)
        fasta_file_as_list = self.thread_safe_general.read_defined_file_to_list(self.mothur_analysis_object.fasta_path)
        taxonomic_screening_fasta_file_path = os.path.join(self.cwd, 'fasta_file_for_tax_screening.fasta')
        self.thread_safe_general.write_list_to_destination(taxonomic_screening_fasta_file_path, fasta_file_as_list)

    def _do_fwd_and_rev_pcr(self):
        try:
            self.mothur_analysis_object.execute_pcr(do_reverse_pcr_as_well=True)
        except RuntimeError as e:
            if str(e) == 'PCR fasta file is blank':
                self.log_qc_error_and_continue(errorreason='No seqs left after PCR')
                raise RuntimeError({'sample_name': self.sample_name})
            self.check_for_error_and_raise_runtime_error(stage_of_qc='PCR')

    def _do_split_abund(self):
        try:
            self.mothur_analysis_object.execute_split_abund()
        except:
            self.check_for_error_and_raise_runtime_error(stage_of_qc='split by abundance')

    def _do_unique_seqs(self):
        try:
            self.mothur_analysis_object.execute_unique_seqs()
        except:
            self.check_for_error_and_raise_runtime_error(stage_of_qc='unique seqs')

    def _do_screen_seqs(self):
        """
        We do two rounds of screen seq. We do the first round of screen seq as a protection again libraries that
        have been made with seq chemistry that is too short to create overlap between the fwd and rev reads.
        We have implemented a penalised mismatch and gapopen score in the make.contigs (both -4). This forces
        mothur to do a better job of identifying overlapping regions. We then screen by minoverlap=30. We do not
        screen by mismatches.
        The second round of screening removes sequences with ambiguous nucleotides.
        This works as a general quality control.
        We do this after PCR to prevent removal of sequenecs that only have ambigous calls outside
        of the primers. Whether this is even possible is probably variable depending of the sequencing technology
        but we keep the check in place as standard.
        """
        try:
            self.mothur_analysis_object.execute_screen_seqs()
            self.mothur_analysis_object.screening_for = 'ambig'
        except RuntimeError as e:
            if self.mothur_analysis_object.screening_for == 'overlap':
                self.check_for_error_and_raise_runtime_error(stage_of_qc='screen seqs overlap', error_summary=str(e))
            elif self.mothur_analysis_object.screening_for == 'ambig':
                self.check_for_error_and_raise_runtime_error(stage_of_qc='screen seqs ambig', error_summary=str(e))

    def _set_unique_and_abs_num_seqs_after_initial_qc(self):
        name_file = self.thread_safe_general.read_defined_file_to_list(self.mothur_analysis_object.name_file_path)

        number_of_contig_seqs_unique = len(name_file)
        self.dss_att_holder.post_qc_unique_num_seqs = number_of_contig_seqs_unique
        sys.stdout.write(
            f'{self.sample_name}: '
            f'data_set_sample_instance_in_q.post_qc_unique_num_seqs = {number_of_contig_seqs_unique}\n')

        abs_count = 0
        for line in name_file:
            abs_count += len(line.split('\t')[1].split(','))

        self.dss_att_holder.post_qc_absolute_num_seqs = abs_count

        sys.stdout.write(
            f'{self.sample_name}: data_set_sample_instance_in_q.post_qc_absolute_num_seqs = {abs_count}\n')

    def check_for_error_and_raise_runtime_error(self, stage_of_qc, error_summary=None):
        if error_summary:
            if error_summary == "no seqs left after overlap seq screening":
                self.log_qc_error_and_continue(
                errorreason=f'No seqs remaining after screen.seqs for minoverlap'
                )
                raise RuntimeError({'sample_name': self.sample_name})
        for stdout_line in self.thread_safe_general.decode_utf8_binary_to_list(
                self.mothur_analysis_object.latest_completed_process_command.stdout
        ):
            if '[WARNING]: Blank fasta name, ignoring read.' in stdout_line:
                self.log_qc_error_and_continue(errorreason=f'Blank fasta name during {stage_of_qc}')
                raise RuntimeError({'sample_name': self.sample_name})
            if 'do not match' in stdout_line:
                self.log_qc_error_and_continue(errorreason=f'error in fastq file during {stage_of_qc}')
                raise RuntimeError({'sample_name': self.sample_name})
            if 'ERROR' in stdout_line:
                self.log_qc_error_and_continue(errorreason=f'error in inital QC during {stage_of_qc}')
                raise RuntimeError({'sample_name': self.sample_name})
        self.log_qc_error_and_continue(errorreason=f'error in inital QC during {stage_of_qc}')
        raise RuntimeError({'sample_name': self.sample_name})

    def log_qc_error_and_continue(self, errorreason):
        print('{}: Error in processing sample'.format(self.sample_name))
        self.dss_att_holder.unique_num_sym_seqs = 0
        self.dss_att_holder.absolute_num_sym_seqs = 0
        self.dss_att_holder.initial_processing_complete = True
        self.dss_att_holder.error_in_processing = True
        self.dss_att_holder.error_reason = errorreason
        self.output_queue_for_attribute_data.put(self.dss_att_holder)


class PotentialSymTaxScreeningHandler:
    """ The purpose of this handler and the executed work is only to get a collection of sequences that will need
    screening against the NCBI database. We also rely on this method to do the blast of our each samples sequences
    against our symclade.fa reference sequences database. We read in the blast.out file in the later functions.
    """
    def __init__(
            self, samples_that_caused_errors_in_qc_list,
            checked_samples_list, list_of_samples_names, num_proc, multiprocess):
        self.multiprocess = multiprocess
        if self.multiprocess:
            self.input_queue = mp_Queue()
            self.manager = Manager()
            self.sub_evalue_sequence_to_num_sampes_found_in_mp_dict = self.manager.dict()
            self.sub_evalue_nucleotide_sequence_to_clade_mp_dict = self.manager.dict()
            self.error_samples_mp_list = self.manager.list(samples_that_caused_errors_in_qc_list)
            self.checked_samples_mp_list = self.manager.list(checked_samples_list)
            self.list_of_sample_names = list_of_samples_names
            self.num_proc = num_proc
            self._load_input_queue()
            self.lock = mp_Lock()
        else:
            self.input_queue = mt_Queue()
            self.sub_evalue_sequence_to_num_sampes_found_in_mp_dict = {}
            self.sub_evalue_nucleotide_sequence_to_clade_mp_dict = {}
            self.error_samples_mp_list = samples_that_caused_errors_in_qc_list
            self.checked_samples_mp_list = checked_samples_list
            self.list_of_sample_names = list_of_samples_names
            self.num_proc = num_proc
            self._load_input_queue()
            self.lock = mt_Lock()

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
        sys.stdout.write('\nPerforming potential sym tax screening QC\n')
        for n in range(self.num_proc):
            if self.multiprocess:
                p = Process(
                    target=self._potential_sym_tax_screening_worker,
                    args=(
                        self.input_queue, 
                        self.error_samples_mp_list,
                        self.checked_samples_mp_list,
                        self.sub_evalue_sequence_to_num_sampes_found_in_mp_dict, 
                        self.sub_evalue_nucleotide_sequence_to_clade_mp_dict,
                        data_loading_temp_working_directory,
                        data_loading_path_to_symclade_db,
                        data_loading_debug, self.lock
                        ))
            else:
                p = Thread(
                target=self._potential_sym_tax_screening_worker,
                args=(
                    self.input_queue, 
                    self.error_samples_mp_list,
                    self.checked_samples_mp_list,
                    self.sub_evalue_sequence_to_num_sampes_found_in_mp_dict, 
                    self.sub_evalue_nucleotide_sequence_to_clade_mp_dict,
                    data_loading_temp_working_directory,
                    data_loading_path_to_symclade_db,
                    data_loading_debug, self.lock
                    ))

            all_processes.append(p)
            p.start()
        for p in all_processes:
            p.join()

    @staticmethod
    def _potential_sym_tax_screening_worker(
            in_q, 
            error_samples_mp_list, 
            checked_samples_mp_list, 
            sub_evalue_sequence_to_num_sampes_found_in_mp_dict, 
            sub_evalue_nucleotide_sequence_to_clade_mp_dict, 
            data_loading_temp_working_directory, 
            data_loading_path_to_symclade_db, 
            data_loading_debug, lock):
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
        for sample_name in iter(in_q.get, 'STOP'):

            # If the sample gave an error during the inital mothur then we don't consider it here.
            if sample_name in error_samples_mp_list:
                continue

            # A sample will be in this list if we have already performed this worker on it and none of its sequences
            # gave matches to the symClade database at below the evalue threshold
            if sample_name in checked_samples_mp_list:
                continue

            taxonomic_screening_worker = PotentialSymTaxScreeningWorker(
                sample_name=sample_name, wkd=data_loading_temp_working_directory,
                path_to_symclade_db=data_loading_path_to_symclade_db, debug=data_loading_debug,
                checked_samples_mp_list=checked_samples_mp_list,
                e_val_collection_mp_dict=sub_evalue_sequence_to_num_sampes_found_in_mp_dict,
                sub_evalue_nucleotide_sequence_to_clade_mp_dict=sub_evalue_nucleotide_sequence_to_clade_mp_dict,
                lock=lock)

            taxonomic_screening_worker.execute_tax_screening()


class PotentialSymTaxScreeningWorker:
    def __init__(
            self, sample_name, wkd, path_to_symclade_db, debug, e_val_collection_mp_dict,
            checked_samples_mp_list, sub_evalue_nucleotide_sequence_to_clade_mp_dict, lock):
        self.thread_safe_general = ThreadSafeGeneral()
        self.sample_name = sample_name
        self.cwd = os.path.join(wkd, self.sample_name)
        self.fasta_file_path = os.path.join(self.cwd, 'fasta_file_for_tax_screening.fasta')
        self.fasta_dict = self.thread_safe_general.create_dict_from_fasta(fasta_path=self.fasta_file_path)
        self.name_file_path = os.path.join(self.cwd, 'name_file_for_tax_screening.names')
        self.name_dict = {
            a.split('\t')[0]: a for a in self.thread_safe_general.read_defined_file_to_list(self.name_file_path)}
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
        # The potential_non_symbiodiniaceae_sequences_list is used to see if there are any samples, that don't have
        # any potential non symbiodnium sequences. i.e. only definite symbiodiniaceae sequences.
        # These samples are added to
        # the checked list and are not checked in following iterations to speed things up.
        self.potential_non_symbiodiniaceae_sequences_list = []
        self.sequence_name_to_clade_dict = None
        self.blast_output_as_list = None
        self.already_processed_blast_seq_result = []
        # this is a managed list that holds the names of samples from which no sequences were thrown out from
        # it will be used in downstream processes.
        self.checked_samples_mp_list = checked_samples_mp_list
        self.lock = lock

    def execute_tax_screening(self):
        sys.stdout.write(f'{self.sample_name}: verifying seqs are Symbiodinium and determining clade\n')

        blastn_analysis = BlastnAnalysis(
            input_file_path=self.fasta_file_path,
            output_file_path=os.path.join(self.cwd, 'blast.out'), db_path=self.path_to_symclade_db,
            output_format_string="6 qseqid sseqid staxids evalue pident qcovs")

        if self.debug:
            blastn_analysis.execute_blastn_analysis(pipe_stdout_sterr=False)
        else:
            blastn_analysis.execute_blastn_analysis(pipe_stdout_sterr=True)

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

        if not self.potential_non_symbiodiniaceae_sequences_list:
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
            # here we are looking for sequences to add to the non_symbiodiniaceae_sequence_list
            # if a sequence fails at any of our if statements it will be added to the non_symbiodiniaceae_sequence_list
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
        The potential_non_symbiodiniaceae_sequences_list is used to see if there are any samples, that don't have
        any potential non symbiodnium sequences. i.e. only definite symbiodiniaceae sequences.
        These samples are added to
        the checked list and are not checked in following iterations to speed things up.
        """

        if identity < 80 or coverage < 95:
            # incorporate the size cutoff here that would normally happen in the further mothur qc later in the code
            self.potential_non_symbiodiniaceae_sequences_list.append(name_of_current_sequence)
            if 184 < len(self.fasta_dict[name_of_current_sequence]) < 310:
                with self.lock:
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
        self.potential_non_symbiodiniaceae_sequences_list.extend(list(sequences_with_no_blast_match_as_set))
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
            data_loading_num_proc, multiprocess, data_loading_dataset_object):
        self.multiprocess = multiprocess
        self.ds_object = data_loading_dataset_object
        if self.multiprocess:
            self.sample_name_mp_input_queue = mp_Queue()
            self.sym_non_sym_mp_manager = Manager()
            self.samples_that_caused_errors_in_qc_mp_list = self.sym_non_sym_mp_manager.list(
                data_loading_samples_that_caused_errors_in_qc_mp_list
            )
            self.sample_attributes_mp_output_queue = mp_Queue()
            self.non_symbiodiniaceae_sequences_list = self.sym_non_sym_mp_manager.list()
            self.num_proc = data_loading_num_proc
            self._populate_input_queue(data_loading_list_of_samples_names)
        else:
            self.sample_name_mp_input_queue = mt_Queue()
            self.samples_that_caused_errors_in_qc_mp_list =  data_loading_samples_that_caused_errors_in_qc_mp_list
            self.sample_attributes_mp_output_queue = mt_Queue()
            self.non_symbiodiniaceae_sequences_list = []
            self.num_proc = data_loading_num_proc
            self._populate_input_queue(data_loading_list_of_samples_names)

    def _populate_input_queue(self, data_loading_list_of_samples_names):
        for sample_name in data_loading_list_of_samples_names:
            self.sample_name_mp_input_queue.put(
                DataSetSample.objects.get(name=sample_name, data_submission_from=self.ds_object))
        for n in range(self.num_proc):
            self.sample_name_mp_input_queue.put('STOP')

    def execute_sym_non_sym_tax_screening(
            self, data_loading_temp_working_directory,
            non_symb_and_size_violation_base_dir_path, data_loading_pre_med_sequence_output_directory_path,
            data_loading_debug):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        sys.stdout.write('\nPerforming sym non sym tax screening QC\n')

        for n in range(self.num_proc):
            if self.multiprocess:
                p = Process(target=self._sym_non_sym_tax_screening_worker, args=(
                    self.sample_name_mp_input_queue, 
                    self.samples_that_caused_errors_in_qc_mp_list, 
                    self.sample_attributes_mp_output_queue,
                    data_loading_temp_working_directory,
                    non_symb_and_size_violation_base_dir_path, 
                    data_loading_pre_med_sequence_output_directory_path,
                    data_loading_debug))
            else:
                p = Thread(target=self._sym_non_sym_tax_screening_worker, args=(
                self.sample_name_mp_input_queue, 
                self.samples_that_caused_errors_in_qc_mp_list, 
                self.sample_attributes_mp_output_queue,
                data_loading_temp_working_directory, 
                non_symb_and_size_violation_base_dir_path,
                data_loading_pre_med_sequence_output_directory_path,
                data_loading_debug))
            all_processes.append(p)
            p.start()

        self._associate_info_to_dss_objects()

        for p in all_processes:
            p.join()

    @staticmethod
    def _sym_non_sym_tax_screening_worker(
        in_q, 
        samples_that_caused_errors_in_qc_mp_list, 
        sample_attributes_mp_output_queue, 
        data_loading_temp_working_directory, 
        data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path,
        data_loading_pre_med_sequence_output_directory_path,
        data_loading_debug):

        for dss in iter(in_q.get, 'STOP'):

            if dss.name in samples_that_caused_errors_in_qc_mp_list:
                continue

            sym_non_sym_tax_screening_worker_object = SymNonSymTaxScreeningWorker(
                data_loading_temp_working_directory=data_loading_temp_working_directory,
                dss=dss,
                data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path=
                data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path,
                data_loading_pre_med_sequence_output_directory_path=data_loading_pre_med_sequence_output_directory_path,
                data_loading_debug=data_loading_debug,
                sample_attributes_mp_output_queue=sample_attributes_mp_output_queue
            )

            try:
                sym_non_sym_tax_screening_worker_object.identify_sym_non_sym_seqs()
            except RuntimeError as e:
                samples_that_caused_errors_in_qc_mp_list.append(e.args[0]['sample_name'])
        sample_attributes_mp_output_queue.put('DONE')

    def _associate_info_to_dss_objects(self):
        # now save the collected data contained in the sample_attributes_holder_mp_dict to the relevant
        # dss_objs.
        done_count = 0

        dss_uid_to_dss_obj_dict = {dss.id: dss for dss in
                                   DataSetSample.objects.filter(data_submission_from=self.ds_object)}
        while done_count < self.num_proc:
            dss = self.sample_attributes_mp_output_queue.get()
            if dss == 'DONE':
                done_count += 1
            else:
                dss.save()


class SymNonSymTaxScreeningWorker:
    def __init__(
            self, data_loading_temp_working_directory, dss,
            data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path,
            data_loading_pre_med_sequence_output_directory_path, data_loading_debug,
            sample_attributes_mp_output_queue
    ):
        # init core objects
        self.dss = dss
        self.cwd = os.path.join(data_loading_temp_working_directory, self.dss.name)
        self.debug = data_loading_debug
        self.sample_attributes_mp_output_queue = sample_attributes_mp_output_queue
        self.thread_safe_general = ThreadSafeGeneral()

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
        self.non_symbiodiniaceae_sequence_name_set_for_sample = set()
        self.sym_size_violation_sequence_name_set_for_sample = set()
        self.sym_no_size_violation_sequence_name_set_for_sample = set()

    def _init_sym_non_size_violation_output_paths(self, data_loading_pre_med_sequence_output_directory_path):
        pre_med_sequence_output_directory_for_sample_path = os.path.join(
            data_loading_pre_med_sequence_output_directory_path, f'{self.dss.id}_{self.dss.name}')
        os.makedirs(pre_med_sequence_output_directory_for_sample_path, exist_ok=True)
        self.pre_med_fasta_path = os.path.join(
            pre_med_sequence_output_directory_for_sample_path, f'pre_med_seqs_{self.dss.name}.fasta'
        )
        self.pre_med_names_file_path = os.path.join(
            pre_med_sequence_output_directory_for_sample_path, f'pre_med_seqs_{self.dss.name}.names'
        )

    def _init_non_sym_and_size_violation_output_paths(
            self, data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path):
        non_symbiodiniaceae_and_size_violation_directory_for_sample_path = os.path.join(
            data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path, self.dss.name)
        os.makedirs(non_symbiodiniaceae_and_size_violation_directory_for_sample_path, exist_ok=True)
        self.non_symbiodiniaceae_seqs_fasta_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.dss.name}_non_symbiodiniaceae_sequences.fasta')
        self.non_symbiodiniaceae_seqs_names_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.dss.name}_non_symbiodiniaceae_sequences.names')
        self.symbiodiniaceae_size_violation_seqs_fasta_output_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.dss.name}_symbiodiniaceae_size_violation_sequences.fasta'
        )
        self.symbiodiniaceae_size_violation_seqs_names_output_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.dss.name}_symbiodiniaceae_size_violation_sequences.names')

    def _init_fasta_name_blast_and_clade_dict_attributes(self):
        fasta_file_path = os.path.join(self.cwd, 'fasta_file_for_tax_screening.fasta')
        self.fasta_dict = self.thread_safe_general.create_dict_from_fasta(fasta_path=fasta_file_path)
        name_file_path = os.path.join(self.cwd, 'name_file_for_tax_screening.names')
        self.name_dict = {
            a.split('\t')[0]: a for a in self.thread_safe_general.read_defined_file_to_list(name_file_path)}
        blast_output_path = os.path.join(self.cwd, 'blast.out')
        self.blast_dict = {blast_line.split('\t')[0]: blast_line for blast_line in
                           self.thread_safe_general.read_defined_file_to_list(blast_output_path)}
        self.sequence_name_to_clade_dict = {
            blast_out_line.split('\t')[0]: blast_out_line.split('\t')[1][-1] for
            blast_out_line in self.blast_dict.values()
        }


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

        self._associate_qc_meta_info_to_dss_objs()

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
                f'seqs_for_med_{self.dss.name}_clade_{clade_of_sequences_to_write_out}.redundant.fasta'
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
                        f.write(f'>{self.dss.name}_{sequence_counter}\n')
                        f.write(f'{self.fasta_dict[sequence_name]}\n')
                        sequence_counter += 1

            self._if_debug_check_length_of_deuniqued_fasta(sample_clade_fasta_path)

    def _if_debug_check_length_of_deuniqued_fasta(self, sample_clade_fasta_path):
        if self.debug:
            deuniqued_fasta = self.thread_safe_general.read_defined_file_to_list(sample_clade_fasta_path)
            if deuniqued_fasta:
                if len(deuniqued_fasta) < 100:
                    print(f'{self.dss.name}: WARNING the dequniqed fasta is < {len(deuniqued_fasta)} lines')
            else:
                print(f'{self.dss.name}: ERROR deuniqued fasta is empty')

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
        for clade_value in set(self.sequence_name_to_clade_dict.values()):
            pre_med_names_path_clade_specific = self.pre_med_names_file_path.replace(
                f'pre_med_seqs_{self.dss.name}', f'pre_med_seqs_{clade_value}_{self.dss.name}')

            with open(pre_med_names_path_clade_specific, 'w') as f:
                for sequence_name in [
                        seq_name for seq_name, clade_val in self.sequence_name_to_clade_dict.items() if
                        clade_val == clade_value and seq_name in
                        self.sym_no_size_violation_sequence_name_set_for_sample]:
                    f.write(f'{self.name_dict[sequence_name]}\n')
                    self.absolute_number_of_sym_no_size_violation_sequences += len(
                        self.name_dict[sequence_name].split('\t')[1].split(','))

    def _write_out_no_size_violation_fasta_to_pre_med_dir(self):
        """ Write out the pre-med seqs into clade separated .fasta and .name pairs"""
        for clade_value in set(self.sequence_name_to_clade_dict.values()):
            pre_med_fasta_path_clade_specific = self.pre_med_fasta_path.replace(
                f'pre_med_seqs_{self.dss.name}', f'pre_med_seqs_{clade_value}_{self.dss.name}')

            with open(pre_med_fasta_path_clade_specific, 'w') as f:
                for sequence_name in [
                    seq_name for seq_name, clade_val in self.sequence_name_to_clade_dict.items() if
                        clade_val == clade_value and seq_name in
                        self.sym_no_size_violation_sequence_name_set_for_sample]:

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
            if query_sequence_name not in self.non_symbiodiniaceae_sequence_name_set_for_sample:
                if 184 < len(self.fasta_dict[query_sequence_name]) < 310:
                    self.sym_no_size_violation_sequence_name_set_for_sample.add(query_sequence_name)
                else:
                    self.sym_size_violation_sequence_name_set_for_sample.add(query_sequence_name)

    def _associate_qc_meta_info_to_dss_objs(self):
        self._associate_non_sym_seq_attributes_to_datasetsample()
        self._associate_sym_seq_no_size_violation_attributes_to_datasetsample()
        self._associate_sym_seq_size_violation_attributes_to_datasetsample()
        self.dss.initial_processing_complete = True
        self.sample_attributes_mp_output_queue.put(self.dss)
        print(f'{self.dss.name}: pre-med QC complete')

    def _associate_sym_seq_size_violation_attributes_to_datasetsample(self):
        # it is important that we set these values from the objects of this class rather
        # from the data_set_sample object attributes as some of these have not been set yet
        self.dss.size_violation_absolute = (
                self.dss.post_qc_absolute_num_seqs -
                self.dss.absolute_num_sym_seqs -
                self.dss.non_sym_absolute_num_seqs
        )
        print(f'{self.dss.name}: size_violation_absolute = {self.dss.size_violation_absolute}')
        self.dss.size_violation_unique = (
                self.dss.post_qc_unique_num_seqs -
                self.dss.unique_num_sym_seqs -
                self.dss.non_sym_unique_num_seqs
        )
        print(f'{self.dss.name}: size_violation_unique = {self.dss.size_violation_unique}')

    def _associate_sym_seq_no_size_violation_attributes_to_datasetsample(self):
        self.dss.unique_num_sym_seqs = len(self.sym_no_size_violation_sequence_name_set_for_sample)
        self.dss.absolute_num_sym_seqs = self.absolute_number_of_sym_no_size_violation_sequences
        print(f'{self.dss}: unique_num_sym_seqs = {len(self.sym_no_size_violation_sequence_name_set_for_sample)}')
        print(f'{self.dss}: absolute_num_sym_seqs = {self.absolute_number_of_sym_no_size_violation_sequences}')

    def _associate_non_sym_seq_attributes_to_datasetsample(self):
        self.dss.non_sym_unique_num_seqs = len(self.non_symbiodiniaceae_sequence_name_set_for_sample)
        self.dss.non_sym_absolute_num_seqs = self.absolute_number_of_non_sym_sequences
        print(
            f'{self.dss.name}: non_sym_unique_num_seqs = {len(self.non_symbiodiniaceae_sequence_name_set_for_sample)}')
        print(f'{self.dss.name}: non_sym_absolute_num_seqs = {self.absolute_number_of_non_sym_sequences}')

    def _log_dataset_attr_and_raise_runtime_error(self):
        # if there are no symbiodiniaceae sequenes then log error and associate meta info
        print(f'{self.dss.name}: QC error.\n No symbiodiniaceae sequences left in sample after pre-med QC.')
        self.dss.non_sym_unique_num_seqs = len(self.non_symbiodiniaceae_sequence_name_set_for_sample)
        self.dss.non_sym_absolute_num_seqs = self.absolute_number_of_non_sym_sequences
        self.dss.size_violation_absolute = self.absolute_number_of_sym_size_violation_sequences
        self.dss.size_violation_unique = len(self.sym_size_violation_sequence_name_set_for_sample)
        self.dss.unique_num_sym_seqs = 0
        self.dss.absolute_num_sym_seqs = 0
        self.dss.initial_processing_complete = True
        self.dss.error_in_processing = True
        self.dss.error_reason = 'No symbiodiniaceae sequences left in sample after pre-med QC'
        self.sample_attributes_mp_output_queue.put(self.dss)
        raise RuntimeError({'sample_name': self.dss.name})

    def _identify_and_write_non_sym_seqs_in_sample(self):
        self._add_seqs_with_no_blast_match_to_non_sym_list()
        self._identify_non_sym_seqs_from_below_match_threshold()
        self._write_out_non_sym_fasta_and_names_files_for_sample()

    def _write_out_non_sym_fasta_and_names_files_for_sample(self):
        if self.non_symbiodiniaceae_sequence_name_set_for_sample:
            self._write_out_non_sym_fasta_for_sample()
            self._write_out_non_sym_names_file_for_sample()

    def _write_out_non_sym_names_file_for_sample(self):
        with open(self.non_symbiodiniaceae_seqs_names_path, 'w') as f:
            for sequence_name in list(self.non_symbiodiniaceae_sequence_name_set_for_sample):
                f.write(f'{self.name_dict[sequence_name]}\n')
                self.absolute_number_of_non_sym_sequences += len(
                    self.name_dict[sequence_name].split('\t')[1].split(','))

    def _write_out_non_sym_fasta_for_sample(self):
        with open(self.non_symbiodiniaceae_seqs_fasta_path, 'w') as f:
            for sequence_name in list(self.non_symbiodiniaceae_sequence_name_set_for_sample):
                f.write(f'>{sequence_name}\n')
                f.write(f'{self.fasta_dict[sequence_name]}\n')

    def _add_seqs_with_no_blast_match_to_non_sym_list(self):
        sequences_with_no_blast_match_as_set = set(self.fasta_dict.keys()) - \
                                               set(self.sequence_name_to_clade_dict.keys())
        self.non_symbiodiniaceae_sequence_name_set_for_sample.update(list(sequences_with_no_blast_match_as_set))
        sys.stdout.write(
            f'{self.dss.name}: {len(sequences_with_no_blast_match_as_set)} sequences thrown out '
            f'initially due to being too divergent from reference sequences\n')

    def _identify_non_sym_seqs_from_below_match_threshold(self):
        """ This method is where the non_symbiodiniaceae sequences for the sample are identified. If they fall below
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
            self.non_symbiodiniaceae_sequence_name_set_for_sample.add(name_of_current_sequence)


class PerformMEDHandler:
    def __init__(self, data_loading_temp_working_directory, data_loading_num_proc, multiprocess):
        # need to get list of the directories in which to perform the MED
        # we want to get a list of the .
        self.multiprocess = multiprocess
        self.temp_working_directory = data_loading_temp_working_directory
        self.num_proc = data_loading_num_proc
        self.list_of_redundant_fasta_paths = []
        self._populate_list_of_redundant_fasta_paths()
        if self.multiprocess:
            self.input_queue_of_redundant_fasta_paths = mp_Queue()
        else:
            self.input_queue_of_redundant_fasta_paths = mt_Queue()
        self._populate_input_queue_of_redundant_fasta_paths()
        self.list_of_med_result_dirs = [
            os.path.join(os.path.dirname(path_to_redundant_fasta), 'MEDOUT') for
            path_to_redundant_fasta in self.list_of_redundant_fasta_paths]
        
    def execute_perform_med_worker(
            self, data_loading_debug, data_loading_path_to_med_padding_executable,
            data_loading_path_to_med_decompose_executable):
        all_processes = []

        for n in range(self.num_proc):
            if self.multiprocess:
                p = Process(target=self._perform_med_worker, args=(
                    self.input_queue_of_redundant_fasta_paths, data_loading_debug,
                    data_loading_path_to_med_padding_executable,
                    data_loading_path_to_med_decompose_executable))
            else:
                p = Thread(target=self._perform_med_worker, args=(
                self.input_queue_of_redundant_fasta_paths, data_loading_debug,
                data_loading_path_to_med_padding_executable,
                data_loading_path_to_med_decompose_executable))
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

    @staticmethod
    def _perform_med_worker(
            in_q, data_loading_debug, data_loading_path_to_med_padding_executable,
            data_loading_path_to_med_decompose_executable):
        for redundant_fata_path in iter(in_q.get, 'STOP'):
            perform_med_worker_instance = PerformMEDWorker(
                redundant_fasta_path=redundant_fata_path, data_loading_debug=data_loading_debug,
                data_loading_path_to_med_padding_executable=data_loading_path_to_med_padding_executable,
                data_loading_path_to_med_decompose_executable=data_loading_path_to_med_decompose_executable)

            perform_med_worker_instance.do_decomposition()


class PerformMEDWorker:
    def __init__(
            self, redundant_fasta_path, data_loading_path_to_med_padding_executable, data_loading_debug,
            data_loading_path_to_med_decompose_executable):
        self.thread_safe_general = ThreadSafeGeneral()
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
        # We cannot add check=True to the subprocess calls as we are expecting some of them to fail
        # when there are too few sequences
        if not self.debug:
            subprocess.run(
                [self.path_to_med_decompose_executable, '-M', str(self.med_m_value), '--skip-gexf-files',
                 '--skip-gen-figures',
                 '--skip-gen-html', '--skip-check-input', '-T', '-o',
                 self.med_output_dir, self.redundant_fasta_path_padded], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elif self.debug:
            subprocess.run(
                [self.path_to_med_decompose_executable, '-M', str(self.med_m_value), '--skip-gexf-files',
                 '--skip-gen-figures',
                 '--skip-gen-html',
                 '--skip-check-input', '-T', '-o',
                 self.med_output_dir, self.redundant_fasta_path_padded])

    def _get_med_m_value(self):
        # Define MED M value dynamically.
        # The M value is a cutoff that looks at the abundance of the most abundant unique sequence in a node
        # if the abundance is lower than M then the node is discarded
        # we have been working recently with an M that equivaltes to 0.4% of 0.004. This was
        # calculated when working with a modelling project where I was subsampling to 1000 sequences. In this
        # scenario the M was set to 4.
        # We should also take care that M doesn't go below 4, so we should use a max choice for the M
        num_of_seqs_to_decompose = len(
            self.thread_safe_general.read_defined_file_to_list(self.redundant_fasta_path_unpadded)) / 2
        return max(4, int(0.004 * num_of_seqs_to_decompose))


class DataSetSampleSequenceCreatorWorker:
    """This class will be responsible for handling a set of med outputs. Objects will be things like the directory,
    the count table, number of samples, number of nodes, these sorts of things."""
    def __init__(self, med_output_directory,
                 data_set_sample_creator_handler_ref_seq_sequence_to_ref_seq_id_dict,
                 data_set_sample_creator_handler_ref_seq_uid_to_ref_seq_name_dict, data_loading_dataset_obj):
        self.thread_safe_general = ThreadSafeGeneral()
        self.output_directory = med_output_directory
        self.sample_name = self.output_directory.split('/')[-3]
        self.clade = self.output_directory.split('/')[-2]
        self.nodes_list_of_nucleotide_sequences = []
        self._populate_nodes_list_of_nucleotide_sequences()
        self.num_med_nodes = len(self.nodes_list_of_nucleotide_sequences)
        self.node_sequence_name_to_ref_seq_id = {}
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
            node_file_as_list = self.thread_safe_general.read_defined_file_to_list(node_file_path)
        except FileNotFoundError:
            raise RuntimeError({'med_output_directory': self.output_directory})

        for i in range(0, len(node_file_as_list), 2):
            node_seq_name = node_file_as_list[i].split('|')[0][1:]
            node_seq_abundance = int(node_file_as_list[i].split('|')[1].split(':')[1])
            node_seq_sequence = node_file_as_list[i+1].replace('-', '')
            self.nodes_list_of_nucleotide_sequences.append(
                NucleotideSequence(name=node_seq_name, abundance=node_seq_abundance, sequence=node_seq_sequence))

    def make_data_set_sample_sequences(self):
        self._associate_med_nodes_to_ref_seq_objs()

        if self._two_or_more_nodes_associated_to_the_same_reference_sequence():
            self._make_associations_and_abundances_in_node_abund_df_unique_again()

        self._update_data_set_sample_med_qc_meta_data_and_clade_totals()

        self._if_more_than_200_seqs_create_clade_collection()

        self._create_data_set_sample_sequences()

    def _associate_med_nodes_to_ref_seq_objs(self):
        for node_nucleotide_sequence_object in self.nodes_list_of_nucleotide_sequences:
            if not self._assign_node_sequence_to_existing_ref_seq(node_nucleotide_sequence_object):
                self._assign_node_sequence_to_new_ref_seq(node_nucleotide_sequence_object)

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
        for dsss_chunk in self.thread_safe_general.chunks(data_set_sample_sequence_list):
            DataSetSampleSequence.objects.bulk_create(dsss_chunk)

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
        for dsss_chunk in self.thread_safe_general.chunks(data_set_sample_sequence_list):
            DataSetSampleSequence.objects.bulk_create(dsss_chunk)
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

            self._update_node_name_abund_in_df_and_nuc_seq_list(node_names_to_be_consolidated,
                                                                summed_abund_of_nodes_of_ref_seq)

    def _get_list_of_node_names_that_need_consolidating(self, non_unique_ref_seq_uid):
        node_names_to_be_consolidated = [
            node_name for node_name in self.node_sequence_name_to_ref_seq_id.keys() if
            self.node_sequence_name_to_ref_seq_id[node_name] == non_unique_ref_seq_uid]
        return node_names_to_be_consolidated

    def _update_node_name_abund_in_df_and_nuc_seq_list(self, node_names_to_be_consolidated,
                                                       summed_abund_of_nodes_of_ref_seq):
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
        This takes into account whether the seq_in_q could be a subset or super set of one of the
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
        self.ref_seq_sequence_to_ref_seq_id_dict[new_ref_seq.sequence] = new_ref_seq.id
        self.node_sequence_name_to_ref_seq_id[node_nucleotide_sequence_object.name] = new_ref_seq.id
        self.ref_seq_uid_to_ref_seq_name_dict[new_ref_seq.id] = str(new_ref_seq)

        sys.stdout.write(f'\r{self.sample_name} clade {self.clade}: '
                         f'Assigning MED node {node_nucleotide_sequence_object.name} '
                         f'to new reference sequence {str(new_ref_seq)}')


class DataSetSampleCreatorHandler:
    """This class will be where we run the code for creating reference sequences, data set sample sequences and
    clade collections."""
    def __init__(self):
        # dictionaries to save us having to do lots of database look ups
        self.ref_seq_uid_to_ref_seq_name_dict = {
            ref_seq.id: str(ref_seq) for ref_seq in ReferenceSequence.objects.all()}
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
                        f'WARNING node file contains only '
                        f'{data_set_sample_sequence_creator_worker.num_med_nodes} sequences.')
            sys.stdout.write(
                f'\n\nPopulating {data_set_sample_sequence_creator_worker.sample_name} with '
                f'clade {data_set_sample_sequence_creator_worker.clade} sequences\n')
            data_set_sample_sequence_creator_worker.make_data_set_sample_sequences()