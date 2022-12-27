#!/usr/bin/env python3
""" SymPortal: a novel analytical framework and platform for coral algal
    symbiont next-generation sequencing ITS2 profiling
    Copyright (C) 2018  Benjamin C C Hume

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see
    https://github.com/didillysquat/SymPortal_framework/tree/master/LICENSE.txt.


    """
__version__ = '0.3.22'


# Django specific settings
import os
from datetime import datetime
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.conf import settings
# ####### Setup Django DB and Models ########
# Ensure settings are read
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
# Your application specific imports
from dbApp.models import DataSet, DataAnalysis, DataSetSample, Study, User, Citation
############################################


import output
import plotting
import sys
import distance
import argparse
import data_loading
import sp_config
if sp_config.system_type == 'remote':
    import textdistance
    from bs4 import BeautifulSoup
    import requests
    import re
import data_analysis
from general import ThreadSafeGeneral
from django_general import CreateStudyAndAssociateUsers
import django_general
from shutil import which
import time
import subprocess
import json
from django.core.exceptions import ObjectDoesNotExist
import logging

class SymPortalWorkFlowManager:
    def __init__(self, custom_args_list=None):
        self.start_time = time.time()
        self.args = self._define_args(custom_args_list)
        # general attributes
        self.thread_safe_general = ThreadSafeGeneral()
        self.symportal_root_directory = os.path.abspath(os.path.dirname(__file__))
        self.dbbackup_dir = os.path.join(self.symportal_root_directory, 'dbBackUp')
        os.makedirs(self.dbbackup_dir, exist_ok=True)
        self.date_time_str = str(datetime.utcnow()).split('.')[0].replace('-','').replace(' ','T').replace(':','')
        self.submitting_user = sp_config.user_name
        self.submitting_user_email = sp_config.user_email
        self.number_of_samples = None
        # for data_loading
        self.data_loading_object = None
        self.data_set_object = None
        self.screen_sub_eval_bool = None
        if sp_config.system_type == 'remote':
            self.screen_sub_eval_bool = True
            self.study = None
        else:
            self.screen_sub_eval_bool = False
        self.reference_db = 'symClade.fa'
        self.output_seq_count_table_obj = None
        # for data analysis
        self.within_clade_cutoff = 0.03
        self.data_analysis_object = None
        self.sp_data_analysis = None
        self.output_type_count_table_obj = None
        self.type_stacked_bar_plotter = None
        # If the shortcut function analyse_next has been used
        # look up the uids of the last analysis and use this plus what has been
        # provided to redefine self.args.analyse
        if self.args.analyse_next:
            self._redefine_arg_analyse()
        # these will be used in all but the data loading. in the dataloading an output dir and html dir
        # are created as part of the dataloading object.
        self.output_dir = None
        self.html_dir = None
        self.js_output_path_dict = {}
        self.js_file_path = None
        # Variables that will hold the distance class objects
        self.unifrac_distance_object = None
        self.braycurtis_distance_object = None

    def _redefine_arg_analyse(self):
        """
        When the user passes the argument --analyse_next then we will find the UIDs
        that were used for the previous analysis and append the passed UIDs to them
        to create a new string. 
        In doing this we will check that all of the dataset IDs exist as it may be that some
        datasets have been deleted since the last analysis
        We will then change self.args.analyse to this value
        and work with that.
        """
        last_analysis = sorted(DataAnalysis.objects.all(), key=lambda x: x.id, reverse=True)[0]
        last_uids = [int(_) for _ in last_analysis.list_of_data_set_uids.split(',')]
        # Then we want to filter these ids for only those that exist
        last_uids = [ds.id for ds in DataSet.objects.filter(id__in=last_uids)]
        new_uids = [int(_) for _ in self.args.analyse_next.split(',')]
        if len(set(last_uids).intersection(set(new_uids))) != 0:
            raise RuntimeError("There appears to be overlap in the uids being provided to the --analyse_next argument "
                               "and the uids provided to the previous "
                               f"analysis:\n\tnew = {new_uids}\n\tlast = {last_uids}")
        last_uids.extend(new_uids)
        new_uid_css = ','.join([str(_) for _ in last_uids])
        self.args.analyse = new_uid_css

    def _define_args(self, custom_args_list=None):
        parser = argparse.ArgumentParser(
            description='Intragenomic analysis of the ITS2 region of the nrDNA',
            epilog='For support email: benjamin.hume@uni-konstanz.de')
        group = parser.add_mutually_exclusive_group(required=True)
        self._define_mutually_exclusive_args(group)
        self._define_additional_args(parser)
        if custom_args_list is not None:
            return parser.parse_args(custom_args_list)
        else:
            return parser.parse_args()

    @staticmethod
    def _define_additional_args(parser):
        parser.add_argument('--num_proc', type=int, help='Number of processors to use', default=1)
        parser.add_argument('--name', help='A name for your input or analysis', default='noName')
        parser.add_argument('--description', help='An optional description', default='No description')
        parser.add_argument('--data_analysis_id', type=int, help='The ID of the data_analysis you wish to output from')
        parser.add_argument('--bootstrap', type=int, help='Number of bootstrap iterations to perform', default=100)
        parser.add_argument(
            '--data_sheet',
            help='An absolute path to the .xlsx file containing the meta-data information for the data_set\'s samples')
        parser.add_argument('--no_figures', action='store_true', help='Skip figure production')
        parser.add_argument('--no_ordinations', action='store_true', help='Skip ordination analysis')
        parser.add_argument('--debug', action='store_true', help='Present additional stdout output', default=False)

        parser.add_argument(
            '--no_output', action='store_true', help='Do no output: count tables, figures, ordinations', default=False)

        parser.add_argument(
            '--distance_method',
            help='Either \'unifrac\', \'braycurtis\', or \'both\' [both]. The method to use when '
                 'calculating distances between its2 type profiles or samples. '
                 '\n\'unifrac\' - ouput only unifrac-derived distance matrices'
                 '\n\'braycurtis\' - output only braycurtis-derived distance matrices'
                 '\n\'all\' - output both unifrac- and braycurtis-derived distance matrices '
                 'mafft and iqtree will be checked for in your PATH. If not found, only braycurtis-derived distances '
                 'will be output', default='both')
        parser.add_argument('--local',
                            help="When passed, only the DataSetSamples of the current output will be used"
                                 " matrices will be calculated using the DIV abundance info from all"
                                 " DataSetSamples the ITS2 type profiles were found in."
                                 " This flag will only have an effect when applied to between ITS2 type profile "
                                 "distances. It will have no effect when calculating between sample distances. "
                                 "[False]",
                            action='store_true', default=False)
        parser.add_argument('--no_pre_med_seqs',
                            help="When passed, DataSetSampleSequencePM objects will not be created"
                                 "[False]", action='store_true', default=False)
        parser.add_argument('--multiprocess', help="When passed, concurrency will be acheived using "
                                                   "multiprocessing rather than multithreading.",
                            action='store_true', default=False)
        parser.add_argument('--force_basal_lineage_separation',
                            help="When passed, cladocopium profiles sequences from the C3, C15 and C1 radiations "
                                 "will not be allowed to occur together in profiles.",
                            action='store_true', default=False)
        # when run as remote
        if sp_config.system_type == 'remote':
            parser.add_argument(
                '--submitting_user_name',
                help='Only for use when running as remote\nallows the association of a different user_name to the '
                     'data_set than the one listed in sp_config', default='not supplied')

            parser.add_argument(
                '--study_user_string',
                help='Only for use when running as remote\nThe comma separated string of the User '
                     'names that should be associated to the Study object.')

            parser.add_argument(
                '--study_name',
                help='Only for use when running as remote\nThe name that will be given to the'
                     'Study associated to the given DataSet object')

            parser.add_argument('--is_cron_loading',
                                help='This is passed only when the loading is being '
                                     'initiated as part of one of the cron jobs.',
                                action='store_true', default=False)

            parser.add_argument(
                '--submitting_user_email',
                help='Only for use when running as remote\nallows the association of a '
                     'different user_email to the data_set '
                     'than the one listed in sp_config', default='not supplied')
    @staticmethod
    def _define_mutually_exclusive_args(group):
        group.add_argument(
            '--load', metavar='path_to_dir',
            help='Run this to load data to the framework\'s database. The first argument to this command must be an '
                 'absolute path to a directory containing  the paired sequencing '
                 'reads in .fastq.gz format. Alternatively, '
                 'this path can point directly to a single compressed file containing the same paired fastq.gz files. '
                 '\nA name must be associated with the data_set using the --name flag. \n'
                 'The number of processes to use '
                 'can also be specified using the --num_proc flag. \nA datasheet can also be uploaded using the '
                 '--data_sheet flag and the full path to the .xlsx data_sheet file (RECOMMENDED). \n'
                 'To skip the generation of figures pass the --no_figures flag.\n To skip the generation of '
                 'ordination files (pairwise distances and PCoA coordinates) pass the --no_ordinations flag')
        group.add_argument(
            '--analyse', metavar='DataSet UIDs',
            help='Analyse one or more data_set objects together. Enter comma separated UIDs of the data_set uids you '
                 'which to analyse. e.g.: 43,44,45. If you wish to use all available dataSubmissions, you may pass '
                 '\'all\' as an argument. To display all data_sets currently submitted to the framework\'s database, '
                 'including their ids, use the \'show_data_sets\' command\nTo skip the generation of figures pass the '
                 '--no_figures flag.\nTo skip the generation of ordination files '
                 '(pairwise distances and PCoA coordinates) '
                 'pass the --no_ordinations flag')
        group.add_argument(
            '--analyse_next', metavar='DataSet UIDs',
            help='This is a convenience function that builds on the --analyse function. Instead of providing a comma '
                 'separated list of the DataSet UIDs to be analysed, the user can provide a comma separated list of '
                 'UIDs that will be analysed IN ADDITION to the DataSets of the last completed analysis. '
                 'In other words, this function saves the user the trouble of having to look up which UIDs made up '
                 'the previous analysis and then adding to this string.')
        group.add_argument(
            '--display_data_sets', action='store_true', help='Display data_sets currently in the framework\'s database')
        group.add_argument(
            '--display_analyses', action='store_true',
            help=' Display data_analysis objects currently stored in the framework\'s database')
        group.add_argument(
            '--print_output_seqs', metavar='DataSet UIDs',
            help='Use this function to output ITS2 sequence count tables for given data_set instances')
        group.add_argument(
            '--print_output_seqs_sample_set', metavar='DataSetSample UIDs',
            help='Use this function to output ITS2 sequence count tables for a collection of DataSetSample instances. '
                 'The input to this function should be a comma separated string of '
                 'the UIDs of the DataSetSample instances '
                 'in question. e.g. 345,346,347,348')
        group.add_argument('--update_citations', help='Check for new articles citing the SymPortal MS.', action='store_true')
        if sp_config.system_type == 'local':
            group.add_argument(
                '--print_output_types', metavar='DataSet UIDs, DataAnalysis UID',
                help='Use this function to output the ITS2 sequence and ITS2 type profile count tables for a given set of '
                     'data_sets that have been run in a given analysis. '
                     'Give the data_set uids that you wish to make outputs '
                     'for as arguments to the --print_output_types flag. To output for multiple data_set objects, '
                     'comma separate the uids of the data_set objects, '
                     'e.g. 44,45,46. Give the ID of the analysis you wish to '
                     'output these from using the --data_analysis_id flag.\nTo skip the generation of figures pass the '
                     '--no_figures flag.')
            group.add_argument(
                '--print_output_types_sample_set', metavar='DataSetSample UIDs, DataAnalysis UID',
                help='Use this function to output the ITS2 sequence and ITS2 type profile count tables for a given set of '
                     'DataSetSample objects that have been run in a given DataAnalysis. Give the DataSetSample '
                     'UIDs that you wish to make outputs from as arguments to the --print_output_types_sample_set flag. '
                     'To output for '
                     'multiple DataSetSample objects, comma separate the UIDs of the DataSetSample objects, '
                     'e.g. 5644,5645,5646. Give the UID of the DataAnalysis you wish to output these from using the '
                     '--data_analysis_id flag.\nTo skip the generation of figures pass the '
                     '--no_figures flag.')
        elif sp_config.system_type == 'remote':
            group.add_argument(
                '--output_study_from_analysis', metavar='Study UID or name',
                help='Use this function to output the ITS2 sequence and ITS2 type profile count tables'
                     ' as well as the associated dissimilarity distances for a given Study object. '
                     'Give the Study UID or name that you wish to output '
                     ' as arguments to the --output_study_from_analysis flag. '
                     'NB. Names should not be numerical. '
                     'Give the ID of the analysis you wish to '
                     'output the Study from using the --data_analysis_id flag.', )
            group.add_argument(
            '--display_studies', action='store_true', help='Display studies currently in the framework\'s database')
        else:
            raise NotImplementedError
        group.add_argument(
            '--between_type_distances', metavar='DataSetSample UIDs, DataAnalysis UID',
            help='Use this function to output UniFrac pairwise distances between ITS2 type profiles clade separated')
        group.add_argument(
            '--between_type_distances_sample_set', metavar='DataSetSample UIDs, DataAnalysis UID',
            help='Use this function to output pairwise distances between ITS2 type profiles clade '
                 'separated from a given collection of DataSetSample objects')
        group.add_argument(
            '--between_type_distances_cct_set', metavar='CladeCollectionType UIDs, DataAnalysis UID',
            help='Use this function to output pairwise distances between a specific set of CladeCollection-AnalysisType'
                 ' associations.')
        group.add_argument(
            '--between_sample_distances', metavar='DataSetSample UIDs',
            help='Use this function to output pairwise distances between samples clade separated from a '
                 'given collection of DataSet objects')
        group.add_argument(
            '--between_sample_distances_sample_set', metavar='DataSetSample UIDs',
            help='Use this function to output pairwise distances between samples clade '
                 'separated from a given collection of DataSetSample objects')
        group.add_argument(
            '--vacuum_database', action='store_true',
            help='Vacuuming the database will free up memory from objects that have been deleted recently')
        group.add_argument('--apply_data_sheet', metavar='DataSet UID',
                           help='Use this function to apply the meta info in a datasheet to '
                                'the DataSetSamples of a given DataSet. Provide the UID of the DataSet to which the '
                                'info should be applied and give the path to the datasheet that contains the '
                                'information using the --data_sheet flag. The sample names in the datasheet must match '
                                'the DataSetSample names exactly. Unpopulated columns in the datasheet will be ignored.'
                                ' I.e. existing meta-information will not be removed from the DataSetSampes if '
                                'information is missing in the datasheet.')
        group.add_argument(
            '--version', '-v', action='store_true',
            help='Output version')

    def start_work_flow(self):
        if self.args.load:
            self.perform_data_loading()
        elif self.args.analyse:
            self._perform_data_analysis()

        # Output data
        elif self.args.print_output_seqs:
            self.perform_stand_alone_sequence_output()
        elif self.args.print_output_seqs_sample_set:
            self.perform_stand_alone_sequence_output()
        elif self._check_for_type_output():
            # Type profile outputs dependent on sp_config.system_type
            return
        # Distances
        elif self.args.between_type_distances:
            self.perform_type_distance_stand_alone()
        elif self.args.between_type_distances_sample_set:
            self.perform_type_distance_stand_alone()
        elif self.args.between_type_distances_cct_set:
            self.perform_type_distance_stand_alone()
        elif self.args.between_sample_distances:
            self._perform_sample_distance_stand_alone()
        elif self.args.between_sample_distances_sample_set:
            self._perform_sample_distance_stand_alone()

        # DB display functions
        elif self._check_for_display_arguments():
            return

        # Apply datasheet
        elif self.args.apply_data_sheet:
            self.apply_datasheet_to_dataset_samples()

        # Finally, if we are on remote check for citation updates
        elif self.args.update_citations:
            if sp_config.system_type == 'remote':
                citation_updater = CitationUpdate()
                try:
                    print('Checking for citation updates')
                    citation_updater.update()
                    print('Complete')
                except Exception as e:
                    print(e)
                    print('ERROR: Citation updating failed.')
            else:
                print('Operation only available on remote system')

    def _check_for_display_arguments(self):
        if self.args.display_data_sets:
            self.perform_display_data_sets()
            return True
        elif self.args.display_analyses:
            self.perform_display_analysis_types()
            return True
        elif self.args.vacuum_database:
            self.perform_vacuum_database()
            return True
        elif self.args.version:
            print(__version__)
            return True
        # Only if we are running as remote, check for study output
        if sp_config.system_type == 'remote':
            if self.args.display_studies:
                self.perform_display_studies()
                return True
        return False

    def _check_for_type_output(self):
        if sp_config.system_type == 'local':
            if self.args.print_output_types:
                self.perform_stand_alone_type_output()
            elif self.args.print_output_types_sample_set:
                self.perform_stand_alone_type_output()
            else:
                return False
        elif sp_config.system_type == 'remote':
            if self.args.output_study_from_analysis:
                # Then we will call one of the below functions after
                # checking that the study exists and whether the Study's
                # collection is set by a number of DataSet or DataSetSample UIDs
                self.output_study_from_analysis()
            else:
                return False
        return True

    # GENERAL
    def _plot_if_not_too_many_samples(self, plotting_function):
        if self.number_of_samples < 1000:
            plotting_function()
        else:
            print(f'Too many samples {self.number_of_samples} to plot.')

    def _set_data_analysis_obj_from_arg_analysis_uid(self):
        self.data_analysis_object = DataAnalysis.objects.get(id=self.args.data_analysis_id)

    def _verify_data_analysis_uid_provided(self):
        if not self.args.data_analysis_id:
            raise RuntimeError(
                'Please provide a data_analysis to ouput from by providing a data_analysis '
                'ID to the --data_analysis_id flag. To see a list of data_analysis objects in the '
                'framework\'s database, use the --display_analyses flag.')

    def _plot_sequence_stacked_bar_from_seq_output_table(self):
        """Plot up the sequence abundances from the output sequence count table. NB this is in the
        case where we have not run an analysis in conjunction, i.e. there are no ITS2 type profiles to consider.
        As such, no ordered list of DataSetSamples should be passed to the plotter."""
        self.seq_stacked_bar_plotter = plotting.SeqStackedBarPlotter(
            output_directory=self.output_seq_count_table_obj.output_dir,
            seq_relative_abund_count_table_path_post_med=self.output_seq_count_table_obj.path_to_seq_output_abund_and_meta_df_absolute,
            no_pre_med_seqs=self.args.no_pre_med_seqs,
            ordered_seq_list=self.output_seq_count_table_obj.clade_abundance_ordered_ref_seq_list,
            date_time_str=self.output_seq_count_table_obj.date_time_str,
            seq_relative_abund_df_pre_med=self.output_seq_count_table_obj.output_df_relative_pre_med)
        self.seq_stacked_bar_plotter.plot_stacked_bar_seqs()

    def _plot_type_stacked_bar_from_type_output_table(self):
        self.type_stacked_bar_plotter = plotting.TypeStackedBarPlotter(
            output_directory=self.output_type_count_table_obj.output_dir,
            type_relative_abund_count_table_path=self.output_type_count_table_obj.path_to_relative_count_table_profiles_abund_and_meta,
            date_time_str=self.output_type_count_table_obj.date_time_str)
        self.type_stacked_bar_plotter.plot_stacked_bar_profiles()

    @staticmethod
    def _plot_type_distances_from_distance_object(distance_object):
        """Search for the paths of the .csv PCoA files and pass these into plotting"""
        sys.stdout.write('\n\nPlotting ITS2 type profile distances\n')
        for pcoa_path in [path for path in distance_object.output_path_list if path.endswith('.csv')]:
            try:
                local_plotter = plotting.DistScatterPlotterTypes(
                    csv_path=pcoa_path, date_time_str=distance_object.date_time_str)
                local_plotter.make_type_dist_scatter_plot()
            except RuntimeError:
                # The error message is printed to stdout at the source
                continue

    @staticmethod
    def _plot_sample_distances_from_distance_object(distance_object):
        """Search for the path of the .csv file that holds the PC coordinates and pass this into plotting"""
        sys.stdout.write('\n\nPlotting sample distances\n')
        for pcoa_path in [path for path in distance_object.output_path_list if path.endswith('.csv')]:
            try:
                local_plotter = plotting.DistScatterPlotterSamples(
                    csv_path=pcoa_path, date_time_str=distance_object.date_time_str)
                local_plotter.make_sample_dist_scatter_plot()
            except RuntimeError:
                # The error message is printed to stdout at the source
                continue

    # DATA ANALYSIS
    def _perform_data_analysis(self):
        self._verify_name_arg_given_analysis()
        self.create_new_data_analysis_obj()
        self.output_dir = os.path.join(
            self.symportal_root_directory, 'outputs', 'analyses', str(self.data_analysis_object.id), self.date_time_str)
        self._set_html_dir_and_js_out_path_from_output_dir()
        self._start_data_analysis()

        if not self.args.no_output:
            os.makedirs(self.html_dir, exist_ok=True)
            self._do_data_analysis_output()
            if not self.args.no_ordinations:
                self._do_data_analysis_ordinations()
            else:
                print('Ordinations skipped at user\'s request')

            # here output the js_output_path item for the DataExplorer
            self._output_js_output_path_dict()
            self._print_analysis_obj_attributes()

        else:
            print('\nOutputs skipped at user\'s request\n')
            self._print_analysis_obj_attributes()

    def _print_analysis_obj_attributes(self):
        try:
            print(f'\n ANALYSIS COMPLETE:\n'
                  f'\tDataAnalysis name: {self.data_analysis_object.name}\n'
                  f'\tDataAnalysis UID: {self.data_analysis_object.id}\n'
                  f'\tStudy name: {self.study.name}\n'
                  f'\tStudy UID: {self.study.name}')
        except AttributeError:
            print(f'\n ANALYSIS COMPLETE: DataAnalysis:\n\tname: '
                  f'{self.data_analysis_object.name}\n\tUID: {self.data_analysis_object.id}\n')
        self.data_analysis_object.analysis_complete_time_stamp = str(
            datetime.utcnow()
        ).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')
        self.data_analysis_object.save()
        print(f'DataSet analysis_complete_time_stamp: '
              f'{self.data_analysis_object.analysis_complete_time_stamp}\n\n\n')

    def _verify_name_arg_given_analysis(self):
        if self.args.name == 'noName':
            sys.exit('Please provide a name using --name')

    def _output_js_output_path_dict(self):
        """Out put the dict that holds the output files so that we can list them in the DataExplorer"""
        # covert the full paths to relative paths and then write out dict
        # https://stackoverflow.com/questions/8693024/how-to-remove-a-path-prefix-in-python
        new_dict = {}
        for k, v in self.js_output_path_dict.items():
            new_dict[k] = os.path.relpath(v, self.output_dir)

        self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
            [{'function_name': 'getDataFilePaths', 'python_obj': new_dict}],
            js_outpath=self.js_file_path)

    def _start_data_analysis(self):
        # Perform the analysis
        self.sp_data_analysis = data_analysis.SPDataAnalysis(
            workflow_manager_parent=self,
            data_analysis_obj=self.data_analysis_object,
            force_basal_lineage_separation=self.args.force_basal_lineage_separation)
        self.sp_data_analysis.analyse_data()

    def _do_data_analysis_output(self):
        self._make_data_analysis_output_type_tables()
        self._make_data_analysis_output_seq_tables()
        self.number_of_samples = len(self.output_type_count_table_obj.sorted_list_of_vdss_uids_to_output)

        if not self.args.no_figures:
            self._plot_if_not_too_many_samples(self._plot_type_stacked_bar_from_type_output_table)

            self._plot_if_not_too_many_samples(self._plot_sequence_stacked_bar_with_ordered_dss_uids_from_type_output)
        else:
            print('\nFigure plotting skipped at user\'s request')

    def _plot_sequence_stacked_bar_with_ordered_dss_uids_from_type_output(self):
        """Plot the sequence abundance info from the output sequence count table ensuring to take in the same
        DataSetSample order as that used in the ITS2 type profile output that was conducted in parallel."""
        self.seq_stacked_bar_plotter = plotting.SeqStackedBarPlotter(
            output_directory=self.output_seq_count_table_obj.output_dir,
            seq_relative_abund_count_table_path_post_med=self.output_seq_count_table_obj.path_to_seq_output_abund_and_meta_df_absolute,
            ordered_sample_uid_list=self.output_type_count_table_obj.sorted_list_of_vdss_uids_to_output,
            no_pre_med_seqs=self.args.no_pre_med_seqs,
            ordered_seq_list=self.output_seq_count_table_obj.clade_abundance_ordered_ref_seq_list,
            date_time_str=self.output_seq_count_table_obj.date_time_str,
            seq_relative_abund_df_pre_med=self.output_seq_count_table_obj.output_df_relative_pre_med)
        self.seq_stacked_bar_plotter.plot_stacked_bar_seqs()

    def _do_data_analysis_ordinations(self):
        self._perform_data_analysis_type_distances()
        if not self.args.no_figures:
            if self.args.distance_method == 'both':
                self._plot_if_not_too_many_samples(
                    lambda: self._plot_type_distances_from_distance_object(self.unifrac_distance_object))
                self._plot_if_not_too_many_samples(
                    lambda: self._plot_type_distances_from_distance_object(self.braycurtis_distance_object))
            elif self.args.distance_method == 'unifrac':
                self._plot_if_not_too_many_samples(
                    lambda: self._plot_type_distances_from_distance_object(self.unifrac_distance_object))
            elif self.args.distance_method == 'braycurtis':
                self._plot_if_not_too_many_samples(
                    lambda: self._plot_type_distances_from_distance_object(self.braycurtis_distance_object))

        self._perform_data_analysis_sample_distances()
        if not self.args.no_figures:
            if self.args.distance_method == 'both':
                self._plot_if_not_too_many_samples(
                    lambda: self._plot_sample_distances_from_distance_object(self.unifrac_distance_object))
                self._plot_if_not_too_many_samples(
                    lambda: self._plot_sample_distances_from_distance_object(self.braycurtis_distance_object))
            elif self.args.distance_method == 'unifrac':
                self._plot_if_not_too_many_samples(
                    lambda: self._plot_sample_distances_from_distance_object(self.unifrac_distance_object))
            elif self.args.distance_method == 'braycurtis':
                self._plot_if_not_too_many_samples(
                    lambda: self._plot_sample_distances_from_distance_object(self.braycurtis_distance_object))

    def _perform_data_analysis_sample_distances(self):
        if self.args.distance_method == 'both':
            if self._check_if_required_packages_found_in_path():
                self._start_analysis_braycurtis_sample_distances()
                self._start_analysis_unifrac_sample_distances()
            else:
                self.args.distance_method = 'braycurtis'
        if self.args.distance_method == 'unifrac':
            self._start_analysis_unifrac_sample_distances()
        elif self.args.distance_method == 'braycurtis':  # braycurtis
            self._start_analysis_braycurtis_sample_distances()

    def _start_analysis_unifrac_sample_distances(self):
        self.unifrac_distance_object = distance.SampleUnifracDistPCoACreator(
            num_processors=self.args.num_proc,
            date_time_str=self.date_time_str,
            data_set_sample_uid_list=self.output_type_count_table_obj.sorted_list_of_vdss_uids_to_output,
            output_dir=self.output_dir,
            html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict)
        self.unifrac_distance_object.compute_unifrac_dists_and_pcoa_coords()

    def _start_analysis_braycurtis_sample_distances(self):
        self.braycurtis_distance_object = distance.SampleBrayCurtisDistPCoACreator(
            date_time_str=self.date_time_str,
            data_set_sample_uid_list=self.output_type_count_table_obj.sorted_list_of_vdss_uids_to_output,
            output_dir=self.output_dir,
            html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict)
        self.braycurtis_distance_object.compute_braycurtis_dists_and_pcoa_coords()

    def _perform_data_analysis_type_distances(self):
        if self.args.distance_method == 'both':
            if self._check_if_required_packages_found_in_path():
                self._start_analysis_braycurtis_type_distances()
                self._start_analysis_unifrac_type_distances()
            else:
                self.args.distance_method = 'braycurtis'
        if self.args.distance_method == 'unifrac':
            self._start_analysis_unifrac_type_distances()
        elif self.args.distance_method == 'braycurtis':  # braycurtis
            self._start_analysis_braycurtis_type_distances()

    def _start_analysis_unifrac_type_distances(self):
        self.unifrac_distance_object = distance.TypeUnifracDistPCoACreator(
            num_processors=self.args.num_proc,
            data_analysis_obj=self.data_analysis_object,
            date_time_str=self.date_time_str,
            data_set_sample_uid_list=self.output_type_count_table_obj.sorted_list_of_vdss_uids_to_output,
            output_dir=self.output_dir,
            local_abunds_only=self.args.local,
            html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict)
        self.unifrac_distance_object.compute_unifrac_dists_and_pcoa_coords()

    def _start_analysis_braycurtis_type_distances(self):
        self.braycurtis_distance_object = distance.TypeBrayCurtisDistPCoACreator(
            data_analysis_obj=self.data_analysis_object,
            date_time_str=self.date_time_str,
            data_set_sample_uid_list=self.output_type_count_table_obj.sorted_list_of_vdss_uids_to_output,
            output_dir=self.output_dir,
            local_abunds_only=self.args.local,
            html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict)
        self.braycurtis_distance_object.compute_braycurtis_dists_and_pcoa_coords()

    def _make_data_analysis_output_seq_tables(self):
        self.output_seq_count_table_obj = output.SequenceCountTableCreator(
            num_proc=self.args.num_proc,
            symportal_root_dir=self.symportal_root_directory,
            no_pre_med_seqs=self.args.no_pre_med_seqs,
            ds_uids_output_str=self.data_analysis_object.list_of_data_set_uids,
            output_dir=self.output_dir,
            sorted_sample_uid_list=self.output_type_count_table_obj.sorted_list_of_vdss_uids_to_output,
            analysis_obj=self.data_analysis_object,
            date_time_str=self.date_time_str,
            html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict, multiprocess=self.args.multiprocess,
            call_type='analysis')
        self.output_seq_count_table_obj.make_seq_output_tables()

    def _make_data_analysis_output_type_tables(self):
        # Write out the AnalysisType count table
        self.output_type_count_table_obj = output.OutputProfileCountTable(
            call_type='analysis', num_proc=self.args.num_proc,
            within_clade_cutoff=self.within_clade_cutoff,
            data_set_uids_to_output=self.sp_data_analysis.list_of_data_set_uids,
            virtual_object_manager=self.sp_data_analysis.virtual_object_manager,
            data_analysis_obj=self.sp_data_analysis.data_analysis_obj,
            output_dir=self.output_dir, html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict,
            date_time_str=self.date_time_str, force_basal_lineage_separation=self.args.force_basal_lineage_separation)
        self.output_type_count_table_obj.output_types()

    def create_new_data_analysis_obj(self):
        self.data_analysis_object = DataAnalysis(
            list_of_data_set_uids=self.args.analyse, within_clade_cutoff=self.within_clade_cutoff,
            name=self.args.name, time_stamp=self.date_time_str,
            submitting_user=self.submitting_user, submitting_user_email=self.submitting_user_email)
        self.data_analysis_object.description = self.args.description
        self.data_analysis_object.save()

    # DATA LOADING
    def perform_data_loading(self):
        self._verify_name_arg_given_load()
        self._execute_data_loading()
        if sp_config.system_type == 'remote' and self.data_loading_object.study and not self.args.no_output:
            self.output_dir = self.data_loading_object.output_directory
            self.study = self.data_loading_object.study
            self._output_study_output_info_items()

    def _execute_data_loading(self):
        if sp_config.system_type == 'remote' and self.args.is_cron_loading:
            self.data_loading_object = data_loading.DataLoading(
                    parent_work_flow_obj=self, datasheet_path=self.args.data_sheet, user_input_path=self.args.load,
                    screen_sub_evalue=self.screen_sub_eval_bool, num_proc=self.args.num_proc, no_fig=self.args.no_figures,
                    no_ord=self.args.no_ordinations, no_output=self.args.no_output,
                    distance_method=self.args.distance_method,
                    no_pre_med_seqs=self.args.no_pre_med_seqs, debug=self.args.debug, multiprocess=self.args.multiprocess,
                    start_time=self.start_time, date_time_str=self.date_time_str,
                    is_cron_loading=True,
                    study_name=self.args.study_name, study_user_string=self.args.study_user_string)
        else:
            self.data_loading_object = data_loading.DataLoading(
                parent_work_flow_obj=self, datasheet_path=self.args.data_sheet, user_input_path=self.args.load,
                screen_sub_evalue=self.screen_sub_eval_bool, num_proc=self.args.num_proc, no_fig=self.args.no_figures,
                no_ord=self.args.no_ordinations, no_output=self.args.no_output,
                distance_method=self.args.distance_method,
                no_pre_med_seqs=self.args.no_pre_med_seqs, debug=self.args.debug, multiprocess=self.args.multiprocess,
                start_time=self.start_time, date_time_str=self.date_time_str,
                is_cron_loading=False)
        
        self.data_loading_object.load_data()

    def _verify_name_arg_given_load(self):
        """If no name arg is provided use the folder name of the self.args.load argument"""
        if self.args.name == 'noName':
            if self.args.load.endswith('/'):
                new_name = self.args.load.split('/')[-2]
            else:
                new_name = self.args.load.split('/')[-1]
            self.args.name = new_name
            print(f'No --name arg provided. Name is being set to {new_name}')

    # STAND_ALONE SEQUENCE OUTPUT
    def perform_stand_alone_sequence_output(self):
        self.output_dir = os.path.abspath(
            os.path.join(self.symportal_root_directory, 'outputs', 'non_analysis', self.date_time_str))
        self._set_html_dir_and_js_out_path_from_output_dir()
        os.makedirs(self.html_dir, exist_ok=True)
        if self.args.print_output_seqs_sample_set:
            self._stand_alone_sequence_output_data_set_sample()
        else:
            self._stand_alone_sequence_output_data_set()
        self.number_of_samples = len(self.output_seq_count_table_obj.sorted_sample_uid_list)
        self._plot_if_not_too_many_samples(self._plot_sequence_stacked_bar_from_seq_output_table)
        self._do_sample_ordination()
        self._output_js_output_path_dict()
        self._print_all_outputs_complete()

    def _do_sample_ordination(self):
        # NB odinations are plot within the below function
        if not self.args.no_ordinations:
            self._do_sample_dist_and_pcoa()

    def _do_sample_dist_and_pcoa(self):
        print('\nCalculating between sample pairwise distances')
        if self.args.distance_method == 'both':
            if self._check_if_required_packages_found_in_path():
                self._do_unifrac_dist_pcoa()
                self._do_braycurtis_dist_pcoa()
            else:
                print('Changing distance method to braycurtis as one or more of the required '
                      'packages could not be found in your PATH')
                self.args.distance_method = 'braycurtis'
        if self.args.distance_method == 'unifrac':
            self._do_unifrac_dist_pcoa()
        elif self.args.distance_method == 'braycurtis':
            self._do_braycurtis_dist_pcoa()

    def _do_braycurtis_dist_pcoa(self):
        """NB this distance output is part of the perform_stand_alone_sequence_output
        and this can be called either with a data set input or a data set sample input.
        The distance PCoACreators have been made so that they can take either sort of input.
        As such, chance the input accordingly"""
        if self.args.print_output_seqs:
            # then we are working with a data set input
            braycurtis_dist_pcoa_creator = distance.SampleBrayCurtisDistPCoACreator(
                date_time_str=self.date_time_str,
                data_set_uid_list=[int(_) for _ in self.args.print_output_seqs.split(',')],
                output_dir=self.output_dir, html_dir=self.html_dir,
                js_output_path_dict=self.js_output_path_dict)
        elif self.args.print_output_seqs_sample_set:
            # then we are working with a data set sample input
            braycurtis_dist_pcoa_creator = distance.SampleBrayCurtisDistPCoACreator(
                date_time_str=self.date_time_str,
                data_set_sample_uid_list=[int(_) for _ in self.args.print_output_seqs_sample_set.split(',')],
                output_dir=self.output_dir, html_dir=self.html_dir,
                js_output_path_dict=self.js_output_path_dict)
        braycurtis_dist_pcoa_creator.compute_braycurtis_dists_and_pcoa_coords()
        if not self.args.no_figures:
            if len(braycurtis_dist_pcoa_creator.clade_col_uid_list) > 1000:
                print(
                    f'Too many samples ({len(braycurtis_dist_pcoa_creator.clade_col_uid_list)}) to generate plots')
            else:
                for output_path in braycurtis_dist_pcoa_creator.output_path_list:
                    if self.this_is_pcoa_path(output_path):
                        clade_of_output = os.path.dirname(output_path).split('/')[-1]
                        sys.stdout.write(f'\n\nGenerating between sample distance plot clade {clade_of_output}\n')
                        try:
                            dist_scatter_plotter_samples = plotting.DistScatterPlotterSamples(
                                csv_path=output_path, date_time_str=self.date_time_str)
                            dist_scatter_plotter_samples.make_sample_dist_scatter_plot()
                        except RuntimeError:
                            # The error message is printed to stdout at the source
                            continue

    def _do_unifrac_dist_pcoa(self):
        unifrac_dict_pcoa_creator = distance.SampleUnifracDistPCoACreator(
            date_time_str=self.date_time_str, output_dir=self.output_dir,
            data_set_uid_list=[int(_) for _ in self.args.print_output_seqs.split(',')],
            num_processors=self.args.num_proc, html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict)
        unifrac_dict_pcoa_creator.compute_unifrac_dists_and_pcoa_coords()
        if not self.args.no_figures:
            if len(unifrac_dict_pcoa_creator.cc_id_to_sample_id.keys()) > 1000:
                print(
                    f'Too many samples ({len(unifrac_dict_pcoa_creator.cc_id_to_sample_id.keys())}) to generate plots')
            else:
                for output_path in unifrac_dict_pcoa_creator.output_path_list:
                    if self.this_is_pcoa_path(output_path):
                        clade_of_output = os.path.dirname(output_path).split('/')[-1]
                        sys.stdout.write(f'\n\nGenerating between sample distance plot clade {clade_of_output}\n')
                        try:
                            dist_scatter_plotter_samples = plotting.DistScatterPlotterSamples(csv_path=output_path,
                                                                                     date_time_str=self.date_time_str)
                            dist_scatter_plotter_samples.make_sample_dist_scatter_plot()
                        except RuntimeError:
                            # The error message is printed to stdout at the source
                            continue

    @staticmethod
    def this_is_pcoa_path(output_path):
        return 'PCoA_coords' in output_path

    def _set_html_dir_and_js_out_path_from_output_dir(self):
        self.html_dir = os.path.join(self.output_dir, 'html')
        self.js_file_path = os.path.join(self.html_dir, 'study_data.js')
        os.makedirs(self.output_dir, exist_ok=True)
        self._set_logging_path()

    def _stand_alone_sequence_output_data_set(self):
        self.output_seq_count_table_obj = output.SequenceCountTableCreator(
            symportal_root_dir=self.symportal_root_directory, call_type='stand_alone',
            no_pre_med_seqs=self.args.no_pre_med_seqs,
            ds_uids_output_str=self.args.print_output_seqs,
            num_proc=self.args.num_proc, output_dir=self.output_dir, date_time_str=self.date_time_str,
            html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict, multiprocess=self.args.multiprocess)
        self.output_seq_count_table_obj.make_seq_output_tables()

    def _stand_alone_sequence_output_data_set_sample(self):
        self.output_seq_count_table_obj = output.SequenceCountTableCreator(
            symportal_root_dir=self.symportal_root_directory, call_type='stand_alone',
            no_pre_med_seqs=self.args.no_pre_med_seqs,
            dss_uids_output_str=self.args.print_output_seqs_sample_set,
            num_proc=self.args.num_proc, html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict,
            output_dir=self.output_dir, date_time_str=self.date_time_str, multiprocess=self.args.multiprocess)
        self.output_seq_count_table_obj.make_seq_output_tables()

    # STAND_ALONE TYPE OUTPUT

    def output_study_from_analysis(self):
        # Check that the provided Study object exists
        # Then run print_output_types_sample_set based on the data_set_samples attribute of the study

        self.study = self._try_to_get_study_object()
        
        # set the data_analysis attribute of the Study
        self._set_data_analysis_obj_from_arg_analysis_uid()
        self.study.data_analysis = self.data_analysis_object
        self.study.save()

        self.args.print_output_types_sample_set = ','.join([str(dss.id) for dss in self.study.data_set_samples.all()])

        # Now rejoin the logic flow for performing a type output as though it were a normal type output
        self.perform_stand_alone_type_output()


    def _try_to_get_study_object(self):
        try:
            return Study.objects.get(name=self.args.output_study_from_analysis)
        except ObjectDoesNotExist:
            pass
        try:
            return Study.objects.get(id=self.args.output_study_from_analysis)
        except ObjectDoesNotExist:
            raise RuntimeError(f'Cannot find Study {self.args.output_study_from_analysis} in the database.')

    def perform_stand_alone_type_output(self):
        self._set_data_analysis_obj_from_arg_analysis_uid()
        self.output_dir = os.path.join(
            self.symportal_root_directory, 'outputs', 'analyses', str(self.data_analysis_object.id), self.date_time_str)
        self._set_html_dir_and_js_out_path_from_output_dir()
        os.makedirs(self.html_dir, exist_ok=True)
        if self.args.print_output_types_sample_set:
            self._stand_alone_type_output_data_set_sample()
            self._stand_alone_seq_output_from_type_output_data_set_sample()
        else:
            self._stand_alone_type_output_data_set()
            self._stand_alone_seq_output_from_type_output_data_set()
        if not self.args.no_figures:
            self.number_of_samples = len(self.output_seq_count_table_obj.sorted_sample_uid_list)
            self._plot_if_not_too_many_samples(self._plot_sequence_stacked_bar_with_ordered_dss_uids_from_type_output)
            self._plot_if_not_too_many_samples(self._plot_type_stacked_bar_from_type_output_table)
        else:
            print('\nFigure plotting skipped at user\'s request')
        if not self.args.no_ordinations:
            self._do_data_analysis_ordinations()
        self._output_js_output_path_dict()

        if sp_config.system_type == 'remote' and self.args.output_study_from_analysis:
            try:
                self._output_study_output_info_items()
            except NotImplementedError as e:
                print(e)
        self._print_all_outputs_complete()

    def _output_study_output_info_items(self):
        """
        Produce the study_output_info.json file in the output directory
        and produce a .bak in the dbBackup directory
        """
        
        bak_path = os.path.join(self.dbbackup_dir, f'symportal_database_backup_{self.date_time_str}.bak')
        study_output_info_path = os.path.join(self.output_dir, 'study_output_info.json')
        # Now output the .json file
        temp_dict = {}
        temp_dict['bak_path'] = bak_path
        temp_dict["time_stamp_str"] = self.date_time_str
        temp_dict["study"] = self.study.name

        # Set the display_online and the data_explorer attribute of the study to True
        # Also set analysis to True
        self.study.display_online = True
        self.study.data_explorer = True
        if self.args.output_study_from_analysis:
            self.study.analysis = True
        elif self.args.load:
            # This is already set as False as default but let's be explicit
            self.study.analysis = False
        self.study.save()

        with open(study_output_info_path, 'w') as f:
            json.dump(obj=temp_dict, fp=f)

        print(f'study_output_info items:\n'
              f'\t{study_output_info_path}\n'
              f'\t{bak_path}')

    def _stand_alone_seq_output_from_type_output_data_set(self):
        self.output_seq_count_table_obj = output.SequenceCountTableCreator(
            call_type='analysis',
            num_proc=self.args.num_proc,
            symportal_root_dir=self.symportal_root_directory,
            no_pre_med_seqs=self.args.no_pre_med_seqs,
            ds_uids_output_str=self.args.print_output_types,
            output_dir=self.output_dir,
            sorted_sample_uid_list=self.output_type_count_table_obj.sorted_list_of_vdss_uids_to_output,
            analysis_obj=self.data_analysis_object,
            date_time_str=self.date_time_str, html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict,
            multiprocess=self.args.multiprocess)
        self.output_seq_count_table_obj.make_seq_output_tables()

    def _stand_alone_seq_output_from_type_output_data_set_sample(self):
        self.output_seq_count_table_obj = output.SequenceCountTableCreator(
            call_type='analysis',
            num_proc=self.args.num_proc,
            symportal_root_dir=self.symportal_root_directory,
            no_pre_med_seqs=self.args.no_pre_med_seqs,
            dss_uids_output_str=self.args.print_output_types_sample_set,
            output_dir=self.output_dir,
            sorted_sample_uid_list=self.output_type_count_table_obj.sorted_list_of_vdss_uids_to_output,
            analysis_obj=self.data_analysis_object,
            date_time_str=self.date_time_str, html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict,
            multiprocess=self.args.multiprocess)
        self.output_seq_count_table_obj.make_seq_output_tables()

    def _stand_alone_type_output_data_set(self):
        ds_uid_list = [int(ds_uid_str) for ds_uid_str in self.args.print_output_types.split(',')]
        self._check_ds_were_part_of_analysis(ds_uid_list)
        self.output_type_count_table_obj = output.OutputProfileCountTable(
            num_proc=self.args.num_proc, within_clade_cutoff=self.within_clade_cutoff,
            call_type='stand_alone', date_time_str=self.date_time_str,
            data_set_uids_to_output=set(ds_uid_list), data_analysis_obj=self.data_analysis_object,
            output_dir=self.output_dir, html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict,
            force_basal_lineage_separation=self.args.force_basal_lineage_separation)
        self.output_type_count_table_obj.output_types()

    def _check_ds_were_part_of_analysis(self, ds_uid_list):
        for ds_uid in ds_uid_list:
            if ds_uid not in [int(uid_str) for uid_str in self.data_analysis_object.list_of_data_set_uids.split(',')]:
                print(f'DataSet UID: {ds_uid} is not part of analysis: {self.data_analysis_object.name}')
                raise RuntimeError

    def _stand_alone_type_output_data_set_sample(self):
        dss_uid_list = [int(dss_uid_str) for dss_uid_str in self.args.print_output_types_sample_set.split(',')]
        self._check_dss_were_part_of_analysis(dss_uid_list)
        self.output_type_count_table_obj = output.OutputProfileCountTable(
            num_proc=self.args.num_proc, within_clade_cutoff=self.within_clade_cutoff,
            call_type='stand_alone', output_dir=self.output_dir, html_dir=self.html_dir,
            js_output_path_dict=self.js_output_path_dict, date_time_str=self.date_time_str,
            data_set_sample_uid_set_to_output=set(dss_uid_list), data_analysis_obj=self.data_analysis_object,
            force_basal_lineage_separation=self.args.force_basal_lineage_separation)
        self.output_type_count_table_obj.output_types()

    def _check_dss_were_part_of_analysis(self, dss_uid_list):
        ds_uid_list_for_query = [int(a) for a in self.data_analysis_object.list_of_data_set_uids.split(',')]
        ds_of_analysis = self._chunk_query_ds_objs_from_ds_uids(ds_uid_list_for_query)
        dss_of_analysis = self._chunk_query_dss_objs_from_ds_uids(ds_of_analysis)
        dss_uids_that_were_part_of_analysis = [dss.id for dss in dss_of_analysis]
        for dss_uid in dss_uid_list:
            if dss_uid not in dss_uids_that_were_part_of_analysis:
                print(f'DataSetSample UID: {dss_uid} was not part of DataAnalysis: {self.data_analysis_object.name}')
                raise RuntimeError

    def _chunk_query_dss_objs_from_ds_uids(self, ds_of_analysis):
        dss_of_analysis = []
        for uid_list in self.thread_safe_general.chunks(ds_of_analysis):
            dss_of_analysis.extend(list(DataSetSample.objects.filter(data_submission_from__in=uid_list)))
        return dss_of_analysis

    def _chunk_query_ds_objs_from_ds_uids(self, ds_uid_list_for_query):
        ds_of_analysis = []
        for uid_list in self.thread_safe_general.chunks(ds_uid_list_for_query):
            ds_of_analysis.extend(list(DataSet.objects.filter(id__in=uid_list)))
        return ds_of_analysis

    # ITS2 TYPE PROFILE STAND_ALONE DISTANCES
    def perform_type_distance_stand_alone(self):
        """Generate the within clade pairwise distances between ITS2 type profiles
        either using a BrayCurtis- or Unifrac-based
        method. Also calculate the PCoA and plot as scatter plot for each."""
        self._verify_data_analysis_uid_provided()
        self._set_data_analysis_obj_from_arg_analysis_uid()
        self.run_type_distances_dependent_on_methods()
        if self.args.distance_method == 'both':
            self._plot_type_distances_from_distance_object(self.braycurtis_distance_object)
            self._plot_type_distances_from_distance_object(self.unifrac_distance_object)
        elif self.args.distance_method == 'unifrac':
            self._plot_type_distances_from_distance_object(self.unifrac_distance_object)
        elif self.args.distance_method == 'braycurtis':
            self._plot_type_distances_from_distance_object(self.braycurtis_distance_object)
        self._output_js_output_path_dict()
        self._print_all_outputs_complete()

    def run_type_distances_dependent_on_methods(self):
        """Start an instance of the correct distance class running."""
        self.output_dir = os.path.join(
                    self.symportal_root_directory, 'outputs', 'ordination', self.date_time_str)
        self._set_html_dir_and_js_out_path_from_output_dir()
        os.makedirs(self.html_dir, exist_ok=True)
        if self.args.distance_method == 'both':
            if self._check_if_required_packages_found_in_path():
                # then do both distance outputs
                if self.args.between_type_distances_sample_set:
                    self._start_type_braycurtis_data_set_samples()
                    self._start_type_unifrac_data_set_samples()
                elif self.args.between_type_distances_cct_set:
                    self._start_type_braycurtis_cct_set()
                    self._start_type_unifrac_cct_set()
                else:
                    self._start_type_braycurtis_data_sets()
                    self._start_type_unifrac_data_sets()
            else:
                print('Changing distance method to braycurtis as one or more of the required '
                      'packages could not be found in your PATH')
                self.args.distance_method = 'braycurtis'

        if self.args.distance_method == 'unifrac':
            if self.args.between_type_distances_sample_set:
                self._start_type_unifrac_data_set_samples()
            elif self.args.between_type_distances_cct_set:
                self._start_type_unifrac_cct_set()
            else:
                self._start_type_unifrac_data_sets()
        elif self.args.distance_method == 'braycurtis':
            if self.args.between_type_distances_sample_set:
                self._start_type_braycurtis_data_set_samples()
            elif self.args.between_type_distances_cct_set:
                self._start_type_braycurtis_cct_set()
            else:
                self._start_type_braycurtis_data_sets()

    @staticmethod
    def _check_if_required_packages_found_in_path():
        """For creating unifrac-derived distances we need
        both iqtree and mafft to be install in the users PATH.
        Here we will check for them. If either of them are not found we will return False"""
        return which('iqtree') and which('mafft')

    # BRAYCURTIS between its2 type profile distance methods
    def _start_type_braycurtis_cct_set(self):
        self.braycurtis_distance_object = distance.TypeBrayCurtisDistPCoACreator(
            data_analysis_obj=self.data_analysis_object,
            date_time_str=self.date_time_str,
            cct_set_uid_list=[int(cct_uid_str) for cct_uid_str in self.args.between_type_distances_cct_set.split(',')],
            local_abunds_only=False, html_dir=self.html_dir,
            output_dir=self.output_dir, js_output_path_dict=self.js_output_path_dict
        )
        self.braycurtis_distance_object.compute_braycurtis_dists_and_pcoa_coords()

    def _start_type_braycurtis_data_sets(self):
        self.braycurtis_distance_object = distance.TypeBrayCurtisDistPCoACreator(
            data_analysis_obj=self.data_analysis_object,
            date_time_str=self.date_time_str,
            data_set_uid_list=[int(ds_uid_str) for ds_uid_str in self.args.between_type_distances.split(',')],
            local_abunds_only=self.args.local, html_dir=self.html_dir,
            output_dir=self.output_dir, js_output_path_dict=self.js_output_path_dict
        )
        self.braycurtis_distance_object.compute_braycurtis_dists_and_pcoa_coords()

    def _start_type_braycurtis_data_set_samples(self):
        self.braycurtis_distance_object = distance.TypeBrayCurtisDistPCoACreator(
            data_analysis_obj=self.data_analysis_object,
            date_time_str=self.date_time_str,
            data_set_sample_uid_list=[int(ds_uid_str) for ds_uid_str in
                                      self.args.between_type_distances_sample_set.split(',')],
            local_abunds_only=self.args.local, html_dir=self.html_dir,
            js_output_path_dict=self.js_output_path_dict, output_dir=self.output_dir
        )
        self.braycurtis_distance_object.compute_braycurtis_dists_and_pcoa_coords()

    # UNIFRAC between its2 type profile distance methods
    def _start_type_unifrac_cct_set(self):
        self.unifrac_distance_object = distance.TypeUnifracDistPCoACreator(
            data_analysis_obj=self.data_analysis_object,
            date_time_str=self.date_time_str,
            num_processors=self.args.num_proc,
            cct_set_uid_list=[int(cct_uid_str) for cct_uid_str in
                              self.args.between_type_distances_cct_set.split(',')],
            local_abunds_only=False, html_dir=self.html_dir,
            output_dir=self.output_dir,
            js_output_path_dict=self.js_output_path_dict
        )
        self.unifrac_distance_object.compute_unifrac_dists_and_pcoa_coords()

    def _start_type_unifrac_data_sets(self):
        self.unifrac_distance_object = distance.TypeUnifracDistPCoACreator(
            data_analysis_obj=self.data_analysis_object,
            date_time_str=self.date_time_str,
            num_processors=self.args.num_proc,
            data_set_uid_list=[int(ds_uid_str) for ds_uid_str in
                               self.args.between_type_distances.split(',')],
            local_abunds_only=self.args.local, output_dir=self.output_dir,
            html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict
        )
        self.unifrac_distance_object.compute_unifrac_dists_and_pcoa_coords()

    def _start_type_unifrac_data_set_samples(self):
        self.unifrac_distance_object = distance.TypeUnifracDistPCoACreator(
            data_analysis_obj=self.data_analysis_object,
            date_time_str=self.date_time_str,
            num_processors=self.args.num_proc,
            data_set_sample_uid_list=[int(ds_uid_str) for ds_uid_str in
                                      self.args.between_type_distances_sample_set.split(',')],
            local_abunds_only=self.args.local, output_dir=self.output_dir,
            html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict
        )
        self.unifrac_distance_object.compute_unifrac_dists_and_pcoa_coords()

    # SAMPLE STAND_ALONE DISTANCES
    def _perform_sample_distance_stand_alone(self):
        self.output_dir = os.path.join(
            self.symportal_root_directory, 'outputs', 'ordination', self.date_time_str)
        self._set_html_dir_and_js_out_path_from_output_dir()
        os.makedirs(self.html_dir, exist_ok=True)
        self._run_sample_distances_dependent_on_methods()
        if self.args.distance_method == 'both':
            self._plot_sample_distances_from_distance_object(self.braycurtis_distance_object)
            self._plot_sample_distances_from_distance_object(self.unifrac_distance_object)
        elif self.args.distance_method == 'unifrac':
            self._plot_sample_distances_from_distance_object(self.unifrac_distance_object)
        elif self.args.distance_method == 'braycurtis':
            self._plot_sample_distances_from_distance_object(self.braycurtis_distance_object)
        self._output_js_output_path_dict()
        self._print_all_outputs_complete()

    def _run_sample_distances_dependent_on_methods(self):
        """Start an instance of the correct distance class running."""
        if self.args.distance_method == 'both':
            if self._check_if_required_packages_found_in_path():
                if self.args.between_sample_distances_sample_set:
                    self._start_sample_braycurtis_data_set_samples()
                    self._start_sample_unifrac_data_set_samples()
                else:
                    self._start_sample_braycurtis_data_sets()
                    self._start_sample_unifrac_data_sets()
            else:
                print('Changing distance method to braycurtis as one or more of the required '
                      'packages could not be found in your PATH')
                self.args.distance_method = 'braycurtis'
        if self.args.distance_method == 'unifrac':
            if self.args.between_sample_distances_sample_set:
                self._start_sample_unifrac_data_set_samples()
            else:
                self._start_sample_unifrac_data_sets()
        elif self.args.distance_method == 'braycurtis':
            if self.args.between_sample_distances_sample_set:
                self._start_sample_braycurtis_data_set_samples()
            else:
                self._start_sample_braycurtis_data_sets()

    @staticmethod
    def _print_all_outputs_complete():
        print('\n\nALL OUTPUTS COMPLETE\n\n')

    def _start_sample_unifrac_data_set_samples(self):
        dss_uid_list = [int(ds_uid_str) for ds_uid_str in self.args.between_sample_distances_sample_set.split(',')]
        self.unifrac_distance_object = distance.SampleUnifracDistPCoACreator(
            data_set_sample_uid_list=dss_uid_list,
            num_processors=self.args.num_proc,
            output_dir=self.output_dir,
            html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict, date_time_str=self.date_time_str)
        self.unifrac_distance_object.compute_unifrac_dists_and_pcoa_coords()

    def _start_sample_unifrac_data_sets(self):
        ds_uid_list = [int(ds_uid_str) for ds_uid_str in self.args.between_sample_distances.split(',')]
        self.unifrac_distance_object = distance.SampleUnifracDistPCoACreator(
            data_set_uid_list=ds_uid_list,
            num_processors=self.args.num_proc,
            output_dir=self.output_dir,
            html_dir=self.html_dir, js_output_path_dict=self.js_output_path_dict, date_time_str=self.date_time_str)
        self.unifrac_distance_object.compute_unifrac_dists_and_pcoa_coords()

    def _start_sample_braycurtis_data_set_samples(self):
        dss_uid_list = [int(ds_uid_str) for ds_uid_str in self.args.between_sample_distances_sample_set.split(',')]
        self.braycurtis_distance_object = distance.SampleBrayCurtisDistPCoACreator(
            date_time_str=self.date_time_str,
            data_set_sample_uid_list=dss_uid_list,
            output_dir=self.output_dir, html_dir=self.html_dir,
            js_output_path_dict=self.js_output_path_dict)
        self.braycurtis_distance_object.compute_braycurtis_dists_and_pcoa_coords()

    def _start_sample_braycurtis_data_sets(self):
        ds_uid_list = [int(ds_uid_str) for ds_uid_str in self.args.between_sample_distances.split(',')]
        self.braycurtis_distance_object = distance.SampleBrayCurtisDistPCoACreator(
            date_time_str=self.date_time_str,
            data_set_uid_list=ds_uid_list,
            output_dir=self.output_dir, html_dir=self.html_dir,
            js_output_path_dict=self.js_output_path_dict)
        self.braycurtis_distance_object.compute_braycurtis_dists_and_pcoa_coords()

    # APPLY DATASHEET TO DATASETSAMPLES
    def apply_datasheet_to_dataset_samples(self):
        try:
            adtdss = django_general.ApplyDatasheetToDataSetSamples(data_set_uid=self.args.apply_data_sheet,
                                                                   data_sheet_path=self.args.data_sheet)
        except RuntimeError as e:
            print(e)
            return
        adtdss.apply_datasheet()

    # VACUUM DB
    def perform_vacuum_database(self):
        print('Vacuuming database')
        self.vacuum_db()
        print('Vacuuming complete')

    @staticmethod
    def vacuum_db():
        from django.db import connection
        cursor = connection.cursor()
        cursor.execute("VACUUM")
        connection.close()

    # DISPLAY DB CONTENTS FUNCTIONS
    @staticmethod
    def perform_display_data_sets():
        data_set_id_to_obj_dict = {ds.id: ds for ds in list(DataSet.objects.all())}
        sorted_list_of_ids = sorted(list(data_set_id_to_obj_dict.keys()))
        for ds_id in sorted_list_of_ids:
            ds_in_q = data_set_id_to_obj_dict[ds_id]
            print(f'{ds_in_q.id}: {ds_in_q.name}\t{ds_in_q.time_stamp}')

    @ staticmethod
    def perform_display_studies():
        ordered_studies = Study.objects.order_by('id')
        for study in ordered_studies:
            print(f'{study.id}\t{study.name}: {",".join([user.name for user in study.user_set.all()])}')
            # print('Users:')
            # for user in study.user_set.all():
            #     print(f'{user.id}\t{user.name}')
            # print('\n\n')

    @staticmethod
    def perform_display_analysis_types():
        data_analysis_id_to_obj_dict = {da.id: da for da in list(DataAnalysis.objects.all())}
        sorted_list_of_ids = sorted(list(data_analysis_id_to_obj_dict.keys()))
        for da_id in sorted_list_of_ids:
            da_in_q = data_analysis_id_to_obj_dict[da_id]
            print(f'{da_in_q.id}: {da_in_q.name}\t{da_in_q.time_stamp}')

    def _set_logging_path(self):
        logging.basicConfig(
            format='%(levelname)s:%(message)s',
            filename=os.path.join(self.output_dir, f'{self.date_time_str}_log.log'),
            filemode='w', level=logging.INFO)
        logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

class CitationUpdate:
    def __init__(self):
        # Try to get the citing pubs live from Google Scholar
        # If this fails fall back to loadinging in the .csv files that have been downloaded manually
        try:
            import compress_pickle
            self.citing_pubs = compress_pickle.load('citing_pubs.p.bz')
            # self.citing_pubs = self._get_citing_pubs()
            # temporarily lets pickle this out so that we can debug without having to do the requrests
            # compress_pickle.dump(self.citing_pubs, 'citing_pubs.p.bz')
        except Exception as e:
            print(e)
            print('Loading citing articles from .csv')
            raise NotImplementedError
        self.current_study_objects_query = Study.objects.all()
        self.current_citation_objects_query = Citation.objects.all()
        self.current_citation_titles = [citation.title for citation in self.current_citation_objects_query]

    def _get_citing_pubs(self):
        """
        Attempt to query google scholar to get a list of the articles that cite the SymPortal paper
        """
        citing_pubs = []
        headers = {
            'accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
            'accept-encoding': 'gzip, deflate, br',
            'accept-language': 'en-GB,en-US;q=0.9,en;q=0.8',
            'cache-control': 'max-age=0',
            'cookie': 'CONSENT=YES+GB.en+20150712-15-0; SEARCH_SAMESITE=CgQIwo8B; ANID=AHWqTUm-NmPN57OwDDERHKOBw9y6CVGj4lBjL8wM6rkfQb-n7712Ab9vofLH9IrY; SID=zgeBmCoAbuLSe5Ye_TQZADF7wxrgPyoajoBT2tYif1Ie7ZBIR3fqzF4yJeXQkrGzQD3PcA.; __Secure-3PSID=zgeBmCoAbuLSe5Ye_TQZADF7wxrgPyoajoBT2tYif1Ie7ZBId5hEurpvAVRB939sr0wK-w.; HSID=AcLCYk7FF8YbFq8YU; SSID=A3lWM-dBOlKYk160r; APISID=jVhpX9WDFTTPhfFk/Aw6hW_bfATMuUOkbm; SAPISID=u4BlUiVIy486neYK/AtnzEUHzoLnfPlVfN; __Secure-HSID=AcLCYk7FF8YbFq8YU; __Secure-SSID=A3lWM-dBOlKYk160r; __Secure-APISID=jVhpX9WDFTTPhfFk/Aw6hW_bfATMuUOkbm; __Secure-3PAPISID=u4BlUiVIy486neYK/AtnzEUHzoLnfPlVfN; __Secure-3PSIDCC=AJi4QfE3goWyj5XBPQfnMcFvaV0-jwyVW1d6uGOr--vNC6Og-FodbE1XZNBb-hdAleu9_YEg; GSP=IN=acfa7a445dbbf298:LD=en:NR=20:LM=1595959706:S=dIrB7oXzhJr2pQWU; NID=204=a0Fv1TbEg7cZ6ADV_FJB5vi6qMg0XL_Cs-Vb-Gqh4aQeGQH04XKhUrbqK5GCKOoLBg7jQ54rguhthZQHSsnqKW8_UlVM3nk0KfaNKN7v3h3upA21XH8HzB8W5Iu8uYs5Di1VOG7mK4YuGxKgP1HcVUQw57NVJOx72WLPYcvfy4pBudQB28P6nKcVCYdjdSZpoxJb6AEabOFC4JZgSREpldfvm0HKPb2tdlBax8HYPY-HmKlcQwagdAyZQZolhJd2ZEaFegppFSSjTRCL86kJkbFqLKAOnZR_soXQFvshLU-g_3BF4-DvEns5__KYRSVR1zdUnFN06H61pub_uipBpw; 1P_JAR=2020-7-29-5; SIDCC=AJi4QfGwhs-FXe-vApcaO2CwX_VFBaw0xVrk1drZxmV4zz5zsV_6VDFD6CjEQ8lLupG-nWRq-M0',
            'referer': 'https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=SymPortal+novel+analytical+framework&btnG=&oq=symportal',
            'sec-fetch-dest': 'document',
            'sec-fetch-mode': 'navigate',
            'sec-fetch-site': 'same-origin',
            'sec-fetch-user': '?1',
            'upgrade-insecure-requests': '1',
            'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/84.0.4147.89 Safari/537.36',
            'x-client-data': 'CKC1yQEIhbbJAQijtskBCMS2yQEIqZ3KAQjyoMoBCI2sygEI/7zKAQjnyMoB'
        }

        current_url = 'https://scholar.google.com/scholar?cites=12056348124469462278&amp;as_sdt=2005&amp;sciodt=0,5&amp;hl=en&amp;oe=ASCII'
        soup = BeautifulSoup(requests.get(
            current_url,
            headers=headers).text, features="html.parser")

        # Find out how many pagntions of articles there are
        # And grab the articles for each page.
        pags = soup.find_all(class_='gs_nma')
        assert(len(pags) != 0)
        for pag in pags:
            if int(pag.text) != 1:
                # Then we need to create a new soup
                # Lets wait for 30s before making the query
                time.sleep(30)
                full_url = 'https://scholar.google.com' + pag['href']
                headers['referer'] = current_url
                current_url = full_url
                soup = BeautifulSoup(requests.get(full_url, headers=headers).text, features="html.parser")

            articles = soup.find_all('div', class_='gs_ri')
            assert(len(articles) != 0)
            year_re = re.compile('20\d{2}')
            for article_div in articles:
                title = article_div.find('h3').find('a').text
                authors = article_div.find(class_='gs_a').text.split('-')[0].rstrip()
                year_search = year_re.search(article_div.find(class_='gs_a').text)
                if year_search is not None:
                    year = year_search.group()
                else:
                    year = ''
                article_url = article_div.find('h3').find('a')['href']
                citing_pubs.append({'title': title, 'authors': authors, 'year': year, 'article_url': article_url})
        # Once we have been through the first set of articles we will need to do another set of articles if there is one
        return citing_pubs

    def update(self):
        """
        For each of the citing_pubs see if there is a matching Study or Citation object
        If there is not, then we will either create a Citation object, associate to an existing Study object
        or create a new study object.
        The Study and Citation objects will then be used by the symportal.org code to populate the published articles
        using SymPortal table. The citation list will be populated directly using the scholarly API.
        """
        for citation in self.citing_pubs:
            while True:
                print(f'Checking for: "{citation["title"]}"')
                print(f'\t{citation["authors"]} {citation["year"]}')
                if citation["title"] in self.current_citation_titles:
                    print('\tFound a Citation object matching the citation')
                else:
                    print('\tCitation does not match a Citation object')
                    self._deal_with_unassociated_citation(citation)

                if self._update_db_objs_and_verify_citation_association(citation=citation):
                    break


    def _update_db_objs_and_verify_citation_association(self, citation):
        # At this point we should have either a Study or citation object created
        # Verify then break
        self.current_study_objects_query = Study.objects.all()
        self.current_citation_objects_query = Citation.objects.all()
        self.current_citation_titles = [citation.title for citation in self.current_citation_objects_query]
        if citation['title'] not in self.current_citation_titles:
            raise RuntimeError(f'Citation {citation.title} still not associated to a Study or Citation object')
        else:
            return True

    def _deal_with_unassociated_citation(self, citation):
        """
        1 = No Study object. Only create Citation object
        2 = Create Study object. Create Citation object if doesn't exist
        3 = Associate to existing study object. Create Citation object if doesn't exist
        :param citation:
        :return:
        """
        option = self._present_user_options(citation)
        if option == '1':
            self._create_new_citation_object(citation)
        elif option == '2':
            self._create_new_study_object(citation)
        elif option == '3':
            self._associate_to_existing_study(citation)
        else:
            raise RuntimeError

    def _present_user_options(self, citation):
        # Present the User with the Study that has the most similar title according to a damerau_levenshtein
        # metric
        study_with_most_similar_title = self._get_most_similar_study(citation)
        print('The most similar Study to the citation in question is: ')
        print(f'< Study: id {study_with_most_similar_title.id}, name {study_with_most_similar_title.name} >')
        print(f'title: {study_with_most_similar_title.title}')
        while True:
            option = input('\n\nPlease select one of the following options: \n'
                           '1 - Citation should not have a Study object associated to it\n'
                           '\t(a Citation object will still be created)\n'
                           '2 - Create a new Study object to associate the citation to\n'
                           '3 - Associate the citation to an existing Study object\n\n[1,2,3]: ')
            if option not in list('123'):
                continue
            else:
                break
        return option

    def _create_new_citation_object(self, citation):
        # Create a citation object
        new_citation = Citation(
            title=citation["title"], article_url=citation['article_url'],
            author_list_string=citation['authors'], year=citation['year']
        )
        new_citation.save()

    def _create_new_study_object(self, citation):
        # Create a new Study object
        csaau = CreateStudyAndAssociateUsers(citation=citation)
        csaau.create_study_and_user_objects()
        self._create_new_citation_object_if_no_exist(citation)

    def _create_new_citation_object_if_no_exist(self, citation):
        if citation['title'] not in self.current_citation_titles:
            print('No Citation object found. Creating new Citation object')
            self._create_new_citation_object(citation=citation)
        else:
            print('Citation object already exists with same title.')

    def _associate_to_existing_study(self, citation):
        while True:
            study_uid = input('Study UID: ')
            try:
                existing_study = Study.objects.get(id=int(study_uid))
                print(f'Study is {existing_study}')
                print(f'\ttitle: {existing_study.title}')
                continue_text = input(f'Study is {existing_study}. Is this correct? [y/n]: ')
                if continue_text == 'y':
                    break
                else:
                    continue
            except ObjectDoesNotExist as e:
                print(e)
                print('Study not found. Please try again.')
        print(f'Updating Study {existing_study}')
        existing_study.title = citation["title"]
        existing_study.author_list_string = citation["authors"]
        existing_study.article_url = citation["article_url"]
        existing_study.is_published = True
        existing_study.display_online = True
        existing_study.save()
        self._create_new_citation_object_if_no_exist(citation=citation)

    def _get_most_similar_study(self, citation):
        max_sim = 0
        max_sim_obj = None
        for study in self.current_study_objects_query:
            sim = textdistance.damerau_levenshtein.normalized_similarity(citation["title"], study.title)
            if sim > max_sim:
                max_sim = sim
                max_sim_obj = study
        return max_sim_obj


if __name__ == "__main__":
    spwfm = SymPortalWorkFlowManager()
    spwfm.start_work_flow()

