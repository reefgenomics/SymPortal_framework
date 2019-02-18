#!/usr/bin/env python3.6
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
    https://github.com/SymPortal/SymPortal_framework/tree/master/LICENSE.txt.
    """

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
from dbApp.models import DataSet, DataAnalysis
############################################

import data_sub_collection_run
import output
import plotting
import sys
import distance
import argparse
import data_loading
import sp_config

class SymPortalWorkFlowManager:
    def __init__(self, custom_args_list=None):
        self.args = self._define_args(custom_args_list)
        # general attributes
        self.symportal_root_directory = os.path.abspath(os.path.dirname(__file__))
        self.date_time_str = str(datetime.now()).replace(' ', '_').replace(':', '-')
        self.submitting_user = sp_config.user_name
        self.submitting_user_email = sp_config.user_email

        # for data_loading
        self.data_set_object = None
        self.screen_sub_eval_bool = None
        if sp_config.system_type == 'remote':
            self.screen_sub_eval_bool = True
        else:
            self.screen_sub_eval_bool = False
        self.reference_db = 'symClade.fa'

        # for data analysis
        self.within_clade_cutoff = 0.03
        self.data_analysis_object = None

        # for dist and pcoa outputs
        self.pcoa_output_path_list = []

    def _define_args(self, custom_args_list=None):
        parser = argparse.ArgumentParser(
            description='Intragenomic analysis of the ITS2 region of the nrDNA',
            epilog='For support email: symportal@gmail.com')
        group = parser.add_mutually_exclusive_group(required=True)
        self.define_mutually_exclusive_args(group)
        self.define_additional_args(group, parser)
        if custom_args_list is not None:
            return parser.parse_args(custom_args_list)
        else:
            return parser.parse_args()

    def define_additional_args(self, group, parser):
        parser.add_argument('--num_proc', type=int, help='Number of processors to use', default=1)
        parser.add_argument('--name', help='A name for your input or analysis', default='noName')
        parser.add_argument('--description', help='An optional description', default='No description')
        parser.add_argument('--data_analysis_id', type=int, help='The ID of the data_analysis you wish to output from')
        group.add_argument(
            '--vacuum_database', action='store_true',
            help='Vacuuming the database will free up memory from objects that have been deleted recently')
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
            help='Either \'unifrac\' or \'braycurtis\', default=braycurtis. The method to use when '
                 'calculating distances between its2 type profiles or samples.', default='braycurtis')
        # when run as remote
        parser.add_argument(
            '--submitting_user_name',
            help='Only for use when running as remote\nallows the association of a different user_name to the '
                 'data_set than the one listed in sp_config', default='not supplied')

        parser.add_argument(
            '--submitting_user_email',
            help='Only for use when running as remote\nallows the association of a different user_email to the data_set '
                 'than the one listed in sp_config', default='not supplied')

    def define_mutually_exclusive_args(self, group):
        group.add_argument(
            '--load', metavar='path_to_dir',
            help='Run this to load data to the framework\'s database. The first argument to this command must be an '
                 'absolute path to a directory containing  the paired sequencing reads in .fastq.gz format. Alternatively, '
                 'this path can point directly to a single compressed file containing the same paired fastq.gz files. '
                 '\nA name must be associated with the data_set using the --name flag. \nThe number of processes to use '
                 'can also be specified using the --num_proc flag. \nA datasheet can also be uploaded using the '
                 '--data_sheet flag and the full path to the .xlsx data_sheet file (RECOMMENDED). \n'
                 'To skip the generation of figures pass the --no_figures flag.\n To skip the generation of '
                 'ordination files (pairwise distances and PCoA coordinates) pass the --no_ordinations flag')
        group.add_argument(
            '--analyse', metavar='data_set uids',
            help='Analyse one or more data_set objects together. Enter comma separated uids of the data_set uids you '
                 'which to analyse. e.g.: 43,44,45. If you wish to use all available dataSubmissions, you may pass '
                 '\'all\' as an argument. To display all data_sets currently submitted to the framework\'s database, '
                 'including their ids, use the \'show_data_sets\' command\nTo skip the generation of figures pass the '
                 '--no_figures flag.\nTo skip the generation of ordination files (pairwise distances and PCoA coordinates) '
                 'pass the --no_ordinations flag')
        group.add_argument(
            '--display_data_sets', action='store_true', help='Display data_sets currently in the framework\'s database')
        group.add_argument(
            '--display_analyses', action='store_true',
            help=' Display data_analysis objects currently stored in the framework\'s database')
        group.add_argument(
            '--print_output_seqs', metavar='data_set uids',
            help='Use this function to output ITS2 sequence count tables for given data_set instances')
        group.add_argument(
            '--print_output_seqs_sample_set', metavar='DataSetSample UIDs',
            help='Use this function to output ITS2 sequence count tables for a collection of DataSetSample instances. '
                 'The input to this function should be a comma separated string of the UIDs of the DataSetSample instances '
                 'in question. e.g. 345,346,347,348')
        group.add_argument(
            '--print_output_types', metavar='data_set uids, analysis ID',
            help='Use this function to output the ITS2 sequence and ITS2 type profile count tables for a given set of '
                 'data_sets that have been run in a given analysis. Give the data_set uids that you wish to make outputs '
                 'for as arguments to the --print_output_types flag. To output for multiple data_set objects, '
                 'comma separate the uids of the data_set objects, e.g. 44,45,46. Give the ID of the analysis you wish to '
                 'output these from using the --data_analysis_id flag.\nTo skip the generation of figures pass the '
                 '--no_figures flag.')
        group.add_argument(
            '--print_output_types_sample_set', metavar='DataSet UIDs, analysis UID',
            help='Use this function to output the ITS2 sequence and ITS2 type profile count tables for a given set of '
                 'DataSetSample objects that have been run in a given DataAnalysis. Give the DataSetSample '
                 'UIDs that you wish to make outputs from as arguments to the --print_output_types flag. To output for '
                 'multiple DataSetSample objects, comma separate the UIDs of the DataSetSample objects, '
                 'e.g. 5644,5645,5646. Give the UID of the DataAnalysis you wish to output these from using the '
                 '--data_analysis_id flag.\nTo skip the generation of figures pass the '
                 '--no_figures flag.')
        group.add_argument(
            '--between_type_distances', metavar='data_set uids, analysis ID',
            help='Use this function to output UniFrac pairwise distances between ITS2 type profiles clade separated')
        group.add_argument(
            '--between_sample_distances', metavar='data_set uids',
            help='Use this function to output UniFrac pairwise distances between samples clade separated from a '
                 'given collection of data_set objects')
        group.add_argument(
            '--between_sample_distances_sample_set', metavar='data_set_sample uids',
            help='Use this function to output UniFrac pairwise distances between samples clade '
                 'separated from a given collection of data_set_sample objects')

    def start_work_flow(self):
        if self.args.load:
            return self.perform_data_loading()
        elif self.args.analyse:
            self.perform_data_analysis()
        elif self.args.print_output_seqs:
            self.perform_sequences_count_table_output()
        elif self.args.print_output_seqs_sample_set:
            self.perform_sequences_count_table_output()
        elif self.args.print_output_types:
            # TODO class abstract and integrate type output by types
            self.perform_type_count_table_output()
        elif self.args.print_output_types_sample_set:
            # TODO class abstract and integrate type output by sample uid
            self.perform_type_count_table_output_sample_uid()
        elif self.args.display_data_sets:
            self.perform_display_data_sets()
        elif self.args.display_analyses:
            self.perform_display_analysis_types()
        elif self.args.between_type_distances:
            #TODO type distances sample uid input
            self.perform_within_clade_type_distance_generation()
        elif self.args.between_sample_distances:
            self.perform_within_clade_sample_distance_generation()
        elif self.args.between_sample_distances_sample_set:
            self.perform_within_clade_sample_distance_generation()
        elif self.args.vacuumDatabase:
            self.perform_vacuum_database()

    def perform_type_count_table_output_sample_uid(self):
        self.verify_data_analysis_uid_provided()
        analysis_object = self.get_data_analysis_by_uid()
        output.output_type_count_tables_data_set_sample_id_input(
            analysisobj=analysis_object, num_processors=self.args.num_proc,
            data_set_sample_ids_to_output_string=self.args.print_output_types_sample_set,
            no_figures=self.args.no_figures,
            output_user=self.submitting_user)

    # data analysis methods
    def perform_data_analysis(self):
        self._verify_name_arg_given()
        self.create_new_data_analysis_obj()
        self.start_data_analysis()

    def create_new_data_analysis_obj(self):
        self.data_analysis_object = DataAnalysis(
            list_of_data_set_uids=self.args.analyse, within_clade_cutoff=self.within_clade_cutoff,
            name=self.args.name, time_stamp=self.date_time_str,
            submitting_user=self.submitting_user, submitting_user_email=self.submitting_user_email)
        self.data_analysis_object.description = self.args.description
        self.data_analysis_object.save()

    def start_data_analysis(self):
        data_sub_collection_run.main(
            data_analysis_object=self.data_analysis_object, num_processors=self.args.num_proc,
            no_figures=self.args.no_figures, no_ordinations=self.args.no_ordinations,
            distance_method=self.args.distance_method, no_output=self.args.no_output,
            debug=self.args.debug)

    # data loading methods
    def perform_data_loading(self):
        self._verify_name_arg_given()
        self.make_new_dataset_object()
        return self._execute_data_loading()

    def _execute_data_loading(self):
        data_loading_object = data_loading.DataLoading(
            data_set_object=self.data_set_object, datasheet_path=self.args.data_sheet, user_input_path=self.args.load,
            screen_sub_evalue=self.screen_sub_eval_bool, num_proc=self.args.num_proc, no_fig=self.args.no_figures,
            no_ord=self.args.no_ordinations, distance_method=self.args.distance_method)
        data_loading_object.load_data()
        return data_loading_object

    def _verify_name_arg_given(self):
        if self.args.name == 'noName':
            sys.exit('Please provide a name using the --name flag. e.g. --name informative_name')

    def make_new_dataset_object(self):
        self.data_set_object = DataSet(
            name=self.args.name, time_stamp=self.date_time_str, reference_fasta_database_used=self.reference_db,
            submitting_user=self.submitting_user, submitting_user_email=self.submitting_user_email)
        self.data_set_object.save()

    # seq output methods
    def perform_sequences_count_table_output(self):
        sequence_count_table_creator = self._execute_seq_output_standalone()

        self._plot_seq_output(sequence_count_table_creator)

    def _plot_seq_output(self, sequence_count_table_creator):
        if len(sequence_count_table_creator.list_of_dss_objects) < 1000:
            seq_stacked_bar_plotter = plotting.SeqStackedBarPlotter(
                output_directory=sequence_count_table_creator.output_dir,
                seq_relative_abund_count_table_path=sequence_count_table_creator.path_to_seq_output_df_relative,
                time_date_str=sequence_count_table_creator.time_date_str)

            seq_stacked_bar_plotter.plot_stacked_bar_seqs()

    def _execute_seq_output_standalone(self):
        if self.args.print_output_seqs:
            sequence_count_table_creator = output.SequenceCountTableCreator(
                symportal_root_dir=self.symportal_root_directory, call_type='stand_alone',
                ds_uids_output_str=self.args.print_output_seqs,
                num_proc=self.args.num_proc, time_date_str=self.date_time_str)
        else:  # self.args.print_output_seqs_sample_set:
            sequence_count_table_creator = output.SequenceCountTableCreator(
                symportal_root_dir=self.symportal_root_directory, call_type='stand_alone',
                ds_uids_output_str=self.args.print_output_seqs_sample_set,
                num_proc=self.args.num_proc, time_date_str=self.date_time_str)
        sequence_count_table_creator.make_output_tables()
        return sequence_count_table_creator

    # type output methods
    def perform_type_count_table_output(self):
        self.verify_data_analysis_uid_provided()

        analysis_object = self.get_data_analysis_by_uid()

        data_sub_collection_run.output_type_count_tables(
            analysisobj=analysis_object, num_processors=self.args.num_proc, call_type='stand_alone',
            datasubstooutput=self.args.print_output_types, no_figures=self.args.no_figures,
            output_user=self.submitting_user)

    def get_data_analysis_by_uid(self):
        return DataAnalysis.objects.get(id=self.data_analysis_object.id)

    def verify_data_analysis_uid_provided(self):
        if not self.args.data_analysis_id:
            raise RuntimeError(
                'Please provide a data_analysis to ouput from by providing a data_analysis '
                'ID to the --data_analysis_id flag. To see a list of data_analysis objects in the '
                'framework\'s database, use the --display_analyses flag.')

    # display db contents functions
    @staticmethod
    def perform_display_data_sets():
        data_set_id_to_obj_dict = {ds.id: ds for ds in list(DataSet.objects.all())}
        sorted_list_of_ids = sorted(list(data_set_id_to_obj_dict.keys()))
        for ds_id in sorted_list_of_ids:
            ds_in_q = data_set_id_to_obj_dict[ds_id]
            print(f'{ds_in_q.id}: {ds_in_q.name}\t{ds_in_q.time_stamp}')

    @staticmethod
    def perform_display_analysis_types():
        data_analysis_id_to_obj_dict = {da.id: da for da in list(DataAnalysis.objects.all())}
        sorted_list_of_ids = sorted(list(data_analysis_id_to_obj_dict.keys()))
        for da_id in sorted_list_of_ids:
            da_in_q = data_analysis_id_to_obj_dict[da_id]
            print(f'{da_in_q.id}: {da_in_q.name}\t{da_in_q.time_stamp}')

    # type distances stand_alone
    def perform_within_clade_type_distance_generation(self):
        self.verify_data_analysis_uid_provided()
        self.run_within_clade_type_distances_dependent_on_methods()
        self._perform_type_distance_plotting()

    def run_within_clade_type_distances_dependent_on_methods(self):
        if self.args.distance_method == 'unifrac':
            self._execute_unifrac_type_dist_pcoa_calc()

        elif self.args.distance_method == 'braycurtis':
            self._execute_braycurtis_type_dist_pcoa_calc()

    def _execute_braycurtis_type_dist_pcoa_calc(self):
        self.pcoa_output_path_list = distance.generate_within_clade_braycurtis_distances_its2_type_profiles(
            data_submission_id_str=self.args.between_type_distances,
            data_analysis_id=self.args.data_analysis_id, call_type='stand_alone',
            date_time_string=self.date_time_str)

    def _execute_unifrac_type_dist_pcoa_calc(self):
        self.pcoa_output_path_list = data_sub_collection_run.generate_within_clade_unifrac_distances_its2_type_profiles(
            data_submission_id_str=self.args.between_type_distances, num_processors=self.args.num_proc,
            data_analysis_id=self.args.data_analysis_id, method='mothur', call_type='stand_alone',
            date_time_string=self.date_time_str, bootstrap_value=self.args.bootstrap)

    def _perform_type_distance_plotting(self):
        for pcoa_path in self.pcoa_output_path_list:
            if 'PCoA_coords' in pcoa_path:
                sys.stdout.write('\nPlotting between its2 type profile distances clade {}\n'.format(
                    os.path.dirname(pcoa_path).split('/')[-1]))
                # then this is a pcoa csv that we should plot
                plotting.plot_between_its2_type_prof_dist_scatter(pcoa_path, date_time_str=self.date_time_str)

    # sample distance stand alone
    def perform_within_clade_sample_distance_generation(self):
        if self.args.distance_method =='unifrac':
            if self.args.between_sample_distances_sample_set:
                raise RuntimeError('DataSetSample distance calculation for custom DataSetSample lists has not yet been '
                                   'implemented using the UniFrac method.\n'
                                   'To output distances between custom DataSetSample lists, please use the braycurtis '
                                   'method.\n'
                                   'e.g: --between_sample_distances_sample_set 1,2,3,4,5,6 --distance_method braycurtis')
            else:
                unifrac_dict_pcoa_creator = self._execute_unifrac_sample_dist_pcoa_calc_ds_uid()
            self._pcoa_samples_plotting(unifrac_dict_pcoa_creator)

        else:
            # braycurtis
            braycurtis_dist_pcoa_creator = self._execute_braycurtis_sample_dist_pcoa_calc()
            self._pcoa_samples_plotting(braycurtis_dist_pcoa_creator)

    def _execute_unifrac_sample_dist_pcoa_calc_ds_uid(self):
        unifrac_dict_pcoa_creator = distance.UnifracDistPCoACreator(
            call_type='stand_alone', date_time_string=self.date_time_str,
            data_set_string=self.args.between_sample_distances, method='mothur', num_processors=self.args.num_proc,
            symportal_root_directory=self.symportal_root_directory)
        unifrac_dict_pcoa_creator.compute_unifrac_dists_and_pcoa_coords()
        return unifrac_dict_pcoa_creator

    def _execute_braycurtis_sample_dist_pcoa_calc(self):
        if self.args.between_sample_distances_sample_set:
            braycurtis_dist_pcoa_creator = distance.BrayCurtisDistPCoACreator(
                symportal_root_directory=self.symportal_root_directory, date_time_string=self.date_time_str,
                data_set_string=self.args.between_sample_distances_sample_set,
                call_type='stand_alone')
        else:  # self.args.between_sample_distances
            braycurtis_dist_pcoa_creator = distance.BrayCurtisDistPCoACreator(
                symportal_root_directory=self.symportal_root_directory, date_time_string=self.date_time_str,
                data_set_string=self.args.between_sample_distances,
                call_type='stand_alone')

        braycurtis_dist_pcoa_creator.compute_braycurtis_dists_and_pcoa_coords()
        return braycurtis_dist_pcoa_creator

    def _pcoa_samples_plotting(self, dist_pcoa_creator):
        for output_path in dist_pcoa_creator.output_file_paths:
            if data_loading.DataLoading.this_is_pcoa_path(output_path):
                clade_of_output = os.path.dirname(output_path).split('/')[-1]
                sys.stdout.write(f'\n\nGenerating between sample distance plot clade {clade_of_output}\n')
                dist_scatter_plotter_samples = plotting.DistScatterPlotterSamples(
                    csv_path=output_path, date_time_str=dist_pcoa_creator.date_time_string)
                dist_scatter_plotter_samples.make_sample_dist_scatter_plot()

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


if __name__ == "__main__":
    spwfm = SymPortalWorkFlowManager()
    spwfm.start_work_flow()
