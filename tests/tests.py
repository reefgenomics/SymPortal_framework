#!/usr/bin/env python3
import os
import main
import django_general
import shutil
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from dbApp.models import DataAnalysis, DataSet
import pandas as pd


class SymPortalTester:
    def __init__(self):
        self.work_flow_manager = None
        self.symportal_testing_root_dir = os.path.abspath(os.path.dirname(__file__))
        self.symportal_root_dir = os.path.abspath(os.path.join(self.symportal_testing_root_dir, '..'))
        self.test_data_dir_path = os.path.join(self.symportal_testing_root_dir, 'data', 'smith_subsampled_data')
        self.data_sheet_file_path = os.path.join(self.test_data_dir_path, 'test_data_submission_input.csv')
        self.num_proc=6
        self.name='testing'
        self.assertion_matching_dir = os.path.join(self.test_data_dir_path, 'assertion_testing')

    def execute_integrated_tests(self):
        self.cleanup_after_previous_tests()
        self._test_data_loading_work_flow()
        self._test_data_analysis_work_flow()
        self.cleanup_after_previous_tests()

    def _test_data_loading_work_flow(self):
        custom_args_list = ['--load', self.test_data_dir_path, '--name', self.name, '--num_proc',
                            str(self.num_proc), '--data_sheet', self.data_sheet_file_path, '--debug']
        self.work_flow_manager = main.SymPortalWorkFlowManager(custom_args_list)
        self.work_flow_manager.start_work_flow()
        # Let's see if we can read in some of the outputs and work with them as assertions
        foo = 'ar'
        # here we should develop some real tests
        # Check the md5sum of the output post_med_seq abund only relative and absolute files
        
        # Post-MED objects to check.
        for path in self.work_flow_manager.data_loading_object.output_path_list:
            if 'seqs.absolute.abund_only' in path:
                absolute_abund_path = path
            if 'seqs.relative.abund_only' in path:
                relative_abund_path = path
            if '.seqs.fasta' in path:
                post_med_fasta_path = path
        
        # Pre-MED objects to check.
        pre_med_absolute_count_table_path = self.work_flow_manager.data_loading_object.sequence_count_table_creator.pre_med_absolute_df_path
        pre_med_relative_count_table_path = self.work_flow_manager.data_loading_object.sequence_count_table_creator.pre_med_relative_df_path
        pre_med_fasta_path = self.work_flow_manager.data_loading_object.sequence_count_table_creator.pre_med_fasta_out_path

        # pre_med absolute
        self._compare_abund_dfs_pre_med(
            assertion_df_path=os.path.join(self.assertion_matching_dir, 'pre_med_absolute_abundance_df.csv'), 
            tested_df_path=pre_med_absolute_count_table_path, absolute=True
            )
        # pre_med relative
        self._compare_abund_dfs_pre_med(
            assertion_df_path=os.path.join(self.assertion_matching_dir, 'pre_med_relative_abundance_df.csv'), 
            tested_df_path=pre_med_relative_count_table_path, absolute=False
            )
        
        # Now compare the DataFrames
        # post_med absolute
        self._compare_abund_dfs_post_med(
            assertion_df_path=os.path.join(self.assertion_matching_dir, 'seqs.absolute.abund_only.txt'), 
            tested_df_path=absolute_abund_path, absolute=True
            )
        # post_med relative
        self._compare_abund_dfs_post_med(
            assertion_df_path=os.path.join(self.assertion_matching_dir, 'seqs.relative.abund_only.txt'), 
            tested_df_path=relative_abund_path, absolute=False
            )
        
    # NB MED produces inconsistent results which means that we have to be careful with how we check the results.
    # It means that the order of the list of sequences and a very small number of the abundances will be different
    def _compare_abund_dfs_post_med(self, assertion_df_path, tested_df_path, absolute):
        a_df = pd.read_table(assertion_df_path, index_col=0)
        t_df = pd.read_table(tested_df_path, index_col=0)
        
        a_tot_abund = a_df.sum().sum()
        t_tot_abund = t_df.sum().sum()
        if a_tot_abund > t_tot_abund:
            if t_tot_abund / a_tot_abund < 0.01:
                raise AssertionError('Sequence abundances do not match for post_med count tables')
        else:
            if a_tot_abund / t_tot_abund < 0.01:
                raise AssertionError('Sequence abundances do not match for post_med count tables')
    
    def _compare_abund_dfs_pre_med(self, assertion_df_path, tested_df_path, absolute=True):
        # We multiply by 1000 and convert to int so that floats can be compared.
        a_df = pd.read_csv(assertion_df_path, index_col=0).drop(columns=['sample_name'])
        t_df = pd.read_csv(tested_df_path, index_col=0).drop(columns=['sample_name'])
        
        # What the sequences are called may differ between databases and even subtly how many sequences are allocated
        # easiest just
        # Check that the overall sums are within 1 % difference.
        a_tot_abund = a_df.sum().sum()
        t_tot_abund = t_df.sum().sum()
        if a_tot_abund > t_tot_abund:
            if t_tot_abund / a_tot_abund < 0.01:
                raise AssertionError('Sequence abundances do not match for pre_med count tables')
        else:
            if a_tot_abund / t_tot_abund < 0.01:
                raise AssertionError('Sequence abundances do not match for pre_med count tables')

    def _test_data_loading_work_flow_no_datasheet(self):
        custom_args_list = ['--load', self.test_data_dir_path, '--name', self.name, '--num_proc',
                            str(self.num_proc), '--debug']
        self.work_flow_manager = main.SymPortalWorkFlowManager(custom_args_list)
        self.work_flow_manager.start_work_flow()

    def _test_data_analysis_work_flow(self):
        custom_args_list = ['--analyse', str(self.work_flow_manager.data_set_object.id), '--name', self.name,
                            '--num_proc', str(self.num_proc)]
        self.work_flow_manager = main.SymPortalWorkFlowManager(custom_args_list)
        self.work_flow_manager.start_work_flow()
        #TODO testing of the profile count table and seq count tables again.

    def cleanup_after_previous_tests(self):
        self.cleanup_after_previously_run_data_analysis_tests()
        self.cleanup_after_previously_run_data_loading_tests()

    def cleanup_after_previously_run_data_loading_tests(self):
        uids_of_testing_data_set_objects = (ds.id for ds in DataSet.objects.filter(name='testing'))
        for ds_uid in uids_of_testing_data_set_objects:
            print(f'Cleaning up after previous data loading test: {ds_uid}')
            self.delete_data_load_output_object_directories(ds_uid)
            django_general.delete_data_set(ds_uid)

    def delete_data_load_output_object_directories(self, data_set_uid):
        directory_to_delete = os.path.abspath(os.path.join(
            self.symportal_root_dir, 'outputs', 'loaded_data_sets', str(data_set_uid)))
        if os.path.exists(directory_to_delete):
            print(f'Deleting {directory_to_delete}')
            shutil.rmtree(directory_to_delete)

    def cleanup_after_previously_run_data_analysis_tests(self):
        uids_of_testing_data_analysis_objects = (da.id for da in DataAnalysis.objects.filter(name='testing'))
        for da_uid in uids_of_testing_data_analysis_objects:
            print(f'Cleaning up after previous data analysis test: {da_uid}')
            self.delete_data_analysis_output_object_directories(da_uid)
            django_general.delete_data_analysis(da_uid)

    def delete_data_analysis_output_object_directories(self, data_analysis_uid):
        directory_to_delete = os.path.abspath(os.path.join(
            self.symportal_root_dir, 'outputs', 'analyses', str(data_analysis_uid)))
        if os.path.exists(directory_to_delete):
            print(f'Deleting {directory_to_delete}')
            shutil.rmtree(directory_to_delete)

if __name__ == "__main__":
    symportal_tester = SymPortalTester()
    symportal_tester.execute_integrated_tests()
