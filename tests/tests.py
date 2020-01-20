#!/usr/bin/env python3.6
import os
import main
import general
import django_general
from pathlib import Path
import shutil
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from dbApp.models import DataAnalysis, DataSet


class SymPortalTester:
    def __init__(self):
        self.work_flow_manager = None
        self.symportal_testing_root_dir = os.path.abspath(os.path.dirname(__file__))
        self.symportal_root_dir = os.path.abspath(os.path.join(self.symportal_testing_root_dir, '..'))
        self.test_data_dir_path = os.path.join(self.symportal_testing_root_dir, 'data', 'smith_subsampled_data')
        self.data_sheet_file_path = os.path.join(self.test_data_dir_path, 'test_data_submission_input.csv')
        self.num_proc=6
        self.name='testing'

    def execute_integrated_tests(self):
        self.cleanup_after_previous_tests()
        self._test_data_loading_work_flow()
        self._test_data_analysis_work_flow()
        self.cleanup_after_previous_tests()

    def _test_data_loading_work_flow(self):
        custom_args_list = ['--load', self.test_data_dir_path, '--name', self.name, '--num_proc', str(self.num_proc), '--data_sheet',
             self.data_sheet_file_path, '--debug']
        self.work_flow_manager = main.SymPortalWorkFlowManager(custom_args_list)
        self.work_flow_manager.start_work_flow()

    def _test_data_loading_work_flow_no_datasheet(self):
        custom_args_list = ['--load', self.test_data_dir_path, '--name', self.name, '--num_proc', str(self.num_proc), '--debug']
        self.work_flow_manager = main.SymPortalWorkFlowManager(custom_args_list)
        self.work_flow_manager.start_work_flow()

    def _test_data_analysis_work_flow(self):
        custom_args_list = ['--analyse', str(self.work_flow_manager.data_set_object.id), '--name', self.name,
                            '--num_proc', str(self.num_proc)]
        self.work_flow_manager = main.SymPortalWorkFlowManager(custom_args_list)
        self.work_flow_manager.start_work_flow()

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
