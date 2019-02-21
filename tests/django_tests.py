#!/usr/bin/env python3.6
from django.test import TransactionTestCase
import os
import main
import shutil

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from dbApp.models import DataSet, DataAnalysis


class SPIntegrativeTesting(TransactionTestCase):
    """This fixture is a dump of a database with contents specifically setup for testing. It contains
    three loaded DataSet objects with all of the corresponding features and a single analysis that contains
    all three of the DataSets.

    The DataSet IDs are:
    1: Smith et al
    7: Anon
    27: Tullia

    The DataAnalysis ID is:
    1: contains DataSets 1, 7 and 27

    NB I was originally running this as a class of TestCase rather than TransactionTestCase, however, when running as
    TestCase the database was being created in memory and so new connections could not be closed and reopened. With
    TransactionTestCase the database is written and so connections can be opened and closed no problem.

    """
    fixtures = ['testing_fixture.json']

    @classmethod
    def setUp(cls):
        cls.symportal_testing_root_dir = os.path.abspath(os.path.dirname(__file__))
        cls.symportal_root_dir = os.path.abspath(os.path.join(cls.symportal_testing_root_dir, '..'))
        cls.test_data_dir_path = os.path.join(cls.symportal_testing_root_dir, 'data', 'smith_subsampled_data')
        cls.data_sheet_file_path = os.path.join(cls.test_data_dir_path, 'test_data_submission_input.csv')
        cls.num_proc = 6
        cls.name = 'testing'

    def test_data_loading(self):
        self._cleanup_after_previously_run_data_loading_tests()
        self._start_data_loading_work_flow()
        self._cleanup_after_previously_run_data_loading_tests()

    def _start_data_loading_work_flow(self):
        custom_args_list = ['--load', self.test_data_dir_path, '--name', self.name, '--num_proc', str(self.num_proc),
                            '--data_sheet',
                            self.data_sheet_file_path]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def _cleanup_after_previously_run_data_loading_tests(self):
        uids_of_testing_data_set_objects = (ds.id for ds in DataSet.objects.filter(name='testing'))
        for ds_uid in uids_of_testing_data_set_objects:
            print(f'Cleaning up after previous data loading test: {ds_uid}')
            self._delete_data_load_output_object_directories(ds_uid)

    def _delete_data_load_output_object_directories(self, data_set_uid):
        directory_to_delete = os.path.abspath(os.path.join(
            self.symportal_root_dir, 'outputs', 'data_set_submissions', str(data_set_uid)))
        if os.path.exists(directory_to_delete):
            shutil.rmtree(directory_to_delete)

    # def test_data_analysis(self):
    #     self._cleanup_after_previously_run_data_analysis_tests()
    #     self._start_data_analysis_work_flow()
    #     self._cleanup_after_previously_run_data_analysis_tests()

    def _start_data_analysis_work_flow(self):
        custom_args_list = ['--analyse', '1,7,27', '--name', self.name,
                            '--num_proc', str(self.num_proc)]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def _cleanup_after_previously_run_data_analysis_tests(self):
        uids_of_testing_data_analysis_objects = (da.id for da in DataAnalysis.objects.filter(name='testing'))
        for da_uid in uids_of_testing_data_analysis_objects:
            print(f'Cleaning up after previous data analysis test: {da_uid}')
            self._delete_data_analysis_output_object_directories(da_uid)

    def _delete_data_analysis_output_object_directories(self, data_analysis_uid):
        directory_to_delete = os.path.abspath(os.path.join(
            self.symportal_root_dir, 'outputs', 'analyses', str(data_analysis_uid)))
        if os.path.exists(directory_to_delete):
            shutil.rmtree(directory_to_delete)