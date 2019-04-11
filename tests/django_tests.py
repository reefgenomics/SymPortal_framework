#!/usr/bin/env python3.6
from django.test import TransactionTestCase
import os
import main
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from dbApp.models import DataSet, DataSetSample


class SPIntegrativeTestingJSONOnly(TransactionTestCase):
    """This class of tests contains tests that can all be run without any pickled binary files.
    The fixture that is loaded contains three datasets and one analysis.
    The analysis is of all three datasets.

    The DataSet IDs are 173, 174 and 175.
    The DataAnalysis ID is 1.

    NB I was originally running this as a class of TestCase rather than TransactionTestCase, however, when running as
    TestCase the database was being created in memory and so new connections could not be closed and reopened. With
    TransactionTestCase the database is written and so connections can be opened and closed no problem.
    """
    fixtures = ['three_dataset_one_analysis_db.json']

    @classmethod
    def setUp(cls):
        cls.symportal_testing_root_dir = os.path.abspath(os.path.dirname(__file__))
        cls.symportal_root_dir = os.path.abspath(os.path.join(cls.symportal_testing_root_dir, '..'))
        cls.test_data_dir_path = os.path.join(cls.symportal_testing_root_dir, 'data', 'smith_subsampled_data')
        cls.data_sheet_file_path = os.path.join(cls.test_data_dir_path, 'test_data_submission_input.csv')
        cls.num_proc = 6
        cls.name = 'testing'

    # TEST DATA LOADING
    def test_data_loading_work_flow(self):
        custom_args_list = ['--load', self.test_data_dir_path, '--name', self.name, '--num_proc', str(self.num_proc),
                            '--data_sheet', self.data_sheet_file_path]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    # TEST DATA ANALYSIS
    def test_data_analysis_work_flow(self):
        custom_args_list = ['--analyse', '173', '--name', self.name,
                            '--num_proc', str(self.num_proc)]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    # TEST STAND_ALONE OUTPUTS
    def test_stand_alone_sequence_outputs_data_set_input(self):
        custom_args_list = ['--print_output_seqs', '173,174', '--num_proc', str(self.num_proc)]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_sequence_outputs_data_set_sample_input(self):
        data_set_samples_173_174 = DataSet.objects.filter(id__in=[173,174])
        data_set_samples_in_173_174 = DataSetSample.objects.filter(data_submission_from__in=data_set_samples_173_174)
        data_set_sample_uids_str = ','.join([str(dss.id) for dss in data_set_samples_in_173_174][2:-2])
        custom_args_list = ['--print_output_seqs_sample_set', data_set_sample_uids_str, '--num_proc', str(self.num_proc)]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_profile_outputs_data_set_input(self):
        custom_args_list = ['--print_output_types', '173', '--data_analysis_id', '1', '--num_proc', str(self.num_proc)]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_profile_outputs_data_set_sample_input(self):
        data_set_sample_173 = DataSet.objects.get(id=173)
        data_set_samples_in_173 = DataSetSample.objects.filter(data_submission_from=data_set_sample_173)
        data_set_sample_uids_str = ','.join([str(dss.id) for dss in data_set_samples_in_173][2:-2])
        custom_args_list = ['--print_output_types_sample_set', data_set_sample_uids_str, '--data_analysis_id', '1', '--num_proc', str(self.num_proc)]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    # TEST STAND_ALONE DISTANCES
    def test_stand_alone_unifrac_type_distances_data_set_uid_input(self):
        custom_args_list = ['--between_type_distances', '173', '--data_analysis_id', '1', '--num_proc', str(self.num_proc), '--distance_method', 'unifrac']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_unifrac_type_distances_data_set_sample_uid_input(self):
        # create a string for dss input that is all but the last 4 dss objects of the DataSet 173
        data_set_173 = DataSet.objects.get(id=173)
        data_set_samples_from_173 = DataSetSample.objects.filter(data_submission_from=data_set_173)
        dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_173][:-4])
        custom_args_list = ['--between_type_distances_sample_set', dss_uids_of_ds_str, '--data_analysis_id', '1', '--num_proc', str(self.num_proc), '--distance_method', 'unifrac']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_braycurtis_type_distances_data_set_uid_input(self):
        custom_args_list = ['--between_type_distances', '173', '--data_analysis_id', '1', '--num_proc', str(self.num_proc), '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_braycurtis_type_distances_data_set_sample_uid_input(self):
        # create a string for dss input that is all but the last 4 dss objects of the DataSet 173
        data_set_173 = DataSet.objects.get(id=173)
        data_set_samples_from_173 = DataSetSample.objects.filter(data_submission_from=data_set_173)
        dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_173][:-4])
        custom_args_list = ['--between_type_distances_sample_set', dss_uids_of_ds_str, '--data_analysis_id', '1', '--num_proc', str(self.num_proc), '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_unifrac_sample_distances_data_set_uid_input(self):
        custom_args_list = ['--between_sample_distances', '173', '--num_proc', str(self.num_proc), '--distance_method', 'unifrac']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_unifrac_sample_distances_data_set_sample_uid_input(self):
        # create a string for dss input that is all but the last 4 dss objects of the DataSet 173
        data_set_173 = DataSet.objects.get(id=173)
        data_set_samples_from_173 = DataSetSample.objects.filter(data_submission_from=data_set_173)
        dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_173][:-4])
        custom_args_list = ['--between_sample_distances_sample_set', dss_uids_of_ds_str, '--num_proc', str(self.num_proc), '--distance_method', 'unifrac']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_braycurtis_sample_distances_data_set_uid_input(self):
        custom_args_list = ['--between_sample_distances', '173', '--num_proc', str(self.num_proc), '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_braycurtis_sample_distances_data_set_sample_uid_input(self):
        # create a string for dss input that is all but the last 4 dss objects of the DataSet 173
        data_set_173 = DataSet.objects.get(id=173)
        data_set_samples_from_173 = DataSetSample.objects.filter(data_submission_from=data_set_173)
        dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_173][:-4])
        custom_args_list = ['--between_sample_distances_sample_set', dss_uids_of_ds_str, '--num_proc', str(self.num_proc), '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()