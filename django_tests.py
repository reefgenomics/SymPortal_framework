#!/usr/bin/env python3.6
from django.test import TransactionTestCase
import os
import main
import shutil
from data_analysis import ArtefactAssessor
import pickle
import plotting
import distance


os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from dbApp.models import DataSet, DataAnalysis, DataSetSample


class SPIntegrativeTestingJSONOnly(TransactionTestCase):
    """This class of tests contains tests that can all be run without any pickled binary files.
    The fixture that is loaded contains three datasets and one analysis..

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
    # def test_data_loading_work_flow(self):
    #     custom_args_list = ['--load', self.test_data_dir_path, '--name', self.name, '--num_proc', str(self.num_proc),
    #                         '--data_sheet', self.data_sheet_file_path]
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()

    # TEST DATA ANALYSIS
    # def test_data_analysis_work_flow(self):
    #     custom_args_list = ['--analyse', '173', '--name', self.name,
    #                         '--num_proc', str(self.num_proc)]
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()

    # TEST STAND_ALONE OUTPUTS
    # def test_stand_alone_sequence_outputs_data_set_input(self):
    #     custom_args_list = ['--print_output_seqs', '173,174', '--num_proc', str(self.num_proc)]
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()

    # def test_stand_alone_sequence_outputs_data_set_sample_input(self):
    #     data_set_samples_173_174 = DataSet.objects.filter(id__in=[173,174])
    #     data_set_samples_in_173_174 = DataSetSample.objects.filter(data_submission_from__in=data_set_samples_173_174)
    #     data_set_sample_uids_str = ','.join([str(dss.id) for dss in data_set_samples_in_173_174][2:-2])
    #     custom_args_list = ['--print_output_seqs_sample_set', data_set_sample_uids_str, '--num_proc', str(self.num_proc)]
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()

    def test_stand_alone_profile_outputs_data_set_input(self):
        custom_args_list = ['--print_output_types', '173', '--data_analysis_id', '1', '--num_proc', str(self.num_proc)]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()
    #
    # def test_stand_alone_profile_outputs_data_set_sample_input(self):
    #     data_set_sample_173 = DataSet.objects.get(id=173)
    #     data_set_samples_in_173 = DataSetSample.objects.filter(data_submission_from=data_set_sample_173)
    #     data_set_sample_uids_str = ','.join([str(dss.id) for dss in data_set_samples_in_173][2:-2])
    #     custom_args_list = ['--print_output_types_sample_set', data_set_sample_uids_str, '--data_analysis_id', '1', '--num_proc', str(self.num_proc)]
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()

    # TEST STAND_ALONE DISTANCES
    # def test_stand_alone_unifrac_type_distances_data_set_uid_input(self):
    #     custom_args_list = ['--between_type_distances', '173', '--data_analysis_id', '1', '--num_proc', str(self.num_proc), '--distance_method', 'unifrac']
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()
    #
    # def test_stand_alone_unifrac_type_distances_data_set_sample_uid_input(self):
    #     # create a string for dss input that is all but the last 4 dss objects of the DataSet 173
    #     data_set_173 = DataSet.objects.get(id=173)
    #     data_set_samples_from_173 = DataSetSample.objects.filter(data_submission_from=data_set_173)
    #     dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_173][:-4])
    #     custom_args_list = ['--between_type_distances_sample_set', dss_uids_of_ds_str, '--data_analysis_id', '1', '--num_proc', str(self.num_proc), '--distance_method', 'unifrac']
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()

    # def test_stand_alone_braycurtis_type_distances_data_set_uid_input(self):
    #     custom_args_list = ['--between_type_distances', '173', '--data_analysis_id', '1', '--num_proc', str(self.num_proc), '--distance_method', 'braycurtis']
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()

    # def test_stand_alone_braycurtis_type_distances_data_set_sample_uid_input(self):
    #     # create a string for dss input that is all but the last 4 dss objects of the DataSet 173
    #     data_set_173 = DataSet.objects.get(id=173)
    #     data_set_samples_from_173 = DataSetSample.objects.filter(data_submission_from=data_set_173)
    #     dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_173][:-4])
    #     custom_args_list = ['--between_type_distances_sample_set', dss_uids_of_ds_str, '--data_analysis_id', '1', '--num_proc', str(self.num_proc), '--distance_method', 'braycurtis']
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()

    # def test_stand_alone_unifrac_sample_distances_data_set_uid_input(self):
    #     custom_args_list = ['--between_sample_distances', '173', '--num_proc', str(self.num_proc), '--distance_method', 'unifrac']
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()
    #
    # def test_stand_alone_unifrac_sample_distances_data_set_sample_uid_input(self):
    #     # create a string for dss input that is all but the last 4 dss objects of the DataSet 173
    #     data_set_173 = DataSet.objects.get(id=173)
    #     data_set_samples_from_173 = DataSetSample.objects.filter(data_submission_from=data_set_173)
    #     dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_173][:-4])
    #     custom_args_list = ['--between_sample_distances_sample_set', dss_uids_of_ds_str, '--num_proc', str(self.num_proc), '--distance_method', 'unifrac']
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()
    #
    # def test_stand_alone_braycurtis_sample_distances_data_set_uid_input(self):
    #     custom_args_list = ['--between_sample_distances', '173', '--num_proc', str(self.num_proc), '--distance_method', 'braycurtis']
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()
    #
    # def test_stand_alone_braycurtis_sample_distances_data_set_sample_uid_input(self):
    #     # create a string for dss input that is all but the last 4 dss objects of the DataSet 173
    #     data_set_173 = DataSet.objects.get(id=173)
    #     data_set_samples_from_173 = DataSetSample.objects.filter(data_submission_from=data_set_173)
    #     dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_173][:-4])
    #     custom_args_list = ['--between_sample_distances_sample_set', dss_uids_of_ds_str, '--num_proc', str(self.num_proc), '--distance_method', 'braycurtis']
    #     test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
    #     test_spwfm.start_work_flow()




class SPIntegrativeTestingBinaries(TransactionTestCase):
    """
    Contains tests that are reliant on pickled objects that are not part of the GitHub repository.
    """

    fixtures = ['three_dataset_one_analysis_db.json']

    @classmethod
    def setUp(cls):
        cls.symportal_testing_root_dir = os.path.abspath(os.path.dirname(__file__))
        cls.symportal_root_dir = os.path.abspath(os.path.join(cls.symportal_testing_root_dir, '..'))
        cls.artefact_assessor_within_clade_cutoff_pickled_object_path = os.path.join(cls.symportal_root_dir, 'tests', 'objects', 'artefact_assessor.p')
        cls.artefact_assessor_reassess_support_pickled_object_path = os.path.join(cls.symportal_root_dir, 'tests', 'objects', 'artefact_assessor_reassess_support.p')
        cls.sp_data_analysis_post_artefact_assess_pickled_object_path = os.path.join(cls.symportal_root_dir, 'tests', 'objects', 'sp_data_analysis_post_artefact_assess.p')
        cls.sp_data_analysis_post_profile_assignment_pickled_object_path = os.path.join(cls.symportal_root_dir, 'tests', 'objects', 'sp_data_analysis_post_profile_assignment.p')
        cls.sp_data_analysis_post_div_naming_pickled_object_path = os.path.join(cls.symportal_root_dir, 'tests', 'objects', 'sp_data_analysis_post_div_naming.p')
        cls.sp_work_flow_post_analysis_pickled_object_path = os.path.join(cls.symportal_root_dir, 'tests', 'objects', 'sp_workflow_post_analysis.p')
        cls.sp_data_analysis_post_del_tempdir_pickled_object_path = os.path.join(cls.symportal_root_dir, 'tests', 'objects',
                                                                          'sp_data_analysis_post_del_tempdir.p')
        cls.sp_data_analysis_post_type_output_pickled_object_path = os.path.join(cls.symportal_root_dir, 'tests',
                                                                                 'objects',
                                                                                 'sp_workflow_post_type_output.p')

        cls.num_proc = 6
        cls.name = 'testing'


    # def test_artefact_assesor_assess_within_clade_cutoff_artefacts(self):
    #     with open(self.artefact_assessor_within_clade_cutoff_pickled_object_path, 'rb') as f:
    #         artefact_assessor = pickle.load(f)
    #         artefact_assessor.assess_within_clade_cutoff_artefacts()
    #         artefact_assessor.reassess_support_of_artefact_div_containing_types()
    #         apples = 'asdf'
    #
    # def test_profile_assignment(self):
    #     with open(self.sp_data_analysis_post_artefact_assess_pickled_object_path, 'rb') as f:
    #         sp_data_analysis = pickle.load(f)
    #         sp_data_analysis._profile_assignment()
    #
    # def test_div_naming(self):
    #     with open(self.sp_data_analysis_post_profile_assignment_pickled_object_path, 'rb') as f:
    #         sp_data_analysis = pickle.load(f)
    #         sp_data_analysis._name_divs()
    #
    # def test_assign_speceis(self):
    #     with open(self.sp_data_analysis_post_div_naming_pickled_object_path, 'rb') as f:
    #         sp_data_analysis = pickle.load(f)
    #         sp_data_analysis._associate_species_designations()
    #
    # def test_make_analysis_type_objects_from_vats(self):
    #     with open(self.sp_data_analysis_post_del_tempdir_pickled_object_path, 'rb') as f:
    #         sp_data_analysis = pickle.load(f)
    #         sp_data_analysis._make_analysis_type_objects_from_vats()
    #
    # def test_analysis_output(self):
    #     with open(self.sp_work_flow_post_analysis_pickled_object_path, 'rb') as f:
    #         sp_workflow = pickle.load(f)
    #         sp_workflow._perform_data_analysis_output_type_tables()
    #         sp_workflow._perform_data_analysis_output_seq_tables()
    #
    # def test_analysis_stacked_bar_plotting(self):
    #     tsbp = plotting.TypeStackedBarPlotter(
    #         output_directory='/Users/humebc/Documents/SymPortal_testing_repo/SymPortal_framework/outputs/analyses/1',
    #         type_relative_abund_count_table_path='/Users/humebc/Documents/SymPortal_testing_repo/SymPortal_framework/outputs/analyses/1/1_testing_2019-03-17_03-33-43.379932.profiles.relative.txt',
    #         time_date_str='2019-03-17_03-33-43.379932')
    #     tsbp.plot_stacked_bar_seqs()
    #
    # def test_perform_data_analysis_plot_seq_bar(self):
    #     with open(os.path.join(self.symportal_root_dir, 'tests', 'objects', 'sample_uid_ordered.p'), 'rb') as f:
    #         sorted_list_of_vdss_uids_to_output = pickle.load(f)
    #     self.seq_stacked_bar_plotter = plotting.SeqStackedBarPlotter(
    #         output_directory='/Users/humebc/Documents/SymPortal_testing_repo/SymPortal_framework/outputs/analyses/1',
    #         seq_relative_abund_count_table_path='/Users/humebc/Documents/SymPortal_testing_repo/SymPortal_framework/outputs/analyses/1/1_testing_2019-03-18_05-11-29.237882.seqs.relative.txt',
    #         ordered_sample_uid_list=sorted_list_of_vdss_uids_to_output,
    #         time_date_str='2019-03-18_05-11-29.237882')
    #     self.seq_stacked_bar_plotter.plot_stacked_bar_seqs()
    #
    # def test_ordination_of_data_sets(self):
    #     bsudpc = distance_dev.BtwnSampleUnifracDistPCoACreator(
    #         call_type='stand_alone', date_time_string=None, method='mothur', num_processors=self.num_proc,
    #         symportal_root_directory=self.symportal_root_dir, data_set_uid_list=[173])
    #     bsudpc.compute_unifrac_dists_and_pcoa_coords()

