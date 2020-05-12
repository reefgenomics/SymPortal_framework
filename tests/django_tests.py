#!/usr/bin/env python3
from django.test import TransactionTestCase
import os
import main
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from dbApp.models import DataSet, DataSetSample, DataAnalysis, CladeCollectionType


class SPIntegrativeTestingJSONOnly(TransactionTestCase):
    """This class of tests contains tests that can all be run without any pickled binary files.
    The fixture that is loaded contains three datasets and one analysis.
    The analysis is of all three datasets.

    The DataSet IDs are 1, 2 and 3.
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
        cls.test_data_dir_path_lite = os.path.join(cls.test_data_dir_path, 'lite')
        cls.data_sheet_file_path_lite = os.path.join(cls.test_data_dir_path_lite, 'test_data_submission_input_lite.csv')
        cls.data_sheet_file_path = os.path.join(cls.test_data_dir_path, 'test_data_submission_input.csv')
        cls.num_proc = 6
        cls.name = 'testing'

    # TEST DATA LOADING
    def test_data_loading_work_flow(self):
        print('\n\nTesting: data_loading_work_flow\n\n')
        custom_args_list = ['--load', self.test_data_dir_path, '--name', self.name, '--num_proc', str(self.num_proc),
                            '--data_sheet', self.data_sheet_file_path]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_data_loading_work_flow_no_pre_med_seqs(self):
        print('\n\nTesting: data_loading_work_flow_no_pre_med_seqs\n\n')
        custom_args_list = ['--load', self.test_data_dir_path_lite, '--name', self.name, '--num_proc',
                            str(self.num_proc), '--data_sheet', self.data_sheet_file_path_lite, '--no_pre_med_seqs',
                            '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()
        custom_args_list = ['--analyse', f'{str(test_spwfm.data_set_object.id)}', '--name', self.name,
                            '--num_proc', str(self.num_proc)]
        test_spwfm_0 = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm_0.start_work_flow()

    def test_data_loading_work_flow_no_output(self):
        print('\n\nTesting: data_loading_work_flow_no_output\n\n')
        custom_args_list = ['--load', self.test_data_dir_path, '--name', self.name, '--num_proc', str(self.num_proc),
                            '--data_sheet', self.data_sheet_file_path, '--no_output', '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_data_loading_work_flow_no_datasheet_zero_byte(self):
        print('\n\nTesting: data_loading_work_flow_no_datasheet\n\n')
        custom_args_list = ['--load', self.test_data_dir_path_lite, '--name', self.name, '--num_proc',
                            str(self.num_proc), '--no_output', '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_data_loading_work_flow_w_datasheet_zero_byte(self):
        print('\n\nTesting: data_loading_work_flow_no_datasheet\n\n')
        custom_args_list = ['--load', self.test_data_dir_path_lite, '--name', self.name, '--data_sheet',
                            self.data_sheet_file_path_lite, '--num_proc', str(self.num_proc), '--no_output',
                            '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    # TEST DATA LOADING only one sample with numerical name
    def test_data_loading_work_flow_number_only(self):
        print('\n\nTesting: data_loading_work_flow_number_only\n\n')
        test_data_dir_path_num_only = os.path.join(self.test_data_dir_path, 'num_only')
        data_sheet_file_path = os.path.join(test_data_dir_path_num_only, 'test_data_submission_input_lite_num_only.csv')
        custom_args_list = ['--load', test_data_dir_path_num_only, '--name', self.name,
                            '--num_proc', str(self.num_proc),
                            '--data_sheet', data_sheet_file_path, '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    # TEST DATA LOADING using a datasheet that contains absolute paths
    def test_data_loading_work_flow_datasheet_w_abs_paths(self):
        """For this test we will use the 'lite' dataset and include the A02 sample
        that has its fwd and rev reads in two seperate directories. We will also use a datasheet
        that has the full paths for both of these reads in it. We have to make this datasheet at run time
        as we don't know the absolute directory of this test dir on the users machine."""
        print('\n\nTesting: data_loading_work_flow_datasheet_w_abs_paths\n\n')
        # Read in the test_data_submission_input_list.csv datasheet and add the A02 sample
        with open(os.path.join(self.test_data_dir_path_lite, 'test_data_submission_input_lite.csv'), 'r') as f:
            lines = [line.split(',') for line in f]
        # add the new sample with the absolute directories in the datasheet
        lines.append(
            ['A02',
             os.path.join(self.test_data_dir_path_lite, 'test_sub_dir_one', 'A02.1_subsampled.fastq'),
             os.path.join(self.test_data_dir_path_lite, 'test_sub_dir_two', 'A02.2_subsampled.fastq'),
             'coral', 'cnidaria', 'anthozoa', 'scleractinia', 'poritidae',
             'porites', 'lobata', '24.339216', '53.059163', '18.06.86', '0-12'])
        new_datasheet_path = os.path.join(self.test_data_dir_path_lite,
                                          'test_data_submission_input_lite_w_abs_paths.csv')
        # write out the new csv
        with open(new_datasheet_path, 'w') as f:
            for line_list in lines:
                csv_line = ",".join(line_list)
                f.write(f'{csv_line}\n')

        custom_args_list = ['--load', self.test_data_dir_path_lite, '--name', self.name, '--num_proc',
                            str(self.num_proc),
                            '--data_sheet', new_datasheet_path, '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    # TEST ANNOTATION OF DATASET WITH DATASHEET
    def test_annotation_of_dataset_with_data_sheet_good(self):
        print('\n\nTesting: test_annotation_of_dataset_with_data_sheet_good')
        dss_uid = 1
        path_to_datasheet = os.path.join(
            self.symportal_testing_root_dir,
            'test_datasheet_for_standalone_application_testing.csv')
        custom_args_list = ['--apply_data_sheet', str(dss_uid), '--data_sheet', path_to_datasheet]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_annotation_of_dataset_with_data_sheet_bad_lat(self):
        print('\n\nTesting: test_annotation_of_dataset_with_data_sheet_bad_lat')
        dss_uid = 1
        path_to_datasheet = os.path.join(
            self.symportal_testing_root_dir,
            'test_datasheet_for_standalone_application_testing_bad_lat.csv')
        custom_args_list = ['--apply_data_sheet', str(dss_uid), '--data_sheet', path_to_datasheet]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_annotation_of_dataset_with_data_sheet_blanks(self):
        print('\n\nTesting: test_annotation_of_dataset_with_data_sheet_blanks')
        dss_uid = 1
        path_to_datasheet = os.path.join(
            self.symportal_testing_root_dir,
            'test_datasheet_for_standalone_application_testing_blanks.csv')
        custom_args_list = ['--apply_data_sheet', str(dss_uid), '--data_sheet', path_to_datasheet]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_annotation_of_dataset_with_data_sheet_wrong_name(self):
        print('\n\nTesting: test_annotation_of_dataset_with_data_sheet_bad_name')
        dss_uid = 1
        path_to_datasheet = os.path.join(
            self.symportal_testing_root_dir,
            'test_datasheet_for_standalone_application_testing_wrong_name.csv')
        custom_args_list = ['--apply_data_sheet', str(dss_uid), '--data_sheet', path_to_datasheet]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    # TEST DATA ANALYSIS
    def test_data_analysis_work_flow(self):
        print('\n\nTesting: data_analysis_work_flow\n\n')
        custom_args_list = ['--analyse', '1', '--name', self.name,
                            '--num_proc', str(self.num_proc)]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_data_analysis_work_flow_force_basal_lineage_separation(self):
        print('\n\nTesting: data_analysis_work_flow_force_basal_lineage_separation\n\n')
        custom_args_list = ['--analyse', '1', '--name', self.name,
                            '--num_proc', str(self.num_proc), '--force_basal_lineage_separation', '--no_output']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()


    def test_data_analysis_work_flow_no_output(self):
        print('\n\nTesting: data_analysis_work_flow\n\n')
        custom_args_list = ['--analyse', '1', '--name', self.name,
                            '--num_proc', str(self.num_proc), '--no_output', '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    # TEST STAND_ALONE OUTPUTS
    def test_stand_alone_sequence_outputs_data_set_input(self):
        print('\n\nTesting: stand_alone_sequence_outputs_data_set_input\n\n')
        custom_args_list = ['--print_output_seqs', '1,2', '--num_proc', str(self.num_proc)]
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_sequence_outputs_data_set_sample_input(self):
        print('\n\nTesting: stand_alone_sequence_outputs_data_set_sample_input\n\n')
        data_set_samples_1_2 = DataSet.objects.filter(id__in=[1, 2])
        data_set_samples_in_1_2 = DataSetSample.objects.filter(data_submission_from__in=data_set_samples_1_2)
        data_set_sample_uids_str = ','.join([str(dss.id) for dss in data_set_samples_in_1_2][2:-2])
        custom_args_list = ['--print_output_seqs_sample_set', data_set_sample_uids_str, '--num_proc',
                            str(self.num_proc), '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_profile_outputs_data_set_input(self):
        print('\n\nTesting: stand_alone_profile_outputs_data_set_input\n\n')
        custom_args_list = ['--print_output_types', '1', '--data_analysis_id', '1', '--num_proc', str(self.num_proc),
                            '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_profile_outputs_data_set_sample_input(self):
        print('\n\nTesting: stand_alone_profile_outputs_data_set_sample_input\n\n')
        data_set_sample_1 = DataSet.objects.get(id=1)
        data_set_samples_in_1 = DataSetSample.objects.filter(data_submission_from=data_set_sample_1)
        data_set_sample_uids_str = ','.join([str(dss.id) for dss in data_set_samples_in_1][2:-2])
        custom_args_list = ['--print_output_types_sample_set', data_set_sample_uids_str,
                            '--data_analysis_id', '1', '--num_proc', str(self.num_proc),
                            '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    # TEST STAND_ALONE DISTANCES
    def test_stand_alone_unifrac_type_distances_data_set_uid_input(self):
        print('\n\nTesting: stand_alone_unifrac_type_distances_data_set_uid_input\n\n')
        custom_args_list = ['--between_type_distances', '1', '--data_analysis_id', '1', '--num_proc',
                            str(self.num_proc), '--distance_method', 'unifrac', '--local']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_unifrac_type_distances_data_set_sample_uid_input(self):
        # create a string for dss input that is all but the last 4 dss objects of the DataSet 1
        print('\n\nTesting: stand_alone_unifrac_type_distances_data_set_sample_uid_input\n\n')
        data_set_1 = DataSet.objects.get(id=1)
        data_set_samples_from_1 = DataSetSample.objects.filter(data_submission_from=data_set_1)
        dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_1][:-4])
        custom_args_list = ['--between_type_distances_sample_set', dss_uids_of_ds_str, '--data_analysis_id',
                            '1', '--num_proc', str(self.num_proc), '--distance_method', 'unifrac']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_braycurtis_type_distances_data_set_uid_input(self):
        print('\n\nTesting: stand_alone_unifrac_type_distances_data_set_uid_input\n\n')
        custom_args_list = ['--between_type_distances', '1', '--data_analysis_id', '1', '--num_proc',
                            str(self.num_proc), '--distance_method', 'braycurtis', '--local']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_braycurtis_type_distances_data_set_sample_uid_input(self):
        # create a string for dss input that is all but the last 4 dss objects of the DataSet 1
        print('\n\nTesting: stand_alone_braycurtis_type_distances_data_set_sample_uid_input\n\n')
        data_set_1 = DataSet.objects.get(id=1)
        data_set_samples_from_1 = DataSetSample.objects.filter(data_submission_from=data_set_1)
        dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_1][:-4])
        custom_args_list = ['--between_type_distances_sample_set', dss_uids_of_ds_str, '--data_analysis_id',
                            '1', '--num_proc', str(self.num_proc), '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_braycurtis_type_distances_cct_uid_input(self):
        # create a string for cct input that is the first half of the ccts associated with dataset 1 and
        # the first half from dataset 2
        print('\n\nTesting: stand_alone_braycurtis_type_distances_cct_uid_input\n\n')
        data_analysis = DataAnalysis.objects.get(id=1)
        cct_one_uid_list = [cct.id for cct in CladeCollectionType.objects.filter(
            clade_collection_found_in__data_set_sample_from__data_submission_from__id=1,
            analysis_type_of__data_analysis_from=data_analysis)]
        cct_two_uid_list = [cct.id for cct in CladeCollectionType.objects.filter(
            clade_collection_found_in__data_set_sample_from__data_submission_from__id=2,
            analysis_type_of__data_analysis_from=data_analysis)]
        one_len = len(cct_one_uid_list)
        two_len = len(cct_two_uid_list)
        cct_uid_list = cct_one_uid_list[:int(-1 * (one_len/2))] + cct_two_uid_list[:int(-1*(two_len/2))]
        print('\n\nTesting: stand_alone_braycurtis_type_distances_cct_uid_input\n\n')
        cct_list_str = ','.join([str(_) for _ in cct_uid_list])
        custom_args_list = ['--between_type_distances_cct_set', cct_list_str, '--data_analysis_id', '1',
                            '--num_proc', str(self.num_proc), '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_unifrac_type_distances_cct_uid_input(self):
        # create a string for cct input that is the first half of the ccts associated with dataset 1 and
        # the first half from dataset 2
        print('\n\nTesting: stand_alone_unifrac_type_distances_cct_uid_input\n\n')
        data_analysis = DataAnalysis.objects.get(id=1)
        cct_one_uid_list = [cct.id for cct in CladeCollectionType.objects.filter(
            clade_collection_found_in__data_set_sample_from__data_submission_from__id=1,
            analysis_type_of__data_analysis_from=data_analysis)]
        cct_two_uid_list = [cct.id for cct in CladeCollectionType.objects.filter(
            clade_collection_found_in__data_set_sample_from__data_submission_from__id=2,
            analysis_type_of__data_analysis_from=data_analysis)]
        one_len = len(cct_one_uid_list)
        two_len = len(cct_two_uid_list)
        cct_uid_list = cct_one_uid_list[:int(-1 * (one_len/2))] + cct_two_uid_list[:int(-1*(two_len/2))]
        print('\n\nTesting: stand_alone_unifrac_type_distances_cct_uid_input\n\n')
        cct_list_str = ','.join([str(_) for _ in cct_uid_list])
        custom_args_list = ['--between_type_distances_cct_set', cct_list_str, '--data_analysis_id', '1',
                            '--num_proc', str(self.num_proc), '--distance_method', 'unifrac']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_unifrac_sample_distances_data_set_uid_input(self):
        print('\n\nTesting: stand_alone_unifrac_sample_distances_data_set_uid_input\n\n')
        custom_args_list = ['--between_sample_distances', '1', '--num_proc', str(self.num_proc),
                            '--distance_method', 'unifrac']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_unifrac_sample_distances_data_set_sample_uid_input(self):
        # create a string for dss input that is all but the last 4 dss objects of the DataSet 1
        print('\n\nTesting: stand_alone_unifrac_sample_distances_data_set_sample_uid_input\n\n')
        data_set_1 = DataSet.objects.get(id=1)
        data_set_samples_from_1 = DataSetSample.objects.filter(data_submission_from=data_set_1)
        dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_1][:-4])
        custom_args_list = ['--between_sample_distances_sample_set', dss_uids_of_ds_str, '--num_proc',
                            str(self.num_proc), '--distance_method', 'unifrac']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_braycurtis_sample_distances_data_set_uid_input(self):
        print('\n\nTesting: stand_alone_braycurtis_sample_distances_data_set_uid_input\n\n')
        custom_args_list = ['--between_sample_distances', '1', '--num_proc', str(self.num_proc),
                            '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()

    def test_stand_alone_braycurtis_sample_distances_data_set_sample_uid_input(self):
        # create a string for dss input that is all but the last 4 dss objects of the DataSet 1
        print('\n\nTesting: stand_alone_braycurtis_sample_distances_data_set_sample_uid_input\n\n')
        data_set_1 = DataSet.objects.get(id=1)
        data_set_samples_from_1 = DataSetSample.objects.filter(data_submission_from=data_set_1)
        dss_uids_of_ds_str = ','.join([str(dss.id) for dss in data_set_samples_from_1][:-4])
        custom_args_list = ['--between_sample_distances_sample_set', dss_uids_of_ds_str, '--num_proc',
                            str(self.num_proc), '--distance_method', 'braycurtis']
        test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        test_spwfm.start_work_flow()