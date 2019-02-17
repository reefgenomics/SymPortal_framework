#!/usr/bin/env python3.6
import os
import main
import json
import general
import sys
from pathlib import Path
import shutil
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from dbApp.models import DataAnalysis, DataSet
import sp_config
import argparse
from output import SequenceCountTableCreator

class SymPortalTester:
    def __init__(self):
        self.symportal_testing_root_dir = os.path.abspath(os.path.dirname(__file__))
        self.system_type = sp_config.system_type
        self.user_name = sp_config.user_name
        self.user_email = sp_config.user_email
        self.data_set_object = None
        self.data_analysis_object = None
        self.argument_parser = self._setup_argument_parser()
        self.list_of_args_for_arg_parser = None
        self.args = None
        self.test_data_dir_path = os.path.join(self.symportal_testing_root_dir, 'data')
        self.data_sheet_file_path = os.path.join(self.test_data_dir_path, 'test_data_submission_input.csv')
        self.completed_data_loading_objects = None
        self.num_proc = 6

    def _setup_argument_parser(self):
        self.argument_parser = argparse.ArgumentParser(
            description='Intragenomic analysis of the ITS2 region of the nrDNA',
            epilog='For support email: symportal@gmail.com')

        group = self.argument_parser.add_mutually_exclusive_group(required=True)

        main.define_mutually_exclusive_args(group)
        main.define_additional_args(group, self.argument_parser)
        return self.argument_parser

    def test_data_loading(self):
        self.args = self.argument_parser.parse_args(
            ['--load', self.test_data_dir_path, '--name', 'testing', '--num_proc', '6', '--data_sheet',
             self.data_sheet_file_path])

        self.completed_data_loading_objects = main.perform_data_loading(self.args)



def test_data_loading_without_command_line(
        test_data_directory, data_sheet_path, num_proc, debug_bool, distance_method_arg,
        no_ord_arg=False, no_fig_arg=False):
    """
    This is the function used to test the loading of data into the database. It will work directly with the
    create_data_submission.py functions rather than by running command line arguments. For the command line version
    please see test_data_loading_to_database_from_command_line above.
    :param test_data_directory: The directory in which the test data is held
    :param data_sheet_path: the data sheet containing the sample names and meta data for the test sample data
    :param num_proc: the number of processors that should be used in the loading
    :param debug_bool: whether the load should be run in debug mode, default=False
    :param distance_method_arg: which distance method should be used (default=braycurtis, alternative=unifrac)
    :param no_ord_arg: bool, whether ordination should not be performed. default = False (ordinations are produced)
    :param no_fig_arg: bool, whether figures should not be produced. default = False (figures are produced)
    :return:
    """

    submitting_user, user_email, system_type = get_system_type_user_name_and_user_email()

    data_set_uid, path_of_output_items_list = perform_load_testing_based_on_system_type(
        data_sheet_path, debug_bool, distance_method_arg, no_fig_arg, no_ord_arg, num_proc,
        system_type, test_data_directory, submitting_user, user_email)

    return data_set_uid, path_of_output_items_list


def perform_load_testing_based_on_system_type(
        data_sheet_path, debug_bool, distance_method_arg, no_fig_arg, no_ord_arg,
        num_proc, system_type, test_data_directory, submitting_user, user_email):
    if system_type == 'remote':
        data_set_uid, path_of_output_items_list = perform_loading_tests_for_remote_system(
            data_sheet_path, debug_bool, distance_method_arg,
            no_fig_arg, no_ord_arg, num_proc, test_data_directory, submitting_user, user_email)

    elif system_type == 'local':
        data_set_uid, path_of_output_items_list = perform_loading_tests_for_local_system(
            data_sheet_path, debug_bool, distance_method_arg,
            no_fig_arg, no_ord_arg, num_proc, test_data_directory, submitting_user, user_email)
    else:
        sys.exit('Unknown system type in settings.py. System type should be set to \'local\' in most instances.')
    return data_set_uid, path_of_output_items_list


def perform_loading_tests_for_remote_system(
        data_sheet_path, debug_bool, distance_method_arg, no_fig_arg,
        no_ord_arg, num_proc, test_data_directory, submitting_user, user_email):

    # perform a load test: 1 - with no sub_e_screening and no data sheet
    perform_load_test_one_remote_system(
        debug_bool, distance_method_arg, no_fig_arg, no_ord_arg, num_proc,
        submitting_user, test_data_directory, user_email)

    # perform a load test: 2 - with sub_e_screening and data sheet
    data_set_uid_two, path_of_output_items_list_two = perform_load_test_two_remote_system(
        data_sheet_path, debug_bool, distance_method_arg, no_fig_arg, no_ord_arg,
        num_proc, submitting_user, test_data_directory, user_email)

    return data_set_uid_two, path_of_output_items_list_two


def perform_load_test_two_remote_system(data_sheet_path, debug_bool, distance_method_arg, no_fig_arg, no_ord_arg,
                                        num_proc, submitting_user, test_data_directory, user_email):
    print('Loading test 2/2:')
    print('Loading test with data_sheet and with screen_sub_evalue_bool == True')
    new_data_set = main.create_new_data_set_object_from_params('testing', submitting_user, user_email)
    data_set_uid_two, path_of_output_items_list_two = perform_loading_test(
        data_sheet_arg=data_sheet_path, debug_bool=debug_bool, distance_method_arg=distance_method_arg,
        input_dir=test_data_directory, new_data_set=new_data_set, no_fig_arg=no_fig_arg, no_ord_arg=no_ord_arg,
        num_proc=num_proc, screen_sub_evalue_bool=True)
    print('Loading test 2/2 SUCCESSFUL')
    return data_set_uid_two, path_of_output_items_list_two


def perform_load_test_one_remote_system(debug_bool, distance_method_arg, no_fig_arg, no_ord_arg, num_proc,
                                        submitting_user, test_data_directory, user_email):
    print('Loading test 1/2:')
    print('Loading test with no data_sheet and with screen_sub_evalue_bool == False')
    new_data_set = main.create_new_data_set_object_from_params('testing', submitting_user, user_email)
    data_set_uid_one, path_of_output_items_list_one = perform_loading_test(
        data_sheet_arg=False, debug_bool=debug_bool, distance_method_arg=distance_method_arg,
        input_dir=test_data_directory, new_data_set=new_data_set, no_fig_arg=no_fig_arg, no_ord_arg=no_ord_arg,
        num_proc=num_proc, screen_sub_evalue_bool=False)
    print('Loading test 1/2 SUCCESSFUL')
    clean_up_after_test_loading(data_set_uid_one)


def perform_loading_tests_for_local_system(data_sheet_path, debug_bool, distance_method_arg, no_fig_arg,
                                           no_ord_arg, num_proc, test_data_directory, submitting_user, user_email):

    # perform a load test: 1 - with no sub_e_screening and no data sheet

    perform_load_test_one_local_system(
        debug_bool, distance_method_arg, no_fig_arg, no_ord_arg,
        num_proc, test_data_directory, submitting_user, user_email)

    # perform a load test: 2 - no sub_e_screening, with data sheet
    data_set_uid_two, path_of_output_items_list_two = perform_load_test_two_local_system(
        data_sheet_path, debug_bool, distance_method_arg, no_fig_arg, no_ord_arg, num_proc,
        test_data_directory, submitting_user, user_email)

    return data_set_uid_two, path_of_output_items_list_two


def perform_load_test_two_local_system(data_sheet_path, debug_bool, distance_method_arg, no_fig_arg,
                                       no_ord_arg, num_proc, test_data_directory, submitting_user,
                                       user_email):
    print('Loading test 2/2:')
    print('Loading test with data_sheet and with screen_sub_evalue_bool == True')
    new_data_set = main.create_new_data_set_object_from_params('testing', submitting_user, user_email)
    data_set_uid_two, path_of_output_items_list_two = perform_loading_test(
        data_sheet_arg=data_sheet_path, debug_bool=debug_bool, distance_method_arg=distance_method_arg,
        input_dir=test_data_directory, new_data_set=new_data_set, no_fig_arg=no_fig_arg, no_ord_arg=no_ord_arg,
        num_proc=num_proc, screen_sub_evalue_bool=False)
    print('Loading test 2/2 SUCCESSFUL')
    return data_set_uid_two, path_of_output_items_list_two


def perform_load_test_one_local_system(debug_bool, distance_method_arg, no_fig_arg, no_ord_arg, num_proc,
                                       test_data_directory, submitting_user, user_email):
    print('Loading test 1/2:')
    print('Loading test with no data_sheet and with screen_sub_evalue_bool == False')
    new_data_set = main.create_new_data_set_object_from_params('testing', submitting_user, user_email)
    data_set_uid_one, path_of_output_items_list_one = perform_loading_test(
        data_sheet_arg=False, debug_bool=debug_bool, distance_method_arg=distance_method_arg,
        input_dir=test_data_directory, new_data_set=new_data_set, no_fig_arg=no_fig_arg, no_ord_arg=no_ord_arg,
        num_proc=num_proc, screen_sub_evalue_bool=False)
    print('Loading test 1/2 SUCCESSFUL')
    clean_up_after_test_loading(data_set_uid_one)


def clean_up_after_test_loading(data_set_uid_one):
    delete_data_load_output_object_directories(data_set_uid_one)
    print('Output directory deleted')
    # delete this data_set
    print('test data_set {} deleted'.format(data_set_uid_one))
    general.delete_data_set(data_set_uid_one)


def clean_up_after_test_analysis(data_analysis_uid):
    delete_data_analysis_output_object_directories(data_analysis_uid)
    print('Output directory deleted')
    # delete this data_set
    print('test data_analysis {} deleted'.format(data_analysis_uid))
    general.delete_data_analysis(data_analysis_uid)


def perform_loading_test(
        data_sheet_arg, debug_bool, distance_method_arg, input_dir, new_data_set, no_fig_arg, no_ord_arg, num_proc,
        screen_sub_evalue_bool):

    data_set_uid, output_path_list = main.load_data(
        data_sheet_arg, debug_bool, distance_method_arg, input_dir,
        new_data_set, no_fig_arg, no_ord_arg, num_proc, screen_sub_evalue_bool
    )
    return data_set_uid, output_path_list


def get_system_type_user_name_and_user_email():
    sp_config_path = os.path.abspath(os.path.join(Path(__file__).parents[1], 'sp_config'))
    with open(sp_config_path, 'r') as f:
        config_dict = json.load(f)
    local_or_remote = config_dict['system_type']
    new_data_set_submitting_user = config_dict['user_name']
    new_data_set_user_email = config_dict['user_email']
    return new_data_set_submitting_user, new_data_set_user_email, local_or_remote


def test_analysis_without_command_line(
        custom_data_set_ids, analysis_name='test_analysis', description_arg='test_analysis',
        distance_method_arg='braycurtis', no_fig_arg=False, no_ord_arg=False, no_output_arg=False, num_proc=6,
        within_clade_cutoff=0.03, debug_bool=False):

    """
    This will run an analysis using the data_set that was generated in the loading tests
    :param custom_data_set_ids: a comma separated string of the data_set uids that will be included in the analysis
    :param analysis_name:
    :param description_arg:
    :param distance_method_arg: the distance matrix to use for the ordination.
    default='braycurtis', alternative='unifrac'
    :param no_ord_arg: bool, whether ordination should not be performed. default = False (ordinations are produced)
    :param no_fig_arg: bool, whether figures should not be produced. default = False (figures are produced)
    :param no_output_arg: bool. if true no figures or ordinations will be produced
    :param num_proc:
    :param within_clade_cutoff: The relative abundance a sequence must be at to be considered in the initial footprints
    :param debug_bool: Whether verbose output should be switched on to aid in debugging.
    :return:
    """

    submitting_user, user_email, system_type = get_system_type_user_name_and_user_email()

    print('Testing analysis 1/1')
    analysis_uid, output_path_list = main.create_analysis_obj_and_run_analysis(
        analysis_name, description_arg, custom_data_set_ids, debug_bool, distance_method_arg, submitting_user,
        user_email, no_fig_arg, no_ord_arg, no_output_arg, num_proc, within_clade_cutoff)
    print('Analysis testing SUCCESSFUL')
    clean_up_after_test_analysis(analysis_uid)
    return analysis_uid, output_path_list


def initiate_test():
    """
    This is the main entry to the testing. It is currently composed of loading a set of data into the database and then
    performing a standalone analysis (i.e. an analysis that only contains the loaded data_set). I have created a set
    of subsampled fastq files from the Smith et al dataset used in the manuscript. Each sample has been subsampled to
    500 sequences. The data is housed in ./tests/data. I have also created a .csv data submission sheet that is housed
    in the same directory. This can be used in the loading of the data.
    :return:
    """

    cleanup_after_previous_tests()

    # run_loading_tests
    test_data_directory = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
    data_sheet_path = os.path.join(test_data_directory, 'test_data_submission_input.csv')

    data_set_uid, path_of_output_items_list = test_data_loading_without_command_line(
        test_data_directory=test_data_directory, data_sheet_path=data_sheet_path,
        num_proc=6, debug_bool=False, distance_method_arg='braycurtis',
        no_ord_arg=False, no_fig_arg=False)

    # run_analysis_tests
    test_analysis_without_command_line(
        custom_data_set_ids=str(data_set_uid), analysis_name='testing',
        description_arg='testing',  distance_method_arg='braycurtis', no_fig_arg=False,
        no_ord_arg=False, no_output_arg=False, num_proc=6, within_clade_cutoff=0.03, debug_bool=False)

    # todo run some checks on the output files to verify there are no subtle bugs

    # todo run the sample by sample ordination

    # todo run all ordinations with unifrac method

    # todo get a list of samples from the submission and run a sample-wise ordination

    # todo run outputs individually for both seqs and for type output

    print('testing complete')


def cleanup_after_previous_tests():
    cleanup_after_previously_run_data_analysis_tests()
    cleanup_after_previously_run_data_loading_tests()


def cleanup_after_previously_run_data_loading_tests():
    uids_of_testing_data_set_objects = (ds.id for ds in DataSet.objects.filter(name='testing'))
    for ds_uid in uids_of_testing_data_set_objects:
        print(f'Cleaning up after previous data loading test: {ds_uid}')
        delete_data_load_output_object_directories(ds_uid)
        general.delete_data_set(ds_uid)


def cleanup_after_previously_run_data_analysis_tests():
    uids_of_testing_data_analysis_objects = (da.id for da in DataAnalysis.objects.filter(name='testing'))
    for da_uid in uids_of_testing_data_analysis_objects:
        print(f'Cleaning up after previous data analysis test: {da_uid}')
        delete_data_analysis_output_object_directories(da_uid)
        general.delete_data_analysis(da_uid)


def delete_data_load_output_object_directories(data_set_uid):
    directory_to_delete = os.path.abspath(os.path.join(
        Path(__file__).parents[1], 'outputs', 'data_set_submissions', str(data_set_uid)))
    if os.path.exists(directory_to_delete):
        shutil.rmtree(directory_to_delete)


def delete_data_analysis_output_object_directories(data_analysis_uid):
    directory_to_delete = os.path.abspath(os.path.join(
        Path(__file__).parents[1], 'outputs', 'analyses', str(data_analysis_uid)))
    if os.path.exists(directory_to_delete):
        shutil.rmtree(directory_to_delete)


if __name__ == "__main__":
    symportal_tester = SymPortalTester()
    symportal_tester.test_data_loading()