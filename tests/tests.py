#!/usr/bin/env python3.6
import os
import subprocess
import main
import json
import general
import sys

def test_data_loading_without_command_line(test_data_directory, data_sheet_path, num_proc, name_for_data_set, debug_bool, distance_method_arg,
                           no_ord_arg=False, no_fig_arg=False):
    """
    This is the function used to test the loading of data into the database. It will work directly with the
    create_data_submission.py functions rather than by running command line arguments. For the command line version
    please see test_data_loading_to_database_from_command_line above.
    :param input_dir: The directory in which the test data is held
    :param data_sheet_arg: the data sheet containing the sample names and meta data for the test sample data
    :param num_proc: the number of processors that should be used in the loading
    :param name_for_data_set: name of the data_set
    :param debug_bool: whether the load should be run in debug mode, default=False
    :param distance_method_arg: which distance method should be used (default=braycurtis, alternative=unifrac)
    :param no_ord_arg: bool, whether ordination should not be performed. default = False (ordinations are produced)
    :param no_fig_arg: bool, whether figures should not be produced. default = False (figures are produced)
    :return:
    """

    submitting_user, user_email, system_type = get_system_type_user_name_and_user_email()

    # generate a data set object that can be used in the submission
    new_data_set = main.create_new_data_set_object_from_params(name_for_data_set, submitting_user, user_email)

    if system_type == 'remote':
        data_set_uid, path_of_output_items_list = perform_loading_tests_for_remote_system(
            data_sheet_path, debug_bool, distance_method_arg,
            new_data_set, no_fig_arg, no_ord_arg, num_proc, test_data_directory)

    elif system_type == 'local':
        data_set_uid, path_of_output_items_list = perform_loading_tests_for_local_system(
            data_sheet_path, debug_bool, distance_method_arg,
            new_data_set, no_fig_arg, no_ord_arg, num_proc, test_data_directory)
    else:
        sys.exit('Unknown system type in settings.py. System type should be set to \'local\' in most instances.')

    return data_set_uid, path_of_output_items_list


def perform_loading_tests_for_remote_system(data_sheet_path, debug_bool, distance_method_arg, new_data_set, no_fig_arg,
                                            no_ord_arg, num_proc, test_data_directory):
    # perform a load test: 1 - with no sub_e_screening and no data sheet
    print('Loading test 1/2:')
    print('Loading test with no data_sheet and with screen_sub_evalue_bool == False')
    data_set_uid_one, path_of_output_items_list_one = perform_loading_test(
        data_sheet_arg=False, debug_bool=debug_bool, distance_method_arg=distance_method_arg,
        input_dir=test_data_directory, new_data_set=new_data_set, no_fig_arg=no_fig_arg, no_ord_arg=no_ord_arg,
        num_proc=num_proc, screen_sub_evalue_bool=False)
    print('Loading test 1/2 SUCCESSFUL')
    # delete this data_set

    # todo see which directories need deleting and write it into a function specific for deleteing files
    # from a data_load

    general.delete_data_set(data_set_uid_one)

    # perform a load test: 2 - with sub_e_screening and data sheet
    print('Loading test 2/2:')
    print('Loading test with data_sheet and with screen_sub_evalue_bool == True')
    data_set_uid_two, path_of_output_items_list_two = perform_loading_test(
        data_sheet_arg=data_sheet_path, debug_bool=debug_bool, distance_method_arg=distance_method_arg,
        input_dir=test_data_directory, new_data_set=new_data_set, no_fig_arg=no_fig_arg, no_ord_arg=no_ord_arg,
        num_proc=num_proc, screen_sub_evalue_bool=True)
    print('Loading test 2/2 SUCCESSFUL')
    return data_set_uid_two, path_of_output_items_list_two


def perform_loading_tests_for_local_system(data_sheet_path, debug_bool, distance_method_arg, new_data_set, no_fig_arg,
                                           no_ord_arg, num_proc, test_data_directory):
    # perform a load test: 1 - with no sub_e_screening and no data sheet
    print('Loading test 1/2:')
    print('Loading test with no data_sheet and with screen_sub_evalue_bool == False')
    data_set_uid_one, path_of_output_items_list_one = perform_loading_test(
        data_sheet_arg=False, debug_bool=debug_bool, distance_method_arg=distance_method_arg,
        input_dir=test_data_directory, new_data_set=new_data_set, no_fig_arg=no_fig_arg, no_ord_arg=no_ord_arg,
        num_proc=num_proc, screen_sub_evalue_bool=False)

    print('Loading test 1/2 SUCCESSFUL')
    # delete this data_set
    # todo see which directories need deleting and write it into a function specific for deleteing files
    # from a data_load
    general.delete_data_set(data_set_uid_one)
    # perform a load test: 2 - no sub_e_screening, with data sheet
    print('Loading test 2/2:')
    print('Loading test with data_sheet and with screen_sub_evalue_bool == True')
    data_set_uid_two, path_of_output_items_list_two = perform_loading_test(
        data_sheet_arg=data_sheet_path, debug_bool=debug_bool, distance_method_arg=distance_method_arg,
        input_dir=test_data_directory, new_data_set=new_data_set, no_fig_arg=no_fig_arg, no_ord_arg=no_ord_arg,
        num_proc=num_proc, screen_sub_evalue_bool=False)

    return data_set_uid_two, path_of_output_items_list_two


def perform_loading_test(
        data_sheet_arg, debug_bool, distance_method_arg, input_dir, new_data_set, no_fig_arg, no_ord_arg, num_proc,
        screen_sub_evalue_bool):

    data_set_uid, output_path_list = main.start_data_submission(
        data_sheet_arg, debug_bool, distance_method_arg, input_dir,
        new_data_set, no_fig_arg, no_ord_arg, num_proc, screen_sub_evalue_bool
    )
    return data_set_uid, output_path_list


def get_system_type_user_name_and_user_email():
    with open('{}/sp_config'.format(os.path.abspath(os.path.dirname(__file__)))) as f:
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

    # run_loading_tests
    test_data_directory = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
    data_sheet_path = test_data_directory + 'test_data_submission_input.csv'

    data_set_uid, path_of_output_items_list = test_data_loading_without_command_line(
        test_data_directory=test_data_directory, data_sheet_path=data_sheet_path,
        num_proc=6, name_for_data_set='testing', debug_bool=False, distance_method_arg='braycurtis',
        no_ord_arg=False, no_fig_arg=False)

    # run_analysis_tests
    analysis_uid, output_path_list = test_analysis_without_command_line(
        custom_data_set_ids=str(data_set_uid), analysis_name='ed_testing_analysis',
        description_arg='ed_testing_analysis',  distance_method_arg='braycurtis', no_fig_arg=False,
        no_ord_arg=False, no_output_arg=False, num_proc=6, within_clade_cutoff=0.03, debug_bool=False)

    # todo create a function that will delete the files from the data_load and the analysis and run here

    # todo run the sample by sample ordination

    # todo run all ordinations with unifrac method

    # todo get a list of samples from the submission and run a sample-wise ordination

    # todo run outputs individually for both seqs and for type output

    print('testing complete')


if __name__ == "__main__":
    initiate_test()

# old code
def test_data_loading_to_database_from_command_line(
        submission_dir='/Users/humebc/Documents/SymPortal_testing_repo/example_data_location',
        data_sheet_path = '/Users/humebc/Documents/SymPortal_testing_repo/'
                          'example_data_location/example_meta_data_input.xlsx',
        has_data_sheet=True, num_proc=6, submission_name='ed_testing'):

    path_to_sp_main = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'main.py'))

    if has_data_sheet:
        subprocess.run([path_to_sp_main, '--submit', submission_dir, '--data_sheet', data_sheet_path,
                        '--num_proc', str(num_proc), '--name', submission_name])
    else:
        subprocess.run([path_to_sp_main, '--submit', submission_dir, '--num_proc', str(num_proc),
                        '--name', submission_name])