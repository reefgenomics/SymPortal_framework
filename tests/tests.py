#!/usr/bin/env python3.6
import os
import subprocess
import main
import json

def test_data_loading_to_database_from_command_line(submission_dir='/Users/humebc/Documents/SymPortal_testing_repo/example_data_location',
                    data_sheet_path = '/Users/humebc/Documents/SymPortal_testing_repo/example_data_location/example_meta_data_input.xlsx',
                    has_data_sheet=True, num_proc=6, submission_name='ed_testing'):
    # writeDir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB')) + '/unnamedRefSeqs.fasta'


    path_to_sp_main = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'main.py'))


    if has_data_sheet:
        subprocess.run([path_to_sp_main, '--submit', submission_dir, '--data_sheet', data_sheet_path, '--num_proc', str(num_proc), '--name', submission_name])
    else:
        subprocess.run([path_to_sp_main, '--submit', submission_dir, '--num_proc', str(num_proc), '--name', submission_name])

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

    # start the main data submission using the paramaters
    # todo if system_type == 'remote' then we should test with screen_sub_evalue_bool set to true and False
    # elif local then we should only test screen_sub_evalue_bool = False

    if system_type == 'remote':

        # perform a load test 1 - with sub_e_screening and data sheet
        print('Loading test 1/4:')
        print('Loading test with data_sheet and with screen_sub_evalue_bool == False')
        data_sub_uid_one = perform_loading_test(data_sheet_arg=data_sheet_path,
                                                debug_bool=debug_bool,
                                                distance_method_arg=distance_method_arg,
                                                input_dir=test_data_directory,
                                                new_data_set=new_data_set,
                                                no_fig_arg=no_fig_arg, no_ord_arg=no_ord_arg, num_proc=num_proc,
                                                screen_sub_evalue_bool=True)
        print('Loading test 1/4 SUCCESSFUL')

        print('Loading test 2/4:')
        print('Loading test with data_sheet and with screen_sub_evalue_bool == False')
        data_sub_id = main.start_data_submission(data_sheet_arg, debug_bool, distance_method_arg, input_dir,
                                                 new_data_set,
                                                 no_fig_arg, no_ord_arg, num_proc, screen_sub_evalue_bool)
        print('Loading test 1/4 SUCCESSFUL')

    return data_sub_id


def perform_loading_test(data_sheet_arg, debug_bool, distance_method_arg, input_dir,
                                                new_data_set, no_fig_arg, no_ord_arg, num_proc, screen_sub_evalue_bool):

    data_sub_id = main.start_data_submission(data_sheet_arg, debug_bool, distance_method_arg, input_dir, new_data_set,
                                             no_fig_arg, no_ord_arg, num_proc, screen_sub_evalue_bool)



def get_system_type_user_name_and_user_email():
    with open('{}/sp_config'.format(os.path.abspath(os.path.dirname(__file__)))) as f:
        config_dict = json.load(f)
    local_or_remote = config_dict['system_type']
    new_data_set_submitting_user = config_dict['user_name']
    new_data_set_user_email = config_dict['user_email']
    return new_data_set_submitting_user, new_data_set_user_email, local_or_remote


def test_analysis_inline(custom_data_set_ids, analysis_name='ed_testing_analysis',
                         description_arg='ed_testing_analysis', distance_method_arg='braycurtis',
                         no_fig_arg=False, no_ord_arg=False, no_output_arg=False, num_proc=6,
                         within_clade_cutoff = 0.03, debug_bool=False):

    this = os.path.abspath(os.path.dirname(__file__))
    print('this = {}'.format(this))

    with open('{}/sp_config'.format(os.path.abspath(os.path.dirname(__file__)))) as f:
        config_dict = json.load(f)
    local_or_remote = config_dict['system_type']

    new_data_set_submitting_user = config_dict['user_name']
    new_data_set_user_email = config_dict['user_email']
    screen_sub_evalue_bool = False

    analysis_id = main.create_analysis_obj_and_run_analysis(analysis_name, description_arg, custom_data_set_ids, debug_bool,
                                         distance_method_arg,
                                         new_data_set_submitting_user, new_data_set_user_email, no_fig_arg,
                                         no_ord_arg, no_output_arg, num_proc, within_clade_cutoff)


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
    data_sub_id = test_data_loading_without_command_line(test_data_directory=test_data_directory,
                                                         data_sheet_path = data_sheet_path,
                                                         num_proc=6, name_for_data_set='testing',
                                                         debug_bool=False, distance_method_arg='braycurtis',
                                                         no_ord_arg=False, no_fig_arg=False)

    # do an analysis test
    analysis_id = test_analysis_inline(custom_data_set_ids=str(data_sub_id), analysis_name='ed_testing_analysis',
                                       description_arg='ed_testing_analysis',  distance_method_arg='braycurtis',
                                       no_fig_arg=False, no_ord_arg=False, no_output_arg=False, num_proc=6,
                                       within_clade_cutoff = 0.03, debug_bool=False)

    #todo run the sample by sample ordination

    #todo run all ordinations with unifrac method

    #todo get a list of samples from the submission and run a sample-wise ordination

    #todo run outputs individually for both seqs and for type output

    #todo generate lightweight fastq files for the testing.


    print('testing complete')
if __name__ == "__main__":
    initiate_test()
