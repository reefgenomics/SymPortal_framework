#!/usr/bin/env python3.6
import os
import subprocess
import main
import json
def test_submission(submission_dir='/Users/humebc/Documents/SymPortal_testing_repo/example_data_location',
                    data_sheet_path = '/Users/humebc/Documents/SymPortal_testing_repo/example_data_location/example_meta_data_input.xlsx',
                    has_data_sheet=True, num_proc=6, submission_name='ed_testing'):
    # writeDir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB')) + '/unnamedRefSeqs.fasta'


    path_to_sp_main = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'main.py'))


    if has_data_sheet:
        subprocess.run([path_to_sp_main, '--submit', submission_dir, '--data_sheet', data_sheet_path, '--num_proc', str(num_proc), '--name', submission_name])
    else:
        subprocess.run([path_to_sp_main, '--submit', submission_dir, '--num_proc', str(num_proc), '--name', submission_name])

def test_submission_inline(input_dir, data_sheet_arg, num_proc, name_for_data_set, debug_bool, distance_method_arg,
                           no_ord_arg=False, no_fig_arg=False):


    # get the system type, user_name and user_email
    with open('{}/sp_config'.format(os.path.abspath(os.path.dirname(__file__)))) as f:
        config_dict = json.load(f)
    local_or_remote = config_dict['system_type']
    new_data_set_submitting_user = config_dict['user_name']
    new_data_set_user_email = config_dict['user_email']

    # set whether to run the sub_evaule determination
    # only perform sub_evalue_screening when working on the remote system
    screen_sub_evalue_bool = False

    # generate a data set object that can be used in the submission
    new_data_set = main.create_new_data_set_object_from_params(name_for_data_set, new_data_set_submitting_user,
                                                          new_data_set_user_email)

    # start the main data submission using the paramaters
    data_sub_id = main.start_data_submission(data_sheet_arg, debug_bool, distance_method_arg, input_dir, new_data_set,
                          no_fig_arg, no_ord_arg, num_proc, screen_sub_evalue_bool)

    return data_sub_id

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

    # do a submission test
    data_sub_id = test_submission_inline(input_dir='/Users/humebc/Documents/SymPortal_testing_repo/example_data_location',
                        data_sheet_arg = '/Users/humebc/Documents/SymPortal_testing_repo/example_data_location/example_meta_data_input.xlsx',
                        num_proc=6, name_for_data_set='ed_testing', debug_bool=False, distance_method_arg='braycurtis',
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
