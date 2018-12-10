from dbApp.models import data_set, reference_sequence, data_set_sample_sequence, analysis_type, analysis_group, \
    data_set_sample, data_analysis, clade_collection, clade_collection_type
import os
import shutil
from multiprocessing import Queue, Process, Manager, current_process
import math
import sys
from plumbum import local
import pandas as pd
from django import db
import subprocess
import re
import numpy as np
from skbio.stats.ordination import pcoa
from general import writeListToDestination, readDefinedFileToList, convert_interleaved_to_sequencial_fasta
import itertools
from scipy.spatial.distance import braycurtis
from datetime import datetime


##### DISTANCE MATRICES #####
def generate_within_clade_UniFrac_distances_ITS2_type_profiles(data_submission_id_str, num_processors, data_analysis_id,
                                                               method, call_type, date_time_string, bootstrap_value=100, noFig=False, output_dir=None):
    ''' This will produce distance matrices between ITS2 type profiles of the same clade.
    It will use exactly the same suite of programs as the generate_within_clade_UniFrac_distances_samples method
    The only difference will be that this will be calculated across ITS2 Type Profiles rather than samples.
    As such, this will require a data_set set and dataAnalysis.
    The average DIV abundances for each of the ITS2 type profiles should first be calculated from which the distances
    can then be calcualated'''

    output_file_paths = []

    if call_type == 'stand_alone':
        wkd = os.path.abspath(os.path.join(os.path.dirname(__file__), 'outputs',
                                           'ordination', '_'.join(str(data_submission_id_str).split(',')),
                                           'between_profiles'))
    else:
        # call_type = 'analysis'
        wkd = '{}/{}'.format(output_dir, 'between_its2_type_profile_distances')

    # get the dataSubmissions
    data_submissions = data_set.objects.filter(id__in=[int(a) for a in str(data_submission_id_str).split(',')])

    # get the dataAnalysis
    data_analysis_obj = data_analysis.objects.get(id=data_analysis_id)

    # go clade by clade
    ITS2_type_profiles_of_data_subs_and_analysis = analysis_type.objects.filter(
        clade_collection_type__cladeCollectionFoundIn__dataSetSampleFrom__dataSubmissionFrom__in=
        data_submissions, dataAnalysisFrom=data_analysis_obj).distinct()

    clade_list = list(set([type_profile.clade for type_profile in ITS2_type_profiles_of_data_subs_and_analysis]))

    # list for storing the paths to the PCoA .csv that will be produced
    PCoA_path_lists = []
    for clade in clade_list:
        # convert query to a list so that we can iterate over twice in exactly the same order
        ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade = list(
            ITS2_type_profiles_of_data_subs_and_analysis.filter(clade=clade))

        if len(ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade) < 2:
            continue

        group_file, name_file, unique_fasta = generate_fasta_name_group_between_profiles(
            ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade)

        # at this point we have the fasta, name and group file
        # write out
        clade_wkd = wkd + '/{}'.format(clade)
        writeListToDestination('{}/unique.fasta'.format(clade_wkd), unique_fasta)
        writeListToDestination('{}/name_file.names'.format(clade_wkd), name_file)
        writeListToDestination('{}/group_file.groups'.format(clade_wkd), group_file)

        # align fasta
        out_file = mafft_align_fasta(clade_wkd, num_proc=num_processors)

        # Generate random data sets
        fseqboot_base = generate_fseqboot_alignments(clade_wkd, bootstrap_value, out_file)

        if method == 'mothur':
            unifrac_dist, unifrac_path = mothur_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file,
                                                                    bootstrap_value, num_processors)
        elif method == 'phylip':
            unifrac_dist, unifrac_path = phylip_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file,
                                                                    bootstrap_value, num_processors)

        PCoA_path = generate_PCoA_coords(clade_wkd, unifrac_dist, date_time_string)

        output_file_paths.append(PCoA_path)
        output_file_paths.append(unifrac_path)

        PCoA_path_lists.append(PCoA_path)
        # Delete the tempDataFolder and contents
        file_to_del = '{}/out_seq_boot_reps'.format(clade_wkd)
        shutil.rmtree(path=file_to_del)

        # now delte all files except for the .csv that holds the coords and the .dist that holds the dists
        list_of_dir = os.listdir(clade_wkd)
        for item in list_of_dir:
            if '.csv' not in item and '.dist' not in item:
                os.remove(os.path.join(clade_wkd, item))

    # output the paths of the new files created
    print('\n\nOutput files:\n')
    for path_of_output_file in output_file_paths:
        print(path_of_output_file)

    return PCoA_path_lists


def generate_fasta_name_group_between_profiles(ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade):
    unique_fasta = []
    unique_seq_id_list = []
    name_file_dict = {}
    ref_seq_id_to_name_file_representative = {}
    group_file = []
    # we will need to get a list of all of the sequences that are found in the above ITS2 type profiles
    # we will then need to align these.
    # we will also need to produce a group file, a name file and a unique fasta file
    # this will be similar to what we have in the generate_within_clade_UniFrac_distances_samples
    # from this point on it should be easy for us to use the same set up as we have in the other method
    # we already have the ratios of each of the divs
    list_of_type_profile_norm_abundance_dicts = []
    for at in ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade:
        sys.stdout.write('\rProcessing {}'.format(at))
        list_of_div_ids = [int(b) for b in at.orderedFootprintList.split(',')]
        foot_print_ratio_array = pd.DataFrame(at.getRatioList())
        normalised_abundance_of_divs_dict = {list_of_div_ids[i]: math.ceil(foot_print_ratio_array[i].mean() * 1000) for
                                             i in range(len(list_of_div_ids))}
        list_of_type_profile_norm_abundance_dicts.append(normalised_abundance_of_divs_dict)
        # go div by div
        for ref_seq_id_key in normalised_abundance_of_divs_dict.keys():
            if ref_seq_id_key in unique_seq_id_list:
                # then this only needs appending to the group file and to the name file
                # the fasta already contains a representative of this.
                temp_name_list = []
                for i in range(normalised_abundance_of_divs_dict[ref_seq_id_key]):
                    integer_seq_name = '{}_id{}_{}'.format(ref_seq_id_key, at.id, i)
                    temp_name_list.append(integer_seq_name)
                    group_file.append('{}\t{}'.format(integer_seq_name, str(at).replace('-', '_').replace('/', '_')))
                name_file_dict[ref_seq_id_to_name_file_representative[ref_seq_id_key]] += temp_name_list
            else:
                # then this needs adding to the fasta
                zero_seq_name = '{}_id{}_{}'.format(ref_seq_id_key, at.id, 0)
                unique_fasta.append('>{}'.format(zero_seq_name))
                # get sequence of ref_seq
                sequence = reference_sequence.objects.get(id=ref_seq_id_key).sequence
                unique_fasta.append(sequence)
                unique_seq_id_list.append(ref_seq_id_key)

                # the name file needs initiating
                # create the temp_name file and add to the master
                # can do group file at same time
                temp_name_list = []
                for i in range(normalised_abundance_of_divs_dict[ref_seq_id_key]):
                    integer_seq_name = '{}_id{}_{}'.format(ref_seq_id_key, at.id, i)
                    temp_name_list.append(integer_seq_name)
                    group_file.append('{}\t{}'.format(integer_seq_name, str(at).replace('-', '_').replace('/', '_')))
                name_file_dict[zero_seq_name] = temp_name_list
                ref_seq_id_to_name_file_representative[ref_seq_id_key] = zero_seq_name

    # now assemble the name file
    name_file = []
    for k, v in name_file_dict.items():
        name_file.append('{}\t{}'.format(k, ','.join(v)))

    return group_file, name_file, unique_fasta

def generate_within_clade_UniFrac_distances_samples_sample_list_input(smpl_id_list_str, num_processors,
                                                    method, bootstrap_value=100, date_time_string=None):



    '''
    This method will be used to generate another .dist file and PCoA coordinates file for a specific set of
    samples. This is super useful when you want to further investigate resolutions according to the distance data.
    TODO We will first implement this in the UniFrac methodology so that we can use this for the SP MS but then
    we should also implement with the BrayCurtis method.
    '''

    output_file_paths = []
    
    smpl_id_list = [int(str_id) for str_id in smpl_id_list_str.split(',')]
    sample_list = data_set_sample.objects.filter(id__in=smpl_id_list)
    clade_collection_list_of_samples = clade_collection.objects.filter(
        dataSetSampleFrom__in=sample_list)

    clades_of_clade_collections = list(set([a.clade for a in clade_collection_list_of_samples]))


    wkd = os.path.abspath(
        os.path.join(os.path.dirname(__file__), 'outputs', 'ordination', 'custom_sample_list',
                     'between_samples', date_time_string))

    # for each clade found in the dataSubmissions' samples
    PCoA_path_lists = []
    for clade_in_question in clades_of_clade_collections:

        clade_wkd = wkd + '/{}'.format(clade_in_question)
        clade_collections_of_clade = clade_collection_list_of_samples.filter(clade=clade_in_question)
        if len(clade_collections_of_clade) < 2:
            continue
        # directly make a name file and unique fasta file

        # NB the MP that we will be writing below is complex and cannot be completed entirely within the worker
        # we will out put to the output queue live and process this in the main thread to create the master
        # fasta and name file

        sys.stdout.write('Creating .name and .fasta files')

        # setup MP
        # Create the queues that will hold the cc information
        taskQueue = Queue()

        # output queue to put the dictionaries that will make the fasta and name files
        outputQueue = Queue()

        # populate the input queue
        for cc in clade_collections_of_clade:
            taskQueue.put(cc)

        # place the stop cues
        for n in range(num_processors):
            taskQueue.put('STOP')

        allProcesses = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        for n in range(num_processors):
            p = Process(target=uni_frac_worker_two, args=(taskQueue, outputQueue))
            allProcesses.append(p)
            p.start()

        # have an output pipe that gets processed as we go along. I.e. before we call wait.
        master_fasta_dict = {}
        master_unique_fasta_seq_list = []
        master_name_dict = {}
        master_group_list = []
        # this dict will link the ref_seq ID to the name of the instance of the ref_seq_id that was used
        master_name_unique_id_dict = {}

        # we will have the worker of each process put a 'STOP' into its output que when it has finished
        # this way we can count how many 'STOP's we have had output and when this equals the number of processors
        # we started then we know that they have all finished and we have processed all of the outputs
        stop_count = 0
        # the fasta_name_set ouput by the worker can simply be a tuple containing
        # (temp_fasta_dict, temp_name_dict, temp_group_list, list_of_IDs, proc_name)
        for fasta_name_set in iter(outputQueue.get, 'STOP'):
            if fasta_name_set == 'EXIT':
                stop_count += 1
                if stop_count == num_processors:
                    break
            else:
                proc_fasta_dict = fasta_name_set[0]
                proc_name_dict = fasta_name_set[1]
                proc_group_list = fasta_name_set[2]
                proc_seq_list = fasta_name_set[3]
                cc_id = fasta_name_set[4]
                cc_name = fasta_name_set[5]
                sys.stdout.write('\rAdding {} to master fasta and name files'.format(cc_name))
                master_group_list += proc_group_list
                for seq_id in proc_seq_list:
                    seq_name = '{}_id{}_0'.format(seq_id, cc_id)
                    if seq_id not in master_unique_fasta_seq_list:
                        master_unique_fasta_seq_list.append(seq_id)
                        # then this ref seq has not yet been added to the fasta and it needs to be
                        # populate the master fasta
                        master_fasta_dict[seq_name] = proc_fasta_dict[seq_name]
                        # update the master_name_unique_id_dict to keep track of which sequence name represents the ref_seq_id
                        master_name_unique_id_dict[seq_id] = seq_name
                        # create a name dict entry
                        master_name_dict[seq_name] = proc_name_dict[seq_name]
                        # extend the master_group_list
                    else:
                        # then this ref seq is already in the fasta so we just need to update the master name
                        # and the group file
                        # look up what the name file key is for the given seq_id
                        master_name_dict[master_name_unique_id_dict[seq_id]] += proc_name_dict[seq_name]

        # process the outputs of the sub processess before we pause to wait for them to complete.
        for p in allProcesses:
            p.join()

        # Now make the fasta file
        fasta_file = []
        for key, value in master_fasta_dict.items():
            fasta_file.extend(['>{}'.format(key), value])

        # Now make the name file
        name_file = []
        for key, value in master_name_dict.items():
            name_file.append('{}\t{}'.format(key, ','.join(value)))

        # output_file_paths.extend(['{}/unique.fasta'.format(clade_wkd), '{}/name_file.names'.format(clade_wkd), '{}/group_file.groups'.format(clade_wkd)])
        writeListToDestination('{}/unique.fasta'.format(clade_wkd), fasta_file)
        writeListToDestination('{}/name_file.names'.format(clade_wkd), name_file)
        writeListToDestination('{}/group_file.groups'.format(clade_wkd), master_group_list)

        out_file = mafft_align_fasta(clade_wkd, num_proc=num_processors)

        ### I am really struglling to get the jmodeltest to run through python.
        # I will therefore work with a fixed model that has been chosen by running jmodeltest on the command line
        # phyml  -i /tmp/jmodeltest7232692498672152838.phy -d nt -n 1 -b 0 --run_id TPM1uf+I -m 012210 -f m -v e -c 1 --no_memory_check -o tlr -s BEST
        # The production of an ML tree is going to be unrealistic with the larger number of sequences
        # it will simply take too long to compute.
        # Instead I will create an NJ consnsus tree.
        # This will involve using the embassy version of the phylip executables

        # First I will have to create random data sets
        fseqboot_base = generate_fseqboot_alignments(clade_wkd, bootstrap_value, out_file)

        if method == 'mothur':
            unifrac_dist, unifrac_path = mothur_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file,
                                                                    bootstrap_value, num_processors, date_time_string)
        elif method == 'phylip':
            unifrac_dist, unifrac_path = phylip_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file,
                                                                    bootstrap_value, num_processors)

        PCoA_path = generate_PCoA_coords(clade_wkd, unifrac_dist, date_time_string)
        PCoA_path_lists.append(PCoA_path)
        # Delete the tempDataFolder and contents
        file_to_del = '{}/out_seq_boot_reps'.format(clade_wkd)
        shutil.rmtree(path=file_to_del)

        output_file_paths.append(PCoA_path)
        output_file_paths.append(unifrac_path)

        # now delte all files except for the .csv that holds the coords and the .dist that holds the dists
        list_of_dir = os.listdir(clade_wkd)
        for item in list_of_dir:
            if '.csv' not in item and '.dist' not in item:
                os.remove(os.path.join(clade_wkd, item))

    # Print output files
    sys.stdout.write('\n\nBetween sample distances output files:\n')
    for path_of_output_file in output_file_paths:
        print(path_of_output_file)

    return PCoA_path_lists

def generate_within_clade_UniFrac_distances_samples(dataSubmission_str, num_processors,
                                                    method, call_type, date_time_string, bootstrap_value=100, output_dir=None):

    # The call_type argument will be used to determine which setting this method is being called from.
    # if it is being called as part of the initial submission call_type='submission', then we will always be working with a single
    # data_set. In this case we should output to the same folder that the submission results were output
    # to. In the case of being output in a standalone manner call_type='stand_alone' then we may be outputting
    # comparisons from several data_sets. As such we cannot rely on being able to put the ordination results
    # into the initial submissions folder. In this case we will use the directory structure that is already
    # in place which will put it in the ordination folder.
    # TODO
    # we can now also colour the ordination plots according to the meta data of the samples, if this is available.

    '''
    This method will generate a distance matrix between samples
    One for each clade.
    I have been giving some thought as to which level we should be generating these distance matrices on.
    If we do them on the all sequence level, i.e. mixture of clades then we are going to end up with the ordination
    separations based primarily on the presence or absence of clades and this will mask the within clade differences
    Biologically it is more important to be able to tell the difference between within clade different types
    than be able to see if they contain different cladal proportions. After all, SP is about the increased resolution
    So... what I propose is that the researcher do a clade level analysis seperately to look at the differnces between
    clades and then they do ordinations clade by clade.

    Now, with the within-clade analyses we have to chose between ITS2 type profile collection based ordinations
    (i.e. only looking at the sequences contained in a clade_collection_type) or working with all of the sequnces that
    are found within a sample of the given clade. With the former, this will require that the sample has gone through
    an analysis, this is not possible for environmental samples. For coral samples this is of course possible, and we
    would then need to choose about whether you would want to have one clade_collection_type / sample pairing per clade
    collection type or whether you would plot a sample with all of its clade_collection_types plotted. I think the
    latter would be better as this would really be telling the difference between the samples rather than the types
    found within the samples. However, I'm now thinking about the cases where samples are missing one or two DIVs and
    maybe SymPortal hasn't done the best job of resolving the true type they contain. In this case the distance
    should still give a true representation of the similarity of this sample to another sample. Then you might start
    to wonder what the point of the ITS2 type profiles are. Well, they would be to tell us WHAT the types are
    whereas the distance matrix would not be able to do this but give us quantification of similarity. So they would
    in fact be nicely complementary.

    So for the time being, given the above logic, we will work on creating cladally separated distance matrices
    for samples of a given data_set or collection of dataSubmissions. For each sample we will use all sequences
    of the given clade found within the sample to calculate the distances. This will output the distance matrix
    in the outputs folder.
    '''

    output_file_paths = []
    data_submissions = data_set.objects.filter(id__in=[int(a) for a in str(dataSubmission_str).split(',')])

    clade_collection_list_of_dataSubmissions = clade_collection.objects.filter(
        dataSetSampleFrom__dataSubmissionFrom__in=data_submissions)

    clades_of_clade_collections = list(set([a.clade for a in clade_collection_list_of_dataSubmissions]))

    if call_type == 'stand_alone':
        wkd = os.path.abspath(
            os.path.join(os.path.dirname(__file__), 'outputs', 'ordination', '_'.join(dataSubmission_str.split(',')),
                         'between_samples'))
    else:
        # call_type == 'submission':
        wkd = output_dir + '/between_sample_distances'
    # for each clade found in the dataSubmissions' samples
    PCoA_path_lists = []
    for clade_in_question in clades_of_clade_collections:

        clade_wkd = wkd + '/{}'.format(clade_in_question)
        clade_collections_of_clade = clade_collection_list_of_dataSubmissions.filter(clade=clade_in_question)
        if len(clade_collections_of_clade) < 2:
            continue
        # directly make a name file and unique fasta file

        # NB the MP that we will be writing below is complex and cannot be completed entirely within the worker
        # we will out put to the output queue live and process this in the main thread to create the master
        # fasta and name file

        sys.stdout.write('Creating .name and .fasta files')

        # setup MP
        # Create the queues that will hold the cc information
        taskQueue = Queue()

        # output queue to put the dictionaries that will make the fasta and name files
        outputQueue = Queue()

        # populate the input queue
        for cc in clade_collections_of_clade:
            taskQueue.put(cc)

        # place the stop cues
        for n in range(num_processors):
            taskQueue.put('STOP')

        allProcesses = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        for n in range(num_processors):
            p = Process(target=uni_frac_worker_two, args=(taskQueue, outputQueue))
            allProcesses.append(p)
            p.start()

        # have an output pipe that gets processed as we go along. I.e. before we call wait.
        master_fasta_dict = {}
        master_unique_fasta_seq_list = []
        master_name_dict = {}
        master_group_list = []
        # this dict will link the ref_seq ID to the name of the instance of the ref_seq_id that was used
        master_name_unique_id_dict = {}

        # we will have the worker of each process put a 'STOP' into its output que when it has finished
        # this way we can count how many 'STOP's we have had output and when this equals the number of processors
        # we started then we know that they have all finished and we have processed all of the outputs
        stop_count = 0
        # the fasta_name_set ouput by the worker can simply be a tuple containing
        # (temp_fasta_dict, temp_name_dict, temp_group_list, list_of_IDs, proc_name)
        for fasta_name_set in iter(outputQueue.get, 'STOP'):
            if fasta_name_set == 'EXIT':
                stop_count += 1
                if stop_count == num_processors:
                    break
            else:
                proc_fasta_dict = fasta_name_set[0]
                proc_name_dict = fasta_name_set[1]
                proc_group_list = fasta_name_set[2]
                proc_seq_list = fasta_name_set[3]
                cc_id = fasta_name_set[4]
                cc_name = fasta_name_set[5]
                sys.stdout.write('\rAdding {} to master fasta and name files'.format(cc_name))
                master_group_list += proc_group_list
                for seq_id in proc_seq_list:
                    seq_name = '{}_id{}_0'.format(seq_id, cc_id)
                    if seq_id not in master_unique_fasta_seq_list:
                        master_unique_fasta_seq_list.append(seq_id)
                        # then this ref seq has not yet been added to the fasta and it needs to be
                        # populate the master fasta
                        master_fasta_dict[seq_name] = proc_fasta_dict[seq_name]
                        # update the master_name_unique_id_dict to keep track of which sequence name represents the ref_seq_id
                        master_name_unique_id_dict[seq_id] = seq_name
                        # create a name dict entry
                        master_name_dict[seq_name] = proc_name_dict[seq_name]
                        # extend the master_group_list
                    else:
                        # then this ref seq is already in the fasta so we just need to update the master name
                        # and the group file
                        # look up what the name file key is for the given seq_id
                        master_name_dict[master_name_unique_id_dict[seq_id]] += proc_name_dict[seq_name]

        # process the outputs of the sub processess before we pause to wait for them to complete.
        for p in allProcesses:
            p.join()

        # Now make the fasta file
        fasta_file = []
        for key, value in master_fasta_dict.items():
            fasta_file.extend(['>{}'.format(key), value])

        # Now make the name file
        name_file = []
        for key, value in master_name_dict.items():
            name_file.append('{}\t{}'.format(key, ','.join(value)))

        # output_file_paths.extend(['{}/unique.fasta'.format(clade_wkd), '{}/name_file.names'.format(clade_wkd), '{}/group_file.groups'.format(clade_wkd)])
        writeListToDestination('{}/unique.fasta'.format(clade_wkd), fasta_file)
        writeListToDestination('{}/name_file.names'.format(clade_wkd), name_file)
        writeListToDestination('{}/group_file.groups'.format(clade_wkd), master_group_list)

        out_file = mafft_align_fasta(clade_wkd, num_proc=num_processors)

        ### I am really struglling to get the jmodeltest to run through python.
        # I will therefore work with a fixed model that has been chosen by running jmodeltest on the command line
        # phyml  -i /tmp/jmodeltest7232692498672152838.phy -d nt -n 1 -b 0 --run_id TPM1uf+I -m 012210 -f m -v e -c 1 --no_memory_check -o tlr -s BEST
        # The production of an ML tree is going to be unrealistic with the larger number of sequences
        # it will simply take too long to compute.
        # Instead I will create an NJ consnsus tree.
        # This will involve using the embassy version of the phylip executables

        # First I will have to create random data sets
        fseqboot_base = generate_fseqboot_alignments(clade_wkd, bootstrap_value, out_file)

        if method == 'mothur':
            unifrac_dist, unifrac_path = mothur_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file,
                                                                    bootstrap_value, num_processors, date_time_string=date_time_string)
        elif method == 'phylip':
            unifrac_dist, unifrac_path = phylip_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file,
                                                                    bootstrap_value, num_processors)

        PCoA_path = generate_PCoA_coords(clade_wkd, unifrac_dist, date_time_string=date_time_string)
        PCoA_path_lists.append(PCoA_path)
        # Delete the tempDataFolder and contents
        file_to_del = '{}/out_seq_boot_reps'.format(clade_wkd)
        shutil.rmtree(path=file_to_del)

        output_file_paths.append(PCoA_path)
        output_file_paths.append(unifrac_path)

        # now delte all files except for the .csv that holds the coords and the .dist that holds the dists
        list_of_dir = os.listdir(clade_wkd)
        for item in list_of_dir:
            if '.csv' not in item and '.dist' not in item:
                os.remove(os.path.join(clade_wkd, item))

    # Print output files
    sys.stdout.write('\n\nBetween sample distances output files:\n')
    for path_of_output_file in output_file_paths:
        print(path_of_output_file)

    return PCoA_path_lists

def generate_within_clade_BrayCurtis_distances_samples_sample_list_input(smpl_id_list_str, date_time_string):
    # The call_type argument will be used to determine which setting this method is being called from.
    # if it is being called as part of the initial submission call_type='submission', then we will always be working with a single
    # data_set. In this case we should output to the same folder that the submission results were output
    # to. In the case of being output in a standalone manner call_type='stand_alone' then we may be outputting
    # comparisons from several data_sets. As such we cannot rely on being able to put the ordination results
    # into the initial submissions folder. In this case we will use the directory structure that is already
    # in place which will put it in the ordination folder.
    # TODO
    # we can now also colour the ordination plots according to the meta data of the samples, if this is available.

    '''
    This method will generate a distance matrix between samples
    One for each clade.
    I have been giving some thought as to which level we should be generating these distance matrices on.
    If we do them on the all sequence level, i.e. mixture of clades then we are going to end up with the ordination
    separations based primarily on the presence or absence of clades and this will mask the within clade differences
    Biologically it is more important to be able to tell the difference between within clade different types
    than be able to see if they contain different cladal proportions. After all, SP is about the increased resolution
    So... what I propose is that the researcher do a clade level analysis seperately to look at the differnces between
    clades and then they do ordinations clade by clade.

    Now, with the within-clade analyses we have to chose between ITS2 type profile collection based ordinations
    (i.e. only looking at the sequences contained in a clade_collection_type) or working with all of the sequnces that
    are found within a sample of the given clade. With the former, this will require that the sample has gone through
    an analysis, this is not possible for environmental samples. For coral samples this is of course possible, and we
    would then need to choose about whether you would want to have one clade_collection_type / sample pairing per clade
    collection type or whether you would plot a sample with all of its clade_collection_types plotted. I think the
    latter would be better as this would really be telling the difference between the samples rather than the types
    found within the samples. However, I'm now thinking about the cases where samples are missing one or two DIVs and
    maybe SymPortal hasn't done the best job of resolving the true type they contain. In this case the distance
    should still give a true representation of the similarity of this sample to another sample. Then you might start
    to wonder what the point of the ITS2 type profiles are. Well, they would be to tell us WHAT the types are
    whereas the distance matrix would not be able to do this but give us quantification of similarity. So they would
    in fact be nicely complementary.

    So for the time being, given the above logic, we will work on creating cladally separated distance matrices
    for samples of a given data_set or collection of dataSubmissions. For each sample we will use all sequences
    of the given clade found within the sample to calculate the distances. This will output the distance matrix
    in the outputs folder.
    '''

    output_file_paths = []

    smpl_id_list = [int(str_id) for str_id in smpl_id_list_str.split(',')]
    sample_list = data_set_sample.objects.filter(id__in=smpl_id_list)
    clade_collection_list_of_samples = clade_collection.objects.filter(
        dataSetSampleFrom__in=sample_list)

    clades_of_clade_collections = list(set([a.clade for a in clade_collection_list_of_samples]))

    wkd = os.path.abspath(
        os.path.join(os.path.dirname(__file__), 'outputs', 'ordination', 'custom_sample_list',
                     'between_samples', date_time_string))


    # for each clade found in the dataSubmissions' samples
    PCoA_path_lists = []
    for clade_in_question in clades_of_clade_collections:

        clade_wkd = wkd + '/{}'.format(clade_in_question)
        # convert to list so that the order is set
        clade_collections_of_clade = list(clade_collection_list_of_samples.filter(clade=clade_in_question))

        if len(clade_collections_of_clade) < 2:
            continue
        # this is where we should start to work with the bray curtis method
        # first thing to do will be to go through each of the clade collections and create a dict
        # that has key as the actual sequence and relative abundance of that sequence
        # we can then store these dict in a dict where the key is the sample ID.
        data_set_samples_seq_rel_abund_of_clade_cols_dict = {}
        for clade_col in clade_collections_of_clade:
            temp_dict = {}
            data_set_sample_sequences_of_clade_col = data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=clade_col)
            total_seqs_ind_clade_col = sum([dsss.abundance for dsss in data_set_sample_sequences_of_clade_col])
            for dsss in data_set_sample_sequences_of_clade_col:
                temp_dict[dsss.referenceSequenceOf.sequence] = dsss.abundance / total_seqs_ind_clade_col
            data_set_samples_seq_rel_abund_of_clade_cols_dict[clade_col.dataSetSampleFrom.id] = temp_dict



        # then we can simply do a pairwise comparison of the clade collections and create distances
        within_clade_distances_dict = {}
        for clade_col_one, clade_col_two in itertools.combinations(list(clade_collections_of_clade), 2):
            # let's work to a virtual subsample of 100 000
            clade_col_one_seq_rel_abundance_dict = data_set_samples_seq_rel_abund_of_clade_cols_dict[clade_col_one.dataSetSampleFrom.id]
            clade_col_two_seq_rel_abundance_dict = data_set_samples_seq_rel_abund_of_clade_cols_dict[
                clade_col_two.dataSetSampleFrom.id]

            # for each comparison. Get a set of all of the sequence and convert this to a list.
            set_of_sequences = set(list(clade_col_one_seq_rel_abundance_dict.keys()))
            set_of_sequences.update(list(clade_col_two_seq_rel_abundance_dict.keys()))
            list_of_sequences = list(set_of_sequences)


            # then iter through the list to get the rel abundances for each of the samples, putting 0 if not found in
            # the sample.
            seq_abundance_list_one = []
            seq_abundance_list_two = []
            for seq in list_of_sequences:
                # populate the abundance list cc one
                if seq in clade_col_one_seq_rel_abundance_dict:
                    seq_abundance_list_one.append(int(100000 * clade_col_one_seq_rel_abundance_dict[seq]))
                else:
                    seq_abundance_list_one.append(0)

                # populate the abundance list cc two
                if seq in clade_col_two_seq_rel_abundance_dict:
                    seq_abundance_list_two.append(int(100000 * clade_col_two_seq_rel_abundance_dict[seq]))
                else:
                    seq_abundance_list_two.append(0)

            distance = braycurtis(seq_abundance_list_one, seq_abundance_list_two)
            # once you have this we should simply be able to crunch the bray-curtis.
            # these distances can be stored in a dictionary by the 'id1/id2' and 'id2/id1'
            within_clade_distances_dict['{}_{}'.format(clade_col_one.dataSetSampleFrom.id, clade_col_two.dataSetSampleFrom.id)] = distance
            within_clade_distances_dict['{}_{}'.format(clade_col_two.dataSetSampleFrom.id, clade_col_one.dataSetSampleFrom.id)] = distance



        # from this dict we can produce the distance file that can be passed into the generate_PCoA_coords method
        distance_out_file = [len(clade_collections_of_clade)]
        for clade_col_outer in clade_collections_of_clade:
            temp_clade_col_string = [clade_col_outer.dataSetSampleFrom.id]

            for clade_col_inner in clade_collections_of_clade:
                if clade_col_outer == clade_col_inner:
                    temp_clade_col_string.append(0)
                else:
                    temp_clade_col_string.append(within_clade_distances_dict['{}_{}'.format(clade_col_outer.dataSetSampleFrom.id, clade_col_inner.dataSetSampleFrom.id)])
            distance_out_file.append('\t'.join([str(distance_item) for distance_item in temp_clade_col_string ]))
        # from here we can hopefully rely on the rest of the methods as they already are. The .dist file should be
        # written out to the clade_wkd.
        os.makedirs(clade_wkd, exist_ok=True)
        dist_out_path = '{}/{}.bray_curtis_within_clade_sample_distances.dist'.format(clade_wkd, date_time_string)

        # for the output version lets also append the sample name to each line so that we can see which sample it is
        # it is important that we otherwise work eith the sample ID as the sample names may not be unique.
        dist_with_sample_name = [distance_out_file[0]]
        list_of_sample_ids = [int(line.split('\t')[0]) for line in distance_out_file[1:]]
        dss_of_outputs = list(data_set_sample.objects.filter(id__in=list_of_sample_ids))
        dict_of_dss_id_to_name = {dss.id:dss.name for dss in dss_of_outputs}
        for line in distance_out_file[1:]:
            temp_list = []
            sample_id = int(line.split('\t')[0])
            sample_name = dict_of_dss_id_to_name[sample_id]
            temp_list.append(sample_name)
            temp_list.extend(line.split('\t'))
            new_line = '\t'.join(temp_list)
            dist_with_sample_name.append(new_line)

        with open(dist_out_path, 'w') as f:
            for line in dist_with_sample_name:
                f.write('{}\n'.format(line))


        PCoA_path = generate_PCoA_coords(clade_wkd, distance_out_file, date_time_string)
        PCoA_path_lists.append(PCoA_path)
        # Delete the tempDataFolder and contents

        output_file_paths.append(PCoA_path)
        output_file_paths.append(dist_out_path)

    # Print output files
    sys.stdout.write('\n\nBetween sample distances output files:\n')
    for path_of_output_file in output_file_paths:
        print(path_of_output_file)

    return PCoA_path_lists

def generate_within_clade_BrayCurtis_distances_samples(dataSubmission_str, call_type, date_time_str, output_dir=None):
    # The call_type argument will be used to determine which setting this method is being called from.
    # if it is being called as part of the initial submission call_type='submission', then we will always be working with a single
    # data_set. In this case we should output to the same folder that the submission results were output
    # to. In the case of being output in a standalone manner call_type='stand_alone' then we may be outputting
    # comparisons from several data_sets. As such we cannot rely on being able to put the ordination results
    # into the initial submissions folder. In this case we will use the directory structure that is already
    # in place which will put it in the ordination folder.
    # TODO
    # we can now also colour the ordination plots according to the meta data of the samples, if this is available.

    '''
    This method will generate a distance matrix between samples
    One for each clade.
    I have been giving some thought as to which level we should be generating these distance matrices on.
    If we do them on the all sequence level, i.e. mixture of clades then we are going to end up with the ordination
    separations based primarily on the presence or absence of clades and this will mask the within clade differences
    Biologically it is more important to be able to tell the difference between within clade different types
    than be able to see if they contain different cladal proportions. After all, SP is about the increased resolution
    So... what I propose is that the researcher do a clade level analysis seperately to look at the differnces between
    clades and then they do ordinations clade by clade.

    Now, with the within-clade analyses we have to chose between ITS2 type profile collection based ordinations
    (i.e. only looking at the sequences contained in a clade_collection_type) or working with all of the sequnces that
    are found within a sample of the given clade. With the former, this will require that the sample has gone through
    an analysis, this is not possible for environmental samples. For coral samples this is of course possible, and we
    would then need to choose about whether you would want to have one clade_collection_type / sample pairing per clade
    collection type or whether you would plot a sample with all of its clade_collection_types plotted. I think the
    latter would be better as this would really be telling the difference between the samples rather than the types
    found within the samples. However, I'm now thinking about the cases where samples are missing one or two DIVs and
    maybe SymPortal hasn't done the best job of resolving the true type they contain. In this case the distance
    should still give a true representation of the similarity of this sample to another sample. Then you might start
    to wonder what the point of the ITS2 type profiles are. Well, they would be to tell us WHAT the types are
    whereas the distance matrix would not be able to do this but give us quantification of similarity. So they would
    in fact be nicely complementary.

    So for the time being, given the above logic, we will work on creating cladally separated distance matrices
    for samples of a given data_set or collection of dataSubmissions. For each sample we will use all sequences
    of the given clade found within the sample to calculate the distances. This will output the distance matrix
    in the outputs folder.
    '''

    output_file_paths = []
    data_submissions = data_set.objects.filter(id__in=[int(a) for a in str(dataSubmission_str).split(',')])

    clade_collection_list_of_dataSubmissions = clade_collection.objects.filter(
        dataSetSampleFrom__dataSubmissionFrom__in=data_submissions)

    clades_of_clade_collections = list(set([a.clade for a in clade_collection_list_of_dataSubmissions]))

    if call_type == 'stand_alone':
        wkd = os.path.abspath(
            os.path.join(os.path.dirname(__file__), 'outputs', 'ordination', '_'.join(dataSubmission_str.split(',')),
                         'between_samples'))
    else:
        # call_type == 'submission':
        wkd = output_dir + '/between_sample_distances'
    # for each clade found in the dataSubmissions' samples
    PCoA_path_lists = []
    for clade_in_question in clades_of_clade_collections:

        clade_wkd = wkd + '/{}'.format(clade_in_question)
        # convert to list so that the order is set
        clade_collections_of_clade = list(clade_collection_list_of_dataSubmissions.filter(clade=clade_in_question))

        if len(clade_collections_of_clade) < 2:
            continue
        # this is where we should start to work with the bray curtis method
        # first thing to do will be to go through each of the clade collections and create a dict
        # that has key as the actual sequence and relative abundance of that sequence
        # we can then store these dict in a dict where the key is the sample ID.
        data_set_samples_seq_rel_abund_of_clade_cols_dict = {}
        for clade_col in clade_collections_of_clade:
            temp_dict = {}
            data_set_sample_sequences_of_clade_col = data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=clade_col)
            total_seqs_ind_clade_col = sum([dsss.abundance for dsss in data_set_sample_sequences_of_clade_col])
            for dsss in data_set_sample_sequences_of_clade_col:
                temp_dict[dsss.referenceSequenceOf.sequence] = dsss.abundance / total_seqs_ind_clade_col
            data_set_samples_seq_rel_abund_of_clade_cols_dict[clade_col.dataSetSampleFrom.id] = temp_dict



        # then we can simply do a pairwise comparison of the clade collections and create distances
        within_clade_distances_dict = {}
        for clade_col_one, clade_col_two in itertools.combinations(list(clade_collections_of_clade), 2):
            # let's work to a virtual subsample of 100 000
            clade_col_one_seq_rel_abundance_dict = data_set_samples_seq_rel_abund_of_clade_cols_dict[clade_col_one.dataSetSampleFrom.id]
            clade_col_two_seq_rel_abundance_dict = data_set_samples_seq_rel_abund_of_clade_cols_dict[
                clade_col_two.dataSetSampleFrom.id]

            # for each comparison. Get a set of all of the sequence and convert this to a list.
            set_of_sequences = set(list(clade_col_one_seq_rel_abundance_dict.keys()))
            set_of_sequences.update(list(clade_col_two_seq_rel_abundance_dict.keys()))
            list_of_sequences = list(set_of_sequences)


            # then iter through the list to get the rel abundances for each of the samples, putting 0 if not found in
            # the sample.
            seq_abundance_list_one = []
            seq_abundance_list_two = []
            for seq in list_of_sequences:
                # populate the abundance list cc one
                if seq in clade_col_one_seq_rel_abundance_dict:
                    seq_abundance_list_one.append(int(100000 * clade_col_one_seq_rel_abundance_dict[seq]))
                else:
                    seq_abundance_list_one.append(0)

                # populate the abundance list cc two
                if seq in clade_col_two_seq_rel_abundance_dict:
                    seq_abundance_list_two.append(int(100000 * clade_col_two_seq_rel_abundance_dict[seq]))
                else:
                    seq_abundance_list_two.append(0)

            distance = braycurtis(seq_abundance_list_one, seq_abundance_list_two)
            # once you have this we should simply be able to crunch the bray-curtis.
            # these distances can be stored in a dictionary by the 'id1/id2' and 'id2/id1'
            within_clade_distances_dict['{}_{}'.format(clade_col_one.dataSetSampleFrom.id, clade_col_two.dataSetSampleFrom.id)] = distance
            within_clade_distances_dict['{}_{}'.format(clade_col_two.dataSetSampleFrom.id, clade_col_one.dataSetSampleFrom.id)] = distance



        # from this dict we can produce the distance file that can be passed into the generate_PCoA_coords method
        distance_out_file = [len(clade_collections_of_clade)]
        for clade_col_outer in clade_collections_of_clade:
            temp_clade_col_string = [clade_col_outer.dataSetSampleFrom.id]

            for clade_col_inner in clade_collections_of_clade:
                if clade_col_outer == clade_col_inner:
                    temp_clade_col_string.append(0)
                else:
                    temp_clade_col_string.append(within_clade_distances_dict['{}_{}'.format(clade_col_outer.dataSetSampleFrom.id, clade_col_inner.dataSetSampleFrom.id)])
            distance_out_file.append('\t'.join([str(distance_item) for distance_item in temp_clade_col_string ]))
        # from here we can hopefully rely on the rest of the methods as they already are. The .dist file should be
        # written out to the clade_wkd.
        os.makedirs(clade_wkd, exist_ok=True)
        dist_out_path = '{}/{}_bray_curtis_within_clade_sample_distances.dist'.format(clade_wkd, date_time_str)

        # for the output version lets also append the sample name to each line so that we can see which sample it is
        # it is important that we otherwise work eith the sample ID as the sample names may not be unique.
        dist_with_sample_name = [distance_out_file[0]]
        list_of_sample_ids = [int(line.split('\t')[0]) for line in distance_out_file[1:]]
        dss_of_outputs = list(data_set_sample.objects.filter(id__in=list_of_sample_ids))
        dict_of_dss_id_to_name = {dss.id:dss.name for dss in dss_of_outputs}
        for line in distance_out_file[1:]:
            temp_list = []
            sample_id = int(line.split('\t')[0])
            sample_name = dict_of_dss_id_to_name[sample_id]
            temp_list.append(sample_name)
            temp_list.extend(line.split('\t'))
            new_line = '\t'.join(temp_list)
            dist_with_sample_name.append(new_line)

        with open(dist_out_path, 'w') as f:
            for line in dist_with_sample_name:
                f.write('{}\n'.format(line))


        PCoA_path = generate_PCoA_coords(clade_wkd, distance_out_file)
        PCoA_path_lists.append(PCoA_path)
        # Delete the tempDataFolder and contents

        output_file_paths.append(PCoA_path)
        output_file_paths.append(dist_out_path)

    # Print output files
    sys.stdout.write('\n\nBetween sample distances output files:\n')
    for path_of_output_file in output_file_paths:
        print(path_of_output_file)

    return PCoA_path_lists

def generate_within_clade_BrayCurtis_distances_ITS2_type_profiles(data_submission_id_str, data_analysis_id, call_type, date_time_string, output_dir=None):
    ''' This will produce distance matrices between ITS2 type profiles of the same clade.
    It will use exactly the same suite of programs as the generate_within_clade_UniFrac_distances_samples method
    The only difference will be that this will be calculated across ITS2 Type Profiles rather than samples.
    As such, this will require a data_set set and dataAnalysis.
    The average DIV abundances for each of the ITS2 type profiles should first be calculated from which the distances
    can then be calcualated'''

    output_file_paths = []

    if call_type == 'stand_alone':
        wkd = os.path.abspath(os.path.join(os.path.dirname(__file__), 'outputs',
                                           'ordination', '_'.join(str(data_submission_id_str).split(',')),
                                           'between_profiles'))
    else:
        # call_type = 'analysis'
        wkd = '{}/{}'.format(output_dir, 'between_its2_type_profile_distances')

    # get the dataSubmissions
    data_submissions = data_set.objects.filter(id__in=[int(a) for a in str(data_submission_id_str).split(',')])

    # get the dataAnalysis
    data_analysis_obj = data_analysis.objects.get(id=data_analysis_id)

    # go clade by clade
    ITS2_type_profiles_of_data_subs_and_analysis = analysis_type.objects.filter(
        clade_collection_type__cladeCollectionFoundIn__dataSetSampleFrom__dataSubmissionFrom__in=
        data_submissions, dataAnalysisFrom=data_analysis_obj).distinct()

    clade_list = list(set([type_profile.clade for type_profile in ITS2_type_profiles_of_data_subs_and_analysis]))

    # list for storing the paths to the PCoA .csv that will be produced
    PCoA_path_lists = []
    for clade in clade_list:
        # convert query to a list so that we can iterate over twice in exactly the same order
        ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade = list(
            ITS2_type_profiles_of_data_subs_and_analysis.filter(clade=clade))

        if len(ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade) < 2:
            continue

        clade_wkd = wkd + '/{}'.format(clade)
        # we will need to get a list of all of the sequences that are found in the above ITS2 type profiles
        # we will then need to align these.
        # we will also need to produce a group file, a name file and a unique fasta file
        # this will be similar to what we have in the generate_within_clade_UniFrac_distances_samples
        # from this point on it should be easy for us to use the same set up as we have in the other method
        # we already have the ratios of each of the divs
        type_profile_rel_abundance_dicts_dict = {}
        for at in ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade:
            sys.stdout.write('\rProcessing {}'.format(at))
            list_of_div_ids = [int(b) for b in at.orderedFootprintList.split(',')]
            foot_print_ratio_array = pd.DataFrame(at.getRatioList())
            rel_abundance_of_divs_dict = {list_of_div_ids[i]: float(foot_print_ratio_array[i].mean())
                                                 for i in range(len(list_of_div_ids))}
            type_profile_rel_abundance_dicts_dict[at.id] = rel_abundance_of_divs_dict

        # here we have a dict for each of the ITS2 type profiles of the clade
        # we can now do pairwise comparisons just like we did for the samples
        # then we can simply do a pairwise comparison of the clade collections and create distances
        within_clade_distances_dict = {}
        for type_one, type_two in itertools.combinations(list(ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade), 2):
            # let's work to a virtual subsample of 100 000
            clade_col_one_seq_rel_abundance_dict = type_profile_rel_abundance_dicts_dict[
                type_one.id]
            clade_col_two_seq_rel_abundance_dict = type_profile_rel_abundance_dicts_dict[
                type_two.id]

            # for each comparison. Get a set of all of the sequence and convert this to a list.
            set_of_sequences = set(list(clade_col_one_seq_rel_abundance_dict.keys()))
            set_of_sequences.update(list(clade_col_two_seq_rel_abundance_dict.keys()))
            list_of_sequences = list(set_of_sequences)

            # then iter through the list to get the rel abundances for each of the samples, putting 0 if not found in
            # the sample.
            seq_abundance_list_one = []
            seq_abundance_list_two = []
            for seq in list_of_sequences:
                # populate the abundance list cc one
                if seq in clade_col_one_seq_rel_abundance_dict:
                    seq_abundance_list_one.append(int(100000 * clade_col_one_seq_rel_abundance_dict[seq]))
                else:
                    seq_abundance_list_one.append(0)

                # populate the abundance list cc two
                if seq in clade_col_two_seq_rel_abundance_dict:
                    seq_abundance_list_two.append(int(100000 * clade_col_two_seq_rel_abundance_dict[seq]))
                else:
                    seq_abundance_list_two.append(0)

            distance = braycurtis(seq_abundance_list_one, seq_abundance_list_two)
            # once you have this we should simply be able to crunch the bray-curtis.
            # these distances can be stored in a dictionary by the 'id1/id2' and 'id2/id1'
            within_clade_distances_dict[
                '{}_{}'.format(type_one.id, type_two.id)] = distance
            within_clade_distances_dict[
                '{}_{}'.format(type_two.id, type_one.id)] = distance

        # from this dict we can produce the distance file that can be passed into the generate_PCoA_coords method
        distance_out_file = [len(ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade)]
        for type_outer in ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade:
            temp_clade_type_string = [type_outer.name]

            for type_inner in ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade:
                if type_outer == type_inner:
                    temp_clade_type_string.append(0)
                else:
                    temp_clade_type_string.append(within_clade_distances_dict[
                                                     '{}_{}'.format(type_outer.id,
                                                                    type_inner.id)])
            distance_out_file.append('\t'.join([str(distance_item) for distance_item in temp_clade_type_string]))
        # from here we can hopefully rely on the rest of the methods as they already are. The .dist file should be
        # written out to the clade_wkd.
        os.makedirs(clade_wkd, exist_ok=True)
        dist_out_path = '{}/{}.bray_curtis_within_clade_sample_distances.dist'.format(clade_wkd, date_time_string)

        with open(dist_out_path, 'w') as f:
            for line in distance_out_file:
                f.write('{}\n'.format(line))


        PCoA_path = generate_PCoA_coords(clade_wkd, distance_out_file, date_time_string)

        output_file_paths.append(PCoA_path)
        output_file_paths.append(dist_out_path)

        PCoA_path_lists.append(PCoA_path)

    # output the paths of the new files created
    print('\n\nOutput files:\n')
    for path_of_output_file in output_file_paths:
        print(path_of_output_file)

    return PCoA_path_lists


def mafft_align_fasta(clade_wkd, num_proc):
    # now mafft align the fasta
    # now align
    # http://plumbum.readthedocs.io/en/latest/local_commands.html#pipelining
    sys.stdout.write('\rAligning sequences')
    mafft = local["mafft"]
    in_file = '{}/unique.fasta'.format(clade_wkd)
    out_file = '{}/unique.aligned.fasta'.format(clade_wkd)

    # NB we were getting an error here because our os.cwd has been in the temp folder which we then deleted
    # when mafft runs it does a look up to get what our cwd is and this was causing a crash.
    # to fix this we will change directory to the outputdir
    os.chdir(clade_wkd)

    # now run mafft including the redirect
    (mafft['--auto', '--thread', num_proc, in_file] > out_file)()
    return out_file


def uni_frac_worker_two(input, output):
    # We can fix this so that it is more multithreadable.
    # first we can fix the having to have a a numerical_name_dict, by just having a counter per
    # process and append the process id to make sure that the sequences are unique
    # we could have an ouput pipe where we could place individual dictionaries which we could process outside
    # Because one thread can process more than one cc and therefore can come across the same ref_seq_id
    # we can't just name them using the range of the abundance as we will end up with sequnce of the same name
    # we will need to have a dict again for each of the processes that keeps track of the count for each of
    # the particular ref_seq_ids
    proc_id = current_process().name

    test_group = []
    # For each clade collection
    for cc in iter(input.get, 'STOP'):
        # set up clade colection
        temp_fasta_dict = {}
        temp_name_dict = {}
        temp_group_list = []
        temp_ref_seq_id_list = []

        sys.stdout.write('\rProcessing cc: {} with {}'.format(cc, proc_id))
        for data_set_sample_seq in data_set_sample_sequence.objects.filter(
                cladeCollectionTwoFoundIn=cc):
            ref_seq_id = data_set_sample_seq.referenceSequenceOf.id
            temp_ref_seq_id_list.append(ref_seq_id)
            unique_seq_name_base = '{}_id{}'.format(ref_seq_id, cc.id)

            smp_name = str(cc.id)
            temp_fasta_dict['{}_{}'.format(unique_seq_name_base, 0)] = data_set_sample_seq.referenceSequenceOf.sequence
            temp_name_list = []

            for i in range(data_set_sample_seq.abundance):
                temp_name_list.append('{}_{}'.format(unique_seq_name_base, i))
                temp_group_list.append('{}\t{}'.format('{}_{}'.format(unique_seq_name_base, i), smp_name))

            temp_name_dict['{}_{}'.format(unique_seq_name_base, 0)] = temp_name_list
        output.put((temp_fasta_dict, temp_name_dict, temp_group_list, temp_ref_seq_id_list, str(cc.id), str(cc)))
    output.put('EXIT')


def concatenate_trees(list_of_tree_paths, out_put_path):
    master_tree_file = []
    for i in range(len(list_of_tree_paths)):
        temp_tree_file = readDefinedFileToList(list_of_tree_paths[i])
        for line in temp_tree_file:
            master_tree_file.append(line)
    writeListToDestination(out_put_path, master_tree_file)


def create_consesnus_tree(clade_wkd, list_of_tree_paths, name_file):
    # Now create consensus tree using sumtrees.py
    concatenate_trees(list_of_tree_paths, '{}/individual_trees'.format(clade_wkd))
    in_file_fconsense = '{}/individual_trees'.format(clade_wkd)

    # # ### DEBUG ###
    # for line in readDefinedFileToList(in_file_fconsense):
    #     if "1 Pro" in line:
    #         bob = 'asdf'
    # #############

    tree_out_file_fconsense_sumtrees = '{}/consensus_tree_sumtrees.newick'.format(clade_wkd)

    completed_consensus = subprocess.run(
        ['sumtrees.py', '-F', 'newick', '--replace', '-o', tree_out_file_fconsense_sumtrees, in_file_fconsense],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if completed_consensus.returncode != 0:
        try:
            if 'recursion' in completed_consensus.stdout.decode('utf-8'):
                sys.exit('There has been a recursion depth error whilst trying to calculate a consensus tree\n'
                         'This occured whilst calculating between sample distances.\n'
                         'This problem occurs when trees are too big for sumtrees.py to process.\n'
                         'This problem can often be solved by increasing the recursion limit used when running sumtrees.py\n'
                         'The dafault value is 1000. This can be increased.\n'
                         'To do this you need to edit the sumtrees.py file.\n'
                         'To locate this file, try typing:\n'
                         ' \'which sumtrees.py\'\n'
                         'This should return the location of the sumtrees.py file\n'
                         'After the line:\n'
                         ' \'import json\'  '
                         '\non a new line add: '
                         '\n\'import sys\' '
                         '\nand then on the next line:\n'
                         'sys.setrecursionlimit(1500)\n'
                         'Save changes and rerun.')

        except:
            sys.exit('There has likely been a recursion depth error whilst trying to calculate a consensus tree\n'
                     'This occured whilst calculating between sample distances.\n'
                     'This problem occurs when trees are too big for sumtrees.py to process.\n'
                     'This problem can often be solved by increasing the recursion limit used when running sumtrees.py\n'
                     'The dafault value is 1000. This can be increased.\n'
                     'To do this you need to edit the sumtrees.py file.\n'
                     'To locate this file, try typing:\n'
                     ' \'which sumtrees.py\'\n'
                     'This should return the location of the sumtrees.py file\n'
                     'After the line:\n'
                     ' \'import json\'  '
                     '\non a new line add: '
                     '\n\'import sys\' '
                     '\nand then on the next line:\n'
                     'sys.setrecursionlimit(1500)\n'
                     'Save changes and rerun.')


    # The consunsus tree output by sumtree is causing problems because it contains metadata.
    # It is important that the tree we use for the unifrac has the branch lengths on it
    # The tree output by consense only has the boostrap values instead of lengths which is not ideal
    # I am planning on removing the metadatafrom the sumtree consensus tree and seeing if that works in mothur's
    # Unifrac. This can be part of the rename_tree code
    rename_tree_two(name_file, tree_out_file_fconsense_sumtrees)
    return tree_out_file_fconsense_sumtrees


def generate_PCoA_coords(clade_wkd, raw_dist_file, date_time_string):
    # simultaneously grab the sample names in the order of the distance matrix and put the matrix into
    # a twoD list and then convert to a numpy array
    temp_two_D_list = []
    sample_names_from_dist_matrix = []
    for line in raw_dist_file[1:]:
        temp_elements = line.split('\t')
        sample_names_from_dist_matrix.append(temp_elements[0].replace(' ', ''))
        temp_two_D_list.append([float(a) for a in temp_elements[1:]])
    uni_frac_dist_array = np.array(temp_two_D_list)
    sys.stdout.write('\rcalculating PCoA coordinates')
    pcoA_full_path = clade_wkd + '/{}.PCoA_coords.csv'.format(date_time_string)
    this = pcoa(uni_frac_dist_array)

    # rename the dataframe index as the sample names
    mapper_dict = {i: j for i, j in enumerate(sample_names_from_dist_matrix)}
    this.samples['sample'] = sample_names_from_dist_matrix
    renamed_dataframe = this.samples.set_index('sample')

    # now add the variance explained as a final row to the renamed_dataframe

    renamed_dataframe = renamed_dataframe.append(this.proportion_explained.rename('proportion_explained'))

    renamed_dataframe.to_csv(pcoA_full_path, index=True, header=True, sep=',')
    return pcoA_full_path


def generate_fseqboot_alignments(clade_wkd, num_reps, out_file):
    # setup fseqboot arguments
    in_file_seqboot = out_file
    out_file_seqboot = in_file_seqboot + '.fseqboot'
    # give the option to install the new phylip suite from their executables
    # or simply download the executables form us and install into the lib/phylipnew/folder
    is_installed = subprocess.call(['which', 'fseqboot'])
    if is_installed == 0:
        fseqboot = local["fseqboot"]
    else:
        fseqboot_path = os.path.join(os.path.dirname(__file__), 'lib/phylipnew/fseqboot')
        if os.path.isfile(fseqboot_path):
            fseqboot = local[fseqboot_path]
        else:
            sys.exit('Cannot find fseqboot in PATH or in local installation at ./lib/phylipnew/fseqboot\n'
                     'For instructions on installing the phylipnew dependencies please visit the SymPortal'
                     'GitHub page: https://github.com/didillysquat/SymPortal_framework/wiki/SymPortal-setup#6-third-party-dependencies')
    # run fseqboot
    sys.stdout.write('\rGenerating multiple datasets')
    (fseqboot['-sequence', in_file_seqboot, '-outfile', out_file_seqboot, '-test', 'b', '-reps', num_reps])()
    # Now divide the fseqboot file up into its 100 different alignments
    fseqboot_file = readDefinedFileToList(out_file_seqboot)
    rep_count = 0
    out_put_holder = []
    out_put_holder.append(fseqboot_file[0])
    fseqboot_base = '{}/out_seq_boot_reps/fseqboot_rep_'.format(clade_wkd)

    # NB below we used to just look for lines that when split had a length of two but this was not specific enough
    # this was causing our system to crash. When we had shorted sequences for example clade A sequences some of the
    # lines only had two groups of 10 bp on them and thus would return a match to line.split() ==2.
    # to fix this we will implent regular expression solution instead
    reg_ex = re.compile('[0-9]+$')
    for line in fseqboot_file[1:]:
        reg_ex_matches_list = reg_ex.findall(line)
        if len(reg_ex_matches_list) == 1:
            writeListToDestination('{}{}'.format(fseqboot_base, rep_count), out_put_holder)
            out_put_holder = []
            out_put_holder.append(line)
            rep_count += 1
        else:
            out_put_holder.append(line)
    writeListToDestination('{}{}'.format(fseqboot_base, rep_count), out_put_holder)
    return fseqboot_base


def perform_unifrac(clade_wkd, tree_path):
    group_file_path = '{}/group_file.groups'.format(clade_wkd)
    name_file_path = '{}/name_file.names'.format(clade_wkd)
    mothur_batch_WU = \
        ['set.dir(input={}, output={})'.format(clade_wkd, clade_wkd),
         'unifrac.weighted(tree={}, group={}, name={}, distance=square, processors=16)'
             .format(tree_path, group_file_path, name_file_path)]
    mothur_batch_path = '{}/mothur_batch_WU'.format(clade_wkd)
    writeListToDestination(mothur_batch_path, mothur_batch_WU)
    # now run the batch file with mothur
    sys.stdout.write('\rcalculating unifrac distances')
    completedProcess = \
        subprocess.run(['mothur', '{}'.format(mothur_batch_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if completedProcess.returncode == 0:
        sys.stdout.write('\rUnifrac successful')
    else:
        sys.stdout.write('\rERROR: {}'.format(completedProcess.sterr.decode('utf-8')))
    return


def rename_tree_two(name_file, tree_out_file_fconsense):
    name_file_reps = []
    for line in name_file:
        name_file_reps.append(line.split('\t')[0])
    seq_re = re.compile('\d+[_ ]id[\d_]+')
    seq_re_meta = re.compile('\[[^\]]*\]')
    sys.stdout.write('\rrenaming tree nodes')
    tree_file = readDefinedFileToList(tree_out_file_fconsense)
    new_tree_file = []
    for line in tree_file:
        # # ### DEBUG ###
        # ### Work out which process name is causing the issue
        #
        #
        # print('processing line {}'.format(tree_file.index(line)))
        # if "1 P" in line or "'4" in line:
        #     buzz = 'asdf'
        #
        #     for i in range(len(line)):
        #         if i < len(line) -1 :
        #             if line[i] + line[i + 1] == "'4":
        #                 if i < 500:
        #                     bobby = line[:i + 500]
        #                 else:
        #                     bobby = line[i-500:i+500]
        # #############
        new_str = line
        examine = list(seq_re.findall(line))

        ### DEBUG ###
        diff = list(set(name_file_reps) - set(examine))
        #############

        # N.B. the sumtrees.py program was causing some very strange behaviour. It was converting '_' to ' '
        # when they were preceeded by a single digit but leaving them as '_' when there were multiple digits before it
        # it took a long time to find what the problem was. It was causing issues in the mothur UniFrac.
        # You will see that I have modified below by having the space function and replacing ' ' with '_' and modifying
        # the regex.
        name_file_rep_match_list = []
        for match_str in examine:
            space = False
            if ' ' in match_str:
                space = True
            if space:

                found = False
                match_str_replced_space = match_str.replace(' ', '_')
                for name_file_rep in name_file_reps:
                    if name_file_rep.startswith(match_str_replced_space):
                        new_str = re.sub("'{}'".format(match_str), name_file_rep, new_str)
                        name_file_rep_match_list.append(name_file_rep)
                        found = True
                        break
            else:
                found = False
                for name_file_rep in name_file_reps:
                    if name_file_rep.startswith(match_str):
                        # then this is a match
                        # now replace the string in the tree file
                        new_str = re.sub('(?<!\d){}'.format(match_str), name_file_rep, new_str)
                        name_file_rep_match_list.append(name_file_rep)
                        found = True
                        break
            if found == False:
                bobbie = 'asdf'

        # now also remove the metadata which is held between square brackets '[]'
        examine_two = list(seq_re_meta.findall(new_str))

        diff = list(set(name_file_reps) - set(name_file_rep_match_list))

        new_str = re.sub('\[[^\]]*\]', '', new_str)

        new_tree_file.append(new_str)

    # here all of the tree_file names should have been replaced. Now write back out.
    sys.stdout.write('\rwriting out tree')
    writeListToDestination(tree_out_file_fconsense, new_tree_file)



# mothur method
def mothur_unifrac_pipeline_MP_worker(input, output, fseqboot_base, clade_wkd):
    for p in iter(input.get, 'STOP'):
        sys.stdout.write('\rProcessing p={} with {}'.format(p, current_process().name))
        # convert the interleaved fasta to sequential fasta
        interleaved_fast = readDefinedFileToList('{}{}'.format(fseqboot_base, p))
        sequen_fast = convert_interleaved_to_sequencial_fasta(interleaved_fast)
        writeListToDestination('{}{}.sequential.fasta'.format(fseqboot_base, p), sequen_fast)
        mothur_batch_dist = \
            ['set.dir(input={}/out_seq_boot_reps/, output={}/out_seq_boot_reps/)'.format(clade_wkd, clade_wkd),
             'dist.seqs(fasta={}, countends=T, output=square)'
                 .format('{}{}.sequential.fasta'.format(fseqboot_base, p))]
        mothur_batch_path = '{}/out_seq_boot_reps/mothur_batch_batch_dist_{}'.format(clade_wkd, p)

        writeListToDestination(mothur_batch_path, mothur_batch_dist)

        # now run the batch file with mothur
        sys.stdout.write('\rCalculating distances...')
        completedProcess = \
            subprocess.run(['mothur', '{}'.format(mothur_batch_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sys.stdout.write('\rDone.')
        # now run in clearcut
        input = '{}{}.sequential.square.dist'.format(fseqboot_base, p)
        mothur_batch_clearcut = \
            ['set.dir(input={}/out_seq_boot_reps/, output={}/out_seq_boot_reps/)'.format(clade_wkd, clade_wkd),
             'clearcut(phylip={}, verbose=t)'
                 .format(input)]

        mothur_batch_path = '{}/out_seq_boot_reps/mothur_batch_batch_clearcut_{}'.format(clade_wkd, p)
        writeListToDestination(mothur_batch_path, mothur_batch_clearcut)
        sys.stdout.write('\rGenerating NJ tree from distances')
        completedProcess = \
            subprocess.run(['mothur', '{}'.format(mothur_batch_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sys.stdout.write('\rDone')
        output.put(input.replace('.dist', '.tre'))
    output.put('kill')


def mothur_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file, num_reps, num_proc, date_time_string):
    # setup MP
    # Create the queues that will hold the cc information
    taskQueue = Queue()

    # output queue to put the dictionaries that will make the fasta and name files
    outputQueue = Queue()

    # populate the input queue
    for p in range(num_reps):
        taskQueue.put(p)

    # place the stop cues
    for n in range(num_proc):
        taskQueue.put('STOP')

    allProcesses = []

    for n in range(num_proc):
        p = Process(target=mothur_unifrac_pipeline_MP_worker, args=(taskQueue, outputQueue, fseqboot_base, clade_wkd))
        allProcesses.append(p)
        p.start()

    # Get list of tree paths from the output queue
    kill_num = 0
    list_of_tree_paths = []
    while 1:
        passedElement = outputQueue.get()
        if passedElement == 'kill':
            kill_num += 1
            if kill_num == num_proc:
                break
        else:
            list_of_tree_paths.append(passedElement)

    for p in allProcesses:
        p.join()

    # Create consensus
    tree_out_file_fconsense_sumtrees = create_consesnus_tree(clade_wkd, list_of_tree_paths, name_file)
    # Perform UniFrac
    perform_unifrac(clade_wkd, tree_out_file_fconsense_sumtrees)

    dist_file_path = '{}/{}'.format(clade_wkd,
                                    tree_out_file_fconsense_sumtrees.split('/')[-1]) + '1.weighted.phylip.dist'

    # here add a date_time_string element to it to make it unique
    dist_file_path_dts = dist_file_path.replace('consensus_tree_sumtrees', '{}.consensus_tree_sumtrees'.format(date_time_string)).replace('1.weighted.phylip', '')

    subprocess.run(['mv', dist_file_path, dist_file_path_dts])

    raw_dist_file = readDefinedFileToList(dist_file_path_dts)
    return raw_dist_file, dist_file_path_dts



# phylip method
def phylip_unifrac_pipeline(clade_wkd, fseqboot_base, name_file, num_reps):
    ### PHYLIP METHOD ###
    # give the option to install the new phylip suite from their executables
    # or simply download the executables form us and install into the ./lib/phylipnew folder
    is_installed = subprocess.call(['which', 'fdnadist'])
    if is_installed == 0:
        fdnadist = local["fdnadist"]
    else:
        fdnadist_path = os.path.join(os.path.dirname(__file__), 'lib/phylipnew/fdnadist')
        if os.path.isfile(fdnadist_path):
            fdnadist = local[fdnadist_path]
        else:
            sys.exit('Cannot find fdnadist in PATH or in local installation at ./lib/phylipnew/fdnadist\n'
                     'For instructions on installing the phylipnew dependencies please visit the SymPortal'
                     'GitHub page: https://github.com/didillysquat/SymPortal_framework/wiki/SymPortal-setup#6-third-party-dependencies')

    for p in range(num_reps):
        # run dnadist
        in_file_dnadist_rep = '{}{}'.format(fseqboot_base, p)
        out_file_dnadist_rep = '{}out_{}'.format(fseqboot_base, p)
        print('calculating distances rep {}'.format(p))
        (fdnadist['-sequence', in_file_dnadist_rep, '-outfile', out_file_dnadist_rep, '-method', 'j'])()

    # give the option to install the new phylip suite from their executables
    # or simply download the executables form us and install into the ./lib/phylipnew folder
    is_installed = subprocess.call(['which', 'fneighbor'])
    if is_installed == 0:
        fneighbor = local["fneighbor"]
    else:
        fneighbor_path = os.path.join(os.path.dirname(__file__), 'lib/phylipnew/fneighbor')
        if os.path.isfile(fneighbor_path):
            fneighbor = local[fneighbor_path]
        else:
            sys.exit('Cannot find fneighbor in PATH or in local installation at ./lib/phylipnew/fneighbor\n'
                     'For instructions on installing the phylipnew dependencies please visit the SymPortal'
                     'GitHub page: https://github.com/didillysquat/SymPortal_framework/wiki/SymPortal-setup#6-third-party-dependencies')


    list_of_tree_paths = []
    for p in range(num_reps):
        print('generating trees rep {}'.format(p))
        fneighbor_in_file_rep = '{}out_{}'.format(fseqboot_base, p)
        fneighbor_out_file_rep = '{}.fneighbor_{}'.format(fneighbor_in_file_rep, p)
        out_tree_file = '{}.nj_tree_{}'.format(fneighbor_in_file_rep, p)
        (fneighbor[
            '-datafile', fneighbor_in_file_rep, '-outfile', fneighbor_out_file_rep, '-jumble', 'Y', '-outtreefile', out_tree_file])()
        list_of_tree_paths.append(out_tree_file)

    tree_out_file_fconsense_sumtrees = create_consesnus_tree(clade_wkd, list_of_tree_paths, name_file)

    perform_unifrac(clade_wkd, tree_out_file_fconsense_sumtrees)

    dist_file_path = '{}/{}'.format(clade_wkd,
                                    tree_out_file_fconsense_sumtrees.split('/')[-1]) + '1.weighted.phylip.dist'

    raw_dist_file = readDefinedFileToList(dist_file_path)

    return raw_dist_file


def phylip_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file, num_reps, num_proc):
    # setup MP
    # Create the queues that will hold the cc information
    taskQueue = Queue()

    # output queue to put the dictionaries that will make the fasta and name files
    outputQueue = Queue()

    # populate the input queue
    for p in range(num_reps):
        taskQueue.put(p)

    # place the stop cues
    for n in range(num_proc):
        taskQueue.put('STOP')

    allProcesses = []

    for n in range(num_proc):
        p = Process(target=phylip_unifrac_pipeline_MP_worker, args=(taskQueue, outputQueue, fseqboot_base))
        allProcesses.append(p)
        p.start()

    # Get list of tree paths from the output queue
    kill_num = 0
    list_of_tree_paths = []
    while 1:
        passedElement = outputQueue.get()
        if passedElement == 'kill':
            kill_num += 1
            if kill_num == num_proc:
                break
        else:
            list_of_tree_paths.append(passedElement)

    for p in allProcesses:
        p.join()

    tree_out_file_fconsense_sumtrees = create_consesnus_tree(clade_wkd, list_of_tree_paths, name_file)

    perform_unifrac(clade_wkd, tree_out_file_fconsense_sumtrees)

    dist_file_path = '{}/{}'.format(clade_wkd,
                                    tree_out_file_fconsense_sumtrees.split('/')[-1]) + '1.weighted.phylip.dist'

    raw_dist_file = readDefinedFileToList(dist_file_path)

    return raw_dist_file, dist_file_path


def phylip_unifrac_pipeline_MP_worker(input, output, fseqboot_base):
    # give the option to install the new phylip suite from their executables
    # or simply download the executables form us and install into the ./lib/phylipnew folder
    is_installed = subprocess.call(['which', 'fdnadist'])
    if is_installed == 0:
        fdnadist = local["fdnadist"]
    else:
        fdnadist_path = os.path.join(os.path.dirname(__file__), 'lib/phylipnew/fdnadist')
        if os.path.isfile(fdnadist_path):
            fdnadist = local[fdnadist_path]
        else:
            sys.exit('Cannot find fdnadist in PATH or in local installation at ./lib/phylipnew/fdnadist\n'
                     'For instructions on installing the phylipnew dependencies please visit the SymPortal'
                     'GitHub page: https://github.com/didillysquat/SymPortal_framework/wiki/SymPortal-setup#6-third-party-dependencies')

    # give the option to install the new phylip suite from their executables
    # or simply download the executables form us and install into the ./lib/phylipnew folder
    is_installed = subprocess.call(['which', 'fneighbor'])
    if is_installed == 0:
        fneighbor = local["fneighbor"]
    else:
        fneighbor_path = os.path.join(os.path.dirname(__file__), 'lib/phylipnew/fneighbor')
        if os.path.isfile(fneighbor_path):
            fneighbor = local[fneighbor_path]
        else:
            sys.exit('Cannot find fneighbor in PATH or in local installation at ./lib/phylipnew/fneighbor\n'
                     'For instructions on installing the phylipnew dependencies please visit the SymPortal'
                     'GitHub page: https://github.com/didillysquat/SymPortal_framework/wiki/SymPortal-setup#6-third-party-dependencies')

    for p in iter(input.get, 'STOP'):
        # run dnadist
        in_file_dnadist_rep = '{}{}'.format(fseqboot_base, p)
        out_file_dnadist_rep = '{}out_{}'.format(fseqboot_base, p)
        print('calculating distances rep {}'.format(p))
        (fdnadist['-sequence', in_file_dnadist_rep, '-outfile', out_file_dnadist_rep, '-method', 'j'])()

        print('generating trees rep {}'.format(p))
        fneighbor_in_file_rep = '{}out_{}'.format(fseqboot_base, p)
        fneighbor_out_file_rep = '{}.fneighbor_{}'.format(fneighbor_in_file_rep, p)
        out_tree_file = '{}.nj_tree_{}'.format(fneighbor_in_file_rep, p)
        (fneighbor[
            '-datafile', fneighbor_in_file_rep, '-outfile', fneighbor_out_file_rep, '-jumble', 'Y', '-outtreefile', out_tree_file])()
        output.put(out_tree_file)
    output.put('kill')

#############################