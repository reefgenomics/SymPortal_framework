from __future__ import division
from dbApp.models import data_set, reference_sequence, data_set_sample_sequence, analysis_type, analysis_group, data_set_sample, data_analysis, clade_collection, clade_collection_type
# import math
import itertools
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from sklearn import manifold
import operator
import subprocess
import os
import re
import json
import string
from collections import defaultdict, Counter
import statistics
import timeit
from multiprocessing import Queue, Process, Manager, current_process
from django import db
import texttable
import pickle
from ctypes import c_int
import shutil
from colorama import init # Fore, Back, Style,
from django.core.exceptions import ObjectDoesNotExist
from datetime import datetime
import pandas as pd
from plumbum import local
from Bio import SeqIO, Nexus
from Bio.Alphabet import IUPAC
from skbio.stats.ordination import pcoa
import sys
import dendropy
import math



if 'PYCHARM_HOSTED' in os.environ:
    convert = False  # in PyCharm, we should disable convert
    strip = False
    print("Hi! You are using PyCharm")
else:
    convert = None
    strip = None

init(convert=convert, strip=strip, autoreset=True)

###### Generic functions ######
def readDefinedFileToList(filename):
    temp_list = []
    with open(filename, mode='r') as reader:
        temp_list = [line.rstrip() for line in reader]
    return temp_list

def writeListToDestination(destination, listToWrite):
    #print('Writing list to ' + destination)
    try:
        os.makedirs(os.path.dirname(destination))
    except FileExistsError:
        pass

    with open(destination, mode='w') as writer:
        i = 0
        while i < len(listToWrite):
            if i != len(listToWrite)-1:
                writer.write(listToWrite[i] + '\n')
            elif i == len(listToWrite)-1:
                writer.write(listToWrite[i])
            i += 1

def convert_interleaved_to_sequencial_fasta(fasta_in):
    fasta_dict = {}
    list_seq_names = []
    list_seq_sequences = []
    num_seqs = int(fasta_in[0].split()[0])
    fasta_cropped = []
    # Get rid of the first line and get rid of the blank lines
    for line in fasta_in[1:]:
        if line != '':
            fasta_cropped.append(line)

    for i in range(len(fasta_cropped)):
            if i < num_seqs:
                # Then we are on one of the inital lines
                list_seq_names.append(fasta_cropped[i].split()[0])
                list_seq_sequences.append(''.join(fasta_cropped[i].split()[1:]))
            else:
                index = i%num_seqs
                list_seq_sequences[index] += ''.join(fasta_cropped[i].split()[1:])

    out_fasta = []
    for name, seq in zip(list_seq_names, list_seq_sequences):
        out_fasta.extend(['>{}'.format(name), seq])

    return out_fasta

###############################


##### DISTANCE MATRICES #####
def generate_within_clade_UniFrac_distances_ITS2_type_profiles(data_submission_id_str, num_processors, data_analysis_id, method, bootstrap_value):
    ''' This will produce distance matrices between ITS2 type profiles of the same clade.
    It will use exactly the same suite of programs as the generate_within_clade_UniFrac_distances_samples method
    The only difference will be that this will be calculated across ITS2 Type Profiles rather than samples.
    As such, this will require a data_set set and dataAnalysis.
    The average DIV abundances for each of the ITS2 type profiles should first be calculated from which the distances
    can then be calcualated'''

    wkd = os.path.abspath(os.path.join(os.path.dirname(__file__), 'outputs',
                                       'ordination','_'.join(data_submission_id_str.split(',')),
                                       'between_profiles', method))

    # get the dataSubmissions
    data_submissions = data_set.objects.filter(id__in=[int(a) for a in data_submission_id_str.split(',')])

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
        ITS2_type_profiles_of_data_subs_and_analysis_list_of_clade = list(ITS2_type_profiles_of_data_subs_and_analysis.filter(clade=clade))

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
        out_file = mafft_align_fasta(clade_wkd)

        # Generate random data sets
        fseqboot_base = generate_fseqboot_alignments(clade_wkd, bootstrap_value, out_file)

        if method == 'fasttree':
            unifrac_dist = run_fast_tree_unifrac_pipeline(clade_wkd, fseqboot_base, name_file, bootstrap_value)
        elif method == 'mothur':
            unifrac_dist = mothur_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file, bootstrap_value, num_processors)
        elif method == 'phylip':
            unifrac_dist = phylip_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file, bootstrap_value, num_processors)

        PCoA_path = generate_PCoA_coords(clade_wkd, unifrac_dist)

        PCoA_path_lists.append(PCoA_path)
        # Delete the tempDataFolder and contents
        file_to_del = '{}/out_seq_boot_reps'.format(clade_wkd)
        shutil.rmtree(path=file_to_del)

        # now delte mothur .logfiles files
        list_of_dir = os.listdir(clade_wkd)
        for item in list_of_dir:
            if item.endswith(".logfile"):
                os.remove(os.path.join(clade_wkd, item))
    apples = 'asdf'
    # get list

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
        print('Processing {}'.format(at))
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
                    group_file.append('{}\t{}'.format(integer_seq_name, str(at).replace('-','_').replace('/','_')))
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
                    group_file.append('{}\t{}'.format(integer_seq_name, str(at).replace('-','_').replace('/','_')))
                name_file_dict[zero_seq_name] = temp_name_list
                ref_seq_id_to_name_file_representative[ref_seq_id_key] = zero_seq_name

    # now assemble the name file
    name_file = []
    for k, v in name_file_dict.items():
        name_file.append('{}\t{}'.format(k, ','.join(v)))

    return group_file, name_file, unique_fasta


def generate_within_clade_UniFrac_distances_samples(dataSubmission_str, num_processors, method, bootstrap_value):
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


    data_submissions = data_set.objects.filter(id__in=[int(a) for a in dataSubmission_str.split(',')])

    clade_collection_list_of_dataSubmissions = clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom__in=data_submissions)

    clades_of_clade_collections = list(set([a.clade for a in clade_collection_list_of_dataSubmissions]))

    wkd = os.path.abspath(os.path.join(os.path.dirname(__file__), 'outputs', 'ordination', '_'.join(dataSubmission_str.split(',')), 'between_samples', method ))
    # for each clade found in the dataSubmissions' samples
    PCoA_path_lists = []
    for clade_in_question in clades_of_clade_collections:
        # ### DEBUG ###
        # if clade_in_question != 'C':
        #     continue
        # ### DEBUG ###
        clade_wkd = wkd + '/{}'.format(clade_in_question)
        clade_collections_of_clade = clade_collection_list_of_dataSubmissions.filter(clade=clade_in_question)
        if len(clade_collections_of_clade) < 2:
            continue
        # directly make a name file and unique fasta file

        #NB the MP that we will be writing below is complex and cannot be completed entirely within the worker
        # we will out put to the output queue live and process this in the main thread to create the master
        # fasta and name file

        print('Creating .name and .fasta files')

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
                print('Adding {} to master fasta and name files'.format(cc_name))
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



        writeListToDestination('{}/unique.fasta'.format(clade_wkd), fasta_file)
        writeListToDestination('{}/name_file.names'.format(clade_wkd), name_file)
        writeListToDestination('{}/group_file.groups'.format(clade_wkd), master_group_list)

        out_file = mafft_align_fasta(clade_wkd)

        ### I am really struglling to get the jmodeltest to run through python.
        # I will therefore work with a fixed model that has been chosen by running jmodeltest on the command line
        # phyml  -i /tmp/jmodeltest7232692498672152838.phy -d nt -n 1 -b 0 --run_id TPM1uf+I -m 012210 -f m -v e -c 1 --no_memory_check -o tlr -s BEST
        # The production of an ML tree is going to be unrealistic with the larger number of sequences
        # it will simply take too long to compute.
        # Instead I will create an NJ consnsus tree.
        # This will involve using the embassy version of the phylip executables



        # First I will have to create random data sets
        fseqboot_base = generate_fseqboot_alignments(clade_wkd, bootstrap_value, out_file)


        if method == 'fasttree':
            unifrac_dist = run_fast_tree_unifrac_pipeline(clade_wkd, fseqboot_base, name_file, bootstrap_value)
        elif method == 'mothur':
            unifrac_dist = mothur_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file, bootstrap_value, num_processors)
        elif method == 'phylip':
            unifrac_dist = phylip_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file, bootstrap_value, num_processors)

        PCoA_path = generate_PCoA_coords(clade_wkd, unifrac_dist)
        PCoA_path_lists.append(PCoA_path)
        # Delete the tempDataFolder and contents
        file_to_del = '{}/out_seq_boot_reps'.format(clade_wkd)
        shutil.rmtree(path=file_to_del)

        # now delte mothur .logfiles files
        list_of_dir = os.listdir(clade_wkd)
        for item in list_of_dir:
            if item.endswith(".logfile"):
                os.remove(os.path.join(clade_wkd, item))
    return PCoA_path_lists


def mafft_align_fasta(clade_wkd):
    # now mafft align the fasta
    # now align
    # http://plumbum.readthedocs.io/en/latest/local_commands.html#pipelining
    print('Aligning sequences')
    mafft = local["mafft"]
    in_file = '{}/unique.fasta'.format(clade_wkd)
    out_file = '{}/unique.aligned.fasta'.format(clade_wkd)
    # now run mafft including the redirect
    (mafft['--auto', '--thread', '16', in_file] > out_file)()
    return out_file

# delete uni_frac_worker
def uni_frac_worker(input, output):
    # We can fix this so that it is more multithreadable.
    # first we can fix the having to have a a numerical_name_dict, by just having a counter per
    # process and append the process id to make sure that the sequences are unique
    # we could have an ouput pipe where we could place individual dictionaries which we could process outside
    # Because one thread can process more than one cc and therefore can come across the same ref_seq_id
    # we can't just name them using the range of the abundance as we will end up with sequnce of the same name
    # we will need to have a dict again for each of the processes that keeps track of the count for each of
    # the particular ref_seq_ids
    proc_id = current_process().name
    ref_seq_to_count_dict = {}
    test_group = []
    for cc in iter(input.get, 'STOP'):

        temp_fasta_dict = {}
        temp_name_dict = {}
        temp_group_list = []
        temp_ref_seq_id_list = []
        print('Processing cc: {} with {}'.format(cc, proc_id))
        for data_set_sample_seq in data_set_sample_sequence.objects.filter(
                cladeCollectionTwoFoundIn=cc):
            ref_seq_id = data_set_sample_seq.referenceSequenceOf.id
            if ref_seq_id in ref_seq_to_count_dict.keys():
                count = ref_seq_to_count_dict[ref_seq_id]
                ref_seq_to_count_dict[ref_seq_id] += data_set_sample_seq.abundance
                new_count = ref_seq_to_count_dict[ref_seq_id]
            else:
                count = 0
                ref_seq_to_count_dict[ref_seq_id] = data_set_sample_seq.abundance
                new_count = ref_seq_to_count_dict[ref_seq_id]
            temp_ref_seq_id_list.append(ref_seq_id)
            unique_seq_name_base = '{}_{}'.format(ref_seq_id, proc_id)

            smp_name = str(cc).replace('-', '_')
            temp_fasta_dict['{}_{}'.format(unique_seq_name_base, 0)] = data_set_sample_seq.referenceSequenceOf.sequence
            temp_name_list = []
            for i in range(count, new_count):
                temp_name_list.append('{}_{}'.format(unique_seq_name_base, i))
                temp_group_list.append('{}\t{}'.format('{}_{}'.format(unique_seq_name_base, i), smp_name))

            temp_name_dict['{}_{}'.format(unique_seq_name_base, 0)] = temp_name_list
        output.put((temp_fasta_dict, temp_name_dict, temp_group_list, temp_ref_seq_id_list, proc_id, str(cc)))
    output.put('EXIT')

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

        print('Processing cc: {} with {}'.format(cc, proc_id))
        for data_set_sample_seq in data_set_sample_sequence.objects.filter(
                cladeCollectionTwoFoundIn=cc):
            ref_seq_id = data_set_sample_seq.referenceSequenceOf.id
            temp_ref_seq_id_list.append(ref_seq_id)
            unique_seq_name_base = '{}_id{}'.format(ref_seq_id, cc.id)

            smp_name = str(cc).replace('-', '_')
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
    # rename_tree(name_file, in_file_fconsense)
    # # ### DEBUG ###
    # for line in readDefinedFileToList(in_file_fconsense):
    #     if "1 Pro" in line:
    #         bob = 'asdf'
    # #############
    # out_file_fconsense = '{}/out_seq_boot_reps/consensus_tree_output_file'.format(clade_wkd)
    # tree_out_file_fconsense = '{}/out_seq_boot_reps/consensus_tree_consense.newick'.format(clade_wkd)
    tree_out_file_fconsense_sumtrees = '{}/consensus_tree_sumtrees.newick'.format(clade_wkd)
    # sumtree = local["sumtrees.py"]
    # output = (sumtree['-F', 'newick', in_file_fconsense] > tree_out_file_fconsense_sumtrees)()

    completed_consensus = subprocess.run(['sumtrees.py', '-F', 'newick', '--replace', '-o', tree_out_file_fconsense_sumtrees, in_file_fconsense], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if completed_consensus.returncode != 0:
        for line in completed_consensus.stdout.decode('utf-8'):
            print(line)
        for line in completed_consensus.sterr.decode('utf-8'):
            print(line)
        apples = 'asdf'

    # The consunsus tree output by sumtree is causing problems because it contains metadata.
    # It is important that the tree we use for the unifrac has the branch lengths on it
    # The tree output by consense only has the boostrap values instead of lengths which is not ideal
    # I am planning on removing the metadatafrom the sumtree consensus tree and seeing if that works in mothur's
    # Unifrac. This can be part of the rename_tree code
    rename_tree_two(name_file, tree_out_file_fconsense_sumtrees)
    return tree_out_file_fconsense_sumtrees

def generate_PCoA_coords(clade_wkd, raw_dist_file):
    # simultaneously grab the sample names in the order of the distance matrix and put the matrix into
    # a twoD list and then convert to a numpy array
    temp_two_D_list = []
    sample_names_from_dist_matrix = []
    for line in raw_dist_file[1:]:
        temp_elements = line.split('\t')
        sample_names_from_dist_matrix.append(temp_elements[0].replace(' ', ''))
        temp_two_D_list.append([float(a) for a in temp_elements[1:]])
    uni_frac_dist_array = np.array(temp_two_D_list)
    print('calculating PCoA coordinates')
    pcoA_full_path = clade_wkd + '/PCoA_coords.csv'
    this = pcoa(uni_frac_dist_array)

    # rename the dataframe index as the sample names
    mapper_dict = {i: j for i, j in enumerate(sample_names_from_dist_matrix)}
    this.samples['sample'] = sample_names_from_dist_matrix
    renamed_dataframe = this.samples.set_index('sample')
    print(renamed_dataframe)
    # now add the variance explained as a final row to the renamed_dataframe

    renamed_dataframe = renamed_dataframe.append(this.proportion_explained.rename('proportion_explained'))

    print('writing out coordinates')
    print(renamed_dataframe)
    renamed_dataframe.to_csv(pcoA_full_path, index=True, header=True, sep=',')
    return pcoA_full_path

def generate_fseqboot_alignments(clade_wkd, num_reps, out_file):
    # setup fseqboot arguments
    in_file_seqboot = out_file
    out_file_seqboot = in_file_seqboot + '.fseqboot'
    fseqboot = local["fseqboot"]
    # run fseqboot
    print('generating multiple datasets')
    (fseqboot['-sequence', in_file_seqboot, '-outfile', out_file_seqboot, '-test', 'b', '-reps', num_reps])()
    # Now divide the fseqboot file up into its 100 different alignments
    fseqboot_file = readDefinedFileToList(out_file_seqboot)
    rep_count = 0
    out_put_holder = []
    out_put_holder.append(fseqboot_file[0])
    fseqboot_base = '{}/out_seq_boot_reps/fseqboot_rep_'.format(clade_wkd)
    for line in fseqboot_file[1:]:
        if len(line.split()) == 2:
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
    print('calculating unifrac distances')
    completedProcess = \
        subprocess.run(['mothur', '{}'.format(mothur_batch_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if completedProcess.returncode == 0:
        print('Unifrac successful')
    else:
        print('ERROR: {}'.format(completedProcess.sterr.decode('utf-8')))
    return

# delete rename_tree
def rename_tree(name_file, tree_out_file_fconsense):
    name_file_reps = []
    for line in name_file:
        name_file_reps.append(line.split('\t')[0])
    seq_re = re.compile('\d+[_ ]Pr\w+-*')
    seq_re_meta = re.compile('\[[^\]]*\]')
    print('renaming tree nodes')
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
    print('writing out tree')
    writeListToDestination(tree_out_file_fconsense, new_tree_file)

def rename_tree_two(name_file, tree_out_file_fconsense):
    name_file_reps = []
    for line in name_file:
        name_file_reps.append(line.split('\t')[0])
    seq_re = re.compile('\d+[_ ]id[\d_]+')
    seq_re_meta = re.compile('\[[^\]]*\]')
    print('renaming tree nodes')
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
    print('writing out tree')
    writeListToDestination(tree_out_file_fconsense, new_tree_file)

def plotPCoA(full_path_to_coords):

    PCoA_coord_array = pd.read_csv(full_path_to_coords)

    # I'm going to assume that the sample order is the same in the PCoA output
    colours = []
    for smp in PCoA_coord_array['sample']:
        if smp[0] == 'Y':
            colours.append('blue')
        elif smp[0] == 'A':
            colours.append('black')
        elif smp[0] == 'C':
            colours.append('green')

    # convert back to numpy array for easy plotting
    # the last item in the series is the proportion explained
    PC1_plot_value = PCoA_coord_array['PC1'].tolist()[:-1]
    PC2_plot_value = PCoA_coord_array['PC2'].tolist()[:-1]

    fig, ax = plt.subplots()

    ax.set_xlabel('PC1 {0:.2f}%'.format(PCoA_coord_array['PC1'].iloc[-1]))
    ax.set_ylabel('PC2 {0:.2f}%'.format(PCoA_coord_array['PC2'].iloc[-1]))


    ax.scatter(PC1_plot_value, PC2_plot_value, color=colours)
    for i, label in enumerate(PCoA_coord_array['sample']):
        ax.annotate(label, (PC1_plot_value[i], PC2_plot_value[i]))

    fig.show()

def plotPCoA_types(full_path_to_coords):

    PCoA_coord_array = pd.read_csv(full_path_to_coords)

    # I'm going to assume that the sample order is the same in the PCoA output
    colours = []
    for smp in PCoA_coord_array['sample']:
        if 'C115c' in smp:
            colours.append('blue')
        elif 'C3aq' in smp:
            colours.append('black')
        elif 'C3am' in smp:
            colours.append('green')
        elif smp =='C3_C3c_C3gulf':
            colours.append('grey')

    # convert back to numpy array for easy plotting
    # the last item in the series is the proportion explained
    PC1_plot_value = PCoA_coord_array['PC1'].tolist()[:-1]
    PC2_plot_value = PCoA_coord_array['PC2'].tolist()[:-1]

    fig, ax = plt.subplots()

    ax.set_xlabel('PC1 {0:.2f}%'.format(PCoA_coord_array['PC1'].iloc[-1]))
    ax.set_ylabel('PC2 {0:.2f}%'.format(PCoA_coord_array['PC2'].iloc[-1]))


    ax.scatter(PC1_plot_value, PC2_plot_value, color=colours)
    for i, label in enumerate(PCoA_coord_array['sample']):
        ax.annotate(label, (PC1_plot_value[i], PC2_plot_value[i]))

    fig.show()

# mothur method
def mothur_unifrac_pipeline_MP_worker(input, output, fseqboot_base, clade_wkd):
    for p in iter(input.get, 'STOP'):
        print('Processing p={} with {}'.format(p, current_process().name))
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
        print('Calculating distances...')
        completedProcess = \
            subprocess.run(['mothur', '{}'.format(mothur_batch_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('Done.')
        # now run in clearcut
        input = '{}{}.sequential.square.dist'.format(fseqboot_base, p)
        mothur_batch_clearcut = \
            ['set.dir(input={}/out_seq_boot_reps/, output={}/out_seq_boot_reps/)'.format(clade_wkd, clade_wkd),
             'clearcut(phylip={}, verbose=t)'
                 .format(input)]

        mothur_batch_path = '{}/out_seq_boot_reps/mothur_batch_batch_clearcut_{}'.format(clade_wkd, p)
        writeListToDestination(mothur_batch_path, mothur_batch_clearcut)
        print('Generating NJ tree from distances')
        completedProcess = \
            subprocess.run(['mothur', '{}'.format(mothur_batch_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('Done')
        output.put(input.replace('.dist', '.tre'))
    output.put('kill')

def mothur_unifrac_pipeline_MP(clade_wkd, fseqboot_base, name_file, num_reps, num_proc):
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

    raw_dist_file = readDefinedFileToList(dist_file_path)
    return raw_dist_file

def mothur_unifrac_pipeline(clade_wkd, fseqboot_base, name_file, num_reps):
    ### MOTHUR TRIAL ###
    list_of_tree_paths = []
    for p in range(num_reps):
        # convert the interleaved fasta to sequential fasta
        interleaved_fast = readDefinedFileToList('{}{}'.format(fseqboot_base, p))
        sequen_fast = convert_interleaved_to_sequencial_fasta(interleaved_fast)
        writeListToDestination('{}{}.sequential.fasta'.format(fseqboot_base, p), sequen_fast)
        mothur_batch_dist = \
            ['set.dir(input={}/out_seq_boot_reps/, output={}/out_seq_boot_reps/)'.format(clade_wkd, clade_wkd),
             'dist.seqs(fasta={}, countends=T, output=square)'
                 .format('{}{}.sequential.fasta'.format(fseqboot_base, p))]
        mothur_batch_path = '{}/out_seq_boot_reps/mothur_batch_batch_dist'.format(clade_wkd)

        writeListToDestination(mothur_batch_path, mothur_batch_dist)

        # now run the batch file with mothur
        print('calculating distances')
        completedProcess = \
            subprocess.run(['mothur', '{}'.format(mothur_batch_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # now run in clearcut
        input = '{}{}.sequential.square.dist'.format(fseqboot_base, p)
        mothur_batch_clearcut = \
            ['set.dir(input={}/out_seq_boot_reps/, output={}/out_seq_boot_reps/)'.format(clade_wkd, clade_wkd),
             'clearcut(phylip={}, verbose=t)'
                 .format(input)]

        mothur_batch_path = '{}/out_seq_boot_reps/mothur_batch_batch_clearcut'.format(clade_wkd)
        writeListToDestination(mothur_batch_path, mothur_batch_clearcut)

        completedProcess = \
            subprocess.run(['mothur', '{}'.format(mothur_batch_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        list_of_tree_paths.append(input.replace('.dist', '.tre'))

    tree_out_file_fconsense_sumtrees = create_consesnus_tree(clade_wkd, list_of_tree_paths, name_file)

    perform_unifrac(clade_wkd, tree_out_file_fconsense_sumtrees)

    dist_file_path = '{}/{}'.format(clade_wkd,
                                    tree_out_file_fconsense_sumtrees.split('/')[-1]) + '1.weighted.phylip.dist'

    raw_dist_file = readDefinedFileToList(dist_file_path)
    return raw_dist_file


# phylip method
def phylip_unifrac_pipeline(clade_wkd, fseqboot_base, name_file, num_reps):
    ### PHYLIP METHOD ###
    fdnadist = local["fdnadist"]
    for p in range(num_reps):
        # run dnadist
        in_file_dnadist_rep = '{}{}'.format(fseqboot_base, p)
        out_file_dnadist_rep = '{}out_{}'.format(fseqboot_base, p)
        print('calculating distances rep {}'.format(p))
        (fdnadist['-sequence', in_file_dnadist_rep, '-outfile', out_file_dnadist_rep, '-method', 'j'])()

    fneighbor = local["fneighbor"]
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

    dist_file_path = '{}/{}'.format(clade_wkd, tree_out_file_fconsense_sumtrees.split('/')[-1]) + '1.weighted.phylip.dist'

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

    return raw_dist_file

def phylip_unifrac_pipeline_MP_worker(input, output, fseqboot_base):
    fdnadist = local["fdnadist"]
    fneighbor = local["fneighbor"]
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


# The fast_tree method of unifraing currently doesn't work for some reason. The output is non-sensical.
# It therefore also hasn't been multiprocessed
def run_fast_tree_unifrac_pipeline(clade_wkd, fseqboot_base, name_file, num_reps):
    #### FastTree Trial #####
    fasttree = local["fasttree"]
    list_of_tree_paths = []
    for p in range(num_reps):
        print('running fast tree rep {}'.format(p))
        interleaved_fast = readDefinedFileToList('{}{}'.format(fseqboot_base, p))
        sequen_fast = convert_interleaved_to_sequencial_fasta(interleaved_fast)
        out = '{}{}.sequential.fasta'.format(fseqboot_base, p)
        writeListToDestination(out, sequen_fast)
        (fasttree['-nt', out] > out.replace('.fasta', '.fasttree.tre'))()
        # list_of_trees_produced.append(out.replace('.fasta', '.fasttree.tre'))
        list_of_tree_paths.append(out.replace('.fasta', '.fasttree.tre'))


    tree_out_file_fconsense_sumtrees = create_consesnus_tree(clade_wkd, list_of_tree_paths, name_file)

    perform_unifrac(clade_wkd, tree_out_file_fconsense_sumtrees)

    dist_file_path = '{}/{}'.format(clade_wkd,
                                    tree_out_file_fconsense_sumtrees.split('/')[-1]) + '1.weighted.phylip.dist'

    raw_dist_file = readDefinedFileToList(dist_file_path)


    return raw_dist_file
#############################






###### Profile Discovery functions ######


def profileDiscovery(nProcessors):


    if not analysisObj.initialTypeDiscoComplete:
    ############# FIND RAW FOOTPRINTS ###############
        cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

        # Get the cladeCollections that are found in the listof datasubmissions that are found in the analysisObj
        cladeCollectionsOfAnalysis = analysisObj.getCladeCollections()


        # List that will hold a dictionary for each clade
        # Each dictionary will hold key = footprint (set of sequences) value = [[] []] where []0 = list of cladeCollections containing given footprint
        # and []1 = list of majority sequence for given sample
        masterCladalListOfFootprintDicts = [{} for clade in cladeList]


        # For each clade Collection add its footprint to the count dict
        # by associating the majority sequence and sample to the corect cladal dict

        # This queue will have the cladeCollections to be processed
        taskQueue = Queue()
        # This will hold the outputted results from the CC being processed
        outputQueue = Queue()
        # This queue will be used to fetch the output dict list from the second listenerworker
        dictResultQueue = Queue()

        for cladecollection in cladeCollectionsOfAnalysis:
            taskQueue.put(cladecollection)


        numProcessors = nProcessors
        for N in range(numProcessors):
            taskQueue.put('STOP')

        allProcesses = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        # Then start the workers
        # workerDiscoveryTwoWorker will process the CCs and pass this info on to the second queue which
        # workerDiscoveryTwoListener will work on
        # Finally workerDiscoveryTwoListener will output its results to the third queue
        for N in range(numProcessors):
            p = Process(target=workerDiscoveryTwoWorker, args=(taskQueue, outputQueue, analysisObj.withinCladeCutOff))
            allProcesses.append(p)
            p.start()

        # Process the output of the multiprocessing
        # We were doing this in a separate worker but it turnsout that queues have a fairly small size limit
        # so we weren't able to pass out the mastercladallistoffootprintdicts after it had been populated
        # As we were only able to do this with one process only it is no loss for us to process this directly in the
        # main process.
        ## http://stackoverflow.com/questions/21641887/python-multiprocessing-process-hangs-on-join-for-large-queue (bob)
        killNum = 0
        while 1:
            passedElement = outputQueue.get()
            if passedElement == 'kill':
                killNum += 1
                if killNum == numProcessors:
                    break
            else:
                # footprintInQ = dictelement[0]
                # cladalDictionaryKey = dictelement[1]
                # CC = dictelement[2]
                # CC.maj() = dictelement[3]
                # TODO 07/12/17 here we are going to start to make changes to the format of the footprintlist
                # called mastercldallistoffootprintdicts
                # For the maj types (passedElement[3]) we are going to put them into their own list rather than
                # have them as items in a single list
                # e.g. a 3d list instead of a 2D list

                # if passedElement[0] not in masterCladalListOfFootprintDicts[passedElement[1]]:
                #     masterCladalListOfFootprintDicts[passedElement[1]][passedElement[0]] = [[passedElement[2]],
                #                                                                             [passedElement[3]]]
                # else:
                #     masterCladalListOfFootprintDicts[passedElement[1]][passedElement[0]][0].append(passedElement[2])
                #     masterCladalListOfFootprintDicts[passedElement[1]][passedElement[0]][1].append(passedElement[3])

                if passedElement[0] not in masterCladalListOfFootprintDicts[passedElement[1]]:
                    masterCladalListOfFootprintDicts[passedElement[1]][passedElement[0]] = [[passedElement[2]],
                                                                                            [[passedElement[3]]]]
                else:
                    masterCladalListOfFootprintDicts[passedElement[1]][passedElement[0]][0].append(passedElement[2])
                    masterCladalListOfFootprintDicts[passedElement[1]][passedElement[0]][1].append([passedElement[3]])

        # First wait for the workers to finish
        for p in allProcesses:
            p.join()
        ################################################


        ########### CHECK RAW FOOTPRINTS FOR SUPPORT AND GENERATE SUPPORTED TYPE PROFILES ###############
        #Now work clade by clade
        for footPrintDict in masterCladalListOfFootprintDicts:
            if footPrintDict: # If there are some clade collections for the given clade
                # The fact that there are this few cladecollections of a clade
                # will be very rare, and in this rare case the Majs will simply be associated to the footprints
                ccCount = sum(len(footPrintDict[akey][0]) for akey in footPrintDict.keys())

                # ##### Debugging (Delete Me) ######
                # for keyItem in footPrintDict.keys():
                #     print(','.join([str(item) for item in keyItem]))
                #     print(','.join([str(cc) for cc in footPrintDict[keyItem][0]]))
                # ##################################
                ######## FIND WHICH FOOTPRINTS ARE SUPPORTED AND WHICH CCs SUPPORT THEM #######
                # collapsedFootPrintDict = collapsePotentialProfiles(footprintList = footPrintDict, reqsupport=max(4, math.ceil(analysisObj.typeSupport*ccCount)), nprocessors=nProcessors)
                # ### DEBUG ###
                # if masterCladalListOfFootprintDicts.index(footPrintDict) == 2:
                #     foo = 'bar'
                # ### DEBUG ###
                collapsedFootPrintDict = collapsePotentialProfiles_initType_objects(footprintList=footPrintDict, reqsupport=4, nprocessors=nProcessors)

                # collapsedFootPrintDict = collapsePotentialProfiles_initType_objects(footprintList=masterCladalListOfFootprintDicts[2],
                #                                                                     reqsupport=4,
                #                                                                     nprocessors=nProcessors)
                ###############################################################################

                ##### Debugging (Delete Me) ######
                if masterCladalListOfFootprintDicts.index(footPrintDict) == 2:
                    # this is Cs
                    apple = 'asdf'


                ## DEBUG
                for initT in collapsedFootPrintDict:
                    print(str(initT.profile))
                    print(','.join([str(cc) for cc in initT.cladeCollection_list]))
                    for cc in initT.cladeCollection_list:
                        if str(cc) in ['July2016-113poc-R1','July2016-59poc-R1','March2016-125poc-R1','March2016-57poc-R1','May2016-113poc-R1','May2016-135poc-R1','May2016-294por-R1','March2016-294poc-R1','May2016-59poc-R1','May2016-125poc-R1']:
                            apples = 'pwers'
                            # thi is a cct that ends up as the C1
                        if str(cc) in ['July2016-159poc-R1','May2016-159poc-R1','March2016-55poc-R1','May2016-55poc-R1','March2016-135poc-R1','July2016-55poc-R1','July2016-135poc-R1','May2016-235poc-R1']:
                            apples = 'asdf'
                            # then this is the C1k types









                ## DEBUG
                ##################################

        ############################################################################################




        ############ CREATE ANALYSIS TYPES BASED ON DISCOVRED FOOTPRINTS ############################
                # TODO 08/12/17 we need to be careful here when we initiate the types as the types we were previously
                ''' generating would represent essentially the majoirty of the ccts sequences. but now some of the types
                will be smaller proportions of the ccts so we should check to see how the initial abundnace of the types
                are calculated. e.g. are they calculated as the proportion of the total seqs in the cct or are we 
                already working on as proportions of the seqs of the type in question. Hopefully it is the latter
                and we were just going with the types that represented the largest number of sequences for the cct.
                '''
                # TODO 08/12/17 I have carefully looked through the type initTypeAtributes method
                ''' Firsly it always works in the context of the sequences found in the type. It produces absoulte
                 counts per sequence in the type for each cladeCollection that the type was supported by.
                 It also produces a count that is relative proportions of each sequence of the type for each CC.
                  Hopefully this is what we are using when we do the second artefact check. I.e. we are looking
                  for the types in the CCs again.
                  For each type, the abslute counts per type sequence per CC are stored in type.footprintSeqAbundances
                  in the order of initalCCs and orderedfootprint list. THere is also the relative version wihich is 
                  stored as type.footprintSeqRatios'''
                # for every footprint that will become an analysis_type
                timeOne = 0
                timeTwo = 0
                print('\n\nCreating analysis types clade {}'.format(cladeList[masterCladalListOfFootprintDicts.index(footPrintDict)]))
                for initialType in collapsedFootPrintDict:
                    footPrintNameString = ','.join([str(refseq) for refseq in initialType.profile])

                    # Work out the corresponding reference_sequence for each Maj of the samples with that corresponding type
                    # Then do len(set()) and see if it is a coDom, i.e. different Maj seqs within the type

                    timeitone = timeit.default_timer()
                    listOfSampSeqs = []



                    # setOfMajSeqs = set(listOfSampSeqs)
                    timeOne += timeit.default_timer() - timeitone
                    timeittwo = timeit.default_timer()

                    if len(initialType.set_of_maj_ref_seqs) > 1: # Then this is a coDom
                        coDom = True

                        # the Counter class (from collections import Counter) may be useful
                        # http://stackoverflow.com/questions/2600191/how-can-i-count-the-occurrences-of-a-list-item-in-python
                        newAnalysisType = analysis_type(coDom=True, dataAnalysisFrom=analysisObj, clade=cladeList[masterCladalListOfFootprintDicts.index(footPrintDict)])

                        newAnalysisType.setMajRefSeqSet(initialType.set_of_maj_ref_seqs)
                        newAnalysisType.initTypeAttributes(initialType.cladeCollection_list, initialType.profile)

                        newAnalysisType.save()
                        print('\rCreating analysis type: {}'.format(newAnalysisType.name), end='')
                    else:
                        coDom = False
                        newAnalysisType = analysis_type(coDom=False, dataAnalysisFrom=analysisObj,
                                           clade=cladeList[masterCladalListOfFootprintDicts.index(footPrintDict)])
                        newAnalysisType.setMajRefSeqSet(initialType.set_of_maj_ref_seqs)
                        newAnalysisType.initTypeAttributes(initialType.cladeCollection_list, initialType.profile)
                        newAnalysisType.save()
                        print('\rCreating analysis type: {}'.format(newAnalysisType.name), end='')

                        apples = 'asdf'
                    timeTwo += timeit.default_timer() - timeitone
                print('\nTimeOne = {}'.format(timeOne))
                print('TimeTne = {}'.format(timeTwo))

                # # DEBUG#
                # for at in analysis_type.objects.filter(dataAnalysisFrom=analysisObj):
                #     print(str(at))
                #     print(','.join([str(cc) for cc in at.getCladeCollectionsFoundInInitially()]))
                #     if 'January2017-280por' in [str(cc) for cc in at.getCladeCollectionsFoundInInitially()]:
                #         appeles = 'pwers'
                # ###
        #############################################################################################

        ############ CHECK FOR ADDITIONAL TYPES NOT FOUND DUE TO INTRAS SPANNING THE WITHINCLADECUTOFF BOUNDARY #########



        print('\n\nChecking for additional artefact types')
        ### CREATION OF SUPER TYPES DUE TO WITHINCLADECOLLECTIONCUTOFF ARTEFACTS ###
        CCToTotalSeqsDict, CCToRefSeqListAndAbundances, typeFootPrintDict, CCToInitialTypeDict = checkForAdditionalArtefactTypes(nProcessors)
        ############################################################################
        testDirPath = os.path.join(os.path.dirname(__file__), 'temp/{}'.format(analysisObj.id))
        os.makedirs(testDirPath, exist_ok=True)
        os.chdir(testDirPath)

        # I think that we must convert the Manager().dict objets to normal dicts before we pickle them.
        # I think this is the reason we are unable to pickle.load the files.
        CCToTotalSeqsDictToDump = dict(CCToTotalSeqsDict)
        pickle.dump(CCToTotalSeqsDictToDump, open("CCToTotalSeqsDict_{}".format(analysisObj.id), "wb"))

        CCToRefSeqListAndAbundancesToDump = dict(CCToRefSeqListAndAbundances)
        pickle.dump(CCToRefSeqListAndAbundancesToDump, open("CCToRefSeqListAndAbundances_{}".format(analysisObj.id), "wb"))

        typeFootPrintDictToDump = dict(typeFootPrintDict)
        pickle.dump(typeFootPrintDictToDump, open("typeFootPrintDict_{}".format(analysisObj.id), "wb"))

        CCToInitialTypeDictToDump = dict(CCToInitialTypeDict)
        pickle.dump(CCToInitialTypeDictToDump, open("CCToInitialTypeDict_{}".format(analysisObj.id), "wb"))
        # also we can perhaps have a save point like analysisTypesDefined below that is associated to the analysisObj
        analysisObj.initialTypeDiscoComplete = True
        analysisObj.save()

        reassessSupportOfArtefactDIVContainingTypes(CCToTotalSeqsDict, CCToRefSeqListAndAbundances, typeFootPrintDict, CCToInitialTypeDict, nProcessors)

    else:
        # print('Checking for additional artefact types')
        # CCToTotalSeqsDict, CCToRefSeqListAndAbundances, typeFootPrintDict, CCToInitialTypeDict = checkForAdditionalArtefactTypes()
        testDirPath = os.path.join(os.path.dirname(__file__), 'temp/{}'.format(analysisObj.id))
        os.chdir(testDirPath)

        CCToTotalSeqsDict = Manager().dict(pickle.load(open("{}/CCToTotalSeqsDict_{}".format(testDirPath, analysisObj.id), "rb")))
        CCToRefSeqListAndAbundances = Manager().dict(pickle.load(open("CCToRefSeqListAndAbundances_{}".format(analysisObj.id), "rb")))
        typeFootPrintDict = Manager().dict(pickle.load(open("typeFootPrintDict_{}".format(analysisObj.id), "rb")))
        CCToInitialTypeDict = Manager().dict(pickle.load(open("CCToInitialTypeDict_{}".format(analysisObj.id), "rb")))

        reassessSupportOfArtefactDIVContainingTypes(CCToTotalSeqsDict, CCToRefSeqListAndAbundances, typeFootPrintDict, CCToInitialTypeDict, nProcessors)

    # # Now print a list of the types
    # allTypesFromDataAnalysis = analysis_type.objects.filter(dataAnalysisFrom=analysisObj)
    # for an in allTypesFromDataAnalysis:
    #     print(an.name)
    #     if an.name == 'D1/D4/D17-D6-D2.2-D2':
    #         apples = 'asdf'
    # apples = 'asdf'
    #################################################################################################################
    analysisObj.analysisTypesDefined = True
    analysisObj.save()
    return

def reassessSupportOfArtefactDIVContainingTypes(CCToTotalSeqsDict, CCToRefSeqListAndAbundances, typeFootPrintDict, CCToInitialTypeDict, cores):
    # TODO 08/12/17 13:41 this is where we're at. We have fixed the cc to initial type dict for the other artefact
    # checking but still need to work on that here as well as the other issues noted below.

    '''
    TODO 08/12/17 This is going to cause some problems with the basal comparisons.
    firstly it is assuming that each cct can only associate with  only one type.
    Instead we will just have to make sure to do checks that mean that each cct can only support one type of each
    basal group.
    We will also have to write checks in at the stages where there is a potential for new types to be created.
    It is important that we are checking types according to relative abundances and ratios. We know that this
    information is stored in the types. For more information on this see the comment that is made whwere
    we first create the inital types.

    Although we are unlocking the DIVs that are found at low abundance in the type Assignement this still does
    not enable a lot of CCs to be assigned to the unlocked types even if they meet the requirements of the unlocked
    DIV. This is because they fall outside the acceptable range of one or more of the other DIVs. This is due to
    the fact that these other ranges are still being defined by the inital CCs that are associated to these types.
    To fix this we must reassess the initally associated CCs to every one of the types that contains an artefact DIV.

    For each type currently found, if it contains an artefact div, go through all CCs to see if they fit the normal
    requeirements e.g. >0.03 for each of the DIVs, BUT allow the 0.005 cutoff for the unlocked DIVs. If you find
    a CC that does fit the type, and it is not that CCs current inital type, then assess to see whether the new
    type represents more seqs than the CCs current inital type. If it does, then add it to a list of CCs that will
    need to be added to this type. There is code from the previous 'checkForAdditionalAretfactTypes' that will
    help with this. E.g. collecting all of the CCs first, then goinging type by type for the CCs that need to be
    removed and checking to see if the types are still viable. If they are not, then we need to delete the type
    and rehome the CCs that are hanging. This was also done in the preivous section. At least no new types will
    be being created.

    By this means, more CCs will be able to be associated to the types that contain these artefact DIVs.
    '''

    # Get list of all types
    allTypesFromDataAnalysis = analysis_type.objects.filter(dataAnalysisFrom=analysisObj)

    # Get list of clades that are represented by the types
    cladeList = set()
    for at in allTypesFromDataAnalysis:
        cladeList.add(at.clade)

    for currentClade in cladeList:
        checked = []
        while 1:
            restart = False
            # Get list of type from scratch as we may be restarting this
            allTypesFromDataAnalysis = analysis_type.objects.filter(dataAnalysisFrom=analysisObj)
            cladalTypesIDs = [at.id for at in allTypesFromDataAnalysis if at.clade == currentClade]
            # For each type
            for typeToCheckID in cladalTypesIDs:
                # if this type has not already been checked
                if typeToCheckID not in checked:
                    print('\n\nChecking {}'.format(analysis_type.objects.get(id = typeToCheckID)))
                    # If this type contains artefact DIVs
                    if typeFootPrintDict[typeToCheckID][2]:
                        # Then this is a type we need to check

                        # This is the list that will hold the CCs that are going to be added to the type's
                        # inital CCs list.
                        supportList = []

                        # This is the info we have to work with from the previous dictionaries
                        artefactDIVIDs = typeFootPrintDict[typeToCheckID][2]
                        nonArtefactDIVIDs = typeFootPrintDict[typeToCheckID][0]
                        refSeqObjsOfTypeList = typeFootPrintDict[typeToCheckID][3]

                        # For each type generate a requirements dict
                        # Dict will be DIV:RequiredRelAbund
                        # Otherwise known as pnt
                        requirementDict = {refSeqObj: (0.03 if refSeqObj.id in nonArtefactDIVIDs else unlockedAbund) for refSeqObj in refSeqObjsOfTypeList}

                        listOfCCsToCheck = [cc for cc in analysisObj.getCladeCollections() if cc.clade == currentClade]

                        taskQueue = Queue()
                        supportListMan = Manager()
                        supportList = supportListMan.list()
                        # outputQueue = Queue()

                        for CC in listOfCCsToCheck:
                            taskQueue.put(CC)

                        for N in range(cores):
                            taskQueue.put('STOP')

                        allProcesses = []

                        # In this case the db.connections were not being recreated and I think this might have something
                        # to do with a connection only being created when there is writing to be done.
                        # as there was only reading operations to be done in the worker task I think that maybe no new connections
                        # were being made.
                        db.connections.close_all()
                        # I think this might be being caused due to the changing of the address of the db in the settings.py file to a relative path
                        # yep this is what was causing the problem so we will need to find some way of defining the location of the database relatively.
                        # DONE. I have changed the settings.py file so that the DB location is now relative to the position of the settings.py file


                        for N in range(cores):
                            p = Process(target=workerArtefactTwo, args=(taskQueue, supportList, CCToInitialTypeDict, typeToCheckID, CCToRefSeqListAndAbundances, refSeqObjsOfTypeList, requirementDict  ))
                            allProcesses.append(p)
                            p.start()

                        for p in allProcesses:
                            p.join()



                        if supportList:
                            restart = True
                            # Then we have found CCs that should be moved to the Type in Q
                            ################# REDISTRIBUTE SUPPORTING TYPES ###############
                            # Now remove the CCs from the types they were previously associated to and update the
                            # CCToInitialTypeDict accordingly. # We will also likely need to reinitiate each of the types
                            # There are likely to be types that have more than one CC being removed from them so the most effective
                            # way to process this is not to go CC by CC but rather type by type
                            # So first get a list of the types that have been effected
                            setOfTypesAffected = set()
                            # at the same time create a dict that tracks which CCs need to be remvoed from each type
                            typeToCCToBeRemovedDict = defaultdict(list)  # we may have to do this manually instead of defaultdict
                            # for CC in supportList:
                            #     if CC.id in CCToInitialTypeDict.keys():
                            #         initialTypeID = CCToInitialTypeDict[CC.id]
                            #         setOfTypesAffected.add(initialTypeID)
                            #         typeToCCToBeRemovedDict[initialTypeID].append(CC)
                            #     else:
                            #         # we don't need to remove this CC from any type but we should add it to the CCToInitialTypeDict
                            #         CCToInitialTypeDict[CC.id] = typeToCheckID


                            for CC in supportList:
                                if CC.id in CCToInitialTypeDict.keys():
                                    initialTypeID = None
                                    if len(CCToInitialTypeDict[CC.id]) == 1:
                                        initialTypeID = CCToInitialTypeDict[CC.id][0]
                                        setOfTypesAffected.add(initialTypeID)
                                        typeToCCToBeRemovedDict[initialTypeID].append(CC)
                                    else:
                                        initialTypeID = find_which_type_is_same_basal_type_as_pnt(
                                            requirementDict, analysis_type.objects.filter(id__in=CCToInitialTypeDict[CC.id])).id
                                        ### DEBUG ###
                                        if not initialTypeID:
                                            foo = 'bar'
                                        ### DEBUG ###
                                        setOfTypesAffected.add(initialTypeID)
                                        typeToCCToBeRemovedDict[initialTypeID].append(CC)
                                else:# we don't need to remove this CC from any type but we should add it to the CCToInitialTypeDict
                                    CCToInitialTypeDict[CC.id] = [typeToCheckID]





                            # Here we have a set of the type IDs for the types affected and a dict associating the CCs to these types
                            # Now go through the types and remove the CCs from the type, change the dict association
                            # and reinitialise the type

                            # THe list to catch and reconsider all of the stranded CCs that might be left over in this process

                            ############## TYPE BY TYPE OF REDISTRIBUTED CCs DELETE OR REINITIATE #################
                            print('Reassessing support of affected types')
                            strandedCCs = []
                            for anType in typeToCCToBeRemovedDict.keys():
                                anTypeInQ = analysis_type.objects.get(id=anType)

                                #################### REMOVE CCs FROM TYPES ######################
                                # Remove CCs from type
                                listOfCCsToBeRemovedStrID = [str(cc.id) for cc in typeToCCToBeRemovedDict[anType]]
                                anTypeInQ.removeCCListFromInitialCladeCollectionList(listOfCCsToBeRemovedStrID)

                                # Change the CC associations in the cctocurrent... dict
                                for ccstrid in listOfCCsToBeRemovedStrID:
                                    CCToInitialTypeDict[int(ccstrid)].remove(anType)
                                    CCToInitialTypeDict[int(ccstrid)].append(typeToCheckID)

                                ################### REASSESS SUPPORT OF TYPE ###################
                                # Now check to see if the typeinQ still has sufficient support
                                # If it does then reinitiate it
                                # BUT if type is only 1 intra in length, then don't require 4 for support
                                # As some of these will not have 4 to start with
                                ################ IF SUPPORTED ##################


                                listOfCCsInType = [cc for cc in clade_collection.objects.filter(id__in=[int(x) for x in anTypeInQ.listOfCladeCollectionsFoundInInitially.split(',') if x != ''])]
                                if listOfCCsInType and len(anTypeInQ.getOrderedFootprintList()) == 1:
                                    print('Short Type {} supported by {} CCs. Reinitiating.'.format(anTypeInQ.name,
                                                                                              len(listOfCCsInType)))
                                    # Then this still has suffiecient support
                                    # reinitiate the type
                                    # listOfCCs = [cc for cc in clade_collection.objects.filter(id__in=[int(x) for x in anTypeInQ.listOfCladeCollectionsFoundInInitially.split(',')])]
                                    anTypeInQ.initTypeAttributes(listOfCC=listOfCCsInType,
                                                                 footprintlistofrefseqs=anTypeInQ.getOrderedFootprintList())
                                elif len(listOfCCsInType) >= 4:
                                    print('Type {} supported by {} CCs. Reinitiating.'.format(anTypeInQ.name, len(listOfCCsInType)))
                                    # Then this still has suffiecient support
                                    # reinitiate the type
                                    # listOfCCs = [cc for cc in clade_collection.objects.filter(id__in=[int(x) for x in anTypeInQ.listOfCladeCollectionsFoundInInitially.split(',')])]
                                    anTypeInQ.initTypeAttributes(listOfCC=listOfCCsInType, footprintlistofrefseqs=anTypeInQ.getOrderedFootprintList())
                                ################ IF UNSUPPORTED #################
                                else:
                                    # THen this antype no longer has sufficient support
                                    # put the CCs into the stranded list and delete this type
                                    # also remove the CCs from the dict of types they are associated to
                                    print('Type {} no longer supported. Deleting. {} CCs stranded.'.format(anTypeInQ.name, str(len(listOfCCsInType))))
                                    del typeFootPrintDict[anTypeInQ.id]
                                    anTypeInQ.delete()

                                    for cc in listOfCCsInType:
                                        if len(CCToInitialTypeDict[cc.id]) > 1:
                                            # Then there are other types associated to this cladeCollection and we should simply
                                            # remove the type in question from the list
                                            CCToInitialTypeDict[cc.id].remove(anType)
                                        else:
                                            # then this only contains one type and we should delte the cc entry in the dict and it will
                                            # become stranded.
                                            del CCToInitialTypeDict[cc.id]
                                            # strandedCCs.extend(listOfCCsInType)
                                            strandedCCs.append(cc)



                            #################### ATTEMPT TO HOME STRANDED TYPES ######################
                            ################### GET COMMON INTRAS #################
                            # Only attempt to make a new type to associate the stranded CCs to if there is sufficient support
                            if len(strandedCCs) >= 4:
                                totalIntraSet = set()
                                for CC in strandedCCs:
                                    totalIntraSet.update(CC.cutOffFootprint(analysisObj.withinCladeCutOff))

                                # Now go through each CC again and remove any intra from the total list that isn't found in the
                                # CCinQ to produce an in common list
                                refSeqsToRemove = set()
                                for CC in strandedCCs:
                                    intrasInCCInQ = CC.cutOffFootprint(analysisObj.withinCladeCutOff)
                                    for refSeq in list(totalIntraSet):
                                        if refSeq not in intrasInCCInQ:
                                            refSeqsToRemove.add(refSeq)

                                # Now create the commons list
                                intrasInCommonList = [refSeq for refSeq in list(totalIntraSet) if
                                                      refSeq not in refSeqsToRemove]

                                exists = False
                                if intrasInCommonList:
                                ############ CHECK IF POTENTIAL NEW TYPE FOOTPRINT ALREADY EXISTS ############
                                # Check to see if the potential newtypes footprint already exists
                                    pntFootprint = set([refSeq.id for refSeq in intrasInCommonList])

                                    typeThatExistsID = 0
                                    for key, footprintdictvalues in typeFootPrintDict.items():
                                        if footprintdictvalues[1] == pntFootprint:
                                            exists = True
                                            typeThatExistsID = key
                                            break

                                ############ IF ALREADY EXISTS; ASSOCIATE STRANDED CCs TO EXISTING TYPE #################
                                if exists:

                                    # If this type already exists then the CCsInQ will have a good type
                                    # to be assigned to and we need do no more
                                    # if this footprint already exists then we should associate these CCs to the type that has this footprint
                                    # so do this in the dictionary but also add these CCs to the cladeCollectionFoundInInitially list and
                                    # reinitiate the type
                                    associateCCsToExistingTypeAndUpdateDicts(cctocurrentinitialtypedict = CCToInitialTypeDict, strandedCCs = strandedCCs, typeThatExistsID=typeThatExistsID, typefootprintdict=typeFootPrintDict)


                                ############ IF DOESN'T EXIST; MAKE NEW TYPE AND ASSOCIATED CCs ###########
                                elif not exists and intrasInCommonList:
                                    # TODO 08/12/17 here we need to make sure that the intrasIncommonList doesn't contain a mixture of
                                    # basal seqs
                                    # will return True if there are multiple basal type in the intras in common list
                                    if not check_if_intrasInCommonList_contains_multiple_basal_seqs(intrasInCommonList):
                                        # He we should create a new type based on the CC collection and footprint above
                                        # and house the CCs in it
                                        makeNewTypeAndAssociateCCsAndUpdateDicts(cctocurrentinitialtypedict = CCToInitialTypeDict, clade = currentClade, intrasInCommonList = intrasInCommonList, strandedCCs = strandedCCs, typefootprintdict = typeFootPrintDict)
                                ############# FIND A NEW HOME IN EXISTING TYPES IF INSUFFICIENT SUPPORT FOR NEW TYPE ##################
                            else:
                                reassociateCCsToExistingTypesAndUpdateDicts(cctocurrentinitialtypedict = CCToInitialTypeDict, cctorefabunddict = CCToRefSeqListAndAbundances, clade = currentClade, strandedCCs = strandedCCs, typefootprintdict = typeFootPrintDict)


                            ############# FINALLY, ADD CCs TO TYPEINQ AND REINITIATE ##############
                            # N.B. The CCs have already had their CCToInitialTypeDict adjusted
                            typeToCheckObj = analysis_type.objects.get(id=typeToCheckID)
                            currentListOfCCsInType = [cc for cc in clade_collection.objects.filter(id__in=[int(x) for x in typeToCheckObj.listOfCladeCollectionsFoundInInitially.split(',') if x != ''])]
                            updatedListOfCCsInType = []
                            updatedListOfCCsInType.extend(currentListOfCCsInType)
                            updatedListOfCCsInType.extend(supportList)

                            typeToCheckObj.initTypeAttributes(listOfCC=updatedListOfCCsInType, footprintlistofrefseqs=typeToCheckObj.getOrderedFootprintList())
                            print('Added {} CCs to {}'.format(len(supportList), typeToCheckObj))
                        checked.append(typeToCheckID)
                        if restart:
                            break
                    else:
                        checked.append(typeToCheckID)
            # If we make it here then we should have been all the way through without having to restart
            # All that remains is to break the While
            if not restart:
                break
        # ### DEBUG ###
        # for at in analysis_type.objects.filter(dataAnalysisFrom=analysisObj, clade=currentClade):
        #     print('{}'.format(at))
        #     print('{}'.format(','.join([str(cct) for cct in at.getCladeCollectionsFoundInInitially()])))
        #     if 'January2017-280por' in [str(cct) for cct in at.getCladeCollectionsFoundInInitially()]:
        #         appeles = 'pwers'
        # #############
    return

def check_if_intrasInCommonList_contains_multiple_basal_seqs(intras_in_common_list):
    basalCount = 0
    C15Found = False
    for rs in intras_in_common_list:
        if 'C3' == rs.name:
            basalCount += 1
        elif 'C1' == rs.name:
            basalCount += 1
        elif 'C15' in rs.name and C15Found == False:
            basalCount += 1
            C15Found = True
    if basalCount > 1:
        return True
    return False

def check_whether_pnt_needs_comparing_to_current_type(requirementDict ,list_of_analysis_types):
    '''This should return the typeobject that has the same basal sequence as the pnt
            either C3, C15, C1 or None'''
    # first get the basal of the pnt
    basalpnt = False
    for rs, abund in requirementDict.items():
        basalpnt = False
        if 'C3' == rs.name:
            basalpnt = 'C3'
            break
        elif 'C1' == rs.name:
            basalpnt = 'C1'
            break
        elif 'C15' in rs.name:
            basalpnt = 'C15'
            break




    # for each at check to see if the result is informative
    for at in list_of_analysis_types:
        basal = False
        for rs in at.getOrderedFootprintList():
            if 'C3' == rs.name:
                basal = 'C3'
                break
            elif 'C1' == rs.name:
                basal = 'C1'
                break
            elif 'C15' in rs.name:
                basal = 'C15'
                break
        if basalpnt == basal:
            # then both types contain the same basal and so this is our type
            return at

    # we will not let CC support types e.g. pnt's that arent the same basal as a type that is already found in the CC's types
    return False

def find_which_type_is_same_basal_type_as_pnt(requirementDict ,list_of_analysis_types):
    '''This should return the typeobject that has the same basal sequence as the pnt
        either C3, C15, C1 or None'''
    # first get the basal of the pnt
    basalpnt = False
    for rs, abund in requirementDict.items():
        basalpnt = False
        if 'C3' == rs.name:
            basalpnt = 'C3'
            break
        elif 'C1' == rs.name:
            basalpnt = 'C1'
            break
        elif 'C15' in rs.name:
            basalpnt = 'C15'
            break

    # here we know if pnt contains a basal seq

    # for each at check to see if the result is informative
    for at in list_of_analysis_types:
        basal = False
        for rs in at.getOrderedFootprintList():
            if 'C3' == rs.name:
                basal = 'C3'
                break
            elif 'C1' == rs.name:
                basal = 'C1'
                break
            elif 'C15' in rs.name:
                basal = 'C15'
                break
        if basalpnt == basal:
            # then both types contain the same basal and so this is our type
            return at
    return False

def workerArtefactTwo(input, supportList, CCToInitialTypeDict, typeToCheckID, CCToRefSeqListAndAbundances, refSeqObjsOfTypeList, requirementDict):
    # Now that we have created all of the type this analysis is going to have, i.e. by doing the artefact checks
    # we now need to check to see if any of the CCs want to support these new types
    # hence here we go through each of the types and each of the CCs in turn for each of these types
    # If we find support then we change the current type of the CC
    for CC in iter(input.get, 'STOP'):
        print('\r{}'.format(CC), end='')
        ######## MP HERE ##########
        # We will need to have a managed list here.
        # We will be able to pass in the managed dicts

        # if str(CC) in ['January2016-148poc', 'January2017-148poc', 'January2017-81poc']:
            # if str(analysis_type.objects.get(id=typeToCheckID)) == 'C42g-C42a-C42.2-C42h-C1-C42b':
            #     apples = 'pears'

        # TODO DEBUG it seems that some of the CCs' ID are not being found in the CCToInitialTypeDict
        # TODO The dict is created earlier on  line 1609. It is created by going through every type
        # in the collecting all of the CCs that it was found in initially and then creating the dict by
        # key = CC id and value = the type.
        # So it is quite possible that if we are going through each CC in the analysis and if a CC didn't have
        # an initial type associated to it then we will have a key error here as the CC would not be incorporated into
        # the CCToInitialTypeDict
        #TODO check to see if we are going through all CCs in this MP
        # It turns out that we are going through the 'list of CCs to check'
        # @ line 609: listOfCCsToCheck = [cc for cc in analysisObj.getCladeCollections() if cc.clade == currentClade]
        # So as you can see it is entirely possible that some of the CCs did not have initial types assigned to them
        # So the question is, why did I make such a rookie error and why did this not break earlier
        # Which leads me to question should each of the CCs be given an initial type, and if so, why are some of the
        # samples not being given initial types.

        # If the CCs current type is the type in question, go to next CC.

        # I am putting in a conditional here to check that the CC.id is in the dict first
        # I see no harm in this. If the CC doesn't already have a type associated to it then we should check
        # to see whether the the current type could fit the cc
        if CC.id in CCToInitialTypeDict.keys():
            if typeToCheckID in CCToInitialTypeDict[CC.id]:
                # then the type we are checking is already associated with the CC in question
                continue
        # else if the CC doesn't currently have a type associated to it then there is no problem continuting to see
        # if it could support the current type in Q

        # Check to see that each of the types DIVs are found in the CC
        CCRelAbundDict = CCToRefSeqListAndAbundances[CC.id]
        refSeqsInCC = CCRelAbundDict.keys()
        #We don't need to do a maj seq check in here because of the check_whether_pnt_needs_comparing_to_current_type
        # function below that will only allow a CC to support a type if it already contains a type of the same basal
        # sequence. This way it is very unlikely that a new pnt will be able to represent more seqs than the current
        # type unless it has the maj seq in it. And even if it does, fair enough, it covers more seqs.
        if set(refSeqObjsOfTypeList).issubset(set(refSeqsInCC)):
            # Then the CC in question contains the intras of the type in question
            # Now check to see if the rel abund requirements are met
            notMet = False
            for refSeqKey in requirementDict.keys():
                if CCRelAbundDict[refSeqKey] < requirementDict[refSeqKey]:
                    # Then this CC does not have the required DIVs at the required rel abund for typeInQ
                    # Move on to check the next CC
                    notMet = True
                    break
            if notMet:
                continue
                # If we have got here then the CC does have the required DIVs at the required
                # rel abundance for the typeInQ.
            # Now must check to see if typeInQ covers more seqs than its initial type

            # Get tot abundance for current inital type
            # TODO 08/12/17
            # So here the question is 1 - does the CC already associate with a type of the basal sequence of the type
            # in question. If yes then we need to compare against this type. If no, then we can carry on with a
            # attempt at finding an association
            if CC.id in CCToInitialTypeDict.keys():
                # this will also check to make sure that the maj of pnt's basal type is found in the CC
                # if it is not then we should not let the CC associate with the potential new type
                currentInitialType = check_whether_pnt_needs_comparing_to_current_type(requirementDict ,analysis_type.objects.filter(id__in=CCToInitialTypeDict[CC.id]))
                # we will only allow a CC to support a potential new type if the pnt is of the same basal as one of the CCs current types
                # This is because we run into trouble in knowing which type to compare the pnt to if they are not the same basal type
                if currentInitialType:
                    # then the clade in question already associated to a type that has the same basal type as the type in question
                    # so we need to check if the new type represents more of the CCs sequences
                    # analysis_type.objects.get(id=CCToInitialTypeDict[CC.id][0])
                    currentTypeSeqRelAbundForCC = []
                    for refSeq in currentInitialType.getOrderedFootprintList():
                        relAbund = CCRelAbundDict[refSeq]
                        currentTypeSeqRelAbundForCC.append(relAbund)

                    # Get tot abundance for type In Q
                    typeInQAbundForCC = []
                    for refSeq in refSeqObjsOfTypeList:
                        relAbund = CCRelAbundDict[refSeq]
                        typeInQAbundForCC.append(relAbund)

                    if sum(typeInQAbundForCC) > sum(currentTypeSeqRelAbundForCC):
                        # Then this CC should be transfered to support the typeInQ
                        # For the time being we will simply hold these CCs in a list
                        # As we can then process all of the CCs from a given type at once
                        # because we will need to update and assess the types that they have come from
                        # to see whether they still have support
                        # Processing the CCs one by one would be far slower
                        supportList.append(CC)
                # else: # We should only allow a CC to support a type if they are differnt basal types.
                #     # supportList.append(CC)
            else:
                # if the CC doesn't currently have a type associated to it then there is no need to see
                # if the current type in Q is a better fit. we can simply add it to the supportList
                # Only allow the CC to give support to the type if the CC's Maj is in the potential type
                if CC.maj().referenceSequenceOf in [rs for rs, abund in requirementDict.items()]:
                    supportList.append(CC)
                # else do not allow CC to support potential new type
    return

def checkTypePairingForArtefactType(typeA, typeB, typefootprintdict, clade, cctocurrentinitialtypedict, cctototdict, cctorefabunddict, cores):
    # NB that the cctocurrentinitialtypedict is based on IDs rather than actual db objects


    ################## CREATE POTENTIAL NEW TYPE PROGRAMATICALLY ###############
    # Create the potential new type that will be tested. This will be made up of all of the combined intras of the two
    # types in question. The abundance requriement for those intras that are not artefact effected will still
    # be 0.03 but for those that are potentially artefact effected it will be lowered to 0.005



    totalListOfIntraIDs = set(typefootprintdict[typeA.id][1])
    totalListOfIntraIDs.update(typefootprintdict[typeB.id][1])


    totalListOfArtefactIntras = set(typefootprintdict[typeA.id][2])
    totalListOfArtefactIntras.update(typefootprintdict[typeB.id][2])


    # We will programatically represent the potential new type (pnt) as a list of tuples, one for each intra
    # first item in tuple will be the refseq of the intra and second item the required abundance for that intra
    pnt = []
    for intraID in totalListOfIntraIDs:
        if intraID in totalListOfArtefactIntras:
            requiredAbundance = unlockedAbund
        else:
            requiredAbundance = 0.03
        pnt.append((reference_sequence.objects.get(id=intraID), requiredAbundance))
    ############################################################################

    # Check to see if the potential newtypes footprint already exists
    pntFootprint = set([item[0].id for item in pnt])
    exists = False
    for key, footprintdictvalues in typefootprintdict.items():
        try:
            if footprintdictvalues[1] == pntFootprint:
                exists = True
                break
        except:
            apples = 'kjhg'
    if exists:
        print('Assessing new type:{}'.format(
            [reference_sequence.objects.get(id=refSeqID).name for refSeqID in pntFootprint]))
        print('Potential new type already exists')
        return False


    # Go through each of the CCs and see if the relative intra abundances meet the requirements of the pnt
    # if they do then add the id of the CC to the supportList
    # This will be used to count if there is sufficient support and to add the CCs to the new type once it is created
    # and to remove these CCs from the old types



    # this is causing a problem as some CCs are not being considered even though they would fit into this type
    # This is because they were not found in either of types A or B initially.
    # I think that really we should be looking through all CCs at this point. At least all CCs that are of the clade
    # The only thing we need to bear in mind is whether we would still enforce the rule of allowing one CC to only
    # support one initial type profile. I think we do need to enforce this rule else, you will have really basic
    # types being supported by CCs. e.g. C3 would get loads of support.
    # This could end up being very expensive so I would consider the following work flow to get this working.
    # Have a dictionary that is CC: current initial typeProfile found in. This will need to be kept upto date.
    # Then for each profile we are testing support for, get list of CCs of the clade, go through each CC as
    # before looking to see if there is support. When you find a CC that matches the requirements, look up the
    # type it was found in initially. Then only give support if current type in consideration represents more of its
    # seqs than its current type. If this is so then add this CC to a support type.
    # Once the list of CCs that support has been made, then check to see if support is great enough. If it is then
    # create the new type as before, but now go through each of the CCs and remove it from the type it was previously
    # associated to. Once you have removed it from the type, check to see if they type still has support.
    # If it is now below support, add the CC to a list of stranded CCs and delete the type
    # Once this is completed we have the new type and a list of stranded CCs. We then need to get the intras in common
    # for the stranded CCs (use 0.03 cutoff). If this profile already exists, then do nothing. else if it doesn't
    # already exisit then we can create a new type from the CC collection and footprint just identified.
    # when doing all of the above, be sure to keep the footprint dict up to date.
    # also the CC to type dict will need to be kept uptodate.

    # Get list of CCs in this analysis that are also of the clade in Q
    #TODO This still needs doing vvvv
    # This is currently a huge ammount of work to do going through all of these CCs
    # Would be good to cut this down



    #I'm making a change. We can include all CCs here and simply modify the workerOne code so that if a CC doesn't have an itnitial type
    # currently it can give support to the PNT if the PNT's DIVs are found in the CC.

    listOfCCsToCheck = [cc for cc in analysisObj.getCladeCollections() if cc.clade == typeA.clade if cc.id in cctocurrentinitialtypedict.keys()]

    # listOfCCsToCheck = [cc for cc in analysisObj.getCladeCollections() if cc.clade == typeA.clade if cc.id in cctocurrentinitialtypedict.keys()]

    supportList = []
    print('Assessing support for potential new type:{}'.format([reference_sequence.objects.get(id=refSeqID).name for refSeqID in pntFootprint]))


    ######################### CHECK EVERY CC OF THE ANALYSISOBJ ################################
    # To see whether the CC supports the pnt
    #TODO multiprocess this. Would need a managed list
    # Support list would be the the managed list


    ############# NEW MP CODE #############

    taskQueue = Queue()
    supportListManager = Manager()
    supportList = supportListManager.list()
    # outputQueue = Queue()


    for CC in listOfCCsToCheck:
        taskQueue.put(CC)

    for N in range(cores):
        taskQueue.put('STOP')

    allProcesses = []


    db.connections.close_all()

    for N in range(cores):
        p = Process(target=workerArtefactOne, args=(taskQueue, supportList, cctorefabunddict, pnt, cctocurrentinitialtypedict))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()


    #######################################



    ####################### IF PNT SUPPORTED CREATE PNT AND REDISTRIBUTE SUPPORTING CCs ############
    # Once the list of CCs that support has been made, then check to see if support is great enough. If it is then
    # create the new type as before, but now go through each of the CCs and remove it from the type it was previously
    # associated to. Once you have removed it from the type, check to see if the type still has support.
    # If it is now below support, add the CC to a list of stranded CCs and delete the type
    if len(supportList) >= 4:

        ##### CREATE NEW TYPE #########
        newAnalysisType = analysis_type(dataAnalysisFrom=analysisObj, clade=clade)
        # listOfCCs = [cc for cc in clade_collection.objects.filter(id__in=supportList)]
        newAnalysisType.initTypeAttributes(supportList, [pntItem[0] for pntItem in pnt])
        newAnalysisType.save()
        print('\nSupport found. Creating new type:{}'.format(newAnalysisType))
        # We need to keep the typefootprintdict upto date when types are created or deleted
        # get list of refseqs in type

        ####### UPDATE TYPEFOOTPRINT ########
        refSeqIDs = set([refSeq.id for refSeq in newAnalysisType.getOrderedFootprintList()])
        artefactIntraIDs = set([int(x) for x in newAnalysisType.artefactIntras.split(',') if x != ''])
        nonArtefactIDs = [id for id in refSeqIDs if id not in artefactIntraIDs]
        footprint = newAnalysisType.getOrderedFootprintList()
        typefootprintdict[newAnalysisType.id] = [nonArtefactIDs, refSeqIDs, artefactIntraIDs, footprint]


        ################# REDISTRIBUTE SUPPORTING TYPES ###############
        # Now remove the CCs from the types they were previously associated to and update the
        # cctocurrentinitialtypedict accordingly. # We will also likely need to reinitiate each of the types
        # There are likely to be types that have more than one CC being removed from them so the most effective
        # way to process this is not to go CC by CC but rather type by type
        # So first get a list of the types that have been effected
        setOfTypesAffected = set()
        # at the same time create a dict that tracks which CCs need to be remvoed from each type
        typeToCCToBeRemovedDict = defaultdict(list) # we may have to do this manually instead of defaultdict
        for CC in supportList:
            if CC.id in cctocurrentinitialtypedict.keys():
                # TODO 08/12/17 again we need to work out which of the types if there are multiples, was effected: DONE but check the new method check_which_type...
                initialTypeID = None
                if len(cctocurrentinitialtypedict[CC.id]) == 1:
                    initialTypeID = cctocurrentinitialtypedict[CC.id][0]
                else:
                    # this will return the analysis_type from the list of analysis types that matches the pnt
                    initial_type = check_which_type_has_same_basal_seq(pnt, analysis_type.objects.filter(id__in=cctocurrentinitialtypedict[CC.id]))
                    if initial_type:
                        initialTypeID = initial_type.id
                    else:
                        # If False then this means that the pnt didn't share any basal seqs in common with the current
                        # init types. In this case, we don't need to remove any of the current initTypes
                        # and we simply need to add the CC to the support of the type
                        initialTypeID = False
                if initialTypeID:
                    setOfTypesAffected.add(initialTypeID)
                    typeToCCToBeRemovedDict[initialTypeID].append(CC)
                else:
                    cctocurrentinitialtypedict[CC.id].append(newAnalysisType.id)
            else:
                # we don't need to remove this CC from any type but we should add it to the cctocurrentinitialtypedict
                cctocurrentinitialtypedict[CC.id] = [newAnalysisType.id]


        # Here we have a set of the type IDs for the types affected and a dict associating the CCs to these types
        # Now go through the types and remove the CCs from the type, change the dict association
        # and reinitialise the type

        # THe list to catch and reconsider all of the stranded CCs that might be left over in this process

        ############## TYPE BY TYPE OF REDISTRIBUTED CCs DELETE OR REINITIATE #################
        print('Reassessing support of affected types')
        strandedCCs = []
        for anType in typeToCCToBeRemovedDict.keys():
            anTypeInQ = analysis_type.objects.get(id=anType)

            #################### REMOVE CCs FROM TYPES ######################
            # Remove CCs from type
            listOfCCsToBeRemovedStrID = [str(cc.id) for cc in typeToCCToBeRemovedDict[anType]]
            anTypeInQ.removeCCListFromInitialCladeCollectionList(listOfCCsToBeRemovedStrID)

            # Change the CC associations in the cctocurrent... dict
            #TODO 08/12/17 remove and add to list DONE.
            for ccstrid in listOfCCsToBeRemovedStrID:
                cctocurrentinitialtypedict[int(ccstrid)].remove(anType)
                cctocurrentinitialtypedict[int(ccstrid)].append(newAnalysisType.id)

            ################### REASSESS SUPPORT OF TYPE ###################
            # Now check to see if the typeinQ still has sufficient support
            # If it does then reinitiate it
            ################ IF SUPPORTED ##################


            listOfCCsInType = [cc for cc in clade_collection.objects.filter(id__in=[int(x) for x in anTypeInQ.listOfCladeCollectionsFoundInInitially.split(',') if x != ''])]
            if len(listOfCCsInType) >= 4:
                print('Type {} supported by {} CCs. Reinitiating.'.format(anTypeInQ.name, len(listOfCCsInType)))
                # Then this still has suffiecient support
                # reinitiate the type
                # listOfCCs = [cc for cc in clade_collection.objects.filter(id__in=[int(x) for x in anTypeInQ.listOfCladeCollectionsFoundInInitially.split(',')])]
                anTypeInQ.initTypeAttributes(listOfCC=listOfCCsInType, footprintlistofrefseqs=anTypeInQ.getOrderedFootprintList())
            ################ IF UNSUPPORTED #################
            else:
                # THen this antype no longer has sufficient support
                # put the CCs into the stranded list and delete this type
                # also remove the CCs from the dict of types they are associated to
                print('Type {} no longer supported. Deleting. {} CCs stranded.'.format(anTypeInQ.name, str(len(listOfCCsInType))))

                del typefootprintdict[anTypeInQ.id]
                anTypeInQ.delete()
                for cc in listOfCCsInType:
                    if len(cctocurrentinitialtypedict[cc.id]) > 1:
                        # Then there are other types associated to this cladeCollection and we should simply
                        # remove the type in question from the list
                        cctocurrentinitialtypedict[cc.id].remove(anType)
                    else:
                        # then this only contains one type and we should delte the cc entry in the dict and it will
                        # become stranded.
                        del cctocurrentinitialtypedict[cc.id]
                        strandedCCs.append(cc)

        # Once this is completed we have the new type and a list of stranded CCs. We then need to get the intras in common
        # for the stranded CCs (use 0.03 cutoff). If this profile already exists, then do nothing. else if it doesn't
        # already exist then we can create a new type from the CC collection and footprint just identified.
        # when doing all of the above, be sure to keep the footprint dict up to date.
        # also the CC to type dict will need to be kept uptodate.

        # Find refSeqs in common
        # First get list of all refseqs above the cutoff in all CCs

        #################### ATTEMPT TO HOME STRANDED TYPES ######################
        # We will first check to see if there are intras in common between the CCs.
        # If there are then we will create a type from that if it doesn't already exist.
        # If there are no intras in common we will assign the CCs to existing types or
        # create maj types to house them if needs be.
        ################### GET COMMON INTRAS #################
        # Only attempt to make a new type to associate the stranded CCs to if there is sufficient support
        if len(strandedCCs) >= 4:
            totalIntraSet = set()
            for CC in strandedCCs:
                totalIntraSet.update(CC.cutOffFootprint(analysisObj.withinCladeCutOff))

            # Now go through each CC again and remove any intra from the total list that isn't found in the
            # CCinQ to produce an in common list
            refSeqsToRemove = set()
            for CC in strandedCCs:
                intrasInCCInQ = CC.cutOffFootprint(analysisObj.withinCladeCutOff)
                for refSeq in list(totalIntraSet):
                    if refSeq not in intrasInCCInQ:
                        refSeqsToRemove.add(refSeq)

            # Now create the commons list
            intrasInCommonList = [refSeq for refSeq in list(totalIntraSet) if refSeq not in refSeqsToRemove]

            exists = False
            if intrasInCommonList:
                ############ CHECK IF POTENTIAL NEW TYPE FOOTPRINT ALREADY EXISTS ############
                # Check to see if the potential newtypes footprint already exists
                pntFootprint = set([refSeq.id for refSeq in intrasInCommonList])

                typeThatExistsID = 0
                for key, footprintdictvalues in typefootprintdict.items():
                    if footprintdictvalues[1] == pntFootprint:
                        exists = True
                        typeThatExistsID = key
                        break

            ############ IF ALREADY EXISTS; ASSOCIATE STRANDED CCs TO EXISTING TYPE #################
            if exists:

                # If this type already exists then the CCsInQ will have a good type
                # to be assigned to and we need do no more
                # if this footprint already exists then we should associate these CCs to the type that has this footprint
                # so do this in the dictionary but also add these CCs to the cladeCollectionFoundInInitially list and
                # reinitiate the type
                associateCCsToExistingTypeAndUpdateDicts(cctocurrentinitialtypedict, strandedCCs, typeThatExistsID, typefootprintdict=typefootprintdict)




            ############ IF DOESN'T EXIST But INTRAS IN COMMONLIST DOES; MAKE NEW TYPE AND ASSOCIATED CCs ###########
            elif not exists and intrasInCommonList:

                # He we should create a new type based on the CC collection and footprint above
                # Create new type
                if not check_if_intrasInCommonList_contains_multiple_basal_seqs(intrasInCommonList):
                    makeNewTypeAndAssociateCCsAndUpdateDicts(cctocurrentinitialtypedict, clade, intrasInCommonList, strandedCCs,
                                               typefootprintdict)
        else:
            reassociateCCsToExistingTypesAndUpdateDicts(cctocurrentinitialtypedict, cctorefabunddict, clade, strandedCCs,
                                          typefootprintdict)
        return True

    ############## IF PNT INSUFFICIENT SUPPORT RETURN FALSE ################
    else:
        print('\nInsufficient support for potential new type')
        return False
    return

def check_which_type_has_same_basal_seq(pnt, list_of_type_objects):
    '''This should return the typeobject that has the same basal sequence as the pnt
    either C3, C15, C1 or None'''
    # first get the basal of the pnt
    basalpnt = False
    for rs, abund in pnt:
        basalpnt = False
        if 'C3' == rs.name:
            basalpnt = 'C3'
            break
        elif 'C1' == rs.name:
            basalpnt = 'C1'
            break
        elif 'C15' in rs.name:
            basalpnt = 'C15'
            break
    # here we know if pnt contains a basal seq

    # for each at check to see if the result is informative
    for at in list_of_type_objects:
        basal = False
        for rs in at.getOrderedFootprintList():
            if 'C3' == rs.name:
                basal = 'C3'
                break
            elif 'C1' == rs.name:
                basal = 'C1'
                break
            elif 'C15' in rs.name:
                basal = 'C15'
                break
        if basalpnt == basal:
            # then both types contain the same basal and so this is our type
            return at
    return False


def workerArtefactOne(input, supportList, cctorefabunddict, pnt, cctocurrentinitialtypedict):
    # pnt is a list of tuples, where 0 is a references sequence and 1 is the required relative abundance

    for CC in iter(input.get, 'STOP'):

        print('\rChecking {} {}'.format(CC, current_process().name), end='')
        ######### COMPARE REF SEQ ABUND INFO FOR CC TO PNT REQUIREMENTS###########
        # We use the cctorefabunddict to look up the relative abundances of the intras in question

        # Here we are checking to see if each CC has the relative abundance required to have the pnt
        # i.e. if each of the pnt's DIVs are found in the CC at suffiecient proportions

        # Holder to see if the requirements were met for this CC and the pnt
        met = True

        # listOfRefSeqsInCC = [dsss.referenceSequenceOf for dsss in data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=CC)]



        # Keep track of what proportion of the CCs sequences this type makes up
        pntSeqRelAbundForCC = []
        for intraReq in pnt:
            if intraReq[0] in cctorefabunddict[CC.id].keys():
                relAbund = cctorefabunddict[CC.id][intraReq[0]]
                pntSeqRelAbundForCC.append(relAbund)
                if relAbund < intraReq[1]:
                    # If we get here then one of the intrasRequirements in the pnt has not been met by
                    # the intra relative abundances in the CCinQ
                    met = False
                    break
            else:
                met = False
                break
        # met = were all requirements met



        ############# IF REQUIREMENTS MET CHECK TO SEE IF PNT SEQS REPRESENT MORE THAN CURRENT TYPE FOR CC #########
        if met: #PNT divs are found in CC
            # The proportion of the CC that the pnt under consideration represents
            pntSeqTotRelAbundForCC = sum(pntSeqRelAbundForCC)

            # When you find a CC that matches the requirements, look up the
            # type it was found in initially. Then only give support if current type in consideration represents more of its
            # seqs than its current type. If this is so then add this CC to a support type.
            # make dicts to speed this up.

            # if the CC in question, into which the PNT fits, already has a type associated to it
            # we need to consider whether this type covers more sequences than the PNT
            # if the PNT covers more, add support, else do not
            if CC.id in cctocurrentinitialtypedict.keys():
                # TODO 08/12/17 here we need a fix as the CCs can now support more than one type
                # we need to get the type that we should be comparing to.
                # If there are multiple types associated to a given CCt then they will be of different basal origins
                # we need to find the one that is the basal origin of the pnt
                # actually this is not true in the case of the C15s. It looks like we have a CC
                # that has multiple C15 types. So...
                # If the check_which_type_has_same_basal_seq returns false, then the non of the currentInitialTypes
                # have a basal in common with the pnt. In this case we can lend our support to the type
                currentInitialTypes = analysis_type.objects.filter(id__in=cctocurrentinitialtypedict[CC.id])
                currentInitialType = None

                if len(currentInitialTypes) == 1:
                    currentInitialType = currentInitialTypes[0]

                elif len(currentInitialTypes) > 1:
                    currentInitialType = check_which_type_has_same_basal_seq(pnt, currentInitialTypes)

                if currentInitialType:
                    currentTypeSeqRelAbundForCC = []
                    for refSeq in currentInitialType.getOrderedFootprintList():
                        relAbund = cctorefabunddict[CC.id][refSeq]
                        currentTypeSeqRelAbundForCC.append(relAbund)
                    currentTypeSeqTotRelAbundForCC = sum(currentTypeSeqRelAbundForCC)
                    # Check to see whether the potential new type covers more of the sequence in the CC.
                    # If so then this CC will now support the new type and no longer support its original initial type

                    ####### PNT SUPPORT #########
                    if pntSeqTotRelAbundForCC > currentTypeSeqTotRelAbundForCC:
                        # Add the CC to the supportList for the pnt
                        supportList.append(CC)
                        # We can't currently change the cctocurrentititialtypedict as we need the pnt types id and we
                        # don't do that until we know it is supported.
                    ####### NO PNT SUPPORT #######
                    else:
                        # If it doesn't then nothing changes
                        pass

                else:
                    # Then the pnt doesn't have a basal seq in common with any of the current intitial types
                    supportList.append(CC)

            # else, if the CCinq doesn't already have a type associated then I see no problem with lending its support
            # to the PNT
            else:
                supportList.append(CC)
    return

def associateCCsToExistingTypeAndUpdateDicts(cctocurrentinitialtypedict, strandedCCs, typeThatExistsID, typefootprintdict):
    for CC in strandedCCs:
        cctocurrentinitialtypedict[CC.id] = [typeThatExistsID]
    try:
        typeThatExists = analysis_type.objects.get(id=typeThatExistsID)
    except:
        apples = 'asdf'
    print('Associating stranded types to existing type {}'.format(typeThatExists.name))
    # add the CCs to the type
    currentListOfInitialCCs = typeThatExists.listOfCladeCollectionsFoundInInitially.split(',')
    currentListOfInitialCCs.extend(str(CC.id) for CC in strandedCCs)
    typeThatExists.listOfCladeCollectionsFoundInInitially = ','.join(currentListOfInitialCCs)
    typeThatExists.save()
    # Re initiate the type in Q to incorporate the new CCs
    listOfCCs = [cc for cc in clade_collection.objects.filter(
        id__in=[int(x) for x in typeThatExists.listOfCladeCollectionsFoundInInitially.split(',')])]
    typeThatExists.initTypeAttributes(listOfCC=listOfCCs,
                                      footprintlistofrefseqs=typeThatExists.getOrderedFootprintList())
    # Update the typefootprintdict with this type
    refSeqIDs = set([refSeq.id for refSeq in typeThatExists.getOrderedFootprintList()])
    artefactIntraIDs = set([int(x) for x in typeThatExists.artefactIntras.split(',') if x != ''])
    nonArtefactIDs = [id for id in refSeqIDs if id not in artefactIntraIDs]
    footprint = typeThatExists.getOrderedFootprintList()
    typefootprintdict[typeThatExists.id] = [nonArtefactIDs, refSeqIDs, artefactIntraIDs, footprint]


def makeNewTypeAndAssociateCCsAndUpdateDicts(cctocurrentinitialtypedict, clade, intrasInCommonList, strandedCCs, typefootprintdict):
    newAnalysisType = analysis_type(dataAnalysisFrom=analysisObj, clade=clade)
    listOfCCs = list(strandedCCs)
    listOfRefSeqs = intrasInCommonList
    newAnalysisType.initTypeAttributes(listOfCCs, listOfRefSeqs)

    # newAnalysisType.save()
    print('Creating new type: {} from {} residual cladeCollections'.format(newAnalysisType, len(strandedCCs)))
    # Update the typefootprintdict with this type
    refSeqIDs = set([refSeq.id for refSeq in newAnalysisType.getOrderedFootprintList()])
    artefactIntraIDs = set([int(x) for x in newAnalysisType.artefactIntras.split(',') if x != ''])
    nonArtefactIDs = [id for id in refSeqIDs if id not in artefactIntraIDs]
    footprint = newAnalysisType.getOrderedFootprintList()
    typefootprintdict[newAnalysisType.id] = [nonArtefactIDs, refSeqIDs, artefactIntraIDs, footprint]
    # Also update the cctocurrentinitialtypedict
    for CC in listOfCCs:
        cctocurrentinitialtypedict[CC.id] = [newAnalysisType.id]
    return newAnalysisType.name

def reassociateCCsToExistingTypesAndUpdateDicts(cctocurrentinitialtypedict, cctorefabunddict, clade, strandedCCs, typefootprintdict):

    # Here we have a number of CCs in the stranded list < 4. We need to work out what to do with them
    # I think at this point we don't want to make any inference with them
    # We will make sure that when we are getting the list of CCs to look for support in that we will not
    # include the CCs that are currently hanging
    # not having CCs reasigned to a type is causing problems.
    # Due to the fact that strict limits are being used to assign types to the discovered types
    # it is meaning that these stranded CCs are getting bad associations
    # We need to reassociate them to the best possible types
    # For each CC we are going to do a sort of mini type assignment

    #DONE because a lot of the single intra types will have been gotten rid of at this point
    # there is a possibility that there will not be a type for the CCs to fit into
    # Also it may be that the CCs now fit into a lesser intra single intra type. e.g. B5c if original type was B5-B5s-B5c.
    # In this case we would obviously rather have the type B5 associated with the CC.
    # So we should check to see if the abundance of the type that the CC has been found in is larger than the abundance of the CCs Maj intra.
    # If it is not then we should simply create a new type of the maj intra and associate the CC to that.
    # Also we need to be mindful of the fact that a CC may not find a match at all, e.g. the bestTypeID will = 'None'.
    # In this case we will also need to make the Maj intra the type.
    if strandedCCs:
        print('Insufficient stranded CCs ({}) to support a new type. Reassociating...'.format(len(strandedCCs)))


        for CC in strandedCCs:
            # Get new list of an types for each CC in case types are made
            listOfAnalysisTypes = [antype for antype in
                                   analysis_type.objects.filter(id__in=list(typefootprintdict.keys()))
                                   if antype.clade == clade]
            highAbund = 0
            bestTypeID = None
            refSeqAbundDictForCC = cctorefabunddict[CC.id]
            for antype in listOfAnalysisTypes:
                footprintDictOfType = typefootprintdict[antype.id]
                # for each refseq in the types footprint
                match = True
                abundOfTypeSeqs = []
                for refSeq in footprintDictOfType[3]:
                    if refSeq in refSeqAbundDictForCC.keys():
                        abundOfRefSeq = refSeqAbundDictForCC[refSeq]
                        if refSeq.id in footprintDictOfType[0]:
                            # Then this is an nonartefact refSeq and so must occur above 0.03
                            if abundOfRefSeq > 0.03:
                                abundOfTypeSeqs.append(abundOfRefSeq)
                            else:
                                # THen this refseq is not found in high enough abundance and so this type cannot be
                                # associated to this CC
                                match = False
                                break
                        else:
                            # THen this is an artefact refSeq and so must occur only above 0.005
                            if abundOfRefSeq > unlockedAbund:
                                abundOfTypeSeqs.append(abundOfRefSeq)
                            else:
                                # Then this refseq is not found in high enough abundance and so this type cannot be
                                # associated to this CC
                                match = False
                                break
                    # If any of the refseqs aren't found in the CC then this type cannot be associated to the CC
                    else:
                        match = False
                        break
                if match:
                    # Then the type fits inside the CC
                    if sum(abundOfTypeSeqs) > highAbund:
                        highAbund = sum(abundOfTypeSeqs)
                        bestTypeID = antype.id
            # When we get here we have looked at all of the antypes for this CC
            # Associate the CC with the type found.

            # Need to check to see if we have found a type that contains more of the CCs seqs than just the Maj Intra
            mostAbundIntraOfCC = max(refSeqAbundDictForCC.items(), key=operator.itemgetter(1))
            abundance = mostAbundIntraOfCC[1]
            refSeq = mostAbundIntraOfCC[0]

            if bestTypeID != None and highAbund >= abundance:
                matchType = analysis_type.objects.get(id=bestTypeID)
                currentListOfInitialCCs = matchType.listOfCladeCollectionsFoundInInitially.split(',')
                currentListOfInitialCCs.append(str(CC.id))
                matchType.listOfCladeCollectionsFoundInInitially = ','.join(currentListOfInitialCCs)
                matchType.save()
                # update type dict
                cctocurrentinitialtypedict[CC.id] = [matchType.id]

                # reinitiate type
                listOfCCs = [cc for cc in
                             clade_collection.objects.filter(id__in=[int(x) for x in currentListOfInitialCCs])]
                matchType.initTypeAttributes(listOfCC=listOfCCs,
                                             footprintlistofrefseqs=matchType.getOrderedFootprintList())

                # update the typeFootPrintDict
                refSeqIDs = set([refSeq.id for refSeq in matchType.getOrderedFootprintList()])
                artefactIntraIDs = set([int(x) for x in matchType.artefactIntras.split(',') if x != ''])
                nonArtefactIDs = [id for id in refSeqIDs if id not in artefactIntraIDs]
                footprint = matchType.getOrderedFootprintList()
                typefootprintdict[matchType.id] = [nonArtefactIDs, refSeqIDs, artefactIntraIDs, footprint]


                print('Stranded CC {} reassociated to {}'.format(CC, matchType.name))
            else:
                # Here we need to create a new type that is simply the Maj intra and associate the CC to it
                newAnalysisTypeName = makeNewTypeAndAssociateCCsAndUpdateDicts(cctocurrentinitialtypedict=cctocurrentinitialtypedict, clade=clade, intrasInCommonList=[refSeq], strandedCCs=[CC], typefootprintdict=typefootprintdict)
                # print('Stranded CC {} associated to new type {}'.format(CC, newAnalysisTypeName))

            #### DELETE ME #####
            # print('Stranded CC {} reassociated to {}'.format(CC, matchType.name))
    return

def createTotalSeqsDictForAllCCs(cores):
    # Generate a dict that simply holds the total number of seqs per CC
    # This will be used when working out relative proportions of seqs in the CC
    # I am going to make sure that this directory already exists.
    # This way I shouldn't have to create it each time

    testDirPath = os.path.join(os.path.dirname(__file__), 'temp/{}'.format(analysisObj.id))
    os.makedirs(testDirPath, exist_ok=True)
    os.chdir(testDirPath)
    try:
        # See if this dict has already been created for this analysis
        # print('We got here')
        CCToTotalSeqsDict = pickle.load(open("CCToTotalSeqsDict_{}".format(analysisObj.id), "rb"))
        print('Loading CCToTotalSeqDict')
        # print('But we also got here. We got here')
    except:
        # if not generate from scratch
        # print('We got here')
        print('Generating CCToTotalSeqDict')

        listOfCCsToCheck = [cc for cc in analysisObj.getCladeCollections()]


        taskQueue = Queue()
        CCToTotalSeqstManager = Manager()
        CCToTotalSeqsDict = CCToTotalSeqstManager.dict()
        #outputQueue = Queue()


        for CC in listOfCCsToCheck:
            taskQueue.put(CC)

        for N in range(cores): # TODO Change back to cores
            taskQueue.put('STOP')

        allProcesses = []

        # In this case the db.connections were not being recreated and I think this might have something
        # to do with a connection only being created when there is writing to be done.
        # as there was only reading operations to be done in the worker task I think that maybe no new connections
        # were being made.
        db.connections.close_all()
        # I think this might be being caused due to the changing of the address of the db in the settings.py file to a relative path
        # yep this is what was causing the problem so we will need to find some way of defining the location of the database relatively.
        # DONE. I have changed the settings.py file so that the DB location is now relative to the position of the settings.py file


        for N in range(cores):
            p = Process(target=workerCCToTotalSeqsDict, args=(taskQueue, CCToTotalSeqsDict))
            allProcesses.append(p)
            p.start()

        for p in allProcesses:
            p.join()
        # For the time being I will not dump here but dump all four dicts at the end
        # pickle.dump(CCToTotalSeqsDict, open("CCToTotalSeqsDict_{}".format(analysisObj.id), "wb"))

        return CCToTotalSeqsDict

def createRefSeqRelAbundsForAllCCsDict(CCToTotalSeqsDict, cores):
    # Generate dict per cc for listof reference sequences and their abundances in the CC
    listOfCCsToCheck = [cc for cc in analysisObj.getCladeCollections()]
    testDirPath = os.path.join(os.path.dirname(__file__), 'temp/{}'.format(analysisObj.id))
    os.chdir(testDirPath)
    try:
        # See if this dict has already been created for this analysis
        CCToRefSeqListAndAbundances = pickle.load(open("CCToRefSeqListAndAbundances_{}".format(analysisObj.id), "rb"))

        return CCToRefSeqListAndAbundances
    except:
        # if not generate from scratch
        # CCToRefSeqListAndAbundances = {}
        # print('Generating CCToRefSeqListAndAbundaces')
        # for CC in listOfCCsToCheck:
        #     print(CC)
        #     listOfDSSSInCC = [dsss for dsss in data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=CC)]
        #     listOfRefSeqsInCC = [dsss.referenceSequenceOf for dsss in listOfDSSSInCC]
        #     listOfAbundances = [dsss.abundance / CCToTotalSeqsDict[CC.id] for dsss in listOfDSSSInCC]
        #     innerDict = {}
        #     for i in range(len(listOfDSSSInCC)):
        #         innerDict[listOfRefSeqsInCC[i]] = listOfAbundances[i]
        #     CCToRefSeqListAndAbundances[CC.id] = innerDict

        #Start the multithreading

        taskQueue = Queue()
        CCToRefSeqListAndAbundancesManager = Manager()
        CCToRefSeqListAndAbundancesDict = CCToRefSeqListAndAbundancesManager.dict()



        for CC in listOfCCsToCheck:
            taskQueue.put(CC)

        for N in range(cores):
            taskQueue.put('STOP')

        allProcesses = []


        db.connections.close_all()


        for N in range(cores):
            p = Process(target=workerCCToRefSeqListAndAbundances, args=(taskQueue, CCToTotalSeqsDict, CCToRefSeqListAndAbundancesDict))
            allProcesses.append(p)
            p.start()

        for p in allProcesses:
            p.join()
        # For the time being I will not pcikle here but instead pickle all four of the dicts when they are done
        # pickle.dump(CCToRefSeqListAndAbundancesDict, open("CCToRefSeqListAndAbundances_{}".format(analysisObj.id), "wb"))

        return CCToRefSeqListAndAbundancesDict

def workerCCToRefSeqListAndAbundances(input, CCToTotalSeqsDict, CCToRefSeqListAndAbundances):
    for CC in iter(input.get, 'STOP'):
        print('\r{} {}'.format(CC, current_process().name), end='')
        listOfDSSSInCC = [dsss for dsss in data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=CC)]
        listOfRefSeqsInCC = [dsss.referenceSequenceOf for dsss in listOfDSSSInCC]
        listOfAbundances = [dsss.abundance / CCToTotalSeqsDict[CC.id] for dsss in listOfDSSSInCC]
        innerDict = {}
        for i in range(len(listOfDSSSInCC)):
            innerDict[listOfRefSeqsInCC[i]] = listOfAbundances[i]
        CCToRefSeqListAndAbundances[CC.id] = innerDict

def createTypeArtefactDivInfoDict(allTypesFromDataAnalysis, cores):
    # Create dict for each type of refSeqs, artefact and non-artefact DIVs, and footprint
    typeFootPrintDict = {}
    for at in allTypesFromDataAnalysis:
        # if str(at) == 'C42g-C42a-C42.2-C42h-C1-C42b':
        #     apples = 'pears'
        # get list of refseqs in type
        refSeqIDs = set([refSeq.id for refSeq in at.getOrderedFootprintList()])
        artefactIntraIDs = set([int(x) for x in at.artefactIntras.split(',') if x != ''])
        # ###DEBUG###
        # for ID in artefactIntraIDs:
        #     print(reference_sequence.objects.get(id=ID).name)
        # #####
        nonArtefactIDs = [id for id in refSeqIDs if id not in artefactIntraIDs]
        footprint = at.getOrderedFootprintList()
        typeFootPrintDict[at.id] = [nonArtefactIDs, refSeqIDs, artefactIntraIDs, footprint]

    return typeFootPrintDict

def createCCtoInitialTypeDict(allTypesFromDataAnalysis):
    #TODO 08/12/17
    # I think I need to change this so that for a single CCID key we have a value that is a list of typeIDs
    # as now we can have multiple inital types supported by one CC
    CCToInitialTypeDict = defaultdict(list)
    CCtoTypeTupleList = []
    for at in allTypesFromDataAnalysis:
        typeID = at.id
        initialCCs = [int(x) for x in at.listOfCladeCollectionsFoundInInitially.split(',')]
        for CCID in initialCCs:
            CCtoTypeTupleList.append((CCID, typeID))

    # Here we have every combination of type and CCID
    for ccid, atid in CCtoTypeTupleList:
        CCToInitialTypeDict[ccid].append(atid)
    #TODO 08/12/17 then check to see if this dict is used else where as if it is then we will need to make adjustments



    return dict(CCToInitialTypeDict)

def check_if_type_pairs_contain_incompatible_basal_seqs(typea, typeb):
    '''The aim of this function is simply to identify types thwhether types contain seqs that are
    incombpatible basals, ie. if a contins C3 and then b contains either c1 or c15 then we don't want
    them to be considered together'''
    # TODO 11/01/18 we can potentially speed this up by keeping the initalT information on basal seqs
    # and having this info be held in the type object.
    refSeqsA = typea.getOrderedFootprintList()
    refSeqsB = typeb.getOrderedFootprintList()

    #first check to see if we have a basal seq in a.
    # If we don't then we can't have an compataility
    basalInA = []
    foundC15A = False
    for rs in refSeqsA:
        if rs.name == 'C3':
            basalInA.append('C3')
        elif rs.name == 'C1':
            basalInA.append('C1')
        elif 'C15' in rs.name and foundC15A == False:
            basalInA.append('C15')
            foundC15A = True

    if len(basalInA) == 0:
        return False

    basalInB = []
    foundC15B = False
    for rs in refSeqsB:
        if rs.name == 'C3':
            basalInB.append('C3')
        elif rs.name == 'C1':
            basalInB.append('C1')
        elif 'C15' in rs.name and foundC15B == False:
            basalInB.append('C15')
            foundC15B = True

    # Now we have the basal seqs found in each of the types - it should only ever be one
    # find non-common items. If there are any, then we know that we have an incompatability
    if len(set(basalInA) ^ set(basalInB)) > 0:
        return True
    return False



def checkForAdditionalArtefactTypes(cores):
    # Generate a dict that simply holds the total number of seqs per CC
    # This will be used when working out relative proportions of seqs in the CC

    CCToTotalSeqsDict = createTotalSeqsDictForAllCCs(cores)


    # Generate dict per cc for listof reference sequences and their abundances in the CC
    CCToRefSeqListAndAbundances = createRefSeqRelAbundsForAllCCsDict(CCToTotalSeqsDict, cores)

    # Get list of all types
    allTypesFromDataAnalysis = analysis_type.objects.filter(dataAnalysisFrom=analysisObj)

    # Get list of clades that are represented by the types
    cladeList = set()
    for at in allTypesFromDataAnalysis:
        cladeList.add(at.clade)

    # Create dict for each type of refSeqs, artefact and non-artefact DIVs, and footprint
    typeFootPrintDict = createTypeArtefactDivInfoDict(allTypesFromDataAnalysis, cores)

    # TODO 08/12/17 this is going to need to be looked at as each CC can now have multiple initial types
    # Create a dict that is CC: inital type found in
    CCToInitialTypeDict = createCCtoInitialTypeDict(allTypesFromDataAnalysis)

    # Do pairwise comparison of types within clades
    #TODO 08/12/17
    # I have gone through all of the remaining code in this function and fixed the cc to initial type dict
    # so that it is compatible with ccts having multiple types.
    # we still need to do the other artefact checking
    for clade in cladeList:
        allTypesFromDataAnalysis = analysis_type.objects.filter(dataAnalysisFrom=analysisObj, clade=clade)

        staticListOfTypes = list(analysis_type.objects.filter(dataAnalysisFrom=analysisObj, clade=clade))
        doneList = []
        while 1:
            restart = False



            # Pairwise comparison of types within clade

            # For implementing multiprocessing in this part of the code
            # We may be able to create a list of all of the pairwise comparisons and then add this to a queue
            # In essence this would mean that all of the code below the a,b iter would go into the worker.
            # However, it might be that the time limiting part of all of this is going to be the checkTypePairingForArtefactType
            # So, perhaps we should look to see whether this component can be multiprocessed.
            # BUT because we may be deleting in and adding types this might be a problem
            # so it might be best if we can try to speed up the internals of the checkTypePa...
            for a, b, in itertools.combinations([at for at in allTypesFromDataAnalysis if at in staticListOfTypes], 2):
                # Check to see if the pairs have already been done as we will be running through this process several
                # times in all likelihood


                if (a.id, b.id) not in doneList:
                    listOfNonArtefactIntrasA = typeFootPrintDict[a.id][0]
                    listOfNonArtefactIntrasB = typeFootPrintDict[b.id][0]

                    # TODO 08/12/17 here we can simply check to see if the pair of types contain incompatible types
                    # if they do then we don't compare them and simply add to the done list
                    # returns False if it is OK to compare, i.e. non incompatabiites
                    # returns True if there are incompatabilities and they cannot be compared
                    if not check_if_type_pairs_contain_incompatible_basal_seqs(a, b):

                        # If all of the nonArtefact intras of A are found in the footprint of B
                        # and if all of the nonArtefact intras of B are found in footprint of A.
                        # i.e. if only artefact intras differentiate the types
                        if set(listOfNonArtefactIntrasA).issubset(typeFootPrintDict[b.id][1]) and set(listOfNonArtefactIntrasB).issubset(typeFootPrintDict[a.id][1]):
                            # Check to see whether either a or b are subsets of each other
                            aSubsetB = typeFootPrintDict[a.id][1].issubset(typeFootPrintDict[b.id][1])
                            bSubsetA = typeFootPrintDict[b.id][1].issubset(typeFootPrintDict[a.id][1])

                            if not aSubsetB and not bSubsetA:
                                listOfArtefactIntrasInBothTypes = list(a.artefactIntras.split(','))
                                listOfArtefactIntrasInBothTypes.extend(b.artefactIntras.split(','))
                                if listOfArtefactIntrasInBothTypes:
                                    # Here we have finally found a type pairing that needs to be checked
                                    print('\nChecking {} and {} for additional artefactual profiles'.format(a, b))
                                    if checkTypePairingForArtefactType(a, b, typeFootPrintDict, clade, CCToInitialTypeDict, CCToTotalSeqsDict, CCToRefSeqListAndAbundances, cores):
                                        # If we did find a new type then we need to start the a,b comparison again
                                        # It should be quick though due to the doneList
                                        restart = True
                                        # We need to recalculate this as we may have delted and created types
                                        # DONE this may be stuck infinitely looping by searching the new types created.
                                        # Maybe only work with the types we started with.
                                        # We could make a list of original types and then just look through the types
                                        # that are found in common between the original types list and an
                                        # updated types list. This way no newly added types will be analysed.
                                        # and we should avoid a run away situation
                                        allTypesFromDataAnalysis = analysis_type.objects.filter(dataAnalysisFrom=analysisObj, clade=clade)
                                        break
                                    else:
                                        doneList.extend([(a.id, b.id), (b.id, a.id)])

                                else:
                                    doneList.extend([(a.id, b.id), (b.id, a.id)])
                            else:
                                doneList.extend([(a.id, b.id), (b.id, a.id)])
                        else:
                            # Whenever we fall out of the series of conditionals we should add both possible combinations
                            # of type a and type b to the done list so that we don't waste time goig through them again
                            doneList.extend([(a.id, b.id), (b.id, a.id)])
                    else:
                        # then they have incompatable basal seqs and must not be compared
                        doneList.extend([(a.id, b.id), (b.id, a.id)])
            # If we make it to here we did a full run through the types without making a new type
            # Time to do the same for the next type
            if not restart:
                break

        # ### DEBUG - See which CCs are associated to which types now ###
        # for at in analysis_type.objects.filter(dataAnalysisFrom=analysisObj, clade=clade):
        #     print('{}'.format(at))
        #     print('{}'.format(','.join([str(cct) for cct in at.getCladeCollectionsFoundInInitially()])))
        #     if 'January2017-280por' in [str(cct) for cct in at.getCladeCollectionsFoundInInitially()]:
        #         appeles = 'pwers'
        # ###############################################################

    return CCToTotalSeqsDict, CCToRefSeqListAndAbundances, typeFootPrintDict, CCToInitialTypeDict

def workerCCToTotalSeqsDict(input, cctototalseqsdict):
    for CC in iter(input.get, 'STOP'):


        print('\r{} {}'.format(CC, current_process().name), end='')
        totSeqsInCC = sum([dsss.abundance for dsss in data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=CC)])
        cctototalseqsdict[CC.id] = totSeqsInCC

def workerDiscoveryTwoWorker(input, output, withncladecutoff):
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    for cladecollection in iter(input.get, 'STOP'):
        clade = cladecollection.clade
        footPrintInQ = cladecollection.cutOffFootprint(withncladecutoff)

        CCMaj = cladecollection.maj()
        # pass queuetwo (footprintInQ, index, CC, CC.maj())
        output.put((footPrintInQ, cladeList.index(clade), cladecollection, CCMaj))
        print('\rFound footprint {}'.format(footPrintInQ), end='')
        # print('Found footprint {}'.format(footPrintInQ))
    output.put('kill')

def check_if_contains_multiple_basal_seqs(footprintdict, keyname):
    ''' The object of this function is to check a given set of referncesequence objects and see if it contains multiple
    basal sequence objects
    if it does return false
    else return true
    just FYI C115 is a C3 derivative

    We will also modify the dictionary. If we find that there are multiple basals in the footprint then we will add a
    maj for each of the basals, for the C3 and C1 this will simply be the C3 or C1 sequence but for the C15 it
    needs to be the most abundnat C15x for each of the samples.
    '''
    basalSeqNames = ['C3', 'C15', 'C1']
    C15Found = False
    count = 0
    for rs in keyname:
        if rs.name == 'C3':
            count += 1
            # Here we need to check that for each of the CCs we have a C3 dsss in the list of maj seqs
            cctList = footprintdict[keyname][0]
            for i in range(len(cctList)):
                c3dsss = data_set_sample_sequence.objects.get(cladeCollectionTwoFoundIn=cctList[i], referenceSequenceOf=rs)

                #TODO check that whether it is a list that we are adding to
                # this will probably fail because maj list is currently just a list.
                # we need to change this to be a 2Dlist we will need to update this when we are first forming the dict
                if c3dsss not in footprintdict[keyname][1][i]:
                    footprintdict[keyname][1][i].append(c3dsss)
        elif rs.name == 'C1':
            count += 1
            # Here we need to check that for each of the CCs we have a C3 dsss in the list of maj seqs
            cctList = footprintdict[keyname][0]
            for i in range(len(cctList)):
                c1dsss = data_set_sample_sequence.objects.get(cladeCollectionTwoFoundIn=cctList[i],
                                                              referenceSequenceOf=rs)
                # TODO check that whether it is a list that we are adding to
                # this will probably fail because maj list is currently just a list.
                # we need to change this to be a 2Dlist we will need to update this when we are first forming the dict
                if c1dsss not in footprintdict[keyname][1][i]:
                    footprintdict[keyname][1][i].append(c1dsss)
        elif 'C15' in rs.name and C15Found == False: # we don't want this to fail just cause we found multiple C15s
            # so only count a C15 once
            count += 1
            C15Found = True
            # This is a little more tricky than the C1 and the C3
            cctList = footprintdict[keyname][0]
            for i in range(len(cctList)):
                # for each cct we need to find the most abundant C15X seq and add this as the majsequence
                listOfSeqsInCCT = list(data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=cctList[i]).order_by('-abundance'))
                for dsss in listOfSeqsInCCT:
                    if 'C15' in dsss.referenceSequenceOf.name:
                        # then this is the seq we want and we want to break after
                        if dsss not in footprintdict[keyname][1][i]:
                            footprintdict[keyname][1][i].append(dsss)
                        break
    if count > 1:
        return False
    else:
        return True

def check_if_contains_multiple_basal_majs(set_of_maj_ref_seqs_large):
    count = 0
    extract = False
    foundC15 = False
    for rs in set_of_maj_ref_seqs_large:
        if rs.name == 'C3':
            count +=1
        elif rs.name == 'C1':
            count +=1
        elif 'C15' in rs.name and not foundC15:
            count += 1
            foundC15 = True
    if count > 1:
        extract = True

    return extract

class initalType():
    # when we collapse an initialType profile we should be careful to remove the maj seqs that are of a reference
    # seq type that are found in the profile of the type we are collapsing to
    def __init__(self, refSeq_set, cladeCollection_list, maj_dsss_list=False):
        self.profile = refSeq_set
        self.profile_length = len(self.profile)
        self.contains_multiple_basal_sequences, \
        self.basalSequence_list = self.check_if_initialType_contains_basal_sequences()
        self.cladeCollection_list = cladeCollection_list
        self.support = len(self.cladeCollection_list)
        # We may move away from using the dsss but for the time being we will use it
        if maj_dsss_list:
            self.majority_sequence_list, self.set_of_maj_ref_seqs = self.create_majority_sequence_list_for_initialType(maj_dsss_list)
        else:
            self.majority_sequence_list, self.set_of_maj_ref_seqs = self.create_majority_sequence_list_for_initalType_from_scratch()

    def __repr__(self):
        return str(self.profile)

    def check_if_initialType_contains_basal_sequences(self):
        '''This function will return two items, firstly a list a bool if there are multiple basal sequences contained
        within the profile_set and secondly it will return a list of the
        I will just check the profile sequence '''
        basal_seq_list = []
        foundC15 = False
        for rs in self.profile:
            if rs.name == 'C3':
                basal_seq_list.append('C3')
            elif rs.name == 'C1':
                basal_seq_list.append('C1')
            elif 'C15' in rs.name and not foundC15:
                basal_seq_list.append('C15')
                foundC15 = True

        if len(basal_seq_list) > 1:
            return True, basal_seq_list
        else:
            return False, basal_seq_list

    def substract_init_type_from_other_init_type(self, other_init_type):
        self.profile = self.profile.difference(other_init_type.profile)
        self.profile_length = len(self.profile)
        self.basalSequence_list = list(set(self.basalSequence_list).difference(set(other_init_type.basalSequence_list)))
        if len(self.basalSequence_list) > 1:
            self.contains_multiple_basal_sequences = True
        else:
            self.contains_multiple_basal_sequences = False
        self.majority_sequence_list, self.set_of_maj_ref_seqs = self.create_majority_sequence_list_for_initalType_from_scratch()

    def create_majority_sequence_list_for_initialType(self, maj_dsss_list):
        # I'm trying to remember what form this takes. I think we'll need to be looking
        # This should be a list of lists. There should be a list for each cladeCollection in the self.cladeCollection_list
        # Within in each of the lists we should have a list of dataSetSampleSequence objects
        # We should already have a list of the dsss's with one dsss for each of the cladeCollections found in the maj_dsss_list
        # we will look to see if there are multiple basal sequences
        # If there are multiple basal sequences then for each cladeCollection within the intial type we will ad
        # the dsss to the list. If there are not multiple basal sequences then we will simply add the dss to a list
        setOfMajRefSeqs = set()
        master_dsss_list = []
        if self.contains_multiple_basal_sequences:
            for clade_collection_obj in self.cladeCollection_list:
                temp_dsss_list = []
                dssstList = list(data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=clade_collection_obj).order_by('-abundance'))
                # for each of the basal seqs in the basal seqs list, find the dsss representative
                for basal_seq in self.basalSequence_list:
                    if basal_seq == 'C15':
                    # Then we just need to find the most abundnant dsss that's name contains the C15

                        for dsss in dssstList:
                            if 'C15' in dsss.referenceSequenceOf.name:
                                temp_dsss_list.append(dsss)
                                # Important to break so that we only add the first and most abundant C15 seq
                                break
                    else:
                        # then we are looking for exact matches
                        for dsss in dssstList:
                            if dsss.referenceSequenceOf.name == basal_seq:
                                temp_dsss_list.append(dsss)
                                break
                # We should also make sure that the original maj sequence is found in the list
                if dssstList[0] not in temp_dsss_list:
                    temp_dsss_list.append(dssstList[0])
                # Make sure that each of the refSeqs for each of the basal or majs are in the maj set
                for dss in temp_dsss_list:
                    setOfMajRefSeqs.add(dss.referenceSequenceOf)
                # Here we should have a list of the dsss instances that represent the basal sequences for the CC in Q
                master_dsss_list.append(temp_dsss_list)
        else:
            # Then there is ony one basal sequence in this initial type and so we simply need to surround the maj
            # with a list.
            for i in range(len(self.cladeCollection_list)):
                master_dsss_list.append(maj_dsss_list[i])
                setOfMajRefSeqs.add(maj_dsss_list[i][0].referenceSequenceOf)


        return master_dsss_list, setOfMajRefSeqs

    def create_majority_sequence_list_for_initalType_from_scratch(self):
        # This will be like above but will not start with the maj_dsss_list
        # we will go through each of the cladeCollections of the type and get the maj sequence for the type

        #if the init type has multiple basal sequences then we will have to find the actual maj and the basal seq dsss
        setOfMajRefSeqs = set()
        master_dsss_list = []
        if self.contains_multiple_basal_sequences:
            for clade_collection_obj in self.cladeCollection_list:
                temp_dsss_list = []
                dsss_in_cc = list(data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=clade_collection_obj).order_by('-abundance'))
                dsss_in_cc_in_profile = [dsss for dsss in dsss_in_cc if dsss.referenceSequenceOf in self.profile]
                # first find the dsss that are the representatives of the basal types
                for basal_seq in self.basalSequence_list:
                    if basal_seq == 'C15':
                    # Then we just need to find the most abundnant dsss that's name contains the C15

                        for dsss in dsss_in_cc_in_profile:
                            if 'C15' in dsss.referenceSequenceOf.name:
                                temp_dsss_list.append(dsss)
                                # Important to break so that we only add the first and most abundant C15 seq
                                break
                    else:
                        # then we are looking for exact matches
                        for dsss in dsss_in_cc_in_profile:
                            if dsss.referenceSequenceOf.name == basal_seq:
                                temp_dsss_list.append(dsss)
                                break
                # now add the actual maj dsss if not one of the basal seqs
                basal_dsss = dsss_in_cc_in_profile[0]
                if basal_dsss not in temp_dsss_list:
                    # Then the actual maj is not already in the list
                    temp_dsss_list.append(basal_dsss)
                for dsss in temp_dsss_list:
                    setOfMajRefSeqs.add(dsss.referenceSequenceOf)
                master_dsss_list.append(temp_dsss_list)

        # else we are just going to be looking for the actual maj dsss
        else:
            for clade_collection_obj in self.cladeCollection_list:
                temp_dsss_list = []
                dsss_in_cc = list(
                    data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=clade_collection_obj).order_by(
                        '-abundance'))
                dsss_in_cc_in_profile = [dsss for dsss in dsss_in_cc if dsss.referenceSequenceOf in self.profile]
                basal_dsss = dsss_in_cc_in_profile[0]
                temp_dsss_list.append(basal_dsss)
                setOfMajRefSeqs.add(basal_dsss.referenceSequenceOf)
                master_dsss_list.append(temp_dsss_list)

        return master_dsss_list, setOfMajRefSeqs


    def absorb_large_init_type(self, large_init_type):
        '''The aim of this function is simply to add the infomation of the large init type to that of the small init
        type'''
        self.cladeCollection_list.extend(large_init_type.cladeCollection_list)
        self.majority_sequence_list.extend(large_init_type.majority_sequence_list)
        self.support = len(self.cladeCollection_list)
        self.set_of_maj_ref_seqs.update(large_init_type.set_of_maj_ref_seqs)

    def extract_support_from_large_initType(self, large_init_type):
        '''The aim of this function differs from above. We are extracting support for this small init_type from
        the large_init type. Once we have extracted the support then we will need to reinitialise the bigtype'''

        #1 - create the list of maj dss lists that will be added to the small init type from the large init type
        # do this by sending over any dss from the big type that is a refseq of the refseqs that the small and
        # large init types have in common
        large_init_type_ref_seqs = large_init_type.set_of_maj_ref_seqs
        small_init_type_ref_seqs = self.set_of_maj_ref_seqs
        ref_seqs_in_common = large_init_type_ref_seqs & small_init_type_ref_seqs
        temp_majdss_list_list = []
        # Keep track of whether some new maj_ref_seqs have been added to the small init type
        # TODO not sure if we need to do this
        new_maj_seq_set = set()

        for i in range(len(large_init_type.majority_sequence_list)):
            list_of_dsss_to_remove_from_large_init_type = []
            temp_dss_list = []
            for j in range(len(large_init_type.majority_sequence_list[i])):
                if large_init_type.majority_sequence_list[i][j].referenceSequenceOf in ref_seqs_in_common:
                    # Then this is one of the dsss that we should remove from the CC in the big init type and
                    # add to the small init type
                    temp_dss_list.append(large_init_type.majority_sequence_list[i][j])
                    #NB important to add the referenceSequenceOf before removing the dsss from the large type
                    new_maj_seq_set.add(large_init_type.majority_sequence_list[i][j].referenceSequenceOf)
                    list_of_dsss_to_remove_from_large_init_type.append(large_init_type.majority_sequence_list[i][j])
            for dsss in list_of_dsss_to_remove_from_large_init_type:
                large_init_type.majority_sequence_list[i].remove(dsss)

            temp_majdss_list_list.append(temp_dss_list)
        # At this point we should have a list of maj dss lists that we can extend the small init type with
        # we have also removed the dss in question from the large init type

        #2 Now extract into the small init type
        self.cladeCollection_list.extend(large_init_type.cladeCollection_list)
        self.majority_sequence_list.extend(temp_majdss_list_list)
        self.set_of_maj_ref_seqs.update(new_maj_seq_set)

        #3 Now modify the large init_type
        # we have already removed the maj seqs from the large init_type
        # need to change, profile, profile length
        # clade_collection_list should not change
        # put through check_if_initial_type_contains_basal_seqs
        # re-initialise set of set_of_maj_ref_seqs

        # new large profile should simply be the ref seqs of the small and large profiles not found in common
        # essentially we extract the small init_types' profile from the large
        large_init_type.profile = large_init_type.profile.difference(self.profile)
        large_init_type.profile_length = len(large_init_type.profile)
        large_init_type.contains_multiple_basal_sequences, \
        large_init_type.basalSequence_list = large_init_type.check_if_initialType_contains_basal_sequences()
        large_init_type.majority_sequence_list, large_init_type.set_of_maj_ref_seqs = large_init_type.create_majority_sequence_list_for_initalType_from_scratch()

    def remove_small_init_type_from_large(self, small_init_type):
        profRefSeq = list(small_init_type.profile)[0]
        self.profile = self.profile.difference(small_init_type.profile)
        self.profile_length = len(self.profile)
        self.contains_multiple_basal_sequences, \
        self.basalSequence_list = self.check_if_initialType_contains_basal_sequences()
        # remove the refSeq and dsss from the majority_sequence_list and set_of_maj_ref_seqs
        # remove refSeq from set_of_maj_ref_seqs
        if profRefSeq in self.set_of_maj_ref_seqs:
            self.set_of_maj_ref_seqs.remove(profRefSeq)
        for i in range(len(small_init_type.cladeCollection_list)): #For each list of dsss
            for j in small_init_type.cladeCollection_list[i]: # for each dsss
                if j.referenceSequenceOf == profRefSeq:
                    del small_init_type.cladeCollection_list[i][j]

        self.majority_sequence_list, self.set_of_maj_ref_seqs = self.create_majority_sequence_list_for_initalType_from_scratch()

    def __str__(self):
        return str(self.profile)

def collapsePotentialProfiles_initType_objects(footprintList, reqsupport, nprocessors):
    #TODO 10/12/17 I have just thought of something else we will have to consider. When you have a profile that
    # has a genuine codom type so there are e.g. MAJs of C3 and C3a, what will happen in this cricumstances.
    # For the time being we will not think too much about this and in the initialType profile classes that we are
    # creating, if we find a multiple basal type situation we will be sure to add the maj type to each of the
    # dsss lists.
    #TODO 10/12/17, I also think we can take the opportunity to move away from working with dataSetSampleSequences
    # and move towards working with referenceSequences, as it seems that we are always converting them to
    # reference sequences any way. Maybe for the time being we will work with dsss and get things up and running
    # before converting to working with referenceSequences.
    #TODO 09/12/17 I want to make sure that the subsets collapsing is still working and taking into account the
    # basal types compatabilites. Actually the subset collapsing is not taking into account basal sequences
    # but this is not a problem. We are essentially just consolidating the larger types into the smaller ones
    # the smaller ones still then have to found as supported using the basal sequences information.

    #TODO 09/12/17 This has become unbarably slow. We have an awful lot of redundancy in here that I would like to fix
    # especially when it comes to checking for basal types. So I am going to create some clases which will only be
    # held in RAM. ie. not in the database.
    # class will be initialType. will contain, profile, profile length, contains_multiple_basal_sequences,
    # basal_sequence_list, clade_collection_list, maj_seq_list, maj_seq_set
    '''This will be our latest battle field 07/12/17. We have a problem where we are getting multiple types
    So the object that we're playing around with here is a dictionary with key that is the set of referencesequences
    objects with a value of a 2d list with containing two lists. the first list is the clade collection object that
    this profile was found in. THe second list is the majority sequence of this set in this clade collection
    so if we want to be able to split up these objects to extract e.g. c15 and c3 sets then we'll need to create another
    dictionary entry. So we were probably relying on their only ever being one cladecollection to every set of sequences
    but if we start to do the extractions then were going to end up with multiple clade collections in the dict
    we will also need to be able to calcualate a new maj sequence for each for each of the extracted set of sequences
    we'll then have to see how this might reek havoc further down stream.
    so we end up creating the analysis types from the information that we have created here so
    I don't think that its going to be a problem that we ahve serveral types for a cgiven clade collection as it simply
    means that one clade collection is going to be found listed in multiple types. So lets give it a go!'''

    # We were having problems with the original.
    # The problem was when you have lots of long footprints that have intras in common i.e. D1 and D1a,
    # these long footprints are collapsed all the way down to D1 or D1a, depending on the maj becuase
    # there aren't any footprints that are JUST D1/D1a. So I have a better approach hopefully that will solve this
    # This approach builds on the original collapseFootprints
    # It starts with an n of the longest type and it tries to collapse each of the types into the n-1 types.
    # However, unlike the original, if it fails to collapse into the n-1 it does not then try to collapse into
    # n-2 types etc. Once we have tested all of the n length types and collapsed those that can be we move those that
    # couldn't into a list of unsupported types.
    # We then get all of the sequences found in the footprints of length n and we do itertools to get all n-1 combinations
    # we then try to fit each of these into the unsupported types. We create new footprints in order of the most supported
    # footprints and assign the associted CCs and Maj seqs. Any n length seqs that still don't have a collapse are
    # kept in the unsupported list. This is iterated dropping through the ns. I think this will be great!
    # 1 - check which supported which not, move supported to supportedlist, move unsupported to unsupported lit
    # 2 - for all unsupported check to see if can be collapsed into fooprints of n-1
    #   if can then collapse into and remove type from not supported list
    # 3 - get all n-1 combos of all seqs found in the n-1 unsupported sequences
    # All the seqs found in a CC can now be accessed from the [2] list (third list) of the foorprintList dict value
    # actually I don't think we'll need to caryy around the [2] list as we can simply use the footprint value
    # to get the list of refseqs - duh!
    # see how many footprints have these types, sort by order
    # create new footprints fot these footprints
    # if supported move to supportedlist
    # if not supported move to the non-supporetd list
    # next n...


    # Lists that will hold the footprints that have and do not have support throughout the collapse
    supportedList = []
    unsupportedList = []


    # convert the footprintList to initalTypes
    initial_types_list = []
    for fpkey, fpvalue in footprintList.items():
        initial_types_list.append(initalType(fpkey, fpvalue[0], fpvalue[1]))

    # ### DEBUG ###
    # for initT in initial_types_list:
    #     for cc in  initT.cladeCollection_list:
    #         if cc.dataSetSampleFrom.name == 'May2016-1poc-R1':
    #             print('{}'.format(initT))
    # ### DEBUG ###

    # start with n at the largest footprint in the footprintlist
    # the lowest n that we will work with is 2
    # when we get to the footprints that are two in length we will do the first type of collapsing
    # i.e. we will collapse into existant 1s if feasible but we will not do the combinations
    # section as this will be moot as each will only be one intra big anyway.
    # once we have done the inital collapse of the 2s then we will sort the 1s into
    # supported or unsupported.
    # At this point all of the footprints that were in the footprintList to begin with will be either in
    # the supported or unsupported list
    # we will then convert these types to the maj types

    # for each length starting at max and dropping by 1 with each increment
    for n in range(max([initT.profile_length for initT in initial_types_list]), 0, -1):
        # This dict will hold the top collapses where bigfootprint = key and smallFootprint = value
        collapseDict = {}
        # populate supported and unsupported list for the next n
        nList = [initT for initT in initial_types_list if initT.profile_length == n]

        # ### DEBUG ###
        for initTThis in unsupportedList:
            if initTThis.profile_length <= n:
                Seomthing = 'asdf'


        if n == 9:
            foo = 'bar'
        #
        # for initT in initial_types_list:
        #     for cc in initT.cladeCollection_list:
        #         if cc.dataSetSampleFrom.name in ['May2016-1poc-R1', 'ZB49']:
        #             debug = True
        #             print('{}'.format(initT))
        #             foo = 'bar'
        # ### DEBUG ###

        for initT in nList:

            # ### DEBUG ###
            # for cc in initT.cladeCollection_list:
            #     if cc.dataSetSampleFrom.name in ['FS1PC9-FS1PC9-N708', 'ZB49'] and cc.clade == 'C':
            #         debug = True
            #         print('{}'.format(initT))
            #         foo = 'bar'
            # ### DEBUG ###


            if initT.support >= reqsupport and not initT.contains_multiple_basal_sequences: # supported
                supportedList.append(initT)
            else: # above support threshold
                unsupportedList.append(initT)

        ############# TRY TO COLLAPSE SIZE n FOOTPRINTS INTO SIZE n-1 FOOTPRINTS ################
        # we will try iterating this as now that we have the potential to find two types in one profile, e.g.
        # a C15 and C3, we may only extract the C3 on the first one but there may still be a C15 in there.
        # whichever had the highest score will have been extracted. Eitherway, along the way, the unsupportedList
        # will have been updated
        repeat = True
        while repeat:
            collapseDict = {}
            repeat = False
            nMinusOneList = [initT for initT in initial_types_list if initT.profile_length == n-1]
            if nMinusOneList:


                for bigFootprint in unsupportedList: # For each big check if small will fit in

                    # ### DEBUG ###
                    # debug = False
                    # for cc in bigFootprint.cladeCollection_list:
                    #     if cc.dataSetSampleFrom.name == 'FS1PC9-FS1PC9-N708' and cc.clade == 'C':
                    #         debug = True
                    #         print('{}'.format(initT))
                    #         foo = 'bar'
                    # ### DEBUG ###


                    print('Assessing discovered footprint {} for supported type'.format('-'.join(str(refseq) for refseq in bigFootprint.profile)), end='\r')
                    # These three are so that we only collapse into footprints that have the maj seqs of the bigfootprints
                    # see code futher down
                    # TODO 07/12/17 this is going to cause some problems for our basal checks as we only have one maj listed
                    # maybe if we identify that footprint contains multiple basals then we can add the most abudant basal maj as well.
                    # We will be performing the check_if_contains... on each of the footprintList keys so when we do this
                    # check we and we find that there are multiple basals, then we can change the dict and add the additional majs

                    # # listOfAnalysesSampleSequencesLarge = footprintList[bigFootprint][1]
                    # listOfAnalysesSampleSequencesLarge = []
                    # for majdssslist in footprintList[bigFootprint][1]:
                    #     for dss in majdssslist:
                    #         listOfAnalysesSampleSequencesLarge.append(dss)
                    # listOfRefSeqsLarge = [smpl.referenceSequenceOf for smpl in listOfAnalysesSampleSequencesLarge]
                    # setOfMajRefSeqsLarge = set(listOfRefSeqsLarge)

                    topScore = 0
                    for smallerFootprint in nMinusOneList:
                        # TODO 08/12/17 we should only consider collapsing into the small footprint if it
                        # doesn't contain multiple basal types. so we should put it through the checker
                        # If the small foot print is a subset of the big footprint consider for collapse
                        # Only collapse if this is the best option i.e. if it give the largest number of support
                        multi_basal = smallerFootprint.contains_multiple_basal_sequences
                        if smallerFootprint.profile.issubset(bigFootprint.profile) and not multi_basal:
                            # Consider this for collapse only if the majsequences of the smaller are a subset of the maj sequences of the larger
                            # e.g. we don't want A B C D being collapsed into B C D when A is the maj of the first and B is a maj of the second

                            # simplest way to check this is to take the setOfMajRefSeqsLarge which is a set of all of the
                            # ref seqs that are majs in the cc that the footprint is found in
                            # and make sure that it is a subset of the smaller footprint in question

                            ######## MAKE SURE THAT THE MAJ SEQS OF THE n FOOTPRINT ARE IN THE n-1 FOOTPRINT
                            # TODO 07/12/17 I think we only need to find one of the large maj refs in the list of smaller maj refs
                            # I don't think we need it to be a complete subset.
                            # if setOfMajRefSeqsLarge.issubset(smallerFootprint):
                            # Just check to see if one of the maj ref seqs is found in common between the small and large footprints
                            # TODO 10/01/18 the above logic is wrong. You need to have all of the set_of_maj_ref_seqs's ref seqs
                            # to be found in the smallerFootprint. If you don't then it means you are collapsing a clade collection
                            # into the smaller footprint where its maj seq may not be found and that is not correct.
                            # This logic is not correct either. e.g. when you have a big footprint that contains multiple basal seqs
                            # e.g. set of maj ref seqs is c15h and c3 then we will never find both of these in a small footprint
                            # st this point as the small footprint would have to contain multiple basals as well
                            # What we actually need to do is ask whether if bigFootprint is not multibasal then, this is fine.
                            # TODO 10/01/18 what we actually need is quite complicated. If the big type is not multi basal,
                            # then we have no problem and we need to find all of the set of maj ref seqs in the small profile
                            # but if the large type is multi basal then it gets a little more complicated
                            # if the large type is multi basal then which of its set of maj ref seqs we need to find in the small profile
                            # is dependent on what the basal seq of the smallfootprint is.
                            # if the small has no basal seqs in it, then we need to find every sequence in the large's set of maj ref seqs
                            # that is NOT a C15x, or the C3 or C1 sequences.
                            # if small basal = C15x then we need to find every one of the large's seqs that isn't C1 or C3
                            # if small basal = C1 then we need to find every on of the large's seqs that isn't C3 or C15x
                            # if small basal = C3 then we need to find every on of the large's seqs that isn't C1 or C15x
                            # we should put this decision into a new function

                            if bigFootprint.contains_multiple_basal_sequences:
                                if does_small_footprint_contain_the_required_ref_seqs_of_the_large_footprint(bigFootprint, smallerFootprint):
                                # if len(bigFootprint.set_of_maj_ref_seqs.intersection(smallerFootprint.profile)) >= 1:
                                    # score = number of samples big was found in plus num samples small was found in
                                    score = bigFootprint.support + smallerFootprint.support
                                    if score > topScore:
                                        topScore = score
                                        repeat = True
                                        collapseDict[bigFootprint] = smallerFootprint

                                        ### DEBUG ####
                                        # if debug:
                                        #     foo = 'bar'
                            else:
                                if bigFootprint.set_of_maj_ref_seqs.issubset(smallerFootprint.profile):
                                    # if len(bigFootprint.set_of_maj_ref_seqs.intersection(smallerFootprint.profile)) >= 1:
                                    # score = number of samples big was found in plus num samples small was found in
                                    score = bigFootprint.support + smallerFootprint.support
                                    if score > topScore:
                                        topScore = score
                                        repeat = True
                                        collapseDict[bigFootprint] = smallerFootprint

                                        # ### DEBUG ####
                                        # if debug:
                                        #     foo = 'bar'

            # Once here we have tried to collapse all of the large seq into all of the smaller footprinst
            # now if there is a collapse do the collapse and then remove it from the unsupported list
            # (we will see if the smaller footprint we added it to is supported next round of n)
            # Or each key of the collapseDict add it to the value.
            # then remove the key from the unsupportedList
            # TODO 08/12/17 here is where we will need to start to implement extraction rather than just deletion for the potentially multiple
            # We will need to be careful that we don't start extracting profiles for non basal mixes, e.g. if we have
            # C3, C3a, C3b, C3d in the big and the small is C3, C3a, C3b we don't want to extract this, we want to use
            # the preivous method of deletion. I guess the best way to tell whether we want to do extraction vs deletion
            # is if there is another basal maj type in the large footprint
            # I'm starting to think that it might be worth having an initial type object where we can store useful
            # information like list of cc, list of dsss and list of basal maj
            # Yep, we now have this and I am doing the initial write now.

            collapse_dict_keys_list = list(collapseDict.keys())
            # TODO 070218 - Go through to make sure that all of the keys are found in the unsupported typelist  -they should be
            ### DEBUG
            count = 0
            for t in range(len(collapse_dict_keys_list)):
                if collapse_dict_keys_list[t] not in unsupportedList:
                    count += 1
            print('{} number of keys not found in unsupported list'.format(count))


            # The collapseDict footprint is now key=large initial type and value = small initial type

            # for footPrintToCollapse in collapseDict.keys():
            for q in range(len(collapse_dict_keys_list)):
                if n == 3:
                    count = 0
                    for t in range(len(collapse_dict_keys_list)):
                        if collapse_dict_keys_list[t] not in unsupportedList:
                            count += 1
                    print('{} number of keys not found in unsupported list'.format(count))
                # ### DEBUG ###
                # for cc in footPrintToCollapse.cladeCollection_list:
                #     if cc.dataSetSampleFrom.name == 'FS1PC9-FS1PC9-N708' and cc.clade == 'C':
                #         debug = True
                #         print('{}'.format(initT))
                #         foo = 'bar'
                # ### DEBUG ###
                # if q == 92 and n == 3:
                #     peasrs = 'asf'
                # if q == 91 and n == 3:
                #     pwer = 'asdf'
                if collapse_dict_keys_list[q] in unsupportedList:
                    appes = 'asdf'



                #### DEBUG
                # returns bool representing whether to extract. If false, delete rather than extract

                if collapse_dict_keys_list[q].support >= reqsupport and not collapse_dict_keys_list[q].contains_multiple_basal_sequences:
                    # Then this type has already had some other leftovers put into it so that it now has the required
                    # support. In this case we can remove the type from the unsupported list and continue to the next
                    unsupportedList.remove(collapse_dict_keys_list[q])
                    continue

                extraction_deletion = collapse_dict_keys_list[q].contains_multiple_basal_sequences
                if not extraction_deletion:
                    # If large type does not contain multiple basal sequences then we do not need to extract and we
                    # can collapse the large into the small.
                    # we need to simply extend the clade collection list of the large type with that of the small
                    # we also need to add the dss lists of the large to the small as well

                    small_init_type = collapseDict[collapse_dict_keys_list[q]]
                    small_init_type.absorb_large_init_type(collapse_dict_keys_list[q])

                    # # below = smaller footprint data = [[smallerCCs += bigger CCs], [smallerMajSeqsForCCs += largerMajSeqsForCCs]]
                    # footprintList[collapseDict[footPrintToCollapse]] = [footprintList[collapseDict[footPrintToCollapse]][0] + footprintList[footPrintToCollapse][0], footprintList[collapseDict[footPrintToCollapse]][1] + footprintList[footPrintToCollapse][1]]

                    # remove the large_init_type from the init_type_list and from the unsupportedList
                    initial_types_list.remove(collapse_dict_keys_list[q])
                    unsupportedList.remove(collapse_dict_keys_list[q])
                else:
                    # Here we need to implement the extraction. It is important that we extract all of the Majs that
                    # should go into the type we are collapsing into
                    # It is not just enough to send over all of the Majs that the big and small majref seqs have in common.
                    # Actually I think it is OK to just send over the Majs that the small and big have in commmon. We need to
                    # remember that we are simply listing ccts for support and the Majs found in those CCts. We can remove
                    # the subsetted sequences from the footprint of the big and put it back in the dict


                    #TODO 171217 consolidating thoughts
                    # so below we want to extend the clade collection list of the small with that of the large
                    # this is not a problem to do and very straight forwards
                    # what will cause us more of an issue is the dsss, that we extend with
                    # for the dssss we want to send over the maj dss list of each of the large CCs if the dss is of
                    # a refseq that is found in the small init types profile or
                    # should all the ref seq basal types of the large be found in the small.
                    # My thoughts below seem to think that's its ok to just add the maj dsss that are of a refseq
                    # that are found in common between the two types

                    # 1 - Collapse big to small
                    # This will all be done within the extract_support_from_large_initType method
                    small_init_type = collapseDict[collapse_dict_keys_list[q]]
                    small_init_type.extract_support_from_large_initType(collapse_dict_keys_list[q])


                    #TODO the large init type should remain in the init_type_list




                    #TODO need to check if the new profile created from the original footPrintToCollapse
                    #  that has now had the small)init_type extracted from it is already shared with another init_type
                    # If so then we need to combine the init_types.
                    # else then we don't need to do anything
                    # If the similar type is also in the collapse dict then we will collapse to that type
                    # else if the type is not in the collapse dict then we will absorb that type.
                    # There should only be a maximum of one initT that has the same
                    match = False
                    for p in range(len(initial_types_list)):
                        if initial_types_list[p].profile == collapse_dict_keys_list[q].profile and initial_types_list[p] != collapse_dict_keys_list[q]:
                            ### DEBUG ###
                            if p == 63:
                                apples = 'asdf'

                            if collapse_dict_keys_list[q] in unsupportedList:
                                appes = 'asdf'
                            ### DEBUG ###
                            # Then we have found an initT that has an exact match for the large initT in Q
                            # In this case we need to create a single initType that is a combination of the two
                            # 070218 we are going to change this and we are going to extract
                            # Let's add the existing initT to the footPinrtToCollapse
                            # This way we can make sure that the footPrintToCoolapse is still in the correct place
                            # i.e. in the unsupported list or not.

                            # If we do find a match then we need to make sure to get rid of the initT that we have
                            # absorbed into the footPrintToCollapse

                            # Should this just be the same as when a small initT absorbs a large initT?
                            # I think so but lets check
                            #TODO check to see that this is appropriate
                            # we need to check specifically if the initial_types[p] is found in the types that
                            # still need to be collapsed. so we need to slice here.
                            match = True
                            if initial_types_list[p] in collapse_dict_keys_list[q+1:]:
                                # then we need to collapse into the found type
                                # because this collapsing can cause the matching type that is also in the collapse
                                # dict to gain support we will also check to if each of the types in the collapse dict
                                # have gained sufficient support to no longer need collapsing. We will do this earlier
                                # in the process, not here.
                                initial_types_list[p].absorb_large_init_type(collapse_dict_keys_list[q])
                                unsupportedList.remove(collapse_dict_keys_list[q])
                                initial_types_list.remove(collapse_dict_keys_list[q])
                                break
                            else:
                                collapse_dict_keys_list[q].absorb_large_init_type(initial_types_list[p])


                                # The initT will also need removing from the unsupportedList if it is in there
                                # This is causing a problem
                                # Some time the initial_types_list[p] that matches the leftover collapse_dict_keys_list[q]
                                # can also be in the collapse_dict_keys_list.
                                # In this case it will be removed from the unsupportedList at this point and won't be removable
                                # later on.
                                # Shouldn't we be checking to see if initial_types_list[p] is smaller than n.
                                if initial_types_list[p] in unsupportedList:
                                    try:
                                        unsupportedList.remove(initial_types_list[p])
                                    except:
                                        pear = 'wer'
                                        # Delete the initT as this has now been absorbed into the footprint to collapse
                                initial_types_list.remove(initial_types_list[p])
                                # If the left over type is less than n then we need to now remove it from the un
                                # supported list as it will be collapsed on another iteration than this one.
                                if collapse_dict_keys_list[q].profile_length < n:
                                    try:
                                        unsupportedList.remove(collapse_dict_keys_list[q])
                                    except:
                                        apples = 'wer'
                                else:
                                    # now we need to check to see if the collapse_dict_keys_list[q] type has support bigger then
                                    # the required. If it does, then it should also be removed from the unsupported list

                                    if collapse_dict_keys_list[q].support >= reqsupport and not collapse_dict_keys_list[q].contains_multiple_basal_sequences:
                                        unsupportedList.remove(collapse_dict_keys_list[q])
                                    else:
                                        # if it doesn't have support then we simply leave it in the unsupportedlist
                                        # and it will go on to be seen if it can be collapsed into one of the insilico
                                        # types that are genearted.
                                        pass
                                break
                    if not match:
                        if collapse_dict_keys_list[q].profile_length < n:
                            try:
                                unsupportedList.remove(collapse_dict_keys_list[q])
                            except:
                                apples = 'wer'




                    #TODO the large init_type does not need removing from the initial type list.
                    # but we still need to do the check to see if the profile length is SMALLER than n. If it is smaller
                    # then n then we should remove it from the unsupported list else we should leave it in the
                    # unsupported list.



        if n > 2:

            # We only need to attempt the further collapsing if there are unsupported types to collapse
            # else move onto the next n
            if len(unsupportedList) > 1:
                # Now we are left with the footprints that are still in the unsupported list

                # ### DEBUG ###
                # if n == 8:
                #     apples = 'pwers'
                #
                # for initT in unsupportedList:
                #
                #
                #     for cc in initT.cladeCollection_list:
                #         if cc.dataSetSampleFrom.name == 'FS1PC9-FS1PC9-N708' and cc.clade == 'C':
                #             debug = True
                #             print('{}'.format(initT))
                #             foo = 'bar'
                #             ### DEBUG ###
                # ### DEBUG

                # we are going to change this.
                # Currently we are generating an enormous generator in the itertools.combinations
                # This is because we are using all of the sequences found in all of the length n types
                # However, I see non reason why we should be using combinatations of sequences found in different profiles.
                # This is not in keeping with the biology.
                # Instead the synthetic types should be made from
                # combinattions of sequences found in each profile
                # This will better represent the biology and also greatly cutdown on computational cost
                # And it will mean that we can multiprocess this
                setOfRefSeqsOfSmpSeqsFoundInCCsOfLengthN = set()
                listOfTypesOfLengthN = [initT for initT in initial_types_list if initT.profile_length == n]
                # Only carry on if we have lengthN footprints to get sequences from
                if listOfTypesOfLengthN:

                    # TODO here is where we should start the new approach.
                    # add the nLength Types into a MP list
                    # then add the stops
                    # then in the worker create a dict that is same as the one below, only the using the nLengthtype instead
                    # Create the queues that will hold the sample information
                    print('\rstarting to generate Synthetic footprints for collapse', end='')

                    taskQueue = Queue()
                    collapseNMerDictManager = Manager()
                    collapseNMerDict = collapseNMerDictManager.dict()
                    #outputQueue = Queue()

                    for nLengthType in listOfTypesOfLengthN:
                        taskQueue.put(nLengthType.profile)

                    numProcessors = nprocessors
                    for N in range(numProcessors):
                        taskQueue.put('STOP')

                    allProcesses = []

                    for N in range(numProcessors):
                        p = Process(target=workerDiscoveryOne, args=(taskQueue, collapseNMerDict, n))
                        allProcesses.append(p)
                        p.start()

                    for p in allProcesses:
                        p.join()
                    # The manager(dict) object was not behaving correctly when I was trying to append items to the list
                    # values. However, when converting to a dict, it does.
                    collapseNMerDict = dict(collapseNMerDict)

                    # outputQueue.put('STOP')

                    # collapseNMerDict = {}
                    # for dictionary in iter(outputQueue.get, 'STOP'):
                    #     collapseNMerDict.update(dictionary)

                    # for nLengthType in listOfTypesOfLengthN:
                    #     # nLengthType should be a frozen set so we will take the refSeqs it contains and add them to the setOfRefSeqs...
                    #     setOfRefSeqsOfSmpSeqsFoundInCCsOfLengthN.update([refSeq for refSeq in nLengthType])

                    #Now go through each of the (n-1)mer footprints and see if they fit into a footprint in the unsuported list
                    #This dict value = the synthetic footprint of length n-1, value, a list of unsupported types that are mastersets of the synthetic footprint

                    # collapseNMerDict = {frozenset(tup): [] for tup in itertools.combinations(setOfRefSeqsOfSmpSeqsFoundInCCsOfLengthN, n-1)}
                    print('\rGenerated {} Synthetic footprints'.format(len(collapseNMerDict)), end='')
                    # We create a dict that is key, type and value is the set of the majrefseqs found in the type
                    # this was being calculated once for every synthetic type. But the number of synthetic types are getting very large
                    # so by having this dict we should speed things up considerably
                    # This info is used to make sure that when collapsing types the maj sequences are shared between the
                    # larger and smaller types

                    #TODO 231217 we should no longer need this
                    # majRefSeqsDict = {}
                    # for nLengthType in unsupportedList:
                    #     listOfAnalysesSampleSequencesLarge = []
                    #     for majdssslist in footprintList[nLengthType][1]:
                    #         for dss in majdssslist:
                    #             listOfAnalysesSampleSequencesLarge.append(dss)
                    #     listOfMajRefSeqsLarge = [smpl.referenceSequenceOf for smpl in listOfAnalysesSampleSequencesLarge]
                    #     setOfMajRefSeqsLarge = set(listOfMajRefSeqsLarge)
                    #     majRefSeqsDict[nLengthType] = setOfMajRefSeqsLarge

                    ## Items for printing out our progress
                    total = len(collapseNMerDict)*len(unsupportedList)
                    count = 0
                    printVal = 0.01
                    print('\rChecking new set of synthetic types', end='')


                    for frTupKey in collapseNMerDict.keys(): # for each of the synthetic footprints
                        # TODO 10/01/18 we need to check that each of the frTupKeys contains multiple basal seqs
                        # If it does then we can't collapse into it.
                        if does_set_of_ref_seqs_contain_multiple_basal_types(frTupKey):
                            continue
                        for nLengthType in unsupportedList: # for each of the unsupported init_types
                            ## Items for printing out progress
                            count += 1
                            percentComplete = count/total
                            if percentComplete > printVal:
                                print("\r%.2f" % percentComplete, end='')
                                printVal = max([printVal + 0.01, percentComplete + 0.01])

                        #For each of the N-1mers see if they fit within the unsupported types
                        # if so then add the footprint into the n-1mers list
                            # we should check for the maj rule again. I.e. we should want all of the Majs of the footprint in Q to be in the kmer footprint
                            # maybe recycle code for this
                            if frTupKey.issubset(nLengthType.profile):
                                # Now check to see that at least one of the set_of_maj_ref_seqs seqs is found in the
                                # frTupKey.
                                if len(nLengthType.set_of_maj_ref_seqs & frTupKey) > 0:
                                    # Then this nLengthType init_type can be collapsed into the frTupKey
                                    collapseNMerDict[frTupKey].append(nLengthType)
                    # Here we have a populated dict.
                    # Order dict
                    # We are only interseted in kmers that were found in the unsupported types
                    # Because alot of the kmers will be found in the sequences they originated from
                    # require the kmer to be associated with more than 1 cladecollection
                    # I am not able to express my logic here perfectly but I am sure that this makes sense
                    ## DEBUG TO INSPECT collapseNMERDICT ##
                    apples = dict(collapseNMerDict)

                    listOfnKmersWithCollapsableFootprints = [kmer for kmer in collapseNMerDict.keys() if len(collapseNMerDict[kmer]) > 1]
                    # No need to continue if there are no footprints that match the nKmers
                    if listOfnKmersWithCollapsableFootprints:
                        # We now go through the kmers according to the number of footprints they were found in
                        orderedListOfPopulatedKmers = sorted(listOfnKmersWithCollapsableFootprints, key=lambda x: len(collapseNMerDict[x]), reverse=True)

                        for kmer in orderedListOfPopulatedKmers:
                            #for each initType in the kmer lists

                        #     ###DEBUG
                        #     for bigFootprintToCollapse in collapseNMerDict[kmer]:
                        #         for cc in bigFootprintToCollapse.cladeCollection_list:
                        #             if cc.dataSetSampleFrom.name == 'FS1PC9-FS1PC9-N708' and cc.clade == 'C' and n == 4:
                        #                 apples = 'asdf'
                        #
                        #
                        #
                        #
                        #
                        # for kmer in orderedListOfPopulatedKmers:
                        #     #for each initType in the kmer lists
                        #     for initT_check in initial_types_list:
                        #         for cc in initT_check.cladeCollection_list:
                        #             if cc.dataSetSampleFrom.name in ['FS1PC9-FS1PC9-N708', 'ZB49'] and cc.clade == 'C' and len(initT_check.profile) > 2:
                        #                 debug = True
                        #                 print('profile = {}, KMER = {}'.format(initT_check.profile, kmer))
                        #
                        #     for bigFootprintToCollapse in collapseNMerDict[kmer]:
                        #         for cc in bigFootprintToCollapse.cladeCollection_list:
                        #             if cc.dataSetSampleFrom.name == 'FS1PC9-FS1PC9-N708' and cc.clade == 'C' and n == 4:
                        #                 apples = 'asdf'
                        #
                        #     ##### DEBUG
                            for k in range(len(collapseNMerDict[kmer])):
                                # chek to see if the footprint is still in the unsupportedlist #TODO!! 11/01/18 also need to check that the bigFootprint hasn't already been collapsed
                                # if it isn't then it has already been associated to a new footprint
                                #   pass it over
                                if collapseNMerDict[kmer][k] in unsupportedList and kmer.issubset(collapseNMerDict[kmer][k].profile): # Then this footprint hasn't been collapsed anywhere yet. Here we also check to make sure that the profile hasn't been changed
                                    # if an init_type already exists with the profile (kmer) in question then add the big init_type to it
                                    #TODO 231217 we will have to check whether the big init_type is a multiple basal seqs
                                    # and therefore whether this is an extraction or a absorption
                                    exists = False
                                    for i in range(len(initial_types_list)):
                                    # for initT_one in initial_types_list:
                                        if initial_types_list[i].profile == kmer:
                                            exists = True
                                            # TODO 231217 we will have to check whether the big init_type is a multiple basal seqs
                                            # and therefore whether this is an extraction or a absorption
                                            # bear in mind that it doesn't matter if the initT we are absorbing or
                                            # extracting into is a multi basal. We will wory about that in the next
                                            # iteration
                                            if collapseNMerDict[kmer][k].contains_multiple_basal_sequences:
                                                # TODO 231317 this still needs writing
                                                # Then we need to extract
                                                initial_types_list[i].extract_support_from_large_initType(collapseNMerDict[kmer][k])

                                                # Once we have extracted into the smaller type
                                                # check to see if the big init Type still contains set_of_maj_ref_seqs
                                                if collapseNMerDict[kmer][k].set_of_maj_ref_seqs:
                                                    #TODO!!!! 10/01/18 we first need to check if the new profile already exists and if it does we need to do as we do above and add it to the similar one
                                                    # if the big init type still exists with maj containing profiles
                                                    # then remove from unsupported list
                                                    for j in range(len(initial_types_list)):
                                                    # for initT_two in initial_types_list:
                                                        #DEBUG
                                                        if j == 22 and k ==2 and i == 553:
                                                            appelse = 'pears'
                                                            if collapseNMerDict[kmer][k] in unsupportedList:
                                                                apples = 'peras'
                                                        wewerehere = False
                                                        if initial_types_list[j].profile == collapseNMerDict[kmer][k].profile and initial_types_list[j] != collapseNMerDict[kmer][k]:
                                                            wewerehere = True
                                                            # Then we have found an initT that has an exact match for the large initT in Q
                                                            # In this case we need to create a single initType that is a combination of the two
                                                            # Let's add the existing initT to the footPinrtToCollapse
                                                            # This way we can make sure that the footPrintToCoolapse is still in the correct place
                                                            # i.e. in the unsupported list or not.

                                                            # If we do find a match then we need to make sure to get rid of the initT that we have
                                                            # absorbed into the footPrintToCollapse

                                                            # Should this just be the same as when a small initT absorbs a large initT?
                                                            # I think so but lets check
                                                            # TODO check to see that this is appropriate
                                                            collapseNMerDict[kmer][k].absorb_large_init_type(initial_types_list[j])

                                                            # the initT will need removing from the inital_types_list and the unsupportedList as it no longer exists.
                                                            if initial_types_list[j] in unsupportedList:
                                                                unsupportedList.remove(initial_types_list[j])
                                                            initial_types_list.remove(initial_types_list[j])
                                                            break

                                                    # We now need to decide if the footprint to collapse should be removed from the unsupportedList.
                                                    # This will depend on if it is longer then n or not.
                                                    if collapseNMerDict[kmer][k].profile_length < n:
                                                        ### DEBUG ###
                                                        try:
                                                            unsupportedList.remove(collapseNMerDict[kmer][k])
                                                        except:
                                                            aples = wewerehere
                                                    # unsupportedList.remove(bigFootprintToCollapse)

                                                else:
                                                    # then the big init_type no longer contains any of its original maj ref seqs and
                                                    # so should be delted.
                                                    initial_types_list.remove(collapseNMerDict[kmer][k])
                                                    unsupportedList.remove(collapseNMerDict[kmer][k])
                                            else:
                                                # TODO 231217 then we need can absorb
                                                initial_types_list[i].absorb_large_init_type(collapseNMerDict[kmer][k])
                                                # make sure that we then get rid of the bigFootprintToCollapse
                                                initial_types_list.remove(collapseNMerDict[kmer][k])
                                                unsupportedList.remove(collapseNMerDict[kmer][k])
                                            break
                                    if not exists: # then the kmer was not already represented by an existing initT
                                        #TODO 270118 We need to check to see if the big type is contains mutiple
                                        # basal if it does then we should extract as above. This should be exactly the
                                        # same code as above but extracting into a new type rather than an existing one
                                        #Try to create a blank type
                                        ########## NEW ###########
                                        if collapseNMerDict[kmer][k].contains_multiple_basal_sequences:
                                            new_blank_initial_type = initalType(refSeq_set=kmer, cladeCollection_list=list(collapseNMerDict[kmer][k].cladeCollection_list))
                                            initial_types_list.append(new_blank_initial_type)

                                            # now remove the above new type's worth of info from the current big footprint
                                            collapseNMerDict[kmer][k].substract_init_type_from_other_init_type(new_blank_initial_type)


                                            if collapseNMerDict[kmer][k].set_of_maj_ref_seqs:
                                                # we first need to check if the new profile already exists and if it does we need to do as we do above and add it to the similar one
                                                # if the big init type is still exists with maj containing profiles
                                                # then remove from unsupported list
                                                for initT in initial_types_list:
                                                    if initT.profile == collapseNMerDict[kmer][k].profile and initT != collapseNMerDict[kmer][k]:
                                                        # Then we have found an initT that has an exact match for the large initT in Q
                                                        # In this case we need to create a single initType that is a combination of the two
                                                        # Let's add the existing initT to the footPinrtToCollapse
                                                        # This way we can make sure that the footPrintToCoolapse is still in the correct place
                                                        # i.e. in the unsupported list or not.

                                                        # If we do find a match then we need to make sure to get rid of the initT that we have
                                                        # absorbed into the footPrintToCollapse

                                                        # Should this just be the same as when a small initT absorbs a large initT?
                                                        # I think so but lets check
                                                        # TODO check to see that this is appropriate
                                                        collapseNMerDict[kmer][k].absorb_large_init_type(initT)

                                                        # the initT will need removing from the inital_types_list and the unsupportedList as it no longer exists.
                                                        initial_types_list.remove(initT)
                                                        if initT in unsupportedList:
                                                            unsupportedList.remove(initT)
                                                        break

                                                # We now need to decide if the footprint to collapse should be removed from the unsupportedList.
                                                # This will depend on if it is longer then n or not.
                                                if collapseNMerDict[kmer][k].profile_length < n:
                                                    unsupportedList.remove(collapseNMerDict[kmer][k])
                                                    # unsupportedList.remove(bigFootprintToCollapse)

                                            else:
                                                # then the big init_type no longer contains any of its original maj ref seqs and
                                                # so should be delted.
                                                initial_types_list.remove(collapseNMerDict[kmer][k])
                                                unsupportedList.remove(collapseNMerDict[kmer][k])
                                        ############## NEW #################
                                        # TODO 231217 in this case we should create a new inital_type
                                        # I am going to have to think about how we will do this
                                        # this still needs to be coded

                                        # the profile will be the kmer
                                        # the cladeCollectionList will be the same as the large init type
                                        # the basalSequence_list can come from the self.check_if_initialType_contains_basal_sequences()

                                        # NB i have made it so that if initalType() doesn't get a majdsss list then it will create one from scratch
                                        ######## ALSO CHANGED ########
                                        else:
                                            new_initial_type = initalType(refSeq_set=kmer, cladeCollection_list=collapseNMerDict[kmer][k].cladeCollection_list)
                                            initial_types_list.append(new_initial_type)
                                        #
                                            # now delete the big init_type from the intial types list and from the unsupportedList
                                            initial_types_list.remove(collapseNMerDict[kmer][k])
                                            if collapseNMerDict[kmer][k] in unsupportedList:
                                                unsupportedList.remove(collapseNMerDict[kmer][k])
                                    # if kmer in footprintList.keys():
                                    #     footprintList[kmer] = [footprintList[kmer][0] + footprintList[bigFootprintToCollapse][0],
                                    #                            footprintList[kmer][1] + footprintList[bigFootprintToCollapse][1]]
                                    #     # then remove the footprint from the unsupported type
                                    #     # Also remove the bigFootprintprintToCollapse from the footprintList as its details are not in another key
                                    #     unsupportedList.remove(bigFootprintToCollapse)
                                    #     del footprintList[bigFootprintToCollapse]
                                    # else:
                                    #     # if the footprint is not already in the footprint list then create newfootprint item and add it to the footprintlist and assign the
                                    #     # CC details to it.
                                    #     # remove the footprint from the unsupported type list
                                    #     # Also remove the bigFootprintprintToCollapse from the footprintList as its details are not in another key
                                    #     footprintList[kmer] = [footprintList[bigFootprintToCollapse][0], footprintList[bigFootprintToCollapse][1]]
                                    #     unsupportedList.remove(bigFootprintToCollapse)
                                    #     del footprintList[bigFootprintToCollapse]
                            # each unsupportedfootprint in the nmer has been collapsed
                        # each of the kmers have been gone through and all unsupportedfootprints that could fit into the kmers have been collapsed
                        # Time for the next n here!
                        # Check to see what is happening once we get to the smaller ns
                        # It may be that we can use the code from below to simply put the maj as the type for footprints that are still
                        # unassigned
        else:
            # At this point we may have some 2s or larger in the unsupported list still
            # we will collapse these to their maj ref seq type
            # Once we have pushed these into the 1s or made new 1s if the maj ref seq type did not already exist
            # then we should collect all of the 1 length footprints and put them in the supported list
            # Or actually in theory so long as we have deleted collapsed footprints from the footprintList we should
            # be able to return the footprint list
            while unsupportedList:
                unsupFp = unsupportedList[0]
                
                # For each clade collection
                #TODO 241217 we need to do each cladeCollection individually here bear in mind.
                for i in range(len(unsupFp.cladeCollection_list)):
                    # check each maj dsss
                    for maj_dss in unsupFp.majority_sequence_list[i]:
                        # Check to see if an initial type with profile of that maj_dss refseq already exists
                        #TODO 241217 it may be worth having a dictionary of the profiles and init types to speed this up.
                        # Let's see how slow it is.
                        found = False
                        for initT in [init for init in initial_types_list if init.profile_length == 1]:
                            if maj_dss.referenceSequenceOf in initT.profile:
                                # Then we have found a profile that has the refseq as a profile
                                # If the big init type is multi basal then we need to extract
                                # if the clade maj dss contains multiple sequences
                                # if unsupFp.contains_multiple_basal_sequences:
                                #     #TODO 241217 check that the large type has had the basal seqs readjusted
                                #     # so that when we do the next maj_dss it may be non multiple basal
                                #     initT.extract_support_from_large_initType(unsupFp)
                                # # else we need to absorb
                                # else:
                                #     initT.absorb_large_init_type(unsupFp)

                                # If found then we simply need to add a clade collection to the initT support
                                # and the appropriate maj dsss
                                initT.cladeCollection_list.append(unsupFp.cladeCollection_list[i])
                                initT.majority_sequence_list.append([maj_dss])

                                found = True
                                break
                        if not found:
                            # If the type init_type doesn't already exist we must create one.
                            new_initial_type = initalType(refSeq_set=frozenset([maj_dss.referenceSequenceOf]), cladeCollection_list=[unsupFp.cladeCollection_list[i]])
                            initial_types_list.append(new_initial_type)
                            # We must then alter the big_init_type
                            # Actually I don't think we do need to alter the big type
                            # We just put it into the small for each of the cladeCcols and maj dsss
                            # once this is done then we delte the big type

                            # unsupFp.remove_small_init_type_from_large(new_initial_type)

                            # If the large init type no longer contains any maj ref seqs then we can delete it
                # Here we have completed collapsing one big_init_type and we can now get rid of the type
                unsupportedList.remove(unsupFp)
                initial_types_list.remove(unsupFp)
                    
            # All unsupported footprints have been associated to their maj and deleted from the footprintList
            return initial_types_list

    return False

def does_set_of_ref_seqs_contain_multiple_basal_types(frozen_set_of_ref_seqs):
    basal_count = 0
    C15_found = False
    for ref_seq in frozen_set_of_ref_seqs:
        if 'C15' in ref_seq.name and not C15_found:
            basal_count += 1
            C15_found = True
            continue
        elif ref_seq.name == 'C3':
            basal_count += 1
            continue
        elif ref_seq.name == 'C1':
            basal_count += 1
            continue
    if basal_count > 1:
        return True
    else:
        return False


def does_small_footprint_contain_the_required_ref_seqs_of_the_large_footprint(big_init_type, small_init_type):

    if small_init_type.basalSequence_list:
        if 'C15' in small_init_type.basalSequence_list[0]:
            # then this is a C15x basal type and we will need to find all sequences that are not C1 or C3
            set_of_seqs_to_find = set()
            ref_seqs_in_big_init_type = list(big_init_type.set_of_maj_ref_seqs)

            for ref_seq in ref_seqs_in_big_init_type:
                if ref_seq.name in ['C1', 'C3']:
                    # then this is a squence we don't need to find
                    pass
                else:
                    set_of_seqs_to_find.add(ref_seq)
            # Here we have the list of the ref_seqs that we need to find in the small_init_type.profile
            if set_of_seqs_to_find.issubset(small_init_type.profile):
                return True
            else:
                return False
        elif small_init_type.basalSequence_list[0] == 'C1':
            # then this is a C1 basal type and we need to find all sequence that are not C15x or C3
            set_of_seqs_to_find = set()
            ref_seqs_in_big_init_type = list(big_init_type.set_of_maj_ref_seqs)

            for ref_seq in ref_seqs_in_big_init_type:
                if 'C15' in ref_seq.name or ref_seq.name == 'C3':
                    # then this is a squence we don't need to find
                    pass
                else:
                    set_of_seqs_to_find.add(ref_seq)
            # Here we have the list of the ref_seqs that we need to find in the small_init_type.profile
            if set_of_seqs_to_find.issubset(small_init_type.profile):
                return True
            else:
                return False
        elif small_init_type.basalSequence_list[0] == 'C3':
            # then this is a C3 basal type and we need to find all sequence that are not C15x or C1
            set_of_seqs_to_find = set()
            ref_seqs_in_big_init_type = list(big_init_type.set_of_maj_ref_seqs)

            for ref_seq in ref_seqs_in_big_init_type:
                if 'C15' in ref_seq.name or ref_seq.name == 'C1':
                    # then this is a squence we don't need to find
                    pass
                else:
                    set_of_seqs_to_find.add(ref_seq)
            # Here we have the list of the ref_seqs that we need to find in the small_init_type.profile
            if set_of_seqs_to_find.issubset(small_init_type.profile):
                return True
            else:
                return False
    else:
        # if the small_init_type doesn't contain a basal sequence sequence, then we need to find all of the seqs
        # in the big_intit_type.set_of_maj_ref_seqs that are not C15x, C1 or C3
        set_of_seqs_to_find = set()
        ref_seqs_in_big_init_type = list(big_init_type.set_of_maj_ref_seqs)

        for ref_seq in ref_seqs_in_big_init_type:
            if 'C15' in ref_seq.name or ref_seq.name in ['C1', 'C3']:
                # then this is a squence we don't need to find
                pass
            else:
                set_of_seqs_to_find.add(ref_seq)
        # Here we have the list of the ref_seqs that we need to find in the small_init_type.profile
        if set_of_seqs_to_find.issubset(small_init_type.profile):
            return True
        else:
            return False

def collapsePotentialProfiles(footprintList, reqsupport, nprocessors):
    #TODO 10/12/17 I have just thought of something else we will have to consider. When you have a profile that
    # has a genuine codom type so there are e.g. MAJs of C3 and C3a, what will happen in this cricumstances.
    # For the time being we will not think too much about this and in the initialType profile classes that we are
    # creating, if we find a multiple basal type situation we will be sure to add the maj type to each of the
    # dsss lists.
    #TODO 10/12/17, I also think we can take the opportunity to move away from working with dataSetSampleSequences
    # and move towards working with referenceSequences, as it seems that we are always converting them to
    # reference sequences any way. Maybe for the time being we will work with dsss and get things up and running
    # before converting to working with referenceSequences.
    #TODO 09/12/17 I want to make sure that the subsets collapsing is still working and taking into account the
    # basal types compatabilites. Actually the subset collapsing is not taking into account basal sequences
    # but this is not a problem. We are essentially just consolidating the larger types into the smaller ones
    # the smaller ones still then have to found as supported using the basal sequences information.

    #TODO 09/12/17 This has become unbarably slow. We have an awful lot of redundancy in here that I would like to fix
    # especially when it comes to checking for basal types. So I am going to create some clases which will only be
    # held in RAM. ie. not in the database.
    # class will be initialType. will contain, profile, profile length, contains_multiple_basal_sequences,
    # basal_sequence_list, clade_collection_list, maj_seq_list, maj_seq_set
    '''This will be our latest battle field 07/12/17. We have a problem where we are getting multiple types
    So the object that we're playing around with here is a dictionary with key that is the set of referencesequences
    objects with a value of a 2d list with containing two lists. the first list is the clade collection object that
    this profile was found in. THe second list is the majority sequence of this set in this clade collection
    so if we want to be able to split up these objects to extract e.g. c15 and c3 sets then we'll need to create another
    dictionary entry. So we were probably relying on their only ever being one cladecollection to every set of sequences
    but if we start to do the extractions then were going to end up with multiple clade collections in the dict
    we will also need to be able to calcualate a new maj sequence for each for each of the extracted set of sequences
    we'll then have to see how this might reek havoc further down stream.
    so we end up creating the analysis types from the information that we have created here so
    I don't think that its going to be a problem that we ahve serveral types for a cgiven clade collection as it simply
    means that one clade collection is going to be found listed in multiple types. So lets give it a go!'''

    # We were having problems with the original.
    # The problem was when you have lots of long footprints that have intras in common i.e. D1 and D1a,
    # these long footprints are collapsed all the way down to D1 or D1a, depending on the maj becuase
    # there aren't any footprints that are JUST D1/D1a. So I have a better approach hopefully that will solve this
    # This approach builds on the original collapseFootprints
    # It starts with an n of the longest type and it tries to collapse each of the types into the n-1 types.
    # However, unlike the original, if it fails to collapse into the n-1 it does not then try to collapse into
    # n-2 types etc. Once we have tested all of the n length types and collapsed those that can be we move those that
    # couldn't into a list of unsupported types.
    # We then get all of the sequences found in the footprints of length n and we do itertools to get all n-1 combinations
    # we then try to fit each of these into the unsupported types. We create new footprints in order of the most supported
    # footprints and assign the associted CCs and Maj seqs. Any n length seqs that still don't have a collapse are
    # kept in the unsupported list. This is iterated dropping through the ns. I think this will be great!
    # 1 - check which supported which not, move supported to supportedlist, move unsupported to unsupported lit
    # 2 - for all unsupported check to see if can be collapsed into fooprints of n-1
    #   if can then collapse into and remove type from not supported list
    # 3 - get all n-1 combos of all seqs found in the n-1 unsupported sequences
    # All the seqs found in a CC can now be accessed from the [2] list (third list) of the foorprintList dict value
    # actually I don't think we'll need to caryy around the [2] list as we can simply use the footprint value
    # to get the list of refseqs - duh!
    # see how many footprints have these types, sort by order
    # create new footprints fot these footprints
    # if supported move to supportedlist
    # if not supported move to the non-supporetd list
    # next n...


    # Lists that will hold the footprints that have and do not have support throughout the collapse
    supportedList = []
    unsupportedList = []




    # start with n at the largest footprint in the footprintlist
    # the lowest n that we will work with is 2
    # when we get to the footprints that are two in length we will do the first type of collapsing
    # i.e. we will collapse into existant 1s if feasible but we will not do the combinations
    # section as this will be moot as each will only be one intra big anyway.
    # once we have done the inital collapse of the 2s then we will sort the 1s into
    # supported or unsupported.
    # At this point all of the footprints that were in the footprintList to begin with will be either in
    # the supported or unsupported list
    # we will then convert these types to the maj types

    # for each length starting at max and dropping by 1 with each increment
    for n in range(max([len(keyName) for keyName in footprintList.keys()]), 0, -1):
        # This dict will hold the top collapses where bigfootprint = key and smallFootprint = value
        collapseDict = {}
        # populate supported and unsupported list for the next n
        nList = [keyName for keyName in footprintList.keys() if len(keyName) == n]
        for keyName in nList:

            # we have to ensure that the check_if_contains... runs so we will assign its return to a variable
            # it may not run if we have it as an additional conditional in the if as the first part of the if will
            # be evaluated first
            check = check_if_contains_multiple_basal_seqs(footprintList, keyName)
            if len(footprintList[keyName][0]) >= reqsupport and check: # supported
                supportedList.append(keyName)
            else: # above support threshold
                unsupportedList.append(keyName)

        ############# TRY TO COLLAPSE SIZE n FOOTPRINTS INTO SIZE n-1 FOOTPRINTS ################
        # we will try iterating this as now that we have the potential to find two types in one profile, e.g.
        # a C15 and C3, we may only extract the C3 on the first one but there may still be a C15 in there.
        # whichever had the highest score will have been extracted. Eitherway, along the way, the unsupportedList
        # will have been updated
        repeat = True
        while repeat:
            collapseDict = {}
            repeat = False
            nMinusOneList = [keyName for keyName in footprintList.keys() if len(keyName) == n-1]
            if nMinusOneList:


                for bigFootprint in unsupportedList: # For each big check if small will fit in
                    print('\nAssessing discovered footprint {} for supported type'.format('-'.join(str(refseq) for refseq in bigFootprint)))
                    # These three are so that we only collapse into footprints that have the maj seqs of the bigfootprints
                    # see code futher down
                    # TODO 07/12/17 this is going to cause some problems for our basal checks as we only have one maj listed
                    # maybe if we identify that footprint contains multiple basals then we can add the most abudant basal maj as well.
                    # We will be performing the check_if_contains... on each of the footprintList keys so when we do this
                    # check we and we find that there are multiple basals, then we can change the dict and add the additional majs

                    # listOfAnalysesSampleSequencesLarge = footprintList[bigFootprint][1]
                    listOfAnalysesSampleSequencesLarge = []
                    for majdssslist in footprintList[bigFootprint][1]:
                        for dss in majdssslist:
                            listOfAnalysesSampleSequencesLarge.append(dss)
                    listOfRefSeqsLarge = [smpl.referenceSequenceOf for smpl in listOfAnalysesSampleSequencesLarge]
                    setOfMajRefSeqsLarge = set(listOfRefSeqsLarge)
                    score = 0
                    topScore = 0
                    for smallerFootprint in nMinusOneList:
                        # TODO 08/12/17 we should only consider collapsing into the small footprint if it
                        # doesn't contain multiple basal types. so we should put it through the checker
                        # If the small foot print is a subset of the big footprint consider for collapse
                        # Only collapse if this is the best option i.e. if it give the largest number of support
                        checked = check_if_contains_multiple_basal_seqs(footprintList, smallerFootprint)
                        if set(smallerFootprint).issubset(set(bigFootprint)) and checked:
                            # Consider this for collapse only if the majsequences of the smaller are a subset of the maj sequences of the larger
                            # e.g. we don't want A B C D being collapsed into B C D when A is the maj of the first and B is a maj of the second

                            # simplest way to check this is to take the setOfMajRefSeqsLarge which is a set of all of the
                            # ref seqs that are majs in the cc that the footprint is found in
                            # and make sure that it is a subset of the smaller footprint in question

                            ######## MAKE SURE THAT THE MAJ SEQS OF THE n FOOTPRINT ARE IN THE n-1 FOOTPRINT
                            # TODO 07/12/17 I think we only need to find one of the large maj refs in the list of smaller maj refs
                            # I don't think we need it to be a complete subset.
                            # if setOfMajRefSeqsLarge.issubset(smallerFootprint):
                            # Just check to see if one of the
                            if len(setOfMajRefSeqsLarge.intersection(smallerFootprint)) >= 1:
                                # score = number of samples big was found in plus num samples small was found in
                                score = len(footprintList[bigFootprint][0]) + len(footprintList[smallerFootprint][0])
                                if score > topScore:
                                    topScore = score
                                    repeat = True
                                    collapseDict[bigFootprint] = smallerFootprint


            # Once here we have tried to collapse all of the large seq into all of the smaller footprinst
            # now if there is a collapse do the collapse and then remove it from the unsupported list
            # (we will see if the smaller footprint we added it to is supported next round of n)
            # Or each key of the collapseDict add it to the value.
            # then remove the key from the unsupportedList
            # TODO 08/12/17 here is where we will need to start to implement extraction rather than just deletion for the potentially multiple
            # We will need to be careful that we don't start extracting profiles for non basal mixes, e.g. if we have
            # C3, C3a, C3b, C3d in the big and the small is C3, C3a, C3b we don't want to extract this, we want to use
            # the preivous method of deletion. I guess the best way to tell whether we want to do extraction vs deletion
            # is if there is another basal maj type in the large footprint
            # I'm starting to think that it might be worth having an initial type object where we can store useful
            # information like list of cc, list of dsss and list of basal maj
            for footPrintToCollapse in collapseDict.keys():
                # returns bool representing whether to extract. If false, delete rather than extract
                extraction_deletion = check_if_contains_multiple_basal_majs(footPrintToCollapse)
                if not extraction_deletion:
                    # below = smaller footprint data = [[smallerCCs += bigger CCs], [smallerMajSeqsForCCs += largerMajSeqsForCCs]]
                    footprintList[collapseDict[footPrintToCollapse]] = [footprintList[collapseDict[footPrintToCollapse]][0] + footprintList[footPrintToCollapse][0], footprintList[collapseDict[footPrintToCollapse]][1] + footprintList[footPrintToCollapse][1]]
                    del footprintList[footPrintToCollapse]
                    unsupportedList.remove(footPrintToCollapse)
                else:
                    # Here we need to implement the extraction. It is important that we extract all of the Majs that
                    # should go into the type we are collapsing into
                    # It is not just enough to send over all of the Majs that the big and small majref seqs have in common.
                    # Actually I think it is OK to just send over the Majs that the small and big have in commmon. We need to
                    # remember that we are simply listing ccts for support and the Majs found in those CCts. We can remove
                    # the subsetted sequences from the footprint of the big and put it back in the dict


                    # 1 - Collapse big to small
                    # Here we have to make a new list of lists for the majSeqs for each cct with only the maj seqs
                    # of the majs found in the small type we are collapsing into

                    # Get list of MajRefSeqs for big and small type
                    listOfAnalysesSampleSequencesLarge = []
                    for majdssslist in footprintList[footPrintToCollapse][1]:
                        for dss in majdssslist:
                            listOfAnalysesSampleSequencesLarge.append(dss)
                    listOfRefSeqsLarge = [smpl.referenceSequenceOf for smpl in listOfAnalysesSampleSequencesLarge]
                    setOfMajRefSeqsLarge = set(listOfRefSeqsLarge)

                    listOfAnalysesSampleSequencesSmall = []
                    for majdssslist in footprintList[collapseDict[footPrintToCollapse]][1]:
                        for dss in majdssslist:
                            listOfAnalysesSampleSequencesSmall.append(dss)
                    listOfRefSeqsSmall = [smpl.referenceSequenceOf for smpl in listOfAnalysesSampleSequencesSmall]
                    setOfMajRefSeqsSmall = set(listOfRefSeqsSmall)

                    listOfMajSeqListsBigFootprint = footprintList[footPrintToCollapse][1]
                    extractedMajSeqs = []

                    refSeqsInCommon = setOfMajRefSeqsLarge.intersection(setOfMajRefSeqsSmall)
                    for i in range(len(listOfMajSeqListsBigFootprint)):
                        temp_listEx = []
                        dsssToRemove = []
                        for j in range(len(listOfMajSeqListsBigFootprint[i])):
                            if listOfMajSeqListsBigFootprint[i][j].referenceSequenceOf in refSeqsInCommon:
                                # Then this is one of the seqs we want
                                temp_listEx.append(listOfMajSeqListsBigFootprint[i][j])
                                # we then need to remove this majSeq from the big footprints list
                                dsssToRemove.append(listOfMajSeqListsBigFootprint[i][j])
                                # footprintList[footPrintToCollapse][1][i].remove(listOfMajSeqListsBigFootprint[i][j])
                        for dsstr in dsssToRemove:
                            footprintList[footPrintToCollapse][1][i].remove(dsstr)
                        # add the list of Maj dss to the list of lits
                        extractedMajSeqs.append(temp_listEx)

                    # Now perform the collapse using the extractedMajSeqs list as the 2D list of maj dsss
                    # below = smaller footprint data = [[smallerCCs += bigger CCs], [smallerMajSeqsForCCs += largerMajSeqsForCCs]]
                    footprintList[collapseDict[footPrintToCollapse]] = [footprintList[collapseDict[footPrintToCollapse]][0] + footprintList[footPrintToCollapse][0], footprintList[collapseDict[footPrintToCollapse]][1] + extractedMajSeqs]

                    # NB if this was an extraction case, then there should still be a maj dss for each of the ccts lef.

                    # 2 Adjust big so that the seqs of the small are removed from its profile
                    # the extracted maj dsss have already been removed
                    # I guess we now need to work out what to do with the big profile post extraction
                    # TODO check that my assumption about the big type being in the unsupported list are true
                    # I'm assuming that the big type was in the unsupported List
                    # The big type will have shrunk in size, I guess that if the big type is still > n then we
                    # can update its key in the unsupportedList. But if it is smaller then n then we should remove it from
                    # the unsupported list.
                    # newKey = smallerseqset not in common with largerseqset
                    updatedKey = collapseDict[footPrintToCollapse] ^ footPrintToCollapse
                    # I suppose we also need to check to see if the new key already exists.
                    if updatedKey in footprintList.keys():
                        # then this already exists and we need to add support to the already existant profile
                        # we should just be able to do simple additions here
                        footprintList[updatedKey] = [footprintList[updatedKey][0] + footprintList[footPrintToCollapse][0], footprintList[updatedKey][1] + footprintList[footPrintToCollapse][1]]

                    else:
                        footprintList[updatedKey] = footprintList[footPrintToCollapse]

                    # Now the updated versin of the bigfootprint has been added to the footprintList, we can delete
                    # old version
                    del footprintList[footPrintToCollapse]

                    # Now remove from unsupported list and place back in depending on the length of the profile compared to n
                    unsupportedList.remove(footPrintToCollapse)
                    if len(updatedKey) >= n:
                        unsupportedList.append(updatedKey)

        if n > 2:
            # We only need to attempt the further collapsing if there are unsupported types to collapse
            # else move onto the next n
            if unsupportedList:
                # Now we are left with the footprints that are still in the unsupported list
                # Get all the seqs that are found in the footprints that are of length n


                # we are going to change this.
                # Currently we are generating an enormous generator in the itertools.combinations
                # This is because we are using all of the sequences found in all of the length n types
                # However, I see non reason why we should be using combinatations of sequences found in different profiles.
                # This is not in keeping with the biology.
                # Instead the synthetic types should be made from
                # combinattions of sequences found in each profile
                # This will better represent the biology and also greatly cutdown on computational cost
                # And it will mean that we can multiprocess this
                setOfRefSeqsOfSmpSeqsFoundInCCsOfLengthN = set()
                listOfTypesOfLengthN = [footprint for footprint in footprintList.keys() if len(footprint) == n]
                # Only carry on if we have lengthN footprints to get sequences from
                if listOfTypesOfLengthN:

                    # TODO here is where we should start the new approach.
                    # add the nLength Types into a MP list
                    # then add the stops
                    # then in the worker create a dict that is same as the one below, only the using the nLengthtype instead
                    # Create the queues that will hold the sample information
                    print('starting to generate Synthetic footprints for collapse')

                    taskQueue = Queue()
                    collapseNMerDictManager = Manager()
                    collapseNMerDict = collapseNMerDictManager.dict()
                    #outputQueue = Queue()

                    for nLengthType in listOfTypesOfLengthN:
                        taskQueue.put(nLengthType)

                    numProcessors = nprocessors
                    for N in range(numProcessors):
                        taskQueue.put('STOP')

                    allProcesses = []

                    for N in range(numProcessors):
                        p = Process(target=workerDiscoveryOne, args=(taskQueue, collapseNMerDict, n))
                        allProcesses.append(p)
                        p.start()

                    for p in allProcesses:
                        p.join()

                    # outputQueue.put('STOP')

                    # collapseNMerDict = {}
                    # for dictionary in iter(outputQueue.get, 'STOP'):
                    #     collapseNMerDict.update(dictionary)

                    # for nLengthType in listOfTypesOfLengthN:
                    #     # nLengthType should be a frozen set so we will take the refSeqs it contains and add them to the setOfRefSeqs...
                    #     setOfRefSeqsOfSmpSeqsFoundInCCsOfLengthN.update([refSeq for refSeq in nLengthType])

                    #Now go through each of the (n-1)mer footprints and see if they fit into a footprint in the unsuported list
                    #This dict value = the synthetic footprint of length n-1, value, a list of unsupported types that are mastersets of the synthetic footprint

                    # collapseNMerDict = {frozenset(tup): [] for tup in itertools.combinations(setOfRefSeqsOfSmpSeqsFoundInCCsOfLengthN, n-1)}
                    print('Generated {} Synthetic footprints'.format(len(collapseNMerDict)))
                    # We create a dict that is key, type and value is the set of the majrefseqs found in the type
                    # this was being calculated once for every synthetic type. But the number of synthetic types are getting very large
                    # so by having this dict we should speed things up considerably
                    # This info is used to make sure that when collapsing types the maj sequences are shared between the
                    # larger and smaller types
                    majRefSeqsDict = {}
                    for nLengthType in unsupportedList:
                        listOfAnalysesSampleSequencesLarge = []
                        for majdssslist in footprintList[nLengthType][1]:
                            for dss in majdssslist:
                                listOfAnalysesSampleSequencesLarge.append(dss)
                        listOfMajRefSeqsLarge = [smpl.referenceSequenceOf for smpl in listOfAnalysesSampleSequencesLarge]
                        setOfMajRefSeqsLarge = set(listOfMajRefSeqsLarge)
                        majRefSeqsDict[nLengthType] = setOfMajRefSeqsLarge

                    ## Items for printing out our progress
                    total = len(collapseNMerDict)*len(unsupportedList)
                    count = 0
                    printVal = 0.01
                    print('Checking new set of synthetic types')


                    for frTupKey in collapseNMerDict.keys(): # for each of the synthetic footprints

                        for nLengthType in unsupportedList:
                            ## Items for printing out progress
                            count += 1
                            percentComplete = count/total
                            if percentComplete > printVal:
                                print("%.2f" % percentComplete)
                                printVal = max([printVal + 0.01, percentComplete + 0.01])

                        #For each of the N-1mers see if they fit within the unsupported types
                        # if so then add the footprint into the n-1mers list
                            # we should check for the maj rule again. I.e. we should want all of the Majs of the footprint in Q to be in the kmer footprint
                            # maybe recycle code for this
                            if frTupKey.issubset(nLengthType):
                                if majRefSeqsDict[nLengthType].issubset(frTupKey):
                                    collapseNMerDict[frTupKey].append(nLengthType)
                    # Here we have a populated dict.
                    # Order dict
                    # We are only interseted in kmers that were found in the unsupported types
                    # Because alot of the kmers will be found in the sequences they originated from
                    # require the kmer to be associated with more than 1 cladecollection
                    # I am not able to express my logic here perfectly but I am sure that this makes sense
                    listOfnKmersWithCollapsableFootprints = [kmer for kmer in collapseNMerDict.keys() if len(collapseNMerDict[kmer]) > 1]
                    # No need to continue if there are no footprints that match the nKmers
                    if listOfnKmersWithCollapsableFootprints:
                        # We now go through the kmers according to the number of footprints they were found in
                        orderedListOfPopulatedKmers = sorted(listOfnKmersWithCollapsableFootprints, key=lambda x: len(collapseNMerDict[x]), reverse=True)

                        for kmer in orderedListOfPopulatedKmers:
                            #for each footprint in their lists
                            for bigFootprintToCollapse in collapseNMerDict[kmer]:
                                # chek to see if the footprint is still in the unsupportedlist
                                # if it isn't then it has already been associated to a new footprint
                                #   pass it over
                                if bigFootprintToCollapse in unsupportedList: # Then this footprint hasn't been collapsed anywhere yet
                                    # if the footprint exists in footprintlist then add the CC details to it
                                    if kmer in footprintList.keys():
                                        footprintList[kmer] = [footprintList[kmer][0] + footprintList[bigFootprintToCollapse][0],
                                                               footprintList[kmer][1] + footprintList[bigFootprintToCollapse][1]]
                                        # then remove the footprint from the unsupported type
                                        # Also remove the bigFootprintprintToCollapse from the footprintList as its details are not in another key
                                        unsupportedList.remove(bigFootprintToCollapse)
                                        del footprintList[bigFootprintToCollapse]
                                    else:
                                        # if the footprint is not already in the footprint list then create newfootprint item and add it to the footprintlist and assign the
                                        # CC details to it.
                                        # remove the footprint from the unsupported type list
                                        # Also remove the bigFootprintprintToCollapse from the footprintList as its details are not in another key
                                        footprintList[kmer] = [footprintList[bigFootprintToCollapse][0], footprintList[bigFootprintToCollapse][1]]
                                        unsupportedList.remove(bigFootprintToCollapse)
                                        del footprintList[bigFootprintToCollapse]
                            # each unsupportedfootprint in the nmer has been collapsed
                        # each of the kmers have been gone through and all unsupportedfootprints that could fit into the kmers have been collapsed
                        # Time for the next n here!
                        # Check to see what is happening once we get to the smaller ns
                        # It may be that we can use the code from below to simply put the maj as the type for footprints that are still
                        # unassigned
        else:
            # At this point we may have some 2s or larger in the unsupported list still
            # we will collapse these to their maj ref seq type
            # Once we have pushed these into the 1s or made new 1s if the maj ref seq type did not already exist
            # then we should collect all of the 1 length footprints and put them in the supported list
            # Or actually in theory so long as we have deleted collapsed footprints from the footprintList we should
            # be able to return the footprint list
            for unsupFp in unsupportedList:
                # for each clade collection check to see if its maj dss is already a profile on its own
                # If it already exists then add to this. Else create new
                for i in range(len(footprintList[unsupFp][0])):
                    # TODO 08/12/17 we can now have multiple majdsss associated with a given profile
                    # So here we should go through each of the maj dsss that are associated with
                    # if the majRefSeq is a footprint add to the existant footprint
                    for majdss in footprintList[unsupFp][1][i]:
                        singleDivProfileInQ = frozenset([majdss.referenceSequenceOf])
                        if singleDivProfileInQ in footprintList.keys():
                            # then we need to add the support to the existant single maj div profile
                            # current single Div in q info = [[current list of CCts] + [cct in question], [current list of dsss lists] + [[majdssofcct]]
                            footprintList[singleDivProfileInQ] = [footprintList[singleDivProfileInQ][0] + [footprintList[unsupFp][0][i]], footprintList[singleDivProfileInQ][1] + [[majdss]]]
                    # if frozenset([footprintList[unsupFp][1][i].referenceSequenceOf]) in footprintList.keys():
                    #     footprintList[frozenset([footprintList[unsupFp][1][i].referenceSequenceOf])] = [footprintList[frozenset([footprintList[unsupFp][1][i].referenceSequenceOf])][0] + [footprintList[unsupFp][0][i]], footprintList[frozenset([footprintList[unsupFp][1][i].referenceSequenceOf])][1] + [footprintList[unsupFp][1][i]] ]
                        else: # create new footprintList entry
                            footprintList[singleDivProfileInQ] = [[footprintList[unsupFp][0][i]], [[majdss]]]
                # remove the unsupportedType from the footprintList
                del footprintList[unsupFp]
            # All unsupported footprints have been associated to their maj and deleted from the footprintList
            return footprintList

    return False

def workerDiscoveryOne(input, collapsenmerdict, n):
    for refSeqFrozenSet in iter(input.get, 'STOP'):
        tempDict = {frozenset(tup): [] for tup in itertools.combinations(refSeqFrozenSet, n - 1)}
        collapsenmerdict.update(tempDict)
        print('\rGenerated iterCombos using {}'.format(current_process().name), end='')
#########################################





###### PRFOILE ASSIGNMENT FUNCTIONS ######
def profileAssignment(nProcessors):
    '''Type assignment cycles through cladeCollections and then types. For each cladeCollections it looks to
    see if each type's footprint can be found. If the foot print is found then it checks to see whether the
    abundances of the type sequences are at ratios within the defining ratios of the type. We use a rule that no two
    types can be assigned to a cladeCollection if they share a DIV. The type which uses up more of the cladeCollections
    sequences is assigned to the cladeCollection. We use DIV proportions are considered in relation to the total
    abundance of the DIV sequences in the CC rather than all seqs in the CC'''


    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    print('Running inferFinalSymbiodiniumTypes()')


    # Get a list of all of the analysisTypes that have been discovered in the type discovery
    querySetOfAnalysisTypes = list(analysis_type.objects.filter(dataAnalysisFrom=analysisObj))
    # Get the IDs of these types
    analysisTypesIDs = [anType.id for anType in querySetOfAnalysisTypes]
    # Get a list of the CCs to check for types within
    querySetOfCCs = analysisObj.getCladeCollections()
    # Get dict of the CC footprint and Type footprints so that we don't have to look them up each time



    ################# FOOTPRINT DICT FOR EACH CC MP #################
    # Make a footprint dict for the CCs
    # This might be an opportunity to use a manager
    # the manager should in theory create a proxy for the object and allow it to be shared between processess
    # http://stackoverflow.com/questions/9436757/how-does-multiprocessing-manager-work-in-python

    manager = Manager()
    CCfootprintDict = manager.dict()

    CCInputQueue = Queue()

    for CC in querySetOfCCs:
        CCInputQueue.put(CC)

    numProcessors = nProcessors
    for N in range(numProcessors):
        CCInputQueue.put('STOP')

    allProcesses = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()
    # Then start the workers
    # workerDiscoveryTwoWorker will process the CCs and pass this info on to the second queue which
    # workerDiscoveryTwoListener will work on
    # Finally workerDiscoveryTwoListener will output its results to the third queue
    for N in range(numProcessors):
        p = Process(target=workerAssignmentOne, args=(CCInputQueue, CCfootprintDict))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

    #################################################################







    ############# INTRA ABUND DICT FOR EACH CC MP ################
    # Make dict of intraAbundances for each CC #
    managerTwo = Manager()
    CCRefSeqAbundDict = managerTwo.dict()

    CCInputQueue = Queue()

    for CC in querySetOfCCs:
        CCInputQueue.put(CC)

    numProcessors = nProcessors
    for N in range(numProcessors):
        CCInputQueue.put('STOP')

    allProcesses = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()

    for N in range(numProcessors):
        p = Process(target=workerAssignmentTwo, args=(CCInputQueue, CCRefSeqAbundDict))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

    apples = 'pausePoint'
    ##############################################################




    ################ TYPE ASSIGNMENT ################
    # TODO  I believe that this is causing problems due to the MP access of the database.
    # This actually makes sense. We go clade collection by cladeCollection,
    # but within each clade collection we go analysis_type by analysis_type.
    # We will be modifying these analysis types by changing their .listOfCladeCollection parameters
    # As this is writing we may well come into conflicts.
    # Two options, hard wire this MP section to only use 1 core, or see if the conflict is cleared up by using a
    # postgresql sever instead of the sqlite server.

    # TO be sure of why this is happening. I will first re-run the MP analysis with iterations and quantify
    # How many mismatches happen between the CC ids attached to each analysis_type vs the number of Cladecollectiontypes
    # there are associated with each at. I will do this with the normal cores being used in this type assignment.
    # I will then do this again with the cores hard coded to 1 (for just this MP) and keep the cores for everyother
    # process. Hopefully this will then confim that it is indeed this part of the code that is causing our problems

    # The main bulk of the type assignment will happen in workerAssignmentThree
    print('\n\nStarting type assignment')

    # Make a footprint dict for the profiles

    antypefootManager = Manager()
    anTypeFootprintDict = antypefootManager.dict({anType.id: frozenset(anType.getOrderedFootprintList()) for anType in querySetOfAnalysisTypes})


    # Two manager objects that will be used in the multiprocessing as outputs

    # This is a list that will hold ID tuples to create clade_collection_type objects in bulk
    managerbulk = Manager()
    bulkCreateList = managerbulk.list()

    # I am no longer using this dict manager it was unstable
    # # This is a dict that will be used to add CCs to the types' cladeCollectionLists all at once
    # typecladecoldictmanager = Manager()
    # typeCladeColDict = typecladecoldictmanager.dict()
    # # Initial population of the managed dict
    # for ID in analysisTypesIDs:
    #     typeCladeColDict[ID] = []

    # Input Queue is CCs again
    profileAssignmentInputQueue = Queue()
    for CC in querySetOfCCs:
        profileAssignmentInputQueue.put(CC)
    for process in range(nProcessors):
        profileAssignmentInputQueue.put('STOP')

    allProcesses = []
    # CCRefSeqAbundDict is a manager dict that was created and populated earlier, as is CCfootprintDict
    # BulkCreateList and typeCladeColDict are managed objects that our outputs go into
    # We will then use these to bulk process at the end of the profile assignment
    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()
    for pro in range(nProcessors):
        p = Process(target=workerAssignmentThree, args=(profileAssignmentInputQueue, CCRefSeqAbundDict, CCfootprintDict, bulkCreateList, analysisTypesIDs, anTypeFootprintDict))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()


    ######################################################
    print('\n\nCreating cladeCollectionTypes in bulk...')
    # Create cladeCollectionTypes in bulk
    cctList = []
    for cct in bulkCreateList:
        cctList.append(clade_collection_type(analysisTypeOf=analysis_type.objects.get(id=cct[0]), cladeCollectionFoundIn=clade_collection.objects.get(id=cct[1])))
    clade_collection_type.objects.bulk_create(cctList)

    # The info in here is the a list of tuples with 0 being the id of the analysis type and 1 the cladeCollection id
    # We can easily do what we were doing with the dictionary manager with the above.
    # Create default  dict
    cladeCollectionTypeMakingDict = defaultdict(list)
    for atID, ccID in bulkCreateList:
        cladeCollectionTypeMakingDict[atID].append(ccID)

    for typeIDKey in cladeCollectionTypeMakingDict.keys():
        analysisTypeInQ = analysis_type.objects.get(id=typeIDKey)
        analysisTypeInQ.addCCListToCladeCollectionList(cladeCollectionTypeMakingDict[typeIDKey])
        analysisTypeInQ.save()

    # #Add the CC's to the profiles
    # #TODO don't use the manager dict it is causing so many problems
    # # instead just get the information required to do the below bit of code from the above
    # # bulk create list
    # for typeKey in typeCladeColDict.keys():
    #     if typeCladeColDict[typeKey]:
    #         analysisTypeInQ = analysis_type.objects.get(id=typeKey)
    #         analysisTypeInQ.addCCListToCladeCollectionList(typeCladeColDict[typeKey])
    #         analysisTypeInQ.save()

    # TODO this is likely not parallelable because it requires writing to the database.
    # updateTypeAttributes updates
    # Update the types' counts and ratios from the up-to-date CC lists
    # See the updateTypeAttributes comment for more details of the use



    for antype in analysis_type.objects.filter(id__in=analysisTypesIDs):
        if antype.listOfCladeCollections:
            antype.updateTypeAttributes()
            antype.save()
            print('\rUpdating type {}'.format(antype.name), end='')


    # TODO this might be parallelable
    multiModalDetection(False)


    analysisObj.analysisTypesAssigned = True
    analysisObj.save()


    return

def workerAssignmentOne(input, dict):
    for CC in iter(input.get, 'STOP'):
        print('\rfootprintdict {}\t{}'.format(CC, current_process()), end='')
        dict[CC] = CC.footprint()

def workerAssignmentTwo(input, dict):
    for CC in iter(input.get, 'STOP'):
        print('\rMaking refSeqAbundDictElement for {}'.format(str(CC)), end='')
        # Make a tempDict that will be returned with the CC as Key in the main CCrefSeqAbundDict
        tempDict = {}
        # get a list of dataSetSampleSequenceTwos and make a refseq to abundance dict
        dssstList = list(data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=CC))
        for dsss in dssstList:
            tempDict[dsss.referenceSequenceOf] = dsss.abundance
        dict[CC] = tempDict

def workerAssignmentThree(input, CCRefSeqAbundDict, CCfootprintDict, bulkCreateList, analysisTypesIDs, anTypeFootprintDict):
    # This is the main assignment process that checks each CC for which types fit within it.
    for CC in iter(input.get, 'STOP'):
        # if str(CC) == 'January2017-280por':
        #     apples = 'pears'


        print('\rProfile assignment {} by process {}'.format(CC, current_process().name), end='')
        # return the dictionaries for this CC for use later on
        # both the abundance of the refseqs in the CC and its footprint
        CCabundDict = CCRefSeqAbundDict[CC]
        CCfootprint = CCfootprintDict[CC]

        # list to keep track of candidate types
        candidateTypesList = []

        # I will cycle through the IDs to prevent having to update the query setofAnalysisTypesevery run
        for i in range(len(analysisTypesIDs)):
            # This is the analysis_type that we are trying to find in the CC
            anTypeInQ = analysis_type.objects.get(id=analysisTypesIDs[i])
            # DEBUG
            if str(anTypeInQ) == 'C15m-C15-C15ae':
                apples = 'pears'
            ###
            anTypeFootprint = anTypeFootprintDict[analysisTypesIDs[i]]

            # Check to see if the type's footprint is found in the CCfootprint.
            if anTypeFootprint.issubset(CCfootprint):
                # If this type contains more than a single seq
                if len(anTypeFootprint) > 1:
                    # This is key code. In ratioCheck we look to see if the refSeqs are found at correct relative
                    # proportions in the CC to be considered part of the type in question.
                    # We are using the abundance of the DIVs in proportion to the total abundance of DIV sequences
                    # in the CCs rather than the total seqs of the CC. This is good a necessary for finding low abundances
                    # of types
                    if ratioCheck(CC, anTypeInQ):
                        # If we pass the ratio check then we continue with the code.
                        pass
                    else:
                        # If we fail the ratio check then this type is not found in the CC
                        continue



                else: # If Type is only one DIV # TODO can't we just look this up with the CCRedSeqAbundDict?
                    # ensure that the seq is found at above 0.05 of the CC
                    seqsInCC = data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=CC)
                    abundanceOfSeqsInCC = sum(aSS.abundance for aSS in seqsInCC)
                    abundanceOfTypeSeqInQ = seqsInCC.get(
                        referenceSequenceOf=list(anTypeFootprintDict[analysisTypesIDs[i]])[0]).abundance
                    if abundanceOfTypeSeqInQ / abundanceOfSeqsInCC < 0.05:
                        continue


                # Here is a type that meets criteria for being included in a sample
                # Previously we were checking to see if [1] the current final types were subsets of the new type,
                # or [2] if the new footprint was a subset of any of the current types
                # However, we have implemented a new control here. The new control is that if two profiles share DIVs,
                # the profile which uses up the most of the samples sequences will be accepted.
                # This new control does away for the need of [1] or [2].



                # List of types to be deleted
                typeToDel = []
                smaller = False
                # candidateTypesList is the list of potentials that are already due to be associated with the CC
                for candidate in candidateTypesList:

                    ##### TODO Delete these they are just for debugging purposes
                    candidateName = analysis_type.objects.get(id=candidate)
                    anTypeName = analysis_type.objects.get(id=analysisTypesIDs[i])
                    #####

                    candidateFootprint = anTypeFootprintDict[candidate]
                    # I think we are comparing sets of refseqs (profiles)
                    # Check to see if the two footprints have any refSeqs in common
                    if bool(anTypeFootprint & candidateFootprint):

                        # then the two candidates share intras and we need to see which one has the most number of sequences in it
                        # we will remove the candidate using up the less sequences from the candidacy for this type
                        totanType = sum([CCabundDict[refSeq] for refSeq in anTypeFootprint])
                        totCan = sum([CCabundDict[refSeq] for refSeq in candidateFootprint])

                        # If the new type has bigger seq total then current candidate, del candidate from candidateTypeList
                        if totanType > totCan:
                            typeToDel.append(candidate)
                        # else if the new type is smaller then do nothing (don't add new type to candidateTypeList)
                        else:
                            smaller = True

                # Remove smaller candidates from candidateTypeList
                for toDel in typeToDel:
                    candidateTypesList.remove(toDel)

                # If all criteria met add new type as candidate
                if smaller == False:
                    candidateTypesList.append(analysisTypesIDs[i])

        # Successful type candidates will be associated with the CC by creating a cladeColletionType that associates
        # between the analysis_type and the cladeCollection


        # To speed up the creation of cladeCollectionTypes, they will be created in bulk. To store the information of
        # which cladeCollectionTypes should be created I will create a tuple which will be the IDs of the analysisTypeOf
        # and cladeCollectionFoundIn paramerters. These will be exported via the bulkCreateList which is a manged list.

        # Previously, we had been associating a cladeCollection to the analysis_type one cladeCollection at a time.
        # We have sped this up by associating multiple cladeCollections to analysisTypes in bulk.
        # To keep track of which CCs need associating to which types within the MPing we have created the managed dict
        # typeCladeColDict. This has the typeID as the key and a list as the value which will contain the CC ids that
        # are to be associated to this type. One MPing is complete the CCs can then be associated to the CC in bulk.

        # if there are some types found in the CC...
        for successfulCandidateID in candidateTypesList:  # These are the types which will need adding to the CC (i.e. create a CCType)
            # Now associate the new types to the CC (i.e. create a CCType per type found in the CC)
            # newCladeCollectionType = clade_collection_type(analysisTypeOf=successfulCandidate, cladeCollectionFoundIn=CC)
            bulkCreateList.append((successfulCandidateID, CC.id))


            # I used to have a manager dict here that I would have a typeID as key and list of CCIDs as value
            # But this really wasn't stable and was causing issues no mater what I tried.
            # I have left the code and comments in below VVVVVVVV

            # temp_list = typeCladeColDict[successfulCandidateID]
            # temp_list.append(CC.id)
            # typeCladeColDict[successfulCandidateID] = temp_list
            #NB This was/is causing all sorts of problems
            # The work around i was using was also not working (commented out above).
            # https://stackoverflow.com/questions/8640367/python-manager-dict-in-multiprocessing
            # I will try implementing this
            # typeCladeColDict[successfulCandidateID] = typeCladeColDict[successfulCandidateID] + [CC.id]
            # This was causing problems and I don't want to have to adjust Manager to include default dict
            # so I am going around the problem with the above, inellegant code.
            # https://stackoverflow.com/questions/26367812/appending-to-list-in-python-dictionary
            #typeCladeColDict[successfulCandidateID].append(CC.id)

    return

def ratioCheck(CC, antype):# allowance of 0.00 and no update
    '''
    This method will assess the type that is currently about to be associated with the cladeCollection to see if the
    abundances in which the type's intras are found make sense.
    For example, if the type is C1/Otu1234 and the average abundances of those two are something like .65 .45 then if we find that the abundances
    of the intras in this sample are like 0.05 and .45 then these intras are likely not due to this type and they will not be included into this sample.
    This function will only take symtypes that have a footprint that contains two or more defining intras.
    I think we will use a set of ratios to define a type. The ratio for each intra will always be the intra in question to the total abundances of the DIVs in the CC.


    '''
    # We will not allow any deviation from the ratio.
    deviationPercentage = 0.00

    # list of ref seqs that make up the type in Q in order of abundance across all
    # cladeCollections the type was found in initially
    orderedFootprintListOfType = antype.getOrderedFootprintList()

    # List that will hold the ratios of the sampleSeqs for the clade collection in q
    footprintInQRatiosList = []

    # For each refseq of type, get the abundance in the CC in Q and sum these. We will use this to calculate
    # what relative abundance each refseq is at relative to only the total abundances of sequences of the CC found in the type in Q
    # In other words, we are not working out the type's refseqs in proportion to all the CCs seqs, just those found in the type (DIVs).
    # Get list of the abundances in order of orderedFoorptintListOfType
    sampleSeqAbundances = []
    for refSeq in orderedFootprintListOfType:
        # Get the abundance of the current refSeq in Q
        sampleSeqAbundances.append(data_set_sample_sequence.objects.get(cladeCollectionTwoFoundIn=CC, referenceSequenceOf=refSeq).abundance)
    total = sum(sampleSeqAbundances)
    footprintInQRatiosList = [x/total for x in sampleSeqAbundances]

    # for each sequence abundance
    # so we now know that for the CC in question we are working out the seq abundanes relative to the total seqs
    # found in the type only, not realtive to all the seqs of the CC - this is good!
    # But now I want to check how we have noted down the MaxMin ratios of the types in the first place. If for each of
    # the DIVs we have divided by the total seqs in the CC then we are comparing two different things.
    # I will check this now.
    # If this is all good then I should be sure to make a note of this in the SymPortal documentations.
    # I have checked this. We are all good. I was mis-remembering. I will aim to make this clearer in the logic doc now.
    currentMaxMinRatioForType = antype.getMaxMinRatios()
    for refSeq in orderedFootprintListOfType[:]:
        currentMax = currentMaxMinRatioForType[orderedFootprintListOfType.index(refSeq)][0]
        currentMin = currentMaxMinRatioForType[orderedFootprintListOfType.index(refSeq)][1]
        maxRatioCutoff = currentMax + deviationPercentage
        minRatioCutoff = currentMin - deviationPercentage

        ratioInQ = footprintInQRatiosList[orderedFootprintListOfType.index(refSeq)]


        if ratioInQ < minRatioCutoff or ratioInQ > maxRatioCutoff:
            return False

    return True

def multiModalDetection(initialRun):


    # If initialRun = True, then this is the first time we are running multiModal Detection and we have not been
    # through the first iteration of assigning types yet. As such we should still be working with the
    # listOfCladeCollectionsFoundInInitially rather than the later listOfCladeCollections
    # This will be important when creating the new types after a successful split candidate is found
    # Here we identify if there is a bionomial distribution in coDom types between the two most abundant intras
    # If we find a binomial distribution which meets our parameters for clearly being two separate distributions
    # then we separate the type into two new types
    if not initialRun: # We only want types that have a cladeCollectionList
        # both empty string and None return False so this is types that have CCs associated to them only
        listOfAllAnalysisTypes = [anObj for anObj in analysis_type.objects.filter(dataAnalysisFrom=analysisObj) if anObj.listOfCladeCollections]
    else: # No need to check if initial run as all types have initial clade collections
        listOfAllAnalysisTypes = [anObj for anObj in analysis_type.objects.all().filter(dataAnalysisFrom=analysisObj)]
    analysisTypesToBeDeleted = []
    # This is a list of the types that have already been successfuly split
    # If this has types in it then we want to re-assess to see if there are further types to split
    # We use this so that we don't have to assess types that have already been found to not be multimodal
    typesSplit = []
    numTypesSplit = 1
    listOfTypesToAnalyse = []
    while numTypesSplit > 0: # Only continue this while we are finding types to split
        numTypesSplit = 0
        if typesSplit:
            listOfTypesToAnalyse = typesSplit
        else:
            listOfTypesToAnalyse = listOfAllAnalysisTypes
        typesSplit = []
        for k in range(len(listOfTypesToAnalyse)):
            orderedFootprintList = listOfTypesToAnalyse[k].getOrderedFootprintList()
            if len(orderedFootprintList) == 1:
                # If there is only one intra then no need to check
                continue
            print('\rAssessing type {} for multimodal characteristics'.format(listOfTypesToAnalyse[k].name), end='')
            # I think we start by just checking the coDom types as these are the most obvious candidates for being
            # bimodal but eventually  we need to check all types
            # if listOfTypesToAnalyse[k].coDom:
            if initialRun:
                cladeCollectionsInType = listOfTypesToAnalyse[k].getCladeCollectionsFoundInInitially()
            else:
                try:
                    cladeCollectionsInType = listOfTypesToAnalyse[k].getCladeCollections()
                except:
                    apples = 'asd'
            # 2D list where each list is cc and each item in list is the relative abundance of the seqinQ in the CC
            # where seqs are in the order of orderedFootprintList
            seqRatiosForType = listOfTypesToAnalyse[k].getRatioList()
            # we can only split a type if we have sufficient instances of that type to look at
            # we will make the cutoff 10 for the time being
            if len(cladeCollectionsInType) > 9:

                # Now that we are using relative abundances rather than ratios we check all intras in the list
                for refSeq in orderedFootprintList:
                    # This index will give us the current intra ratio we are checking
                    index = orderedFootprintList.index(refSeq)
                    # This list will be a list of the abundances of the CC list for the refSeq in Q
                    # These are the values that will be plotted virtually in the histogram
                    listOfRatios = [ratioVal[index] for ratioVal in seqRatiosForType]
                    # linspace (start, stop, numvalues)
                    x_grid = np.linspace(min(listOfRatios)-1, max(listOfRatios)+1, 2000)
                    try:
                        # This was causing errors when we were trying to calculate it for a list of 0s.
                        # As such I have put in a harmless continue as we be sure that the ref_seq in question is
                        # not going to be a reason to split the type as it is all 0s.
                        kde = gaussian_kde(listOfRatios)
                    except:
                        continue
                    pdf = kde.evaluate(x_grid)

                    # returns index of the local max's in pdf
                    # Using these index on x_grid will give you x of maxs, use on pdf will give you y
                    c = list((np.diff(np.sign(np.diff(pdf))) < 0).nonzero()[0] + 1)
                    modes = len(c)

                    # UNCOMMENT for visualisation of types
                    # plotHists(pdf, x_grid, listOfRatios, listOfTypesToAnalyse[k].name)

                    # If this appears to be a bimodal distribution

                    if modes == 2:
                        # Must be sufficient separation between the peaks in x axis
                        xDiffValid = False
                        if x_grid[c[1]] - x_grid[c[0]] > 0.7:
                            xDiffValid = True
                        # plotHists(pdf, x_grid, listOfRatios, listOfTypesToAnalyse[k].name)
                        # Must also be sufficient diff between minima y and small peak y
                        # This represents the x spread and overlap of the two peaks
                        d = list((np.diff(np.sign(np.diff(pdf))) != 0).nonzero()[0] + 1)  # max and min indices

                        if min([pdf[d[0]], pdf[d[2]]]) == 0:
                            xDiffValid = False
                        else:
                            if pdf[d[1]] / min([pdf[d[0]], pdf[d[2]]]) > 0.85:  # Insufficient separation of peaks
                                xDiffValid = False
                                # Then here we have a candidate for separating into two types
                                # Seperate the samples either side of the minima x axis
                                # Find the minimum x value and then use this ratio as the separator for the two modes
                                # Because the ratio information was put in the same order as the samplesFoundInAsFinal
                                # we can work out which samples are which side of the ratio
                                # Given that this relies on the order of the ratio list being the same as that of the cladeCollection
                                # We need to check with the new HumeFst code that this order is the same or find a more definitive way to keep track of which
                                # clade collection the samples came from.
                                # I have checked this and the ratios list should indeed be in the same order as the clades list.
                                # So the strategy below is robust.
                        if xDiffValid:
                            if initialRun:
                                orderedListOfCCsInType = listOfTypesToAnalyse[k].getCladeCollectionsFoundInInitially()
                            else:
                                orderedListOfCCsInType = listOfTypesToAnalyse[k].getCladeCollections()
                            CCsForTypeA = []
                            CCsForTypeB = []
                            # returns the index of local max and mins in pdf
                            # index 1 is the min
                            minX = x_grid[list(((np.diff(np.sign(np.diff(pdf))) != 0).nonzero()[0] + 1))[1]]
                            # for each sample assign ratio to one of the two new types
                            for i in range(len(listOfRatios)):
                                if listOfRatios[i] < minX:
                                    try:
                                        CCsForTypeA.append(orderedListOfCCsInType[i])
                                    except:
                                        apples = 'gotworms'
                                else:
                                    CCsForTypeB.append(orderedListOfCCsInType[i])

                            # Only create new types if each type is supported by three samples

                            # get list of cladCollectionTypes that are associated with this analysis_type
                            # create two new analysisTypes according to the cc split
                            # (we will need to claculate the setofMajSeqs inorder to initiate the new analysis type)
                            # for each of the cladeCollection in the analysisTypes create a new clade_collection_type
                            # Delete the old analysis_type -- NB this should delete the cladCollectionTypes as well.
                            # We should probably do a try: except: to try to delete cladeCollectionTypes that were
                            # associated with the analysis types to make sure that they are deleted




                            if len(CCsForTypeA) >= 3 and len(CCsForTypeB) >= 3:
                                if listOfTypesToAnalyse[k].name == '21899/21901/21905-21916-21934':
                                    apples = 'asdf'
                                print('Splitting type {}'.format(listOfTypesToAnalyse[k].name))
                                # create two new analysisTypes according to the cc split
                                # (we will need to claculate the setofMajSeqs inorder to initiate the new analysis type)
                                # Get set of representative ref seqs of the maj seqs from each cc collection

                                #TODO check my logic here but,  the maj in the clade collection is not necessarily in the type
                                # so therefore it won't count towards whether a type is codom or not.
                                # we need to find the most abundant refseq of the CC that is in the type
                                # for each of the CCsForTypeX lists, go through ordered list generated and if found in the initial types footprint then add
                                # to the majSeqRefSeqSetX
                                # I have fixed this now
                                # TODO this is currently very slow and can likely be imporved. One option would be
                                # just to check a small number of the seqs first rather than all of them and then
                                # if we don't find a sequence from the CC that is in the type then we can go and check
                                # them all.
                                # When all is false we only search the top 10 seqs of the clade collection
                                # do this first then if needs be search all
                                all = False
                                majSeqRefSeqSetA = set()
                                for i in range(2):# Do twice, once with only top 10, second if needed, with all
                                    for cc in CCsForTypeA:
                                        orderedFootprintOfRefSeqsByAbundance = cc.orderedListOfRefSeqs(all)
                                        for refSequence in orderedFootprintOfRefSeqsByAbundance:
                                            if refSequence in orderedFootprintList:
                                                majSeqRefSeqSetA.add(refSequence)
                                                all = False
                                                break
                                            else:
                                                all = True
                                    if all == False: #we found all the refseqs we need
                                        break

                                all = False
                                majSeqRefSeqSetB = set()
                                for i in range(2):  # Do twice, once with only top 10, second if needed, with all
                                    for cc in CCsForTypeB:
                                        orderedFootprintOfRefSeqsByAbundance = cc.orderedListOfRefSeqs(all)
                                        for refSequence in orderedFootprintOfRefSeqsByAbundance:
                                            if refSequence in orderedFootprintList:
                                                majSeqRefSeqSetB.add(refSequence)
                                                all = False
                                                break
                                            else:
                                                all = True
                                    if all == False:  # we found all the refseqs we need
                                        break

                                #majSeqRefSeqSetA = frozenset([smpSeq.referenceSequenceOf for smpSeq in [cc.maj() for cc in CCsForTypeA]])
                                #majSeqRefSeqSetB = frozenset([smpSeq.referenceSequenceOf for smpSeq in [cc.maj() for cc in CCsForTypeB]])


                                # Initiate analysisTypeA first
                                if len(majSeqRefSeqSetA) > 1: # Then coDom = True
                                    newAnalysisTypeA = analysis_type(dataAnalysisFrom=analysisObj, coDom=True, clade=CCsForTypeA[0].clade)
                                    newAnalysisTypeA.save()
                                else:
                                    newAnalysisTypeA = analysis_type(dataAnalysisFrom=analysisObj, coDom=False, clade=CCsForTypeA[0].clade)
                                    newAnalysisTypeA.save()
                                newAnalysisTypeA.setMajRefSeqSet(majSeqRefSeqSetA)
                                # If this is the first run of the multiModal detection then we are still
                                # working with the listOfCladeCollectionsFoundInInitially rather than the later listOfCladeCollections
                                # equally we will be calling initTypeAtttributes rather than updateTypeAttributes

                                if initialRun:
                                    newAnalysisTypeA.initTypeAttributes(listOfCC=CCsForTypeA, footprintlistofrefseqs=orderedFootprintList)
                                else:
                                    # before calling updateTypeAttributes we must set self.listOfCladeCollections and self.orderedFootprintList
                                    # NB that self.orderedFootprintList does not have to be ordered
                                    newAnalysisTypeA.listOfCladeCollections = ','.join([str(cc.id) for cc in CCsForTypeA])
                                    newAnalysisTypeA.orderedFootprintList = ','.join([str(refSeq.id) for refSeq in orderedFootprintList])

                                    newAnalysisTypeA.save()
                                    if newAnalysisTypeA.orderedFootprintList == '':
                                        apples = 'adsf'
                                    newAnalysisTypeA.updateTypeAttributes()
                                newAnalysisTypeA.save()
                                print('Creating type {} for multimodal characteristics'.format(newAnalysisTypeA.name))
                                # for each of the cladeCollection in the analysisTypes create a new clade_collection_type
                                # This was causing a lot of trouble because I was creating cladeCollectionTypes on the
                                # initial run which was bad. adding the cct but adding the CCs to the listOfCladeCollectionsFoundInInitially
                                # was causing problems during the type assignment.
                                # I will put in this if so that the ccts are only created for non-initial runs
                                if not initialRun:
                                    listOfCladeCollectionTypes = []
                                    for CC in CCsForTypeA:
                                        newcladeCollectionType = clade_collection_type(analysisTypeOf=newAnalysisTypeA, cladeCollectionFoundIn = CC)
                                        listOfCladeCollectionTypes.append(newcladeCollectionType)
                                    clade_collection_type.objects.bulk_create(listOfCladeCollectionTypes)


                                # Initiate analysisTypeB second
                                if len(majSeqRefSeqSetB) > 1:  # Then coDom = True
                                    newAnalysisTypeB = analysis_type(dataAnalysisFrom=analysisObj, coDom=True, clade=CCsForTypeB[0].clade)
                                    newAnalysisTypeB.save()
                                else:
                                    newAnalysisTypeB = analysis_type(dataAnalysisFrom=analysisObj, coDom=False, clade=CCsForTypeB[0].clade)
                                    newAnalysisTypeB.save()
                                newAnalysisTypeB.setMajRefSeqSet(majSeqRefSeqSetB)
                                if initialRun:
                                    newAnalysisTypeB.initTypeAttributes(listOfCC=CCsForTypeB, footprintlistofrefseqs=orderedFootprintList)

                                else:
                                # before calling updateTypeAttributes we must set self.listOfCladeCollections and self.orderedFootprintList
                                # NB that self.orderedFootprintList does not have to be ordered
                                # If this is the first run of the multiModal detection then we are still
                                # working with the listOfCladeCollectionsFoundInInitially rather than the later listOfCladeCollections
                                # equally we will be calling initTypeAtttributes rather than updateTypeAttributes
                                    newAnalysisTypeB.listOfCladeCollections = ','.join([str(cc.id) for cc in CCsForTypeB])
                                    newAnalysisTypeB.orderedFootprintList = ','.join([str(refSeq.id) for refSeq in orderedFootprintList])
                                    newAnalysisTypeB.save()
                                    if newAnalysisTypeB.orderedFootprintList == '':
                                        apples = 'adsf'
                                    newAnalysisTypeB.updateTypeAttributes()
                                newAnalysisTypeB.save()
                                print('Creating type {} for multimodal characteristics'.format(newAnalysisTypeB.name))
                                # for each of the cladeCollection in the analysisTypes create a new clade_collection_type
                                if not initialRun:
                                    listOfCladeCollectionTypes = []
                                    for CC in CCsForTypeB:
                                        newcladeCollectionType = clade_collection_type(analysisTypeOf=newAnalysisTypeB, cladeCollectionFoundIn=CC)
                                        listOfCladeCollectionTypes.append(newcladeCollectionType)
                                    clade_collection_type.objects.bulk_create(listOfCladeCollectionTypes)
                                # Both analysis types are initialized here
                                # Delete the old analysis_type -- NB this should delete the cladCollectionTypes as well.
                                # The actual deleting should be done outside of the analysis_type for loop (k loop)
                                # To do this we will keep a list with the types to be deleted
                                analysisTypesToBeDeleted.append(listOfTypesToAnalyse[k])
                                typesSplit.extend([newAnalysisTypeA, newAnalysisTypeB])
                                numTypesSplit += 1
                                #when we are looking at doing more than just the first intra
                                # we should break from the for loop passing through all of the intras as we have already
                                # found a type that should be split
                                break
    # if binomials were found then it might be worth passing back through the types incase some of the
    # new types have further splits in them
    # to do this make a list of types that were split and then only check these
    # also keep track of numtypes split. then do the modal detection while, types split >0
    # you could do if listOfTypesSplit, to see whether there is a list to work off
    # else just do all types.

    # Delete the redundant analysisTypes
    #cladeCollectionsFirst
    cctToDelete = clade_collection_type.objects.all().filter(analysisTypeOf__in=analysisTypesToBeDeleted)
    for cc in cctToDelete:
        cc.delete()
    for anlType in analysisTypesToBeDeleted:
        anlType.delete()


    #NB We will need to be careful of the potential for two types to now have the same name.
    # I am not going to fix this now but we will need to fix this before we do the user
    # output as this is likely what they will use as their unique identifier
    # For the sake of the code though we will deal with the type instance or the
    # id so they will always be unique
    return

def plotHists(pdf, x_grid, newlist, typename):
    plt.interactive(False)
    fig, ax = plt.subplots(2, sharex=True)
    ax[0].hist(newlist, 100)
    ax[0].set_title(typename)
    ax[0].set_xlim([-2, 4])
    ax[1].plot(x_grid, pdf, color='blue', alpha=0.5, lw=3)
    plt.show()
    return
###############################





###### PROFILE COLLAPSE FUNCTIONS ######
def profileCollapse():


    # create the analsysisGroup instances and assign the analysisTypes to the groups
    #UNCOMMENT
    assignProfilesToGroups()


    #removeSubsetTypes()
    # Create a list of the groups so that we can maintain order
    listOfGroups = [anGrp for anGrp in analysis_group.objects.filter(dataAnalysisFrom=analysisObj)]
    i = 0
    initial = True
    listNotToCollapse = []
    while i < len(listOfGroups): # For each group
        ##### FOR DEBUG

        # groupinQ = analysis_group.objects.get(id=832)

        ######
        skip = False
        print('Assessing collapse within a group that needs a name')
        anlTypesInGroup = analysis_type.objects.filter(analysisGroupOf=listOfGroups[i])

        #### We will have the first thing which makes sure that types are sufficiently differnt in their footprint
        ##### but then we will have a second which will collapse if they share a critical number of types
        # Its hard to know what this number should be but
        # I think we should ass up all of the CC that the two types are found in and then count the num
        # that they have in common. If the num in common / total num > 0.5 or if num in common > 6 then we collapse
        # hopfully we can just reuse the code below
        # anlTypesInGroup = analysis_type.objects.filter(analysisGroupOf=listOfGroups[i])
        bob = [str(atg) for atg in anlTypesInGroup]
        if '22118-21978-22150' in bob:
            apples = 'asdf'

        # namedList = [an.name for an in anlTypesInGroup]
        # if 'D1/D4-D6-D17-D2.2-D2' in namedList:

        if len(anlTypesInGroup) > 1:





        # if len(anlTypesInGroup) > 1:
            # initial true is for first time looking at this group generate all the defRefSeqIntra info from scratch
            # but when we come back to here on an iteration looking at the same group after having collapsed
            # two groups to form a new one, we don't want to have to define all from scratch
            # we only want to update the intra info for the newly created type (the collapsed types will have
            # been removed)
            if initial == True:
                # Key = analysis_type, value = list where each item is the average ratio of the intra in question to the maj
                # these values are in order of abundance of defining intra across all of the types in the group
                # if listOfGroups[i].id == 832:
                #     apples = 'asdf'
                infoListDict, orderedListOfDefiningRefSeqs = getDefRefSeqIntraAbundanceForAllProfiles(analysisGrporAnType=listOfGroups[i], isinitial=initial)

                # Create the CCRefSeqAbundDict which is a dict that has a CC as key and the value is a dict that
                # is a referenceSeqeuce as key to its abundance in the type for the CCinQ
                listOfCladeCollectionsWithTypesFromThisGroup = clade_collection.objects.filter(clade_collection_type__analysisTypeOf__analysisGroupOf=listOfGroups[i])
                CCRefSeqAbundDict = createCCRefSeqAbundDict(listOfCladeCollectionsWithTypesFromThisGroup)

                # infoListDict, orderedListOfDefiningRefSeqs = getDefRefSeqIntraAbundanceForAllTypes(analysisGrporAnType=groupinQ, isinitial=initial)
                initial = False # set to false so that from now on, unless i is increased, we update the infoList
            # This list will hold those types which pass the collapse test but had an empty consesnusfootprint
            # and as such should not be collapsed

            for a, b in itertools.combinations(anlTypesInGroup, 2):
                # if {a.name, b.name} == {'D1/D4-D6-D17-D2.2-D2', 'D4/D1/D6-D2.2-D2'}:
                #     apples = 'asdf'
                # if frozenset([a,b]) in listNotToCollapse:
                #     continue



                collapseStatus = assessProfileCollapse(a, b, infoListDict, orderedListOfDefiningRefSeqs)


                if collapseStatus == 'typeProfileCollapse':

                    # want to make sure that the new list of CladeCollections hold each CC only once
                    newCCColection = set()
                    newCCColection.update(a.listOfCladeCollections.split(','), b.listOfCladeCollections.split(','))
                    newCCColection = [int(x) for x in list(newCCColection)]

                    # TODO it may be that we can imporve this
                    # we could make a list of the def intras in both the foot prints, make a set of these
                    # and then see which of these intras are found in all of the CCs
                    # have a think about this with a fresh brain

                    ########
                    # I am going to try to implement the above improvement here
                    # make a list of the def intras in both the footprints, make a set of these
                    allDefIntras = set()
                    allDefIntras.update(a.getOrderedFootprintList(), b.getOrderedFootprintList())
                    listOfAllCCs = clade_collection.objects.filter(id__in=newCCColection)
                    consensusFootprint = allDefIntras
                    for CC in listOfAllCCs:
                        ccfootprint = CC.footprint()
                        consensusFootprint = consensusFootprint & ccfootprint
                        # at this point consensusFootprint should be a footprint that is made up of the defining intras
                        # from type a and b that are are found in all of the cladeCollections
                        # TODO this is causing a problem, it is possible to have two types that do not share a maj. e.g. 22098 and 21978.
                        # in this case they should not be collapsed. Because although they share similar ratios of the majority of their
                        # intras (not necessarily their defining intras) if they share no common maj then this means that
                        # all samples they were found in have the one maj in the first type and a different single maj in the second
                        # this would be a strange situation to collapse into one type.
                        ########
                    # so here we have a conditional where if the the consensusFoorprint is empty we abort the collapse
                    if not consensusFootprint:
                        listNotToCollapse.append(frozenset([a,b]))
                        continue
                    newMajSeqRefSeqSet = set()
                    newMajSeqRefSeqSet.update(a.getMajRefSeqSet(), b.getMajRefSeqSet())
                    if len(newMajSeqRefSeqSet) > 1:# codom
                        newAnalysisType = analysis_type(dataAnalysisFrom=analysisObj, coDom=True, clade=a.clade, analysisGroupOf=listOfGroups[i])
                    else:
                        newAnalysisType = analysis_type(dataAnalysisFrom=analysisObj, coDom=False, clade=a.clade, analysisGroupOf=listOfGroups[i])
                    print('Collapsing types {} and {}'.format(str(a), str(b)))
                    newAnalysisType.setMajRefSeqSet(newMajSeqRefSeqSet)




                    #newFootPrintList = list(set(a.getOrderedFootprintList()) & set(b.getOrderedFootprintList()))
                    newFootPrintList = list(consensusFootprint)
                    newAnalysisType.listOfCladeCollections = ','.join([str(x) for x in newCCColection])
                    # New anType will have footprint that is the shared defining intras of the two types being collapsed
                    # This should cover all cases. e.g. when one big and one small, big will collapse to small
                    # when same size but different footprints, only common intras that are defining in both will stay as new footprint
                    # In this way we maintain the golden rul that any sample with a give type will definitely contain
                    # all of the sequences that define the type
                    # I have updated the above by doing the conesnsus Footprint

                    #NB orderedFootprintlist doesn't actually have to be ordered. That will
                    # happen when we update type attributes
                    newAnalysisType.orderedFootprintList = ','.join([str(refSeq.id) for refSeq in newFootPrintList])
                    if newAnalysisType.orderedFootprintList == '':
                        apples = 'adsf'
                    newAnalysisType.save()
                    newAnalysisType.updateTypeAttributes()
                    newAnalysisType.save()
                    print('Creating type {}'.format(str(newAnalysisType)))
                    # Delete the cladeCollectionTypes that were associated with the collapsing types
                    clade_collection_type.objects.filter(analysisTypeOf__in=[a, b]).delete()
                    # Get query set of cladeCollections associted with newType and create
                    # clade_collection_type for each
                    newCCs = clade_collection.objects.filter(id__in=newCCColection)
                    for CC in newCCs:
                        newcladeCollectionType = clade_collection_type(analysisTypeOf=newAnalysisType, cladeCollectionFoundIn=CC)
                        newcladeCollectionType.save()

                    # Now delete the types that are being collapsed
                    # But before we do that, delete the types from the group we are currently in and append the new group
                    # The types are automatically removed from the analysis by deleting the types
                    # as the info of which types are held in which grops are held in the types themselves
                    # Also the new type is already added to the analysis_group by the type attribute analysisGroupOf
                    del infoListDict[a]
                    del infoListDict[b]
                    a.delete()
                    b.delete()


                    # update the infoList
                    # set initial to false so that the infoList isn't generated from scratch
                    initial = False
                    newDictelement = getDefRefSeqIntraAbundanceForAllProfiles(analysisGrporAnType=newAnalysisType, isinitial=initial, orderedlistofdefiningrefseqs=orderedListOfDefiningRefSeqs)
                    infoListDict.update(newDictelement)

                    # Check to see that none of the CCs associted with this group have types that have defining intras in common
                    # For those CCs that do, delete the type representing the least seqs in that CC
                    # by deleting the CCtype and removing the CC from the type's cladeCollectionList
                    # if a type is removed then need to recalculate its infolist
                    # so the below method will output a list of type IDs that need the infolist regenerating
                    typesAffectedIds = removeSubsetTypes(listOfGroups[i], CCRefSeqAbundDict)
                    # regenerate the infoList element for each of the affected types
                    initial = False
                    for typeID in typesAffectedIds:
                        newDictelementSub = getDefRefSeqIntraAbundanceForAllProfiles(analysisGrporAnType=analysis_type.objects.get(id=typeID), isinitial=initial, orderedlistofdefiningrefseqs=orderedListOfDefiningRefSeqs)
                        infoListDict.update(newDictelementSub)

                    # Now break from the for loop that is comparing the types within the group we are in
                    # we haven't increased i at this point so we will analyse the same group again
                    # but this time with the new type and without the collapsed types
                    # if we make it through the type comparisons for a group and we don't find a collapse then
                    # we will make it to the i increase and we will move onto the next group
                    skip = True
                    break
                elif collapseStatus == 'CCCollapse':
                    # Here we are going to create a new type that contains the CCs that these two type have in common.
                    # Simple, but brilliant solution to co-occurence. The only issues I see occuring are what
                    # to do with the remaining two types, as one of them may have only very few CCs left in it
                    # we will also need to update the remaining two types (if they have CCs left in them)
                    # with the CCs for the new types removed. Hopefully we can use the current update type code for that
                    typeACCs = [int(x) for x in a.listOfCladeCollections.split(',')]
                    typeBCCs = [int(x) for x in b.listOfCladeCollections.split(',')]
                    newCCColection = list(set(typeACCs).intersection(typeBCCs))
                    print('Removing {} co-occurances CCs from {} and {}'.format(len(newCCColection), a, b))
                    allCCsList = set()
                    allCCsList.update(typeACCs)
                    allCCsList.update(typeBCCs)
                    # First create new type

                    # Work out the new footprint
                    # This will simply be a combination of the two type's footprints as we know they are all
                    # found in the same samples so must all comtains the same intras
                    consensusFootprint = set()
                    consensusFootprint.update(a.getOrderedFootprintList(), b.getOrderedFootprintList())
                    consensusFootprint = list(consensusFootprint)

                    # Work out MajRefSeqSet
                    newAnalysisType = analysis_type(dataAnalysisFrom=analysisObj, clade=a.clade, analysisGroupOf=listOfGroups[i])
                    newAnalysisType.listOfCladeCollections = ','.join([str(x) for x in newCCColection])
                    # NB orderedFootprintlist doesn't actually have to be ordered. That will
                    # happen when we update type attributes
                    newAnalysisType.orderedFootprintList = ','.join([str(refSeq.id) for refSeq in consensusFootprint])
                    newAnalysisType.save()
                    # This updateType will also set the type.coDom and the type.MajRefSeqSet
                    newAnalysisType.updateTypeAttributes()
                    newAnalysisType.save()

                    print('Creating type {}'.format(str(newAnalysisType)))

                    # Delete the cladeCollectionTypes that were associated with the old types
                    cctsToDelete = list(clade_collection_type.objects.filter(analysisTypeOf__in=[a, b], cladeCollectionFoundIn__in=clade_collection.objects.filter(id__in=newCCColection)))
                    for cct in cctsToDelete:
                        cct.delete()

                    # Create clade collection type for each CC in the new type
                    newCCs = clade_collection.objects.filter(id__in=newCCColection)
                    for CC in newCCs:
                        newcladeCollectionType = clade_collection_type(analysisTypeOf=newAnalysisType, cladeCollectionFoundIn=CC)
                        newcladeCollectionType.save()

                    # Delete the intraInfo for old types, we will remake this
                    del infoListDict[a]
                    del infoListDict[b]

                    # Remove the CCs of the new type from the old types
                    aDeleted = False
                    bDeleted = False
                    CCsRemainingInA = list(set(typeACCs).difference(newCCColection))
                    CCsRemainingInB = list(set(typeBCCs).difference(newCCColection))
                    if CCsRemainingInA:
                        a.listOfCladeCollections = ','.join([str(x) for x in CCsRemainingInA])
                        a.save()
                    else:
                        aDeleted = True
                        a.delete()
                    if CCsRemainingInB:
                        b.listOfCladeCollections = ','.join([str(x) for x in CCsRemainingInB])
                        b.save()
                    else:
                        bDeleted = True
                        b.delete()

                    # Update the types
                    if not aDeleted:
                        a.updateTypeAttributes()
                        a.save()
                    if not bDeleted:
                        b.updateTypeAttributes()
                        b.save()

                    # update the infoList
                    # set initial to false so that the infoList isn't generated from scratch
                    initial = False
                    for antype, deleted in zip([newAnalysisType, a, b], [False, aDeleted, bDeleted]):
                        if not deleted:
                            newDictelement = getDefRefSeqIntraAbundanceForAllProfiles(analysisGrporAnType=antype,
                                                                                   isinitial=initial,
                                                                                   orderedlistofdefiningrefseqs=orderedListOfDefiningRefSeqs)
                            infoListDict.update(newDictelement)

                    # Check to see that none of the CCs associted with this group have types that have defining intras in common
                    # For those CCs that do, delete the type representing the least seqs in that CC
                    # by deleting the CCtype and removing the CC from the type's cladeCollectionList
                    # if a type is removed then need to recalculate its infolist
                    # so the below method will output a list of type IDs that need the infolist regenerating
                    typesAffectedIds = removeSubsetTypes(listOfGroups[i], CCRefSeqAbundDict)
                    # regenerate the infoList element for each of the affected types
                    initial = False
                    for typeID in typesAffectedIds:
                        newDictelementSub = getDefRefSeqIntraAbundanceForAllProfiles(analysisGrporAnType=analysis_type.objects.get(id=typeID), isinitial=initial, orderedlistofdefiningrefseqs=orderedListOfDefiningRefSeqs)
                        infoListDict.update(newDictelementSub)


                    # Now break from the for loop that is comparing the types within the group we are in
                    # we haven't increased i at this point so we will analyse the same group again
                    # but this time with the intra infos for the old types and the new type
                    # if we make it through the type comparisons for a group and we don't find a collapse then
                    # we will make it to the i increase and we will move onto the next group

                    skip = True
                    break


        else:
            i += 1
            listNotToCollapse = []
            initial = True
            continue
        if not skip:
            listNotToCollapse = []
            initial = True
            i += 1



    analysisObj.analysisTypesCollapsed = True
    analysisObj.save()
    return

def createCCRefSeqAbundDict(querySetOfCCs):
    ############# INTRA ABUND DICT FOR EACH CC MP ################
    # Make dict of intraAbundances for each CC #TODO multiprocess this
    CCCol = clade_collection.objects.filter(id__in=querySetOfCCs)
    managerTwo = Manager()
    CCRefSeqAbundDict = managerTwo.dict()

    CCInputQueue = Queue()

    for CC in CCCol:
        CCInputQueue.put(CC)

    numProcessors = nProcessors
    for N in range(numProcessors):
        CCInputQueue.put('STOP')

    allProcesses = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()

    for N in range(numProcessors):
        p = Process(target=workerAssignmentTwo, args=(CCInputQueue, CCRefSeqAbundDict))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

    return CCRefSeqAbundDict
    ##############################################################

def checkForSharedIntraTypes(cctwo, ccrefseqabunddict, antypefootprintdict, groupid):
    # Keep a list of which types have been effected so that intraInfoList can be generated
    # this will be done by returning setOfTypeIdsAffected
    setOfTypeIdsAffected = set()
    # Work out which types are associted with the CC in Q
    CCinQ = clade_collection.objects.get(id=cctwo)
    listOfCCTypesOfCCTwo = clade_collection_type.objects.filter(cladeCollectionFoundIn=CCinQ,
                                                              analysisTypeOf__dataAnalysisFrom=analysisObj)

    # List of types in the CC that are of the group in Q
    groupInQ = analysis_group.objects.get(id=groupid)

    listOfTypesInCCTwo = list(analysis_type.objects.filter(clade_collection_type__cladeCollectionFoundIn=CCinQ, dataAnalysisFrom=analysisObj, analysisGroupOf=groupInQ))
    # listOfTypesInCCTwo = [cctype.analysisTypeOf for cctype in listOfCCTypesOfCCTwo]
    CCabundDict = ccrefseqabunddict[CCinQ]
    # If any two types share intras
    # Work out which of the types represents the larger number of seqs for the CC in Q
    typeToDel = []
    for a, b in itertools.combinations(listOfTypesInCCTwo, 2):

        # This was causing problems because some of e.g. when a type has D1/D4... and D5 and D5 was collapsed with D5/D4
        # Then the new D5/D4 type shared D4 with D1/D4... original type and had more seqs so binned it.
        # So instead I think the simplest action will work best here. Bin the shorter type
        aFootprint = antypefootprintdict[a.id]
        bFootprint = antypefootprintdict[b.id]
        if bool(aFootprint & bFootprint):
            aLength = len(a.getOrderedFootprintList())
            bLength = len(b.getOrderedFootprintList())

            # If the new type has bigger seq total then current candidate, del candidate from
            if aLength > bLength:
                # Then add id of b to typeToDel. Then done with the comparisons, del the CCtype and remove the CC from the type in Q
                typeToDel.append(b.id)
            else:
                typeToDel.append(a.id)
    # For each of the type in the del list, delete CCtype and remove CC from type
    for typeID in typeToDel:
        analysisTypeInQ = analysis_type.objects.get(id=typeID)
        listOfCCTypesOfCCTwo.get(analysisTypeOf=analysisTypeInQ).delete()
        analysisTypeInQ.removeCCFromCladeCollectionList(CCinQ)
        analysisTypeInQ.save()
        if len(analysisTypeInQ.listOfCladeCollections) < 2:
            # THen this type no longer has any CCs associated with it, so it won't have any
            # cctypes either. Therefore, delete.
            analysisTypeInQ.delete()
        else:
            analysisTypeInQ.updateTypeAttributes()
            analysisTypeInQ.save()
            setOfTypeIdsAffected.add(analysisTypeInQ.id)

    return setOfTypeIdsAffected

def removeSubsetTypes( group, ccrefseqabunddict,):
    # This will check to see that cladeCollections don't contain types that are perfect subsets of each other
    # This needs to be updated so that it makes sure that a given CC cannot contain two types that contain the same
    # an intra in common. If this is the case then the type which makes up for the largest number of sequences must be
    # kept.
    # There is code for this in the type assignment section. Bear in mind that if we find a type that needs to be
    # disassociated from a CC then we will need to recalculate the intra info for that type as well as the new type
    # created.

    # Keep track of which analysisTypes have been affected by this check by updating this set
    analysisTypeIdSet = set()




    listOfCladeCollectionsWithTypesFromThisGroupIDs = set()
    listOfTypes = analysis_type.objects.filter(analysisGroupOf = group)
    anTypeFootprintDict = {anType.id: frozenset(anType.getOrderedFootprintList()) for anType in listOfTypes}
    cladeCollectionTypes = clade_collection_type.objects.filter(analysisTypeOf__in = listOfTypes).distinct()
    for cctype in cladeCollectionTypes:
        listOfCladeCollectionsWithTypesFromThisGroupIDs.add(cctype.cladeCollectionFoundIn.id)



    for cctwo in listOfCladeCollectionsWithTypesFromThisGroupIDs:
        analysisTypeIdSet.update(checkForSharedIntraTypes(cctwo, ccrefseqabunddict, anTypeFootprintDict, group.id))





    return analysisTypeIdSet

def assignProfilesToGroups():
    '''

    Place each of the types into their respective groups
    Within groups assess each type pair for collapse
    Make note of which types collapse
    Where either of the types in a collapsing pair also collapse with alternative types but the other inital
    collapsing pair doesn't collapse with this then collapse the type with the maximum collapsing possibilities
    to the type with the closest Fst.
    Once a collapse has occured... lets get to here first and see what it looks like
    :return:
    '''

    # groupnames are the Maj types found in types
    # but if coDom types exist the groups contain all types that have majs that are found in the codoms
    # Get list of all types that have cladeCollections associated with them
    listOfAnalysisTypes = [anlType for anlType in analysis_type.objects.filter(dataAnalysisFrom=analysisObj) if anlType.listOfCladeCollections]



    # Create groupings of Types that have majs in common
    # groupNamesList contains one list per group. Each group list contains multiple types that share majs
    # This is basically creating a network of types that have Majs in common
    # e.g. C5, C3, C5/C3 would all be in the same group/. If C5/C3 were not present, then C5 and C3 would be seperate groups
    groupNamesList = []
    for anlType in listOfAnalysisTypes:
        # Assess whether to add typeinQ to a new group or to an already existent group
        if str(anlType) == 'C3a':
            apples = 'asdf'
        try:
            print('\rAssigning group to {}'.format(anlType), end='')
        except:
            apples = 'asdf'
        groupNamesList = createGroups(groupNamesList, anlType)

    # Initialize groups
    # for each group link the types it contains to the analysis_group instance
    for group in groupNamesList:
        print('Creating new group without name')
        newAnalysisGroup = analysis_group(dataAnalysisFrom = analysisObj)
        # Name each of the groups that were just created
        # newAnalysisGroup.generateName()
        newAnalysisGroup.save()
        for anlType in group: # for each of the analysisTypes in the group
            anlType.analysisGroupOf = newAnalysisGroup
            anlType.save()

    analysisObj.analysisTypesCollapsed = True
    analysisObj.save()

    return

def createGroups(groupNamesList, newAnlType):
    #Identify if there is a group that already contains Majs that are found in the anltype in question
    # If identified add the anlType to this grouping else make a new grouping

    groupList = groupNamesList
    indexMatch = set()
    for tempGroup in groupList: # for each grouping of anlTypes which is a list
        listOfMajRefSeqsFoundInTypesInTempGroup = []
        for aType in tempGroup:
            listOfMajRefSeqsFoundInTypesInTempGroup.extend(list(aType.getMajRefSeqSet()))
        for majRefSeq in newAnlType.getMajRefSeqSet():
            if majRefSeq in listOfMajRefSeqsFoundInTypesInTempGroup:
                # Then this type shares a majseq with one of the types in this group
                # Get index so that these two groups can be merged
                indexMatch.add(groupList.index(tempGroup))

    indexMatch = list(indexMatch)
    if indexMatch:
        if len(indexMatch) > 2:
            apples = 'asdf'
        # if we do have an index in the index Match then we have several types of one or more groups that
        # need to be made into a new group that include our new type
        # Get a list of all of the types that will be in the new group
        listOfTypesInGroups = []
        listOfGroups = [groupList[i] for i in indexMatch]
        for group in listOfGroups:
            for anType in group:
                listOfTypesInGroups.append(anType)

        # Also add our new type
        listOfTypesInGroups.append(newAnlType)
        # Now delete the groups that shared majs with our new type
        # Do this by creating a new  list without the sequences to delete
        newGroupNames = [groupList[index] for index in range(len(groupList)) if index not in indexMatch]
        # Append the new group of types to the newGroupNames and return
        # The list of types is a 2D list so must add each type from each list
        newGroupNames.append(listOfTypesInGroups)
        return newGroupNames

    else:
        # If at this point the indexMatch is empty them we have not found a group of types that share majs with our new
        # type. so, we make a new group for it and add it to the groupNamesList and return the list
        # in preparation for the next type to be assessed
        groupList.append([newAnlType])
        return groupList

    return

def getDefRefSeqIntraAbundanceForAllProfiles(analysisGrporAnType, isinitial, **kwargs):
    '''
    If isinital == True then this is the first run through of the profile collapse
    we will produce a dictionary where key is analysis_type and value is a list that contains the average relative
    abundance of each defining intra from all of the types found in the group
    when true, we must first generate a listOfDefiningRefSeqsAllTypes to work out which intras are found across
    all of the types in the gruop
    if isinitial == False, then we are sort of topping up the previously created dictionary with info for the types
    that have just been created. The output will only be a single type key to infolist value that will then be used
    to update the current infolistdict outside of this function
    This has been updated from the previous ratio system to the relabundance system
    '''

    if isinitial:
        listOfTypesInGroup = analysis_type.objects.filter(analysisGroupOf=analysisGrporAnType)
        # Get a list of all definingRefSeqs in the analysisTypes found in all of the types in the group
        definingRefSeqs = []
        for antype in listOfTypesInGroup:
            if str(antype) == '22094' or str(antype) == '22094/22523' or str(antype) == '22557':
                apples = 'asdf'
            definingRefSeqs += antype.getOrderedFootprintList()
        listOfDefiningRefSeqsAllTypes = list(set(definingRefSeqs))
        # Get all analysedSampleSequences that are of these referenceSequences and from the analysisTypes in question
        # https://docs.djangoproject.com/en/1.10/topics/db/queries/
        # 'Lookups that span relationships'
        # Turns out the the DJANGO documentation on this reverse span lookup is awful and wrong.
        # I figured it out using http://stackoverflow.com/questions/3831480/how-to-span-lookup-relationship-in-django
        # However probably best to use my example below as example

        anSmpSeqsOfTypes = data_set_sample_sequence.objects.filter(referenceSequenceOf__in=listOfDefiningRefSeqsAllTypes, cladeCollectionTwoFoundIn__clade_collection_type__analysisTypeOf__in=listOfTypesInGroup)

        bob = list(anSmpSeqsOfTypes)
        # for each of the refSeqs add up the abundances of the analysedSampleSequences that are of the refSeq type.
        # Keep these sums in a dict and use them to get the order of most abundance seqs.# This will be the oder that the seqs are compared in. e.g. ordered over all seqs for all types in the group
        refSeqCountDict = {}
        for refSeq in listOfDefiningRefSeqsAllTypes:
            refSeqCountDict[refSeq] = sum(anSmpSeq.abundance for anSmpSeq in anSmpSeqsOfTypes.filter(referenceSequenceOf=refSeq))
        # We should be careful here to check that the most abundant seq found in the orderedListOfDefiningRefSeqs is
        # found in all of the analysis type footprints. We could have a rare case where it is not


        #A dictionary that will hold an info list per type
        infoListDict = {}


        # Some of the managers can be created outside of this below loop
        managerOrderedListOfDefiningRefSeqs = Manager()
        orderedListOfDefiningRefSeqs = managerOrderedListOfDefiningRefSeqs.list([a[0] for a in sorted(refSeqCountDict.items(), key=operator.itemgetter(1), reverse=True)])

        for antype in listOfTypesInGroup:
            # Create the infolist

            print('calculating intra info for {}'.format(antype.name))

            #TODO multiprocess here
            # I think that this is the rate limiting step and therefore should be multiprocessed at this level
            # When we have the very slow info lists to make it is for the types which have a lot of clade collections
            # in them.
            # Input will be list of clade collections from the listOfCladeCollectionsOfType
            # Managed objects will be listOfCladeCollectionsOfType, orderedListOfDefiningRefSeqs, infoList
            # This is a lot of objects to be transfering to the processess, hopefully the overhead will not be too large

            # Create Managers (some managers are created outside this loop just above)
            managerListOfCladeCollectionsOfType = Manager()
            listOfCladeCollectionsOfType = managerListOfCladeCollectionsOfType.list(list(clade_collection.objects.filter(clade_collection_type__analysisTypeOf=antype)))

            managerInfoList = Manager()
            infoList = managerInfoList.list([[] for cc in listOfCladeCollectionsOfType])

            # managerNewInfoList = Manager()
            # newInfoList = managerNewInfoList.list([0 for refseq in orderedListOfDefiningRefSeqs])

            CCInputQueue = Queue()

            for CC in listOfCladeCollectionsOfType:
                CCInputQueue.put(CC)

            numProcessors = 6
            for N in range(numProcessors):
                CCInputQueue.put('STOP')

            allProcesses = []

            # close all connections to the db so that they are automatically recreated for each process
            # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
            db.connections.close_all()

            for N in range(numProcessors):
                p = Process(target=workerIntraAbundList, args=(CCInputQueue, orderedListOfDefiningRefSeqs, listOfCladeCollectionsOfType, infoList))
                allProcesses.append(p)
                p.start()

            for p in allProcesses:
                p.join()


            newInfoList = [0 for a in range(len(infoList[0]))]
            for I in range(len(infoList)):
                total = sum(infoList[I])
                for j in range(len(infoList[I])):# For each intra in order
                    newInfoList[j] += infoList[I][j] / total


            # create a dict where the anType is the key and the respective infoList is the value
            newInfoList = [x/len(infoList) for x in newInfoList]
            if sum(newInfoList) != 1:
                apples = 'asdf'
                print(str(sum(newInfoList)))

            infoListDict[antype] = newInfoList
            apples = 'pears'
        return infoListDict, orderedListOfDefiningRefSeqs
    else:
        # If here args passed in =
        # analysisGrporAnType = newAnalysisType,
        # isinitial = initial,
        # orderedlistofdefiningrefseqs = orderedListOfDefiningRefSeqs

        orderedListOfDefiningRefSeqs = kwargs.get('orderedlistofdefiningrefseqs', None)


        print('calculating intra info for {}'.format(analysisGrporAnType.name))



        managerListOfCladeCollectionsOfType = Manager()
        listOfCladeCollectionsOfType = managerListOfCladeCollectionsOfType.list(
            list(clade_collection.objects.filter(clade_collection_type__analysisTypeOf=analysisGrporAnType)))

        managerInfoList = Manager()
        infoList = managerInfoList.list([[] for cc in listOfCladeCollectionsOfType])

        CCInputQueue = Queue()

        for CC in listOfCladeCollectionsOfType:
            CCInputQueue.put(CC)

        numProcessors = 10
        for N in range(numProcessors):
            CCInputQueue.put('STOP')

        allProcesses = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        for N in range(numProcessors):
            p = Process(target=workerIntraAbundList,
                        args=(CCInputQueue, orderedListOfDefiningRefSeqs, listOfCladeCollectionsOfType, infoList))
            allProcesses.append(p)
            p.start()

        for p in allProcesses:
            p.join()

        newInfoList = [0 for a in range(len(infoList[0]))]
        for I in range(len(infoList)):
            total = sum(infoList[I])
            for j in range(len(infoList[I])):  # For each intra in order
                newInfoList[j] += infoList[I][j] / total

        # create a dict where the anType is the key and the respective infoList is the value
        newInfoList = [x / len(infoList) for x in newInfoList]


        newdictelement = {analysisGrporAnType: newInfoList}
        return newdictelement

    return

def workerIntraAbundList(CCInputQueue, orderedListOfDefiningRefSeqs, listOfCladeCollectionsOfType, infoList):

    for CC in iter(CCInputQueue.get, 'STOP'):
        print('{} is creating intraAbundList for {}'.format(current_process().name, CC))
        CCindex = listOfCladeCollectionsOfType.index(CC)

        infoListInQ = infoList[CCindex]
        for k in range(len(orderedListOfDefiningRefSeqs)):  # For each refSeq
            #TODO this is currently causing problems. The sqlite db is not handling concurrency well
            # in some instances we are not adding anything to the list
            try:
                #https://stackoverflow.com/questions/29392650/python-using-a-multidimensional-multiprocessing-manager-list

                infoListInQ.append(data_set_sample_sequence.objects.get(cladeCollectionTwoFoundIn=CC, referenceSequenceOf=orderedListOfDefiningRefSeqs[k]).abundance)

            except data_set_sample_sequence.DoesNotExist:
                infoListInQ.append(0)

        infoList[CCindex] = infoListInQ

def assessProfileCollapse(typeOne, typeTwo, infoListDict, orderOfRefSeqsInInfoListDictItems):
    # TODO collapsing types that are only one intra in length is continually causing us problems.
    # Given that we don't gain a lot by collapsing them, and that having them as singles makes it easy to identify when
    # types haven't been well called. And if we often we start to create run away types that get bigger and bigger
    # we are going to ignore collapse criteria if a type is only 1 intra long.
    #firstly check to see if the types co-occur significantly

    # If either of the types is only 1 in length
    if len(typeOne.getOrderedFootprintList()) == 1 or len(typeTwo.getOrderedFootprintList()) == 1:
        return 'noCollapse'

    if CCsInCommon(typeOne, typeTwo):
        return 'CCCollapse'

    #Then check to see if the ITS2 profiles of each of the types are significantly differnt
    if assessProfilesHaveSigDiff(typeOne=typeOne, typeTwo=typeTwo, infoListDict=infoListDict, orderOfRefSeqsInInfoListDictItems=orderOfRefSeqsInInfoListDictItems):
        return 'typeProfileCollapse'

    return 'noCollapse'

def CCsInCommon(typea, typeb):
    # get list of CC for each type
    typeACCs = [int(a) for a in typea.listOfCladeCollections.split(',')]
    typeBCCs = [int(a) for a in typeb.listOfCladeCollections.split(',')]
    CCsInCommonList = list(set(typeACCs).intersection(typeBCCs))
    allCCsList = set()
    allCCsList.update(typeACCs)
    allCCsList.update(typeBCCs)

    if len(CCsInCommonList) > 5:
        print('{} has {} CCs'.format(typea, len(typeACCs)))
        print('{} has {} CCs'.format(typeb, len(typeBCCs)))
        print('{} in common CCs'.format(len(CCsInCommonList)))
        return True
    if len(CCsInCommonList)/len(allCCsList) > 0.5:
        print('{} has {} CCs'.format(typea, len(typeACCs)))
        print('{} has {} CCs'.format(typeb, len(typeBCCs)))
        print('{} in common CCs'.format(len(CCsInCommonList)))
        return True
    return False

def assessProfilesHaveSigDiff(typeOne, typeTwo, infoListDict, orderOfRefSeqsInInfoListDictItems):

    ### The below essentially checks to see that the two types are sufficiently unique in their definitions.
    # i.e. we only want to have types that have genuine difference is their defining ITS2 seqs
    # A separate method will check to see if the two types co-occur in a significant number of sequencs.
    # These types will also be collapsed
    # assess each of the average relative abundances for each intra between the two types. If we can't find an intra pair that
    # has a ratio smaller than 0.5 then this is grounds for collapse. I.e. we have to demonstrate a significant
    # difference in intra abundances in order to not collapse.




    # Work out the set of intras that are found in both of the types
    setOfIntras = typeOne.getOrderedFootprintList()
    setOfIntras.extend(typeTwo.getOrderedFootprintList())
    setOfIntras = list(set(setOfIntras))

    # Work through each of the intras doing the comparisons
    # For each intra check to see if the average relative abundace of that intra in the type is within one half of the other
    # If all of the ratios are within the 0.5 then we will return true and collapse, else false and no collapse
    for i in range(len(setOfIntras)):
        numerator = min(infoListDict[typeOne][orderOfRefSeqsInInfoListDictItems.index(setOfIntras[i])], infoListDict[typeTwo][orderOfRefSeqsInInfoListDictItems.index(setOfIntras[i])])
        denominator = max(infoListDict[typeOne][orderOfRefSeqsInInfoListDictItems.index(setOfIntras[i])], infoListDict[typeTwo][orderOfRefSeqsInInfoListDictItems.index(setOfIntras[i])])
        #### We are going to shake things up a bit. We are going to ensist that
        # one of the numerator or denominator must be above 5%, else we won't count
        # the division towards not collapsing
        if numerator > 0.03 or denominator > 0.03:
            if denominator != 0:
                if numerator / denominator < 0.5:
                    return False

    return True

def createNewRefSeqName(closestMatch, listofseqnamesthatalreadyexist, listofrefseqdbnames):
    # This code happens when we have a seq that needs a name
    # We know the seq is not an exact match so we make a derivative name
    # closestMatch[0] = closestMatch name
    # closestMatch[1] = match %
    # closetsMatch[2] = coverage
    matchObj = re.match("^[A-I]{1}[0-9]{1,3}", closestMatch)
    baseName = matchObj.group(0)


    # If we have got here then it is time to derive a name
    # This will be the base name plus the next available alpha concatenation
    # This is somewhat ugly code but goes through first one set of 26 letters
    # Then two sets and then three, i.e. aab aac aad ... aba abb abc etc.
    alphaList = string.ascii_letters[0:26]
    for alpha in alphaList:
        if baseName + alpha not in listofseqnamesthatalreadyexist and baseName + alpha not in listofrefseqdbnames:
            # TODO we also need to check if the new name is already in the reference database.
            # if it is we can't associate to this name as later on we might find a seq that does
            # have a 100% coverage and match ot that name
            return baseName + alpha
    for alpha in alphaList:
        for alphaTwo in alphaList:
            if baseName + alpha + alphaTwo not in listofseqnamesthatalreadyexist and baseName + alpha + alphaTwo not in listofrefseqdbnames:
                return baseName + alpha + alphaTwo
    for alpha in alphaList:
        for alphaTwo in alphaList:
            for alphaThree in alphaList:
                if baseName + alpha + alphaTwo + alphaThree not in listofseqnamesthatalreadyexist and baseName + alpha + alphaTwo + alphaThree not in listofrefseqdbnames:
                    return baseName + alpha + alphaTwo + alphaThree

    return False
###########################





###### OUTPUT FORMAT FUNCTIONS ######
def assignSpecies():

    at = analysis_type.objects.filter(dataAnalysisFrom=analysisObj)
    IDs = [att.id for att in at]

    # For each analysis type check and assign the associated species
    for m in range(len(IDs)):
        assignedSpecies = []
        att = analysis_type.objects.get(id=IDs[m])
        # This is a 2D list of sequence abundances.
        # Each list is the seq abundances found in a clade collection
        # Each value within a list is the abundance of the defining sequences (in order of att.orderedFootprintList)
        refSeqOrder = att.orderedFootprintList
        refSeqIDsList = [int(a) for a in refSeqOrder.split(',')]
        footprintAbunInfo = json.loads(att.footprintSeqAbundances)
        clade = att.clade
        # A list that will contain the average abundance of the def ref seq in question
        averageAbundanceList = [0 for refSeq in refSeqIDsList]
        for j in range(len(footprintAbunInfo)):
            for i in range(len(refSeqIDsList)):
                averageAbundanceList[i] += int(footprintAbunInfo[j][i])/sum([int(a) for a in footprintAbunInfo[j]])
        # Here we have the abundances summated for each seq in the averageAbundanceList
        # Just have to divide by the num cladeCollections to get averages
        for i in range(len(averageAbundanceList)):
            averageAbundanceList[i] = averageAbundanceList[i]/len(footprintAbunInfo)
        # Now we have the averge abundance proportions we make dict to the refseq name
        listOfRefSeqs = []

        for i in range(len(refSeqIDsList)):
            listOfRefSeqs.append(reference_sequence.objects.get(id=refSeqIDsList[i]))
        listOfAllDefSeqNames = [refseq.name for refseq in listOfRefSeqs]
        proportionDict = {L: K for L, K in zip(listOfAllDefSeqNames, averageAbundanceList)}
        majSeqs = [refSeq.name for refSeq in reference_sequence.objects.filter(id__in=att.MajRefSeqSet.split(','))]
        # Now we have all of the info we need to work out which species to associate to the type
        if clade == 'A':
            if 'A1' in majSeqs:
                assignedSpecies.append('S. microadriaticum')
            if 'A2' in majSeqs:
                assignedSpecies.append('S. pilosum')
            if 'A3' in majSeqs:
                assignedSpecies.extend(['S. natans', 'S. tridacnidorum'])
            if 'A4' in majSeqs:
                assignedSpecies.append('S. linucheae')
        elif clade == 'B':
            if 'B1' in majSeqs:
                assignedSpecies.extend(['S. minutum', 'S. antillogorgium', 'S. pseudominutum'])
            if 'B2' in majSeqs:
                assignedSpecies.append('S. psygmophilum')
            if 'B4' in majSeqs:
                assignedSpecies.append('S. muscatinei')
            if 'B7' in majSeqs or 'B13' in majSeqs:
                assignedSpecies.append('S. endomadracis')
            if 'B2a' in majSeqs:
                assignedSpecies.append('S. aenigmaticum')
        elif clade == 'C':
            if 'C1' in majSeqs:
                assignedSpecies.append('S. goreaui')
            if 'C3' in listOfAllDefSeqNames and 'C3gulf' in listOfAllDefSeqNames:
                assignedSpecies.append('S. thermophilum')
        elif clade == 'D':
            # I have decided that we are not going to take into account the abundance of non-maj intragenomic
            # defining sequences. e.g. D6 when calling associated species.
            # This is because there can be a large difference sample to sample in the abundance of the sequences
            # Rather we will assign both clade D species if the required sequences are present
            # We are also giving the researcher the average abundances and SDs for each output type
            if 'D1' in majSeqs:
                if 'D4' not in listOfAllDefSeqNames:
                    assignedSpecies.append('S. glynnii')
                else:# There is a significant abundance of D4
                    if 'D6' in listOfAllDefSeqNames:
                        # Then there is a significant amount of D6
                        assignedSpecies.extend(['S. glynnii', 'S. trenchii'])
                    else:
                        # THere is D1, D4 but not D6
                        assignedSpecies.append('S. trenchii')
            if 'D8' in majSeqs or 'D12' in majSeqs or 'D13' in majSeqs:
                assignedSpecies.append('S. eurythalpos')
            if 'D15' in majSeqs:
                assignedSpecies.append('S. boreum')
        elif clade == 'E':
            assignedSpecies.append('S. voratum')
        elif clade == 'F':
            if 'F1' in majSeqs:
                assignedSpecies.append('S. kawagutii')
        elif clade == 'G':
            pass
        elif clade == 'H':
            pass
        elif clade == 'I':
            pass

        if not assignedSpecies:#If no suggested species have been associated
            att.species = 'None'
            att.save()
        else:
            att.species = ','.join(assignedSpecies)
            att.save()

    analysisObj.speciesAssociated = True
    analysisObj.save()
    return

def sort_list_of_types_by_clade_collections_in_current_output(list_of_analysis_types, clade_collection_types_in_current_output):
    # create a list of tupples that is the type ID and the number of this outputs CCs that it was found in.

    tuple_list = []
    clade_collections_in_current_output = [clade_collection_type_obj.cladeCollectionFoundIn.id for clade_collection_type_obj in clade_collection_types_in_current_output]
    for at in list_of_analysis_types:
        list_of_clade_collections_found_in = [int(x) for x in at.listOfCladeCollections.split(',')]
        num_clade_collections_of_output = list(set(clade_collections_in_current_output).intersection(list_of_clade_collections_found_in))
        tuple_list.append((at.id, len(num_clade_collections_of_output)))

    type_IDs_sorted = sorted(tuple_list, key=lambda x: x[1], reverse=True)

    return [analysis_type.objects.get(id=x[0]) for x in type_IDs_sorted]

def formatOutput(analysisobj, numProcessors, datasubstooutput):

    analysisObj = analysisobj
    # This is one of the last things to do before we can use our first dataset
    # The table will have types as columns and rows as samples
    # its rows will give various data about each type as well as how much they were present in each sample

    # It will produce four output files, for DIV abundances and proportions and Type abundances and proportions.
    # found in the given clade collection.
    # Types will be listed first by clade and then by the number of clade collections the types were found in
    # The species groups will be numbered consecutively as they appear
    # Each type's ID will also be given as a UID.
    # The date the database was accessed and the version should also be noted
    # The formal species descriptions which correspond to the found ITS2 type will be noted for each type
    # Finally the AccessionNumber of each of the defining reference species will also be noted
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    '''
    32: VoolstraSet
    40: EmilyData
    47: EdData
    51: ed_24_07_17_singapore
    52: UTS_ONE_25_07_17
    53: UTS_TWO_25_07_17
    54: 20171008_jitErn
    55: 20171009_roberto
    56: 20171016_alejandro
    57: 20171016_Pappas
    58: 20171017_roberto
    59: 20171017_tullia
    60: 20170719_alejandro
    61: 20170713_alejandro
    62: 20170506_roberto
    63: ryan_first
    '''


    dataSubmissionsToOutput = [int(a) for a in datasubstooutput.split(',')]

    # Get collection of types that are specific for the dataSubmissions we are looking at
    querySetOfDataSubmissions = data_set.objects.filter(id__in=dataSubmissionsToOutput)

    cladeCollectionsFromDSs = clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom__in=querySetOfDataSubmissions)
    cladeCollectionTypesFromThisAnalysisObjAndDataSubmission = clade_collection_type.objects.filter(cladeCollectionFoundIn__in=cladeCollectionsFromDSs, analysisTypeOf__dataAnalysisFrom=analysisObj).distinct()

    # cladeCollectionsFromDSs = list(clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom__in=dataSubmissionsToOutput))
    # cladeCollectionTypesFromThisAnalysisObjAndDataSubmission = list(clade_collection_type.objects.filter(cladeCollectionFoundIn__in=cladeCollectionsFromDSs, analysisTypeOf__dataAnalysisFrom=analysisObj).distinct())

    at = set()
    for cct in cladeCollectionTypesFromThisAnalysisObjAndDataSubmission:
        at.add(cct.analysisTypeOf)
    at = list(at)

    # at = analysis_type.objects.filter(dataAnalysisFrom=analysisObj)
    IDs = [att.id for att in at]
    outputTableOne = []
    outputTableTwo = []
    # Need to get a list of Samples that are part of the dataAnalysis
    #listOfDataSetSamples = list(data_set_sample.objects.filter(dataSubmissionFrom__in=data_set.objects.filter(id__in=[int(a) for a in analysisObj.listOfDataSubmissions.split(',')])))
    listOfDataSetSamples = list(data_set_sample.objects.filter(dataSubmissionFrom__in=data_set.objects.filter(id__in=dataSubmissionsToOutput)))
    # The header row for the table
    # outputTableOne.append('ITS2 type profile UID\tClade\tSpeciesGrp\tType_Profile\tMajority ITS2 Seq\tAssociated Species\tType abund local\tType abund DB\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))
    # outputTableTwo.append('ITS2 type profile UID\tClade\tSpeciesGrp\tType_Profile\tMajority ITS2 Seq\tAssociated Species\tType abund local\tType abund DB\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))

    # I am making some changes to the output format for the SP manuscript
    outputTableOne.append('ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))
    outputTableTwo.append('ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))

    # Now go through the types, firstly by clade and then by the number of cladeCollections they were found in
    # Populate a 2D list of types with a list per clade
    typesCladalList = [[] for clade in cladeList]
    for att in at:
        try:
            if len(att.listOfCladeCollections) > 0:
                typesCladalList[cladeList.index(att.clade)].append(att)
        except:
            pass


    # Go clade by clade
    for i in range(len(typesCladalList)):
        if typesCladalList[i]:

            ###### CALCULATE REL ABUND AND SD OF DEF INTRAS FOR THIS TYPE ###############
            #     # For each type we want to calculate the average proportions of the defining seqs in that type
            #     # We will also calculate an SD.
            #     # We will do both of these calculations using the footprintSeqAbundances, listOfCladeCollections
            #     # and orderedFoorprintList attributes of each type

            ########## MAKE GROUP COUNTER AND DICT ###########
            # These are now going to be managed items for use in the MP
            # We want to name the groups that types are found in sequencially
            # To get the next number to name a group we will use the groupCount
            # To look up what number has been assigned to groups that have already been printed we will use the groupDict
            # groupCount = 0
            # groupDict = {}
            ###################################################

            # This is currently being sorted according to the number of listOfCladeCollections the type was found in
            # which translates to the number of clade collections across the whole data analysis, i.e. samples from this
            #  out put also from any other samples inlcuded in the dataanalysis we want to try to output according to
            # the number of samples of the output it was found in instead
            # get a sorted list of the types by most commonly found first
            sortedListOfTypes = sorted(typesCladalList[i], key=lambda x: len(x.listOfCladeCollections.split(',')), reverse=True)

            # We need to produce a list of analysisTypes according to the number of cladecollections of the output
            # they are found in. This shouldn't be too tricky.
            # We will pass in typesCladalList[i] and we will simply look to see which of the .listOfCladeCollections
            # of each type are found in the cladeCollectionTypesFromThisAnalysisObjAndDataSubmission

            sortedListOfTypes = sort_list_of_types_by_clade_collections_in_current_output(typesCladalList[i], cladeCollectionTypesFromThisAnalysisObjAndDataSubmission)

            # Here we will MP at the type level
            # i.e. each type will be processed on a differnt core. In order for this to work we should have a managed
            # dictionary where the key can be the types ID and the value can be the datalist
            # In order to get the types output in the correct order we should use the sortedListOfTypes to resort the
            # data once the MPing has been done.

            # This dict will hold all of the output rows that the MPing has created.
            typeOutputManager = Manager()
            typeOutputManagedDict = typeOutputManager.dict({type: None for type in sortedListOfTypes})
            # This dict will hold details of which groups have already been written
            groupDictManager = Manager()
            groupDictManagedDict = groupDictManager.dict()

            groupCountManager = Manager()
            groupCountManagedInt = groupCountManager.Value(c_int, 0)

            listOfDataSetSamplesManager = Manager()
            listOfDataSetSamplesManagedList = listOfDataSetSamplesManager.list(listOfDataSetSamples)

            typeInputQueue = Queue()

            for anType in sortedListOfTypes:
                typeInputQueue.put(anType)


            for N in range(numProcessors):
                typeInputQueue.put('STOP')

            allProcesses = []

            # close all connections to the db so that they are automatically recreated for each process
            # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
            db.connections.close_all()

            for N in range(numProcessors):
                p = Process(target=outputWorkerOne, args=(typeInputQueue, groupDictManagedDict, groupCountManagedInt, listOfDataSetSamplesManagedList, cladeList[i], typeOutputManagedDict))
                allProcesses.append(p)
                p.start()

            for p in allProcesses:
                p.join()


            # # now go type by type
            # for j in range(len(sortedListOfTypes)):
            #
            #     ###### CALCULATE REL ABUND AND SD OF DEF INTRAS FOR THIS TYPE ###############
            #     # For each type we want to calculate the average proportions of the defining seqs in that type
            #     # We will also calculate an SD.
            #     # We will do both of these calculations using the footprintSeqAbundances, listOfCladeCollections
            #     # and orderedFoorprintList attributes of each type
            #     footprintAbundances = json.loads(sortedListOfTypes[j].footprintSeqAbundances)
            #     orderedFootPrintlist = sortedListOfTypes[j].orderedFootprintList.split(',')
            #     # We have to make a decision as ot whether this average should represent all the findings of this type
            #     # or whether we should only represent the averages of samples found in this dataset.
            #     # I think it needs to be a global average. Because the type is defined based on all samples in the
            #     # SymPortal db.
            #     averageList = [[] for defseq in orderedFootPrintlist]
            #     for a in range(len(footprintAbundances)):
            #         for b in range(len(orderedFootPrintlist)):
            #             averageList[b].append(float(int(footprintAbundances[a][b])/sum(footprintAbundances[a])))
            #     # Here we have a list for each def seq (as proportions), this can now have totals, SDs and averages calculated on it
            #     # then convert to proportions by dividing by total
            #     SDList = []
            #     TotList = []
            #
            #     for a in range(len(averageList)):
            #         SDList.append(statistics.pstdev(averageList[a]))
            #         TotList.append(sum(averageList[a]))
            #     tot = sum(TotList)
            #     for a in range(len(TotList)):
            #         TotList[a] = TotList[a]/tot
            #
            #     # The TotList now contains the proportions of each def seq
            #     # The SDList now contains the SDs for each of the def seqs proportions
            #     ###########################################################################
            #
            #
            #
            #     # Counter that will increment with every sample type is found in
            #     # This is counter will end up being the type abund local value
            #     abundanceCount = 0
            #
            #     typeInQ = sortedListOfTypes[j]
            #     clade = typeInQ.clade
            #
            #     ###### GET GROUP NAME INFO ########
            #     groupID = typeInQ.analysisGroupOf.id
            #     if groupID not in groupDict.keys():
            #         groupCount += 1
            #         groupDict[groupID] = groupCount
            #     ###################################
            #
            #
            #     # For each type create a holder that will hold 0s for each sample until populated below
            #     dataRowRaw = [0 for smpl in listOfDataSetSamples]
            #     dataRowProp = [0 for smpl in listOfDataSetSamples]
            #
            #
            #     typeCCIDs = [int(a) for a in typeInQ.listOfCladeCollections.split(',')]
            #
            #
            #     # We also need the type abund db value. We can get this from the type cladeCollections
            #     globalCount = len(typeCCIDs)
            #
            #
            #
            #     for k in range(len(listOfDataSetSamples)):# Within each type go through each of the samples
            #         print('Outputing data for\ttype: {} \tsample: {}'.format(sortedListOfTypes[j], listOfDataSetSamples[k]))
            #         # when looking at sample, get the clade collections found in the sample
            #         CCsInSample = clade_collection.objects.filter(dataSetSampleFrom = listOfDataSetSamples[k])
            #         # because we are within a clade here the type can only have one cladecollection in common with the sample
            #         # so check to see if there are any items in common between the type's CC list and the samples CC
            #         ccInSampleIDs = [a.id for a in CCsInSample]
            #
            #         IDsInCommon = list(set(ccInSampleIDs).intersection(typeCCIDs))
            #         if IDsInCommon:# Then this type was found in this sample
            #             abundanceCount += 1
            #             # Need to work out how many seqs were found from the sample for this type
            #             # Also need to work this out as a proportion of all of the Symbiodinium seqs found in the sample
            #             ccInQ = clade_collection.objects.get(id=IDsInCommon[0])
            #             # The number of sequences that were returned for the sample in question
            #             totSeqNum = data_set_sample.objects.get(clade_collection=ccInQ).finalTotSeqNum
            #             # The number of sequences that make up the type in q
            #             # The index of the CC in the type's listOfCladeCollections
            #             CCindexInType = typeCCIDs.index(ccInQ.id)
            #             # List containing the abundances of each of the ref seqs that make up the type in the given CC
            #             seqAbundanceInfoForCCandTypeInQ = json.loads(typeInQ.footprintSeqAbundances)[CCindexInType]
            #             # total of the above list, i.e. seqs in type
            #             sumOfDefRefSeqAbundInType = sum(seqAbundanceInfoForCCandTypeInQ)
            #             # type abundance as proportion of the total seqs found in the sample
            #             typeProp = sumOfDefRefSeqAbundInType/totSeqNum
            #             # Now populate the dataRow with the sumOfDefRefSeqAbundInType and the typeProp
            #             dataRowRaw[k] = sumOfDefRefSeqAbundInType
            #             dataRowProp[k] = typeProp
            #
            #     # Type Profile
            #     typeID = sortedListOfTypes[j].id
            #     speciesGrp = '{}_Grp_{}'.format(cladeList[i], groupDict[groupID])
            #     typeProfile = typeInQ.name
            #     species = typeInQ.species
            #     typeAbundance = abundanceCount
            #     # This is currently putting out the Majs in an odd order due to the fact that the MajRefSeqSet is
            #     # a set. Instead we should get the Majs from the name.
            #     # majITS2 = '/'.join([str(refseq) for refseq in sortedListOfTypes[j].getMajRefSeqSet()])
            #     majITS2, majList = getMajList(sortedListOfTypes[j])
            #     typeAbundAndSDString = getAbundStr(TotList, SDList, majList)
            #
            #     sequenceAccession = typeInQ.generateName(accession=True)
            #     rowOne = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(typeID, clade, speciesGrp, typeProfile, majITS2,  species, str(typeAbundance), str(globalCount), '\t'.join([str(a) for a in dataRowRaw]), sequenceAccession, typeAbundAndSDString)
            #     rowTwo = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(typeID, clade, speciesGrp, typeProfile, majITS2,  species, str(typeAbundance), str(globalCount), '\t'.join(["{:.3f}".format(prop) for prop in dataRowProp]), sequenceAccession, typeAbundAndSDString)
            #     outputTableOne.append(rowOne)
            #     outputTableTwo.append(rowTwo)

            # Here the MP is completed and we now need to add the lists that are held in the managed dict to
            # the outputTableOne and outputTableTwo

            for antype in sortedListOfTypes:
                outputTableOne.append(typeOutputManagedDict[antype][0])
                outputTableTwo.append(typeOutputManagedDict[antype][1])

    # At this point we should have gone through each of the clades in order and for each clade we should have
    # created an analysis_type row in order of the sortedListOfTypes.

    # The current out put has samples on the Y and types on the X.
    # We want to reverse this as there will often be more samples than types
    # https://stackoverflow.com/questions/6473679/transpose-list-of-lists

    # Put into list of lists format to transpose
    newOutOne = [a.split('\t') for a in outputTableOne]
    newOutTwo = [a.split('\t') for a in outputTableTwo]

    # Transpose
    outputTableOne = list(map(list, zip(*newOutOne)))
    outputTableTwo = list(map(list, zip(*newOutTwo)))

    # Put back into tab delim format
    outputTableOne = ['\t'.join(a) for a in outputTableOne]
    outputTableTwo = ['\t'.join(a) for a in outputTableTwo]

    outputDir = os.path.join(os.path.dirname(__file__), 'outputs/{}'.format(analysisObj.id))
    os.makedirs(outputDir, exist_ok=True)
    os.chdir(outputDir)

    writeListToDestination('{}/{}_{}.profiles.counts.txt'.format(outputDir, analysisObj.id, '_'.join([str(a) for a in dataSubmissionsToOutput])), outputTableOne)
    writeListToDestination('{}/{}_{}.profiles.props.txt'.format(outputDir, analysisObj.id, '_'.join([str(a) for a in dataSubmissionsToOutput])), outputTableTwo)

    ########################## ITS2 INTRA ABUND COUNT TABLE ################################
    #TODO this is where we're going to have to work with the sequences that aren't part of a type.
    # Esentially we want to end up with the noName sequences divieded up cladally.
    # So at the moment this will mean that we will divide up the current no names but also
    # add information about the other cladal sequences that didn't make it into a cladeCollection


    # Now we need to output the ITS2 intragenomic abundances that the types are based on for each CC
    # Already have querySetOfDataSubmissions and cladeCollectionsFromDSs

    ################ GET ORDERED LIST OF INTRAS BY CLADE THEN ABUNDANCE #################
    # This is the master list that will hold the eventually ordered defining referenceSequences
    cladeAbundanceOrderedRefSeqList = []
    refSeqsInDSs = reference_sequence.objects.filter(data_set_sample_sequence__cladeCollectionTwoFoundIn__in = cladeCollectionsFromDSs).distinct()

    # Get list of clades found
    cladeSet = set()
    for refSeq in refSeqsInDSs:
        cladeSet.add(refSeq.clade)
    # Order the list of clades alphabetically
    cladeList = [a for a in cladeList if a in cladeSet]


    for i in range(len(cladeList)):
        intraCounter = defaultdict(int)

        # refSeqs of clade
        refSeqsOfClade = refSeqsInDSs.filter(clade=cladeList[i])

        # dataSetSampleSequencesOfClade
        dataSetSampleSequencesOfClade = data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn__in = cladeCollectionsFromDSs, referenceSequenceOf__in = refSeqsOfClade)
        #TODO This takes a bit of time and could in theory be multiprocessed
        # Popualate the counter dict
        # TODO MP here ############# DEBUG
        intraCounterManager = Manager()
        intraCounter = intraCounterManager.dict()

        dsssQueue = Queue()

        for dsss in dataSetSampleSequencesOfClade:
            dsssQueue.put(dsss)

        for N in range(numProcessors):
            dsssQueue.put('STOP')

        allProcesses = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        print('Creating refSeq counts for clade {}'.format(cladeList[i]))
        for N in range(numProcessors):
            p = Process(target=outputWorkerTwo, args=(dsssQueue, intraCounter))
            allProcesses.append(p)
            p.start()

        for p in allProcesses:
            p.join()

        # for dsss in dataSetSampleSequencesOfClade:
        #     intraCounter[dsss.referenceSequenceOf.name] += dsss.abundance

        # Sort counter dict by abundance
        sortedIntraCounter = sorted(intraCounter.items(), key=operator.itemgetter(1), reverse=True)

        # Add the refSeqNames to the master cladeAbundanceOrderedRefSeqList
        # Do not add the noName items. These will be summed across the clades and entered into the list
        # once their final abundances have been calculated
        cladeAbundanceOrderedRefSeqList.extend([a[0] for a in sortedIntraCounter if a[0] != 'noName'])

    # At this point we have the cladeAbundanceOrderedRefSeqList
    # Now we just need to add the noName refSeq to the list

    # cladeAbundanceOrderedRefSeqList.append('noName')


    # Create the two output tables as lists
    intraAbundCountTable = []
    intraAbundPropTable = []
    # Create the header line one for prop one for count
    # add the cladal noName categories directly to the header
    # do not add them to the cladeAbundanceOrderedRefSeqList so that we don't append 0s to them
    # we will add the abundances in the outputWorkerThree method below

    headerPre = 'Samples\t{}'.format('\t'.join(cladeAbundanceOrderedRefSeqList))
    noNameStrings = '\t'.join(['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']])
    header = '\t'.join([headerPre, noNameStrings])


    intraAbundCountTable.append(header)
    intraAbundPropTable.append(header)
    ######################################################################################

    ############ POPULATE TABLE WITH CC DATA #############
    # For every sample
    sampleList = data_set_sample.objects.filter(dataSubmissionFrom__in=querySetOfDataSubmissions)


    # In order to MP this we will have to pay attention to the order. As we can't keep the order as we work with the
    # MPs we will do as we did above for the profies outputs and have an output dict that will have the sample as key
    # and its corresponding data row as the value. Then when we have finished the MPing we will go through the
    # sampleList order and output the data in this order.
    # So we will need a managed output dict.

    sampleOutputDictManager = Manager()
    managedSampleOutputDict = sampleOutputDictManager.dict()

    cladeAbundanceOrderedRefSeqListManager = Manager()
    cladeAbundanceOrderedRefSeqList = cladeAbundanceOrderedRefSeqListManager.list(cladeAbundanceOrderedRefSeqList)



    dssQueue = Queue()

    for dss in sampleList:
        dssQueue.put(dss)

    for N in range(numProcessors):
        dssQueue.put('STOP')

    allProcesses = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()

    for N in range(numProcessors):
        p = Process(target=outputWorkerThree, args=(dssQueue, managedSampleOutputDict, cladeAbundanceOrderedRefSeqList, querySetOfDataSubmissions))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

    # for i in range(len(sampleList)):
    #     print('Outputting DIV data for {}'.format(sampleList[i].name))
    #     # List that will hold the row
    #     sampleRowDataCounts = []
    #     sampleRowDataProps = []
    #
    #     # CCs of sample
    #     sampleCCs = cladeCollectionsFromDSs.filter(dataSetSampleFrom = sampleList[i])
    #
    #     # Get number of seqs found in Sample
    #     dsssInSample = data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn__in=sampleCCs)
    #     sampleSeqTot = sum([dsss.abundance for dsss in dsssInSample])
    #
    #     # create dataSetSampleSeqDict for sample
    #     dataSetSampleSeqDict = defaultdict(int)
    #     for dsss in dsssInSample:
    #         dataSetSampleSeqDict[dsss.referenceSequenceOf.name] += dsss.abundance
    #
    #     # populate row in order of the master sorted intras
    #     for name in cladeAbundanceOrderedRefSeqList:
    #         if name in dataSetSampleSeqDict.keys():
    #             counts = dataSetSampleSeqDict[name]
    #             sampleRowDataCounts.append(dataSetSampleSeqDict[name])
    #             sampleRowDataProps.append(counts/sampleSeqTot)
    #         else:
    #             sampleRowDataCounts.append(0)
    #             sampleRowDataProps.append(0)
    #
    #     # Here we have all of the data populated
    #     # Now write to table
    #     intraAbundCountTable.append('{}\t{}'.format(sampleList[i].name, '\t'.join(str(x) for x in sampleRowDataCounts)))
    #     intraAbundPropTable.append('{}\t{}'.format(sampleList[i].name, '\t'.join("{:.3f}".format(x) for x in sampleRowDataProps)))
    #     # ["{:.3f}".format(prop) for prop in dataRowProp]

    # Now before writing out we need to add to the intraAbund objects from the managedSampleOutputDict in order of
    # sampleList
    for dss in sampleList:
        intraAbundCountTable.append(managedSampleOutputDict[dss][0])
        intraAbundPropTable.append(managedSampleOutputDict[dss][1])



    # Now add the accesion number / UID for each of the DIVs
    accession_list = []
    accession_list.append('DIV_accession')
    for ref_seq_name in cladeAbundanceOrderedRefSeqList:
        ref_seq = reference_sequence.objects.get(name=ref_seq_name)
        ref_seq_accession = ''
        if ref_seq.accession != 'None':
            ref_seq_accession = ref_seq.accession
        else:
            ref_seq_accession = str(ref_seq.id)
        accession_list.append(ref_seq_accession)

    intraAbundCountTable.append('\t'.join(accession_list))
    intraAbundPropTable.append('\t'.join(accession_list))

    # Here we have the tables populated and ready to output
    writeListToDestination('{}/{}_{}.DIVs.counts.txt'.format(outputDir, analysisObj.id,
                                                             '_'.join([str(a) for a in dataSubmissionsToOutput])),
                           intraAbundCountTable)
    writeListToDestination('{}/{}_{}.DIVs.props.txt'.format(outputDir, analysisObj.id,
                                                            '_'.join([str(a) for a in dataSubmissionsToOutput])),
                           intraAbundPropTable)

    return

def div_output_pre_analysis(datasubstooutput, numProcessors, output_dir, sorted_sample_list=None, analysis_obj_id=None):
    # This function will take in a list of dataSubmissions from which to output the DIV information of the samples
    # it will work very similarly to the second half of normal formatOutput_ord function. I will try to reuse as
    # much code as possible

    ########################## ITS2 INTRA ABUND COUNT TABLE ################################
    # This is where we're going to have to work with the sequences that aren't part of a type.
    # Esentially we want to end up with the noName sequences divieded up cladally.
    # So at the moment this will mean that we will divide up the current no names but also
    # add information about the other cladal sequences that didn't make it into a cladeCollection


    # Now we need to output the ITS2 intragenomic abundances that the types are based on for each CC
    # Already have querySetOfDataSubmissions and cladeCollectionsFromDSs

    # list to hold the paths of the outputted files
    output_path_list = []

    ################ GET ORDERED LIST OF INTRAS BY CLADE THEN ABUNDANCE #################
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    dataSubmissionsToOutput = [int(a) for a in datasubstooutput.split(',')]

    # outputDir = os.path.join(os.path.dirname(__file__), 'outputs/non_analysis/')
    # Get collection of types that are specific for the dataSubmissions we are looking at
    querySetOfDataSubmissions = data_set.objects.filter(id__in=dataSubmissionsToOutput)

    cladeCollectionsFromDSs = clade_collection.objects.filter(
        dataSetSampleFrom__dataSubmissionFrom__in=querySetOfDataSubmissions)

    cladeAbundanceOrderedRefSeqList = []
    refSeqsInDSs = reference_sequence.objects.filter(
        data_set_sample_sequence__cladeCollectionTwoFoundIn__in=cladeCollectionsFromDSs).distinct()

    # Get list of clades found
    cladeSet = set()
    for refSeq in refSeqsInDSs:
        cladeSet.add(refSeq.clade)
    # Order the list of clades alphabetically
    sub_clade_list = [a for a in clade_list if a in cladeSet]

    for i in range(len(sub_clade_list)):
        intraCounter = defaultdict(int)

        # refSeqs of clade
        refSeqsOfClade = refSeqsInDSs.filter(clade=sub_clade_list[i])

        # dataSetSampleSequencesOfClade
        dataSetSampleSequencesOfClade = data_set_sample_sequence.objects.filter(
            cladeCollectionTwoFoundIn__in=cladeCollectionsFromDSs, referenceSequenceOf__in=refSeqsOfClade)

        # Popualate the counter dict

        intraCounterManager = Manager()
        intraCounter = intraCounterManager.dict()

        dsssQueue = Queue()

        for dsss in dataSetSampleSequencesOfClade:
            dsssQueue.put(dsss)

        for N in range(numProcessors):
            dsssQueue.put('STOP')

        allProcesses = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        print('Creating refSeq counts for clade {}'.format(sub_clade_list[i]))
        for N in range(numProcessors):
            p = Process(target=outputWorkerTwo_pre_analysis, args=(dsssQueue, intraCounter))
            allProcesses.append(p)
            p.start()

        for p in allProcesses:
            p.join()

        # for dsss in dataSetSampleSequencesOfClade:
        #     intraCounter[dsss.referenceSequenceOf.name] += dsss.abundance

        # Sort counter dict by abundance
        sortedIntraCounter = sorted(intraCounter.items(), key=operator.itemgetter(1), reverse=True)

        # Add the refSeqNames to the master cladeAbundanceOrderedRefSeqList
        # Do not add the noName items. These will be summed across the clades and entered into the list
        # once their final abundances have been calculated
        cladeAbundanceOrderedRefSeqList.extend([a[0] for a in sortedIntraCounter if a[0] != 'noName'])



    # Create the two output tables as lists
    intraAbundCountTable = []
    intraAbundPropTable = []
    # Create the header line one for prop one for count
    # add the cladal noName categories directly to the header
    # do not add them to the cladeAbundanceOrderedRefSeqList so that we don't append 0s to them
    # we will add the abundances in the outputWorkerThree method below

    headerPre = 'Samples\t{}'.format('\t'.join([a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]))
    noNameStrings = '\t'.join(['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']])
    qc_stats = '\t'.join(['raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                          'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                          'post_taxa_id_absolute_non_symbiodinium_seqs',
                          'post_taxa_id_unique_non_symbiodinium_seqs'])
    # append the noName sequences as individual sequence abundances
    break_down_strings = '\t'.join([a for a in cladeAbundanceOrderedRefSeqList if '_' in a])
    header = '\t'.join([headerPre, noNameStrings, qc_stats, break_down_strings])

    intraAbundCountTable.append(header)
    intraAbundPropTable.append(header)
    ######################################################################################

    ############ POPULATE TABLE WITH CC DATA #############
    # For every sample
    sampleList = data_set_sample.objects.filter(dataSubmissionFrom__in=querySetOfDataSubmissions)

    # In order to MP this we will have to pay attention to the order. As we can't keep the order as we work with the
    # MPs we will do as we did above for the profies outputs and have an output dict that will have the sample as key
    # and its corresponding data row as the value. Then when we have finished the MPing we will go through the
    # sampleList order and output the data in this order.
    # So we will need a managed output dict.

    sampleOutputDictManager = Manager()
    managedSampleOutputDict = sampleOutputDictManager.dict()

    cladeAbundanceOrderedRefSeqListManager = Manager()
    cladeAbundanceOrderedRefSeqList = cladeAbundanceOrderedRefSeqListManager.list(cladeAbundanceOrderedRefSeqList)

    dssQueue = Queue()

    for dss in sampleList:
        dssQueue.put(dss)

    for N in range(numProcessors):
        dssQueue.put('STOP')

    allProcesses = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()

    for N in range(numProcessors):
        p = Process(target=outputWorkerThree_pre_analysis, args=(
        dssQueue, managedSampleOutputDict, cladeAbundanceOrderedRefSeqList, querySetOfDataSubmissions))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

    # Now before writing out we need to add to the intraAbund objects from the managedSampleOutputDict in order of
    # sampleList
    #### DEVELOPMENT ####
    # So that we can inspect the managedSampleOutputDict
    managedSampleOutputDict_dict = dict(managedSampleOutputDict)

    #####################

    # If there is a sorted sample list, make sure that it matches the samples that we are outputting

    if sorted_sample_list:
        sorted_sample_list_names = [a.name for a in sorted_sample_list]

        if len(sorted_sample_list_names) != len(sampleList):
            sys.exit('Name in sorted_sample_list do not match those to be outputted!')

        provided_name_to_samp_name_dict = {}
        for provided_name in sorted_sample_list_names:
            match_list = []
            for smp in sampleList:
                if provided_name in smp.name:
                    match_list.append(smp.name)
                    provided_name_to_samp_name_dict[provided_name] = smp
            if len(match_list) > 1:
                sys.exit('Sample name {} matches more than one output sample ({}).'.format(provided_name, '\t'.join(match_list)))
            if len(match_list) == 0:
                sys.exit('Sample name {} does not match any output sample.'.format(provided_name))

        # if we got to hear then the sorted_sample_list looks good

        for provided_name in sorted_sample_list_names:
            dss = provided_name_to_samp_name_dict[provided_name]
            intraAbundCountTable.append(managedSampleOutputDict[dss][0])
            intraAbundPropTable.append(managedSampleOutputDict[dss][1])
    else:
        for dss in sampleList:
            intraAbundCountTable.append(managedSampleOutputDict[dss][0])
            intraAbundPropTable.append(managedSampleOutputDict[dss][1])

    # Now add the accesion number / UID for each of the DIVs
    accession_list = []
    accession_list.append('DIV_accession')
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
        ref_seq = reference_sequence.objects.get(name=ref_seq_name)
        if ref_seq.accession and ref_seq.accession != 'None':
            ref_seq_accession = ref_seq.accession
        else:
            ref_seq_accession = str(ref_seq.id)
        accession_list.append(ref_seq_accession)

    accession_str = '\t'.join(accession_list)

    # now add the blanks for the noName cladal columns
    accession_str += ''.join(['\t' for cld in clade_list])
    # add blanks for the qc stats
    accession_str += ''.join(['\t' for i in range(7)])

    # now add the accession numbers for the noName break down seqs
    accession_list = []
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
        accession_list.append(ref_seq_name.split('_')[0])

    accession_str = '\t'.join([accession_str, '\t'.join(accession_list)])

    intraAbundCountTable.append(accession_str)
    intraAbundPropTable.append(accession_str)

    # Here we have the tables populated and ready to output
    if analysis_obj_id:
        path_to_div_absolute = '{}/{}_{}.DIVs.absolute.txt'.format(output_dir, analysis_obj_id,
                                                                '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        path_to_div_absolute = '{}/{}.DIVs.absolute.txt'.format(output_dir,
                                                                   '_'.join([str(a) for a in dataSubmissionsToOutput]))
    writeListToDestination(path_to_div_absolute, intraAbundCountTable)
    output_path_list.append(path_to_div_absolute)

    if analysis_obj_id:
        path_to_div_relative = '{}/{}_{}.DIVs.relative.txt'.format(output_dir, analysis_obj_id,
                                                                '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        path_to_div_relative = '{}/{}.DIVs.relative.txt'.format(output_dir, '_'.join([str(a) for a in dataSubmissionsToOutput]))

    writeListToDestination(path_to_div_relative, intraAbundPropTable)
    output_path_list.append(path_to_div_relative)

    # Output a .fasta for of all of the sequences found in the analysis
    fasta_output_list = []
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
        ref_seq = reference_sequence.objects.get(name=ref_seq_name)
        fasta_output_list.append('>{}'.format(ref_seq_name))
        fasta_output_list.append(ref_seq.sequence)
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
        ref_seq = reference_sequence.objects.get(id=int(ref_seq_name.split('_')[0]))
        fasta_output_list.append('>{}'.format(ref_seq_name))
        fasta_output_list.append(ref_seq.sequence)

    if analysis_obj_id:
        fasta_path = '{}/{}_{}.DIVs.fasta'.format(output_dir, analysis_obj_id, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        fasta_path = '{}/{}.DIVs.fasta'.format(output_dir, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    writeListToDestination(fasta_path, fasta_output_list)
    output_path_list.append(fasta_path)

    print('DIV table output files:')
    for path_item in output_path_list:
        print(path_item)

    return

def div_output_pre_analysis_new_meta(datasubstooutput, numProcessors, output_dir, sorted_sample_list=None, analysis_obj_id=None):
    # This function will take in a list of dataSubmissions from which to output the DIV information of the samples
    # it will work very similarly to the second half of normal formatOutput_ord function. I will try to reuse as
    # much code as possible

    ########################## ITS2 INTRA ABUND COUNT TABLE ################################
    # This is where we're going to have to work with the sequences that aren't part of a type.
    # Esentially we want to end up with the noName sequences divieded up cladally.
    # So at the moment this will mean that we will divide up the current no names but also
    # add information about the other cladal sequences that didn't make it into a cladeCollection


    # Now we need to output the ITS2 intragenomic abundances that the types are based on for each CC
    # Already have querySetOfDataSubmissions and cladeCollectionsFromDSs

    # list to hold the paths of the outputted files
    output_path_list = []

    ################ GET ORDERED LIST OF INTRAS BY CLADE THEN ABUNDANCE #################
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    dataSubmissionsToOutput = [int(a) for a in datasubstooutput.split(',')]

    # outputDir = os.path.join(os.path.dirname(__file__), 'outputs/non_analysis/')
    # Get collection of types that are specific for the dataSubmissions we are looking at
    querySetOfDataSubmissions = data_set.objects.filter(id__in=dataSubmissionsToOutput)

    cladeCollectionsFromDSs = clade_collection.objects.filter(
        dataSetSampleFrom__dataSubmissionFrom__in=querySetOfDataSubmissions)

    cladeAbundanceOrderedRefSeqList = []
    refSeqsInDSs = reference_sequence.objects.filter(
        data_set_sample_sequence__cladeCollectionTwoFoundIn__in=cladeCollectionsFromDSs).distinct()

    # Get list of clades found
    cladeSet = set()
    for refSeq in refSeqsInDSs:
        cladeSet.add(refSeq.clade)
    # Order the list of clades alphabetically
    sub_clade_list = [a for a in clade_list if a in cladeSet]

    for i in range(len(sub_clade_list)):
        intraCounter = defaultdict(int)

        # refSeqs of clade
        refSeqsOfClade = refSeqsInDSs.filter(clade=sub_clade_list[i])

        # dataSetSampleSequencesOfClade
        dataSetSampleSequencesOfClade = data_set_sample_sequence.objects.filter(
            cladeCollectionTwoFoundIn__in=cladeCollectionsFromDSs, referenceSequenceOf__in=refSeqsOfClade)

        # Popualate the counter dict

        intraCounterManager = Manager()
        intraCounter = intraCounterManager.dict()

        dsssQueue = Queue()

        for dsss in dataSetSampleSequencesOfClade:
            dsssQueue.put(dsss)

        for N in range(numProcessors):
            dsssQueue.put('STOP')

        allProcesses = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        print('Creating refSeq counts for clade {}'.format(sub_clade_list[i]))
        for N in range(numProcessors):
            p = Process(target=outputWorkerTwo_pre_analysis, args=(dsssQueue, intraCounter))
            allProcesses.append(p)
            p.start()

        for p in allProcesses:
            p.join()

        # for dsss in dataSetSampleSequencesOfClade:
        #     intraCounter[dsss.referenceSequenceOf.name] += dsss.abundance

        # Sort counter dict by abundance
        sortedIntraCounter = sorted(intraCounter.items(), key=operator.itemgetter(1), reverse=True)

        # Add the refSeqNames to the master cladeAbundanceOrderedRefSeqList
        # Do not add the noName items. These will be summed across the clades and entered into the list
        # once their final abundances have been calculated
        cladeAbundanceOrderedRefSeqList.extend([a[0] for a in sortedIntraCounter if a[0] != 'noName'])



    # Create the two output tables as lists
    intraAbundCountTable = []
    intraAbundPropTable = []
    # Create the header line one for prop one for count
    # add the cladal noName categories directly to the header
    # do not add them to the cladeAbundanceOrderedRefSeqList so that we don't append 0s to them
    # we will add the abundances in the outputWorkerThree method below

    headerPre = 'Samples\t{}'.format('\t'.join([a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]))
    noNameStrings = '\t'.join(['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']])
    qc_stats = '\t'.join(['raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                          'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                          'size_screening_violation_absolute', 'size_screening_violation_uniue',
                          'post_taxa_id_absolute_non_symbiodinium_seqs',
                          'post_taxa_id_unique_non_symbiodinium_seqs', 'post_med_absolute', 'post_med_unique'])
    # append the noName sequences as individual sequence abundances
    break_down_strings = '\t'.join([a for a in cladeAbundanceOrderedRefSeqList if '_' in a])
    header = '\t'.join([headerPre, noNameStrings, qc_stats, break_down_strings])

    intraAbundCountTable.append(header)
    intraAbundPropTable.append(header)
    ######################################################################################

    ############ POPULATE TABLE WITH CC DATA #############
    # For every sample
    sampleList = data_set_sample.objects.filter(dataSubmissionFrom__in=querySetOfDataSubmissions)

    # In order to MP this we will have to pay attention to the order. As we can't keep the order as we work with the
    # MPs we will do as we did above for the profies outputs and have an output dict that will have the sample as key
    # and its corresponding data row as the value. Then when we have finished the MPing we will go through the
    # sampleList order and output the data in this order.
    # So we will need a managed output dict.

    sampleOutputDictManager = Manager()
    managedSampleOutputDict = sampleOutputDictManager.dict()

    cladeAbundanceOrderedRefSeqListManager = Manager()
    cladeAbundanceOrderedRefSeqList = cladeAbundanceOrderedRefSeqListManager.list(cladeAbundanceOrderedRefSeqList)

    dssQueue = Queue()

    for dss in sampleList:
        dssQueue.put(dss)

    for N in range(numProcessors):
        dssQueue.put('STOP')

    allProcesses = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()

    for N in range(numProcessors):
        p = Process(target=outputWorkerThree_pre_analysis, args=(
        dssQueue, managedSampleOutputDict, cladeAbundanceOrderedRefSeqList, querySetOfDataSubmissions))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

    # Now before writing out we need to add to the intraAbund objects from the managedSampleOutputDict in order of
    # sampleList
    #### DEVELOPMENT ####
    # So that we can inspect the managedSampleOutputDict
    managedSampleOutputDict_dict = dict(managedSampleOutputDict)

    #####################

    # If there is a sorted sample list, make sure that it matches the samples that we are outputting

    if sorted_sample_list:
        sorted_sample_list_names = [a.name for a in sorted_sample_list]

        if len(sorted_sample_list_names) != len(sampleList):
            sys.exit('Name in sorted_sample_list do not match those to be outputted!')

        provided_name_to_samp_name_dict = {}
        for provided_name in sorted_sample_list_names:
            match_list = []
            for smp in sampleList:
                if provided_name in smp.name:
                    match_list.append(smp.name)
                    provided_name_to_samp_name_dict[provided_name] = smp
            if len(match_list) > 1:
                sys.exit('Sample name {} matches more than one output sample ({}).'.format(provided_name, '\t'.join(match_list)))
            if len(match_list) == 0:
                sys.exit('Sample name {} does not match any output sample.'.format(provided_name))

        # if we got to hear then the sorted_sample_list looks good

        for provided_name in sorted_sample_list_names:
            dss = provided_name_to_samp_name_dict[provided_name]
            intraAbundCountTable.append(managedSampleOutputDict[dss][0])
            intraAbundPropTable.append(managedSampleOutputDict[dss][1])
    else:
        for dss in sampleList:
            intraAbundCountTable.append(managedSampleOutputDict[dss][0])
            intraAbundPropTable.append(managedSampleOutputDict[dss][1])

    # Now add the accesion number / UID for each of the DIVs
    accession_list = []
    accession_list.append('DIV_accession')
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
        ref_seq = reference_sequence.objects.get(name=ref_seq_name)
        if ref_seq.accession and ref_seq.accession != 'None':
            ref_seq_accession = ref_seq.accession
        else:
            ref_seq_accession = str(ref_seq.id)
        accession_list.append(ref_seq_accession)

    accession_str = '\t'.join(accession_list)

    # now add the blanks for the noName cladal columns
    accession_str += ''.join(['\t' for cld in clade_list])
    # add blanks for the qc stats
    accession_str += ''.join(['\t' for i in range(11)])

    # now add the accession numbers for the noName break down seqs
    accession_list = []
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
        accession_list.append(ref_seq_name.split('_')[0])

    accession_str = '\t'.join([accession_str, '\t'.join(accession_list)])

    intraAbundCountTable.append(accession_str)
    intraAbundPropTable.append(accession_str)

    # Here we have the tables populated and ready to output
    if analysis_obj_id:
        path_to_div_absolute = '{}/{}_{}.DIVs.absolute.txt'.format(output_dir, analysis_obj_id,
                                                                '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        path_to_div_absolute = '{}/{}.DIVs.absolute.txt'.format(output_dir,
                                                                   '_'.join([str(a) for a in dataSubmissionsToOutput]))
    writeListToDestination(path_to_div_absolute, intraAbundCountTable)
    output_path_list.append(path_to_div_absolute)

    if analysis_obj_id:
        path_to_div_relative = '{}/{}_{}.DIVs.relative.txt'.format(output_dir, analysis_obj_id,
                                                                '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        path_to_div_relative = '{}/{}.DIVs.relative.txt'.format(output_dir, '_'.join([str(a) for a in dataSubmissionsToOutput]))

    writeListToDestination(path_to_div_relative, intraAbundPropTable)
    output_path_list.append(path_to_div_relative)

    # Output a .fasta for of all of the sequences found in the analysis
    fasta_output_list = []
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
        ref_seq = reference_sequence.objects.get(name=ref_seq_name)
        fasta_output_list.append('>{}'.format(ref_seq_name))
        fasta_output_list.append(ref_seq.sequence)
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
        ref_seq = reference_sequence.objects.get(id=int(ref_seq_name.split('_')[0]))
        fasta_output_list.append('>{}'.format(ref_seq_name))
        fasta_output_list.append(ref_seq.sequence)

    if analysis_obj_id:
        fasta_path = '{}/{}_{}.DIVs.fasta'.format(output_dir, analysis_obj_id, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        fasta_path = '{}/{}.DIVs.fasta'.format(output_dir, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    writeListToDestination(fasta_path, fasta_output_list)
    output_path_list.append(fasta_path)

    print('DIV table output files:')
    for path_item in output_path_list:
        print(path_item)

    return

def div_output_pre_analysis_new_meta_and_new_dss_structure(datasubstooutput, numProcessors, output_dir, sorted_sample_list=None, analysis_obj_id=None):
    # This function will take in a list of dataSubmissions from which to output the DIV information of the samples
    # it will work very similarly to the second half of normal formatOutput_ord function. I will try to reuse as
    # much code as possible

    ########################## ITS2 INTRA ABUND COUNT TABLE ################################
    # This is where we're going to have to work with the sequences that aren't part of a type.
    # Esentially we want to end up with the noName sequences divieded up cladally.
    # So at the moment this will mean that we will divide up the current no names but also
    # add information about the other cladal sequences that didn't make it into a cladeCollection


    # Now we need to output the ITS2 intragenomic abundances that the types are based on for each CC
    # Already have querySetOfDataSubmissions and cladeCollectionsFromDSs

    # list to hold the paths of the outputted files
    output_path_list = []

    ################ GET ORDERED LIST OF INTRAS BY CLADE THEN ABUNDANCE #################
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    dataSubmissionsToOutput = [int(a) for a in datasubstooutput.split(',')]

    # outputDir = os.path.join(os.path.dirname(__file__), 'outputs/non_analysis/')
    # Get collection of types that are specific for the dataSubmissions we are looking at
    querySetOfDataSubmissions = data_set.objects.filter(id__in=dataSubmissionsToOutput)

    # cladeCollectionsFromDSs = clade_collection.objects.filter(
    #     dataSetSampleFrom__dataSubmissionFrom__in=querySetOfDataSubmissions)

    cladeAbundanceOrderedRefSeqList = []
    refSeqsInDSs = reference_sequence.objects.filter(
        data_set_sample_sequence__data_set_sample_from__dataSubmissionFrom__in=querySetOfDataSubmissions).distinct()

    # Get list of clades found
    cladeSet = set()
    for refSeq in refSeqsInDSs:
        cladeSet.add(refSeq.clade)
    # Order the list of clades alphabetically
    sub_clade_list = [a for a in clade_list if a in cladeSet]

    for i in range(len(sub_clade_list)):
        intraCounter = defaultdict(int)

        # refSeqs of clade
        refSeqsOfClade = refSeqsInDSs.filter(clade=sub_clade_list[i])

        # dataSetSampleSequencesOfClade
        dataSetSampleSequencesOfClade = data_set_sample_sequence.objects.filter(
            data_set_sample_from__dataSubmissionFrom__in=querySetOfDataSubmissions, referenceSequenceOf__in=refSeqsOfClade)

        # Popualate the counter dict

        intraCounterManager = Manager()
        intraCounter = intraCounterManager.dict()

        dsssQueue = Queue()

        for dsss in dataSetSampleSequencesOfClade:
            dsssQueue.put(dsss)

        for N in range(numProcessors):
            dsssQueue.put('STOP')

        allProcesses = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        print('Creating refSeq counts for clade {}'.format(sub_clade_list[i]))
        for N in range(numProcessors):
            p = Process(target=outputWorkerTwo_pre_analysis, args=(dsssQueue, intraCounter))
            allProcesses.append(p)
            p.start()

        for p in allProcesses:
            p.join()

        # for dsss in dataSetSampleSequencesOfClade:
        #     intraCounter[dsss.referenceSequenceOf.name] += dsss.abundance

        # Sort counter dict by abundance
        sortedIntraCounter = sorted(intraCounter.items(), key=operator.itemgetter(1), reverse=True)

        # Add the refSeqNames to the master cladeAbundanceOrderedRefSeqList
        # Do not add the noName items. These will be summed across the clades and entered into the list
        # once their final abundances have been calculated
        cladeAbundanceOrderedRefSeqList.extend([a[0] for a in sortedIntraCounter if a[0] != 'noName'])



    # Create the two output tables as lists
    intraAbundCountTable = []
    intraAbundPropTable = []
    # Create the header line one for prop one for count
    # add the cladal noName categories directly to the header
    # do not add them to the cladeAbundanceOrderedRefSeqList so that we don't append 0s to them
    # we will add the abundances in the outputWorkerThree method below

    headerPre = 'Samples\t{}'.format('\t'.join([a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]))
    noNameStrings = '\t'.join(['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']])
    qc_stats = '\t'.join(['raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                          'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                          'size_screening_violation_absolute', 'size_screening_violation_uniue',
                          'post_taxa_id_absolute_non_symbiodinium_seqs',
                          'post_taxa_id_unique_non_symbiodinium_seqs', 'post_med_absolute', 'post_med_unique'])
    # append the noName sequences as individual sequence abundances
    break_down_strings = '\t'.join([a for a in cladeAbundanceOrderedRefSeqList if '_' in a])
    header = '\t'.join([headerPre, noNameStrings, qc_stats, break_down_strings])

    intraAbundCountTable.append(header)
    intraAbundPropTable.append(header)
    ######################################################################################

    ############ POPULATE TABLE WITH CC DATA #############
    # For every sample
    sampleList = data_set_sample.objects.filter(dataSubmissionFrom__in=querySetOfDataSubmissions)

    # In order to MP this we will have to pay attention to the order. As we can't keep the order as we work with the
    # MPs we will do as we did above for the profies outputs and have an output dict that will have the sample as key
    # and its corresponding data row as the value. Then when we have finished the MPing we will go through the
    # sampleList order and output the data in this order.
    # So we will need a managed output dict.

    sampleOutputDictManager = Manager()
    managedSampleOutputDict = sampleOutputDictManager.dict()

    cladeAbundanceOrderedRefSeqListManager = Manager()
    cladeAbundanceOrderedRefSeqList = cladeAbundanceOrderedRefSeqListManager.list(cladeAbundanceOrderedRefSeqList)

    dssQueue = Queue()

    for dss in sampleList:
        dssQueue.put(dss)

    for N in range(numProcessors):
        dssQueue.put('STOP')

    allProcesses = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()

    for N in range(numProcessors):
        p = Process(target=outputWorkerThree_pre_analysis_new_dss_structure, args=(
        dssQueue, managedSampleOutputDict, cladeAbundanceOrderedRefSeqList, querySetOfDataSubmissions))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

    # Now before writing out we need to add to the intraAbund objects from the managedSampleOutputDict in order of
    # sampleList
    #### DEVELOPMENT ####
    # So that we can inspect the managedSampleOutputDict
    managedSampleOutputDict_dict = dict(managedSampleOutputDict)

    #####################

    # If there is a sorted sample list, make sure that it matches the samples that we are outputting

    if sorted_sample_list:
        sorted_sample_list_names = [a.name for a in sorted_sample_list]

        if len(sorted_sample_list_names) != len(sampleList):
            sys.exit('Name in sorted_sample_list do not match those to be outputted!')

        provided_name_to_samp_name_dict = {}
        for provided_name in sorted_sample_list_names:
            match_list = []
            for smp in sampleList:
                if provided_name in smp.name:
                    match_list.append(smp.name)
                    provided_name_to_samp_name_dict[provided_name] = smp
            if len(match_list) > 1:
                sys.exit('Sample name {} matches more than one output sample ({}).'.format(provided_name, '\t'.join(match_list)))
            if len(match_list) == 0:
                sys.exit('Sample name {} does not match any output sample.'.format(provided_name))

        # if we got to hear then the sorted_sample_list looks good

        for provided_name in sorted_sample_list_names:
            dss = provided_name_to_samp_name_dict[provided_name]
            intraAbundCountTable.append(managedSampleOutputDict[dss][0])
            intraAbundPropTable.append(managedSampleOutputDict[dss][1])
    else:
        for dss in sampleList:
            intraAbundCountTable.append(managedSampleOutputDict[dss][0])
            intraAbundPropTable.append(managedSampleOutputDict[dss][1])

    # Now add the accesion number / UID for each of the DIVs
    accession_list = []
    accession_list.append('DIV_accession')
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
        ref_seq = reference_sequence.objects.get(name=ref_seq_name)
        if ref_seq.accession and ref_seq.accession != 'None':
            ref_seq_accession = ref_seq.accession
        else:
            ref_seq_accession = str(ref_seq.id)
        accession_list.append(ref_seq_accession)

    accession_str = '\t'.join(accession_list)

    # now add the blanks for the noName cladal columns
    accession_str += ''.join(['\t' for cld in clade_list])
    # add blanks for the qc stats
    accession_str += ''.join(['\t' for i in range(11)])

    # now add the accession numbers for the noName break down seqs
    accession_list = []
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
        accession_list.append(ref_seq_name.split('_')[0])

    accession_str = '\t'.join([accession_str, '\t'.join(accession_list)])

    intraAbundCountTable.append(accession_str)
    intraAbundPropTable.append(accession_str)

    # Here we have the tables populated and ready to output
    if analysis_obj_id:
        path_to_div_absolute = '{}/{}_{}.DIVs.absolute.txt'.format(output_dir, analysis_obj_id,
                                                                '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        path_to_div_absolute = '{}/{}.DIVs.absolute.txt'.format(output_dir,
                                                                   '_'.join([str(a) for a in dataSubmissionsToOutput]))
    writeListToDestination(path_to_div_absolute, intraAbundCountTable)
    output_path_list.append(path_to_div_absolute)

    if analysis_obj_id:
        path_to_div_relative = '{}/{}_{}.DIVs.relative.txt'.format(output_dir, analysis_obj_id,
                                                                '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        path_to_div_relative = '{}/{}.DIVs.relative.txt'.format(output_dir, '_'.join([str(a) for a in dataSubmissionsToOutput]))

    writeListToDestination(path_to_div_relative, intraAbundPropTable)
    output_path_list.append(path_to_div_relative)

    # Output a .fasta for of all of the sequences found in the analysis
    fasta_output_list = []
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
        ref_seq = reference_sequence.objects.get(name=ref_seq_name)
        fasta_output_list.append('>{}'.format(ref_seq_name))
        fasta_output_list.append(ref_seq.sequence)
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
        ref_seq = reference_sequence.objects.get(id=int(ref_seq_name.split('_')[0]))
        fasta_output_list.append('>{}'.format(ref_seq_name))
        fasta_output_list.append(ref_seq.sequence)

    if analysis_obj_id:
        fasta_path = '{}/{}_{}.DIVs.fasta'.format(output_dir, analysis_obj_id, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        fasta_path = '{}/{}.DIVs.fasta'.format(output_dir, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    writeListToDestination(fasta_path, fasta_output_list)
    output_path_list.append(fasta_path)

    print('DIV table output files:')
    for path_item in output_path_list:
        print(path_item)

    return

def formatOutput_ord(analysisobj, numProcessors, datasubstooutput):

    analysisObj = analysisobj
    # This is one of the last things to do before we can use our first dataset
    # The table will have types as columns and rows as samples
    # its rows will give various data about each type as well as how much they were present in each sample

    # It will produce four output files, for DIV abundances and proportions and Type abundances and proportions.
    # found in the given clade collection.
    # Types will be listed first by clade and then by the number of clade collections the types were found in
    # The species groups will be numbered consecutively as they appear
    # Each type's ID will also be given as a UID.
    # The date the database was accessed and the version should also be noted
    # The formal species descriptions which correspond to the found ITS2 type will be noted for each type
    # Finally the AccessionNumber of each of the defining reference species will also be noted
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    '''
    32: VoolstraSet
    40: EmilyData
    47: EdData
    51: ed_24_07_17_singapore
    52: UTS_ONE_25_07_17
    53: UTS_TWO_25_07_17
    54: 20171008_jitErn
    55: 20171009_roberto
    56: 20171016_alejandro
    57: 20171016_Pappas
    58: 20171017_roberto
    59: 20171017_tullia
    60: 20170719_alejandro
    61: 20170713_alejandro
    62: 20170506_roberto
    63: ryan_first
    '''


    # List of the paths to the files that have been output
    output_files_list = []

    dataSubmissionsToOutput = [int(a) for a in datasubstooutput.split(',')]

    # Get collection of types that are specific for the dataSubmissions we are looking at
    querySetOfDataSubmissions = data_set.objects.filter(id__in=dataSubmissionsToOutput)

    cladeCollectionsFromDSs = clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom__in=querySetOfDataSubmissions)
    cladeCollectionTypesFromThisAnalysisObjAndDataSubmission = clade_collection_type.objects.filter(cladeCollectionFoundIn__in=cladeCollectionsFromDSs, analysisTypeOf__dataAnalysisFrom=analysisObj).distinct()

    # cladeCollectionsFromDSs = list(clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom__in=dataSubmissionsToOutput))
    # cladeCollectionTypesFromThisAnalysisObjAndDataSubmission = list(clade_collection_type.objects.filter(cladeCollectionFoundIn__in=cladeCollectionsFromDSs, analysisTypeOf__dataAnalysisFrom=analysisObj).distinct())

    at = set()
    for cct in cladeCollectionTypesFromThisAnalysisObjAndDataSubmission:
        at.add(cct.analysisTypeOf)
    at = list(at)

    # at = analysis_type.objects.filter(dataAnalysisFrom=analysisObj)
    IDs = [att.id for att in at]
    outputTableOne = []
    outputTableTwo = []
    # Need to get a list of Samples that are part of the dataAnalysis
    #listOfDataSetSamples = list(data_set_sample.objects.filter(dataSubmissionFrom__in=data_set.objects.filter(id__in=[int(a) for a in analysisObj.listOfDataSubmissions.split(',')])))
    listOfDataSetSamples = list(data_set_sample.objects.filter(dataSubmissionFrom__in=data_set.objects.filter(id__in=dataSubmissionsToOutput)))
    # The header row for the table
    # outputTableOne.append('ITS2 type profile UID\tClade\tSpeciesGrp\tType_Profile\tMajority ITS2 Seq\tAssociated Species\tType abund local\tType abund DB\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))
    # outputTableTwo.append('ITS2 type profile UID\tClade\tSpeciesGrp\tType_Profile\tMajority ITS2 Seq\tAssociated Species\tType abund local\tType abund DB\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))

    # I am making some changes to the output format for the SP manuscript
    # outputTableOne.append('ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))
    # outputTableTwo.append('ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))



    # Now go through the types, firstly by clade and then by the number of cladeCollections they were found in
    # Populate a 2D list of types with a list per clade
    typesCladalList = [[] for clade in cladeList]
    for att in at:
        try:
            if len(att.listOfCladeCollections) > 0:
                typesCladalList[cladeList.index(att.clade)].append(att)
        except:
            pass


    # Items for for creating the new samples sorted output
    across_clade_type_sample_abund_dict = dict()
    across_clade_sorted_type_order = []
    # Go clade by clade
    for i in range(len(typesCladalList)):
        if typesCladalList[i]:

            ###### CALCULATE REL ABUND AND SD OF DEF INTRAS FOR THIS TYPE ###############
            #     # For each type we want to calculate the average proportions of the defining seqs in that type
            #     # We will also calculate an SD.
            #     # We will do both of these calculations using the footprintSeqAbundances, listOfCladeCollections
            #     # and orderedFoorprintList attributes of each type

            ########## MAKE GROUP COUNTER AND DICT ###########
            # These are now going to be managed items for use in the MP
            # We want to name the groups that types are found in sequencially
            # To get the next number to name a group we will use the groupCount
            # To look up what number has been assigned to groups that have already been printed we will use the groupDict
            # groupCount = 0
            # groupDict = {}
            ###################################################

            # This is currently being sorted according to the number of listOfCladeCollections the type was found in
            # which translates to the number of clade collections across the whole data analysis, i.e. samples from this
            #  out put also from any other samples inlcuded in the dataanalysis we want to try to output according to
            # the number of samples of the output it was found in instead
            # get a sorted list of the types by most commonly found first
            sortedListOfTypes = sorted(typesCladalList[i], key=lambda x: len(x.listOfCladeCollections.split(',')), reverse=True)

            # We need to produce a list of analysisTypes according to the number of cladecollections of the output
            # they are found in. This shouldn't be too tricky.
            # We will pass in typesCladalList[i] and we will simply look to see which of the .listOfCladeCollections
            # of each type are found in the cladeCollectionTypesFromThisAnalysisObjAndDataSubmission

            sortedListOfTypes = sort_list_of_types_by_clade_collections_in_current_output(typesCladalList[i], cladeCollectionTypesFromThisAnalysisObjAndDataSubmission)

            # Here we will MP at the type level
            # i.e. each type will be processed on a differnt core. In order for this to work we should have a managed
            # dictionary where the key can be the types ID and the value can be the datalist
            # In order to get the types output in the correct order we should use the sortedListOfTypes to resort the
            # data once the MPing has been done.

            # This dict will hold all of the output rows that the MPing has created.
            typeOutputManager = Manager()
            typeOutputManagedDict = typeOutputManager.dict({type: None for type in sortedListOfTypes})
            # This dict will hold details of which groups have already been written
            groupDictManager = Manager()
            groupDictManagedDict = groupDictManager.dict()

            groupCountManager = Manager()
            groupCountManagedInt = groupCountManager.Value(c_int, 0)

            listOfDataSetSamplesManager = Manager()
            listOfDataSetSamplesManagedList = listOfDataSetSamplesManager.list(listOfDataSetSamples)

            typeInputQueue = Queue()

            for anType in sortedListOfTypes:
                typeInputQueue.put(anType)


            for N in range(numProcessors):
                typeInputQueue.put('STOP')

            allProcesses = []

            # close all connections to the db so that they are automatically recreated for each process
            # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
            db.connections.close_all()

            for N in range(numProcessors):
                p = Process(target=outputWorkerOne, args=(typeInputQueue, groupDictManagedDict, groupCountManagedInt, listOfDataSetSamplesManagedList, cladeList[i], typeOutputManagedDict))
                allProcesses.append(p)
                p.start()

            for p in allProcesses:
                p.join()




            # Here the MP is completed and we now need to add the lists that are held in the managed dict to
            # the outputTableOne and outputTableTwo

            #### Development ####
            # So we have the current list of samples stored in a mangaed list that is the listOfDataSetSamplesManaged
            # the dictionary that we have below is key = analysis_type and the value is two lists
            # the first list is the raw counts and the second is the prop counts
            # In each of these lists the string is a string of the row values for each of the types
            # we can extract the seq abundances from these values.
            # we are currently doing this by clade so we will have to extract this info for each clade
            # I think we can simply add the items for each of the within clade dicts to each other and end up with
            # an across clades dict then we should be able to work with this to get the order we are looking for

            across_clade_type_sample_abund_dict.update(typeOutputManagedDict)
            across_clade_sorted_type_order.extend(sortedListOfTypes)
            typeOutputManagedDict_dict = dict(typeOutputManagedDict)
            #####################


            # for antype in sortedListOfTypes:
            #     outputTableOne.append(typeOutputManagedDict[antype][0])
            #     outputTableTwo.append(typeOutputManagedDict[antype][1])


   # At this stage we have the across_clade_type_sample_abund_dict that contains the type to row values
    # Below is the pseudo code for getting the sorted list of samples
    # get a list of the types for the output
    # Ideally we want to have this list of types sorted according to the most abundant in the output first
    # we currently don't have this except separated by clade
    # We can get this from the dict though
    type_to_abund_list = []
    for analysis_type_obj in across_clade_type_sample_abund_dict.keys():
        type_to_abund_list.append((analysis_type_obj, int(across_clade_type_sample_abund_dict[analysis_type_obj][0].split('\t')[4])))

    # now we can sort this list according to the abundance and this will give us the order of types that we want
    sorted_analysis_type_abundance_list = [a[0] for a in sorted(type_to_abund_list, key=lambda x: x[1], reverse=True)]

    # now that we have the sorted list of types
    # we want to extract some information
    # For each sample, we want to have the abundance of each of the types in that sample in the order of the types
    dataSetSample_to_type_abundance_dict = dict()
    for i in range(len(listOfDataSetSamples)):
        temp_abundance_list = []
        for analysis_type_obj in sorted_analysis_type_abundance_list:
            temp_abundance = int(across_clade_type_sample_abund_dict[analysis_type_obj][0].split('\t')[7:-1][i])
            temp_abundance_list.append(temp_abundance)
        dataSetSample_to_type_abundance_dict[listOfDataSetSamples[i]] = temp_abundance_list

    # list of samples that still need to be sorted
    samples_still_to_be_sorted = listOfDataSetSamples

    # list of samples that have been added to the sorted list
    samples_that_have_been_sorted = []

    # for each analysis type in order of the most abundant types first
    for j in range(len(sorted_analysis_type_abundance_list)):
        print('Sorting {}'.format(sorted_analysis_type_abundance_list[j]))
        # go through all of the samples and if a sample has the type as its most abundant type
        # then add it and its relative abundance of the type to to the below list
        list_of_samples_containing_analysis_type = []
        for data_set_sample_obj in dataSetSample_to_type_abundance_dict.keys():
            if data_set_sample_obj not in samples_that_have_been_sorted:
                # see if the type in question is the most abundant type in the sample in question
                list_of_abundances = dataSetSample_to_type_abundance_dict[data_set_sample_obj]
                abundance_of_type_in_sample = list_of_abundances[j]
                if abundance_of_type_in_sample != 0:
                    # then we need to see if this type is the most abundant in the sample
                    if all(i <= abundance_of_type_in_sample for i in list_of_abundances):
                        # then this type is the most abundant type in the sample and it should be added to the
                        # list_of_samples_containing_analysis_type
                        list_of_samples_containing_analysis_type.append((data_set_sample_obj, abundance_of_type_in_sample/sum(list_of_abundances)))
        # now that we've been through all of the samples for the given type
        # we need to sort the list_of_samples_containing_analysis_type by the second value of the tuple
        # we then need to add the samples in this order to the sample_that_have_been_sorted and remove them
        # from the samples_still_to_be_sorted
        sorted_list_of_samples_containing_analysis_type = sorted(list_of_samples_containing_analysis_type, key=lambda x: x[1], reverse=True)
        sorted_list_of_samples_for_this_type = [a[0] for a in sorted_list_of_samples_containing_analysis_type]

        for data_set_sample_obj in sorted_list_of_samples_for_this_type:
            samples_that_have_been_sorted.append(data_set_sample_obj)

    # here the samples_that_have_been_sorted should be populated in the sorted order that we want

    # Here the samples we should be left with should be samples that don't have a type associated to them
    # negatives etc.
    samples_not_ordered = list(set(listOfDataSetSamples) - set(samples_that_have_been_sorted))

    # add these samples to the list
    samples_that_have_been_sorted.extend(samples_not_ordered)

    # these are the tables that the ITS2 type data will be stored in.
    outputTableOne.append(
        'ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format(
            '\t'.join([dataSamp.name for dataSamp in samples_that_have_been_sorted])))
    outputTableTwo.append(
        'ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format(
            '\t'.join([dataSamp.name for dataSamp in samples_that_have_been_sorted])))
    
    # at this point we need to rearrange the data that is held in the current row strings
    for analysis_type_obj in across_clade_sorted_type_order:
        current_abundance_list_counts = across_clade_type_sample_abund_dict[analysis_type_obj][0].split('\t')[7:-1]
        current_abundance_list_props = across_clade_type_sample_abund_dict[analysis_type_obj][1].split('\t')[7:-1]
        new_abundance_list_counts = []
        new_abundance_list_props = []
        for new_sample in samples_that_have_been_sorted:
            new_abundance_list_counts.append(current_abundance_list_counts[listOfDataSetSamples.index(new_sample)])
            new_abundance_list_props.append(current_abundance_list_props[listOfDataSetSamples.index(new_sample)])

        new_string_counts = '{}\t{}\t{}'.format('\t'.join(across_clade_type_sample_abund_dict[analysis_type_obj][0].split('\t')[:7]), '\t'.join(new_abundance_list_counts), across_clade_type_sample_abund_dict[analysis_type_obj][0].split('\t')[-1])
        new_sting_props = '{}\t{}\t{}'.format('\t'.join(across_clade_type_sample_abund_dict[analysis_type_obj][1].split('\t')[:7]), '\t'.join(new_abundance_list_props), across_clade_type_sample_abund_dict[analysis_type_obj][1].split('\t')[-1])

        outputTableOne.append(new_string_counts)
        outputTableTwo.append(new_sting_props)





    # At this point we should have gone through each of the clades in order and for each clade we should have
    # created an analysis_type row in order of the sortedListOfTypes.

    # The current out put has samples on the Y and types on the X.
    # We want to reverse this as there will often be more samples than types
    # https://stackoverflow.com/questions/6473679/transpose-list-of-lists

    # Put into list of lists format to transpose
    newOutOne = [a.split('\t') for a in outputTableOne]
    newOutTwo = [a.split('\t') for a in outputTableTwo]

    # Transpose
    outputTableOne = list(map(list, zip(*newOutOne)))
    outputTableTwo = list(map(list, zip(*newOutTwo)))

    # Put back into tab delim format
    outputTableOne = ['\t'.join(a) for a in outputTableOne]
    outputTableTwo = ['\t'.join(a) for a in outputTableTwo]

    outputDir = os.path.join(os.path.dirname(__file__), 'outputs/{}'.format(analysisObj.id))
    os.makedirs(outputDir, exist_ok=True)
    os.chdir(outputDir)

    # Finally append the species references to the tables
    species_ref_dict = {
        'S. microadriaticum' : 'Freudenthal, H. D. (1962). Symbiodinium gen. nov. and Symbiodinium microadriaticum '
                               'sp. nov., a Zooxanthella: Taxonomy, Life Cycle, and Morphology. The Journal of '
                               'Protozoology 9(1): 45-52',
        'S. pilosum' : 'Trench, R. (2000). Validation of some currently used invalid names of dinoflagellates. '
                       'Journal of Phycology 36(5): 972-972.\tTrench, R. K. and R. J. Blank (1987). '
                       'Symbiodinium microadriaticum Freudenthal, S. goreauii sp. nov., S. kawagutii sp. nov. and '
                       'S. pilosum sp. nov.: Gymnodinioid dinoflagellate symbionts of marine invertebrates. '
                       'Journal of Phycology 23(3): 469-481.',
        'S. natans' : 'Hansen, G. and N. Daugbjerg (2009). Symbiodinium natans sp. nob.: A free-living '
                      'dinoflagellate from Tenerife (northeast Atlantic Ocean). Journal of Phycology 45(1): 251-263.',
        'S. tridacnidorum' : 'Lee, S. Y., H. J. Jeong, N. S. Kang, T. Y. Jang, S. H. Jang and T. C. Lajeunesse (2015). '
                             'Symbiodinium tridacnidorum sp. nov., a dinoflagellate common to Indo-Pacific giant clams,'
                             ' and a revised morphological description of Symbiodinium microadriaticum Freudenthal, '
                             'emended Trench & Blank. European Journal of Phycology 50(2): 155-172.',
        'S. linucheae' : 'Trench, R. K. and L.-v. Thinh (1995). Gymnodinium linucheae sp. nov.: The dinoflagellate '
                         'symbiont of the jellyfish Linuche unguiculata. European Journal of Phycology 30(2): 149-154.',
        'S. minutum' : 'Lajeunesse, T. C., J. E. Parkinson and J. D. Reimer (2012). A genetics-based description of '
                       'Symbiodinium minutum sp. nov. and S. psygmophilum sp. nov. (dinophyceae), two dinoflagellates '
                       'symbiotic with cnidaria. Journal of Phycology 48(6): 1380-1391.',
        'S. antillogorgium' : 'Parkinson, J. E., M. A. Coffroth and T. C. LaJeunesse (2015). "New species of Clade B '
                              'Symbiodinium (Dinophyceae) from the greater Caribbean belong to different functional '
                              'guilds: S. aenigmaticum sp. nov., S. antillogorgium sp. nov., S. endomadracis sp. nov., '
                              'and S. pseudominutum sp. nov." Journal of phycology 51(5): 850-858.',
        'S. pseudominutum' : 'Parkinson, J. E., M. A. Coffroth and T. C. LaJeunesse (2015). "New species of Clade B '
                              'Symbiodinium (Dinophyceae) from the greater Caribbean belong to different functional '
                              'guilds: S. aenigmaticum sp. nov., S. antillogorgium sp. nov., S. endomadracis sp. nov., '
                              'and S. pseudominutum sp. nov." Journal of phycology 51(5): 850-858.',
        'S. psygmophilum' : 'Lajeunesse, T. C., J. E. Parkinson and J. D. Reimer (2012). A genetics-based description of '
                       'Symbiodinium minutum sp. nov. and S. psygmophilum sp. nov. (dinophyceae), two dinoflagellates '
                       'symbiotic with cnidaria. Journal of Phycology 48(6): 1380-1391.',
        'S. muscatinei' : 'No reference available',
        'S. endomadracis' : 'Parkinson, J. E., M. A. Coffroth and T. C. LaJeunesse (2015). "New species of Clade B '
                              'Symbiodinium (Dinophyceae) from the greater Caribbean belong to different functional '
                              'guilds: S. aenigmaticum sp. nov., S. antillogorgium sp. nov., S. endomadracis sp. nov., '
                              'and S. pseudominutum sp. nov." Journal of phycology 51(5): 850-858.',
        'S. aenigmaticum' : 'Parkinson, J. E., M. A. Coffroth and T. C. LaJeunesse (2015). "New species of Clade B '
                              'Symbiodinium (Dinophyceae) from the greater Caribbean belong to different functional '
                              'guilds: S. aenigmaticum sp. nov., S. antillogorgium sp. nov., S. endomadracis sp. nov., '
                              'and S. pseudominutum sp. nov." Journal of phycology 51(5): 850-858.',
        'S. goreaui' : 'Trench, R. (2000). Validation of some currently used invalid names of dinoflagellates. '
                       'Journal of Phycology 36(5): 972-972.\tTrench, R. K. and R. J. Blank (1987). '
                       'Symbiodinium microadriaticum Freudenthal, S. goreauii sp. nov., S. kawagutii sp. nov. and '
                       'S. pilosum sp. nov.: Gymnodinioid dinoflagellate symbionts of marine invertebrates. '
                       'Journal of Phycology 23(3): 469-481.',
        'S. thermophilum' : 'Hume, B. C. C., C. D`Angelo, E. G. Smith, J. R. Stevens, J. Burt and J. Wiedenmann (2015).'
                            ' Symbiodinium thermophilum sp. nov., a thermotolerant symbiotic alga prevalent in corals '
                            'of the world`s hottest sea, the Persian/Arabian Gulf. Sci. Rep. 5.',
        'S. glynnii' : 'LaJeunesse, T. C., D. T. Pettay, E. M. Sampayo, N. Phongsuwan, B. Brown, D. O. Obura, O. '
                       'Hoegh-Guldberg and W. K. Fitt (2010). Long-standing environmental conditions, geographic '
                       'isolation and host-symbiont specificity influence the relative ecological dominance and '
                       'genetic diversification of coral endosymbionts in the genus Symbiodinium. Journal of '
                       'Biogeography 37(5): 785-800.',
        'S. trenchii' : 'LaJeunesse, T. C., D. C. Wham, D. T. Pettay, J. E. Parkinson, S. Keshavmurthy and C. A. Chen '
                        '(2014). Ecologically differentiated stress-tolerant endosymbionts in the dinoflagellate genus'
                        ' Symbiodinium (Dinophyceae) Clade D are different species. Phycologia 53(4): 305-319.',
        'S. eurythalpos' : 'LaJeunesse, T. C., D. C. Wham, D. T. Pettay, J. E. Parkinson, S. Keshavmurthy and C. A. Chen '
                        '(2014). Ecologically differentiated stress-tolerant endosymbionts in the dinoflagellate genus'
                        ' Symbiodinium (Dinophyceae) Clade D are different species. Phycologia 53(4): 305-319.',
        'S. boreum' : 'LaJeunesse, T. C., D. C. Wham, D. T. Pettay, J. E. Parkinson, S. Keshavmurthy and C. A. Chen '
                        '(2014). "Ecologically differentiated stress-tolerant endosymbionts in the dinoflagellate genus'
                        ' Symbiodinium (Dinophyceae) Clade D are different species." Phycologia 53(4): 305-319.',
        'S. voratum' : 'Jeong, H. J., S. Y. Lee, N. S. Kang, Y. D. Yoo, A. S. Lim, M. J. Lee, H. S. Kim, W. Yih, H. '
                       'Yamashita and T. C. LaJeunesse (2014). Genetics and Morphology Characterize the Dinoflagellate'
                       ' Symbiodinium voratum, n. sp., (Dinophyceae) as the Sole Representative of Symbiodinium Clade E'
                       '. Journal of Eukaryotic Microbiology 61(1): 75-94.',
        'S. kawagutii' : 'Trench, R. (2000). Validation of some currently used invalid names of dinoflagellates. '
                       'Journal of Phycology 36(5): 972-972.\tTrench, R. K. and R. J. Blank (1987). '
                       'Symbiodinium microadriaticum Freudenthal, S. goreauii sp. nov., S. kawagutii sp. nov. and '
                       'S. pilosum sp. nov.: Gymnodinioid dinoflagellate symbionts of marine invertebrates. '
                       'Journal of Phycology 23(3): 469-481.'
    }

    species_set = set()
    for analysis_type_obj in across_clade_sorted_type_order:
        if analysis_type_obj.species != 'None':
            species_set.update(analysis_type_obj.species.split(','))

    # put in the species reference title
    outputTableOne.append('Species references')
    outputTableTwo.append('Species references')

    # now add the references for each of the associated species
    for species in species_set:
        if species in species_ref_dict.keys():
            outputTableOne.append(species + '\t' + species_ref_dict[species])
            outputTableTwo.append(species + '\t' + species_ref_dict[species])


    path_to_profiles_absolute = '{}/{}_{}.profiles.absolute.txt'.format(outputDir, analysisObj.id, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    writeListToDestination(path_to_profiles_absolute, outputTableOne)
    output_files_list.append(path_to_profiles_absolute)

    path_to_profiles_rel = '{}/{}_{}.profiles.relative.txt'.format(outputDir, analysisObj.id, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    writeListToDestination(path_to_profiles_rel, outputTableTwo)
    output_files_list.append(path_to_profiles_rel)

    # ########################## ITS2 INTRA ABUND COUNT TABLE ################################
    div_output_pre_analysis_new_meta_and_new_dss_structure(datasubstooutput=datasubstooutput, numProcessors=numProcessors, output_dir=outputDir, sorted_sample_list=samples_that_have_been_sorted, analysis_obj_id=analysisobj.id)

    print('ITS2 type profile output files:')
    for output_file in output_files_list:
        print(output_file)
    #
    # # This is where we're going to have to work with the sequences that aren't part of a type.
    # # Esentially we want to end up with the noName sequences divieded up cladally.
    # # So at the moment this will mean that we will divide up the current no names but also
    # # add information about the other cladal sequences that didn't make it into a cladeCollection
    #
    #
    # # Now we need to output the ITS2 intragenomic abundances that the types are based on for each CC
    # # Already have querySetOfDataSubmissions and cladeCollectionsFromDSs
    #
    # ################ GET ORDERED LIST OF INTRAS BY CLADE THEN ABUNDANCE #################
    # # This is the master list that will hold the eventually ordered defining referenceSequences
    # cladeAbundanceOrderedRefSeqList = []
    # refSeqsInDSs = reference_sequence.objects.filter(data_set_sample_sequence__cladeCollectionTwoFoundIn__in = cladeCollectionsFromDSs).distinct()
    #
    # # Get list of clades found
    # cladeSet = set()
    # for refSeq in refSeqsInDSs:
    #     cladeSet.add(refSeq.clade)
    # # Order the list of clades alphabetically
    # cladeList = [a for a in cladeList if a in cladeSet]
    #
    #
    # for i in range(len(cladeList)):
    #     intraCounter = defaultdict(int)
    #
    #     # refSeqs of clade
    #     refSeqsOfClade = refSeqsInDSs.filter(clade=cladeList[i])
    #
    #     # dataSetSampleSequencesOfClade
    #     dataSetSampleSequencesOfClade = data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn__in = cladeCollectionsFromDSs, referenceSequenceOf__in = refSeqsOfClade)
    #
    #     # Popualate the counter dict
    #
    #     intraCounterManager = Manager()
    #     intraCounter = intraCounterManager.dict()
    #
    #     dsssQueue = Queue()
    #
    #     for dsss in dataSetSampleSequencesOfClade:
    #         dsssQueue.put(dsss)
    #
    #     for N in range(numProcessors):
    #         dsssQueue.put('STOP')
    #
    #     allProcesses = []
    #
    #     # close all connections to the db so that they are automatically recreated for each process
    #     # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    #     db.connections.close_all()
    #
    #     print('Creating refSeq counts for clade {}'.format(cladeList[i]))
    #     for N in range(numProcessors):
    #         p = Process(target=outputWorkerTwo, args=(dsssQueue, intraCounter))
    #         allProcesses.append(p)
    #         p.start()
    #
    #     for p in allProcesses:
    #         p.join()
    #
    #     # for dsss in dataSetSampleSequencesOfClade:
    #     #     intraCounter[dsss.referenceSequenceOf.name] += dsss.abundance
    #
    #     # Sort counter dict by abundance
    #     sortedIntraCounter = sorted(intraCounter.items(), key=operator.itemgetter(1), reverse=True)
    #
    #     # Add the refSeqNames to the master cladeAbundanceOrderedRefSeqList
    #     # Do not add the noName items. These will be summed across the clades and entered into the list
    #     # once their final abundances have been calculated
    #     cladeAbundanceOrderedRefSeqList.extend([a[0] for a in sortedIntraCounter if a[0] != 'noName'])
    #
    # # At this point we have the cladeAbundanceOrderedRefSeqList
    # # Now we just need to add the noName refSeq to the list
    #
    # # cladeAbundanceOrderedRefSeqList.append('noName')
    #
    #
    # # Create the two output tables as lists
    # intraAbundCountTable = []
    # intraAbundPropTable = []
    # # Create the header line one for prop one for count
    # # add the cladal noName categories directly to the header
    # # do not add them to the cladeAbundanceOrderedRefSeqList so that we don't append 0s to them
    # # we will add the abundances in the outputWorkerThree method below
    #
    # headerPre = 'Samples\t{}'.format('\t'.join(cladeAbundanceOrderedRefSeqList))
    # noNameStrings = '\t'.join(['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']])
    # header = '\t'.join([headerPre, noNameStrings])
    #
    #
    # intraAbundCountTable.append(header)
    # intraAbundPropTable.append(header)
    # ######################################################################################
    #
    # ############ POPULATE TABLE WITH CC DATA #############
    # # For every sample
    # sampleList = data_set_sample.objects.filter(dataSubmissionFrom__in=querySetOfDataSubmissions)
    #
    #
    # # In order to MP this we will have to pay attention to the order. As we can't keep the order as we work with the
    # # MPs we will do as we did above for the profies outputs and have an output dict that will have the sample as key
    # # and its corresponding data row as the value. Then when we have finished the MPing we will go through the
    # # sampleList order and output the data in this order.
    # # So we will need a managed output dict.
    #
    # sampleOutputDictManager = Manager()
    # managedSampleOutputDict = sampleOutputDictManager.dict()
    #
    # cladeAbundanceOrderedRefSeqListManager = Manager()
    # cladeAbundanceOrderedRefSeqList = cladeAbundanceOrderedRefSeqListManager.list(cladeAbundanceOrderedRefSeqList)
    #
    #
    #
    # dssQueue = Queue()
    #
    # for dss in sampleList:
    #     dssQueue.put(dss)
    #
    # for N in range(numProcessors):
    #     dssQueue.put('STOP')
    #
    # allProcesses = []
    #
    # # close all connections to the db so that they are automatically recreated for each process
    # # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    # db.connections.close_all()
    #
    # for N in range(numProcessors):
    #     p = Process(target=outputWorkerThree, args=(dssQueue, managedSampleOutputDict, cladeAbundanceOrderedRefSeqList, querySetOfDataSubmissions))
    #     allProcesses.append(p)
    #     p.start()
    #
    # for p in allProcesses:
    #     p.join()
    #
    #
    # # Now before writing out we need to add to the intraAbund objects from the managedSampleOutputDict in order of
    # # sampleList
    # #### DEVELOPMENT ####
    # # So that we can inspect the managedSampleOutputDict
    # managedSampleOutputDict_dict = dict(managedSampleOutputDict)
    #
    #
    # #####################
    #
    #
    # for dss in samples_that_have_been_sorted:
    #     intraAbundCountTable.append(managedSampleOutputDict[dss][0])
    #     intraAbundPropTable.append(managedSampleOutputDict[dss][1])
    #
    #
    #
    # # Now add the accesion number / UID for each of the DIVs
    # accession_list = []
    # accession_list.append('DIV_accession')
    # for ref_seq_name in cladeAbundanceOrderedRefSeqList:
    #     ref_seq = reference_sequence.objects.get(name=ref_seq_name)
    #     ref_seq_accession = ''
    #     if ref_seq.accession and ref_seq.accession != 'None':
    #         ref_seq_accession = ref_seq.accession
    #     else:
    #         ref_seq_accession = str(ref_seq.id)
    #     accession_list.append(ref_seq_accession)
    #
    # intraAbundCountTable.append('\t'.join(accession_list))
    # intraAbundPropTable.append('\t'.join(accession_list))
    #
    # # Here we have the tables populated and ready to output
    # path_to_DIV_absolute = '{}/{}_{}.DIVs.absolute.txt'.format(outputDir, analysisObj.id, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    # writeListToDestination(path_to_DIV_absolute, intraAbundCountTable)
    # output_files_list.append(path_to_DIV_absolute)
    #
    # path_to_DIV_rel = '{}/{}_{}.DIVs.relative.txt'.format(outputDir, analysisObj.id, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    # writeListToDestination(path_to_DIV_rel, intraAbundPropTable)
    # output_files_list.append(path_to_DIV_rel)
    #
    # # Output a .fasta for of all of the sequences found in the analysis
    # fasta_output_list = []
    # for ref_seq_name in cladeAbundanceOrderedRefSeqList:
    #     ref_seq = reference_sequence.objects.get(name=ref_seq_name)
    #     fasta_output_list.append('>{}'.format(ref_seq_name))
    #     fasta_output_list.append(ref_seq.sequence)
    #
    # path_to_fasta = '{}/{}_{}.DIVs.fasta'.format(outputDir, analysisObj.id, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    # writeListToDestination(path_to_fasta, fasta_output_list)
    # output_files_list.append(path_to_fasta)
    #
    # path_to_sorted_list = '{}/{}_{}_sorted_sample_list.txt'.format(outputDir, analysisObj.id, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    # writeListToDestination(path_to_sorted_list, [str(a.id) for a in samples_that_have_been_sorted])
    # output_files_list.append(path_to_sorted_list)
    #
    # print('Output files:')
    # for output_file in output_files_list:
    #     print(output_file)

    return

def outputWorkerThree(input, outDict, cladeAbundanceOrderedRefSeqList, querySetOfDataSubmissions):
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    for dss in iter(input.get, 'STOP'):


        print('Outputting DIV data for {} with {}'.format(dss.name, current_process().name))
        # List that will hold the row
        sampleRowDataCounts = []
        sampleRowDataProps = []
        cladalAbundances = [int(a) for a in json.loads(dss.cladalSeqTotals)]
        sampleSeqTot = sum(cladalAbundances)

        # TODO investigate. THe particular case was a sample called T1P1 from one of alejandros datasets
        # For some reason a sample got through all of the error catching with 0 sequence in it.
        # I'm not sure how this is happening but for the time being I will just be sure to write it out as
        # 0s across the board.
        if dss.errorInProcessing or sampleSeqTot == 0:
            # Then this sample had a problem in the sequencing and we need to just output 0s across the board
            for name in cladeAbundanceOrderedRefSeqList:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            for cl in cladeList:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # Here we need to add the string to the outputDict rather than the intraAbund table objects
            countData = '{}\t{}'.format(dss.name, '\t'.join(str(x) for x in sampleRowDataCounts))
            propData = '{}\t{}'.format(dss.name, '\t'.join("{:.3f}".format(x) for x in sampleRowDataProps))
            outDict[dss] = [countData, propData]
            continue

        # CCs of sample
        cladeCollectionsFromDSs = clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom__in=querySetOfDataSubmissions)
        sampleCCs = cladeCollectionsFromDSs.filter(dataSetSampleFrom=dss)

        # Get number of seqs found in Sample
        dsssInSample = data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn__in=sampleCCs)







        # sampleSeqTot = sum([dsss.abundance for dsss in dsssInSample])

        # create dataSetSampleSeqDict for sample

        dataSetSampleSeqDict = defaultdict(int)
        # this dict will hold noName seqs abundances separated by clade.
        # when we come to populate the info on the sequences that were not in clade collections then we will
        # check to see if there is info for the clade in question in this dict first.
        # if the clade in question isn't in this dict then we will use the information from the
        # samples cladalSeqTotals parameter
        noNameCladalHolder = defaultdict(int)
        for dsss in dsssInSample:
            if dsss.referenceSequenceOf.name != 'noName':
                dataSetSampleSeqDict[dsss.referenceSequenceOf.name] += dsss.abundance
            elif dsss.referenceSequenceOf.name == 'noName':
                noNameCladalHolder[dsss.referenceSequenceOf.clade] += dsss.abundance
        # populate row in order of the master sorted intras
        # Do not populate any of the noName categories

        # #### DEBUG ####
        # # We seem to be creating some rows of data that don't add up to one so I want to verify if this is being caused
        # # due to the sampleSeqTot mismatching with the sum of the dsss.abundances
        # this = sum(dict(dataSetSampleSeqDict).values()) + sum(dict(noNameCladalHolder).values())
        # nonCC_cladalTot = 0
        # for clade in cladeList:
        #     if clade not in [a.clade for a in sampleCCs]:
        #         nonCC_cladalTot += cladalAbundances[cladeList.index(clade)]
        #
        # if sampleSeqTot != (sum(dataSetSampleSeqDict.values()) + sum(noNameCladalHolder.values()) + nonCC_cladalTot):
        #
        #     apples = 'peras'
        #     bob = list(dsssInSample)
        #     harry = list(sampleCCs)
        #     holder = []
        #     for d in dsssInSample:
        #         holder.append((d.referenceSequenceOf.name, d.abundance, d.referenceSequenceOf.clade))

        for name in cladeAbundanceOrderedRefSeqList:
            if name in dataSetSampleSeqDict.keys():
                counts = dataSetSampleSeqDict[name]
                sampleRowDataCounts.append(counts)
                sampleRowDataProps.append(counts / sampleSeqTot)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)

        # here we can make the data for the noName cladal divisions
        # We were having a problem here
        list_of_clades_in_CCs = [a.clade for a in sampleCCs]
        for clade in cladeList:
            if clade in noNameCladalHolder.keys():
                # then there was a cladeCollection of this clade AND there was noName seqs
                # so we can add this info
                sampleRowDataCounts.append(noNameCladalHolder[clade])
                sampleRowDataProps.append(noNameCladalHolder[clade]/sampleSeqTot)
            elif clade not in list_of_clades_in_CCs:
                # then there were no cladeCollections of this clade. There may be no sequences of this clade found
                # in the sample. Either way the seqs of this clade found will be in the cladalSeqTotals parameter
                # of the sample
                sampleRowDataCounts.append(cladalAbundances[cladeList.index(clade)])
                sampleRowDataProps.append(cladalAbundances[cladeList.index(clade)]/sampleSeqTot)
            elif clade in list_of_clades_in_CCs:
                # we need this eventuality for when we have a clade collection but no 'noName' seqs. Only
                # seqs that are named. In this case we will not enter either of the above ifs and we will need
                # to enter a 0 here.
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)


        # Here we need to add the string to the outputDict rather than the intraAbund table objects
        countData = '{}\t{}'.format(dss.name, '\t'.join(str(x) for x in sampleRowDataCounts))
        propData = '{}\t{}'.format(dss.name, '\t'.join("{:.3f}".format(x) for x in sampleRowDataProps))
        outDict[dss] = [countData, propData]


def outputWorkerThree_pre_analysis(input, outDict, cladeAbundanceOrderedRefSeqList, querySetOfDataSubmissions):
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    for dss in iter(input.get, 'STOP'):

        print('Outputting DIV data for {} with {}'.format(dss.name, current_process().name))
        # List that will hold the row
        sampleRowDataCounts = []
        sampleRowDataProps = []
        cladalAbundances = [int(a) for a in json.loads(dss.cladalSeqTotals)]
        sampleSeqTot = sum(cladalAbundances)

        if dss.errorInProcessing or sampleSeqTot == 0:
            # Then this sample had a problem in the sequencing and we need to just output 0s across the board
            # Named sequences get 0s
            for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # Clade totals get 0s
            for cl in cladeList:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # Add in the qc totals if possible
            # For the proportions we will have to add zeros as we cannot do proportions
            # CONTIGS
            # This is the absolute number of sequences after make.contigs
            if dss.initialTotSeqNum:
                contig_num = dss.initialTotSeqNum
                sampleRowDataCounts.append(contig_num)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # POST-QC
            # store the aboslute number of sequences after sequencing QC at this stage
            if dss.post_seq_qc_absolute_num_seqs:
                post_qc_absolute = dss.post_seq_qc_absolute_num_seqs
                sampleRowDataCounts.append(post_qc_absolute)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # This is the unique number of sequences after the sequencing QC
            if dss.initialUniqueSeqNum:
                post_qc_unique = dss.initialUniqueSeqNum
                sampleRowDataCounts.append(post_qc_unique)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # POST TAXA-ID
            # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
            if dss.finalTotSeqNum:
                tax_id_symbiodinium_absolute = dss.finalTotSeqNum
                sampleRowDataCounts.append(tax_id_symbiodinium_absolute)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # Same as above but the number of unique seqs
            if dss.finalUniqueSeqNum:
                tax_id_symbiodinium_unique = dss.finalUniqueSeqNum
                sampleRowDataCounts.append(tax_id_symbiodinium_unique)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # size violation absolute
            if dss.size_violation_absolute:
                size_viol_ab = dss.size_violation_absolute
                sampleRowDataCounts.append(size_viol_ab)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # size violation unique
            if dss.size_violation_unique:
                size_viol_uni = dss.size_violation_unique
                sampleRowDataCounts.append(size_viol_uni)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # store the abosolute number of sequenes that were not considered Symbiodinium
            if dss.non_sym_absolute_num_seqs:
                tax_id_non_symbiodinum_abosulte = dss.non_sym_absolute_num_seqs
                sampleRowDataCounts.append(tax_id_non_symbiodinum_abosulte)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # This is the number of unique sequences that were not considered Symbiodinium
            if dss.nonSymSeqsNum:
                tax_id_non_symbiodinium_unique = dss.nonSymSeqsNum
                sampleRowDataCounts.append(tax_id_non_symbiodinium_unique)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # post-med absolute
            if dss.post_med_absolute:
                post_med_abs = dss.post_med_absolute
                sampleRowDataCounts.append(post_med_abs)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # post-med absolute
            if dss.post_med_unique:
                post_med_uni = dss.post_med_unique
                sampleRowDataCounts.append(post_med_uni)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)


            # Breakdown details get 0s
            for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # Here we need to add the string to the outputDict rather than the intraAbund table objects
            countData = '{}\t{}'.format(dss.name, '\t'.join(str(x) for x in sampleRowDataCounts))
            propData = '{}\t{}'.format(dss.name, '\t'.join("{:.3f}".format(x) for x in sampleRowDataProps))
            outDict[dss] = [countData, propData]
            continue

        # CCs of sample
        cladeCollectionsFromDSs = clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom__in=querySetOfDataSubmissions)
        sampleCCs = cladeCollectionsFromDSs.filter(dataSetSampleFrom=dss)

        # Get number of seqs found in Sample
        dsssInSample = data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn__in=sampleCCs)


        # sampleSeqTot = sum([dsss.abundance for dsss in dsssInSample])

        # create dataSetSampleSeqDict for sample

        dataSetSampleSeqDict = defaultdict(int)
        # this dict will hold noName seqs abundances separated by clade.
        # when we come to populate the info on the sequences that were not in clade collections then we will
        # check to see if there is info for the clade in question in this dict first.
        # if the clade in question isn't in this dict then we will use the information from the
        # samples cladalSeqTotals parameter
        noNameCladalHolder = defaultdict(int)
        no_name_break_down_dict = defaultdict(int)
        for dsss in dsssInSample:
            ref_seq_obj = dsss.referenceSequenceOf
            ref_seq_obj_name = ref_seq_obj.name
            ref_seq_obj_clade = ref_seq_obj.clade
            if ref_seq_obj_name != 'noName':
                dataSetSampleSeqDict[ref_seq_obj_name] += dsss.abundance
            elif ref_seq_obj_name == 'noName':
                noNameCladalHolder[ref_seq_obj_clade] += dsss.abundance
                no_name_break_down_dict[str(ref_seq_obj.id) + '_{}'.format(ref_seq_obj_clade)] += dsss.abundance        # populate row in order of the master sorted intras
        # Do not populate any of the noName categories


        for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
            if name in dataSetSampleSeqDict.keys():
                counts = dataSetSampleSeqDict[name]
                sampleRowDataCounts.append(counts)
                sampleRowDataProps.append(counts / sampleSeqTot)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)

        # here we can make the data for the noName cladal divisions
        # We were having a problem here
        list_of_clades_in_CCs = [a.clade for a in sampleCCs]
        for clade in cladeList:
            if clade in noNameCladalHolder.keys():
                # then there was a cladeCollection of this clade AND there was noName seqs
                # so we can add this info
                sampleRowDataCounts.append(noNameCladalHolder[clade])
                sampleRowDataProps.append(noNameCladalHolder[clade]/sampleSeqTot)
            elif clade not in list_of_clades_in_CCs:
                # then there were no cladeCollections of this clade. There may be no sequences of this clade found
                # in the sample. Either way the seqs of this clade found will be in the cladalSeqTotals parameter
                # of the sample
                sampleRowDataCounts.append(cladalAbundances[cladeList.index(clade)])
                sampleRowDataProps.append(cladalAbundances[cladeList.index(clade)]/sampleSeqTot)
            elif clade in list_of_clades_in_CCs:
                # we need this eventuality for when we have a clade collection but no 'noName' seqs. Only
                # seqs that are named. In this case we will not enter either of the above ifs and we will need
                # to enter a 0 here.
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)

        # Here we add in the post qc and post-taxa id counts
        # For the absolute counts we will report the absolute seq number
        # For the realteive counts we will report these as proportions of the sampleSeqTot.
        # I.e. we will have numbers larger than 1 for many of the values and the symbiodinium seqs should be 1

        # CONTIGS
        # This is the absolute number of sequences after make.contigs
        contig_num = dss.initialTotSeqNum
        sampleRowDataCounts.append(contig_num)
        sampleRowDataProps.append(contig_num/sampleSeqTot)

        # POST-QC
        # store the aboslute number of sequences after sequencing QC at this stage
        post_qc_absolute = dss.post_seq_qc_absolute_num_seqs
        sampleRowDataCounts.append(post_qc_absolute)
        sampleRowDataProps.append(post_qc_absolute / sampleSeqTot)
        # This is the unique number of sequences after the sequencing QC
        post_qc_unique = dss.initialUniqueSeqNum
        sampleRowDataCounts.append(post_qc_unique)
        sampleRowDataProps.append(post_qc_unique / sampleSeqTot)

        # POST TAXA-ID
        # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
        tax_id_symbiodinium_absolute = dss.finalTotSeqNum
        sampleRowDataCounts.append(tax_id_symbiodinium_absolute)
        sampleRowDataProps.append(tax_id_symbiodinium_absolute / sampleSeqTot)
        # Same as above but the number of unique seqs
        tax_id_symbiodinium_unique = dss.finalUniqueSeqNum
        sampleRowDataCounts.append(tax_id_symbiodinium_unique)
        sampleRowDataProps.append(tax_id_symbiodinium_unique / sampleSeqTot)
        # store the absolute number of sequences lost to size cutoff violations
        size_violation_aboslute = dss.size_violation_absolute
        sampleRowDataCounts.append(size_violation_aboslute)
        sampleRowDataProps.append(size_violation_aboslute / sampleSeqTot)
        # store the unique size cutoff violations
        size_violation_unique = dss.size_violation_unique
        sampleRowDataCounts.append(size_violation_unique)
        sampleRowDataProps.append(size_violation_unique / sampleSeqTot)
        # store the abosolute number of sequenes that were not considered Symbiodinium
        tax_id_non_symbiodinum_abosulte = dss.non_sym_absolute_num_seqs
        sampleRowDataCounts.append(tax_id_non_symbiodinum_abosulte)
        sampleRowDataProps.append(tax_id_non_symbiodinum_abosulte / sampleSeqTot)
        # This is the number of unique sequences that were not considered Symbiodinium
        tax_id_non_symbiodinium_unique = dss.nonSymSeqsNum
        sampleRowDataCounts.append(tax_id_non_symbiodinium_unique)
        sampleRowDataProps.append(tax_id_non_symbiodinium_unique / sampleSeqTot)

        # Post MED absolute
        post_med_absolute = dss.post_med_absolute
        sampleRowDataCounts.append(post_med_absolute)
        sampleRowDataProps.append(post_med_absolute / sampleSeqTot)
        # Post MED unique
        post_med_unique = dss.post_med_unique
        sampleRowDataCounts.append(post_med_unique)
        sampleRowDataProps.append(post_med_unique / sampleSeqTot)

        # now add the noName clade breakdown sequence counts
        for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
            if name in no_name_break_down_dict.keys():
                counts = no_name_break_down_dict[name]
                sampleRowDataCounts.append(counts)
                sampleRowDataProps.append(counts / sampleSeqTot)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)



        # Here we need to add the string to the outputDict rather than the intraAbund table objects
        countData = '{}\t{}'.format(dss.name, '\t'.join(str(x) for x in sampleRowDataCounts))
        propData = '{}\t{}'.format(dss.name, '\t'.join("{:.3f}".format(x) for x in sampleRowDataProps))
        outDict[dss] = [countData, propData]


def outputWorkerThree_pre_analysis_new_dss_structure(input, outDict, cladeAbundanceOrderedRefSeqList, querySetOfDataSubmissions):
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    for dss in iter(input.get, 'STOP'):

        print('Outputting DIV data for {} with {}'.format(dss.name, current_process().name))
        # List that will hold the row
        sampleRowDataCounts = []
        sampleRowDataProps = []
        cladalAbundances = [int(a) for a in json.loads(dss.cladalSeqTotals)]
        sampleSeqTot = sum(cladalAbundances)

        if dss.errorInProcessing or sampleSeqTot == 0:
            # Then this sample had a problem in the sequencing and we need to just output 0s across the board
            # Named sequences get 0s
            for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # Clade totals get 0s
            for cl in cladeList:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # Add in the qc totals if possible
            # For the proportions we will have to add zeros as we cannot do proportions
            # CONTIGS
            # This is the absolute number of sequences after make.contigs
            if dss.initialTotSeqNum:
                contig_num = dss.initialTotSeqNum
                sampleRowDataCounts.append(contig_num)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # POST-QC
            # store the aboslute number of sequences after sequencing QC at this stage
            if dss.post_seq_qc_absolute_num_seqs:
                post_qc_absolute = dss.post_seq_qc_absolute_num_seqs
                sampleRowDataCounts.append(post_qc_absolute)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # This is the unique number of sequences after the sequencing QC
            if dss.initialUniqueSeqNum:
                post_qc_unique = dss.initialUniqueSeqNum
                sampleRowDataCounts.append(post_qc_unique)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # POST TAXA-ID
            # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
            if dss.finalTotSeqNum:
                tax_id_symbiodinium_absolute = dss.finalTotSeqNum
                sampleRowDataCounts.append(tax_id_symbiodinium_absolute)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # Same as above but the number of unique seqs
            if dss.finalUniqueSeqNum:
                tax_id_symbiodinium_unique = dss.finalUniqueSeqNum
                sampleRowDataCounts.append(tax_id_symbiodinium_unique)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # size violation absolute
            if dss.size_violation_absolute:
                size_viol_ab = dss.size_violation_absolute
                sampleRowDataCounts.append(size_viol_ab)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # size violation unique
            if dss.size_violation_unique:
                size_viol_uni = dss.size_violation_unique
                sampleRowDataCounts.append(size_viol_uni)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # store the abosolute number of sequenes that were not considered Symbiodinium
            if dss.non_sym_absolute_num_seqs:
                tax_id_non_symbiodinum_abosulte = dss.non_sym_absolute_num_seqs
                sampleRowDataCounts.append(tax_id_non_symbiodinum_abosulte)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # This is the number of unique sequences that were not considered Symbiodinium
            if dss.nonSymSeqsNum:
                tax_id_non_symbiodinium_unique = dss.nonSymSeqsNum
                sampleRowDataCounts.append(tax_id_non_symbiodinium_unique)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # post-med absolute
            if dss.post_med_absolute:
                post_med_abs = dss.post_med_absolute
                sampleRowDataCounts.append(post_med_abs)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # post-med absolute
            if dss.post_med_unique:
                post_med_uni = dss.post_med_unique
                sampleRowDataCounts.append(post_med_uni)
                sampleRowDataProps.append(0)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)

            # Breakdown details get 0s
            for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # Here we need to add the string to the outputDict rather than the intraAbund table objects
            countData = '{}\t{}'.format(dss.name, '\t'.join(str(x) for x in sampleRowDataCounts))
            propData = '{}\t{}'.format(dss.name, '\t'.join("{:.3f}".format(x) for x in sampleRowDataProps))
            outDict[dss] = [countData, propData]
            continue

        # CCs of sample
        # cladeCollectionsFromDSs = clade_collection.objects.filter(
        #     dataSetSampleFrom__dataSubmissionFrom__in=querySetOfDataSubmissions)
        # sampleCCs = cladeCollectionsFromDSs.filter(dataSetSampleFrom=dss)

        # Get number of seqs found in Sample
        dsssInSample = data_set_sample_sequence.objects.filter(data_set_sample_from=dss)

        # sampleSeqTot = sum([dsss.abundance for dsss in dsssInSample])

        # create dataSetSampleSeqDict for sample

        dataSetSampleSeqDict = defaultdict(int)
        # this dict will hold noName seqs abundances separated by clade.
        # when we come to populate the info on the sequences that were not in clade collections then we will
        # check to see if there is info for the clade in question in this dict first.
        # if the clade in question isn't in this dict then we will use the information from the
        # samples cladalSeqTotals parameter
        noNameCladalHolder = defaultdict(int)
        no_name_break_down_dict = defaultdict(int)
        for dsss in dsssInSample:
            ref_seq_obj = dsss.referenceSequenceOf
            ref_seq_obj_name = ref_seq_obj.name
            ref_seq_obj_clade = ref_seq_obj.clade
            if ref_seq_obj_name != 'noName':
                dataSetSampleSeqDict[ref_seq_obj_name] += dsss.abundance
            elif ref_seq_obj_name == 'noName':
                noNameCladalHolder[ref_seq_obj_clade] += dsss.abundance
                no_name_break_down_dict[str(ref_seq_obj.id) + '_{}'.format(
                    ref_seq_obj_clade)] += dsss.abundance  # populate row in order of the master sorted intras
        # Do not populate any of the noName categories


        for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
            if name in dataSetSampleSeqDict.keys():
                counts = dataSetSampleSeqDict[name]
                sampleRowDataCounts.append(counts)
                sampleRowDataProps.append(counts / sampleSeqTot)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)

        # here we can make the data for the noName cladal divisions
        # We were having a problem here

        # for each clade in cladeList, just look to see if there is a noName.
        # sample_CCs = clade_collection.objects.filter(dataSetSampleFrom=dss)
        # list_of_clades_in_CCs = [a.clade for a in sample_CCs]
        for clade in cladeList:
            if clade in noNameCladalHolder.keys():
                sampleRowDataCounts.append(noNameCladalHolder[clade])
                sampleRowDataProps.append(noNameCladalHolder[clade] / sampleSeqTot)
            else:
                # then there were no no_names of this clade
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # if clade in noNameCladalHolder.keys():
            #     # then there was a cladeCollection of this clade AND there was noName seqs
            #     # so we can add this info
            #     sampleRowDataCounts.append(noNameCladalHolder[clade])
            #     sampleRowDataProps.append(noNameCladalHolder[clade] / sampleSeqTot)
            # elif clade not in list_of_clades_in_CCs:
            #     # then there were no cladeCollections of this clade. There may be no sequences of this clade found
            #     # in the sample. Either way the seqs of this clade found will be in the cladalSeqTotals parameter
            #     # of the sample
            #     sampleRowDataCounts.append(cladalAbundances[cladeList.index(clade)])
            #     sampleRowDataProps.append(cladalAbundances[cladeList.index(clade)] / sampleSeqTot)
            # elif clade in list_of_clades_in_CCs:
            #     # we need this eventuality for when we have a clade collection but no 'noName' seqs. Only
            #     # seqs that are named. In this case we will not enter either of the above ifs and we will need
            #     # to enter a 0 here.
            #     sampleRowDataCounts.append(0)
            #     sampleRowDataProps.append(0)

        # Here we add in the post qc and post-taxa id counts
        # For the absolute counts we will report the absolute seq number
        # For the realteive counts we will report these as proportions of the sampleSeqTot.
        # I.e. we will have numbers larger than 1 for many of the values and the symbiodinium seqs should be 1

        # CONTIGS
        # This is the absolute number of sequences after make.contigs
        contig_num = dss.initialTotSeqNum
        sampleRowDataCounts.append(contig_num)
        sampleRowDataProps.append(contig_num / sampleSeqTot)

        # POST-QC
        # store the aboslute number of sequences after sequencing QC at this stage
        post_qc_absolute = dss.post_seq_qc_absolute_num_seqs
        sampleRowDataCounts.append(post_qc_absolute)
        sampleRowDataProps.append(post_qc_absolute / sampleSeqTot)
        # This is the unique number of sequences after the sequencing QC
        post_qc_unique = dss.initialUniqueSeqNum
        sampleRowDataCounts.append(post_qc_unique)
        sampleRowDataProps.append(post_qc_unique / sampleSeqTot)

        # POST TAXA-ID
        # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
        tax_id_symbiodinium_absolute = dss.finalTotSeqNum
        sampleRowDataCounts.append(tax_id_symbiodinium_absolute)
        sampleRowDataProps.append(tax_id_symbiodinium_absolute / sampleSeqTot)
        # Same as above but the number of unique seqs
        tax_id_symbiodinium_unique = dss.finalUniqueSeqNum
        sampleRowDataCounts.append(tax_id_symbiodinium_unique)
        sampleRowDataProps.append(tax_id_symbiodinium_unique / sampleSeqTot)
        # store the absolute number of sequences lost to size cutoff violations
        size_violation_aboslute = dss.size_violation_absolute
        sampleRowDataCounts.append(size_violation_aboslute)
        sampleRowDataProps.append(size_violation_aboslute / sampleSeqTot)
        # store the unique size cutoff violations
        size_violation_unique = dss.size_violation_unique
        sampleRowDataCounts.append(size_violation_unique)
        sampleRowDataProps.append(size_violation_unique / sampleSeqTot)
        # store the abosolute number of sequenes that were not considered Symbiodinium
        tax_id_non_symbiodinum_abosulte = dss.non_sym_absolute_num_seqs
        sampleRowDataCounts.append(tax_id_non_symbiodinum_abosulte)
        sampleRowDataProps.append(tax_id_non_symbiodinum_abosulte / sampleSeqTot)
        # This is the number of unique sequences that were not considered Symbiodinium
        tax_id_non_symbiodinium_unique = dss.nonSymSeqsNum
        sampleRowDataCounts.append(tax_id_non_symbiodinium_unique)
        sampleRowDataProps.append(tax_id_non_symbiodinium_unique / sampleSeqTot)

        # Post MED absolute
        post_med_absolute = dss.post_med_absolute
        sampleRowDataCounts.append(post_med_absolute)
        sampleRowDataProps.append(post_med_absolute / sampleSeqTot)
        # Post MED unique
        post_med_unique = dss.post_med_unique
        sampleRowDataCounts.append(post_med_unique)
        sampleRowDataProps.append(post_med_unique / sampleSeqTot)

        # now add the noName clade breakdown sequence counts
        for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
            if name in no_name_break_down_dict.keys():
                counts = no_name_break_down_dict[name]
                sampleRowDataCounts.append(counts)
                sampleRowDataProps.append(counts / sampleSeqTot)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)

        # Here we need to add the string to the outputDict rather than the intraAbund table objects
        countData = '{}\t{}'.format(dss.name, '\t'.join(str(x) for x in sampleRowDataCounts))
        propData = '{}\t{}'.format(dss.name, '\t'.join("{:.3f}".format(x) for x in sampleRowDataProps))
        outDict[dss] = [countData, propData]

def outputWorkerTwo_pre_analysis(input, outDict):
    for dsss in iter(input.get, 'STOP'):
        if dsss.referenceSequenceOf.name == 'noName':
            name_unit = str(dsss.referenceSequenceOf.id) + '_{}'.format(dsss.referenceSequenceOf.clade)
        else:
            name_unit = dsss.referenceSequenceOf.name
        if name_unit in outDict.keys():
            outDict[name_unit] = outDict[name_unit] + dsss.abundance
        else:
            outDict[name_unit] = dsss.abundance


def outputWorkerTwo(input, outDict):
    for dsss in iter(input.get, 'STOP'):
        if dsss.referenceSequenceOf.name in outDict.keys():
            outDict[dsss.referenceSequenceOf.name] = outDict[dsss.referenceSequenceOf.name] + dsss.abundance
        else:
            outDict[dsss.referenceSequenceOf.name] = dsss.abundance



def outputWorkerOne(input, groupDict, groupCount, listOfDataSetSamples, clade, outputDict):
    for anType in iter(input.get, 'STOP'):  # Within each type go through each of the samples
        print('Outputing data for\ttype: {} with {}'.format(anType, current_process().name))

        ###### CALCULATE REL ABUND AND SD OF DEF INTRAS FOR THIS TYPE ###############
        # For each type we want to calculate the average proportions of the defining seqs in that type
        # We will also calculate an SD.
        # We will do both of these calculations using the footprintSeqAbundances, listOfCladeCollections
        # and orderedFoorprintList attributes of each type
        footprintAbundances = json.loads(anType.footprintSeqAbundances)
        orderedFootPrintlist = anType.orderedFootprintList.split(',')
        # We have to make a decision as to whether this average should represent all the findings of this type
        # or whether we should only represent the averages of samples found in this dataset.
        # I think it needs to be a global average. Because the type is defined based on all samples in the
        # SymPortal db.
        averageList = [[] for defseq in orderedFootPrintlist]
        for a in range(len(footprintAbundances)):
            for b in range(len(orderedFootPrintlist)):
                averageList[b].append(float(int(footprintAbundances[a][b]) / sum(footprintAbundances[a])))
        # Here we have a list for each def seq (as proportions), this can now have totals, SDs and averages calculated on it
        # then convert to proportions by dividing by total
        SDList = []
        TotList = []

        for a in range(len(averageList)):
            SDList.append(statistics.pstdev(averageList[a]))
            TotList.append(sum(averageList[a]))
        tot = sum(TotList)
        for a in range(len(TotList)):
            TotList[a] = TotList[a] / tot

        # The TotList now contains the proportions of each def seq
        # The SDList now contains the SDs for each of the def seqs proportions
        ###########################################################################



        # Counter that will increment with every sample type is found in
        # This is counter will end up being the type abund local value
        abundanceCount = 0

        typeInQ = anType
        clade = typeInQ.clade

        ###### GET GROUP NAME INFO ########
        groupID = typeInQ.analysisGroupOf.id
        if groupID not in groupDict.keys():
            # with groupCount.get_lock():
            groupCount.value += 1

            groupDict[groupID] = groupCount.value
        ###################################


        # For each type create a holder that will hold 0s for each sample until populated below
        dataRowRaw = [0 for smpl in listOfDataSetSamples]
        dataRowProp = [0 for smpl in listOfDataSetSamples]
        typeCCIDs = [int(a) for a in typeInQ.listOfCladeCollections.split(',')]

        # We also need the type abund db value. We can get this from the type cladeCollections
        globalCount = len(typeCCIDs)

        for k in range(len(listOfDataSetSamples)):  # Within each type go through each of the samples
            # when looking at sample, get the clade collections found in the sample
            CCsInSample = clade_collection.objects.filter(dataSetSampleFrom=listOfDataSetSamples[k])
            # because we are within a clade here the type can only have one cladecollection in common with the sample
            # so check to see if there are any items in common between the type's CC list and the samples CC
            ccInSampleIDs = [a.id for a in CCsInSample]

            IDsInCommon = list(set(ccInSampleIDs).intersection(typeCCIDs))
            if IDsInCommon:  # Then this type was found in this sample
                abundanceCount += 1
                # Need to work out how many seqs were found from the sample for this type
                # Also need to work this out as a proportion of all of the Symbiodinium seqs found in the sample
                ccInQ = clade_collection.objects.get(id=IDsInCommon[0])
                # The number of sequences that were returned for the sample in question
                totSeqNum = data_set_sample.objects.get(clade_collection=ccInQ).finalTotSeqNum
                # The number of sequences that make up the type in q
                # The index of the CC in the type's listOfCladeCollections
                CCindexInType = typeCCIDs.index(ccInQ.id)
                # List containing the abundances of each of the ref seqs that make up the type in the given CC
                seqAbundanceInfoForCCandTypeInQ = json.loads(typeInQ.footprintSeqAbundances)[CCindexInType]
                # total of the above list, i.e. seqs in type
                sumOfDefRefSeqAbundInType = sum(seqAbundanceInfoForCCandTypeInQ)
                # type abundance as proportion of the total seqs found in the sample
                typeProp = sumOfDefRefSeqAbundInType / totSeqNum
                # Now populate the dataRow with the sumOfDefRefSeqAbundInType and the typeProp
                dataRowRaw[k] = sumOfDefRefSeqAbundInType
                dataRowProp[k] = typeProp

        # Type Profile
        typeID = anType.id
        speciesGrp = '{}_Grp_{}'.format(clade, groupDict[groupID])
        typeProfile = typeInQ.name
        species = typeInQ.species
        typeAbundance = abundanceCount
        # This is currently putting out the Majs in an odd order due to the fact that the MajRefSeqSet is
        # a set. Instead we should get the Majs from the name.
        # majITS2 = '/'.join([str(refseq) for refseq in sortedListOfTypes[j].getMajRefSeqSet()])
        majITS2, majList = getMajList(anType)
        typeAbundAndSDString = getAbundStr(TotList, SDList, majList)

        sequenceAccession = typeInQ.generateName(accession=True)
        rowOne = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(typeID, clade,  majITS2,
                                                                     species, str(typeAbundance), str(globalCount), typeProfile,
                                                                     '\t'.join([str(a) for a in dataRowRaw]),
                                                                     sequenceAccession, typeAbundAndSDString)
        rowTwo = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(typeID, clade, majITS2,
                                                                     species, str(typeAbundance), str(globalCount), typeProfile,
                                                                     '\t'.join(["{:.3f}".format(prop) for prop in
                                                                                dataRowProp]), sequenceAccession,
                                                                     typeAbundAndSDString)
        # Here we will use the outputDict instead of adding it to the outputTableOne
        # We will then get these elements in the dict into the right order with the types we will then
        # write out each of the type elements in the order of the orderedTypeLists

        outputDict[anType]  = [rowOne, rowTwo]

        # outputTableOne.append(rowOne)
        # outputTableTwo.append(rowTwo)

def getMajList(atype):
    if str(atype) in ['D1/D2-D4', 'D1/D6/D2.2-D4', 'B5/B5e-B5a-B5b']:
        apples = 'asdf'
    # This is a little tricky. Have to take into account that the Maj seqs are not always more abundant
    # than the non-Maj seqs.
    #e.g. type B5/B5e-B5a-B5b, B5e is acutally the least abundant of the seqs
    # Therefore we cannot simply grab the orderedFooprintList seqs in order assuming they are the Majs

    name = atype.name
    count = name.count('/')
    majList = []
    seqsInOrderOfAbunIDs = atype.orderedFootprintList.split(',')
    majSeqsIDs = atype.MajRefSeqSet.split(',')
    for index in range(count + 1):
        for item in range(len(seqsInOrderOfAbunIDs)):
            if seqsInOrderOfAbunIDs[item] in majSeqsIDs:
                majList.append(reference_sequence.objects.get(id=int(seqsInOrderOfAbunIDs[item])).name)
                del seqsInOrderOfAbunIDs[item]
                break
    majStringOutput = '/'.join(majList)
    return majStringOutput, majList

def getAbundStr(totlist, sdlist, majlist):
    majComp = []
    lessComp = []
    for i in range(len(totlist)):
        totString = "{0:.3f}".format(totlist[i])
        sdString = "{0:.3f}".format(sdlist[i])
        if i in range(len(majlist)):
            majComp.append('{}[{}]'.format(totString, sdString))
        else:
            lessComp.append('{}[{}]'.format(totString, sdString))
    majCompStr = '/'.join(majComp)
    lessCompStr = '-'.join(lessComp)
    if lessCompStr:
        abundOutputStr = '-'.join([majCompStr, lessCompStr])
    else:
        abundOutputStr = majCompStr
    return abundOutputStr

def dbProfileOverviewFig():
    cladeList = ['A', 'B', 'C', 'D']

    analysisTypes = analysis_type.objects.filter(dataAnalysisFrom=analysisObj)
    cladalTypeList = [[] for clade in cladeList]
    for ant in analysisTypes:
        clade = ant.clade
        print(clade)
        if clade in cladeList:
            cladalTypeList[cladeList.index(clade)].append(ant)
        else:
            continue
    tw=0
    # Here we have the analysistypes split by clades A B C D
    plt.figure(1)
    for clt in cladalTypeList:
        defDict = defaultdict(int)
        for ant in clt:#for analysistypes of given clade
            # count how many cladeCollections the type was found in
            # i.e. count cladeCollectionTypes of clade
            typeAbundance = len(clade_collection_type.objects.filter(analysisTypeOf = ant))
            defDict[ant] += typeAbundance
        # At this stage we can work out the N most abundant types for the clade
        # We can then work out how many times each of the analysistypes was found per data_set

        # 32: VoolstraSet
        # 40: EmilyData
        # 43: Ryan
        # 47: EdData
        # 48: Roberto
        # 49: Alejandro
        countDict = Counter(defDict)
        # IDs of the respective data_set for use in filtering db filtering below
        voolID = 32
        emID = 40
        ryID = 43


        # TODO put this into for loop format so that we don't to have to write out for each new additin to the DBV
        # we need to end up with a list of lists [[numberofCCsForDS for type in countDict.most_common(20)] for data_set in DBV]
        dataSubmissionsInDBV = list(data_set.objects.filter(id__in=[int(x) for x in analysisObj.listOfDataSubmissions.split(',')]))
        dataForPlottingListOfLists = [[] for dS in dataSubmissionsInDBV]

        typeAbundanceData = [abundance for item, abundance in countDict.most_common(20)]
        typeNames = [item for item, abundance in countDict.most_common(20)]
        for i in range(len(dataForPlottingListOfLists)):
            for item, abundance in countDict.most_common(20):
                CCsWithTypeInDS = list(clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom=dataSubmissionsInDBV[i], clade_collection_type__analysisTypeOf=item).distinct())
                numCCsWithTypeInDS = len(CCsWithTypeInDS)
                dataForPlottingListOfLists[i].append(numCCsWithTypeInDS)
        # Here we have the list populated

        # VVV OLD CODE VVV #
        # for item, abundance in countDict.most_common(20):
        #     typeAbundanceData.append(abundance)
        #     typeNames.append(item)
        #     # Now for each type find out how many cladeCollections they were in for each given data_set
        #     # and append to the relavent list.
        #     #Voolstra first
        #     voolDsSpecificCCsWithType = list(clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom=data_set.objects.get(id=32), clade_collection_type__analysisTypeOf=item).distinct())
        #     numVoolDsSpecificCCsWithType = len(voolDsSpecificCCsWithType)
        #     voolstraAbundanceData.append(numVoolDsSpecificCCsWithType)
        #     # Emily next
        #     emDsSpecificCCsWithType = list(clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom=data_set.objects.get(id=40), clade_collection_type__analysisTypeOf=item).distinct())
        #     numEmDsSpecificCCsWithType = len(emDsSpecificCCsWithType)
        #     emilyAbundanceData.append(numEmDsSpecificCCsWithType)
        #     # Ryan final
        #     ryDsSpecificCCsWithType = list(clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom=data_set.objects.get(id=43), clade_collection_type__analysisTypeOf=item).distinct())
        #     numRyDsSpecificCCsWithType = len(ryDsSpecificCCsWithType)
        #     ryanAbundanceData.append(numRyDsSpecificCCsWithType)
        # ^^^ OLD CODE ^^^ #

        plotnumber = cladalTypeList.index(clt) + 1
        # to adust for the summary plot

        # a list of color hex values that we can cycle through
        colorHexList = os.path.abspath(os.path.dirname(__file__)) + '/newColourList'


        ax = plt.subplot(2, len(cladeList), plotnumber)
        ind = np.arange(len(countDict.most_common(20)))
        width = 0.35

        listOfPs = []
        # First calculate the propdata for each dataSet
        propDataList = []
        for i in range(len(dataForPlottingListOfLists)):
            propDataList.append([a / abund for a, abund in zip(dataForPlottingListOfLists[i], typeAbundanceData)])

        for i in range(len(dataForPlottingListOfLists)):
            propData = propDataList[i]
            newBottom = [0 for anType in typeAbundanceData]
            # This gorgeous little nest calculates the new bottom for each bar
            for j in range(len(typeAbundanceData)): # for each abundance in the newBottomList
                for k in range(0,i):
                    newBottom[j] += propDataList[k][j]
            p = ax.bar(ind, propData, width, bottom=newBottom, color=colorHexList[i])
            listOfPs.append(p)

        # VVV OLD CODE VVV #
        # voolPropData = [a/abund for a, abund in zip(voolstraAbundanceData, typeAbundanceData)]
        # emPropData = [a/abund for a, abund in zip(emilyAbundanceData, typeAbundanceData)]
        # ryPropData = [a/abund for a, abund in zip(ryanAbundanceData, typeAbundanceData)]
        # p0 = ax.bar(ind, voolPropData, width, color='#ffffff')
        # p1 = ax.bar(ind, emPropData, width, bottom=voolPropData, color='#6600cc')
        # newBottom = [L + M for L, M in zip(voolPropData, emPropData)]
        # p2 = ax.bar(ind, ryPropData, width, bottom=newBottom, color='#cc0000')
        # ^^^ OLD CODE ^^^ #

        if tw == 0:
            tw +=1
            # plt.figlegend(handles=(p0, p1, p2), labels=('Voolstra et al.', 'Howells et al.', 'McMinds et al.'), loc='lower left')
            dsNames = [ds.name for ds in dataSubmissionsInDBV]
            plt.figlegend(handles=(listOfPs), labels=(dsNames), loc='lower left')

        ax.get_xaxis().set_visible(False)
        # ax.set_xticks(ind)
        # ax.set_xticklabels(['{}:  {}'.format(typeNm, str(abun)) for typeNm, abun in zip(typeNames, typeAbundanceData)], rotation='vertical')
        ax.set_title('Clade {}'.format(cladeList[cladalTypeList.index(clt)]), fontsize=16)
        if plotnumber == 1:
            ax.set_ylabel('Type-dataset proportions', fontsize=12)



        ax = plt.subplot(2, len(cladeList), plotnumber + len(cladeList))
        p3 = ax.bar(ind, typeAbundanceData, width, color='#000000')
        ax.set_xticks([a + width for a in ind])
        ax.set_xticklabels(['{}:  {}'.format(typeNm, str(abun)) for typeNm, abun in zip(typeNames, typeAbundanceData)], rotation='vertical')
        if plotnumber == 1:
            ax.set_ylabel('Type abundance (num samples)', fontsize=12)
        #ax.set_title('Clade {}'.format(cladeList[cladalTypeList.index(clt)]), fontsize=16)
    plt.suptitle('SymPortal Database Overview', fontsize=18)
    plt.tight_layout()
    # tight_layout doesn't take into account the title and so we have to move the subplots down manually
    plt.subplots_adjust(top=0.92)
    plt.show()
    return

def outputFasta():
    namedRefSeqCol = reference_sequence.objects.filter(hasName=True)
    fasta = []
    for refSeq in namedRefSeqCol:
        fasta.extend(['>{}_{}'.format(refSeq.name, refSeq.id), '{}'.format(refSeq.sequence)])
    outputDir = os.path.join(os.path.dirname(__file__), 'outputs/{}'.format(analysisObj.id))
    os.makedirs(outputDir, exist_ok=True)

    writeListToDestination('{}/SymPortalDBSequences.fasta'.format(outputDir), fasta)
    return

def namingRefSeqsUsedInDefs():
    listOfSeqNamesThatAlreadyExist = [refSeq.name for refSeq in reference_sequence.objects.filter(hasName=True)]
    listOfSeqNamesThatAlreadyExist.append('D1a')

    # The full path of the refSeqDB.fa
    refFasta = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB/referenceSeqDB')) + '/refSeqDB.fa'

    # List of the reference sequence names
    listOfRefSeqDBNames = [a[1:] for a in refFasta if a[0] == '>']

    # Get list of referenceSeqs that are found in the analysis and are used in the type descriptions
    # Get list of analysisTypes
    listOfDefiningRefSeqIDs = set()
    at = analysis_type.objects.filter(dataAnalysisFrom=analysisObj)
    for att in at:
        listOfDefiningRefSeqIDs.update([int(a) for a in att.orderedFootprintList.split(',')])

    # Here we have a list of the refseq IDs that are used in definitions of types and therefore need names
    # Get collection of the actual refseqs
    listOfDefiningRefSeqIDs = list(listOfDefiningRefSeqIDs)

    # we only need to name the refSeqs that don't already have a name
    listOfDefiningRefSeqs = list(reference_sequence.objects.filter(id__in=listOfDefiningRefSeqIDs, hasName=False))

    ### DEBUG ###
    # Perform sanity check to see that all of the seqs are independent of each other e.g. none are subsets or super sets
    for a, b in itertools.combinations(listOfDefiningRefSeqs,2):
        if a.sequence in b.sequence:
            apples = 'asdf'
        elif b.sequence in a.sequence:
            apples = 'asdf'
    ### DEBUG ###

    # When we made the code to assign sequences to refseqs or create new refseqs we were smart in that
    # if we had a new sequence that was bigger than a refseq we collapsed the big seq to the refseq.
    # this means that we should just be able to blast the refseqs against the refSeqDB.fa and call the 100% matches
    # the names in the refSeqDB.
    if listOfDefiningRefSeqs:
        # Perform the blast of refSeq
        # create a fasta file to blast against
        fastaToBlast = []
        for rs in listOfDefiningRefSeqs:
            fastaToBlast.extend(['>{}'.format(rs.id), '{}'.format(rs.sequence)])
        writeDir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB')) + '/unnamedRefSeqs.fasta'
        writeListToDestination(writeDir, fastaToBlast)

        blastOutputPath = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB')) + '/blast.out'

        outputFmt = '6 qseqid sseqid evalue pident qcovs'
        inputPath = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB')) + '/unnamedRefSeqs.fasta'

        os.chdir(os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB')))

        # Run local blast
        # completedProcess = subprocess.run([blastnPath, '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'symbiodinium.fa', '-max_target_seqs', '1', '-num_threads', '1'])
        completedProcess = subprocess.run(
            ['blastn', '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'refSeqDB.fa',
             '-max_target_seqs', '1', '-num_threads', '3'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Read in blast output
        blastOutputFile = readDefinedFileToList(blastOutputPath)

        # First go through and assign names to the 100% coverage and matches
        for bo in blastOutputFile:
            splitEl = bo.split('\t')
            if int(float(splitEl[3])) == 100 and int(splitEl[4]) == 100:
                potentialName = splitEl[1]
                if potentialName not in listOfSeqNamesThatAlreadyExist:
                    # Then we can assign this name to the refSeq
                    refSeqInQ = reference_sequence.objects.get(id=int(splitEl[0]))
                    refSeqInQ.name = potentialName
                    refSeqInQ.hasName = True
                    refSeqInQ.save()
                    listOfSeqNamesThatAlreadyExist.append(potentialName)
                else:
                    apples = 'asdf' # We have an issue if we get here because we have a 100% match to a refseq
                    # but that ref seqname has aleady been associated with a different seq
        # todo we can comment out the below to prevent new names from being created
        # Now assign names to those that aren't exact matches
        for bo in blastOutputFile:
            splitEl=bo.split('\t')
            refSeqInQ = reference_sequence.objects.get(id=int(splitEl[0]))
            if not refSeqInQ.hasName:
                newName = createNewRefSeqName(splitEl[1], listOfSeqNamesThatAlreadyExist, listOfRefSeqDBNames)
                refSeqInQ.name = newName
                refSeqInQ.hasName = True
                refSeqInQ.save()
                listOfSeqNamesThatAlreadyExist.append(newName)
        #Finally update the type names
        IDs = [att.id for att in at]
        for i in range(len(IDs)):
            typeInQ = analysis_type.objects.get(id=IDs[i])
            typeInQ.name = typeInQ.generateName()
            typeInQ.save()

    # update the refSeqDB blast database so that it holds all named sequences
    path_for_refSeqDB = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB', 'refSeqDB.fa'))
    # make a fasta that is all of the named seqs
    fasta_to_write = []
    for rs in reference_sequence.objects.filter(hasName=True):
        fasta_to_write.extend(['>{}'.format(rs.name), '{}'.format(rs.sequence)])

    # here we have the new fasta in a list. Now write it out
    writeListToDestination(path_for_refSeqDB, fasta_to_write)

    # now create blast db from the fasta
    completed_process = subprocess.run(
        ['makeblastdb', '-in', path_for_refSeqDB, '-dbtype', 'nucl', '-title',
         'refSeqDB'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    analysisObj.refSeqsNamed = True
    analysisObj.save()
    return
#####################################



###### MAIN ######
def main(dataanalysistwoobject, cores):

    ##### CLEAN UP tempData FOLDER ####
    # Delete the tempDataFolder and contents
    shutil.rmtree(os.path.abspath(os.path.join(os.path.dirname(__file__), 'temp')))
    # recreate the tempDataFolder
    os.makedirs(os.path.abspath(os.path.join(os.path.dirname(__file__), 'temp')),
                exist_ok=True)

    # Set the dataAnlysisTwoObject in Q as a global object analysisObj for easy reference
    global analysisObj
    analysisObj = dataanalysistwoobject
    global nProcessors
    nProcessors = cores
    #TODO eventually I would like to attach this parameter to the data_set object
    global unlockedAbund
    unlockedAbund = 0.0001




    start_time = timeit.default_timer()
    ####### PROFILE DISCOVERY ##########
    print('Profile discovery')
    if not analysisObj.analysisTypesDefined:
        profileDiscovery(nProcessors)
    ####################################
    elapsed_time = timeit.default_timer() - start_time
    print('\n\nPROFILE DISCOVERY {}'.format(elapsed_time))

    start_time = timeit.default_timer()
    ####### PROFILE ASSIGNMENT #########
    print('\n\nProfile Assignment')
    if not analysisObj.analysisTypesAssigned:
        # This includes the multiModal analyses
        profileAssignment(cores)
    ########################    ############
    elapsed_time = timeit.default_timer() - start_time
    print('PROFILE ASSIGNMENT {}'.format(elapsed_time))


    ####### PROFILE COLLAPSE #########
    start_time = timeit.default_timer()
    print('GROUP ASSIGNMENT')
    if not analysisObj.analysisTypesCollapsed:
        assignProfilesToGroups()
        # profileCollapse()
    elapsed_time = timeit.default_timer() - start_time
    print('PROFILE COLLAPSE {}'.format(elapsed_time))
    ####################################


    ####### SEQUENCE NAMING #########
    print('Naming defining reference Sequences')
    if not analysisObj.refSeqsNamed:
        namingRefSeqsUsedInDefs()
    ####################################


    ####### SPECIES ASSIGNMENT #########
    print('Assigning Species')
    if not analysisObj.speciesAssociated:
        assignSpecies()
    ####################################



    # ####### Outputs #########
    # print('Formatting Output')
    # formatOutput(cores)
    # # dbProfileOverviewFig()
    # outputFasta()
    # ####################################


    ##### CLEAN UP tempData FOLDER ####
    # Delete the tempDataFolder and contents
    shutil.rmtree(os.path.abspath(os.path.join(os.path.dirname(__file__), 'temp')))
    # recreate the tempDataFolder
    os.makedirs(os.path.abspath(os.path.join(os.path.dirname(__file__), 'temp')),
                exist_ok=True)
#################################################


########### Ancillary Function ##################
# A collection of functions that support some of the less major functions of SymPortal e.g. producing figures
# or checking the differentiation of types
##### Type Comparison ######
def compareIntraAbundOfTypes(analysisIdsOrNames, datID):
    '''
    This function will take either a list of IDs of types or a list of type names
    it will produce a text table output showing the relative abundnces of each of the intras
    so that the similarity of the types can be compared
    beacuase different analyses can have types of the same names, the function will also require the dataAnalysisTwoID
    '''

    isNames = False
    try:
        int(analysisIdsOrNames[0])
    except ValueError:
        isNames = True

    if isNames:# Then we are dealing with names rather than id's
        listOfTypes = analysis_type.objects.filter(name__in=analysisIdsOrNames, dataAnalysisFrom= data_analysis.objects.get(id=datID))
    else:
        listOfTypes = analysis_type.objects.filter(id__in=analysisIdsOrNames)



    ############## CREATE INFOLISTS FOR EACH TYPE ################
    # Get a list of all definingRefSeqs in the analysisTypes found in all of the types in the group
    definingRefSeqs = []
    for antype in listOfTypes:
        if str(antype) == '22094' or str(antype) == '22094/22523' or str(antype) == '22557':
            apples = 'asdf'
        definingRefSeqs += antype.getOrderedFootprintList()
    listOfDefiningRefSeqsAllTypes = list(set(definingRefSeqs))
    # Get all analysedSampleSequences that are of these referenceSequences and from the analysisTypes in question
    # https://docs.djangoproject.com/en/1.10/topics/db/queries/
    # 'Lookups that span relationships'
    # Turns out the the DJANGO documentation on this reverse span lookup is awful and wrong.
    # I figured it out using http://stackoverflow.com/questions/3831480/how-to-span-lookup-relationship-in-django
    # However probably best to use my example below as example

    anSmpSeqsOfTypes = data_set_sample_sequence.objects.filter(referenceSequenceOf__in=listOfDefiningRefSeqsAllTypes, cladeCollectionTwoFoundIn__clade_collection_type__analysisTypeOf__in=listOfTypes)


    # for each of the refSeqs add up the abundances of the analysedSampleSequences that are of the refSeq type.
    # Keep these sums in a dict and use them to get the order of most abundance seqs.# This will be the oder that the seqs are compared in. e.g. ordered over all seqs for all types in the group
    refSeqCountDict = {}
    for refSeq in listOfDefiningRefSeqsAllTypes:
        refSeqCountDict[refSeq] = sum(anSmpSeq.abundance for anSmpSeq in anSmpSeqsOfTypes.filter(referenceSequenceOf=refSeq))
    # We should be careful here to check that the most abundant seq found in the orderedListOfDefiningRefSeqs is
    # found in all of the analysis type footprints. We could have a rare case where it is not


    #A dictionary that will hold an info list per type
    infoListDict = {}


    # Some of the managers can be created outside of this below loop
    managerOrderedListOfDefiningRefSeqs = Manager()
    orderedListOfDefiningRefSeqs = managerOrderedListOfDefiningRefSeqs.list([a[0] for a in sorted(refSeqCountDict.items(), key=operator.itemgetter(1), reverse=True)])

    for antype in listOfTypes:
        # Create the infolist

        print('calculating intra info for {}'.format(antype.name))

        #TODO multiprocess here
        # I think that this is the rate limiting step and therefore should be multiprocessed at this level
        # When we have the very slow info lists to make it is for the types which have a lot of clade collections
        # in them.
        # Input will be list of clade collections from the listOfCladeCollectionsOfType
        # Managed objects will be listOfCladeCollectionsOfType, orderedListOfDefiningRefSeqs, infoList
        # This is a lot of objects to be transfering to the processess, hopefully the overhead will not be too large

        # Create Managers (some managers are created outside this loop just above)
        managerListOfCladeCollectionsOfType = Manager()
        listOfCladeCollectionsOfType = managerListOfCladeCollectionsOfType.list(list(clade_collection.objects.filter(clade_collection_type__analysisTypeOf=antype)))

        managerInfoList = Manager()
        infoList = managerInfoList.list([[] for cc in listOfCladeCollectionsOfType])

        # managerNewInfoList = Manager()
        # newInfoList = managerNewInfoList.list([0 for refseq in orderedListOfDefiningRefSeqs])

        CCInputQueue = Queue()

        for CC in listOfCladeCollectionsOfType:
            CCInputQueue.put(CC)

        numProcessors = 20
        for N in range(numProcessors):
            CCInputQueue.put('STOP')

        allProcesses = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        for N in range(numProcessors):
            p = Process(target=workerIntraAbundList, args=(CCInputQueue, orderedListOfDefiningRefSeqs, listOfCladeCollectionsOfType, infoList))
            allProcesses.append(p)
            p.start()

        for p in allProcesses:
            p.join()


        newInfoList = [0 for a in range(len(infoList[0]))]
        for I in range(len(infoList)):
            total = sum(infoList[I])
            for j in range(len(infoList[I])):# For each intra in order
                newInfoList[j] += infoList[I][j] / total


        # create a dict where the anType is the key and the respective infoList is the value
        newInfoList = [x/len(infoList) for x in newInfoList]
        if sum(newInfoList) != 1:
            apples = 'asdf'
            print(str(sum(newInfoList)))

        infoListDict[antype] = newInfoList
        apples = 'pears'
    ########################################################

    # Here we have info lists for each type that is going to be compared
    # For each type the relative abundance of each intra (infoListDict) listed in the order of orderedListOfDefiningRefSeqs

    # Create table object
    x = texttable.Texttable()
    tableList = []

    # One off header row
    header = []
    header.append('Type')
    header.extend(str(refSeq) for refSeq in orderedListOfDefiningRefSeqs)
    tableList.append(header)

    # Add row for each type in comparison
    for antypes in infoListDict.keys():
        rowList = []
        rowList.append(str(antypes))
        rowList.extend([str(a)[:5] for a in infoListDict[antypes]])
        tableList.append(rowList)

    x.add_rows(tableList)

    print(x.draw())

    apples = 'asdf'
    return