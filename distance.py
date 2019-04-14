from dbApp.models import (DataSet, ReferenceSequence, DataSetSampleSequence,
                          AnalysisType, DataSetSample, DataAnalysis, CladeCollection)
import os
import shutil
from multiprocessing import Queue, Process, current_process
import math
import sys
from plumbum import local
import pandas as pd
from django import db
import subprocess
import re
import numpy as np
from skbio.stats.ordination import pcoa
import general
from general import write_list_to_destination, read_defined_file_to_list, convert_interleaved_to_sequencial_fasta_first_line_removal, sqrt_transform_abundance_df
import itertools
from scipy.spatial.distance import braycurtis
from symportal_utils import MothurAnalysis, SequenceCollection
from datetime import datetime

# General methods
class SequenceCollectionComplete(Exception):
    pass


class FseqbootAlignmentGenerator:
    """This class is to be used in the UnifracDistanceCreator to handle the creation of the fseqboot alignments
    for a given clade. It produces a directory (self.output_dir) that contains a number of alignment replicates."""

    def __init__(self, parent_unifrac_dist_creator, clade_in_question, num_reps):
        self.parent_unifrac_dist_creator = parent_unifrac_dist_creator
        self.clade = clade_in_question
        self.num_reps = num_reps
        self.input_fasta_path = self.parent_unifrac_dist_creator.clade_master_fasta_file_aligned_path
        self.output_seqboot_file_path = self.input_fasta_path + '.fseqboot'
        self.fseqboot_local = None
        self._setup_fseqboot_plum_bum_local()
        self.output_dir = os.path.join(self.parent_unifrac_dist_creator.output_dir, self.clade)
        # this is the base of the path that will be used to build the path of each of the alignment replicates
        # to build the complete path a rep_count string (str(count)) will be appended to this base.
        self.fseqboot_clade_root_dir = os.path.join(self.output_dir, 'out_seq_boot_reps')
        os.makedirs(self.fseqboot_clade_root_dir, exist_ok=True)
        self.fseqboot_base = os.path.join(self.fseqboot_clade_root_dir, 'fseqboot_rep_')
        self.fseqboot_file_as_list = None
        self.rep_count = 0
        self.fseqboot_individual_alignment_as_list = None

    def do_fseqboot_alignment_generation(self):
        self._execute_fseqboot()

        self.fseqboot_file_as_list = read_defined_file_to_list(self.output_seqboot_file_path)

        self._divide_fseqboot_file_into_indi_alignments_and_write_out()

    def _divide_fseqboot_file_into_indi_alignments_and_write_out(self):
        # Now divide the fseqboot file up into its 100 different alignments
        # here we use the re to match if its the end of an alignment. More often it will not be and so we will just
        # add the next line to the alignment.
        # when we do find the end of an alignmnet then we write it out. And grab the next line to start the next
        # alignment
        self.fseqboot_individual_alignment_as_list = [self.fseqboot_file_as_list[0]]
        reg_ex = re.compile('[0-9]+$')
        for line in self.fseqboot_file_as_list[1:]:
            reg_ex_matches_list = reg_ex.findall(line)
            if len(reg_ex_matches_list) == 1:
                write_list_to_destination('{}{}'.format(self.fseqboot_base, self.rep_count),
                                          self.fseqboot_individual_alignment_as_list)
                self.fseqboot_individual_alignment_as_list = [line]
                self.rep_count += 1
            else:
                self.fseqboot_individual_alignment_as_list.append(line)
        write_list_to_destination(f'{self.fseqboot_base}{self.rep_count}', self.fseqboot_individual_alignment_as_list)

    def _execute_fseqboot(self):
        sys.stdout.write('\rGenerating multiple fseqboot alignments')
        (self.fseqboot_local[
            '-sequence', self.input_fasta_path, '-outfile', self.output_seqboot_file_path,
            '-test', 'b', '-reps', self.num_reps])()

    def _setup_fseqboot_plum_bum_local(self):
        is_installed = subprocess.call(['which', 'fseqboot'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if is_installed == 0:
            self.fseqboot_local = local["fseqboot"]
        else:
            fseqboot_path = os.path.join(
                self.parent_unifrac_dist_creator.symportal_root_dir, 'lib', 'phylipnew', 'fseqboot')
            if os.path.isfile(fseqboot_path):
                self.fseqboot_local = local[fseqboot_path]
            else:
                raise RuntimeError('Cannot find fseqboot in PATH or in local installation at ./lib/phylipnew/fseqboot\n'
                                   'For instructions on installing the phylipnew dependencies '
                                   'please visit the SymPortal '
                                   'GitHub page: https://github.com/didillysquat/SymPortal_framework/'
                                   'wiki/SymPortal-setup#6-third-party-dependencies')


class UnifracSubcladeHandler:
    """The output of this is a list that contains a list of paths to trees, one for each replicate fasta"""

    def __init__(self, parent_unifrac_creator, clade_in_question, fseqbootbase):
        self.clade = clade_in_question
        self.parent_unifrac_creator = parent_unifrac_creator
        self.fseqbootbase = fseqbootbase
        self.input_queue_of_rep_numbers = Queue()
        self._populate_input_queue()
        self.output_queue_of_paths_to_trees = Queue()
        self.output_dir = os.path.join(self.parent_unifrac_creator.output_dir, self.clade)
        self.list_of_output_tree_paths = []
        self.concat_tree_file_path = os.path.join(self.output_dir, 'concatenated_tree_file')
        self.consensus_tree_file_path = os.path.join(self.output_dir, 'consensus_tree_sumtrees.newick')
        self.unifrac_dist_file_path = None

    def perform_unifrac(self):
        mothur_analysis = MothurAnalysis.init_for_weighted_unifrac(
            tree_path=self.consensus_tree_file_path,
            group_file_path=self.parent_unifrac_creator.clade_master_group_file_path,
            name_file_path=self.parent_unifrac_creator.clade_master_names_file_path, name=self.clade)
        mothur_analysis.execute_weighted_unifrac()
        self.unifrac_dist_file_path = mothur_analysis.dist_file_path

    def _populate_input_queue(self):
        # populate the input queue
        for rep_number in range(self.parent_unifrac_creator.bootstrap_value):
            self.input_queue_of_rep_numbers.put(rep_number)

        # place the stop cues
        for n in range(self.parent_unifrac_creator.num_proc):
            self.input_queue_of_rep_numbers.put('STOP')

    def execute_unifrac_mothur_worker(self):
        """This worker simply runs two mothur analyses. The first makes a .dist phylip distance matrix. The second
        makes a clear cut tree from this distance matrix."""
        all_processes = []

        for n in range(self.parent_unifrac_creator.num_proc):
            p = Process(target=self._mothur_unifrac_pipeline_mp_worker,
                        args=())
            all_processes.append(p)
            p.start()

        # Get list of tree paths from the output queue
        kill_num = 0
        while 1:
            passed_element = self.output_queue_of_paths_to_trees.get()
            if passed_element == 'kill':
                kill_num += 1
                if kill_num == self.parent_unifrac_creator.num_proc:
                    break
            else:
                self.list_of_output_tree_paths.append(passed_element)

        for p in all_processes:
            p.join()

    def _mothur_unifrac_pipeline_mp_worker(self):
        for rep_num in iter(self.input_queue_of_rep_numbers.get, 'STOP'):
            unifrac_mothur_worker = UnifracMothurWorker(rep_num=rep_num, fseqbootbase=self.fseqbootbase)
            unifrac_mothur_worker.make_trees()
            self.output_queue_of_paths_to_trees.put(unifrac_mothur_worker.output_tree_path)
        self.output_queue_of_paths_to_trees.put('kill')

    def create_consensus_tree(self):
        self._concatenate_trees()
        self._create_consensus_tree_with_meta_data_to_be_removed()
        self._rename_tree_nodes_and_remove_meta_data()

    def _rename_tree_nodes_and_remove_meta_data(self):
        """This will use the name file to rename the tree nodes and also remove metadata from the treefile.
        """
        name_file = read_defined_file_to_list(self.parent_unifrac_creator.clade_master_names_file_path)
        name_file_reps = []
        for line in name_file:
            name_file_reps.append(line.split('\t')[0])

        seq_re = re.compile('\d+[_ ]id[\d_]+')
        sys.stdout.write('\rrenaming tree nodes')
        tree_file = read_defined_file_to_list(self.consensus_tree_file_path)

        new_tree_file = []
        for line in tree_file:
            new_str = line
            examine = list(seq_re.findall(line))

            # N.B. the sumtrees.py program was causing some very strange behaviour. It was converting '_' to ' '
            # when they were preceeded by a single digit but leaving them as '_'
            # when there were multiple digits before it
            # it took a long time to find what the problem was. It was causing issues in the mothur UniFrac.
            # You will see that I have modified below by having the space function and replacing ' '
            # with '_' and modifying
            # the regex.
            name_file_rep_match_list = []
            for match_str in examine:
                space = False
                if ' ' in match_str:
                    space = True
                if space:

                    match_str_replced_space = match_str.replace(' ', '_')
                    for name_file_rep in name_file_reps:
                        if name_file_rep.startswith(match_str_replced_space):
                            new_str = re.sub("'{}'".format(match_str), name_file_rep, new_str)
                            name_file_rep_match_list.append(name_file_rep)
                            found = True
                            break
                else:
                    for name_file_rep in name_file_reps:
                        if name_file_rep.startswith(match_str):
                            # then this is a match
                            # now replace the string in the tree file
                            new_str = re.sub('(?<!\d){}'.format(match_str), name_file_rep, new_str)
                            name_file_rep_match_list.append(name_file_rep)
                            break

            # now also remove the metadata which is held between square brackets '[]'
            new_str = re.sub('\[[^\]]*\]', '', new_str)
            new_tree_file.append(new_str)

        # here all of the tree_file names should have been replaced. Now write back out.
        sys.stdout.write('\rwriting out tree')
        write_list_to_destination(self.consensus_tree_file_path, new_tree_file)

    def _create_consensus_tree_with_meta_data_to_be_removed(self):
        """run sumtrees.py on the concatenated tree file to create a
        consensus tree that will need annotating with the distances.
        """
        completed_consensus = subprocess.run([
            'sumtrees.py', '-F', 'newick', '--replace', '-o',
            self.consensus_tree_file_path, self.concat_tree_file_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if completed_consensus.returncode != 0:
            try:
                if 'recursion' in completed_consensus.stdout.decode('ISO-8859-1'):
                    sys.exit('There has been a recursion depth error whilst trying to calculate a consensus tree\n'
                             'This occured whilst calculating between sample distances.\n'
                             'This problem occurs when trees are too big for sumtrees.py to process.\n'
                             'This problem can often be solved by increasing the recursion '
                             'limit used when running sumtrees.py\n'
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
                         'This problem can often be solved by increasing the recursion '
                         'limit used when running sumtrees.py\n'
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

    def _concatenate_trees(self):
        master_tree_file = []
        for i in range(len(self.list_of_output_tree_paths)):
            temp_tree_file = read_defined_file_to_list(self.list_of_output_tree_paths[i])
            for line in temp_tree_file:
                master_tree_file.append(line)
        write_list_to_destination(self.concat_tree_file_path, master_tree_file)


class UnifracMothurWorker:
    def __init__(self, rep_num, fseqbootbase):
        self.rep_num = rep_num
        self.fseqbootbase = fseqbootbase
        self.fasta_as_list = None
        self.sequential_fasta_path = f'{self.fseqbootbase}{self.rep_num}.sequential.fasta'
        self.output_tree_path = None

    def make_trees(self):
        sys.stdout.write(f'\rProcessing rep number {self.rep_num} with {current_process().name}')

        # convert the interleaved fasta to sequential fasta
        self._convert_interleaved_fasta_to_sequential_and_write_out()

        rep_sequence_collection = SequenceCollection(path_to_file=self.sequential_fasta_path, name=self.rep_num)
        mothur_analysis = MothurAnalysis.init_from_sequence_collection(sequence_collection=rep_sequence_collection)
        # produce the distance file that can then be used in the clearcut analysis tree making
        mothur_analysis.execute_dist_seqs()
        mothur_analysis.execute_clearcut()
        self.output_tree_path = mothur_analysis.tree_file_path

    def _convert_interleaved_fasta_to_sequential_and_write_out(self):
        self.fasta_as_list = read_defined_file_to_list('{}{}'.format(self.fseqbootbase, self.rep_num))
        self.fasta_as_list = self._convert_interleaved_to_sequencial_fasta_first_line_removal()
        write_list_to_destination(self.sequential_fasta_path, self.fasta_as_list)

    def _convert_interleaved_to_sequencial_fasta_first_line_removal(self):
        list_seq_names = []
        list_seq_sequences = []
        num_seqs = int(self.fasta_as_list[0].split()[0])
        fasta_cropped = []
        # Get rid of the first line and get rid of the blank lines
        for line in self.fasta_as_list[1:]:
            if line != '':
                fasta_cropped.append(line)

        for i in range(len(fasta_cropped)):
            if i < num_seqs:
                # Then we are on one of the inital lines
                list_seq_names.append(fasta_cropped[i].split()[0])
                list_seq_sequences.append(''.join(fasta_cropped[i].split()[1:]))
            else:
                index = i % num_seqs
                list_seq_sequences[index] += ''.join(fasta_cropped[i].split()[1:])

        out_fasta = []
        for name, seq in zip(list_seq_names, list_seq_sequences):
            out_fasta.extend(['>{}'.format(name), seq])

        return out_fasta


# UniFrac classes
class TypeUnifracSeqAbundanceMPCollection:
    """The purpose of this class is to be used during the BtwnTypeUnifracDistanceCreator as a tidy holder for the
     sequences information that is collected from each AnalysisType that is part of the output."""

    def __init__(self, analysis_type, proc_id, is_sqrt_transf):
        self.fasta_dict = {}
        self.name_dict = {}
        self.group_list = []
        self.ref_seq_id_list = []
        self.analysis_type_or_clade_collection = analysis_type
        self.proc_id = proc_id
        self.is_sqrt_transf = is_sqrt_transf

    def collect_seq_info(self):
        """This function returns information used to build the master fasta and names files that are required to make
        the unifrac. Currently, this information is collected on a type for type basis for every CladeCollection that
        the type was found in regardless of wether the CladeCollections are part of the output or not.
        For example, if a type was found in samples 1, 3, 4, and 7. and 1, 3, 4 are part of this output, the
        average relative abundances of the type in question will still be calculated using the relative abundance
        information from all four of the CladeCoolection.
        TODO However, in some circumstances it will be useful to be
        able to limit which instances of the type are used in defining the average relative abundances of the DIVs.
        To do this we will provide a flag that is --local. In the above example, if --local was applied then the local
        abundances would be calculated only from the datasetsamples in question."""
        sys.stdout.write('\rProcessing AnalysisType: {} with {}'.format(self.analysis_type_or_clade_collection, self.proc_id))
        ref_seq_uids_of_analysis_type = [int(b) for b in self.analysis_type_or_clade_collection.ordered_footprint_list.split(',')]

        df = pd.DataFrame(self.analysis_type_or_clade_collection.get_ratio_list())

        if self.is_sqrt_transf:
            df = general.sqrt_transform_abundance_df(df)

        normalised_abundance_of_divs_dict = {
            ref_seq_uids_of_analysis_type[i]: math.ceil(df[i].mean() * 10000) for
            i in range(len(ref_seq_uids_of_analysis_type))}

        for ref_seq in ReferenceSequence.objects.filter(id__in=ref_seq_uids_of_analysis_type):
            self.ref_seq_id_list.append(ref_seq.id)
            unique_seq_name_base = '{}_id{}'.format(ref_seq.id, self.analysis_type_or_clade_collection.id)

            analysis_uid_str = str(self.analysis_type_or_clade_collection.id)

            self.fasta_dict[
                '{}_{}'.format(unique_seq_name_base, 0)] = ref_seq.sequence
            temp_name_list = []

            for i in range(normalised_abundance_of_divs_dict[ref_seq.id]):
                temp_name_list.append('{}_{}'.format(unique_seq_name_base, i))
                self.group_list.append('{}\t{}'.format('{}_{}'.format(unique_seq_name_base, i), analysis_uid_str))

            self.name_dict['{}_{}'.format(unique_seq_name_base, 0)] = temp_name_list


class SampleUnifracSeqAbundanceMPCollection:
    """The purpose of this class is to be used during the BtwnSampleUnifracDistanceCreator as a tidy holder for the
     sequences information that is collected from each CladeCollection that is part of the output."""

    def __init__(self, clade_collection, proc_id, is_sqrt_transf):
        self.fasta_dict = {}
        self.name_dict = {}
        self.group_list = []
        self.ref_seq_id_list = []
        self.analysis_type_or_clade_collection = clade_collection
        self.proc_id = proc_id
        self.is_sqrt_transf = is_sqrt_transf

    def collect_seq_info(self):
        sys.stdout.write('\rProcessing cc: {} with {}'.format(self.analysis_type_or_clade_collection, self.proc_id))
        list_of_dsss_in_cc = list(DataSetSampleSequence.objects.filter(
                clade_collection_found_in=self.analysis_type_or_clade_collection))
        normalisation_sequencing_depth = 10000
        if self.is_sqrt_transf:
            normalised_abund_dict = self._make_norm_abund_dict_sqrt(list_of_dsss_in_cc, normalisation_sequencing_depth)
        else:
            normalised_abund_dict = self._make_norm_abund_dict_no_sqrt(list_of_dsss_in_cc,
                                                                       normalisation_sequencing_depth)

        for data_set_sample_seq in list_of_dsss_in_cc:
            ref_seq_id = data_set_sample_seq.reference_sequence_of.id
            self.ref_seq_id_list.append(ref_seq_id)
            unique_seq_name_base = '{}_id{}'.format(ref_seq_id, self.analysis_type_or_clade_collection.id)

            smp_name = str(self.analysis_type_or_clade_collection.id)
            self.fasta_dict[
                '{}_{}'.format(unique_seq_name_base, 0)] = data_set_sample_seq.reference_sequence_of.sequence
            temp_name_list = []

            for i in range(normalised_abund_dict[data_set_sample_seq.id]):
                temp_name_list.append('{}_{}'.format(unique_seq_name_base, i))
                self.group_list.append('{}\t{}'.format('{}_{}'.format(unique_seq_name_base, i), smp_name))

            self.name_dict['{}_{}'.format(unique_seq_name_base, 0)] = temp_name_list

    def _make_norm_abund_dict_no_sqrt(self, list_of_dsss_in_cc, normalisation_sequencing_depth):
        total_seqs_of_cc = sum([dss.abundance for dss in list_of_dsss_in_cc])
        normalised_abund_dict = {dsss.id: int((dsss.abundance / total_seqs_of_cc) * normalisation_sequencing_depth) for
                                 dsss in list_of_dsss_in_cc}
        return normalised_abund_dict

    def _make_norm_abund_dict_sqrt(self, list_of_dsss_in_cc, normalisation_sequencing_depth):
        sqrt_abundances_dict = {
            dsss.id: math.sqrt(dsss.abundance) for dsss in
            list_of_dsss_in_cc}
        total_seqs_of_cc = sum(sqrt_abundances_dict.values())
        normalised_abund_dict = {
        dsss.id: int((sqrt_abundances_dict[dsss.id] / total_seqs_of_cc) * normalisation_sequencing_depth)
        for dsss in list_of_dsss_in_cc}
        return normalised_abund_dict


class BaseUnifracDistanceCreatorHandlerOne:
    """This is the set of classes (one base, two derived (one for samples one for profiles)) that generate
    the master fasta, names and group files per clade that are used in calculating the UniFrac distances.
    """
    def __init__(self, parent_unifrac_dist_creator, clade_in_question):
        self.parent_unifrac_dist_creator = parent_unifrac_dist_creator
        self.clade = clade_in_question
        self.output_dir = os.path.join(self.parent_unifrac_dist_creator.output_dir, self.clade)
        os.makedirs(self.output_dir, exist_ok=True)
        self.output_unifrac_seq_abund_mp_collection_queue = Queue()
        # variables for processing and consolidating the sequence collection outputs
        self.master_fasta_dict = {}
        self.master_unique_fasta_seq_list = []
        self.master_name_dict = {}
        self.master_group_file_as_list = []
        # this dict will link the ref_seq ID to the name of the instance of the ref_seq_id that was used
        self.master_name_unique_id_dict = {}
        # we will have the worker of each process put a 'EXIT' into its output que when it has finished
        # this way we can count how many 'EXIT's we have had output and when this equals the number of processors
        # we started then we know that they have all finished and we have processed all of the outputs
        self.stop_count = 0
        self.master_names_file_as_list = []
        self.master_fasta_file_as_list = []
        self.master_names_file_path = os.path.join(self.output_dir, 'unifrac_master_names_file.names')
        self.master_fasta_file_path = os.path.join(self.output_dir, 'unifrac_master_fasta_file.fasta')
        self.master_group_file_path = os.path.join(self.output_dir, 'unifrac_master_group_file.fasta')

    def _populate_the_master_seq_collection_objects(self):
        """This method collects the output of the sequence collection mp methods and populates the objects that
        represnt the master names, fasta and group files that we are creating to calculate the UniFrac"""
        for seq_collection_object in iter(self.output_unifrac_seq_abund_mp_collection_queue.get, 'STOP'):
            try:
                should_continue = self._check_if_all_seq_collection_objects_have_been_processed(seq_collection_object)
                if should_continue:
                    continue
                sys.stdout.write(
                    f'\rAdding {str(seq_collection_object.analysis_type_or_clade_collection)} to master fasta and name files')

                self._update_master_group_list(seq_collection_object)

                for seq_id in seq_collection_object.ref_seq_id_list:
                    seq_name = f'{seq_id}_id{seq_collection_object.analysis_type_or_clade_collection.id}_0'
                    if seq_id not in self.master_unique_fasta_seq_list:
                        self._populate_new_seq_info_in_master_dict_objects(seq_collection_object, seq_id, seq_name)
                    else:
                        self._populate_existing_seq_info_in_master_dict_objects(seq_collection_object, seq_id, seq_name)
            except SequenceCollectionComplete:
                break

    def _populate_existing_seq_info_in_master_dict_objects(self, seq_collection_object, seq_id, seq_name):
        # then this ref seq is already in the fasta so we just need to update the master name dict
        self.master_name_dict[self.master_name_unique_id_dict[seq_id]] += seq_collection_object.name_dict[seq_name]

    def _populate_new_seq_info_in_master_dict_objects(self, seq_collection_object, seq_id, seq_name):
        self.master_unique_fasta_seq_list.append(seq_id)
        # then this ref seq has not yet been added to the fasta and it needs to be
        # populate the master fasta
        self.master_fasta_dict[seq_name] = seq_collection_object.fasta_dict[seq_name]
        # update the master_name_unique_id_dict to keep track of
        # which sequence name represents the ref_seq_id
        self.master_name_unique_id_dict[seq_id] = seq_name
        # create a name dict entry
        self.master_name_dict[seq_name] = seq_collection_object.name_dict[seq_name]

    def _update_master_group_list(self, seq_collection_object):
        self.master_group_file_as_list += seq_collection_object.group_list

    def _check_if_all_seq_collection_objects_have_been_processed(self, seq_collection_object):
        if seq_collection_object == 'EXIT':
            self.stop_count += 1
            if self.stop_count == self.parent_unifrac_dist_creator.num_proc:
                raise SequenceCollectionComplete
            return True
        return False

    def _convert_fasta_and_name_dict_to_lists_for_writing(self):
        # convert names_dict to names file as list
        self._name_dict_to_file_as_list()

        # convert fasta dict to list
        self._fasta_dict_to_file_as_list()

    def _name_dict_to_file_as_list(self):
        for k, v in self.master_name_dict.items():
            name_list_string = ','.join(v)
            self.master_names_file_as_list.append(f'{k}\t{name_list_string}')

    def _fasta_dict_to_file_as_list(self):
        for k, v in self.master_fasta_dict.items():
            self.master_fasta_file_as_list.extend([f'>{k}', v])

    def _write_out_master_fasta_names_and_group_files(self):
        write_list_to_destination(self.master_fasta_file_path, self.master_fasta_file_as_list)
        write_list_to_destination(self.master_names_file_path, self.master_names_file_as_list)
        write_list_to_destination(self.master_group_file_path, self.master_group_file_as_list)

    @staticmethod
    def _once_complete_wait_for_processes_to_complete(all_processes):
        # process the outputs of the sub processess before we pause to wait for them to complete.
        for p in all_processes:
            p.join()


class BtwnTypeUnifracDistanceCreatorHandlerOne(BaseUnifracDistanceCreatorHandlerOne):
    """Derived class specifically for the between ITS2 type profile distance calculations"""

    def __init__(self, parent_unifrac_dist_creator, clade_in_question):
        super().__init__(
            parent_unifrac_dist_creator=parent_unifrac_dist_creator,
            clade_in_question=clade_in_question)

        self.analysis_types_of_clade = self.parent_unifrac_dist_creator.analysis_types_from_data_set_samples.filter(
            clade=self.clade)
        self._raise_runtime_error_if_not_enough_analysis_types()
        self.input_analysis_type_queue = Queue()
        self.output_unifrac_seq_abund_mp_collection_queue = Queue()
        self._populate_input_queue()


    def execute_unifrac_distance_creator_worker_one(self):
        all_processes = self._start_sequence_collection_running()

        self._populate_the_master_seq_collection_objects()

        self._once_complete_wait_for_processes_to_complete(all_processes)

        self._convert_fasta_and_name_dict_to_lists_for_writing()

        self._write_out_master_fasta_names_and_group_files()

    def _start_sequence_collection_running(self):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        for n in range(self.parent_unifrac_dist_creator.num_proc):
            p = Process(target=self._unifrac_distance_creator_worker_one, args=())
            all_processes.append(p)
            p.start()
        return all_processes

    def _unifrac_distance_creator_worker_one(self):
        proc_id = current_process().name
        for at in iter(self.input_analysis_type_queue.get, 'STOP'):
            unifrac_seq_abundance_mp_collection = TypeUnifracSeqAbundanceMPCollection(analysis_type=at, proc_id=proc_id, is_sqrt_transf=self.parent_unifrac_dist_creator.is_sqrt_transf)
            unifrac_seq_abundance_mp_collection.collect_seq_info()
            self.output_unifrac_seq_abund_mp_collection_queue.put(unifrac_seq_abundance_mp_collection)
        self.output_unifrac_seq_abund_mp_collection_queue.put('EXIT')

    def _populate_input_queue(self):
        for at in self.analysis_types_of_clade:
            self.input_analysis_type_queue.put(at)

        for n in range(self.parent_unifrac_dist_creator.num_proc):
            self.input_analysis_type_queue.put('STOP')

    def _raise_runtime_error_if_not_enough_analysis_types(self):
        if len(self.analysis_types_of_clade) < 2:
            raise RuntimeWarning({'message': 'insufficient objects of clade'})


class BtwnSampleUnifracDistanceCreatorHandlerOne(BaseUnifracDistanceCreatorHandlerOne):
    """Derived class specifically for the between ample distance calculations"""

    def __init__(self, parent_unifrac_dist_creator, clade_in_question):
        super().__init__(
            parent_unifrac_dist_creator=parent_unifrac_dist_creator,
            clade_in_question=clade_in_question)
        self.clade_collections_of_clade = self.parent_unifrac_dist_creator.clade_collections_from_data_set_samples.filter(
            clade=self.clade)
        self._raise_runtime_error_if_not_enough_clade_collections()
        self.input_clade_collection_queue = Queue()

        self._populate_input_queue()

    def execute_unifrac_distance_creator_worker_one(self):
        all_processes = self._start_sequence_collection_running()

        self._populate_the_master_seq_collection_objects()

        self._once_complete_wait_for_processes_to_complete(all_processes)

        self._convert_fasta_and_name_dict_to_lists_for_writing()

        self._write_out_master_fasta_names_and_group_files()

    def _unifrac_distance_creator_worker_one(self):
        proc_id = current_process().name
        for cc in iter(self.input_clade_collection_queue.get, 'STOP'):
            unifrac_seq_abundance_mp_collection = SampleUnifracSeqAbundanceMPCollection(
                clade_collection=cc, proc_id=proc_id, is_sqrt_transf=self.parent_unifrac_dist_creator.is_sqrt_transf)
            unifrac_seq_abundance_mp_collection.collect_seq_info()
            self.output_unifrac_seq_abund_mp_collection_queue.put(unifrac_seq_abundance_mp_collection)
        self.output_unifrac_seq_abund_mp_collection_queue.put('EXIT')

    def _start_sequence_collection_running(self):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        for n in range(self.parent_unifrac_dist_creator.num_proc):
            p = Process(target=self._unifrac_distance_creator_worker_one, args=())
            all_processes.append(p)
            p.start()
        return all_processes

    def _populate_input_queue(self):
        for cc in self.clade_collections_of_clade:
            self.input_clade_collection_queue.put(cc)

        for n in range(self.parent_unifrac_dist_creator.num_proc):
            self.input_clade_collection_queue.put('STOP')

    def _raise_runtime_error_if_not_enough_clade_collections(self):
        if len(self.clade_collections_of_clade) < 2:
            raise RuntimeWarning({'message': 'insufficient objects of clade'})


class BaseUnifracDistPCoACreator:
    """Base class for BtwnTypeUnifracDistPCoACreator and BtwnSampleUnifracDistPCoACreator.
    These classes are used for generating UniFrac distances between either ITS2 type profiles
    or Samples."""
    def __init__(
            self, num_proc, bootstrap_val, output_dir, data_set_uid_list, data_set_sample_uid_list,
            symportal_root_directory, call_type, date_time_string, profiles_or_samples, is_sqrt_transf):
        # 'profiles' or 'samples'
        self.profiles_or_samples = profiles_or_samples
        self.num_proc = num_proc
        self.bootstrap_value = bootstrap_val
        self.output_dir = output_dir
        self.data_set_sample_uid_list = self._set_data_set_sample_uid_list(data_set_sample_uid_list, data_set_uid_list)
        self.output_file_paths = []

        # Clade based files
        # the paths output from the creation of the master fasta, names and group files
        # these will be updated with every clade_in_question
        self.clade_master_names_file_path = None
        self.clade_master_fasta_file_unaligned_path = None
        self.clade_master_fasta_file_aligned_path = None
        self.clade_master_group_file_path = None
        # this is the folder in which all of the replicate work was done. We will use this to delete this dir later
        self.clade_fseqboot_root_dir = None
        # this is the base of the path that will be used to build the path of each of the alignment replicates
        # to build the complete path a rep_count string (str(count)) will be appended to this base.
        self.clade_fseqboot_base = None
        # the output list of trees from the mothur portion of this analysis
        self.clade_tree_path_list = None

        self.symportal_root_dir = symportal_root_directory
        self.call_type = call_type
        if date_time_string:
            self.date_time_string = date_time_string
        else:
            self.date_time_string = str(datetime.now()).replace(' ', '_').replace(':', '-')
        self.output_file_paths = []
        self.clade_output_dir = None
        # path to the .csv file that will hold the PCoA coordinates
        self.clade_pcoa_coord_file_path = None
        # path to the .dist file that holds the unifrac or braycurtis derived sample clade-separated paired distances
        self.clade_dist_file_path = None
        self.clade_dist_file_as_list = None
        self.is_sqrt_transf = is_sqrt_transf

    def _set_data_set_sample_uid_list(self, data_set_sample_uid_list, data_set_uid_list):
        if data_set_sample_uid_list:
            return data_set_sample_uid_list
        else:
            return [dss.id for dss in DataSetSample.objects.filter(data_submission_from__in=data_set_uid_list)]

    def _write_output_paths_to_stdout(self):
        print('UniFrac and PCoA computation complete. Ouput files:')
        for output_path in self.output_file_paths:
            print(output_path)

    def _append_output_files_to_output_list(self):
        self.output_file_paths.extend([self.clade_dist_file_path, self.clade_pcoa_coord_file_path])

    def _clean_up_temp_files(self):
        if os.path.exists(self.clade_fseqboot_root_dir):
            shutil.rmtree(path=self.clade_fseqboot_root_dir)

        # now delte all files except for the .csv that holds the coords and the .dist that holds the dists
        list_of_dir = os.listdir(self.clade_output_dir)
        for item in list_of_dir:
            if '.csv' not in item and '.dist' not in item:
                os.remove(os.path.join(self.clade_output_dir, item))

    def _compute_unifrac_distances(self, clade_in_question):
        unifrac_subclade_handler = UnifracSubcladeHandler(
            clade_in_question=clade_in_question,
            fseqbootbase=self.clade_fseqboot_base,
            parent_unifrac_creator=self)
        unifrac_subclade_handler.execute_unifrac_mothur_worker()
        unifrac_subclade_handler.create_consensus_tree()
        unifrac_subclade_handler.perform_unifrac()
        self.clade_dist_file_path = unifrac_subclade_handler.unifrac_dist_file_path
        self._append_date_time_string_to_unifrac_dist_path(clade_in_question)
        self.clade_dist_file_as_list = [line.replace(' ', '') for line in
                                        read_defined_file_to_list(self.clade_dist_file_path)]

    def _make_fseqboot_replicate_alignments(self, clade_in_question):
        fseqboot_alignment_generator = FseqbootAlignmentGenerator(clade_in_question=clade_in_question,
                                                                  parent_unifrac_dist_creator=self,
                                                                  num_reps=self.bootstrap_value)
        fseqboot_alignment_generator.do_fseqboot_alignment_generation()
        self.clade_fseqboot_base = fseqboot_alignment_generator.fseqboot_base
        self.clade_fseqboot_root_dir = fseqboot_alignment_generator.fseqboot_clade_root_dir

    def _append_date_time_string_to_unifrac_dist_path(self, clade_in_question):
        # here add a date_time_string element to it to make it unique
        old_dist_path = self.clade_dist_file_path
        directory_of_path = os.path.dirname(old_dist_path)
        new_file_name = f'{self.date_time_string}_unifrac_{self.profiles_or_samples}_distances_{clade_in_question}.dist'
        self.clade_dist_file_path = os.path.join(directory_of_path, new_file_name)
        os.rename(old_dist_path, self.clade_dist_file_path)

    def _align_master_fasta(self):
        self.clade_master_fasta_file_aligned_path = self.clade_master_fasta_file_unaligned_path.replace(
            '.fasta', '.aligned.fasta')
        general.mafft_align_fasta(
            input_path=self.clade_master_fasta_file_unaligned_path,
            output_path=self.clade_master_fasta_file_aligned_path,
            num_proc=self.num_proc)

    def _compute_pcoa_coords(self, clade, profiles_or_samples_string):
        # simultaneously grab the sample names in the order of the distance matrix and put the matrix into
        # a twoD list and then convert to a numpy array
        self.clade_pcoa_coord_file_path = os.path.join(
            self.clade_output_dir, f'{self.date_time_string}.unifrac_{profiles_or_samples_string}_PCoA_coords_{clade}.csv')
        raw_dist_file = read_defined_file_to_list(self.clade_dist_file_path)

        temp_two_d_list = []
        sample_names_from_dist_matrix = []
        sample_ids_from_dist_matrix = []
        for line in raw_dist_file[1:]:
            temp_elements = line.split('\t')
            sample_names_from_dist_matrix.append(temp_elements[0].replace(' ', ''))
            sample_ids_from_dist_matrix.append(int(temp_elements[1]))
            temp_two_d_list.append([float(a) for a in temp_elements[2:]])

        dist_as_np_array = np.array(temp_two_d_list)

        sys.stdout.write('\rcalculating PCoA coordinates')

        pcoa_output = pcoa(dist_as_np_array)

        # rename the pcoa dataframe index as the sample names
        pcoa_output.samples['sample'] = sample_names_from_dist_matrix
        renamed_pcoa_dataframe = pcoa_output.samples.set_index('sample')

        # now add the variance explained as a final row to the renamed_dataframe
        renamed_pcoa_dataframe = renamed_pcoa_dataframe.append(
            pcoa_output.proportion_explained.rename('proportion_explained'))

        sample_ids_from_dist_matrix.append(0)
        renamed_pcoa_dataframe.insert(loc=0, column='sample_uid', value=sample_ids_from_dist_matrix)
        renamed_pcoa_dataframe['sample_uid'] = renamed_pcoa_dataframe['sample_uid'].astype('int')

        renamed_pcoa_dataframe.to_csv(self.clade_pcoa_coord_file_path, index=True, header=True, sep=',')


class TypeUnifracDistPCoACreator(BaseUnifracDistPCoACreator):
    """Class for calculating the UniFrac distances between ITS2 type profiles, producing PCoA coordinates from
    these distances and plotting the resultant PCoA coordinates. These are all done on a clade by clade basis.
    TODO implement the --local flag. Also implement a --cct_uids flag which will allow users calculate distances
     for specific sets of its2 type profiles where the distances are calculated using specific sets of CladeCollections
     to work out the average relative abundances of the DIVs.
    """
    def __init__(
            self, symportal_root_directory, num_processors, call_type, data_analysis_obj, date_time_string=None,
            bootstrap_value=100, output_dir=None, data_set_uid_list=None, data_set_sample_uid_list=None, is_sqrt_transf=False):

        super().__init__(
            num_proc=num_processors,bootstrap_val=bootstrap_value, output_dir=output_dir,
            data_set_uid_list=data_set_uid_list, data_set_sample_uid_list=data_set_sample_uid_list, symportal_root_directory=symportal_root_directory, call_type=call_type,
            date_time_string=date_time_string, profiles_or_samples='profiles', is_sqrt_transf=is_sqrt_transf)

        self.data_analysis_obj = data_analysis_obj
        self.analysis_types_from_data_set_samples = AnalysisType.objects.filter(
            data_analysis_from=self.data_analysis_obj,
            cladecollectiontype__clade_collection_found_in__data_set_sample_from__in=self.data_set_sample_uid_list)
        self.clades_for_dist_calcs = list(set([at.clade for at in self.analysis_types_from_data_set_samples]))
        self.output_dir = self._setup_output_dir(call_type, output_dir)

    def compute_unifrac_dists_and_pcoa_coords(self):
        for clade_in_question in self.clades_for_dist_calcs:
            self.clade_output_dir = os.path.join(self.output_dir, clade_in_question)
            try:
                self._create_and_write_out_master_fasta_names_and_group_files(clade_in_question)
            except RuntimeWarning as w:
                if str(w) == 'insufficient objects of clade':
                    continue

            self._align_master_fasta()

            # The production of an ML tree is going to be unrealistic with the larger number of sequences
            # it will simply take too long to compute.
            # Instead I will create an NJ consensus tree.
            # This will involve using the embassy version of the phylip executables

            # First I will have to create random data sets
            self._make_fseqboot_replicate_alignments(clade_in_question)

            self._compute_unifrac_distances(clade_in_question)

            self._add_sample_uids_to_dist_file_and_write()

            self._compute_pcoa_coords(clade=clade_in_question, profiles_or_samples_string='profiles')

            self._clean_up_temp_files()

            self._append_output_files_to_output_list()

        self._write_output_paths_to_stdout()

    def _add_sample_uids_to_dist_file_and_write(self):
        # for the output version lets also append the sample name to each line so that we can see which sample it is
        # it is important that we otherwise work eith the sample ID as the sample names may not be unique.
        dist_with_sample_name = [self.clade_dist_file_as_list[0]]
        list_of_at_ids = [int(line.split('\t')[0]) for line in self.clade_dist_file_as_list[1:]]
        at_of_outputs = list(AnalysisType.objects.filter(id__in=list_of_at_ids))
        dict_of_at_id_to_analysis_type_name = {at.id: at.name for at in at_of_outputs}
        for line in self.clade_dist_file_as_list[1:]:
            temp_list = []
            at_id = int(line.split('\t')[0])
            sample_name = dict_of_at_id_to_analysis_type_name[at_id]
            temp_list.append(sample_name)
            temp_list.extend(line.split('\t'))
            new_line = '\t'.join(temp_list)
            dist_with_sample_name.append(new_line)
        self.clade_dist_file_as_list = dist_with_sample_name
        write_list_to_destination(self.clade_dist_file_path, self.clade_dist_file_as_list)

    def _create_and_write_out_master_fasta_names_and_group_files(self, clade_in_question):
        sys.stdout.write('Creating master .name and .fasta files for UniFrac')
        try:
            unifrac_dist_creator_handler_one = BtwnTypeUnifracDistanceCreatorHandlerOne(
                clade_in_question=clade_in_question,
                parent_unifrac_dist_creator=self)
        except RuntimeWarning as w:
            if w.args[0]['message'] == 'insufficient objects of clade':
                raise RuntimeWarning('insufficient objects of clade')
            else:
                raise RuntimeError(f'Unknown error in {BtwnSampleUnifracDistanceCreatorHandlerOne.__name__} init')
        unifrac_dist_creator_handler_one.execute_unifrac_distance_creator_worker_one()
        self.clade_master_names_file_path = unifrac_dist_creator_handler_one.master_names_file_path
        self.clade_master_fasta_file_unaligned_path = unifrac_dist_creator_handler_one.master_fasta_file_path
        self.clade_master_group_file_path = unifrac_dist_creator_handler_one.master_group_file_path

    def _setup_output_dir(self, call_type, output_dir):
        if call_type == 'stand_alone':
            return os.path.abspath(
                os.path.join(
                    self.symportal_root_dir, 'outputs', 'ordination', self.date_time_string, 'between_profiles'))
        else:
            # call_type == 'analysis':
            return os.path.join(output_dir, 'between_profile_distances')


class SampleUnifracDistPCoACreator(BaseUnifracDistPCoACreator):
    """

    This method will generate a distance matrix between samples using the UniFrac method.
    One for each clade.
    It will also perform a PCoA for each distance matrix. a .dist file and a .csv with the pcoa coords will be output

    The call_type argument will be used to determine which setting this method is being called from.
    if it is being called as part of the initial submission call_type='submission',
    then we will always be working with a single
    data_set. In this case we should output to the same folder that the submission results were output
    to. In the case of being output in a standalone manner (call_type='stand_alone') then we may be outputting
    comparisons from several data_sets. As such we cannot rely on being able to put the ordination results
    into the initial submissions folder. In this case we will use the directory structure that is already
    in place which will put it in the ordination folder.

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
    """

    def __init__(
            self, symportal_root_directory, num_processors, call_type, date_time_string=None,
            bootstrap_value=100, output_dir=None, data_set_uid_list=None, data_set_sample_uid_list=None, is_sqrt_transf=False):
        super().__init__(
            num_proc=num_processors, bootstrap_val=bootstrap_value, output_dir=output_dir,
            data_set_uid_list=data_set_uid_list, data_set_sample_uid_list=data_set_sample_uid_list,
            symportal_root_directory=symportal_root_directory, call_type=call_type,
            date_time_string=date_time_string, profiles_or_samples='samples', is_sqrt_transf=is_sqrt_transf)

        self.clade_collections_from_data_set_samples = CladeCollection.objects.filter(
            data_set_sample_from__in=self.data_set_sample_uid_list)
        self.clades_for_dist_calcs = list(set([a.clade for a in self.clade_collections_from_data_set_samples]))
        self.output_dir = self._setup_output_dir(call_type, output_dir)


    def compute_unifrac_dists_and_pcoa_coords(self):
        for clade_in_question in self.clades_for_dist_calcs:

            self.clade_output_dir = os.path.join(self.output_dir, clade_in_question)

            try:
                self._create_and_write_out_master_fasta_names_and_group_files(clade_in_question)
            except RuntimeWarning as w:
                if str(w) == 'insufficient objects of clade':
                    continue

            self._align_master_fasta()

            # The production of an ML tree is going to be unrealistic with the larger number of sequences
            # it will simply take too long to compute.
            # Instead I will create an NJ consensus tree.
            # This will involve using the embassy version of the phylip executables

            # First I will have to create random data sets
            self._make_fseqboot_replicate_alignments(clade_in_question)

            self._compute_unifrac_distances(clade_in_question)

            self._add_sample_uids_to_dist_file_and_write()

            self._compute_pcoa_coords(clade=clade_in_question, profiles_or_samples_string='samples')

            self._clean_up_temp_files()

            self._append_output_files_to_output_list()

        self._write_output_paths_to_stdout()

    def _add_sample_uids_to_dist_file_and_write(self):
        # for the output version lets also append the sample name to each line so that we can see which sample it is
        # it is important that we otherwise work eith the sample ID as the sample names may not be unique.
        dist_with_sample_name = [self.clade_dist_file_as_list[0]]
        list_of_cc_ids = [int(line.split('\t')[0]) for line in self.clade_dist_file_as_list[1:]]
        cc_of_outputs = list(CladeCollection.objects.filter(id__in=list_of_cc_ids))
        dict_of_cc_id_to_sample_name = {cc.id: cc.data_set_sample_from.name for cc in cc_of_outputs}
        for line in self.clade_dist_file_as_list[1:]:
            temp_list = []
            cc_id = int(line.split('\t')[0])
            sample_name = dict_of_cc_id_to_sample_name[cc_id]
            temp_list.append(sample_name)
            temp_list.extend(line.split('\t'))
            new_line = '\t'.join(temp_list)
            dist_with_sample_name.append(new_line)
        self.clade_dist_file_as_list = dist_with_sample_name
        write_list_to_destination(self.clade_dist_file_path, self.clade_dist_file_as_list)

    def _create_and_write_out_master_fasta_names_and_group_files(self, clade_in_question):
        sys.stdout.write('Creating master .name and .fasta files for UniFrac')
        try:
            unifrac_dist_creator_handler_one = BtwnSampleUnifracDistanceCreatorHandlerOne(
                clade_in_question=clade_in_question,
                parent_unifrac_dist_creator=self)
        except RuntimeWarning as w:
            if w.args[0]['message'] == 'insufficient objects of clade':
                raise RuntimeWarning('insufficient objects of clade')
            else:
                raise RuntimeError(f'Unknown error in {BtwnSampleUnifracDistanceCreatorHandlerOne.__name__} init')
        unifrac_dist_creator_handler_one.execute_unifrac_distance_creator_worker_one()
        self.clade_master_names_file_path = unifrac_dist_creator_handler_one.master_names_file_path
        self.clade_master_fasta_file_unaligned_path = unifrac_dist_creator_handler_one.master_fasta_file_path
        self.clade_master_group_file_path = unifrac_dist_creator_handler_one.master_group_file_path

    def _setup_output_dir(self, call_type, output_dir):
        if call_type == 'stand_alone':
            return os.path.abspath(
                os.path.join(
                    self.symportal_root_dir, 'outputs', 'ordination', self.date_time_string.replace('.','_'), 'between_samples'))
        else:
            # call_type == 'submission' or 'analysis':
            return os.path.join(output_dir, 'between_sample_distances')


# BrayCurtis classes
class BaseBrayCurtisDistPCoACreator:
    def __init__(self, symportal_root_directory, call_type, date_time_string, profiles_or_samples):
        self.symportal_root_dir = symportal_root_directory
        self.call_type = call_type
        if date_time_string:
            self.date_time_string = date_time_string
        else:
            self.date_time_string = str(datetime.now()).replace(' ', '_').replace(':', '-')
        self.output_file_paths = []
        self.clade_output_dir = None
        # path to the .csv file that will hold the PCoA coordinates
        self.clade_pcoa_coord_file_path = None
        # path to the .dist file that holds the unifrac or braycurtis derived sample clade-separated paired distances
        self.clade_dist_file_path = None
        self.clade_dist_file_as_list = None
        # either 'profiles' or 'samples'
        self.profiles_or_samples = profiles_or_samples

        # clade specific attributes. Will be updated for every clade processe
        self.objs_of_clade = None
        self.clade_rs_uid_to_normalised_abund_clade_dict = {}
        self.clade_within_clade_distances_dict = {}

    def _compute_pcoa_coords(self, clade):
        # simultaneously grab the sample names in the order of the distance matrix and put the matrix into
        # a twoD list and then convert to a numpy array
        self.clade_pcoa_coord_file_path = os.path.join(
            self.clade_output_dir, f'{self.date_time_string}.bray_curtis_{self.profiles_or_samples}_PCoA_coords_{clade}.csv')
        raw_dist_file = read_defined_file_to_list(self.clade_dist_file_path)

        temp_two_d_list = []
        object_names_from_dist_matrix = []
        object_ids_from_dist_matrix = []
        for line in raw_dist_file[1:]:
            temp_elements = line.split('\t')
            object_names_from_dist_matrix.append(temp_elements[0].replace(' ', ''))
            object_ids_from_dist_matrix.append(int(temp_elements[1]))
            temp_two_d_list.append([float(a) for a in temp_elements[2:]])

        dist_as_np_array = np.array(temp_two_d_list)

        sys.stdout.write('\rcalculating PCoA coordinates')

        pcoa_output = pcoa(dist_as_np_array)

        # rename the pcoa dataframe index as the sample names
        pcoa_output.samples['sample'] = object_names_from_dist_matrix
        renamed_pcoa_dataframe = pcoa_output.samples.set_index('sample')

        # now add the variance explained as a final row to the renamed_dataframe
        renamed_pcoa_dataframe = renamed_pcoa_dataframe.append(
            pcoa_output.proportion_explained.rename('proportion_explained'))

        object_ids_from_dist_matrix.append(0)
        if self.profiles_or_samples == 'samples':
            renamed_pcoa_dataframe.insert(loc=0, column='sample_uid', value=object_ids_from_dist_matrix)
            renamed_pcoa_dataframe['sample_uid'] = renamed_pcoa_dataframe['sample_uid'].astype('int')
        else:  # 'profiles'
            renamed_pcoa_dataframe.insert(loc=0, column='analysis_type_uid', value=object_ids_from_dist_matrix)
            renamed_pcoa_dataframe['analysis_type_uid'] = renamed_pcoa_dataframe['analysis_type_uid'].astype('int')

        renamed_pcoa_dataframe.to_csv(self.clade_pcoa_coord_file_path, index=True, header=True, sep=',')

    def _set_data_set_sample_uid_list(self, data_set_sample_uid_list, data_set_uid_list):
        if data_set_sample_uid_list:
            return data_set_sample_uid_list
        else:
            return [dss.id for dss in DataSetSample.objects.filter(data_submission_from__in=data_set_uid_list)]

    def _compute_braycurtis_btwn_obj_pairs(self):
        for obj_one, obj_two in itertools.combinations(list(self.objs_of_clade), 2):
            # let's work to a virtual subsample of 100 000
            obj_one_seq_rel_abundance_dict = self.clade_rs_uid_to_normalised_abund_clade_dict[
                obj_one.id]
            obj_two_seq_rel_abundance_dict = self.clade_rs_uid_to_normalised_abund_clade_dict[
                obj_two.id]

            # for each comparison. Get a set of all of the sequence and convert this to a list.
            set_of_rs_uids = set(list(obj_one_seq_rel_abundance_dict.keys()))
            set_of_rs_uids.update(list(obj_two_seq_rel_abundance_dict.keys()))
            list_of_rs_uids = list(set_of_rs_uids)

            # then iter through the list to get the rel abundances for each of the samples, putting 0 if not found in
            # the sample.
            rs_abundance_list_one = []
            rs_abundance_list_two = []
            for rs_uid in list_of_rs_uids:
                # populate the abundance list cc one
                if rs_uid in obj_one_seq_rel_abundance_dict:
                    rs_abundance_list_one.append(obj_one_seq_rel_abundance_dict[rs_uid])
                else:
                    rs_abundance_list_one.append(0)

                # populate the abundance list cc two
                if rs_uid in obj_two_seq_rel_abundance_dict:
                    rs_abundance_list_two.append(obj_two_seq_rel_abundance_dict[rs_uid])
                else:
                    rs_abundance_list_two.append(0)

            # once you have this we should simply be able to crunch the bray-curtis.
            distance = braycurtis(rs_abundance_list_one, rs_abundance_list_two)

            # these distances can be stored in a dictionary by the 'id1/id2' and 'id2/id1'
            self.clade_within_clade_distances_dict[frozenset({obj_one.id, obj_two.id})] = distance

    def _generate_distance_file(self):
        self.clade_dist_file_as_list = [len(self.objs_of_clade)]
        for obj_outer in self.objs_of_clade:
            temp_at_string = [obj_outer.id]

            for obj_inner in self.objs_of_clade:
                if obj_outer == obj_inner:
                    temp_at_string.append(0)
                else:
                    temp_at_string.append(
                        self.clade_within_clade_distances_dict[frozenset({obj_outer.id, obj_inner.id})])
            self.clade_dist_file_as_list.append(
                '\t'.join([str(distance_item) for distance_item in temp_at_string]))

    def _add_obj_uids_to_dist_file_and_write(self):
        # for the output version lets also append the sample name to each line so that we can see which sample it is
        # it is important that we otherwise work eith the sample ID as the sample names may not be unique.
        dist_with_obj_name = [self.clade_dist_file_as_list[0]]
        list_of_obj_uids = [int(line.split('\t')[0]) for line in self.clade_dist_file_as_list[1:]]
        if self.profiles_or_samples == 'samples':
            objs_of_outputs = list(CladeCollection.objects.filter(id__in=list_of_obj_uids))
            dict_of_obj_id_to_obj_name = {obj.id: obj.data_set_sample_from.name for obj in objs_of_outputs}
        else:  # 'profiles'
            objs_of_outputs = list(AnalysisType.objects.filter(id__in=list_of_obj_uids))
            dict_of_obj_id_to_obj_name = {obj.id: obj.name for obj in objs_of_outputs}

        for line in self.clade_dist_file_as_list[1:]:
            temp_list = []
            obj_id = int(line.split('\t')[0])
            obj_name = dict_of_obj_id_to_obj_name[obj_id]
            temp_list.append(obj_name)
            temp_list.extend(line.split('\t'))
            new_line = '\t'.join(temp_list)
            dist_with_obj_name.append(new_line)
        self.clade_dist_file_as_list = dist_with_obj_name
        write_list_to_destination(self.clade_dist_file_path, self.clade_dist_file_as_list)

    def _write_out_dist_file(self):
        write_list_to_destination(self.clade_dist_file_path, self.clade_dist_file_as_list)

    def _write_output_paths_to_stdout(self):
        print('\n\nBrayCurtis distances and PCoA computation complete. Output files:')
        for output_path in self.output_file_paths:
            print(output_path)

    def _append_output_files_to_output_list(self):
        self.output_file_paths.extend([self.clade_dist_file_path, self.clade_pcoa_coord_file_path])


class SampleBrayCurtisDistPCoACreator(BaseBrayCurtisDistPCoACreator):
    """This will produce braycurtis based clade separated distance matrices and compute PCoAs from these
    matrices. It will allow input of data set ids or as sample uid lists that are comma separated.
    It may be that we can have a bas class here which will allow us to have a bray curtis that is specific
    to types and samples.

    Abbreviations
    -------------
    ds = DataSet
    dss = DataSetSample
    dsss = DataSetSampleSequence
    cc = CladeCollection
    """

    def __init__(
            self, symportal_root_directory, date_time_string=None, data_set_sample_uid_list=None,
            data_set_uid_list=None, call_type=None, output_dir=None, is_sqrt_transf=False):
        super().__init__(
            symportal_root_directory=symportal_root_directory, call_type=call_type, date_time_string=date_time_string, profiles_or_samples='samples')

        self.data_set_sample_uid_list = self._set_data_set_sample_uid_list(
            data_set_sample_uid_list=data_set_sample_uid_list, data_set_uid_list=data_set_uid_list)

        self.cc_list_for_output = CladeCollection.objects.filter(
            data_set_sample_from__in=self.data_set_sample_uid_list)
        self.clades_of_ccs = list(set([a.clade for a in self.cc_list_for_output]))
        self.output_dir = self._set_output_dir(output_dir=output_dir)
        self.is_sqrt_transf = is_sqrt_transf


    def _set_output_dir(self, output_dir):
        if self.call_type == 'stand_alone':
            new_output_dir = os.path.join(
                self.symportal_root_dir, 'outputs', 'ordination', self.date_time_string.replace('.','_'),
                'between_samples')
        else:
            # call_type == 'submission':
            new_output_dir = os.path.join(output_dir, 'between_sample_distances')
        return new_output_dir

    def compute_braycurtis_dists_and_pcoa_coords(self):
        print('\n\nComputing sample pairwise distances and PCoA coordinates using the BrayCurtis method\n')
        for clade_in_question in self.clades_of_ccs:
            self._init_clade_dirs_and_paths(clade_in_question)

            self.objs_of_clade = list(self.cc_list_for_output.filter(clade=clade_in_question))
            if len(self.objs_of_clade) < 2:
                continue
            self._create_rs_uid_to_normalised_abund_dict_for_each_obj_samples()
            self._compute_braycurtis_btwn_obj_pairs()
            self._generate_distance_file()
            self._add_obj_uids_to_dist_file_and_write()
            self._write_out_dist_file()
            self._compute_pcoa_coords(clade=clade_in_question)
            self._append_output_files_to_output_list()
        self._write_output_paths_to_stdout()

    def _init_clade_dirs_and_paths(self, clade_in_question):
        self.clade_output_dir = os.path.join(self.output_dir, clade_in_question)
        self.clade_dist_file_path = os.path.join(
            self.clade_output_dir, f'{self.date_time_string}.bray_curtis_sample_distances_{clade_in_question}.dist')
        os.makedirs(self.clade_output_dir, exist_ok=True)

    def _create_rs_uid_to_normalised_abund_dict_for_each_obj_samples(self):
        # Go through each of the clade collections and create a dict
        # that has key as ref_seq_uid and relative abundance of that sequence
        # we can then store these dict in a dict where the key is the sample ID.
        for clade_col in self.objs_of_clade:
            temp_dict = {}
            list_of_dss_in_cc = list(DataSetSampleSequence.objects.filter(
                clade_collection_found_in=clade_col))

            if self.is_sqrt_transf:
                sqrt_abundances_dict = {
                    dsss.id: math.sqrt(dsss.abundance) for dsss in
                    list_of_dss_in_cc}
                total_seqs_ind_clade_col = sum(sqrt_abundances_dict.values())
                for dsss in list_of_dss_in_cc:
                    temp_dict[dsss.reference_sequence_of.id] = (sqrt_abundances_dict[dsss.id] /
                                                                total_seqs_ind_clade_col) * 10000
            else:
                total_seqs_ind_clade_col = sum([dsss.abundance for dsss in list_of_dss_in_cc])
                for dsss in list_of_dss_in_cc:
                    temp_dict[dsss.reference_sequence_of.id] = (dsss.abundance / total_seqs_ind_clade_col)*10000

            self.clade_rs_uid_to_normalised_abund_clade_dict[clade_col.id] = temp_dict

    @staticmethod
    def _infer_is_dataset_of_datasetsample(smpl_id_list_str, data_set_string):
        if smpl_id_list_str is not None and data_set_string is not None:
            raise RuntimeError('Please input smpl_id_list_str OR data_set_string')
        elif smpl_id_list_str and data_set_string:
            raise RuntimeError('Please input smpl_id_list_str OR data_set_string')
        else:
            if smpl_id_list_str:
                return True
            else:
                return False


class TypeBrayCurtisDistPCoACreator(BaseBrayCurtisDistPCoACreator):
    """This will produce braycurtis based clade separated distance matrices and compute PCoAs from these
    matrices. It will allow input of data set ids or as sample uid lists that are comma separated.
    It may be that we can have a bas class here which will allow us to have a bray curtis that is specific
    to types and samples.

    Abbreviations
    -------------
    ds = DataSet
    dss = DataSetSample
    dsss = DataSetSampleSequence
    cc = CladeCollection
    """

    def __init__(
            self, symportal_root_directory, data_analysis_obj, date_time_string=None, data_set_sample_uid_list=None,
            data_set_uid_list=None, call_type=None, output_dir=None, is_sqrt_transf=False):
        super().__init__(
            symportal_root_directory=symportal_root_directory, call_type=call_type, date_time_string=date_time_string, profiles_or_samples='profiles')

        self.data_set_sample_uid_list = self._set_data_set_sample_uid_list(
            data_set_sample_uid_list=data_set_sample_uid_list, data_set_uid_list=data_set_uid_list)
        self.data_analysis_obj = data_analysis_obj
        self.at_list_for_output = AnalysisType.objects.filter(
            data_analysis_from=self.data_analysis_obj,
            cladecollectiontype__clade_collection_found_in__data_set_sample_from__in=self.data_set_sample_uid_list).distinct()
        self.clades_of_ats = list(set([at.clade for at in self.at_list_for_output]))
        self.output_dir = self._set_output_dir(output_dir=output_dir)
        self.is_sqrt_transf = is_sqrt_transf


    def _set_output_dir(self, output_dir):
        if self.call_type == 'stand_alone':
            new_output_dir = os.path.join(
                self.symportal_root_dir, 'outputs', 'ordination', self.date_time_string.replace('.','_'),
                'between_profiles')
        else:
            # call_type == 'analysis':
            new_output_dir = os.path.join(output_dir, 'between_profile_distances')
        return new_output_dir



    def compute_braycurtis_dists_and_pcoa_coords(self):
        print('\n\nComputing ITS2 type profile pairwise distances and PCoA coordinates using the BrayCurtis method\n')
        for clade_in_question in self.clades_of_ats:
            self._init_clade_dirs_and_paths(clade_in_question)

            self.objs_of_clade = list(self.at_list_for_output.filter(clade=clade_in_question))
            if len(self.objs_of_clade) < 2:
                continue
            self._create_rs_uid_to_normalised_abund_dict_for_each_obj_profiles()
            self._compute_braycurtis_btwn_obj_pairs()
            self._generate_distance_file()
            self._add_obj_uids_to_dist_file_and_write()
            self._write_out_dist_file()
            self._compute_pcoa_coords(clade=clade_in_question)
            self._append_output_files_to_output_list()
        self._write_output_paths_to_stdout()

    def _init_clade_dirs_and_paths(self, clade_in_question):
        self.clade_output_dir = os.path.join(self.output_dir, clade_in_question)
        self.clade_dist_file_path = os.path.join(
            self.clade_output_dir,
            f'{self.date_time_string}.bray_curtis_within_clade_profile_distances_{clade_in_question}.dist')
        os.makedirs(self.clade_output_dir, exist_ok=True)

    def _create_rs_uid_to_normalised_abund_dict_for_each_obj_profiles(self):
        # Go through each of the clade collections and create a dict
        # that has key as ref_seq_uid and relative abundance of that sequence
        # we can then store these dict in a dict where the key is the sample ID.
        for at in self.objs_of_clade:
            ref_seq_uids_of_analysis_type = [int(b) for b in at.ordered_footprint_list.split(',')]

            df = pd.DataFrame(at.get_ratio_list())

            if self.is_sqrt_transf:
                df = general.sqrt_transform_abundance_df(df)

            normalised_abundance_of_divs_dict = {
                ref_seq_uids_of_analysis_type[i]: math.ceil(df[i].mean() * 100000) for
                i in range(len(ref_seq_uids_of_analysis_type))}

            self.clade_rs_uid_to_normalised_abund_clade_dict[at.id] = normalised_abundance_of_divs_dict

