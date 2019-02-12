import os
import subprocess
import sys
from plumbum import local
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from dbApp.models import DataSet, DataSetSample, DataSetSampleSequence, CladeCollection, ReferenceSequence, DataAnalysis
from general import (
    decode_utf8_binary_to_list, create_dict_from_fasta, create_seq_name_to_abundance_dict_from_name_file,
    combine_two_fasta_files, remove_primer_mismatch_annotations_from_fasta, write_list_to_destination,
    read_defined_file_to_list, mafft_align_fasta)
from pickle import dump, load
from multiprocessing import Queue, Manager, Process, current_process
from django import db
import json
from datetime import datetime
from matplotlib.patches import Rectangle  # Rectangle is used despite it being greyed out in pycharm
from matplotlib.collections import PatchCollection
# https://stackoverflow.com/questions/21784641/installation-issue-with-matplotlib-python
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.ioff()
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import re
from skbio.stats.ordination import pcoa
import shutil
import itertools
from scipy.spatial.distance import braycurtis

class BlastnAnalysis:
    def __init__(
            self, input_file_path, output_file_path,
            db_path='/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload/nt', max_target_seqs=1,
            num_threads=1, output_format_string="6 qseqid sseqid staxids evalue pident qcovs staxid stitle ssciname",
            blastn_exec_path='blastn'
    ):

        self.input_file_path = input_file_path
        self.output_file_path = output_file_path
        self.db_path = db_path
        self.output_format_string = output_format_string
        self.max_target_seqs = max_target_seqs
        self.num_threads = num_threads
        self.blastn_exec_path = blastn_exec_path

    def execute(self, pipe_stdout_sterr=True):
        if pipe_stdout_sterr:
            completedProcess = subprocess.run([
                self.blastn_exec_path, '-out', self.output_file_path, '-outfmt', self.output_format_string, '-query',
                self.input_file_path, '-db', self.db_path,
                '-max_target_seqs', f'{self.max_target_seqs}', '-num_threads', f'{self.num_threads}'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            completedProcess = subprocess.run([
                self.blastn_exec_path, '-out', self.output_file_path, '-outfmt', self.output_format_string, '-query',
                self.input_file_path, '-db', self.db_path,
                '-max_target_seqs', f'{self.max_target_seqs}', '-num_threads', f'{self.num_threads}'])
        return completedProcess

    def return_blast_output_as_list(self):
        return read_defined_file_to_list(self.output_format_string)

    def return_blast_results_dict(self):
        blast_output_file_as_list = self.return_blast_output_as_list()
        blast_output_dict = defaultdict(list)
        for line in blast_output_file_as_list:
            blast_output_dict[line.split('\t')[0]].append('\t'.join(line.split('\t')[1:]))
        return blast_output_file_as_list


class MothurAnalysis:

    def __init__(
            self, sequence_collection=None,  input_dir=None, output_dir=None, name=None,
            fastq_gz_fwd_path=None, fastq_gz_rev_path=None,
             name_file_path=None, mothur_execution_path='mothur', auto_convert_fastq_to_fasta=True,
            pcr_fwd_primer=None, pcr_rev_primer=None, pcr_oligo_file_path=None,
            pcr_fwd_primer_mismatch=2, pcr_rev_primer_mismatch=2, pcr_analysis_name=None, num_processors=10,
            stdout_and_sterr_to_pipe=True, tree_file_path=None, group_file_path=None, is_unifrac_analysis=False
            ):

        self._setup_core_attributes(auto_convert_fastq_to_fasta, fastq_gz_fwd_path, fastq_gz_rev_path,
                                    input_dir, mothur_execution_path, name, name_file_path, output_dir,
                                    sequence_collection, num_processors, stdout_and_sterr_to_pipe, is_unifrac_analysis)


        self._setup_pcr_analysis_attributes(pcr_analysis_name, pcr_fwd_primer, pcr_fwd_primer_mismatch,
                                            pcr_oligo_file_path, pcr_rev_primer, pcr_rev_primer_mismatch)




    def perform_weighted_unifrac(self):
        apples = 'asdf'
    def _setup_core_attributes(self, auto_convert_fastq_to_fasta, fastq_gz_fwd_path, fastq_gz_rev_path,
                               input_dir, mothur_execution_path, name, name_file_path, output_dir,
                               sequence_collection, num_processors, stdout_and_sterr_to_pipe, is_unifrac_analysis,
                               tree_path, group_file_path):



        if is_unifrac_analysis:
            self._setup_unifrac_attributes(
                tree_path, group_file_path, name_file_path, input_dir, output_dir,
                mothur_execution_path, num_processors, stdout_and_sterr_to_pipe)
        else:
            self._verify_that_is_either_sequence_collection_or_fastq_pair(fastq_gz_fwd_path, fastq_gz_rev_path,
                                                                          sequence_collection)
            if sequence_collection is not None:
                self._setup_sequence_collection_attribute(auto_convert_fastq_to_fasta, name, sequence_collection)
            elif sequence_collection is None:
                self._setup_fastq_attributes(fastq_gz_fwd_path, fastq_gz_rev_path)
            self._setup_remainder_of_core_attributes(input_dir, mothur_execution_path, name_file_path,
                                                     output_dir, sequence_collection, num_processors, stdout_and_sterr_to_pipe)

    def _setup_unifrac_attributes(self, tree_path, group_file_path, name_file_path, input_dir, output_dir,
                mothur_execution_path, num_processors, stdout_and_sterr_to_pipe):
        self.tree_file_path = tree_path
        self.group_file_path = group_file_path
        self.name_file_path = name_file_path
        if input_dir is None:
            self.input_dir = os.path.dirname(tree_path)
        else:
            self.input_dir = input_dir
        if output_dir is None:
            self.output_dir = os.path.dirname(tree_path)
        else:
            self.output_dir = output_dir
        self.exec_path = mothur_execution_path
        self.processors = num_processors
        self.stdout_and_sterr_to_pipe = stdout_and_sterr_to_pipe



    def _setup_remainder_of_core_attributes(self, input_dir, mothur_execution_path, name_file_path,
                                            output_dir, sequence_collection, num_processors, stdout_and_sterr_to_pipe):
        self.exec_path = mothur_execution_path
        if input_dir is None:
            self.input_dir = os.path.dirname(sequence_collection.file_path)
        else:
            self.input_dir = input_dir
        if output_dir is None:
            self.output_dir = os.path.dirname(sequence_collection.file_path)
        else:
            self.output_dir = input_dir

        self.name_file_path = name_file_path
        self.mothur_batch_file_path = None
        self.processors = num_processors
        # we need to have seperate latest completed process objects for the actual commands and for the summaries
        # this is so that we can still extract useful information housed in the stdout from running the command
        # once the execute... function has been completed. Else, this information is lost due to it being replaced
        # by the completed_process of the summary that is automatically run after each command.
        self.latest_completed_process_command = None
        self.latest_completed_process_summary = None
        self.latest_summary_output_as_list = None
        self.latest_summary_path = None
        self.stdout_and_sterr_to_pipe = stdout_and_sterr_to_pipe
        self.dist_file_path = None
        self.tree_file_path = None

    def _setup_fastq_attributes(self, fastq_gz_fwd_path, fastq_gz_rev_path):
        self.fastq_gz_fwd_path = fastq_gz_fwd_path
        self.fastq_gz_rev_path = fastq_gz_rev_path
        self.sequence_collection = None
        self.fasta_path = None

    def _setup_sequence_collection_attribute(self, auto_convert_fastq_to_fasta, name, sequence_collection):
        self.fastq_gz_fwd_path = None
        self.fastq_gz_rev_path = None
        if sequence_collection.file_type == 'fastq':
            self._convert_to_fasta_or_raise_value_error(auto_convert_fastq_to_fasta, sequence_collection)
        if name is None:
            self.name = sequence_collection.name
        else:
            self.name = name
        self.sequence_collection = sequence_collection
        self.fasta_path = self.sequence_collection.file_path

    def _convert_to_fasta_or_raise_value_error(self, auto_convert_fastq_to_fasta, sequence_collection):
        if auto_convert_fastq_to_fasta:
            print('SequenceCollection must be of type fasta\n. Running SeqeunceCollection.convert_to_fasta.\n')
            sequence_collection.convert_to_fasta()
        else:
            ValueError('SequenceCollection must be of type fasta. You can use the SequenceCollection')

    def _verify_that_is_either_sequence_collection_or_fastq_pair(self, fastq_gz_fwd_path, fastq_gz_rev_path,
                                                                 sequence_collection):
        # if a group_file was given then
        if sequence_collection and (fastq_gz_fwd_path or fastq_gz_rev_path):
            raise ValueError(
                'Please create a MothurAnalysis from either a sequence_collection OR a pair of fastq_gz files.\n'
                'MothurAnalysis.from_pair_of_fastq_gz_files or MothurAnalysis.from_sequence_collection')

    def _setup_pcr_analysis_attributes(self, pcr_analysis_name, pcr_fwd_primer, pcr_fwd_primer_mismatch,
                                       pcr_oligo_file_path, pcr_rev_primer, pcr_rev_primer_mismatch):
        if pcr_analysis_name:
            if pcr_analysis_name.lower() in ['symvar', 'sym_var']:
                self.pcr_fwd_primer = 'GAATTGCAGAACTCCGTGAACC'
                self.rev_primer = 'GAATTGCAGAACTCCGTGAACC',

            elif pcr_analysis_name.lower() in ['laj', 'lajeunesse']:
                self.pcr_fwd_primer = 'GAATTGCAGAACTCCGTG'
                self.pcr_rev_primer = 'CGGGTTCWCTTGTYTGACTTCATGC'
            else:
                raise ValueError(
                    'pcr_analysis_name \'{}\' is not recognised.\nOptions are \'symvar\' or \'lajeunesse\'.'
                )
        else:
            self.pcr_fwd_primer = pcr_fwd_primer
            self.pcr_rev_primer = pcr_rev_primer
        self.pcr_fwd_primer_mismatch = pcr_fwd_primer_mismatch
        self.pcr_rev_primer_mismatch = pcr_rev_primer_mismatch
        self.pcr_oligo_file_path = pcr_oligo_file_path

    # init class methods for the MothurAnalysis
    @classmethod
    def init_from_pair_of_fastq_gz_files(cls, name, fastq_gz_fwd_path, fastq_gz_rev_path,
                                         output_dir=None, mothur_execution_string='mothur', num_processors=10,
                                         stdout_and_sterr_to_pipe=True
                                         ):
        return cls(name=name, sequence_collection=None, mothur_execution_path=mothur_execution_string,
                   input_dir=os.path.dirname(os.path.abspath(fastq_gz_fwd_path)), output_dir=output_dir,
                   fastq_gz_fwd_path=fastq_gz_fwd_path, fastq_gz_rev_path=fastq_gz_rev_path,
                   name_file_path=None, num_processors=num_processors,
                   stdout_and_sterr_to_pipe=stdout_and_sterr_to_pipe)

    @classmethod
    def init_from_sequence_collection(cls, sequence_collection, name=None, input_dir=None,
                                      output_dir=None, mothur_execution_path='mothur',
                                      pcr_fwd_primer=None, pcr_rev_primer=None, pcr_oligo_file_path=None,
                                      pcr_fwd_primer_mismatch=2, pcr_rev_primer_mismatch=2, pcr_analysis_name=None,
                                      num_processors=10, stdout_and_sterr_to_pipe=True):
        return cls(
            name=name, sequence_collection=sequence_collection, input_dir=input_dir,
            output_dir=output_dir, mothur_execution_path=mothur_execution_path, pcr_fwd_primer=pcr_fwd_primer,
            pcr_rev_primer=pcr_rev_primer, pcr_oligo_file_path=pcr_oligo_file_path,
            pcr_fwd_primer_mismatch=pcr_fwd_primer_mismatch, pcr_rev_primer_mismatch=pcr_rev_primer_mismatch,
            pcr_analysis_name=pcr_analysis_name, num_processors=num_processors,
            stdout_and_sterr_to_pipe=stdout_and_sterr_to_pipe
        )

    @classmethod
    def init_for_weighted_unifrac(
            cls, tree_path, group_file_path, name_file_path, input_dir=None,
            output_dir=None, mothur_execution_string='mothur', num_processors=10,
            stdout_and_sterr_to_pipe=True, is_unifrac_analysis=True):

        return cls(
            mothur_execution_path=mothur_execution_string, input_dir=input_dir, output_dir=output_dir,
            name_file_path=name_file_path, num_processors=num_processors,
            stdout_and_sterr_to_pipe=stdout_and_sterr_to_pipe,
            tree_file_path=tree_path, group_file_path=group_file_path, is_unifrac_analysis=is_unifrac_analysis
        )

    # ########################################

    # main mothur commands
    def execute_screen_seqs(self, argument_dictionary):
        """This will perform a mothur screen.seqs.
        Because there are so many arguments taht the screen seqs command can take we will use a dictionary
        to determine how the mothur batch file should be made. The dictionary is very simple: the key should
        be the argument and the value should be the value. e.g.
        argument_dictionary = {'max_length':'500', 'min_length':'100'}
        """
        self._screen_seqs_make_and_write_mothur_batch_file(argument_dictionary)
        self._run_mothur_batch_file_command()
        good_fasta_path = self._screen_seqs_extract_good_output_path()
        self.fasta_path = good_fasta_path
        self._update_sequence_collection_from_fasta_file()
        self.__execute_summary()

    def execute_pcr(self, do_reverse_pcr_as_well=False):
        """This will perform a mothur pcr.seqs analysis.
        if do_reverse_pcr__as_well is true then we will also reverse complement the fasta a perform the
        """

        self._pcr_validate_attributes_are_set()

        self._pcr_make_and_write_oligo_file_if_doesnt_exist()

        self._pcr_make_and_write_mothur_batch_file()

        self._run_mothur_batch_file_command()

        fwd_output_scrapped_fasta_path, fwd_output_good_fasta_path = self._pcr_extract_good_and_scrap_output_paths()

        remove_primer_mismatch_annotations_from_fasta(fwd_output_scrapped_fasta_path)
        remove_primer_mismatch_annotations_from_fasta(fwd_output_good_fasta_path)


        # then we should clean up the output_bad_fasta
        # then reverse complement it
        # then do a pcr on it again using the same oligo set as the first run
        # we should then get the output from that pcr and add it to the previous run
        if do_reverse_pcr_as_well:
            self.fasta_path = fwd_output_scrapped_fasta_path
            self._rev_comp_make_and_write_mothur_batch_file()
            self._run_mothur_batch_file_command()
            self.fasta_path = self._extract_output_path_first_line()
            self._pcr_make_and_write_mothur_batch_file()
            self._run_mothur_batch_file_command()
            rev_output_good_fasta_path = self._pcr_extract_good_and_scrap_output_paths()[1]
            remove_primer_mismatch_annotations_from_fasta(rev_output_good_fasta_path)
            self._make_new_fasta_path_for_fwd_rev_combined(rev_output_good_fasta_path)
            # now create a fasta that is the good fasta from both of the pcrs. this will become the new mothuranalysis fasta.

            combine_two_fasta_files(
                path_one=fwd_output_good_fasta_path,
                path_two=rev_output_good_fasta_path,
                path_for_combined=self.fasta_path
            )
        else:
            self.fasta_path = fwd_output_good_fasta_path
        if self.name_file_path:
            self._update_sequence_collection_from_fasta_name_pair()
        else:
            self._update_sequence_collection_from_fasta_file()

    def execute_make_contigs(self):
        """
        This will use the fastq_gz_fwd_path and fastq_gz_rev_paths to make a .file file that will be used
        as input to mothurs make.contigs command.
        N.B. Although in theory we can use the fastq_gz_fwd_path and the rev path directly as arguments to the mothur.contigs
        there appears to be a bug that doesn't allow this to work. Using a .file file is fine though. The .file file
        is in the format "path_to_file_1 path_to_file_2" i.e the paths only separated by a space.
        :return:
        """
        # create .file file for the fwd fastq pair and the reverse fastq pair
        dot_file_file_path = self._make_contig_make_and_write_out_dot_file()

        self._make_contig_make_and_write_mothur_batch(dot_file_file_path)

        self._run_mothur_batch_file_command()

        self.fasta_path = self._extract_output_path_first_line()

        self._update_sequence_collection_from_fasta_file()

        self.__execute_summary()

    def execute_unique_seqs(self):

        self._unique_seqs_make_and_write_mothur_batch()

        self._run_mothur_batch_file_command()

        self.name_file_path, self.fasta_path = self._extract_output_path_two_lines()

        self._update_sequence_collection_from_fasta_name_pair()

        self.__execute_summary()

    def execute_split_abund(self, abund_cutoff=2):

        self._split_abund_make_and_write_mothur_batch(abund_cutoff)

        self._run_mothur_batch_file_command()

        self.name_file_path, self.fasta_path = self._split_abund_extract_output_path_name_and_fasta()

        self._update_sequence_collection_from_fasta_name_pair()

        self.__execute_summary()

    def execute_dist_seqs(self):
        self._dist_seqs_make_and_write_mothur_batch()

        self._run_mothur_batch_file_command()

        self.dist_file_path = self._extract_output_path_first_line()

    def execute_clearcut(self):
        self._validate_dist_file()
        self._clearcut_make_and_write_mothur_batch()
        self._run_mothur_batch_file_command()
        self.tree_file_path = self._extract_output_path_first_line()

    def execute_weighted_unifrac(self):
        self._weighted_unifrac_make_and_write_mothur_batch()
        self._run_mothur_batch_file_command()
        if self.latest_completed_process_command.returncode == 0:
            sys.stdout.write('\rUnifrac successful')
        else:
            sys.stdout.write('\rERROR: {}'.format(self.latest_completed_process_command.sterr.decode('utf-8')))
        # TODO verify that this is the appropriate method to collect the output
        self.dist_file_path = self._extract_output_path_first_line()


    def __execute_summary(self):
        self._summarise_make_and_write_mothur_batch()
        self._run_mothur_batch_file_command()
        self.latest_summary_path = self._extract_output_path_first_line()
        self.latest_summary_output_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)

    # #####################
    def _weighted_unifrac_make_and_write_mothur_batch(self):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'unifrac.weighted(tree={self.tree_file_path}, group={self.group_file_path}, name={self.name_file_path},'
            f' distance=square, processors={self.processors})'
        ]

        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')

        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _validate_dist_file(self):
        if self.dist_file_path is None:
            raise RuntimeError('A .dist file must exist to run the clearcut analysis')

    def _clearcut_make_and_write_mothur_batch(self):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'clearcut(phylip={self.dist_file_path}, verbose=t)'
        ]

        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')

        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _dist_seqs_make_and_write_mothur_batch(self):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'dist.seqs(fasta={self.fasta_path}, countends=T, output=square)'
        ]

        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')

        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _split_abund_extract_output_path_name_and_fasta(self):
        stdout_string_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            print(stdout_string_as_list[i])
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i + 2], stdout_string_as_list[i + 4]

    def _split_abund_make_and_write_mothur_batch(self, abund_cutoff):
        if self.name_file_path:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'split.abund(fasta={self.fasta_path}, name={self.name_file_path}, cutoff={abund_cutoff})'
            ]
        else:
            raise RuntimeError(
                'Non name_file_path present. '
                'A name file is necessary to be able to assess the abundances of sequences in the .fasta file'
            )
        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _unique_seqs_make_and_write_mothur_batch(self):
        if self.name_file_path:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'unique.seqs(fasta={self.fasta_path}, name={self.name_file_path})'
            ]
        else:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'unique.seqs(fasta={self.fasta_path})'
            ]
        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _summarise_make_and_write_mothur_batch(self):
        if self.name_file_path:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'summary.seqs(fasta={self.fasta_path})'
            ]
        else:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'summary.seqs(fasta={self.fasta_path}, name={self.name_file_path})'
            ]
        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)


    def _make_contig_make_and_write_mothur_batch(self, dot_file_file_path):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'make.contigs(file={dot_file_file_path})'
        ]
        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _make_contig_make_and_write_out_dot_file(self):
        dot_file_file = [f'{self.fastq_gz_fwd_path} {self.fastq_gz_rev_path}']
        dot_file_file_path = os.path.join(self.input_dir, 'fastq_pair.file')
        write_list_to_destination(dot_file_file_path, dot_file_file)
        return dot_file_file_path


    def _make_new_fasta_path_for_fwd_rev_combined(self, rev_output_good_fasta_path):
        self.fasta_path = rev_output_good_fasta_path.replace('.scrap.pcr.rc.pcr', '.pcr.combined')

    def _update_sequence_collection_from_fasta_file(self):
        self.sequence_collection.set_list_of_nucleotide_sequences_from_fasta_or_fastq(self.fasta_path)

    def _update_sequence_collection_from_fasta_name_pair(self):
        self.sequence_collection.generate_sequence_collection_from_fasta_name_pair(self.fasta_path)

    def _rev_comp_make_and_write_mothur_batch_file(self):
        mothur_batch_file = self._make_rev_complement_mothur_batch_file()
        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _make_rev_complement_mothur_batch_file(self):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'reverse.seqs(fasta={self.fasta_path})'
        ]
        return mothur_batch_file

    def _extract_output_path_first_line(self):
        stdout_string_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            print(stdout_string_as_list[i])
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i+1]

    def _extract_output_path_two_lines(self):
        stdout_string_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            print(stdout_string_as_list[i])
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i + 1], stdout_string_as_list[i + 2]

    def _pcr_extract_good_and_scrap_output_paths(self):
        stdout_string_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            print(stdout_string_as_list[i])
            if 'Output File Names' in stdout_string_as_list[i]:
                output_good_fasta_path = stdout_string_as_list[i + 1]
                output_scrapped_fasta_path = stdout_string_as_list[i + 3]
                return output_scrapped_fasta_path, output_good_fasta_path

    def _screen_seqs_extract_good_output_path(self):
        stdout_string_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            print(stdout_string_as_list[i])
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i + 1]

    def _run_mothur_batch_file_command(self):
        if self.stdout_and_sterr_to_pipe:
            self.latest_completed_process_command = subprocess.run(
                [self.exec_path, self.mothur_batch_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        else:
            self.latest_completed_process_command = subprocess.run(
                [self.exec_path, self.mothur_batch_file_path])

    def _run_mothur_batch_file_summary(self, stdout_and_sterr_to_pipe=False):
        if stdout_and_sterr_to_pipe:
            self.latest_completed_process_summary = subprocess.run(
                [self.exec_path, self.mothur_batch_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        else:
            self.latest_completed_process_summary = subprocess.run(
                [self.exec_path, self.mothur_batch_file_path])

    def _pcr_make_and_write_mothur_batch_file(self):
        mothur_batch_file = self._pcr_make_mothur_batch_file()
        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _screen_seqs_make_and_write_mothur_batch_file(self, argument_dictionary):
        mothur_batch_file = self._screen_seqs_make_mothur_batch_file(argument_dictionary)
        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _screen_seqs_create_additional_arguments_string(self, argument_dict):
        individual_argument_strings = []
        for k, v in argument_dict.items():
            if v is not None:
                individual_argument_strings.append(f'{k}={v}')
        return ', '.join(individual_argument_strings)

    def _screen_seqs_make_mothur_batch_file(self, argument_dict):
        additional_arguments_string = self._screen_seqs_create_additional_arguments_string(argument_dict)
        if self.name_file_path:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'screen.seqs(fasta={self.fasta_path}, name={self.name_file_path}, {additional_arguments_string})'
            ]

        else:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'screen.seqs(fasta={self.fasta_path}, {additional_arguments_string})'
            ]
        return mothur_batch_file

    def _pcr_make_mothur_batch_file(self):
        if self.name_file_path:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'pcr.seqs(fasta={self.fasta_path}, name={self.name_file_path}, oligos={self.pcr_oligo_file_path}, '
                f'pdiffs={self.pcr_fwd_primer_mismatch}, rdiffs={self.pcr_rev_primer_mismatch}, processors={self.processors})'
            ]

        else:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'pcr.seqs(fasta={self.fasta_path}, oligos={self.pcr_oligo_file_path}, '
                f'pdiffs={self.pcr_fwd_primer_mismatch}, rdiffs={self.pcr_rev_primer_mismatch}, processors={self.processors})'
            ]
        return mothur_batch_file

    def _pcr_make_and_write_oligo_file_if_doesnt_exist(self):
        if self.pcr_oligo_file_path is None:
            oligo_file = [
                f'forward\t{self.pcr_fwd_primer}',
                f'reverse\t{self.pcr_rev_primer}'
            ]
            self.pcr_oligo_file_path = os.path.join(self.input_dir, 'oligo_file.oligo')
            write_list_to_destination(self.pcr_oligo_file_path, oligo_file)

    def _pcr_validate_attributes_are_set(self):
        sys.stdout.write(f'\nValidating PCR attributes are set\n')
        if self.fasta_path is None:
            raise RuntimeError('Fasta_path is None. A valid fasta_path is required to perform the pcr method.')
        if self.pcr_fwd_primer is None or self.pcr_rev_primer is None:
            if self.pcr_fwd_primer is None and self.pcr_rev_primer is None:
                raise RuntimeError('Please set fwd_primer and rev_primer: ')
            elif self.pcr_fwd_primer is None:
                raise RuntimeError('Please set fwd_primer.')
            elif self.pcr_rev_primer is None:
                raise RuntimeError('Please set fwd_primer.')
        sys.stdout.write(f'\nPCR attributes: OK\n')


class SequenceCollection:
    """ A sequence collection is a set of sequences either generated from a fastq file or from a fasta file.
    It cannot be created directly from binary files or from paired files. As such, to generate a SequenceCollection
    for example from a pair of fastaq.gz files, you would first have to run a mothur contig analysis and create
    the SeqeunceCollection from the resultant fasta file that is generated."""
    def __init__(self, name, path_to_file=None, auto_convert_to_fasta=True):
        self.name = name
        self.file_path = path_to_file
        # self.file_as_list = read_defined_file_to_list(self.file_path)
        self.file_type = self.infer_file_type()
        self.list_of_nucleotide_sequences = None
        self.set_list_of_nucleotide_sequences_from_fasta_or_fastq()
        if auto_convert_to_fasta:
            self.convert_to_fasta()


    def convert_to_fasta(self):
        self.file_path = self.write_out_as_fasta()
        self.file_type = 'fasta'

    def __len__(self):
        return(len(self.list_of_nucleotide_sequences))

    def write_out_as_fasta(self, path_for_fasta_file = None):
        if self.file_type == 'fasta':
            print(f'SequenceCollection is already of type fasta and a fasta file already exists: {self.file_path}')
            return
        if self.file_type == 'fastq':
            if path_for_fasta_file is None:
                fasta_path = self.infer_fasta_path_from_current_fastq_path()
                write_list_to_destination(destination=fasta_path, list_to_write=self.as_fasta())
            else:
                fasta_path = path_for_fasta_file
                write_list_to_destination(destination=fasta_path, list_to_write=self.as_fasta())
            return fasta_path


    def infer_fasta_path_from_current_fastq_path(self):
        return self.file_path.replace('fastq', 'fasta')

    def set_list_of_nucleotide_sequences_from_fasta_or_fastq(self, alt_fasta_path=None):
        """This will generate a list of NucleotideSequence objects.
        It will do this with a fasta or fastq file as the sole input.
        As such no abundance data will be collected for each of the NucleotideSequence objects."""
        if self.file_type == 'fasta':
            self.parse_fasta_file_and_extract_nucleotide_sequence_objects(alternative_fasta_file_path=alt_fasta_path)
        elif self.file_type == 'fastq':
            self.parse_fastq_file_and_extract_nucleotide_sequence_objects()

    def generate_sequence_collection_from_fasta_name_pair(self, name_file_path, fasta_file_path):
        """This will generate a list of NucleotideSequence objects.
        It will do this with a fasta or fastq file as the sole input.
        As such no abundance data will be collected for each of the NucleotideSequence objects."""
        list_of_nucleotide_sequence_objects = []
        fasta_dict = create_dict_from_fasta(fasta_path=fasta_file_path)
        seq_name_to_abundace_dict = create_seq_name_to_abundance_dict_from_name_file(name_file_path)
        for seq_name, seq_sequence in fasta_dict.items():
            list_of_nucleotide_sequence_objects.append(
                NucleotideSequence(sequence=seq_sequence, name=seq_name, abundance=seq_name_to_abundace_dict[seq_name])
            )
        self.list_of_nucleotide_sequences = list_of_nucleotide_sequence_objects

    def parse_fasta_file_and_extract_nucleotide_sequence_objects(self, alternative_fasta_file_path=None):
        list_of_nucleotide_sequence_objects = []
        if alternative_fasta_file_path:
            self.file_path = alternative_fasta_file_path
        fasta_file = read_defined_file_to_list(self.file_path)
        for i in range(0, len(fasta_file), 2):
            list_of_nucleotide_sequence_objects.append(
                NucleotideSequence(sequence=fasta_file[i+1], name=fasta_file[i][1:])
            )
        self.list_of_nucleotide_sequences = list_of_nucleotide_sequence_objects

    def parse_fastq_file_and_extract_nucleotide_sequence_objects(self):
        list_of_nuleotide_sequence_objects = []
        fastq_file_as_list = read_defined_file_to_list(self.file_path)
        for i in range(len(fastq_file_as_list)):
            if i < len(fastq_file_as_list) - 2:
                if self.is_fastq_defline(fastq_file_as_list, i):
                    self.create_new_nuc_seq_object_and_add_to_list(fastq_file_as_list, i, list_of_nuleotide_sequence_objects)
        self.list_of_nucleotide_sequences = list_of_nuleotide_sequence_objects

    def is_fastq_defline(self, fastsq_file, index_value):
        if fastsq_file[index_value].startswith('@') and fastsq_file[index_value + 2][0] == '+':
            return True

    def create_new_nuc_seq_object_and_add_to_list(self, fastq_file_as_list, index_val, list_of_nuleotide_sequence_objects):
        name, sequence = self.get_single_fastq_info_from_fastq_file_by_index(fastq_file_as_list, index_val)
        list_of_nuleotide_sequence_objects.append(NucleotideSequence(sequence=sequence, name=name))

    def get_single_fastq_info_from_fastq_file_by_index(self, fastq_file_as_list, index_val):
        name = fastq_file_as_list[index_val][1:].split(' ')[0]
        sequence = fastq_file_as_list[index_val + 1]
        return name, sequence

    def infer_file_type(self):
        if 'fasta' in self.file_path:
            return 'fasta'
        elif 'fastq' in self.file_path:
            return 'fastq'
        else:
            raise ValueError('Input file used to create the SequenceCollection must be either fasta or fastq')

    def as_fasta(self):
        fasta_file = []
        for seq_obj in self.list_of_nucleotide_sequences:
            fasta_file.extend([f'>{seq_obj.name}', f'{seq_obj.sequence}'])
        return fasta_file


class NucleotideSequence:
    def __init__(self, sequence, name=None, abundance=None):
        self.sequence = sequence
        self.length = len(sequence)
        self.name = name
        self.abundance = abundance


class InitialMothurHandler:
    def __init__(self, sample_fastq_pairs, num_proc):
        self.input_queue_containing_pairs_of_fastq_file_paths = Queue()
        self.worker_manager = Manager()
        self.samples_that_caused_errors_in_qc_mp_list = self.worker_manager.list()
        self.num_proc = num_proc
        self._populate_input_queue(sample_fastq_pairs)

    def _populate_input_queue(self, sample_fastq_pairs):
        for fastq_path_pair in sample_fastq_pairs:
            self.input_queue_containing_pairs_of_fastq_file_paths.put(fastq_path_pair)

        for n in range(self.num_proc):
            self.input_queue_containing_pairs_of_fastq_file_paths.put('STOP')

    def execute_worker_initial_mothur(
            self, data_loading_dataset_object, data_loading_temp_working_directory, data_loading_debug
    ):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        sys.stdout.write('\nPerforming QC\n')
        for n in range(self.num_proc):
            p = Process(target=self._worker_initial_mothur,
                        args=(data_loading_dataset_object, data_loading_temp_working_directory, data_loading_debug)
                        )

            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

    def _worker_initial_mothur(
            self, data_loading_dataset_object, data_loading_temp_working_directory, data_loading_debug
    ):
        """
        This worker performs the pre-MED processing that is primarily mothur-based.
        This QC includes making contigs, screening for ambigous calls (0 allowed) and homopolymer (maxhomop=5)
        Discarding singletons and doublets, in silico PCR. It also checks whether sequences are rev compliment.
        This is all done through the use of an InitialMothurWorker class which in turn makes use of the MothurAnalysis
        class that does the heavy lifting of running the mothur commands in sequence.
        """

        for contigpair in iter(self.input_queue_containing_pairs_of_fastq_file_paths.get, 'STOP'):

            initial_morthur_worker = InitialMothurWorker(
                contig_pair=contigpair,
                data_set_object=data_loading_dataset_object,
                temp_wkd=data_loading_temp_working_directory,
                debug=data_loading_debug
            )

            try:
                initial_morthur_worker.execute()
            except RuntimeError as e:
                self.samples_that_caused_errors_in_qc_mp_list.append(e.sample_name)
        return


class InitialMothurWorker:
    def __init__(self, contig_pair, data_set_object, temp_wkd, debug):
        self.sample_name = contig_pair.split('\t')[0].replace('[dS]', '-')
        self.data_set_sample = DataSetSample.objects.get(
                name=self.sample_name, data_submission_from=data_set_object
            )
        self.cwd = os.path.join(temp_wkd, self.sample_name)
        os.makedirs(self.cwd, exist_ok=True)
        self.pre_med_seq_dump_dir = self.cwd.replace('tempData', 'pre_MED_seqs')
        os.makedirs(self.pre_med_seq_dump_dir, exist_ok=True)
        self.mothur_analysis_object = MothurAnalysis(
            pcr_analysis_name='symvar',
            input_dir=self.cwd,
            output_dir=self.cwd,
            fastq_gz_fwd_path=contig_pair.split('\t')[1],
            fastq_gz_rev_path=contig_pair.split('\t')[2],
            stdout_and_sterr_to_pipe=debug)
        self.debug = debug

    def execute(self):
        sys.stdout.write(f'{self.sample_name}: QC started\n')

        self._do_make_contigs()

        self._set_absolute_num_seqs_after_make_contigs()

        self._do_screen_seqs()

        self._do_unique_seqs()

        self._do_split_abund()

        self._do_unique_seqs()

        self._do_fwd_and_rev_pcr()

        self._do_unique_seqs()

        self._set_unique_num_seqs_after_initial_qc()

        self._set_absolute_num_seqs_after_inital_qc()

        self._save_changes_to_data_set_sample()

        sys.stdout.write(f'{self.sample_name}: Initial mothur complete\n')

        self._write_out_final_name_and_fasta_for_tax_screening()


    def _write_out_final_name_and_fasta_for_tax_screening(self):
        name_file_as_list = read_defined_file_to_list(self.mothur_analysis_object.name_file_path)
        taxonomic_screening_name_file_path = os.path.join(self.cwd, 'name_file_for_tax_screening.names')
        write_list_to_destination(taxonomic_screening_name_file_path, name_file_as_list)
        fasta_file_as_list = read_defined_file_to_list(self.mothur_analysis_object.fasta_path)
        taxonomic_screening_fasta_file_path = os.path.join(self.cwd, 'fasta_file_for_tax_screening.fasta')
        write_list_to_destination(taxonomic_screening_fasta_file_path, fasta_file_as_list)

    def _save_changes_to_data_set_sample(self):
        self.data_set_sample.save()

    def _do_fwd_and_rev_pcr(self):
        self.mothur_analysis_object.execute_pcr(do_reverse_pcr_as_well=True)
        self.check_for_no_seqs_after_pcr_and_raise_runtime_error()

    def _do_split_abund(self):
        self.mothur_analysis_object.execute_split_abund(abund_cutoff=2)
        self.check_for_error_and_raise_runtime_error()

    def _do_unique_seqs(self):
        self.mothur_analysis_object.execute_unique_seqs()
        self.check_for_error_and_raise_runtime_error()

    def _do_screen_seqs(self):
        self.mothur_analysis_object.execute_screen_seqs(argument_dictionary={'maxambig': 0, 'maxhomop': 5})
        self.check_for_error_and_raise_runtime_error()

    def _do_make_contigs(self):
        self.mothur_analysis_object.execute_make_contigs()
        self.check_for_error_and_raise_runtime_error()

    def _set_absolute_num_seqs_after_make_contigs(self):
        number_of_contig_seqs_absolute = len(
            read_defined_file_to_list(self.mothur_analysis_object.latest_summary_path)
        ) - 1
        self.data_set_sample.num_contigs = number_of_contig_seqs_absolute
        sys.stdout.write(
            f'{self.sample_name}: data_set_sample_instance_in_q.num_contigs = {number_of_contig_seqs_absolute}\n')

    def _set_unique_num_seqs_after_initial_qc(self):
        number_of_contig_seqs_unique = len(
            read_defined_file_to_list(self.mothur_analysis_object.latest_summary_path)) - 1
        self.data_set_sample.post_qc_unique_num_seqs = number_of_contig_seqs_unique
        sys.stdout.write(
            f'{self.sample_name}: '
            f'data_set_sample_instance_in_q.post_qc_unique_num_seqs = {number_of_contig_seqs_unique}\n')

    def _set_absolute_num_seqs_after_inital_qc(self):
        last_summary = read_defined_file_to_list(self.mothur_analysis_object.latest_summary_path)
        absolute_count = 0
        for line in last_summary[1:]:
            absolute_count += int(line.split('\t')[6])
        self.data_set_sample.post_qc_absolute_num_seqs = absolute_count
        sys.stdout.write(
            f'{self.sample_name}: data_set_sample_instance_in_q.post_qc_absolute_num_seqs = {absolute_count}\n')

    def check_for_no_seqs_after_pcr_and_raise_runtime_error(self):
        if len(self.mothur_analysis_object.sequence_collection) == 0:
            self.log_qc_error_and_continue(errorreason='No seqs left after PCR')
            raise RuntimeError(sample_name=self.sample_name)

    def check_for_error_and_raise_runtime_error(self):
        for stdout_line in decode_utf8_binary_to_list(
                self.mothur_analysis_object.latest_completed_process_command.stdout
        ):
            if '[WARNING]: Blank fasta name, ignoring read.' in stdout_line:
                self.log_qc_error_and_continue(errorreason='Blank fasta name')
                raise RuntimeError(sample_name=self.sample_name)
            if 'ERROR' in stdout_line:
                self.log_qc_error_and_continue(errorreason='error in inital QC')
                raise RuntimeError(sample_name=self.sample_name)

    def log_qc_error_and_continue(self, errorreason):
        print('Error in processing sample: {}'.format(self.sample_name))
        self.data_set_sample.unique_num_sym_seqs = 0
        self.data_set_sample.absolute_num_sym_seqs = 0
        self.data_set_sample.initial_processing_complete = True
        self.data_set_sample.error_in_processing = True
        self.data_set_sample.error_reason = errorreason
        self._save_changes_to_data_set_sample()


class PotentialSymTaxScreeningHandler:
    """ The purpose of this handler and the executed work is only to get a collection of sequences that will need
    screening against the NCBI database. We also rely on this method to do the blast of our each samples sequences
    against our symclade.fa reference sequences database. We read in the blast.out file in the later functions.
    """
    def __init__(self, samples_that_caused_errors_in_qc_list, checked_samples_list, list_of_samples_names, num_proc):
        self.input_queue = Queue()
        self.manager = Manager()
        self.sub_evalue_sequence_to_num_sampes_found_in_mp_dict = self.manager.dict()
        self.sub_evalue_nucleotide_sequence_to_clade_mp_dict = self.manager.dict()
        self.error_samples_mp_list = self.manager.list(samples_that_caused_errors_in_qc_list)
        self.checked_samples_mp_list = self.manager.list(checked_samples_list)
        self.list_of_sample_names = list_of_samples_names
        self.num_proc = num_proc
        self._load_input_queue()

    def _load_input_queue(self):
        # load up the input q
        for sample_name in self.list_of_sample_names:
            self.input_queue.put(sample_name)
        # load in the STOPs
        for n in range(self.num_proc):
            self.input_queue.put('STOP')

    def execute_potential_sym_tax_screening(self, data_loading_temp_working_directory, data_loading_path_to_symclade_db, data_loading_debug):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        sys.stdout.write('\nPerforming QC\n')
        for n in range(self.num_proc):
            new_process = Process(
                target=self._potential_sym_tax_screening_worker,
                args=(data_loading_temp_working_directory, data_loading_path_to_symclade_db, data_loading_debug))

            all_processes.append(new_process)
            new_process.start()
        for process in all_processes:
            process.join()

    def _potential_sym_tax_screening_worker(
            self, data_loading_temp_working_directory, data_loading_path_to_symclade_db, data_loading_debug):
        """
        :param input_q: The multiprocessing queue that holds a list of the sample names
        :param e_val_collection_dict: This is a managed dictionary where key is a nucleotide sequence that has:
        1 - provided a match in the blast analysis
        2 - is of suitable size
        3 - but has an evalue match below the cuttof
        the value of the dict is an int that represents how many samples this nucleotide sequence was found in
        :param err_smpl_list: This is a managed list containing sample names of the samples that had errors during the
        initial mothur qc and therefore don't require taxonomic screening performed on them
        :param checked_samples: This is a list of sample names for samples that were found to contain only
        Symbiodinium sequences or have already had all potential Symbiodinium sequences screened and so don't
        require any further taxonomic screening
        :return: Whilst this doesn't return anything a number of objects are picked out in each of the local
        working directories for use in the workers that follow this one.
        """
        for sample_name in iter(self.input_queue.get, 'STOP'):

            # If the sample gave an error during the inital mothur then we don't consider it here.
            if sample_name in self.error_samples_mp_list:
                continue

            # A sample will be in this list if we have already performed this worker on it and none of its sequences
            # gave matches to the symClade database at below the evalue threshold
            if sample_name in self.checked_samples_mp_list:
                continue

            taxonomic_screening_worker = PotentialSymTaxScreeningWorker(
                sample_name=sample_name, wkd=data_loading_temp_working_directory,
                path_to_symclade_db=data_loading_path_to_symclade_db, debug=data_loading_debug,
                checked_samples_mp_list=self.checked_samples_mp_list,
                e_val_collection_mp_dict=self.sub_evalue_sequence_to_num_sampes_found_in_mp_dict,
                sub_evalue_nucleotide_sequence_to_clade_mp_dict=self.sub_evalue_nucleotide_sequence_to_clade_mp_dict)

            taxonomic_screening_worker.execute()


class PotentialSymTaxScreeningWorker:
    def __init__(
            self, sample_name, wkd, path_to_symclade_db, debug, e_val_collection_mp_dict,
            checked_samples_mp_list, sub_evalue_nucleotide_sequence_to_clade_mp_dict):
        self.sample_name = sample_name
        self.cwd = os.path.join(wkd, self.sample_name)
        self.fasta_file_path = os.path.join(self.cwd, 'fasta_file_for_tax_screening.fasta')
        self.fasta_dict = create_dict_from_fasta(fasta_path=self.fasta_file_path)
        self.name_file_path = os.path.join(self.cwd, 'name_file_for_tax_screening.names')
        self.name_dict = {a.split('\t')[0]: a for a in read_defined_file_to_list(self.name_file_path)}
        self.path_to_symclade_db = path_to_symclade_db
        self.debug = debug
        # This is a managed dictionary where key is a nucleotide sequence that has:
        # 1 - provided a match in the blast analysis
        # 2 - is of suitable size
        # 3 - but has an evalue match below the cuttof
        # the value of the dict is an int that represents how many samples this nucleotide sequence was found in
        self.e_val_collection_mp_dict = e_val_collection_mp_dict
        # This dictionary will be used outside of the multiprocessing to append the clade of a given sequences
        # that is being added to the symClade reference database
        self.sub_evalue_nucleotide_sequence_to_clade_mp_dict = sub_evalue_nucleotide_sequence_to_clade_mp_dict
        # The potential_non_symbiodinium_sequences_list is used to see if there are any samples, that don't have
        # any potential non symbiodnium sequences. i.e. only definite symbiodinium sequences. These samples are added to
        # the checked list and are not checked in following iterations to speed things up.
        self.potential_non_symbiodinium_sequences_list = []
        self.sequence_name_to_clade_dict = None
        self.blast_output_as_list = None
        self.already_processed_blast_seq_result = []
        # this is a managed list that holds the names of samples from which no sequences were thrown out from
        # it will be used in downstream processes.
        self.checked_samples_mp_list = checked_samples_mp_list

    def execute(self):
        sys.stdout.write(f'{self.sample_name}: verifying seqs are Symbiodinium and determining clade\n')

        blastn_analysis = BlastnAnalysis(
            input_file_path=self.fasta_file_path,
            output_file_path=os.path.join(self.cwd, 'blast.out'), db_path=self.path_to_symclade_db,
            output_format_string="6 qseqid sseqid staxids evalue pident qcovs")

        if self.debug:
            blastn_analysis.execute(pipe_stdout_sterr=False)
        else:
            blastn_analysis.execute(pipe_stdout_sterr=True)

        sys.stdout.write(f'{self.sample_name}: BLAST complete\n')

        self.blast_output_as_list = blastn_analysis.return_blast_output_as_list()

        self._if_debug_warn_if_blast_out_empty_or_low_seqs()

        self.sequence_name_to_clade_dict = {
            blast_out_line.split('\t')[0]: blast_out_line.split('\t')[1][-1] for blast_out_line in self.blast_output_as_list}

        self._add_seqs_with_no_blast_match_to_non_sym_list()

        # NB blast results sometimes return several matches for the same seq.
        # as such we will use the already_processed_blast_seq_resulst to make sure that we only
        # process each sequence once.
        self._identify_and_allocate_non_sym_and_sub_e_seqs()

        if not self.potential_non_symbiodinium_sequences_list:
            self.checked_samples_mp_list.append(self.sample_name)

    def _identify_and_allocate_non_sym_and_sub_e_seqs(self):
        for line in self.blast_output_as_list:
            name_of_current_sequence = line.split('\t')[0]
            if name_of_current_sequence in self.already_processed_blast_seq_result:
                continue
            self.already_processed_blast_seq_result.append(name_of_current_sequence)
            identity = float(line.split('\t')[4])
            coverage = float(line.split('\t')[5])

            # noinspection PyPep8,PyBroadException
            # here we are looking for sequences to add to the non_symbiodinium_sequence_list
            # if a sequence fails at any of our if statements it will be added to the non_symbiodinium_sequence_list
            try:
                evalue_power = int(line.split('\t')[3].split('-')[1])
                # With the smallest sequences i.e. 185bp it is impossible to get above the 100 threshold
                # even if there is an exact match. As such, we sould also look at the match identity and coverage
                if evalue_power < 100:  # evalue cut off, collect sequences that don't make the cut
                    self._if_ident_cov_size_good_add_seq_to_non_sym_list_and_eval_dict(
                        coverage, identity, name_of_current_sequence
                    )
            except:
                # here we weren't able to extract the evalue_power for some reason.
                self._if_ident_cov_size_good_add_seq_to_non_sym_list_and_eval_dict(
                    coverage, identity, name_of_current_sequence
                )

    def _if_ident_cov_size_good_add_seq_to_non_sym_list_and_eval_dict(self, coverage, identity, name_of_current_sequence):
        """This method will add nucleotide sequences that gave blast matches but that were below the evalue
        or indentity and coverage thresholds to the evalue collection dict and
        record how many samples that sequence was found in.
        It will also take into account the size thresholds that would normally
        happen later in the code during further mothur qc.
        Finally, it will also populate a nucleotide sequence to clade dictionary that will be used outside of
        the MPing to append the clade of a given sequences that is being added to the symClade reference database.
        The potential_non_symbiodinium_sequences_list is used to see if there are any samples, that don't have
        any potential non symbiodnium sequences. i.e. only definite symbiodinium sequences. These samples are added to
        the checked list and are not checked in following iterations to speed things up.
        """

        if identity < 80 or coverage < 95:
            # incorporate the size cutoff here that would normally happen in the further mothur qc later in the code
            self.potential_non_symbiodinium_sequences_list.append(name_of_current_sequence)
            if 184 < len(self.fasta_dict[name_of_current_sequence]) < 310:
                if self.fasta_dict[name_of_current_sequence] in self.e_val_collection_mp_dict.keys():
                    self.e_val_collection_mp_dict[self.fasta_dict[name_of_current_sequence]] += 1
                else:
                    self.e_val_collection_mp_dict[self.fasta_dict[name_of_current_sequence]] = 1
                    self.sub_evalue_nucleotide_sequence_to_clade_mp_dict[
                        self.fasta_dict[name_of_current_sequence]
                    ] = self.sequence_name_to_clade_dict[name_of_current_sequence]

    def _add_seqs_with_no_blast_match_to_non_sym_list(self):
        sequences_with_no_blast_match_as_set = set(self.fasta_dict.keys()) - set(self.sequence_name_to_clade_dict.keys())
        self.potential_non_symbiodinium_sequences_list.extend(list(sequences_with_no_blast_match_as_set))
        sys.stdout.write(
            f'{self.sample_name}: {len(sequences_with_no_blast_match_as_set)} sequences thrown out '
            f'initially due to being too divergent from reference sequences\n')

    def _if_debug_warn_if_blast_out_empty_or_low_seqs(self):
        if self.debug:
            if not self.blast_output_as_list:
                print(f'WARNING blast output file is empty for {self.sample_name}')
            else:
                if len(self.blast_output_as_list) < 10:
                    print(
                        f'WARNING blast output file for {self.sample_name} '
                        f'is only {len(self.blast_output_as_list)} lines long')


class SymNonSymTaxScreeningHandler:
    def __init__(
            self, data_loading_samples_that_caused_errors_in_qc_mp_list, data_loading_list_of_samples_names,
            data_loading_num_proc):
        self.sample_name_mp_input_queue = Queue()
        self.sym_non_sym_mp_manager = Manager()
        self.samples_that_caused_errors_in_qc_mp_list = self.sym_non_sym_mp_manager.list(
            data_loading_samples_that_caused_errors_in_qc_mp_list
        )
        self.non_symbiodinium_sequences_list = self.sym_non_sym_mp_manager.list()
        self.num_proc = data_loading_num_proc
        self._populate_input_queue(data_loading_list_of_samples_names)

    def _populate_input_queue(self, data_loading_list_of_samples_names):
        for sample_name in data_loading_list_of_samples_names:
            self.sample_name_mp_input_queue.put(sample_name)
        for n in range(self.num_proc):
            self.sample_name_mp_input_queue.put('STOP')

    def execute_sym_non_sym_tax_screening(self, data_loading_temp_working_directory, data_loading_dataset_object, non_symb_and_size_violation_base_dir_path, data_loading_pre_med_sequence_output_directory_path, data_loading_debug):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        sys.stdout.write('\nPerforming QC\n')

        for n in range(self.num_proc):
            p = Process(target=self._sym_non_sym_tax_screening_worker, args=(data_loading_temp_working_directory, data_loading_dataset_object, non_symb_and_size_violation_base_dir_path, data_loading_pre_med_sequence_output_directory_path, data_loading_debug))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

    def _sym_non_sym_tax_screening_worker(self, data_loading_temp_working_directory, data_loading_dataset_object, data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path, data_loading_pre_med_sequence_output_directory_path, data_loading_debug):

        for sample_name in iter(self.sample_name_mp_input_queue.get, 'STOP'):
            if sample_name in self.samples_that_caused_errors_in_qc_mp_list:
                continue

            sym_non_sym_tax_screening_worker_object = SymNonSymTaxScreeningWorker(
                data_loading_temp_working_directory=data_loading_temp_working_directory,
                data_loading_dataset_object=data_loading_dataset_object, sample_name=sample_name,
                data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path=
                data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path,
                data_loading_pre_med_sequence_output_directory_path=data_loading_pre_med_sequence_output_directory_path,
                data_loading_debug=data_loading_debug
            )

            try:
                sym_non_sym_tax_screening_worker_object.execute()
            except RuntimeError as e:
                self.samples_that_caused_errors_in_qc_mp_list.append(e.sample_name)


class SymNonSymTaxScreeningWorker:
    def __init__(
            self, data_loading_temp_working_directory, sample_name, data_loading_dataset_object,
            data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path,
            data_loading_pre_med_sequence_output_directory_path, data_loading_debug
    ):
        self._init_core_class_attributes(
            data_loading_dataset_object, data_loading_temp_working_directory, sample_name, data_loading_debug)
        self._init_fasta_name_blast_and_clade_dict_attributes()
        self._init_non_sym_and_size_violation_output_paths(
            data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path)
        self._init_sym_non_size_violation_output_paths(data_loading_pre_med_sequence_output_directory_path)
        self._init_sets_for_categorizing_sequences()
        self._init_qc_meta_info_counters()

    def _init_qc_meta_info_counters(self):
        self.absolute_number_of_non_sym_sequences = 0
        self.absolute_number_of_sym_no_size_violation_sequences = 0
        self.absolute_number_of_sym_size_violation_sequences = 0

    def _init_sets_for_categorizing_sequences(self):
        self.non_symbiodinium_sequence_name_set_for_sample = set()
        self.sym_size_violation_sequence_name_set_for_sample = set()
        self.sym_no_size_violation_sequence_name_set_for_sample = set()

    def _init_sym_non_size_violation_output_paths(self, data_loading_pre_med_sequence_output_directory_path):
        pre_med_sequence_output_directory_for_sample_path = os.path.join(
            data_loading_pre_med_sequence_output_directory_path, self.sample_name)
        os.makedirs(pre_med_sequence_output_directory_for_sample_path, exist_ok=True)
        self.pre_med_fasta_path = os.path.join(
            pre_med_sequence_output_directory_for_sample_path, f'pre_med_seqs_{self.sample_name}.fasta'
        )
        self.pre_med_names_file_path = os.path.join(
            pre_med_sequence_output_directory_for_sample_path, f'pre_med_seqs_{self.sample_name}.names'
        )

    def _init_non_sym_and_size_violation_output_paths(
            self, data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path):
        non_symbiodiniaceae_and_size_violation_directory_for_sample_path = os.path.join(
            data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path, self.sample_name)
        os.makedirs(non_symbiodiniaceae_and_size_violation_directory_for_sample_path, exist_ok=True)
        self.non_symbiodiniaceae_seqs_fasta_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.sample_name}_non_symbiodiniaceae_sequences.fasta')
        self.non_symbiodiniaceae_seqs_names_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.sample_name}_non_symbiodiniaceae_sequences.names')
        self.symbiodiniaceae_size_violation_seqs_fasta_output_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.sample_name}_symbiodiniaceae_size_violation_sequences.fasta'
        )
        self.symbiodiniaceae_size_violation_seqs_names_output_path = os.path.join(
            non_symbiodiniaceae_and_size_violation_directory_for_sample_path,
            f'{self.sample_name}_symbiodiniaceae_size_violation_sequences.names')

    def _init_fasta_name_blast_and_clade_dict_attributes(self):
        fasta_file_path = os.path.join(self.cwd, 'fasta_file_for_tax_screening.fasta')
        self.fasta_dict = create_dict_from_fasta(fasta_path=fasta_file_path)
        name_file_path = os.path.join(self.cwd, 'name_file_for_tax_screening.names')
        self.name_dict = {a.split('\t')[0]: a for a in read_defined_file_to_list(name_file_path)}
        blast_output_path = os.path.join(self.cwd, 'blast.out')
        self.blast_dict = {blast_line.split('\t')[0]: blast_line for blast_line in
                           read_defined_file_to_list(blast_output_path)}
        self.sequence_name_to_clade_dict = {
            blast_out_line.split('\t')[0]: blast_out_line.split('\t')[1][-1] for
            blast_out_line in self.blast_dict.values()
        }

    def _init_core_class_attributes(
            self, data_loading_dataset_object, data_loading_temp_working_directory, sample_name, data_loading_debug):
        self.sample_name = sample_name
        self.cwd = os.path.join(data_loading_temp_working_directory, sample_name)
        self.dataset_object = DataSetSample.objects.get(
            name=sample_name, data_submission_from=data_loading_dataset_object
        )
        self.debug = data_loading_debug

    def execute(self):
        """This method completes the pre-med quality control.
        Having gone through the initial mothur qc, and having screened for potential symbiodiniaceae sequences
        (if running on the remote system), this method now identifies
        1 - the non symbiodiniaceae sequences and writes them out to the output dir
        2 - identifies the symbiodiniaceae sequences that violate our size range thresholds and writes them out
        also to the output dir
        3 - and finally the symbiodiniaceae sequences that do not violate our size range thresholds (sequences that
        will be carried through into med decomposition). These are also written out both as clade seperated
        (one redundant fasta for each clade) in the temp working directory (for med processing), and as a non clade separated .name
        .file set (pre-med-seqs) in the output directory.
        This method also populates all of the dataset qc metadata attributes accordingly.
        """
        self._identify_and_write_non_sym_seqs_in_sample()

        self._get_size_violation_and_non_size_violations_seq_sets()

        self._write_out_size_violation_seqs()

        if self.sym_no_size_violation_sequence_name_set_for_sample:
            self._write_out_no_size_violation_seqs()
        else:
            self._log_dataset_attr_and_raise_runtime_error()

        self._associate_qc_meta_info_to_data_sheet()

    def _write_out_no_size_violation_seqs(self):
        self._write_out_no_size_violation_seqs_to_pre_med_dirs()
        clades_of_non_violation_seqs = self._get_set_of_clades_represented_by_no_size_violation_seqs()
        self._write_out_no_size_violation_seqs_redundant_fasta_clade_separated(clades_of_non_violation_seqs)

    def _write_out_no_size_violation_seqs_redundant_fasta_clade_separated(self, clades_of_non_violation_seqs):
        for clade_of_sequences_to_write_out in clades_of_non_violation_seqs:
            sequence_names_of_clade = self._get_sequence_names_of_clade_for_no_size_violation_sequences(
                clade_of_sequences_to_write_out
            )

            sample_clade_fasta_path = os.path.join(
                self.cwd,
                clade_of_sequences_to_write_out,
                f'seqs_for_med_{self.sample_name}_clade_{clade_of_sequences_to_write_out}.redundant.fasta'
            )
            with open(sample_clade_fasta_path, 'w') as f:
                for sequence_name in sequence_names_of_clade:
                    sequence_counter = 0
                    for i in range(len(self.name_dict[sequence_name].split('\t')[1].split(','))):
                        f.write(f'>{self.sample_name}_{sequence_name}_{sequence_counter}\n')
                        f.write(f'{self.fasta_dict[sequence_name]}\n')
                        sequence_counter += 1

            self._if_debug_check_length_of_deuniqued_fasta(sample_clade_fasta_path)

    def _if_debug_check_length_of_deuniqued_fasta(self, sample_clade_fasta_path):
        if self.debug:
            deuniqued_fasta = read_defined_file_to_list(sample_clade_fasta_path)
            if deuniqued_fasta:
                if len(deuniqued_fasta) < 100:
                    print(f'{self.sample_name}: WARNING the dequniqed fasta is < {len(deuniqued_fasta)} lines')
            else:
                print(f'{self.sample_name}: ERROR deuniqued fasta is empty')

    def _get_sequence_names_of_clade_for_no_size_violation_sequences(self, clade_of_sequences_to_write_out):
        sequence_names_of_clade = [
            sequence_name for sequence_name in
            self.sym_no_size_violation_sequence_name_set_for_sample if
            self.sequence_name_to_clade_dict[sequence_name] == clade_of_sequences_to_write_out
        ]
        return sequence_names_of_clade

    def _get_set_of_clades_represented_by_no_size_violation_seqs(self):
        clades_of_non_violation_seqs = set(
            [
                clade_value for sequence_name, clade_value in self.blast_dict.items()
                if sequence_name in self.sym_no_size_violation_sequence_name_set_for_sample
            ]
        )
        return clades_of_non_violation_seqs

    def _write_out_no_size_violation_seqs_to_pre_med_dirs(self):
        self._write_out_no_size_violation_fasta_to_pre_med_dir()
        self._write_out_no_size_violation_names_file_to_pre_med_dir()

    def _write_out_no_size_violation_names_file_to_pre_med_dir(self):
        with open(self.pre_med_names_file_path, 'w') as f:
            for sequence_name in list(self.sym_no_size_violation_sequence_name_set_for_sample):
                f.write(f'{self.name_dict[sequence_name]}\n')
                self.absolute_number_of_sym_no_size_violation_sequences += len(
                    self.name_dict[sequence_name].split('\t')[1].split(','))

    def _write_out_no_size_violation_fasta_to_pre_med_dir(self):
        with open(self.pre_med_fasta_path, 'w') as f:
            for sequence_name in list(self.sym_no_size_violation_sequence_name_set_for_sample):
                f.write(f'>{sequence_name}\n')
                f.write(f'{self.fasta_dict[sequence_name]}\n')

    def _write_out_size_violation_seqs(self):
        if self.sym_size_violation_sequence_name_set_for_sample:
            self._write_out_size_violation_fasta()

            self._write_out_size_violation_names_file()

    def _write_out_size_violation_names_file(self):
        with open(self.symbiodiniaceae_size_violation_seqs_names_output_path, 'w') as f:
            for sequence_name in list(self.sym_size_violation_sequence_name_set_for_sample):
                f.write(f'{self.name_dict[sequence_name]}\n')
                self.absolute_number_of_sym_size_violation_sequences += len(
                    self.name_dict[sequence_name].split('\t')[1].split(','))

    def _write_out_size_violation_fasta(self):
        with open(self.symbiodiniaceae_size_violation_seqs_fasta_output_path, 'w') as f:
            for sequence_name in list(self.sym_size_violation_sequence_name_set_for_sample):
                f.write(f'>{sequence_name}\n')
                f.write(f'{self.fasta_dict[sequence_name]}\n')

    def _get_size_violation_and_non_size_violations_seq_sets(self):
        for query_sequence_name, blast_line in self.blast_dict.items():
            if query_sequence_name not in self.non_symbiodinium_sequence_name_set_for_sample:
                if 184 < len(self.fasta_dict[query_sequence_name]) < 310:
                    self.sym_no_size_violation_sequence_name_set_for_sample.add(query_sequence_name)
                else:
                    self.sym_size_violation_sequence_name_set_for_sample.add(query_sequence_name)

    def _associate_qc_meta_info_to_data_sheet(self):
        self._associate_non_sym_seq_attributes_to_dataset()
        self._associate_sym_seq_no_size_violation_attributes_to_dataset()
        self._associate_sym_seq_size_violation_attributes_to_dataset()
        self.dataset_object.initial_processing_complete = True
        self.dataset_object.save()
        print(f'{self.sample_name}: pre-med QC complete')

    def _associate_sym_seq_size_violation_attributes_to_dataset(self):
        self.dataset_object.size_violation_absolute = (
                self.dataset_object.post_qc_absolute_num_seqs -
                self.dataset_object.absolute_num_sym_seqs -
                self.dataset_object.non_sym_absolute_num_seqs
        )
        print(f'{self.sample_name}: size_violation_absolute = {self.dataset_object.size_violation_absolute}')
        self.dataset_object.size_violation_unique = (
                self.dataset_object.post_qc_unique_num_seqs -
                self.dataset_object.unique_num_sym_seqs -
                self.dataset_object.non_sym_unique_num_seqs
        )
        print(f'{self.sample_name}: size_violation_unique = {self.dataset_object.size_violation_unique}')

    def _associate_sym_seq_no_size_violation_attributes_to_dataset(self):
        self.dataset_object.unique_num_sym_seqs = len(self.sym_no_size_violation_sequence_name_set_for_sample)
        self.dataset_object.absolute_num_sym_seqs = self.absolute_number_of_sym_no_size_violation_sequences
        print(f'{self.sample_name}: unique_num_sym_seqs = {self.dataset_object.unique_num_sym_seqs}')
        print(f'{self.sample_name}: absolute_num_sym_seqs = {self.absolute_number_of_sym_no_size_violation_sequences}')

    def _associate_non_sym_seq_attributes_to_dataset(self):
        self.dataset_object.non_sym_unique_num_seqs = len(self.non_symbiodinium_sequence_name_set_for_sample)
        self.dataset_object.non_sym_absolute_num_seqs = self.absolute_number_of_non_sym_sequences
        print(f'{self.sample_name}: non_sym_unique_num_seqs = {self.dataset_object.non_sym_unique_num_seqs}')
        print(f'{self.sample_name}: non_sym_absolute_num_seqs = {self.absolute_number_of_non_sym_sequences}')

    def _log_dataset_attr_and_raise_runtime_error(self):
        # if there are no symbiodiniaceae sequenes then log error and associate meta info
        print(f'{self.sample_name}: QC error.\n No symbiodiniaceae sequences left in sample after pre-med QC.')
        self.dataset_object.non_sym_unique_num_seqs = len(self.non_symbiodinium_sequence_name_set_for_sample)
        self.dataset_object.non_sym_absolute_num_seqs = self.absolute_number_of_non_sym_sequences
        self.dataset_object.size_violation_absolute = self.absolute_number_of_sym_size_violation_sequences
        self.dataset_object.size_violation_unique = len(self.sym_size_violation_sequence_name_set_for_sample)
        self.dataset_object.unique_num_sym_seqs = 0
        self.dataset_object.absolute_num_sym_seqs = 0
        self.dataset_object.initial_processing_complete = True
        self.dataset_object.error_in_processing = True
        self.dataset_object.error_reason = 'No symbiodiniaceae sequences left in sample after pre-med QC'
        self.dataset_object.save()
        raise RuntimeError(sample_name=self.sample_name)

    def _identify_and_write_non_sym_seqs_in_sample(self):
        self._add_seqs_with_no_blast_match_to_non_sym_list()
        self._identify_non_sym_seqs_from_below_match_threshold()
        self._write_out_non_sym_fasta_and_names_files_for_sample()

    def _write_out_non_sym_fasta_and_names_files_for_sample(self):
        if self.non_symbiodinium_sequence_name_set_for_sample:
            self._write_out_non_sym_fasta_for_sample()
            self._write_out_non_sym_names_file_for_sample()

    def _write_out_non_sym_names_file_for_sample(self):
        with open(self.non_symbiodiniaceae_seqs_names_path, 'w') as f:
            for sequence_name in list(self.non_symbiodinium_sequence_name_set_for_sample):
                f.write(f'{self.name_dict[sequence_name]}\n')
                self.absolute_number_of_non_sym_sequences += len(self.name_dict[sequence_name].split('\t')[1].split(','))

    def _write_out_non_sym_fasta_for_sample(self):
        with open(self.non_symbiodiniaceae_seqs_fasta_path, 'w') as f:
            for sequence_name in list(self.non_symbiodinium_sequence_name_set_for_sample):
                f.write(f'>{sequence_name}\n')
                f.write(f'{self.fasta_dict[sequence_name]}\n')

    def _add_seqs_with_no_blast_match_to_non_sym_list(self):
        sequences_with_no_blast_match_as_set = set(self.fasta_dict.keys()) - set(self.sequence_name_to_clade_dict.keys())
        self.non_symbiodinium_sequence_name_set_for_sample.update(list(sequences_with_no_blast_match_as_set))
        sys.stdout.write(
            f'{self.sample_name}: {len(sequences_with_no_blast_match_as_set)} sequences thrown out '
            f'initially due to being too divergent from reference sequences\n')

    def _identify_non_sym_seqs_from_below_match_threshold(self):
        """ This method is where the non_symbiodinium sequences for the sample are identified. If they fall below
        the given identity, coverage and evalues for their match to the symClade db then they are conisdered
        non-symbiodinum. This does not take into size at all. We will screen and report size seperately.
        """
        for blast_sequence_name, blast_line in self.blast_dict.items():
            # noinspection PyPep8,PyBroadException
            try:
                evalue_power = int(blast_line.split('\t')[3].split('-')[1])
                # With the smallest sequences i.e. 185bp it is impossible to get above the 100 threshold
                # even if there is an exact match. As such, we sould also look at the match identity and coverage
                if evalue_power < 100:  # evalue cut off, collect sequences that don't make the cut
                    self._if_low_identity_and_coverage_add_to_non_symb_set(
                        name_of_current_sequence=blast_sequence_name, blast_line=blast_line
                    )
            except:
                # here we weren't able to extract the evalue_power for some reason.
                self._if_low_identity_and_coverage_add_to_non_symb_set(
                    name_of_current_sequence=blast_sequence_name, blast_line=blast_line
                )

    def _if_low_identity_and_coverage_add_to_non_symb_set(self, name_of_current_sequence, blast_line):
        percentage_identity = float(blast_line.split('\t')[4])
        percentage_coverage = float(blast_line.split('\t')[5])
        if percentage_identity < 80 or percentage_coverage < 95:
            # incorporate the size cutoff here that would normally happen in the further mothur qc later in the code
            self.non_symbiodinium_sequence_name_set_for_sample.add(name_of_current_sequence)


class PerformMEDHandler:
    def __init__(self, data_loading_list_of_samples_names, data_loading_temp_working_directory, data_loading_num_proc):
        # need to get list of the directories in which to perform the MED
        # we want to get a list of the .
        self.list_of_samples = data_loading_list_of_samples_names
        self.temp_working_directory = data_loading_temp_working_directory
        self.num_proc = data_loading_num_proc
        self.list_of_redundant_fasta_paths = []
        self._populate_list_of_redundant_fasta_paths()
        self.input_queue_of_redundant_fasta_paths = Queue()
        self._populate_input_queue_of_redundant_fasta_paths()
        self.list_of_med_result_dirs = [
            os.path.join(os.path.dirname(path_to_redundant_fasta), 'MEDOUT') for
            path_to_redundant_fasta in self.list_of_redundant_fasta_paths]

    def execute_perform_med_worker(self, data_loading_debug, data_loading_path_to_med_padding_executable, data_loading_path_to_med_decompoase_executable):
        all_processes = []

        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        for n in range(self.num_proc):
            p = Process(target=self._deunique_worker, args=(
                data_loading_debug, data_loading_path_to_med_padding_executable, data_loading_path_to_med_decompoase_executable))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

    def _populate_list_of_redundant_fasta_paths(self):
        for dirpath, dirnames, files in os.walk(self.temp_working_directory):
            for file_name in files:
                if file_name.endswith('redundant.fasta'):
                    self.list_of_redundant_fasta_paths.append(os.path.join(dirpath, file_name))

    def _populate_input_queue_of_redundant_fasta_paths(self):
        for redundant_fasta_path in self.list_of_redundant_fasta_paths:
            self.input_queue_of_redundant_fasta_paths.put(redundant_fasta_path)
        for n in range(self.num_proc):
            self.input_queue_of_redundant_fasta_paths.put('STOP')

    def _deunique_worker(self, data_loading_debug, data_loading_path_to_med_padding_executable, data_loading_path_to_med_decompose_executable):
        for redundant_fata_path in iter(self.input_queue_of_redundant_fasta_paths.get, 'STOP'):

            perform_med_worker_instance = PerformMEDWorker(redundant_fata_path, data_loading_debug, data_loading_path_to_med_padding_executable, data_loading_path_to_med_decompose_executable)

            perform_med_worker_instance.execute()


class PerformMEDWorker:
    def __init__(self, redundant_fasta_path, data_loading_path_to_med_padding_executable, data_loading_debug, data_loading_path_to_med_decompose_executable):
        self.redundant_fasta_path_unpadded = redundant_fasta_path
        self.redundant_fasta_path_padded = self.redundant_fasta_path_unpadded.replace('.fasta', '.padded.fasta')
        self.cwd = os.path.dirname(self.redundant_fasta_path_unpadded)
        self.sample_name = self.cwd.split('/')[-2]
        self.debug = data_loading_debug
        self.path_to_med_padding_executable = data_loading_path_to_med_padding_executable
        self.path_to_med_decompose_executable = data_loading_path_to_med_decompose_executable
        self.med_output_dir = os.path.join(os.path.dirname(self.redundant_fasta_path_unpadded), 'MEDOUT')
        os.makedirs(self.med_output_dir, exist_ok=True)
        self.med_m_value = self._get_med_m_value()

    def execute(self):
        # TODO check this business about how MED determines the sample name, i.e. whether the underscores matter.
        sys.stdout.write(f'{self.sample_name}: starting MED analysis\n')
        sys.stdout.write(f'{self.sample_name}: padding sequences\n')
        subprocess.run([
            self.path_to_med_padding_executable,
            '-o', self.redundant_fasta_path_padded,
            self.redundant_fasta_path_unpadded], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sys.stdout.write(f'{self.sample_name}: decomposing\n')
        self._do_the_decomposition()
        sys.stdout.write(f'{self.sample_name}: MED analysis complete\n')

    def _do_the_decomposition(self):
        if not self.debug:
            subprocess.run(
                [self.path_to_med_decompose_executable, '-M', str(self.med_m_value), '--skip-gexf-files',
                 '--skip-gen-figures',
                 '--skip-gen-html', '--skip-check-input', '-o',
                 self.med_output_dir, self.redundant_fasta_path_padded], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elif self.debug:
            subprocess.run(
                [self.path_to_med_decompose_executable, '-M', str(self.med_m_value), '--skip-gexf-files',
                 '--skip-gen-figures',
                 '--skip-gen-html',
                 '--skip-check-input', '-o',
                 self.med_output_dir, self.redundant_fasta_path_padded])

    def _get_med_m_value(self):
        # Define MED M value dynamically.
        # The M value is a cutoff that looks at the abundance of the most abundant unique sequence in a node
        # if the abundance is lower than M then the node is discarded
        # we have been working recently with an M that equivaltes to 0.4% of 0.004. This was
        # calculated when working with a modelling project where I was subsampling to 1000 sequences. In this
        # scenario the M was set to 4.
        # We should also take care that M doesn't go below 4, so we should use a max choice for the M
        num_of_seqs_to_decompose = len(read_defined_file_to_list(self.redundant_fasta_path_padded)) / 2
        return max(4, int(0.004 * num_of_seqs_to_decompose))


class DataSetSampleCreatorWorker:
    """This class will be responsible for handling a set of med outputs. Objects will be things like the directory,
    the count table, number of samples, number of nodes, these sorts of things."""
    def __init__(self, med_output_directory,
                 data_set_sample_creator_handler_ref_seq_sequence_to_ref_seq_id_dict,
                 data_set_sample_creator_handler_ref_seq_uid_to_ref_seq_name_dict, data_loading_dataset_obj):
        self.output_directory = med_output_directory
        self.sample_name = self.output_directory.split('/')[-4]
        self.clade = self.output_directory.split('/')[-3]
        self.nodes_list_of_nucleotide_sequences = []
        self._populate_nodes_list_of_nucleotide_sequences()
        self.num_med_nodes = len(self.nodes_list_of_nucleotide_sequences)
        self.node_sequence_name_to_ref_seq_id = {}
        self.ref_seq_sequence_to_ref_seq_id_dict = data_set_sample_creator_handler_ref_seq_sequence_to_ref_seq_id_dict
        self.ref_seq_uid_to_ref_seq_name_dict = data_set_sample_creator_handler_ref_seq_uid_to_ref_seq_name_dict
        self.node_abundance_df = pd.read_csv(
            os.path.join(self.output_directory, 'MATRIX-COUNT.txt'), delimiter='\t', header=0)
        self.total_num_sequences = sum(self.node_abundance_df.iloc[0])
        self.dataset_sample_object = DataSetSample.objects.get(
            data_submission_from=data_loading_dataset_obj, name=self.sample_name)
        self.clade_collection_object = None

    def _populate_nodes_list_of_nucleotide_sequences(self):
        node_file_path = os.path.join(self.output_directory, 'NODE-REPRESENTATIVES.fasta')
        try:
            node_file_as_list = read_defined_file_to_list(node_file_path)
        except FileNotFoundError:
            raise RuntimeError(med_output_directory=self.output_directory)

        for i in range(0, len(node_file_as_list), 2):
            node_seq_name = node_file_as_list[i].split('|')[0][1:]
            node_seq_abundance = int(node_file_as_list[i].split('|')[1].split(':')[1])
            node_seq_sequence = node_file_as_list[i+1].replace('-','')
            self.nodes_list_of_nucleotide_sequences.append(
                NucleotideSequence(name=node_seq_name, abundance=node_seq_abundance, sequence=node_seq_sequence))

    def execute(self):
        for node_nucleotide_sequence_object in self.nodes_list_of_nucleotide_sequences:
            if not self._assign_node_sequence_to_existing_ref_seq(node_nucleotide_sequence_object):
                self._assign_node_sequence_to_new_ref_seq(node_nucleotide_sequence_object)

        if self._two_or_more_nodes_associated_to_the_same_reference_sequence():
            self._make_associations_and_abundances_in_node_abund_df_unique_again()

        self._update_data_set_sample_med_qc_meta_data_and_clade_totals()

        self._if_more_than_200_seqs_create_clade_collection()

        self._create_data_set_sample_sequences()

    def _create_data_set_sample_sequences(self):
        if self._we_made_a_clade_collection():
            self._create_data_set_sample_sequences_with_clade_collection()
        else:
            self._create_data_set_sample_sequences_without_clade_collection()

    def _create_data_set_sample_sequences_without_clade_collection(self):
        data_set_sample_sequence_list = []
        for node_nucleotide_sequence_object in self.nodes_list_of_nucleotide_sequences:
            associated_ref_seq_id = self.node_sequence_name_to_ref_seq_id[node_nucleotide_sequence_object.name]
            associated_ref_seq_object = ReferenceSequence.objects.get(id=associated_ref_seq_id)
            df_index_label = self.node_abundance_df.index.values.tolist()[0]
            dss = DataSetSampleSequence(reference_sequence_of=associated_ref_seq_object,
                                        abundance=self.node_abundance_df.at[
                                            df_index_label, node_nucleotide_sequence_object.name],
                                        data_set_sample_from=self.dataset_sample_object)
            data_set_sample_sequence_list.append(dss)
        DataSetSampleSequence.objects.bulk_create(data_set_sample_sequence_list)

    def _create_data_set_sample_sequences_with_clade_collection(self):
        data_set_sample_sequence_list = []
        associated_ref_seq_uid_as_str_list = []
        for node_nucleotide_sequence_object in self.nodes_list_of_nucleotide_sequences:
            associated_ref_seq_id = self.node_sequence_name_to_ref_seq_id[node_nucleotide_sequence_object.name]
            associated_ref_seq_uid_as_str_list.append(str(associated_ref_seq_id))
            associated_ref_seq_object = ReferenceSequence.objects.get(id=associated_ref_seq_id)
            df_index_label = self.node_abundance_df.index.values.tolist()[0]
            dss = DataSetSampleSequence(
                reference_sequence_of=associated_ref_seq_object,
                clade_collection_found_in=self.clade_collection_object,
                abundance=self.node_abundance_df.at[df_index_label, node_nucleotide_sequence_object.name],
                data_set_sample_from=self.dataset_sample_object)
            data_set_sample_sequence_list.append(dss)
        # Save all of the newly created dss
        DataSetSampleSequence.objects.bulk_create(data_set_sample_sequence_list)
        self.clade_collection_object.footprint = ','.join(associated_ref_seq_uid_as_str_list)
        self.clade_collection_object.save()

    def _we_made_a_clade_collection(self):
        return self.total_num_sequences > 200

    def _if_more_than_200_seqs_create_clade_collection(self):
        if self.total_num_sequences > 200:
            sys.stdout.write(
                f'\n{self.sample_name} clade {self.clade}: '
                f'{self.total_num_sequences} sequences. Creating CladeCollection_object\n')
            new_cc = CladeCollection(clade=self.clade, data_set_sample_from=self.dataset_sample_object)
            new_cc.save()
            self.clade_collection_object = new_cc
        else:
            sys.stdout.write(
                f'\n{self.sample_name} clade {self.clade}: {self.total_num_sequences} sequences. '
                f'Insufficient sequence to create a CladeCollection_object\n')

    def _update_data_set_sample_med_qc_meta_data_and_clade_totals(self):
        self.dataset_sample_object.post_med_absolute += self.total_num_sequences
        self.dataset_sample_object.post_med_unique += len(self.node_abundance_df.iloc[0])
        # Update the cladal_seq_totals
        cladal_seq_abundance_counter = [int(a) for a in json.loads(self.dataset_sample_object.cladal_seq_totals)]
        clade_index = list('ABCDEFGHI').index(self.clade)
        cladal_seq_abundance_counter[clade_index] = self.total_num_sequences
        self.dataset_sample_object.cladal_seq_totals = json.dumps([str(a) for a in cladal_seq_abundance_counter])
        self.dataset_sample_object.save()

    def _two_or_more_nodes_associated_to_the_same_reference_sequence(self):
        """Multiple nodes may be assigned to the same reference sequence. We only want to create one
        data set sample sequence object per reference sequence. As such we need to consolidate abundances
        """
        return len(set(self.node_sequence_name_to_ref_seq_id.values())) != len(self.node_sequence_name_to_ref_seq_id.keys())

    def _make_associations_and_abundances_in_node_abund_df_unique_again(self):
        list_of_non_unique_ref_seq_uids = [
            ref_seq_uid for ref_seq_uid, count in
            Counter(self.node_sequence_name_to_ref_seq_id.values()).items() if count > 1]
        for non_unique_ref_seq_uid in list_of_non_unique_ref_seq_uids:

            node_names_to_be_consolidated = self._get_list_of_node_names_that_need_consolidating(non_unique_ref_seq_uid)

            summed_abund_of_nodes_of_ref_seq = self._get_summed_abundances_of_the_nodes_to_be_consolidated(
                node_names_to_be_consolidated)

            self._del_all_but_first_of_non_unique_nodes_from_df_and_node_to_ref_seq_dict(node_names_to_be_consolidated)

            self._update_node_name_abund_in_df(node_names_to_be_consolidated, summed_abund_of_nodes_of_ref_seq)

    def _get_list_of_node_names_that_need_consolidating(self, non_unique_ref_seq_uid):
        node_names_to_be_consolidated = [
            node_name for node_name in self.node_sequence_name_to_ref_seq_id.keys() if
            self.node_sequence_name_to_ref_seq_id[node_name] == non_unique_ref_seq_uid]
        return node_names_to_be_consolidated

    def _update_node_name_abund_in_df(self, node_names_to_be_consolidated, summed_abund_of_nodes_of_ref_seq):
        df_index_name = self.node_abundance_df.index.values.tolist()[0]
        self.node_abundance_df.at[df_index_name, node_names_to_be_consolidated[0]] = summed_abund_of_nodes_of_ref_seq

    def _del_all_but_first_of_non_unique_nodes_from_df_and_node_to_ref_seq_dict(self, node_names_to_be_consolidated):
        # now delete columns for all but the first of the node_names_to_be_consolidated from df
        self.node_abundance_df = self.node_abundance_df.drop(columns=node_names_to_be_consolidated[1:])
        # and from the node_names_to_ref_seq_id dict
        for node_name in node_names_to_be_consolidated[1:]:
            del self.node_sequence_name_to_ref_seq_id[node_name]

    def _get_summed_abundances_of_the_nodes_to_be_consolidated(self, node_names_to_be_consolidated):
        summed_abund_of_nodes_of_ref_seq = sum(self.node_abundance_df[node_names_to_be_consolidated].iloc[0])
        return summed_abund_of_nodes_of_ref_seq

    def _assign_node_sequence_to_existing_ref_seq(self, node_nucleotide_sequence_object):
        """ We use this to look to see if there is an equivalent ref_seq Sequence for the sequence in question
        This take into account whether the seq_in_q could be a subset or super set of one of the
        ref_seq.sequences.
        Will return false if no ref_seq match is found
        """
        if self._node_sequence_exactly_matches_reference_sequence_sequence(node_nucleotide_sequence_object):
            return self._associate_node_seq_to_ref_seq_by_exact_match_and_return_true(node_nucleotide_sequence_object)
        elif self._node_sequence_matches_reference_sequence_sequence_plus_adenine(node_nucleotide_sequence_object):  # This was a seq shorter than refseq but we can associate
            return self._associate_node_seq_to_ref_seq_by_adenine_match_and_return_true(node_nucleotide_sequence_object)
        else:
            return self._search_for_super_set_match_and_associate_if_found_else_return_false(
                node_nucleotide_sequence_object)

    def _node_sequence_exactly_matches_reference_sequence_sequence(self, node_nucleotide_sequence_object):
        return node_nucleotide_sequence_object.sequence in self.ref_seq_sequence_to_ref_seq_id_dict

    def _node_sequence_matches_reference_sequence_sequence_plus_adenine(self, node_nucleotide_sequence_object):
        return 'A' + node_nucleotide_sequence_object.sequence in self.ref_seq_sequence_to_ref_seq_id_dict

    def _search_for_super_set_match_and_associate_if_found_else_return_false(self, node_nucleotide_sequence_object):
        # or if the seq in question is bigger than a refseq sequence and is a super set of it
        # In either of these cases we should consider this a match and use the refseq matched to.
        # This might be very coputationally expensive but lets give it a go
        for ref_seq_sequence in self.ref_seq_sequence_to_ref_seq_id_dict.keys():
            if node_nucleotide_sequence_object.sequence in ref_seq_sequence or \
                    ref_seq_sequence in node_nucleotide_sequence_object.sequence:
                # Then this is a match
                self.node_sequence_name_to_ref_seq_id[node_nucleotide_sequence_object.name] = \
                self.ref_seq_sequence_to_ref_seq_id_dict[ref_seq_sequence]
                name_of_reference_sequence = self.ref_seq_uid_to_ref_seq_name_dict[
                    self.ref_seq_sequence_to_ref_seq_id_dict[ref_seq_sequence]]
                self._print_succesful_association_details_to_stdout(node_nucleotide_sequence_object,
                                                                    name_of_reference_sequence)
                return True
        return False

    def _associate_node_seq_to_ref_seq_by_adenine_match_and_return_true(self, node_nucleotide_sequence_object):
        self.node_sequence_name_to_ref_seq_id[
            node_nucleotide_sequence_object.name] = self.ref_seq_sequence_to_ref_seq_id_dict[
            'A' + node_nucleotide_sequence_object.sequence]
        name_of_reference_sequence = self.ref_seq_uid_to_ref_seq_name_dict[
            self.ref_seq_sequence_to_ref_seq_id_dict['A' + node_nucleotide_sequence_object.sequence]]
        self._print_succesful_association_details_to_stdout(node_nucleotide_sequence_object, name_of_reference_sequence)
        return True

    def _associate_node_seq_to_ref_seq_by_exact_match_and_return_true(self, node_nucleotide_sequence_object):
        self.node_sequence_name_to_ref_seq_id[
            node_nucleotide_sequence_object.name] = self.ref_seq_sequence_to_ref_seq_id_dict[
            node_nucleotide_sequence_object.sequence]
        name_of_reference_sequence = self.ref_seq_uid_to_ref_seq_name_dict[
            self.ref_seq_sequence_to_ref_seq_id_dict[node_nucleotide_sequence_object.sequence]]
        self._print_succesful_association_details_to_stdout(node_nucleotide_sequence_object, name_of_reference_sequence)
        return True

    def _print_succesful_association_details_to_stdout(self, node_nucleotide_sequence_object, name_of_reference_sequence):
        sys.stdout.write(f'\r{self.sample_name} clade {self.clade}:\n'
                         f'Assigning MED node {node_nucleotide_sequence_object.name} '
                         f'to existing reference sequence {name_of_reference_sequence}')

    def _assign_node_sequence_to_new_ref_seq(self, node_nucleotide_sequence_object):
        new_ref_seq = ReferenceSequence(clade=self.clade, sequence=node_nucleotide_sequence_object.sequence)
        new_ref_seq.save()
        new_ref_seq.name = str(new_ref_seq.id)
        new_ref_seq.save()
        self.ref_seq_sequence_to_ref_seq_id_dict[new_ref_seq.sequence] = new_ref_seq.id
        self.node_sequence_name_to_ref_seq_id[node_nucleotide_sequence_object.name] = new_ref_seq.id
        self.ref_seq_uid_to_ref_seq_name_dict[new_ref_seq.id] = new_ref_seq.name

        sys.stdout.write(f'\r{self.sample_name} clade {self.clade}:\n'
                         f'Assigning MED node {node_nucleotide_sequence_object.name} '
                         f'to new reference sequence {new_ref_seq.name}')


class DataSetSampleCreatorHandler:
    """This class will be where we run the code for creating reference sequences, data set sample sequences and
    clade collections."""
    def __init__(self):
        # dictionaries to save us having to do lots of database look ups
        self.ref_seq_uid_to_ref_seq_name_dict = {
            ref_seq.id: ref_seq.name for ref_seq in ReferenceSequence.objects.all()}
        self.ref_seq_sequence_to_ref_seq_id_dict = {
            ref_seq.sequence: ref_seq.id for ref_seq in ReferenceSequence.objects.all()}


    def execute_data_set_sample_creation(self, data_loading_list_of_med_output_directories, data_loading_debug, data_loading_dataset_object):
        for med_output_directory in data_loading_list_of_med_output_directories:
            try:
                med_output_object = DataSetSampleCreatorWorker(
                    med_output_directory=med_output_directory,
                    data_loading_dataset_obj=data_loading_dataset_object,
                    data_set_sample_creator_handler_ref_seq_sequence_to_ref_seq_id_dict=
                    self.ref_seq_sequence_to_ref_seq_id_dict,
                    data_set_sample_creator_handler_ref_seq_uid_to_ref_seq_name_dict=
                    self.ref_seq_uid_to_ref_seq_name_dict)
            except RuntimeError as e:
                print(f'{e.med_output_directory}: File not found during DataSetSample creation.')
                continue
            if data_loading_debug:
                if med_output_object.num_med_nodes < 10:
                    print(f'{med_output_directory}: WARNING node file contains only {med_output_object.num_med_nodes} sequences.')
            sys.stdout.write(
                f'\n\nPopulating {med_output_object.sample_name} with clade {med_output_object.clade} sequences\n')


class SequenceCountTableCreator:
    """ This is essentially broken into two parts. The first part goes through all of the DataSetSamples from
    the DataSets of the output and collects abundance information. The second part then puts this abundance
    information into a dataframe for both the absoulte and the relative abundance.

    """
    def __init__(self, call_type, output_dir, data_set_uids_to_output_as_comma_sep_string, num_proc, sorted_sample_uid_list=None, analysis_obj_id=None, time_date_str=None, output_user=None):
        self._init_core_vars(
            analysis_obj_id, call_type, data_set_uids_to_output_as_comma_sep_string, num_proc,
            output_dir, output_user, sorted_sample_uid_list, time_date_str)
        self._init_seq_abundance_collection_objects()
        self._init_vars_for_putting_together_the_dfs()
        self._init_output_paths()

    def _init_core_vars(self, analysis_obj_id, call_type, data_set_uids_to_output_as_comma_sep_string, num_proc,
                        output_dir, output_user, sorted_sample_uid_list, time_date_str):
        self.num_proc = num_proc
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.sorted_sample_uid_list = sorted_sample_uid_list
        self.analysis_obj_id = analysis_obj_id
        if time_date_str:
            self.time_date_str = time_date_str
        else:
            self.time_date_str = str(datetime.now()).replace(' ', '_').replace(':', '-')
        self.call_type = call_type
        self.output_user = output_user
        self.clade_list = list('ABCDEFGHI')
        uids_of_data_sets_to_output = [int(a) for a in data_set_uids_to_output_as_comma_sep_string.split(',')]
        self.data_set_objects_to_output = DataSet.objects.filter(id__in=uids_of_data_sets_to_output)
        self.ref_seqs_in_datasets = ReferenceSequence.objects.filter(
            datasetsamplesequence__data_set_sample_from__data_submission_from__in=self.data_set_objects_to_output).distinct()
        set_of_clades_found = {ref_seq.clade for ref_seq in self.ref_seqs_in_datasets}
        self.ordered_list_of_clades_found = [clade for clade in self.clade_list if clade in set_of_clades_found]
        self.list_of_data_set_sample_objects = DataSetSample.objects.filter(
            data_submission_from__in=self.data_set_objects_to_output)

    def _init_seq_abundance_collection_objects(self):
        """Output objects from first worker to be used by second worker"""
        self.dss_id_to_list_of_dsss_objects_dict_mp_dict = None
        self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict = None
        self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict = None
        # this is the list that we will use the self.annotated_dss_name_to_cummulative_rel_abund_mp_dict to create
        # it is a list of the ref_seqs_ordered first by clade then by abundance.
        self.clade_abundance_ordered_ref_seq_list = []

    def _init_vars_for_putting_together_the_dfs(self):
        # variables concerned with putting together the dataframes
        self.dss_id_to_pandas_series_results_list_dict = None
        self.output_df_absolute = None
        self.output_df_relative = None
        self.output_seqs_fasta_as_list = []

    def _init_output_paths(self):
        self.output_paths_list = []
        if self.analysis_obj_id:
            data_analysis_obj = DataAnalysis.objects.get(id=self.analysis_obj_id)
            self.path_to_seq_output_df_absolute = os.path.join(
                self.output_dir,
                f'{self.analysis_obj_id}_{data_analysis_obj.name}_{self.time_date_str}.seqs.absolute.txt')
            self.path_to_seq_output_df_relative = os.path.join(
                self.output_dir,
                f'{self.analysis_obj_id}_{data_analysis_obj.name}_{self.time_date_str}.seqs.relative.txt')

            self.output_fasta_path = os.path.join(self.output_dir,
                                                  f'{self.analysis_obj_id}_{data_analysis_obj.name}_{self.time_date_str}.seqs.fasta')

        else:
            self.path_to_seq_output_df_absolute = os.path.join(self.output_dir,
                                                               f'{self.time_date_str}.seqs.absolute.txt')
            self.path_to_seq_output_df_relative = os.path.join(self.output_dir,
                                                               f'{self.time_date_str}.seqs.relative.txt')
            self.output_fasta_path = os.path.join(self.output_dir, f'{self.time_date_str}.seqs.fasta')

    def execute_output(self):
        self._collect_abundances_for_creating_the_output()

        self._generate_sample_output_series()

        self._create_ordered_output_dfs_from_series()

        self._add_uids_for_seqs_to_dfs()

        self._append_meta_info_to_df()

        self._write_out_dfs_and_fasta()

    def _write_out_dfs_and_fasta(self):
        self.output_df_absolute.to_csv(self.path_to_seq_output_df_absolute, sep="\t")
        self.output_paths_list.append(self.path_to_seq_output_df_absolute)
        self.output_df_relative.to_csv(self.path_to_seq_output_df_relative, sep="\t")
        self.output_paths_list.append(self.path_to_seq_output_df_relative)
        # we created the fasta above.
        write_list_to_destination(self.output_fasta_path, self.output_seqs_fasta_as_list)
        self.output_paths_list.append(self.output_fasta_path)
        print('\nITS2 sequence output files:')
        for path_item in self.output_path_list:
            print(path_item)

    def _append_meta_info_to_df(self):
        # Now append the meta infromation for each of the data_sets that make up the output contents
        # this is information like the submitting user, what the uids of the datasets are etc.
        # There are several ways that this can be called.
        # it can be called as part of the submission: call_type = submission
        # part of an analysis output: call_type = analysis
        # or stand alone: call_type = 'stand_alone'
        # we should have an output for each scenario
        if self.call_type == 'submission':
            self._append_meta_info_to_df_submission()
        elif self.call_type == 'analysis':
            self._append_meta_info_to_df_analysis()
        else:
            # call_type=='stand_alone'
            self._append_meta_info_to_df_stand_alone()

    def _append_meta_info_to_df_submission(self):
        data_set_object = self.data_set_objects_to_output[0]
        # there will only be one data_set object
        meta_info_string_items = [
            f'Output as part of data_set submission ID: {data_set_object.id}; '
            f'submitting_user: {data_set_object.submitting_user}; '
            f'time_stamp: {data_set_object.time_stamp}']
        temp_series = pd.Series(meta_info_string_items, index=[list(self.output_df_absolute)[0]],
                                name='meta_info_summary')
        self.output_df_absolute = self.output_df_absolute.append(temp_series)
        self.output_df_relative = self.output_df_relative.append(temp_series)

    def _append_meta_info_to_df_analysis(self):
        data_analysis_obj = DataAnalysis.objects.get(id=self.analysis_obj_id)
        num_data_set_objects_as_part_of_analysis = len(data_analysis_obj.list_of_data_set_uids.split(','))
        meta_info_string_items = [
            f'Output as part of data_analysis ID: {data_analysis_obj.id}; '
            f'Number of data_set objects as part of analysis = {num_data_set_objects_as_part_of_analysis}; '
            f'submitting_user: {data_analysis_obj.submitting_user}; time_stamp: {data_analysis_obj.time_stamp}']
        temp_series = pd.Series(meta_info_string_items, index=[list(self.output_df_absolute)[0]],
                                name='meta_info_summary')
        self.output_df_absolute = self.output_df_absolute.append(temp_series)
        self.output_df_relative = self.output_df_relative.append(temp_series)
        for data_set_object in self.data_set_objects_to_output:
            data_set_meta_list = [
                f'Data_set ID: {data_set_object.id}; '
                f'Data_set name: {data_set_object.name}; '
                f'submitting_user: {data_set_object.submitting_user}; '
                f'time_stamp: {data_set_object.time_stamp}']

            temp_series = pd.Series(data_set_meta_list, index=[list(self.output_df_absolute)[0]], name='data_set_info')
            self.output_df_absolute = self.output_df_absolute.append(temp_series)
            self.output_df_relative = self.output_df_relative.append(temp_series)

    def _append_meta_info_to_df_stand_alone(self):
        meta_info_string_items = [
            f'Stand_alone output by {self.output_user} on {self.time_date_str}; '
            f'Number of data_set objects as part of output = {len(self.data_set_objects_to_output)}']
        temp_series = pd.Series(meta_info_string_items, index=[list(self.output_df_absolute)[0]],
                                name='meta_info_summary')
        self.output_df_absolute = self.output_df_absolute.append(temp_series)
        self.output_df_relative = self.output_df_relative.append(temp_series)
        for data_set_object in self.data_set_objects_to_output:
            data_set_meta_list = [
                f'Data_set ID: {data_set_object.id}; '
                f'Data_set name: {data_set_object.name}; '
                f'submitting_user: {data_set_object.submitting_user}; '
                f'time_stamp: {data_set_object.time_stamp}']
            temp_series = pd.Series(data_set_meta_list, index=[list(self.output_df_absolute)[0]], name='data_set_info')
            self.output_df_absolute = self.output_df_absolute.append(temp_series)
            self.output_df_relative = self.output_df_relative.append(temp_series)

    def _add_uids_for_seqs_to_dfs(self):
        """Now add the UID for each of the sequences"""
        sys.stdout.write('\nGenerating accession and fasta\n')
        reference_sequences_in_data_sets_no_name = ReferenceSequence.objects.filter(
            datasetsamplesequence__data_set_sample_from__data_submission_from__in=self.list_of_data_set_sample_objects,
            has_name=False).distinct()
        reference_sequences_in_data_sets_has_name = ReferenceSequence.objects.filter(
            datasetsamplesequence__data_set_sample_from__data_submission_from__in=self.list_of_data_set_sample_objects,
            has_name=True).distinct()
        no_name_dict = {rs.id: rs.sequence for rs in reference_sequences_in_data_sets_no_name}
        has_name_dict = {rs.name: (rs.id, rs.sequence) for rs in reference_sequences_in_data_sets_has_name}
        accession_list = []
        num_cols = len(list(self.output_df_relative))
        for i, col_name in enumerate(list(self.output_df_relative)):
            sys.stdout.write('\rAppending accession info and creating fasta {}: {}/{}'.format(col_name, i, num_cols))
            if col_name in self.clade_abundance_ordered_ref_seq_list:
                if col_name[-2] == '_':
                    col_name_id = int(col_name[:-2])
                    accession_list.append(str(col_name_id))
                    self.output_seqs_fasta_as_list.append('>{}'.format(col_name))
                    self.output_seqs_fasta_as_list.append(no_name_dict[col_name_id])
                else:
                    col_name_tup = has_name_dict[col_name]
                    accession_list.append(str(col_name_tup[0]))
                    self.output_seqs_fasta_as_list.append('>{}'.format(col_name))
                    self.output_seqs_fasta_as_list.append(col_name_tup[1])
            else:
                accession_list.append(np.nan)
        temp_series = pd.Series(accession_list, name='seq_accession', index=list(self.output_df_relative))
        self.output_df_absolute = self.output_df_absolute.append(temp_series)
        self.output_df_relative = self.output_df_relative.append(temp_series)

    def _create_ordered_output_dfs_from_series(self):
        """Put together the pandas series that hold sequences abundance outputs for each sample in order of the samples
        either according to a predefined ordered list or by an order that will be generated below."""
        if self.sorted_sample_uid_list:
            sys.stdout.write('\nValidating sorted sample list and ordering dataframe accordingly\n')
            self._check_sorted_sample_list_is_valid()

            self._create_ordered_output_dfs_from_series_with_sorted_sample_list()

        else:
            sys.stdout.write('\nGenerating ordered sample list and ordering dataframe accordingly\n')
            self.sorted_sample_uid_list = self._generate_ordered_sample_list()

            self._create_ordered_output_dfs_from_series_with_sorted_sample_list()

    def _generate_ordered_sample_list(self):
        """ Returns a list which is simply the ids of the samples ordered
        This will order the samples according to which sequence is their most abundant.
        I.e. samples found to have the sequence which is most abundant in the largest number of sequences
        will be first. Within each maj sequence, the samples will be sorted by the abundance of that sequence
        in the sample.
        At the moment we are also ordering by clade just so that you see samples with the A's at the top
        of the output so that we minimise the number of 0's in the top left of the output
        honestly I think we could perhaps get rid of this and just use the over all abundance of the sequences
        discounting clade. This is what we do for the clade order when plotting.
        """
        output_df_relative = self._make_raw_relative_abund_df_from_series()
        ordered_sample_list = self._get_sample_order_from_rel_seq_abund_df(output_df_relative)
        return ordered_sample_list

    def _make_raw_relative_abund_df_from_series(self):
        output_df_relative = pd.concat(
            [list_of_series[1] for list_of_series in self.dss_id_to_pandas_series_results_list_dict.values()],
            axis=1)
        output_df_relative = output_df_relative.T
        # now remove the rest of the non abundance columns
        non_seq_columns = [
            'sample_name', 'raw_contigs', 'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_qc_absolute_seqs',
            'post_qc_unique_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
            'post_taxa_id_absolute_symbiodinium_seqs',
            'post_taxa_id_unique_symbiodinium_seqs', 'post_med_absolute', 'post_med_unique',
            'size_screening_violation_absolute', 'size_screening_violation_unique']
        no_name_seq_columns = ['noName Clade {}'.format(clade) for clade in list('ABCDEFGHI')]
        cols_to_drop = non_seq_columns + no_name_seq_columns
        output_df_relative.drop(columns=cols_to_drop, inplace=True)
        return output_df_relative

    def _get_sample_order_from_rel_seq_abund_df(self, sequence_only_df_relative):

        max_seq_ddict, no_maj_samps, seq_to_samp_ddict = self._generate_most_abundant_sequence_dictionaries(
            sequence_only_df_relative)

        return self._generate_ordered_sample_list_from_most_abund_seq_dicts(max_seq_ddict, no_maj_samps,
                                                                            seq_to_samp_ddict)

    def _generate_ordered_sample_list_from_most_abund_seq_dicts(self, max_seq_ddict, no_maj_samps, seq_to_samp_ddict):
        # then once we have compelted this for all sequences go clade by clade
        # and generate the sample order
        ordered_sample_list_by_uid = []
        sys.stdout.write('\nGoing clade by clade sorting by abundance\n')
        for clade in list('ABCDEFGHI'):
            sys.stdout.write('\rGetting clade {} seqs'.format(clade))
            tup_list_of_clade = []
            # get the clade specific list of the max_seq_ddict
            for k, v in max_seq_ddict.items():
                sys.stdout.write('\r{}'.format(k))
                if k.startswith(clade) or k[-2:] == '_{}'.format(clade):
                    tup_list_of_clade.append((k, v))

            if not tup_list_of_clade:
                continue
            # now get an ordered list of the sequences for this clade
            sys.stdout.write('\rOrdering clade {} seqs'.format(clade))

            ordered_sequence_of_clade_list = [x[0] for x in sorted(tup_list_of_clade, key=lambda x: x[1], reverse=True)]

            for seq_to_order_samples_by in ordered_sequence_of_clade_list:
                sys.stdout.write('\r{}'.format(seq_to_order_samples_by))
                tup_list_of_samples_that_had_sequence_as_most_abund = seq_to_samp_ddict[seq_to_order_samples_by]
                ordered_list_of_samples_for_seq_ordered = \
                    [x[0] for x in
                     sorted(tup_list_of_samples_that_had_sequence_as_most_abund, key=lambda x: x[1], reverse=True)]
                ordered_sample_list_by_uid.extend(ordered_list_of_samples_for_seq_ordered)
        # finally add in the samples that didn't have a maj sequence
        ordered_sample_list_by_uid.extend(no_maj_samps)
        return ordered_sample_list_by_uid

    def _generate_most_abundant_sequence_dictionaries(self, sequence_only_df_relative):
        # {sequence_name_found_to_be_most_abund_in_sample: num_samples_it_was_found_to_be_most_abund_in}
        max_seq_ddict = defaultdict(int)
        # {most_abundant_seq_name: [(dss.id, rel_abund_of_most_abund_seq) for samples with that seq as most abund]}
        seq_to_samp_ddict = defaultdict(list)
        # a list to hold the names of samples in which there was no most abundant sequence identified
        no_maj_samps = []
        for sample_to_sort_uid in sequence_only_df_relative.index.values.tolist():
            sys.stdout.write(f'\r{sample_to_sort_uid}: Getting maj seq for sample')
            sample_series_as_float = self._get_sample_seq_abund_info_as_pd_series_float_type(
                sample_to_sort_uid, sequence_only_df_relative)
            max_abund_seq = self._get_name_of_most_abundant_seq(sample_series_as_float)
            max_rel_abund = self._get_rel_abund_of_most_abund_seq(sample_series_as_float)
            if not max_rel_abund > 0:
                no_maj_samps.append(sample_to_sort_uid)
            else:
                # add a tup of sample name and rel abund of seq to the seq_to_samp_dict
                seq_to_samp_ddict[max_abund_seq].append((sample_to_sort_uid, max_rel_abund))
                # add this to the ddict count
                max_seq_ddict[max_abund_seq] += 1
        return max_seq_ddict, no_maj_samps, seq_to_samp_ddict

    def _get_sample_seq_abund_info_as_pd_series_float_type(self, sample_to_sort_uid, sequence_only_df_relative):
        return sequence_only_df_relative.loc[sample_to_sort_uid].astype('float')

    def _get_rel_abund_of_most_abund_seq(self, sample_series_as_float):
        return sample_series_as_float.max()

    def _get_name_of_most_abundant_seq(self, sample_series_as_float):
        max_abund_seq = sample_series_as_float.idxmax()
        return max_abund_seq

    def _create_ordered_output_dfs_from_series_with_sorted_sample_list(self):
        # NB I was originally performing the concat directly on the managedSampleOutputDict (i.e. the mp dict)
        # but this was starting to produce errors. Starting to work on the dss_id_to_pandas_series_results_list_dict
        #  (i.e. normal, not mp, dict) seems to not produce these errors.
        sys.stdout.write('\rPopulating the absolute dataframe with series. This could take a while...')
        output_df_absolute = pd.concat(
            [list_of_series[0] for list_of_series in
             self.dss_id_to_pandas_series_results_list_dict.values()], axis=1)
        sys.stdout.write('\rPopulating the relative dataframe with series. This could take a while...')
        output_df_relative = pd.concat(
            [list_of_series[1] for list_of_series in
             self.dss_id_to_pandas_series_results_list_dict.values()], axis=1)
        # now transpose
        output_df_absolute = output_df_absolute.T
        output_df_relative = output_df_relative.T
        # now make sure that the order is correct.
        self.output_df_absolute = output_df_absolute.reindex(self.sorted_sample_uid_list)
        self.output_df_relative = output_df_relative.reindex(self.sorted_sample_uid_list)

    def _check_sorted_sample_list_is_valid(self):
        if len(self.sorted_sample_uid_list) != len(self.sample_list):
            raise RuntimeError({'message': 'Number of items in sorted_sample_list do not match those to be outputted!'})
        if self._smpls_in_sorted_smpl_list_not_in_list_of_samples():
            raise RuntimeError(
                {'message': 'Sample list passed in does not match sample list from db query'})

    def _smpls_in_sorted_smpl_list_not_in_list_of_samples(self):
        return list(set(self.sorted_sample_uid_list).difference(set([dss.id for dss in self.list_of_data_set_sample_objects])))

    def _generate_sample_output_series(self):
        """This generate a pandas series for each of the samples. It uses the ordered ReferenceSequence list created
         in the previous method as well as the other two dictionaries made.
         One df for absolute abundances and one for relative abundances. These series will be put together
         and ordered to construct the output data frames that will be written out for the user.
        """
        seq_count_table_output_series_generator_handler = SequenceCountTableOutputSeriesGeneratorHandler(
            clade_abundance_ordered_ref_seq_list=self.clade_abundance_ordered_ref_seq_list,
            dss_list=self.list_of_data_set_sample_objects,
            num_proc=self.num_proc)
        seq_count_table_output_series_generator_handler.execute_sequence_count_table_dataframe_contructor_handler(
            dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict=
            self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict,
            dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict=
            self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict)
        self.dss_id_to_pandas_series_results_list_dict = \
            dict(seq_count_table_output_series_generator_handler.dss_id_to_pandas_series_results_list_mp_dict)

    def _collect_abundances_for_creating_the_output(self):
        sequence_count_table_ordered_seqs_handler_instance = SequenceCountTableCollectAbundanceHandler(
            self.list_of_data_set_sample_objects,
            self.ref_seqs_in_datasets,
            self.num_proc,
            self.ordered_list_of_clades_found)
        sequence_count_table_ordered_seqs_handler_instance.execute_sequence_count_table_ordered_seqs_worker()
        # update the dictionaries that will be used in the second worker from the first worker
        self.update_dicts_for_the_second_worker_from_first_worker(sequence_count_table_ordered_seqs_handler_instance)

    def update_dicts_for_the_second_worker_from_first_worker(self, sequence_count_table_ordered_seqs_handler_instance):
        self.dss_id_to_list_of_dsss_objects_dict_mp_dict = \
            sequence_count_table_ordered_seqs_handler_instance.dss_id_to_list_of_dsss_objects_mp_dict

        self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict = \
            sequence_count_table_ordered_seqs_handler_instance.\
                dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict

        self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict = \
            sequence_count_table_ordered_seqs_handler_instance.\
                dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict

        self.clade_abundance_ordered_ref_seq_list = \
            sequence_count_table_ordered_seqs_handler_instance.clade_abundance_ordered_ref_seq_list


class SequenceCountTableCollectAbundanceHandler:
    """The purpose of this handler and the associated worker is to populate three dictionaries that will be used
    in making the count table output.
    1 - dict(ref_seq_name : cumulative relative abundance for each sequence across all samples)
    2 - sample_id : list(
                         dict(ref_seq_of_sample_name:absolute_abundance_of_dsss_in_sample),
                         dict(ref_seq_of_sample_name:relative_abundance_of_dsss_in_sample)
                         )
    3 - sample_id : list(
                         dict(clade:total_abund_of_no_name_seqs_of_clade_in_q_),
                         dict(clade:relative_abund_of_no_name_seqs_of_clade_in_q_)
                         )
    Abbreviations:
    ds = DataSet
    dss = DataSetSample
    dsss = DataSetSampleSequence
    ref_seq = ReferenceSeqeunce
    The end product of this method will be returned to the count table creator. The first dict will be used to create a
    list of the ReferenceSequence objects of this output ordered first by clade and then by cumulative relative
    abundance across all samples in the output.
    """
    def __init__(self, seq_count_table_creator_list_of_data_set_sample_objects, seq_count_table_creator_ref_seqs_in_datasets, seq_count_table_creator_num_procesors, seq_count_table_creator_ordered_list_of_clades_found):

        self.list_of_data_set_sample_objects = seq_count_table_creator_list_of_data_set_sample_objects
        self.mp_manager = Manager()
        self.input_dss_mp_queue = Queue()
        self._populate_input_dss_mp_queue()

        self.ref_seq_names_clade_annotated = [
        ref_seq.name if ref_seq.has_name else
        str(ref_seq.id) + '_{}'.format(ref_seq.clade) for
            ref_seq in seq_count_table_creator_ref_seqs_in_datasets]
        self.ordered_list_of_clades_found = seq_count_table_creator_ordered_list_of_clades_found
        self.num_proc = seq_count_table_creator_num_procesors
        #TODO we were previously creating an MP dictionary for every proc used. We were then collecting them afterwards
        # I'm not sure if there was a good reason for doing this, but I don't see any comments to the contrary.
        # it should not be necessary to have a dict for every proc. Instead we can just have on mp dict.
        # we should check that this is still working as expected.
        # self.list_of_dictionaries_for_processes = self._generate_list_of_dicts_for_processes()
        self.dss_id_to_list_of_dsss_objects_mp_dict = self.mp_manager.dict()
        self._populate_dss_id_to_list_of_dsss_objects()
        self.annotated_dss_name_to_cummulative_rel_abund_mp_dict = self.mp_manager.dict(
            {refSeq_name: 0 for refSeq_name in self.ref_seq_names_clade_annotated})
        self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict = self.mp_manager.dict(),
        self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict = self.mp_manager.dict()

        # this is the list that we will use the self.annotated_dss_name_to_cummulative_rel_abund_mp_dict to create
        # it is a list of the ref_seqs_ordered first by clade then by abundance.
        self.clade_abundance_ordered_ref_seq_list = []

    def execute_sequence_count_table_ordered_seqs_worker(self):
        all_processes = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        for n in range(self.num_proc):
            p = Process(target=self._sequence_count_table_ordered_seqs_worker, args=())
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

        self._generate_clade_abundance_ordered_ref_seq_list_from_seq_name_abund_dict()

    def _generate_clade_abundance_ordered_ref_seq_list_from_seq_name_abund_dict(self):
        for i in range(len(self.ordered_list_of_clades_found)):
            temp_within_clade_list_for_sorting = []
            for seq_name, abund_val in self.annotated_dss_name_to_cummulative_rel_abund_mp_dict.items():
                if seq_name.startswith(self.ordered_list_of_clades_found[i]) or seq_name[
                                                                                -2:] == f'_{self.ordered_list_of_clades_found[i]}':
                    # then this is a seq of the clade in Q and we should add to the temp list
                    temp_within_clade_list_for_sorting.append((seq_name, abund_val))
            # now sort the temp_within_clade_list_for_sorting and add to the cladeAbundanceOrderedRefSeqList
            sorted_within_clade = [
                a[0] for a in sorted(temp_within_clade_list_for_sorting, key=lambda x: x[1], reverse=True)]

            self.clade_abundance_ordered_ref_seq_list.extend(sorted_within_clade)

    def _sequence_count_table_ordered_seqs_worker(self):

        for dss in iter(self.input_dss_mp_queue.get, 'STOP'):
            sys.stdout.write(f'\r{dss.name}: collecting seq abundances')
            sequence_count_table_ordered_seqs_worker_instance = SequenceCountTableCollectAbundanceWorker(
                dss.id,
                self.annotated_dss_name_to_cummulative_rel_abund_mp_dict,
                self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict,
                self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict,
                self.ref_seq_names_clade_annotated, self.dss_id_to_list_of_dsss_objects_mp_dict)
            sequence_count_table_ordered_seqs_worker_instance.execute()



    def _populate_input_dss_mp_queue(self):
        for dss in self.list_of_data_set_sample_objects:
            self.input_dss_mp_queue.put(dss)

        for N in range(self.num_proc):
            self.input_dss_mp_queue.put('STOP')

    def _populate_dss_id_to_list_of_dsss_objects(self):
        for dss in self.list_of_data_set_sample_objects:
            sys.stdout.write('\r{}'.format(dss.name))
            self.dss_id_to_list_of_dsss_objects_mp_dict[dss.id] = list(
                DataSetSampleSequence.objects.filter(data_set_sample_from=dss))


class SequenceCountTableCollectAbundanceWorker:
    def __init__(self, dss,
                annotated_dss_name_to_cummulative_rel_abund_mp_dict,
                dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts,
                dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict,
                ref_seq_names_clade_annotated,
                 dss_id_to_list_of_dsss_objects_dict):

        self.dss = dss
        self.total_sequence_abundance_for_sample = sum([int(a) for a in json.loads(dss.cladal_seq_totals)])
        self.annotated_dss_name_to_cummulative_rel_abund_mp_dict = annotated_dss_name_to_cummulative_rel_abund_mp_dict
        self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts = dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts
        self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict = dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict
        self.ref_seq_names_clade_annotated = ref_seq_names_clade_annotated
        self.dss_id_to_list_of_dsss_objects_dict = dss_id_to_list_of_dsss_objects_dict
        cladal_abundances = [int(a) for a in json.loads(self.dss.cladal_seq_totals)]
        self.total_abundance_of_sequences_in_sample = sum(cladal_abundances)

    def execute(self):

        clade_summary_absolute_dict, clade_summary_relative_dict = \
            self._generate_empty_noname_seq_abund_summary_by_clade_dicts()

        smple_seq_count_aboslute_dict, smple_seq_count_relative_dict = self._generate_empty_seq_name_to_abund_dicts()

        dsss_in_sample = self.dss_id_to_list_of_dsss_objects_dict[self.dss.id]

        for dsss in dsss_in_sample:
            # determine what the name of the seq will be in the output
            name_unit = self._determine_output_name_of_dsss_and_pop_noname_clade_dicts(
                clade_summary_absolute_dict, clade_summary_relative_dict, dsss)

            self._populate_abs_and_rel_abundances_for_dsss(dsss, name_unit, smple_seq_count_aboslute_dict,
                                                           smple_seq_count_relative_dict)

        self._associate_sample_abundances_to_mp_dicts(
            clade_summary_absolute_dict,
            clade_summary_relative_dict,
            smple_seq_count_aboslute_dict,
            smple_seq_count_relative_dict)

    def _associate_sample_abundances_to_mp_dicts(self, clade_summary_absolute_dict, clade_summary_relative_dict,
                                                 smple_seq_count_aboslute_dict, smple_seq_count_relative_dict):
        self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts[self.dss.id] = [smple_seq_count_aboslute_dict,
                                                                                         smple_seq_count_relative_dict]
        self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict[self.dss.id] = [
            clade_summary_absolute_dict, clade_summary_relative_dict]

    def _populate_abs_and_rel_abundances_for_dsss(self, dsss, name_unit, smple_seq_count_aboslute_dict,
                                                  smple_seq_count_relative_dict):
        rel_abund_of_dsss = dsss.abundance / self.total_abundance_of_sequences_in_sample
        self.annotated_dss_name_to_cummulative_rel_abund_mp_dict[name_unit] += rel_abund_of_dsss
        smple_seq_count_aboslute_dict[name_unit] += dsss.abundance
        smple_seq_count_relative_dict[name_unit] += rel_abund_of_dsss

    def _determine_output_name_of_dsss_and_pop_noname_clade_dicts(
            self, clade_summary_absolute_dict, clade_summary_relative_dict, dsss):
        if not dsss.reference_sequence_of.has_name:
            name_unit = str(dsss.reference_sequence_of.id) + '_{}'.format(dsss.reference_sequence_of.clade)
            # the clade summries are only for the noName seqs
            clade_summary_absolute_dict[dsss.reference_sequence_of.clade] += dsss.abundance
            clade_summary_relative_dict[
                dsss.reference_sequence_of.clade] += dsss.abundance / self.total_abundance_of_sequences_in_sample
        else:
            name_unit = dsss.reference_sequence_of.name
        return name_unit

    def _generate_empty_seq_name_to_abund_dicts(self):
        smple_seq_count_aboslute_dict = {seq_name: 0 for seq_name in self.ref_seq_names_clade_annotated}
        smple_seq_count_relative_dict = {seq_name: 0 for seq_name in self.ref_seq_names_clade_annotated}
        return smple_seq_count_aboslute_dict, smple_seq_count_relative_dict

    def _generate_empty_noname_seq_abund_summary_by_clade_dicts(self):
        clade_summary_absolute_dict = {clade: 0 for clade in list('ABCDEFGHI')}
        clade_summary_relative_dict = {clade: 0 for clade in list('ABCDEFGHI')}
        return clade_summary_absolute_dict, clade_summary_relative_dict


class SequenceCountTableOutputSeriesGeneratorHandler:
    def __init__(self, clade_abundance_ordered_ref_seq_list, dss_list, num_proc):
        self.num_proc = num_proc
        self.dss_list = dss_list
        self.clade_abundance_ordered_ref_seq_list = clade_abundance_ordered_ref_seq_list
        self.output_df_header = self._create_output_df_header()
        self.worker_manager = Manager()
        # dss.id : [pandas_series_for_absolute_abundace, pandas_series_for_absolute_abundace]
        self.dss_id_to_pandas_series_results_list_mp_dict = self.worker_manager.dict()
        self.dss_input_queue = Queue()
        self._populate_dss_input_queue()


    def execute_sequence_count_table_dataframe_contructor_handler(
            self,
            dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict,
            dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict):
        all_processes = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        sys.stdout.write('\n\nOutputting seq data\n')
        for N in range(self.num_proc):
            p = Process(target=self._output_df_contructor_worker, args=(
                dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict,
                dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict
            ))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

    def _output_df_contructor_worker(self, sequence_count_table_creator_dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict, sequence_count_table_creator_dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict):

        for dss in iter(self.dss_input_queue.get, 'STOP'):
            seq_count_table_df_contructor_worker_instance = SequenceCountTableOutputSeriesGeneratorWorker(
                dss,
                self.dss_id_to_pandas_series_results_list_mp_dict,
                self.output_df_header,
                self.clade_abundance_ordered_ref_seq_list,
                sequence_count_table_creator_dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict,
                sequence_count_table_creator_dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict)

            seq_count_table_df_contructor_worker_instance.execute()


    def _populate_dss_input_queue(self):
        for dss in self.dss_list:
            self.dss_input_queue.put(dss)

        for N in range(self.num_proc):
            self.dss_input_queue.put('STOP')

    def _create_output_df_header(self):
        header_pre = self.clade_abundance_ordered_ref_seq_list
        no_name_summary_strings = ['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']]
        qc_stats = [
            'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs', 'post_taxa_id_absolute_symbiodinium_seqs',
            'post_taxa_id_unique_symbiodinium_seqs', 'size_screening_violation_absolute',
            'size_screening_violation_unique',
            'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
            'post_med_absolute',
            'post_med_unique']

        # append the noName sequences as individual sequence abundances
        return ['sample_name'] + qc_stats + no_name_summary_strings + header_pre


class SequenceCountTableOutputSeriesGeneratorWorker:
    def __init__(
            self, dss,
            dss_id_to_pandas_series_results_list_mp_dict, output_df_header, clade_abundance_ordered_ref_seq_list,
                sequence_count_table_creator_dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict,
                sequence_count_table_creator_dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict
    ):
        self.dss = dss
        self.dss_id_to_pandas_series_results_list_mp_dict = dss_id_to_pandas_series_results_list_mp_dict
        self.output_df_header = output_df_header
        self.clade_abundance_ordered_ref_seq_list = clade_abundance_ordered_ref_seq_list
        # dss.id : [{dsss:absolute abundance in dss}, {dsss:relative abundance in dss}]
        self.dss_abundance_mp_dict = \
            sequence_count_table_creator_dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict
        # dss.id : [{clade:total absolute abundance of no name seqs from that clade},
        #           {clade:total relative abundance of no name seqs from that clade}
        #          ]
        self.dss_noname_clade_summary_mp_dict = \
            sequence_count_table_creator_dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict
        self.sample_row_data_absolute = []
        self.sample_row_data_relative = []
        self.sample_seq_tot = sum([int(a) for a in json.loads(dss.cladal_seq_totals)])

    def execute(self):
        sys.stdout.write(f'\r{self.dss.name}: Creating data ouput row')
        if self._dss_had_problem_in_processing():
            self.sample_row_data_absolute.append(self.dss.name)
            self.sample_row_data_relative.append(self.dss.name)

            self._populate_quality_control_data_of_failed_sample()

            self._output_the_failed_sample_pandas_series()
            return

        self._populate_quality_control_data_of_successful_sample()

        self._output_the_successful_sample_pandas_series()

    def _output_the_successful_sample_pandas_series(self):
        sample_series_absolute = pd.Series(self.sample_row_data_absolute, index=self.output_df_header, name=self.dss.id)
        sample_series_relative = pd.Series(self.sample_row_data_relative, index=self.output_df_header, name=self.dss.id)
        self.dss_id_to_pandas_series_results_list_mp_dict[self.dss.id] = [sample_series_absolute, sample_series_relative]

    def _populate_quality_control_data_of_successful_sample(self):
        # Here we add in the post qc and post-taxa id counts
        # For the absolute counts we will report the absolute seq number
        # For the relative counts we will report these as proportions of the sampleSeqTot.
        # I.e. we will have numbers larger than 1 for many of the values and the symbiodinium seqs should be 1
        self.sample_row_data_absolute.append(self.dss.name)
        self.sample_row_data_relative.append(self.dss.name)

        # CONTIGS
        # This is the absolute number of sequences after make.contigs
        contig_num = self.dss.num_contigs
        self.sample_row_data_absolute.append(contig_num)
        self.sample_row_data_relative.append(contig_num / self.sample_seq_tot)
        # POST-QC
        # store the aboslute number of sequences after sequencing QC at this stage
        post_qc_absolute = self.dss.post_qc_absolute_num_seqs
        self.sample_row_data_absolute.append(post_qc_absolute)
        self.sample_row_data_relative.append(post_qc_absolute / self.sample_seq_tot)
        # This is the unique number of sequences after the sequencing QC
        post_qc_unique = self.dss.post_qc_unique_num_seqs
        self.sample_row_data_absolute.append(post_qc_unique)
        self.sample_row_data_relative.append(post_qc_unique / self.sample_seq_tot)
        # POST TAXA-ID
        # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
        tax_id_symbiodinium_absolute = self.dss.absolute_num_sym_seqs
        self.sample_row_data_absolute.append(tax_id_symbiodinium_absolute)
        self.sample_row_data_relative.append(tax_id_symbiodinium_absolute / self.sample_seq_tot)
        # Same as above but the number of unique seqs
        tax_id_symbiodinium_unique = self.dss.unique_num_sym_seqs
        self.sample_row_data_absolute.append(tax_id_symbiodinium_unique)
        self.sample_row_data_relative.append(tax_id_symbiodinium_unique / self.sample_seq_tot)
        # store the absolute number of sequences lost to size cutoff violations
        size_violation_aboslute = self.dss.size_violation_absolute
        self.sample_row_data_absolute.append(size_violation_aboslute)
        self.sample_row_data_relative.append(size_violation_aboslute / self.sample_seq_tot)
        # store the unique size cutoff violations
        size_violation_unique = self.dss.size_violation_unique
        self.sample_row_data_absolute.append(size_violation_unique)
        self.sample_row_data_relative.append(size_violation_unique / self.sample_seq_tot)
        # store the abosolute number of sequenes that were not considered Symbiodinium
        tax_id_non_symbiodinum_abosulte = self.dss.non_sym_absolute_num_seqs
        self.sample_row_data_absolute.append(tax_id_non_symbiodinum_abosulte)
        self.sample_row_data_relative.append(tax_id_non_symbiodinum_abosulte / self.sample_seq_tot)
        # This is the number of unique sequences that were not considered Symbiodinium
        tax_id_non_symbiodinium_unique = self.dss.non_sym_unique_num_seqs
        self.sample_row_data_absolute.append(tax_id_non_symbiodinium_unique)
        self.sample_row_data_relative.append(tax_id_non_symbiodinium_unique / self.sample_seq_tot)
        # Post MED absolute
        post_med_absolute = self.dss.post_med_absolute
        self.sample_row_data_absolute.append(post_med_absolute)
        self.sample_row_data_relative.append(post_med_absolute / self.sample_seq_tot)
        # Post MED unique
        post_med_unique = self.dss.post_med_unique
        self.sample_row_data_absolute.append(post_med_unique)
        self.sample_row_data_relative.append(post_med_unique / self.sample_seq_tot)

        # now add the clade divided summaries of the clades
        for clade in list('ABCDEFGHI'):
            self.sample_row_data_absolute.append(self.dss_noname_clade_summary_mp_dict[self.dss.id][0][clade])
            self.sample_row_data_relative.append(self.dss_noname_clade_summary_mp_dict[self.dss.id][1][clade])

        # and append these abundances in order of cladeAbundanceOrderedRefSeqList to
        # the sampleRowDataCounts and the sampleRowDataProps
        for seq_name in self.clade_abundance_ordered_ref_seq_list:
            sys.stdout.write('\rOutputting seq data for {}: sequence {}'.format(self.dss.name, seq_name))
            self.sample_row_data_absolute.append(self.dss_abundance_mp_dict[self.dss.id][0][seq_name])
            self.sample_row_data_relative.append(self.dss_abundance_mp_dict[self.dss.id][1][seq_name])

    def _output_the_failed_sample_pandas_series(self):
        sample_series_absolute = pd.Series(self.sample_row_data_absolute, index=self.output_df_header, name=self.dss.id)
        sample_series_relative = pd.Series(self.sample_row_data_relative, index=self.output_df_header, name=self.dss.id)
        self.dss_id_to_pandas_series_results_list_mp_dict[self.dss.id] = [sample_series_absolute,
                                                                          sample_series_relative]

    def _populate_quality_control_data_of_failed_sample(self):
        # Add in the qc totals if possible
        # For the proportions we will have to add zeros as we cannot do proportions
        # CONTIGS
        # This is the absolute number of sequences after make.contigs

        if self.dss.num_contigs:
            contig_num = self.dss.num_contigs
            self.sample_row_data_absolute.append(contig_num)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # POST-QC
        # store the aboslute number of sequences after sequencing QC at this stage
        if self.dss.post_qc_absolute_num_seqs:
            post_qc_absolute = self.dss.post_qc_absolute_num_seqs
            self.sample_row_data_absolute.append(post_qc_absolute)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # This is the unique number of sequences after the sequencing QC
        if self.dss.post_qc_unique_num_seqs:
            post_qc_unique = self.dss.post_qc_unique_num_seqs
            self.sample_row_data_absolute.append(post_qc_unique)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # POST TAXA-ID
        # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
        if self.dss.absolute_num_sym_seqs:
            tax_id_symbiodinium_absolute = self.dss.absolute_num_sym_seqs
            self.sample_row_data_absolute.append(tax_id_symbiodinium_absolute)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # Same as above but the number of unique seqs
        if self.dss.unique_num_sym_seqs:
            tax_id_symbiodinium_unique = self.dss.unique_num_sym_seqs
            self.sample_row_data_absolute.append(tax_id_symbiodinium_unique)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # size violation absolute
        if self.dss.size_violation_absolute:
            size_viol_ab = self.dss.size_violation_absolute
            self.sample_row_data_absolute.append(size_viol_ab)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # size violation unique
        if self.dss.size_violation_unique:
            size_viol_uni = self.dss.size_violation_unique
            self.sample_row_data_absolute.append(size_viol_uni)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # store the abosolute number of sequenes that were not considered Symbiodinium
        if self.dss.non_sym_absolute_num_seqs:
            tax_id_non_symbiodinum_abosulte = self.dss.non_sym_absolute_num_seqs
            self.sample_row_data_absolute.append(tax_id_non_symbiodinum_abosulte)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # This is the number of unique sequences that were not considered Symbiodinium
        if self.dss.non_sym_unique_num_seqs:
            tax_id_non_symbiodinium_unique = self.dss.non_sym_unique_num_seqs
            self.sample_row_data_absolute.append(tax_id_non_symbiodinium_unique)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # post-med absolute
        if self.dss.post_med_absolute:
            post_med_abs = self.dss.post_med_absolute
            self.sample_row_data_absolute.append(post_med_abs)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # post-med absolute
        if self.dss.post_med_unique:
            post_med_uni = self.dss.post_med_unique
            self.sample_row_data_absolute.append(post_med_uni)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)

        # no name clade summaries get 0.
        for _ in list('ABCDEFGHI'):
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)

        # All sequences get 0s
        for _ in self.clade_abundance_ordered_ref_seq_list:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)

        # All sequences get 0s
        for _ in self.clade_abundance_ordered_ref_seq_list:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)

    def _dss_had_problem_in_processing(self):
        return self.dss.error_in_processing or self.sample_seq_tot == 0


class SubPlotter:
    """A class that can be used for the sub plots of the SeqStackedBarPlotter and the TypeStackedBarPlotter.
    """
    def __init__(self, parent_plotter_instance, index_of_this_subplot):
        self.parent_plotter = parent_plotter_instance
        self.patches_list = []
        self.x_tick_label_list = []
        self.x_index_for_plot = 0
        self.colour_list = []
        self.index_of_this_subplot = index_of_this_subplot
        self.subplot_axes = self.parent_plotter.axarr[self.index_of_this_subplot]
        self.end_slice = self._get_end_index_for_slicing_plotting_data()
        self.num_samples_in_this_subplot = len(
                self.parent_plotter.output_count_table_as_df.index.values.tolist()[
                self.index_of_this_subplot * self.parent_plotter.samples_per_subplot:self.end_slice])
        # Make a custom colour map
        # https://matplotlib.org/api/_as_gen/matplotlib.colors.ListedColormap.html
        self.listed_colour_map = None
        self.patches_collection = None


    def plot_seq_subplot(self):
        self._create_rect_patches_and_populate_colour_list()

        self._make_listed_colour_map()

        self._make_patches_collection()

        self._draw_patches_on_axes()

        self._format_axes()



    def _format_axes(self):
        # make it so that the x axes is constant length that will be the num of samples per subplot
        self.subplot_axes.set_xlim(0 - 0.5, self.parent_plotter.samples_per_subplot - 0.5)
        self.subplot_axes.set_ylim(0, 1)
        self.subplot_axes.set_xticks(range(self.num_samples_in_this_subplot))
        self.subplot_axes.set_xticklabels(self.x_tick_label_list, rotation='vertical', fontsize=6)
        self.subplot_axes.spines['right'].set_visible(False)
        self.subplot_axes.spines['top'].set_visible(False)
        # as well as getting rid of the top and right axis splines
        # I'd also like to restrict the bottom spine to where there are samples plotted but also
        # maintain the width of the samples
        # I think the easiest way to do this is to hack a bit by setting the x axis spines to invisible
        # and then drawing on a line at y = 0 between the smallest and largest ind (+- 0.5)
        self.subplot_axes.spines['bottom'].set_visible(False)
        self.subplot_axes.add_line(
            Line2D((0 - 0.5, self.num_samples_in_this_subplot - 0.5), (0, 0), linewidth=2, color='black'))

    def _draw_patches_on_axes(self):
        self.subplot_axes.add_collection(self.patches_collection)
        self.subplot_axes.autoscale_view()
        self.subplot_axes.figure.canvas.draw()

    def _make_patches_collection(self):
        patches_collection = PatchCollection(self.patches_list, cmap=self.listed_colour_map)
        patches_collection.set_array(np.arange(len(self.patches_list)))

    def _make_listed_colour_map(self):
        self.listed_colour_map = ListedColormap(self.colour_list)

    def _get_end_index_for_slicing_plotting_data(self):
        if self.index_of_this_subplot == self.parent_plotter.number_of_subplots - 1:
            end_slice = self.parent_plotter.num_samples
        else:
            end_slice = self.parent_plotter.samples_per_subplot * (self.index_of_this_subplot + 1)
        return end_slice

    def _create_rect_patches_and_populate_colour_list(self):
        for sample in self.parent_plotter.output_count_table_as_df.index.values.tolist()[
                      self.index_of_this_subplot * self.parent_plotter.samples_per_subplot:self.end_slice]:
            sys.stdout.write(f'\rPlotting sample: {sample}')
            self._add_sample_names_to_tick_label_list(sample)
            # for each sample we will start at 0 for the y and then add the height of each bar to this
            bottom = 0
            # for each sequence, create a rect patch
            # the rect will be 1 in width and centered about the ind value.
            for seq in list(self.parent_plotter.output_count_table_as_df):
                # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
                rel_abund = self.parent_plotter.output_count_table_as_df.loc[sample, seq]
                if rel_abund > 0:
                    self.patches_list.append(Rectangle(
                        (self.x_index_for_plot - 0.5, bottom),
                        1,
                        rel_abund, color=self.parent_plotter.colour_dict[seq]))
                    self.colour_list.append(self.parent_plotter.colour_dict[seq])
                    bottom += rel_abund
            self.x_index_for_plot += 1

    def _add_sample_names_to_tick_label_list(self, sample):
        sample_name = self.parent_plotter.smpl_id_to_smp_name_dict[int(sample)]
        if len(sample_name) < 20:
            self.x_tick_label_list.append(self.parent_plotter.smpl_id_to_smp_name_dict[int(sample)])
        else:
            self.x_tick_label_list.append(f'uid_{int(sample)}')


class LegendPlotter:
    """This class can be used by the SeqStackedBarPlotter and the TypeStackedBarPlotter to handle
    the plotting of the legend subplot.
    """
    def __init__(self, parent_plotter):
        self.parent_plotter = parent_plotter
        self.ax_to_plot_on = self.parent_plotter.axarr[-1]
        # legend setup parameters
        self.y_coord_increments = 100 / self.parent_plotter.max_n_rows
        self.leg_box_depth = 2 / 3 * self.y_coord_increments

        self.x_coord_increments = 100 / self.parent_plotter.max_n_cols
        self.leg_box_width = self.x_coord_increments / 3
        self._set_n_rows_and_last_row_len()
        self.sequence_count = 0

    def plot_legend_seqs(self):
        self._set_ylim_and_x_lim_and_invert_y_axis()

        self._plot_legend_rows()

        self._remove_frames_from_axis()

    def _plot_legend_rows(self):
        sys.stdout.write(
            f'\nGenerating figure legend for {str(self.parent_plotter.num_leg_cells)} most common sequences\n')

        for row_increment in range(min(self.n_rows, self.parent_plotter.max_n_rows)):

            if self._this_is_last_row_of_legend(row_increment=row_increment):
                for col_increment in range(self.parent_plotter.max_n_cols):
                    self._plot_legend_row(row_increment=row_increment, col_increment=col_increment)

                    self.sequence_count += 1
            else:
                for col_increment in range(self.last_row_len):
                    self._plot_legend_row(row_increment=row_increment, col_increment=col_increment)

                    self.sequence_count += 1

    def _set_ylim_and_x_lim_and_invert_y_axis(self):
        # Once we know the number of rows, we can also adjust the y axis limits
        self.ax_to_plot_on.set_xlim(0, 100)
        self.ax_to_plot_on.set_ylim(0, ((self.n_rows - 1) * self.y_coord_increments) + self.leg_box_depth)
        self.ax_to_plot_on.invert_yaxis()

    def _set_n_rows_and_last_row_len(self):
        if len(self.parent_plotter.ordered_list_of_seqs_names) < self.parent_plotter.num_leg_cells:
            if len(self.parent_plotter.ordered_list_of_seqs_names) % self.parent_plotter.max_n_cols != 0:
                self.n_rows = int(
                    len(self.parent_plotter.ordered_list_of_seqs_names) / self.parent_plotter.max_n_cols) + 1
                self.last_row_len = len(self.parent_plotter.ordered_list_of_seqs_names) % self.parent_plotter.max_n_cols
            else:
                self.n_rows = int(len(self.parent_plotter.ordered_list_of_seqs_names) / self.parent_plotter.max_n_cols)
                self.last_row_len = self.parent_plotter.max_n_cols
        else:
            self.n_rows = self.parent_plotter.max_n_rows
            self.last_row_len = self.parent_plotter.max_n_cols


    def _this_is_last_row_of_legend(self, row_increment):
        return (row_increment + 1) != self.n_rows

    def _remove_frames_from_axis(self):
        self.ax_to_plot_on.set_frame_on(False)
        self.ax_to_plot_on.get_xaxis().set_visible(False)
        self.ax_to_plot_on.get_yaxis().set_visible(False)

    def _plot_legend_row(self, row_increment, col_increment):
        leg_box_x, leg_box_y = self._add_legend_rect(col_increment=col_increment, row_increment=row_increment)
        self._add_legend_text(leg_box_x, leg_box_y)

    def _add_legend_text(self, leg_box_x, leg_box_y):
        text_x = leg_box_x + self.leg_box_width + (0.2 * self.leg_box_width)
        text_y = leg_box_y + (0.5 * self.leg_box_depth)
        self.ax_to_plot_on.text(
            text_x, text_y, self.parent_plotter.ordered_list_of_seqs_names[self.sequence_count],
            verticalalignment='center', fontsize=8)

    def _add_legend_rect(self, col_increment, row_increment):
        leg_box_x = col_increment * self.x_coord_increments
        leg_box_y = row_increment * self.y_coord_increments
        self.ax_to_plot_on.add_patch(Rectangle((leg_box_x, leg_box_y),
                                           width=self.leg_box_width, height=self.leg_box_depth,
                                           color=self.parent_plotter.colour_dict[self.parent_plotter.ordered_list_of_seqs_names[self.sequence_count]]))
        return leg_box_x, leg_box_y


class SeqStackedBarPlotter:
    def __init__(self, seq_relative_abund_count_table_path, output_directory, time_date_str=None, ordered_sample_uid_list=None):
        self.seq_relative_abund_count_table_path = seq_relative_abund_count_table_path
        self.output_directory = output_directory
        if time_date_str:
            self.time_date_str = time_date_str
        else:
            self.time_date_str = str(datetime.now()).replace(' ', '_').replace(':', '-')
        self.fig_output_base = os.path.join(self.output_directory, f'{self.time_date_str}')
        self.smpl_id_to_smp_name_dict = None
        self.output_count_table_as_df = self._create_output_df_and_populate_smpl_id_to_smp_name_dict()
        self.ordered_list_of_seqs_names = self._set_ordered_list_of_seqs_names()
        # legend parameters and vars
        self.max_n_cols = 8
        self.max_n_rows = 7
        self.num_leg_cells = self.max_n_rows * self.max_n_cols
        self.colour_dict = self._set_colour_dict()
        # plotting vars
        self.ordered_sample_uid_list = self._set_ordered_sample_uid_list_and_reorder_df(ordered_sample_uid_list)
        self.num_samples = len(self.output_count_table_as_df.index.values.tolist())
        self.samples_per_subplot = 50
        self.number_of_subplots = self._infer_number_of_subplots()
        # we add  1 to the n_subplots here for the legend at the bottom
        self.f, self.axarr = plt.subplots(self.number_of_subplots + 1, 1, figsize=(10, 3 * self.number_of_subplots))
        self.output_path_list = []


    def plot_stacked_bar_seqs(self):
        for sub_plot_index in range(self.number_of_subplots):
            sub_plotter = SubPlotter(index_of_this_subplot=sub_plot_index, parent_plotter_instance=self)
            sub_plotter.plot_seq_subplot()

        self._plot_legend()

        plt.tight_layout()

        self._write_out_plot()

        self.output_path_list.extend(
            [
                f'{self.fig_output_base}_seq_abundance_stacked_bar_plot.svg',
                f'{self.fig_output_base}_seq_abundance_stacked_bar_plot.png'
            ])

    def _write_out_plot(self):
        sys.stdout.write('\nsaving as .svg\n')
        plt.savefig('{}_seq_abundance_stacked_bar_plot.svg'.format(self.fig_output_base))
        sys.stdout.write('\nsaving as .png\n')
        plt.savefig('{}_seq_abundance_stacked_bar_plot.png'.format(self.fig_output_base))

    def _plot_legend(self):
        legend_plotter = LegendPlotter(parent_plotter=self)
        legend_plotter.plot_legend_seqs()

    def _get_end_index_for_slicing_plotting_data(self, i):
        if i == self.number_of_subplots - 1:
            end_slice = self.num_samples
        else:
            end_slice = self.samples_per_subplot * (i + 1)
        return end_slice

    def _add_sample_names_to_tick_label_list(self, sample, x_tick_label_list):
        sample_name = self.smpl_id_to_smp_name_dict[int(sample)]
        if len(sample_name) < 20:
            x_tick_label_list.append(self.smpl_id_to_smp_name_dict[int(sample)])
        else:
            x_tick_label_list.append(f'uid_{int(sample)}')

    def _infer_number_of_subplots(self):
        if (self.num_samples % self.samples_per_subplot) != 0:
            number_of_subplots = int(self.num_samples / self.samples_per_subplot) + 1
        else:
            number_of_subplots = int(self.num_samples / self.samples_per_subplot)
        return number_of_subplots

    def _set_ordered_sample_uid_list_and_reorder_df(self, ordered_sample_uid_list):
        """If we are plotting this in companion with an ITS2 type profile output then we will be passed a
        ordered_sample_uid_list. It is very useful to have the ITS2 type profile output figure and the seq figure
        in the same sample order for direct comparison.
        If this output is not associated with an ITS2 type profile output then we will need to
        generate the sample order from scratch.
        In theory the output should already be somewhat ordered in that the samples should be in order of similarity.
        However, these have the artifical clade ordering (sorted by clade and then by abundance of seqs) so for the
        plotting it will probably be better to get a new
        order for the samples that is not constrained to the order of the clades. For this we should order as usual
        according to the most common majority sequences and then within this grouping we should order according to the
        the abundance of these sequences within the samples.
        """
        if not ordered_sample_uid_list:


            self.ordered_sample_list = self._generate_sample_order_de_novo()
        else:
            self.ordered_sample_list = ordered_sample_uid_list

        self._reorder_df_by_new_sample_and_seq_order()

    def _reorder_df_by_new_sample_and_seq_order(self):
        self.output_count_table_as_df = self.output_count_table_as_df[self.ordered_list_of_seqs_names]
        self.output_count_table_as_df = self.output_count_table_as_df.reindex(
            [int(a) for a in self.ordered_sample_list])

    def _generate_sample_order_de_novo(self):
        """At this stage we have the ordered list of seqs we now need to order the samples
        this method will return us the names of the samples in order that they should be plotted.
        """
        # {sequence_name_found_to_be_most_abund_in_sample: num_samples_it_was_found_to_be_most_abund_in}
        max_seq_ddict = defaultdict(int)
        # {most_abundant_seq_name: [(dss.id, rel_abund_of_most_abund_seq) for samples with that seq as most abund]}
        seq_to_samp_ddict = defaultdict(list)

        # for each sample get the columns name of the max value of a div not including the columns in the following:
        no_maj_seq = []
        for sample_id_to_sort in self.output_count_table_as_df.index.values.tolist():
            smp_series = self.output_count_table_as_df.loc[sample_id_to_sort].astype('float')
            max_abund_seq_name = smp_series.idxmax()
            rel_abund_of_max_abund_seq_name = smp_series.max()
            if not rel_abund_of_max_abund_seq_name > 0:
                no_maj_seq.append(sample_id_to_sort)
            else:
                # add a tup of sample name and rel abund of seq to the seq_to_samp_dict
                seq_to_samp_ddict[max_abund_seq_name].append((sample_id_to_sort, rel_abund_of_max_abund_seq_name))
                # add this to the ddict count
                max_seq_ddict[max_abund_seq_name] += 1

        # then once we have compelted this for all sequences
        # generate the sample order according to the sequence order
        ordered_sample_list = []
        # get an ordered list of the sequencs according to the max_seq_ddict
        ordered_list_of_sequences = [x[0] for x in sorted(max_seq_ddict.items(), key=lambda x: x[1], reverse=True)]

        for seq_to_order_samples_by in ordered_list_of_sequences:
            tup_list_of_samples_that_had_sequence_as_most_abund = seq_to_samp_ddict[seq_to_order_samples_by]
            ordered_list_of_samples_for_seq_ordered = [
                x[0] for x in sorted(tup_list_of_samples_that_had_sequence_as_most_abund,
                                     key=lambda x: x[1], reverse=True)]
            ordered_sample_list.extend(ordered_list_of_samples_for_seq_ordered)

        ordered_sample_list.extend(no_maj_seq)

        return ordered_sample_list

    def _set_colour_dict(self):
        """Create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
        to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
        If we are only going to have a legend that is cols x rows as shown below, then we should only use
        that many colours in the plotting."""
        colour_palette, grey_palette = self._get_colour_lists()

        temp_colour_dict = {}
        for i in range(len(self.ordered_list_of_seqs_names)):
            if i < self.num_leg_cells:
                temp_colour_dict[self.ordered_list_of_seqs_names[i]] = colour_palette[i]
            else:
                grey_index = i % len(grey_palette)
                temp_colour_dict[self.ordered_list_of_seqs_names[i]] = grey_palette[grey_index]
        return temp_colour_dict

    def _get_colour_lists(self):
        colour_palette = self._get_colour_list()
        grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
        return colour_palette, grey_palette

    def _set_ordered_list_of_seqs_names(self):
        """Get a list of the sequences in order of their abundance and use this list to create the colour dict
        The abundances are got by simply summing up the columns
        """
        abundance_dict = {}
        for col in list(self.output_count_table_as_df):
            abundance_dict[col] = sum(self.output_count_table_as_df[col])

        # get the names of the sequences sorted according to their totalled abundance
        self.ordered_list_of_seqs_names = [x[0] for x in sorted(abundance_dict.items(), key=lambda x: x[1], reverse=True)]

    def _create_output_df_and_populate_smpl_id_to_smp_name_dict(self):
        """Drop the QC columns from the SP output df and also drop the clade summation columns
        we will be left with just columns for each one of the sequences found in the samples
        we need to drop the rows first before we can make the smp_id_to_smp_name_dict else
        we will have the final row names in the index which are not convertable to int
        need to make the smp_id_to_smp_name_dict before dropping the sample_name col"""

        sp_output_df = pd.read_csv(self.seq_relative_abund_count_table_path, sep='\t', lineterminator='\n', header=0, index_col=0)

        meta_index_to_cut_from = self._get_df_index_to_drop_from(sp_output_df)

        self._drop_meta_info_rows_from_df(meta_index_to_cut_from, sp_output_df)

        self._populate_smpl_id_to_smp_name_dict(sp_output_df)

        self._drop_non_seq_abund_cols_and_set_df_types(sp_output_df)

    def _drop_non_seq_abund_cols_and_set_df_types(self, sp_output_df):
        sp_output_df.drop(
            columns=['sample_name', 'noName Clade A', 'noName Clade B', 'noName Clade C', 'noName Clade D',
                     'noName Clade E', 'noName Clade F', 'noName Clade G', 'noName Clade H',
                     'noName Clade I', 'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                     'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                     'post_taxa_id_absolute_non_symbiodinium_seqs',
                     'post_taxa_id_unique_non_symbiodinium_seqs',
                     'size_screening_violation_absolute', 'size_screening_violation_unique',
                     'post_med_absolute', 'post_med_unique'
                     ], inplace=True)
        sp_output_df = sp_output_df.astype('float')
        sp_output_df.index = sp_output_df.index.astype('int')

    def _populate_smpl_id_to_smp_name_dict(self, sp_output_df):
        self.smpl_id_to_smp_name_dict = {
            int(uid): smp_name for uid, smp_name in
            zip(sp_output_df.index.values.tolist(), sp_output_df['sample_name'].values.tolist())}

    def _drop_meta_info_rows_from_df(self, meta_index_to_cut_from, sp_output_df):
        sp_output_df.drop(index=sp_output_df.index[range(meta_index_to_cut_from, 0, 1)], inplace=True)

    def _get_df_index_to_drop_from(self, sp_output_df):
        # In order to be able to drop the DIV row at the end and the meta information rows, we should
        # drop all rows that are after the DIV column. We will pass in an index value to the .drop
        # that is called here. To do this we need to work out which index we are working with
        meta_index_to_cut_from = None
        index_values_as_list = sp_output_df.index.values.tolist()
        for i in range(-1, -(len(index_values_as_list)), -1):
            if index_values_as_list[i].startswith('seq'):
                # then this is the index (in negative notation) that we need to cut from
                meta_index_to_cut_from = i
                break
        return meta_index_to_cut_from



    @staticmethod
    def _get_colour_list():
        colour_list = [
            "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#FFDBE5", "#7A4900",
            "#0000A6",
            "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400",
            "#4FC601",
            "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
            "#B903AA",
            "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
            "#372101",
            "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
            "#0CBD66",
            "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
            "#BEC459",
            "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9",
            "#FF913F",
            "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
            "#7900D7",
            "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
            "#201625",
            "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
            "#CB7E98",
            "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489",
            "#806C66",
            "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
            "#C895C5",
            "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58", "#7A7BFF", "#D68E01",
            "#353339",
            "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A",
            "#001325",
            "#02525F", "#0AA3F7", "#E98176", "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75",
            "#8D8546",
            "#9695C5", "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
            "#00005F",
            "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058",
            "#A45B02",
            "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
            "#F4D749",
            "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE",
            "#C6DC99",
            "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
            "#C6005A",
            "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183", "#AA5199", "#B5D6C3",
            "#A38469",
            "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433", "#789EC9", "#6D80BA", "#953F00", "#5EFF03",
            "#E4FFFC",
            "#1BE177", "#BCB1E5", "#76912F", "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213",
            "#A76F42",
            "#89412E", "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
            "#BDC9D2",
            "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71",
            "#868E7E",
            "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C",
            "#00B57F",
            "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
        return colour_list


class UnifracSeqAbundanceMPCollection:
    """The purpose of this class is to be used during the UnifracDistanceCreator as a tidy holder for the
     sequences information that is collected from each CladeCollection that is part of the output."""
    def __init__(self, clade_collection, proc_id):
        self.fasta_dict = {}
        self.name_dict = {}
        self.group_list = []
        self.ref_seq_id_list = []
        self.clade_collection = clade_collection
        self.proc_id = proc_id

    def collect_seq_info(self):
        sys.stdout.write('\rProcessing cc: {} with {}'.format(self.clade_collection, self.proc_id))
        for data_set_sample_seq in DataSetSampleSequence.objects.filter(
                clade_collection_found_in=self.clade_collection):
            ref_seq_id = data_set_sample_seq.reference_sequence_of.id
            self.ref_seq_id_list.append(ref_seq_id)
            unique_seq_name_base = '{}_id{}'.format(ref_seq_id, self.clade_collection.id)

            smp_name = str(self.clade_collection.id)
            self.fasta_dict['{}_{}'.format(unique_seq_name_base, 0)] = data_set_sample_seq.reference_sequence_of.sequence
            temp_name_list = []

            for i in range(data_set_sample_seq.abundance):
                temp_name_list.append('{}_{}'.format(unique_seq_name_base, i))
                self.group_list.append('{}\t{}'.format('{}_{}'.format(unique_seq_name_base, i), smp_name))

            self.name_dict['{}_{}'.format(unique_seq_name_base, 0)] = temp_name_list


class SequenceCollectionComplete(Exception): pass


class UnifracDistanceCreatorHandlerOne:
    """The purpose of this handler will be to setup the sequence collection for a given clade"""
    def __init__(self, parent_unifrac_dist_creator, clade_in_question):
        self.parent_unifrac_dist_creator = parent_unifrac_dist_creator
        self.clade = clade_in_question
        self.output_dir = os.path.join(self.parent_unifrac_dist_creator.output_dir, self.clade)
        self.clade_collections_of_clade = self.parent_unifrac_dist_creator.clade_collections_from_data_sets.filter(
            clade=self.clade)
        self._raise_runtime_error_if_not_enough_clade_collections()
        self.input_clade_collection_queue = Queue()
        self.output_unifrac_seq_abund_mp_collection_queue = Queue()
        self._populate_input_queue()
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

    def _unifrac_distance_creator_worker_one(self):
        proc_id = current_process().name
        for cc in iter(self.input_clade_collection_queue.get, 'STOP'):
            unifrac_seq_abundance_mp_collection = UnifracSeqAbundanceMPCollection(clade_collection=cc, proc_id=proc_id)
            unifrac_seq_abundance_mp_collection.collect_seq_info()
            self.output_unifrac_seq_abund_mp_collection_queue.put(unifrac_seq_abundance_mp_collection)
        self.output_unifrac_seq_abund_mp_collection_queue.put('EXIT')

    def execute_unifrac_distance_creator_worker_one(self):
        all_processes = self._start_sequence_collection_running()

        self._populate_the_master_seq_collection_objects()

        self._once_complete_wait_for_processes_to_complete(all_processes)

        self._write_out_master_fasta_names_and_group_files()

    def _write_out_master_fasta_names_and_group_files(self):
        write_list_to_destination(self.master_fasta_file_path, self.master_fasta_file_as_list)
        write_list_to_destination(self.master_names_file_path, self.master_names_file_as_list)
        write_list_to_destination(self.master_group_file_path, self.master_group_file_as_list)

    def _once_complete_wait_for_processes_to_complete(self, all_processes):
        # process the outputs of the sub processess before we pause to wait for them to complete.
        for p in all_processes:
            p.join()

    def _populate_the_master_seq_collection_objects(self):
        """This method collects the output of the sequence collection mp methods and populates the objects that
        represnt the master names, fasta and group files that we are creating to calculate the UniFrac"""
        for seq_collection_object in iter(self.output_unifrac_seq_abund_mp_collection_queue.get, 'STOP'):
            try:
                self._check_if_all_seq_collection_objects_have_been_processed(seq_collection_object)

                sys.stdout.write(
                    f'\rAdding {seq_collection_object.clade_collection.name} to master fasta and name files')

                self._update_master_group_list(seq_collection_object)

                for seq_id in seq_collection_object.ref_seq_id_list:
                    seq_name = f'{seq_id}_id{seq_collection_object.clade_collection.id}_0'
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
            self.output_unifrac_seq_abund_mp_collection_queue.put(cc)

        for n in range(self.parent_unifrac_dist_creator.num_proc):
            self.output_unifrac_seq_abund_mp_collection_queue.put('STOP')


    def _raise_runtime_error_if_not_enough_clade_collections(self):
        if len(self.clade_collections_of_clade) < 2:
            raise RuntimeWarning({'message' : 'insufficient clade collections'})


class FseqbootAlignmentGenerator:
    """This class is to be used in the UnifracDistanceCreator to handle the creation of the fseqboot alignments
    for a given clade. It produces a directory (self.output_dir) that contains a number of alignment replicates."""
    def __init__(self, parent_unifrac_dist_creator, clade_in_question, num_reps):
        self.parent_unifrac_dist_creator = parent_unifrac_dist_creator
        self.clade = clade_in_question
        self.num_reps = num_reps
        self.input_fasta_path = self.parent_unifrac_dist_creator.master_fasta_file_aligned_path
        self.output_seqboot_file_path = self.input_fasta_path + '.fseqboot'
        self.fseqboot_local = None
        self._setup_fseqboot_plum_bum_local()
        self.output_dir = os.path.join(self.parent_unifrac_dist_creator.output_dir, self.clade)
        # this is the base of the path that will be used to build the path of each of the alignment replicates
        # to build the complete path a rep_count string (str(count)) will be appended to this base.
        self.fseqboot_clade_root_dir = os.path.join(self.output_dir, 'out_seq_boot_reps')
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
        write_list_to_destination(f'{fseqboot_base}{rep_count}', self.fseqboot_individual_alignment_as_list)

    def _execute_fseqboot(self):
        sys.stdout.write('\rGenerating multiple fseqboot alignments')
        (self.fseqboot_local[
            '-sequence', self.input_fasta_path, '-outfile', self.output_seqboot_file_path,
            '-test', 'b', '-reps', self.num_reps])()

    def _setup_fseqboot_plum_bum_local(self):
        is_installed = subprocess.call(['which', 'fseqboot'])
        if is_installed == 0:
            self.fseqboot_local = local["fseqboot"]
        else:
            fseqboot_path = os.path.join(self.parent_unifrac_dist_creator.symportal_root_dir, 'lib', 'phylipnew', 'fseqboot')
            if os.path.isfile(fseqboot_path):
                self.fseqboot_local = local[fseqboot_path]
            else:
                raise RuntimeError('Cannot find fseqboot in PATH or in local installation at ./lib/phylipnew/fseqboot\n'
                         'For instructions on installing the phylipnew dependencies please visit the SymPortal'
                         'GitHub page: https://github.com/didillysquat/SymPortal_framework/'
                         'wiki/SymPortal-setup#6-third-party-dependencies')

class UnifracSubcladeHandler:
    """The output of this is a list that contains a list of paths to trees, one for each replicate fasta"""
    def __init__(self, parent_unifrac_creator, clade_in_question, fseqbootbase):
        self.clade = clade_in_question
        self.parent_unifrac_creator = parent_unifrac_creator
        self.fseqbootbase = fseqbootbase
        self.input_queue_of_rep_numbers = Queue()
        self.output_queue_of_paths_to_trees = Queue()
        self.output_dir = os.path.join(self.parent_unifrac_creator.output_dir, self.clade)
        self.list_of_output_tree_paths = None
        self.concat_tree_file_path = os.path.join(self.output_dir, 'concatenated_tree_file')
        self.consensus_tree_file_path = os.path.join(self.output_dir, 'consensus_tree_sumtrees.newick')
        self.unifrac_dist_file_path = None

    def perform_unifrac(self):
        mothur_analysis = MothurAnalysis.init_for_weighted_unifrac(
            tree_path=self.consensus_tree_file_path,
            group_file_path=self.parent_unifrac_creator.clade_master_group_file_path,
            name_file_path=self.parent_unifrac_creator.clade_master_name_file_path)
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
            unifrac_mothur_worker = UnifracMothurWorker(rep_num=rep_num, fseqbootbase = self.fseqbootbase)
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
        name_file = read_defined_file_to_list(self.parent_unifrac_creator.master_names_file_path)
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
                if 'recursion' in completed_consensus.stdout.decode('utf-8'):
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




class UnifracDistPCoACreator:
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
    def __init__(self, symportal_root_directory, data_set_string, num_processors, method, call_type, date_time_string, bootstrap_value=100, output_dir=None):
        self.symportal_root_dir = symportal_root_directory
        self.call_type = call_type
        self.data_sets_to_output = DataSet.objects.filter(id__in=[int(a) for a in str(data_set_string).split(',')])
        self.num_proc = num_processors
        self.method = method
        if date_time_string:
            self.date_time_string = date_time_string
        else:
            self.date_time_string = str(datetime.now()).replace(' ', '_').replace(':', '-')
        self.bootstrap_value = bootstrap_value
        self.output_file_paths = []
        self.clade_collections_from_data_sets = CladeCollection.objects.filter(
            data_set_sample_from__data_submission_from__in=self.data_sets_to_output)
        self.clades_of_clade_collections_list = list(set([a.clade for a in self.clade_collections_from_data_sets]))
        self.output_dir = self._setup_output_dir(call_type, data_set_string, output_dir)
        self.output_path_list = []
        # Clade based files
        self.clade_output_dir = None
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
        # the unifrac distances
        self.clade_unifrac_dist_file_path = None
        self.clade_pcoa_coord_file_path = None

    def compute_unifrac_dists_and_pcoa_coords(self):
        for clade_in_question in self.clade_collections_from_data_sets:
            self.clade_output_dir = os.path.join(self.output_dir, clade_in_question)

            try:
                self._create_and_write_out_master_fasta_names_and_group_files(clade_in_question)
            except RuntimeWarning as w:
                if w.message == 'insufficient clade collections':
                    continue

            self._align_master_fasta()

            # The production of an ML tree is going to be unrealistic with the larger number of sequences
            # it will simply take too long to compute.
            # Instead I will create an NJ consensus tree.
            # This will involve using the embassy version of the phylip executables

            # First I will have to create random data sets
            self._make_fseqboot_replicate_alignments(clade_in_question)

            self._compute_unifrac_distances(clade_in_question)

            self._compute_pcoa_coords()

            self._clean_up_temp_files()

            self._append_output_files_to_output_list()

    def _append_output_files_to_output_list(self):
        self.output_path_list.extend([self.clade_unifrac_dist_file_path, self.clade_pcoa_coord_file_path])

    def _clean_up_temp_files(self):
        if os.path.exists(self.clade_fseqboot_root_dir):
            shutil.rmtree(path=self.clade_fseqboot_root_dir)

        # now delte all files except for the .csv that holds the coords and the .dist that holds the dists
        list_of_dir = os.listdir(self.clade_output_dir)
        for item in list_of_dir:
            if '.csv' not in item and '.dist' not in item:
                os.remove(os.path.join(self.clade_output_dir, item))

    def _compute_pcoa_coords(self):
        # simultaneously grab the sample names in the order of the distance matrix and put the matrix into
        # a twoD list and then convert to a numpy array
        self.clade_pcoa_coord_file_path = os.path.join(
            self.clade_output_dir, f'{self.date_time_string}.PCoA_coords.csv')
        raw_dist_file = read_defined_file_to_list(self.clade_unifrac_dist_file_path)

        temp_two_d_list = []
        sample_names_from_dist_matrix = []
        for line in raw_dist_file[1:]:
            temp_elements = line.split('\t')
            sample_names_from_dist_matrix.append(temp_elements[0].replace(' ', ''))
            temp_two_d_list.append([float(a) for a in temp_elements[1:]])

        uni_frac_dist_as_np_array = np.array(temp_two_d_list)

        sys.stdout.write('\rcalculating PCoA coordinates')

        pcoa_output = pcoa(uni_frac_dist_as_np_array)

        # rename the pcoa dataframe index as the sample names
        pcoa_output.samples['sample'] = sample_names_from_dist_matrix
        renamed_pcoa_dataframe = pcoa_output.samples.set_index('sample')

        # now add the variance explained as a final row to the renamed_dataframe
        renamed_pcoa_dataframe = renamed_pcoa_dataframe.append(pcoa_output.proportion_explained.rename('proportion_explained'))

        renamed_pcoa_dataframe.to_csv(self.clade_pcoa_coord_file_path, index=True, header=True, sep=',')

    def _compute_unifrac_distances(self, clade_in_question):
        unifrac_subclade_handler = UnifracSubcladeHandler(
            clade_in_question=clade_in_question,
            fseqbootbase=self.clade_fseqboot_base,
            parent_unifrac_creator=self)
        unifrac_subclade_handler.execute_unifrac_mothur_worker()
        unifrac_subclade_handler.create_consensus_tree()
        unifrac_subclade_handler.perform_unifrac()
        self.clade_unifrac_dist_file_path = unifrac_subclade_handler.unifrac_dist_file_path
        self._append_date_time_string_to_unifrac_dist_path()

    def _make_fseqboot_replicate_alignments(self, clade_in_question):
        fseqboot_alignment_generator = FseqbootAlignmentGenerator(clade_in_question=clade_in_question,
                                                                  parent_unifrac_dist_creator=self,
                                                                  num_reps=self.bootstrap_value)
        fseqboot_alignment_generator.do_fseqboot_alignment_generation()
        self.clade_fseqboot_base = fseqboot_alignment_generator.fseqboot_base
        self.clade_fseqboot_root_dir = fseqboot_alignment_generator.fseqboot_clade_root_dir

    def _append_date_time_string_to_unifrac_dist_path(self):
        # here add a date_time_string element to it to make it unique
        old_dist_path = self.clade_unifrac_dist_file_path
        dist_path_extension = self.clade_unifrac_dist_file_path.split('.')[-1]
        self.clade_unifrac_dist_file_path = self.clade_unifrac_dist_file_path.replace(f'.{dist_path_extension}', f'.{self.date_time_string}.{dist_path_extension}')
        os.rename(old_dist_path, self.clade_unifrac_dist_file_path)

    def _align_master_fasta(self):
        self.clade_master_fasta_file_aligned_path = self.clade_master_fasta_file_unaligned_path.replace('.fasta', '.aligned.fasta')
        mafft_align_fasta(
            input_path=self.clade_master_fasta_file_unaligned_path,
            output_path=self.clade_master_fasta_file_aligned_path,
            num_proc=self.num_proc)

    def _create_and_write_out_master_fasta_names_and_group_files(self, clade_in_question):
        sys.stdout.write('Creating master .name and .fasta files for UniFrac')
        try:
            unifrac_dist_creator_handler_one = UnifracDistanceCreatorHandlerOne(
                clade_in_question=clade_in_question,
                parent_unifrac_dist_creator=self)
        except RuntimeWarning as w:
            if w.message == 'insufficient clade collections':
                raise RuntimeWarning({'message': 'insufficient clade collections'})
            else:
                raise RuntimeError(f'Unknown error in {UnifracDistanceCreatorHandlerOne.__name__} init')
        unifrac_dist_creator_handler_one.execute_unifrac_distance_creator_worker_one()
        self.clade_master_names_file_path = unifrac_dist_creator_handler_one.master_names_file_path
        self.clade_master_fasta_file_unaligned_path = unifrac_dist_creator_handler_one.master_fasta_file_path
        self.clade_master_group_file_path = unifrac_dist_creator_handler_one.master_group_file_path

    def _setup_output_dir(self, call_type, data_set_string, output_dir):
            if call_type == 'stand_alone':
                output_dir = os.path.abspath(
                    os.path.join(self.symportal_root_dir, 'outputs', 'ordination', '_'.join(data_set_string.split(',')),
                                 'between_samples'))
            else:
                # call_type == 'submission':
                output_dir = output_dir + '/between_sample_distances'
            return output_dir

class GenericDistanceCreator:
    def __init__(self):

class BrayCurtisDisPCoACreator:
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

    def __init__(self, date_time_string, symportal_root_dir, smpl_id_list_str=None, data_set_string=None, call_type=None, output_dir=None):
        self.is_datasetsample_not_dataset_based = self._infer_is_dataset_of_datasetsample(
            smpl_id_list_str, data_set_string)
        if self.is_datasetsample_not_dataset_based:
            self._init_datasetsample_based_attributes(date_time_string, smpl_id_list_str, symportal_root_dir)
        else:
            self._init_dataset_based_attributes(call_type, data_set_string, output_dir, symportal_root_dir)
        self.output_paths_list = []
        # clade specific attributes. Will be updated for every clade processed
        self.clade_output_dir = None
        self.ccs_of_clade = None
        self.clade_dsss_seq_to_rel_abund_for_ccs_of_clade_dict = {}
        self.clade_within_clade_distances_dict = {}
        self.clade_distance_file_as_list = None
        self.clade_dist_file_output_path = None
        self.date_time_string = date_time_string

    def compute_unifrac_dists_and_pcoa_coords(self):
        for clade_in_question in self.clades_of_ccs:
            self._init_clade_dirs_and_paths(clade_in_question)

            self.ccs_of_clade = list(self.cc_list_for_output.filter(clade=clade_in_question))
            if len(self.ccs_of_clade) < 2:
                continue
            self._create_dsss_sequence_to_rel_abund_dict_for_each_cc()
            self._compute_braycurtis_btwn_cc_pairs()
            self._generate_distance_file()
            self._add_sample_uids_to_dist_file()
            self._write_out_dist_file()

    def _write_out_dist_file(self):
        write_list_to_destination(self.clade_dist_file_output_path, self.clade_distance_file_as_list)

    def _add_sample_uids_to_dist_file(self):
        # for the output version lets also append the sample name to each line so that we can see which sample it is
        # it is important that we otherwise work eith the sample ID as the sample names may not be unique.
        dist_with_sample_name = [self.clade_distance_file_as_list[0]]
        list_of_cc_ids = [int(line.split('\t')[0]) for line in self.clade_distance_file_as_list[1:]]
        cc_of_outputs = list(CladeCollection.objects.filter(id__in=list_of_cc_ids))
        dict_of_cc_id_to_sample_name = {cc.id: cc.data_set_sample_from.name for cc in cc_of_outputs}
        for line in self.clade_distance_file_as_list[1:]:
            temp_list = []
            cc_id = int(line.split('\t')[0])
            sample_name = dict_of_cc_id_to_sample_name[cc_id]
            temp_list.append(sample_name)
            temp_list.extend(line.split('\t'))
            new_line = '\t'.join(temp_list)
            dist_with_sample_name.append(new_line)
        self.clade_distance_file_as_list = dist_with_sample_name

    def _init_clade_dirs_and_paths(self, clade_in_question):
        self.clade_output_dir = os.path.join(self.output_dir, clade_in_question)
        self.clade_dist_file_output_path = os.path.join(
            self.clade_output_dir, f'{self.date_time_string}.bray_curtis_within_clade_sample_distances.dist')
        os.makedirs(self.clade_output_dir, exist_ok=True)

    def _generate_distance_file(self):
        self.clade_distance_file_as_list = [len(self.ccs_of_clade)]
        for clade_col_outer in self.ccs_of_clade:
            temp_clade_col_string = [clade_col_outer.id]

            for clade_col_inner in self.ccs_of_clade:
                if clade_col_outer == clade_col_inner:
                    temp_clade_col_string.append(0)
                else:
                    temp_clade_col_string.append(
                        self.clade_within_clade_distances_dict['{}_{}'.format(clade_col_outer.id, clade_col_inner.id)])
            self.clade_distance_file_as_list.append(
                '\t'.join([str(distance_item) for distance_item in temp_clade_col_string]))

    def _compute_braycurtis_btwn_cc_pairs(self):
        for clade_col_one, clade_col_two in itertools.combinations(list(self.ccs_of_clade), 2):
            # let's work to a virtual subsample of 100 000
            clade_col_one_seq_rel_abundance_dict = self.clade_dsss_seq_to_rel_abund_for_ccs_of_clade_dict[
                clade_col_one.id]
            clade_col_two_seq_rel_abundance_dict = self.clade_dsss_seq_to_rel_abund_for_ccs_of_clade_dict[
                clade_col_two.id]

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

            # once you have this we should simply be able to crunch the bray-curtis.
            distance = braycurtis(seq_abundance_list_one, seq_abundance_list_two)

            # these distances can be stored in a dictionary by the 'id1/id2' and 'id2/id1'
            self.clade_within_clade_distances_dict['{}_{}'.format(clade_col_one.id, clade_col_two.id)] = distance
            self.clade_within_clade_distances_dict['{}_{}'.format(clade_col_two.id, clade_col_one.id)] = distance

    def _create_dsss_sequence_to_rel_abund_dict_for_each_cc(self):
        # Go through each of the clade collections and create a dict
        # that has key as the actual sequence and relative abundance of that sequence
        # we can then store these dict in a dict where the key is the sample ID.
        for clade_col in self.ccs_of_clade:
            temp_dict = {}
            data_set_sample_sequences_of_clade_col = DataSetSampleSequence.objects.filter(
                clade_collection_found_in=clade_col)
            total_seqs_ind_clade_col = sum([dsss.abundance for dsss in data_set_sample_sequences_of_clade_col])
            for dsss in data_set_sample_sequences_of_clade_col:
                temp_dict[dsss.reference_sequence_of.sequence] = dsss.abundance / total_seqs_ind_clade_col
            self.clade_dsss_seq_to_rel_abund_for_ccs_of_clade_dict[clade_col.id] = temp_dict

    def _init_datasetsample_based_attributes(self, date_time_string, smpl_id_list_str, symportal_root_dir):
        self.dss_list = DataSetSample.objects.filter(id__in=[int(str_id) for str_id in smpl_id_list_str.split(',')])
        self.cc_list_for_output = CladeCollection.objects.filter(
            data_set_sample_from__in=self.dss_list)
        self.clades_of_ccs = list(set([a.clade for a in self.cc_list_for_output]))
        self.output_dir = os.path.abspath(
            os.path.join(
                symportal_root_dir, 'outputs', 'ordination',
                'custom_sample_list', 'between_samples', date_time_string))

    def _init_dataset_based_attributes(self, call_type, data_set_string, output_dir, symportal_root_dir):
        self.ds_list = DataSet.objects.filter(id__in=[int(a) for a in str(data_set_string).split(',')])
        self.cc_list_for_output = CladeCollection.objects.filter(
            data_set_sample_from__data_submission_from__in=self.ds_list)
        self.clades_of_ccs = list(set([a.clade for a in self.cc_list_for_output]))
        if call_type == 'stand_alone':
            self.output_dir = os.path.abspath(os.path.join(symportal_root_dir, 'outputs', 'ordination',
                                                           '_'.join(data_set_string.split(',')),
                                                           'between_samples'))
        else:
            # call_type == 'submission':
            self.output_dir = os.path.join(output_dir + 'between_sample_distances')

    def _infer_is_dataset_of_datasetsample(self, smpl_id_list_str, data_set_string):
        if smpl_id_list_str is not None and data_set_string is not None:
            raise RuntimeError('Please input smpl_id_list_str OR data_set_string')
        elif smpl_id_list_str and data_set_string:
            raise RuntimeError('Please input smpl_id_list_str OR data_set_string')
        else:
            if smpl_id_list_str:
                return True
            else:
                return False


















