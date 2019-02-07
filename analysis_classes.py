import os
import subprocess
import sys
from plumbum import local
import pandas as pd
from collections import defaultdict
from dbApp.models import DataSetSample
from general import (
    decode_utf8_binary_to_list, create_dict_from_fasta, create_seq_name_to_abundance_dict_from_name_file,
    combine_two_fasta_files, remove_primer_mismatch_annotations_from_fasta, write_list_to_destination,
    read_defined_file_to_list)
from pickle import dump, load
from multiprocessing import Queue, Manager, Process
from django import db

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
            stdout_and_sterr_to_pipe=True
            ):

        self._setup_core_attributes(auto_convert_fastq_to_fasta, fastq_gz_fwd_path, fastq_gz_rev_path,
                                    input_dir, mothur_execution_path, name, name_file_path, output_dir,
                                    sequence_collection, num_processors, stdout_and_sterr_to_pipe)


        self._setup_pcr_analysis_attributes(pcr_analysis_name, pcr_fwd_primer, pcr_fwd_primer_mismatch,
                                            pcr_oligo_file_path, pcr_rev_primer, pcr_rev_primer_mismatch)

    def _setup_core_attributes(self, auto_convert_fastq_to_fasta, fastq_gz_fwd_path, fastq_gz_rev_path,
                               input_dir, mothur_execution_path, name, name_file_path, output_dir,
                               sequence_collection, num_processors, stdout_and_sterr_to_pipe):

        self._verify_that_is_either_sequence_collection_or_fastq_pair(fastq_gz_fwd_path, fastq_gz_rev_path,
                                                                      sequence_collection)

        if sequence_collection is not None:
            self._setup_sequence_collection_attribute(auto_convert_fastq_to_fasta, name, sequence_collection)
        elif sequence_collection is None:
            self._setup_fastq_attributes(fastq_gz_fwd_path, fastq_gz_rev_path)

        self._setup_remainder_of_core_attributes(input_dir, mothur_execution_path, name_file_path,
                                                 output_dir, sequence_collection, num_processors, stdout_and_sterr_to_pipe)

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
    def from_pair_of_fastq_gz_files(cls, name, fastq_gz_fwd_path, fastq_gz_rev_path,
                                    output_dir=None, mothur_execution_string='mothur', num_processors=10,
                                    stdout_and_sterr_to_pipe=True
                                    ):
        return cls(name=name, sequence_collection=None, mothur_execution_path=mothur_execution_string,
                   input_dir=os.path.dirname(os.path.abspath(fastq_gz_fwd_path)), output_dir=output_dir,
                   fastq_gz_fwd_path=fastq_gz_fwd_path, fastq_gz_rev_path=fastq_gz_rev_path,
                   name_file_path=None, num_processors=num_processors,
                   stdout_and_sterr_to_pipe=stdout_and_sterr_to_pipe)

    @classmethod
    def from_sequence_collection(cls, sequence_collection, name=None, input_dir=None,
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

    def __execute_summary(self):
        self._summarise_make_and_write_mothur_batch()
        self._run_mothur_batch_file_command()
        self.latest_summary_path = self._extract_output_path_first_line()
        self.latest_summary_output_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)

    # #####################

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


class SymNonSymTaxScreeningWorker:
    def __init__(
            self, data_loading_temp_working_directory, sample_name, data_loading_dataset_object,
            data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path,
            data_loading_pre_med_sequence_output_directory_path
    ):
        self._init_core_class_attributes(data_loading_dataset_object, data_loading_temp_working_directory, sample_name)
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
            self, data_loading_dataset_object, data_loading_temp_working_directory, sample_name):
        self.sample_name = sample_name
        self.cwd = os.path.join(data_loading_temp_working_directory, sample_name)
        self.dataset_object = DataSetSample.objects.get(
            name=sample_name, data_submission_from=data_loading_dataset_object
        )

    def execute(self):
        """This method completes the pre-med quality control.
        Having gone through the initial mothur qc, and having screened for potential symbiodiniaceae sequences
        (if running on the remote system), this method now identifies
        1 - the non symbiodiniaceae sequences and writes them out to the output dir
        2 - identifies the symbiodiniaceae sequences that violate our size range thresholds and writes them out
        also to the output dir
        3 - and finally the symbiodiniaceae sequences that do not violate our size range thresholds (sequences that
        will be carried through into med decomposition). These are also written out both as clade seperated
        .fasta and .name files in the temp working directory (for med processing), and as a non clade separated .name
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
        self._write_out_no_size_violation_seqs_clade_separated(clades_of_non_violation_seqs)

    def _write_out_no_size_violation_seqs_clade_separated(self, clades_of_non_violation_seqs):
        for clade_of_sequences_to_write_out in clades_of_non_violation_seqs:
            sequence_names_of_clade = self._get_sequence_names_of_clade_for_no_size_violation_sequences(
                clade_of_sequences_to_write_out
            )

            self._write_out_no_size_violation_clade_specific_fasta(
                clade_of_sequences_to_write_out, sequence_names_of_clade
            )

            self._write_out_no_size_violation_clade_specific_names_file(
                clade_of_sequences_to_write_out, sequence_names_of_clade
            )

    def _write_out_no_size_violation_clade_specific_names_file(
            self, clade_of_sequences_to_write_out, sequence_names_of_clade):
        sample_clade_names_file_path = os.path.join(
            self.cwd,
            clade_of_sequences_to_write_out,
            f'seqs_for_med_{self.sample_name}_clade_{clade_of_sequences_to_write_out}.names'
        )
        with open(sample_clade_names_file_path, 'w') as f:
            for sequence_name in sequence_names_of_clade:
                f.write(f'{self.name_dict[sequence_name]}\n')

    def _write_out_no_size_violation_clade_specific_fasta(
            self, clade_of_sequences_to_write_out, sequence_names_of_clade):
        sample_clade_fasta_path = os.path.join(
            self.cwd,
            clade_of_sequences_to_write_out,
            f'seqs_for_med_{self.sample_name}_clade_{clade_of_sequences_to_write_out}.fasta'
        )
        with open(sample_clade_fasta_path, 'w') as f:
            for sequence_name in sequence_names_of_clade:
                f.write(f'>{sequence_name}\n')
                f.write(f'{self.fasta_dict[sequence_name]}\n')

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

    def execute_sym_non_sym_tax_screening(self, data_loading_temp_working_directory, data_loading_dataset_object, data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path, data_loading_pre_med_sequence_output_directory_path):
        all_processes = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        sys.stdout.write('\nPerforming QC\n')

        for n in range(self.num_proc):
            p = Process(target=self.__sym_non_sym_tax_screening_worker, args=(data_loading_temp_working_directory, data_loading_dataset_object, data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path, data_loading_pre_med_sequence_output_directory_path))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

    def __sym_non_sym_tax_screening_worker(self, data_loading_temp_working_directory, data_loading_dataset_object, data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path, data_loading_pre_med_sequence_output_directory_path):

        for sample_name in iter(self.sample_name_mp_input_queue.get, 'STOP'):
            if sample_name in self.samples_that_caused_errors_in_qc_mp_list:
                continue

            sym_non_sym_tax_screening_worker_object = SymNonSymTaxScreeningWorker(
                data_loading_temp_working_directory=data_loading_temp_working_directory,
                data_loading_dataset_object=data_loading_dataset_object, sample_name=sample_name,
                data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path=
                data_loading_non_symbiodiniaceae_and_size_violation_base_directory_path,
                data_loading_pre_med_sequence_output_directory_path=data_loading_pre_med_sequence_output_directory_path
            )

            try:
                sym_non_sym_tax_screening_worker_object.execute()
            except RuntimeError as e:
                self.samples_that_caused_errors_in_qc_mp_list.append(e.sample_name)
