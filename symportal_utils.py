from collections import defaultdict
import subprocess
import os
import sys
from general import (read_defined_file_to_list, write_list_to_destination, create_dict_from_fasta,
                     combine_two_fasta_files, decode_utf8_binary_to_list,
                     create_seq_name_to_abundance_dict_from_name_file)

class BlastnAnalysis:
    def __init__(
            self, input_file_path, output_file_path,
            db_path='/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload/nt', max_target_seqs=1,
            num_threads=1, output_format_string="6 qseqid sseqid staxids evalue pident qcovs staxid stitle ssciname",
            blastn_exec_path='blastn', pipe_stdout_sterr=True
    ):

        self.input_file_path = input_file_path
        self.output_file_path = output_file_path
        self.db_path = db_path
        self.output_format_string = output_format_string
        self.max_target_seqs = max_target_seqs
        self.num_threads = num_threads
        self.blastn_exec_path = blastn_exec_path
        self.pipe_stdout_sterr = pipe_stdout_sterr

    def execute_blastn_analysis(self, pipe_stdout_sterr=True):
        if not pipe_stdout_sterr:
            completed_process = subprocess.run([
                self.blastn_exec_path, '-out', self.output_file_path, '-outfmt', self.output_format_string, '-query',
                self.input_file_path, '-db', self.db_path,
                '-max_target_seqs', f'{self.max_target_seqs}', '-num_threads', f'{self.num_threads}'])

        elif not self.pipe_stdout_sterr:
            completed_process = subprocess.run([
                self.blastn_exec_path, '-out', self.output_file_path, '-outfmt', self.output_format_string, '-query',
                self.input_file_path, '-db', self.db_path,
                '-max_target_seqs', f'{self.max_target_seqs}', '-num_threads', f'{self.num_threads}'])
        else:
            completed_process = subprocess.run([
                self.blastn_exec_path, '-out', self.output_file_path, '-outfmt', self.output_format_string, '-query',
                self.input_file_path, '-db', self.db_path,
                '-max_target_seqs', f'{self.max_target_seqs}', '-num_threads', f'{self.num_threads}'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return completed_process

    def return_blast_output_as_list(self):
        return read_defined_file_to_list(self.output_file_path)

    def return_blast_results_dict(self):
        blast_output_file_as_list = self.return_blast_output_as_list()
        blast_output_dict = defaultdict(list)
        for line in blast_output_file_as_list:
            blast_output_dict[line.split('\t')[0]].append('\t'.join(line.split('\t')[1:]))
        return blast_output_dict

    def make_db(self, title_for_db):
        if self.pipe_stdout_sterr:
            subprocess.run(
                ['makeblastdb', '-in', self.db_path, '-dbtype', 'nucl', '-title',
                 title_for_db], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            subprocess.run(
                ['makeblastdb', '-in', self.db_path, '-dbtype', 'nucl', '-title',
                 title_for_db])


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
                                    sequence_collection, num_processors, stdout_and_sterr_to_pipe, is_unifrac_analysis,
                                    tree_file_path, group_file_path)

        self._setup_pcr_analysis_attributes(pcr_analysis_name, pcr_fwd_primer, pcr_fwd_primer_mismatch,
                                            pcr_oligo_file_path, pcr_rev_primer, pcr_rev_primer_mismatch)

    def _setup_core_attributes(self, auto_convert_fastq_to_fasta, fastq_gz_fwd_path, fastq_gz_rev_path,
                               input_dir, mothur_execution_path, name, name_file_path, output_dir,
                               sequence_collection, num_processors, stdout_and_sterr_to_pipe, is_unifrac_analysis,
                               tree_file_path, group_file_path):

        if is_unifrac_analysis:
            self._setup_unifrac_attributes(
                tree_file_path, group_file_path, name_file_path, input_dir, output_dir,
                mothur_execution_path, num_processors, stdout_and_sterr_to_pipe, name)
        else:
            self._verify_that_is_either_sequence_collection_or_fastq_pair(fastq_gz_fwd_path, fastq_gz_rev_path,
                                                                          sequence_collection)
            if sequence_collection is not None:
                self._setup_sequence_collection_attribute(auto_convert_fastq_to_fasta, name, sequence_collection)
            elif sequence_collection is None:
                self._setup_fastq_attributes(fastq_gz_fwd_path, fastq_gz_rev_path, name)
            self._setup_remainder_of_core_attributes(
                input_dir, mothur_execution_path, name_file_path, output_dir, sequence_collection,
                num_processors, stdout_and_sterr_to_pipe)

    def _setup_unifrac_attributes(
            self, tree_path, group_file_path, name_file_path, input_dir, output_dir,
            mothur_execution_path, num_processors, stdout_and_sterr_to_pipe, name):
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
        self.fasta_path = None
        self.name = name

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
            self.output_dir = output_dir

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

    def _setup_fastq_attributes(self, fastq_gz_fwd_path, fastq_gz_rev_path, name):
        self.fastq_gz_fwd_path = fastq_gz_fwd_path
        self.fastq_gz_rev_path = fastq_gz_rev_path
        self.sequence_collection = None
        self.fasta_path = None
        self.name = name

    def _setup_sequence_collection_attribute(self, auto_convert_fastq_to_fasta, name, sequence_collection):
        self.fastq_gz_fwd_path = None
        self.fastq_gz_rev_path = None
        if sequence_collection.file_type == 'fastq':
            self._convert_to_fasta_or_raise_value_error(auto_convert_fastq_to_fasta, sequence_collection)
        if name is None:
            if sequence_collection.name is not None:
                self.name = sequence_collection.name
            else:
                raise RuntimeError('A name must be provided for the MothurAnalysis as one was not provided for the passed SequenceCollection')
        else:
            self.name = name
        self.sequence_collection = sequence_collection
        self.fasta_path = self.sequence_collection.file_path

    @staticmethod
    def _convert_to_fasta_or_raise_value_error(auto_convert_fastq_to_fasta, sequence_collection):
        if auto_convert_fastq_to_fasta:
            print('SequenceCollection must be of type fasta\n. Running SeqeunceCollection.convert_to_fasta.\n')
            sequence_collection.convert_to_fasta()
        else:
            ValueError('SequenceCollection must be of type fasta. You can use the SequenceCollection')

    @staticmethod
    def _verify_that_is_either_sequence_collection_or_fastq_pair(fastq_gz_fwd_path, fastq_gz_rev_path,
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
                self.pcr_rev_primer = 'CGGGTTCWCTTGTYTGACTTCATGC',

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
    def init_from_pair_of_fastq_gz_files(
            cls, name, fastq_gz_fwd_path, fastq_gz_rev_path, output_dir=None, mothur_execution_string='mothur',
            num_processors=1, stdout_and_sterr_to_pipe=True, pcr_fwd_primer=None, pcr_rev_primer=None,
            pcr_oligo_file_path=None, pcr_fwd_primer_mismatch=2, pcr_rev_primer_mismatch=2, pcr_analysis_name=None, input_dir=None):

        return cls(name=name, sequence_collection=None, mothur_execution_path=mothur_execution_string,
                   input_dir=input_dir, output_dir=output_dir,
                   fastq_gz_fwd_path=fastq_gz_fwd_path, fastq_gz_rev_path=fastq_gz_rev_path,
                   name_file_path=None, num_processors=num_processors,
                   stdout_and_sterr_to_pipe=stdout_and_sterr_to_pipe, pcr_fwd_primer=pcr_fwd_primer,
                   pcr_rev_primer=pcr_rev_primer, pcr_oligo_file_path=pcr_oligo_file_path,
                   pcr_fwd_primer_mismatch=pcr_fwd_primer_mismatch, pcr_rev_primer_mismatch=pcr_rev_primer_mismatch,
                   pcr_analysis_name=pcr_analysis_name)

    @classmethod
    def init_from_sequence_collection(cls, sequence_collection, name=None, input_dir=None,
                                      output_dir=None, mothur_execution_path='mothur',
                                      pcr_fwd_primer=None, pcr_rev_primer=None, pcr_oligo_file_path=None,
                                      pcr_fwd_primer_mismatch=2, pcr_rev_primer_mismatch=2, pcr_analysis_name=None,
                                      num_processors=1, stdout_and_sterr_to_pipe=True):
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
            cls, tree_path, group_file_path, name_file_path, name, input_dir=None,
            output_dir=None, mothur_execution_string='mothur', num_processors=1,
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
        good_fasta_path = self._extract_output_path_first_line_command()
        if good_fasta_path is None:
            raise RuntimeError
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

        remove_primer_mismatch_annotations_from_fasta(fwd_output_good_fasta_path)

        # In some uncommon cases, all amplicons gave a PCR match and there is no scrapped fastsa to do a rev PCR on
        # then we should clean up the output_bad_fasta
        # then reverse complement it
        # then do a pcr on it again using the same oligo set as the first run
        # we should then get the output from that pcr and add it to the previous run
        if do_reverse_pcr_as_well and self._if_scrap_fasta_exists_clean_and_write_out(fwd_output_scrapped_fasta_path):
            remove_primer_mismatch_annotations_from_fasta(fwd_output_scrapped_fasta_path)
            self.fasta_path = fwd_output_scrapped_fasta_path
            self._rev_comp_make_and_write_mothur_batch_file()
            self._run_mothur_batch_file_command()
            self.fasta_path = self._extract_output_path_first_line_command()
            if self.fasta_path is None:
                raise RuntimeError
            self._pcr_make_and_write_mothur_batch_file()
            self._run_mothur_batch_file_command()
            rev_output_good_fasta_path = self._pcr_extract_good_and_scrap_output_paths()[1]
            remove_primer_mismatch_annotations_from_fasta(rev_output_good_fasta_path)
            self._make_new_fasta_path_for_fwd_rev_combined(rev_output_good_fasta_path)

            # now create a fasta that is the good fasta from both of the pcrs.
            # this will become the new mothuranalysis fasta.
            combine_two_fasta_files(
                path_one=fwd_output_good_fasta_path,
                path_two=rev_output_good_fasta_path,
                path_for_combined=self.fasta_path
            )
        else:
            self.fasta_path = fwd_output_good_fasta_path
        if len(read_defined_file_to_list(self.fasta_path)) == 0:
            raise RuntimeError('PCR fasta file is blank')
        if self.name_file_path:
            self._update_sequence_collection_from_fasta_name_pair()
        else:
            self._update_sequence_collection_from_fasta_file()

    def execute_make_contigs(self):
        """
        This will use the fastq_gz_fwd_path and fastq_gz_rev_paths to make a .file file that will be used
        as input to mothurs make.contigs command.
        N.B. Although in theory we can use the fastq_gz_fwd_path and the
        rev path directly as arguments to the mothur.contigs
        there appears to be a bug that doesn't allow this to work. Using a .file file is fine though. The .file file
        is in the format "path_to_file_1 path_to_file_2" i.e the paths only separated by a space.
        :return:
        """
        # create .file file for the fwd fastq pair and the reverse fastq pair
        dot_file_file_path = self._make_contig_make_and_write_out_dot_file()

        self._make_contig_make_and_write_mothur_batch(dot_file_file_path)

        self._run_mothur_batch_file_command()

        self.fasta_path = self._extract_output_path_first_line_command()
        if self.fasta_path is None:
            raise RuntimeError

        self._update_sequence_collection_from_fasta_file()

        self.__execute_summary()

    def execute_unique_seqs(self):

        self._unique_seqs_make_and_write_mothur_batch()

        self._run_mothur_batch_file_command()

        self.name_file_path, self.fasta_path = self._extract_output_path_two_lines()
        if self.name_file_path is None or self.fasta_path is None:
            raise RuntimeError

        self._update_sequence_collection_from_fasta_name_pair()

        self.__execute_summary()

    def execute_split_abund(self, abund_cutoff=2):

        self._split_abund_make_and_write_mothur_batch(abund_cutoff)

        self._run_mothur_batch_file_command()

        self.name_file_path, self.fasta_path = self._split_abund_extract_output_path_name_and_fasta()
        if self.name_file_path is None or self.fasta_path is None:
            raise RuntimeError

        self._update_sequence_collection_from_fasta_name_pair()

        self.__execute_summary()

    def execute_dist_seqs(self):
        self._dist_seqs_make_and_write_mothur_batch()

        self._run_mothur_batch_file_command()

        self.dist_file_path = self._extract_output_path_first_line_command()
        if self.dist_file_path is None:
            raise RuntimeError('Dist file is None')


    def execute_clearcut(self):
        self._validate_dist_file()
        self._clearcut_make_and_write_mothur_batch()
        self._run_mothur_batch_file_command()
        self.tree_file_path = self._extract_output_path_first_line_command()
        if self.tree_file_path is None:
            raise RuntimeError('Tree file path is None')


    def execute_weighted_unifrac(self):
        self._weighted_unifrac_make_and_write_mothur_batch()
        self._run_mothur_batch_file_command()
        if self.latest_completed_process_command.returncode == 0:
            sys.stdout.write('\rUnifrac successful')
        else:
            sys.stdout.write('\rERROR: {}'.format(self.latest_completed_process_command.sterr.decode('ISO-8859-1')))

        self.dist_file_path = self._extract_output_path_second_line_command()

    def __execute_summary(self):
        self._summarise_make_and_write_mothur_batch()
        self._run_mothur_batch_file_summary()
        self.latest_summary_path = self._extract_output_path_first_line_summary()
        self.latest_summary_output_as_list = decode_utf8_binary_to_list(self.latest_completed_process_summary.stdout)

    # #####################
    def _if_scrap_fasta_exists_clean_and_write_out(self, fwd_output_scrapped_fasta_path):
        # NB Mothur will not output a scrap fasta file if there are no scrap fasta. Also NB that mothur will output
        # sequence names with no sequence for sequences that have multiple matches for a given primer.
        # we should screen for these and remove them.
        if fwd_output_scrapped_fasta_path == '':
            return False
        else:
            scrapped_fasta_as_list = read_defined_file_to_list(fwd_output_scrapped_fasta_path)
            if scrapped_fasta_as_list:
                new_scrapped_fasta = self._make_new_fasta_no_multi_match_lines(scrapped_fasta_as_list)
                if new_scrapped_fasta:
                    write_list_to_destination(fwd_output_scrapped_fasta_path, new_scrapped_fasta)
                    return True
                else:
                    return False
            else:
                return False

    def _make_new_fasta_no_multi_match_lines(self, scrapped_fasta_as_list):
        new_scrapped_fasta = []
        for i in range(0, len(scrapped_fasta_as_list), 2):
            if not 'multipleMatches' in scrapped_fasta_as_list[i]:
                new_scrapped_fasta.extend([scrapped_fasta_as_list[i], scrapped_fasta_as_list[i + 1]])
        return new_scrapped_fasta

    def _weighted_unifrac_make_and_write_mothur_batch(self):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'unifrac.weighted(tree={self.tree_file_path}, group={self.group_file_path}, name={self.name_file_path},'
            f' distance=square, processors={self.processors})'
        ]

        self._set_mothur_batch_file_path()

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

        self._set_mothur_batch_file_path()

        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _dist_seqs_make_and_write_mothur_batch(self):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'dist.seqs(fasta={self.fasta_path}, countends=T, output=square)'
        ]

        self._set_mothur_batch_file_path()

        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _split_abund_extract_output_path_name_and_fasta(self):
        stdout_string_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):

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
        self._set_mothur_batch_file_path()
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
        self._set_mothur_batch_file_path()
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _set_mothur_batch_file_path(self):
        self.mothur_batch_file_path = os.path.join(self.input_dir, f'{self.name}_mothur_batch_file')

    def _summarise_make_and_write_mothur_batch(self):
        if self.name_file_path:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'summary.seqs(fasta={self.fasta_path}, name={self.name_file_path})'
            ]
        else:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'summary.seqs(fasta={self.fasta_path})'
            ]
        self._set_mothur_batch_file_path()
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _make_contig_make_and_write_mothur_batch(self, dot_file_file_path):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'make.contigs(file={dot_file_file_path})'
        ]
        self._set_mothur_batch_file_path()
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _make_contig_make_and_write_out_dot_file(self):
        dot_file_file = [f'{self.fastq_gz_fwd_path} {self.fastq_gz_rev_path}']
        dot_file_file_path = os.path.join(self.input_dir, f'{self.name}_fastq_pair.file')
        write_list_to_destination(dot_file_file_path, dot_file_file)
        return dot_file_file_path

    def _make_new_fasta_path_for_fwd_rev_combined(self, rev_output_good_fasta_path):
        self.fasta_path = rev_output_good_fasta_path.replace('.scrap.pcr.rc.pcr', '.pcr.combined')

    def _update_sequence_collection_from_fasta_file(self):
        if self.sequence_collection is None:
            # then we need to create a sequence collection from the new fasta path
            self.sequence_collection = SequenceCollection(path_to_file=self.fasta_path)
        else:
            self.sequence_collection.set_list_of_nucleotide_sequences_from_fasta_or_fastq(self.fasta_path)

    def _update_sequence_collection_from_fasta_name_pair(self):
        self.sequence_collection.generate_sequence_collection_from_fasta_name_pair(fasta_file_path=self.fasta_path, name_file_path=self.name_file_path)

    def _rev_comp_make_and_write_mothur_batch_file(self):
        mothur_batch_file = self._make_rev_complement_mothur_batch_file()
        self._set_mothur_batch_file_path()
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _make_rev_complement_mothur_batch_file(self):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'reverse.seqs(fasta={self.fasta_path})'
        ]
        return mothur_batch_file

    def _extract_output_path_first_line_command(self):
        stdout_string_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i+1]

    def _extract_output_path_second_line_command(self):
        stdout_string_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i+2]

    def _extract_output_path_first_line_summary(self):
        stdout_string_as_list = decode_utf8_binary_to_list(self.latest_completed_process_summary.stdout)
        for i in range(len(stdout_string_as_list)):
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i+1]

    def _extract_output_path_two_lines(self):
        stdout_string_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i + 1], stdout_string_as_list[i + 2]

    def _pcr_extract_good_and_scrap_output_paths(self):
        stdout_string_as_list = decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            if 'Output File Names' in stdout_string_as_list[i]:
                output_good_fasta_path = stdout_string_as_list[i + 1]
                output_scrapped_fasta_path = stdout_string_as_list[i + 3]
                return output_scrapped_fasta_path, output_good_fasta_path

    def _run_mothur_batch_file_command(self):
        """Run the mothur batch file that does make.contigs. NB that with some dodgee fastq pairs, mothur gets stuck
        in a loop printing out warnings. We will therefore read through the stdout and look for warnings and kill
        the process if these warnings get too high."""
        if self.stdout_and_sterr_to_pipe:
            warning_count = 0
            with subprocess.Popen([self.exec_path, self.mothur_batch_file_path],stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
                for byte_line in proc.stdout:
                    if 'WARNING' in byte_line.decode('ISO-8859-1'):
                        warning_count += 1
                        if warning_count > 100:
                            proc.kill()
                            raise RuntimeError('bad fastq, mothur stuck in loop')

        else:
            warning_count = 0
            with subprocess.Popen([self.exec_path, self.mothur_batch_file_path], stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE) as proc:
                for byte_line in proc.stdout:
                    line_str = byte_line.decode('ISO-8859-1')
                    print(line_str)
                    if 'WARNING' in line_str:
                        warning_count += 1
                        if warning_count > 100:
                            proc.kill()
                            raise RuntimeError('bad fastq, mothur stuck in loop')

    def _run_mothur_batch_file_summary(self):
        if self.stdout_and_sterr_to_pipe:
            self.latest_completed_process_summary = subprocess.run(
                [self.exec_path, self.mothur_batch_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        else:
            self.latest_completed_process_summary = subprocess.run(
                [self.exec_path, self.mothur_batch_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            for line in decode_utf8_binary_to_list(self.latest_completed_process_summary.stdout):
                print(line)

    def _pcr_make_and_write_mothur_batch_file(self):
        mothur_batch_file = self._pcr_make_mothur_batch_file()
        self._set_mothur_batch_file_path()
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _screen_seqs_make_and_write_mothur_batch_file(self, argument_dictionary):
        mothur_batch_file = self._screen_seqs_make_mothur_batch_file(argument_dictionary)
        self._set_mothur_batch_file_path()
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    @staticmethod
    def _screen_seqs_create_additional_arguments_string(argument_dict):
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
                f'pdiffs={self.pcr_fwd_primer_mismatch}, rdiffs={self.pcr_rev_primer_mismatch}, '
                f'processors={self.processors})'
            ]

        else:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'pcr.seqs(fasta={self.fasta_path}, oligos={self.pcr_oligo_file_path}, '
                f'pdiffs={self.pcr_fwd_primer_mismatch}, rdiffs={self.pcr_rev_primer_mismatch}, '
                f'processors={self.processors})'
            ]
        return mothur_batch_file

    def _pcr_make_and_write_oligo_file_if_doesnt_exist(self):
        if self.pcr_oligo_file_path is None:
            oligo_file = [
                f'forward\t{self.pcr_fwd_primer}',
                f'reverse\t{self.pcr_rev_primer}'
            ]
            self.pcr_oligo_file_path = os.path.join(self.input_dir, f'{self.name}_oligo_file.oligo')
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
    def __init__(self, path_to_file, name=None, auto_convert_to_fasta=True):
        self.name = name
        self.file_path = path_to_file
        # self.file_as_list = read_defined_file_to_list(self.file_path)
        self.file_type = self.infer_file_type()
        self.list_of_nucleotide_sequences = None
        self.set_list_of_nucleotide_sequences_from_fasta_or_fastq()
        if auto_convert_to_fasta:
            if self.file_type != 'fasta':
                self.convert_to_fasta()

    def convert_to_fasta(self):
        self.file_path = self.write_out_as_fasta()
        self.file_type = 'fasta'

    def __len__(self):
        return len(self.list_of_nucleotide_sequences)

    def write_out_as_fasta(self):
        self.file_path = self.infer_fasta_path_from_current_fastq_path()
        write_list_to_destination(destination=self.file_path, list_to_write=self.as_fasta())

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
        seq_name_to_abundace_dict = create_seq_name_to_abundance_dict_from_name_file(name_file_path=name_file_path)
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
                    self.create_new_nuc_seq_object_and_add_to_list(
                        fastq_file_as_list, i, list_of_nuleotide_sequence_objects)
        self.list_of_nucleotide_sequences = list_of_nuleotide_sequence_objects

    @staticmethod
    def is_fastq_defline(fastsq_file, index_value):
        if fastsq_file[index_value].startswith('@') and fastsq_file[index_value + 2][0] == '+':
            return True

    def create_new_nuc_seq_object_and_add_to_list(
            self, fastq_file_as_list, index_val, list_of_nuleotide_sequence_objects):
        name, sequence = self.get_single_fastq_info_from_fastq_file_by_index(fastq_file_as_list, index_val)
        list_of_nuleotide_sequence_objects.append(NucleotideSequence(sequence=sequence, name=name))

    @staticmethod
    def get_single_fastq_info_from_fastq_file_by_index(fastq_file_as_list, index_val):
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

def remove_primer_mismatch_annotations_from_fasta(fasta_path):
    temp_fasta = []
    fasta_to_clean = read_defined_file_to_list(fasta_path)
    for i in range(len(fasta_to_clean) - 1):
        if fasta_to_clean[i]:
            if fasta_to_clean[i][0] == '>' and fasta_to_clean[i + 1]:
                if '|' in fasta_to_clean[i]:
                    temp_fasta.extend([fasta_to_clean[i].split('|')[0], fasta_to_clean[i + 1]])
                else:
                    temp_fasta.extend([fasta_to_clean[i].split('\t')[0], fasta_to_clean[i + 1]])
    write_list_to_destination(fasta_path, temp_fasta)