from collections import defaultdict
import subprocess
import os
from general import ThreadSafeGeneral

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
        self.thread_safe_general = ThreadSafeGeneral()

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
        return self.thread_safe_general.read_defined_file_to_list(self.output_file_path)

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
            self, input_dir, output_dir, name,
            fastq_gz_fwd_path, fastq_gz_rev_path, stdout_and_sterr_to_pipe):

        self.name = name
        self.fastq_gz_fwd_path = fastq_gz_fwd_path
        self.fastq_gz_rev_path = fastq_gz_rev_path
        self.fasta_path = None
        # The path to the .report file output from make.contigs
        # This is required when screening for overlap
        self.report_path = None
        # There are two screening stages, first for sequence overlap and mismatch
        # then for ambiguous sequences. We use the same methods to run the screening and pass in the self.screening_for
        # to modify the parameter that we are screening for. Once overlap screening is complete, self.screening_for
        # will be updated to 'ambig'
        self.screening_for = 'overlap'
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.name_file_path = None
        self.mothur_batch_file_path = None
        self.mothur_batch_file = None

        # we need to have seperate latest completed process objects for the actual commands and for the summaries
        # this is so that we can still extract useful information housed in the stdout from running the command
        # once the execute... function has been completed. Else, this information is lost due to it being replaced
        # by the completed_process of the summary that is automatically run after each command.
        self.latest_completed_process_command = None
        self.latest_completed_process_summary = None
        self.latest_summary_output_as_list = None
        self.latest_summary_path = None
        self.stdout_and_sterr_to_pipe = stdout_and_sterr_to_pipe
        # The SymVar primers
        self.pcr_fwd_primer = 'GAATTGCAGAACTCCGTGAACC'
        self.pcr_rev_primer = 'CGGGTTCWCTTGTYTGACTTCATGC'
        self.pcr_fwd_primer_mismatch = 2
        self.pcr_rev_primer_mismatch = 2
        self.pcr_oligo_file_path = None
        self.dot_file_file_path = None
        self.mothur_batch_file_path = os.path.join(self.input_dir, f'{self.name}_mothur_batch_file')
        # if there were error messages raised during the make.contigs
        self.make_contigs_error = False
        self.stdout_as_list = None
        self.thread_safe_general = ThreadSafeGeneral()


    # ########################################

    # main mothur commands
    def execute_make_contigs(self):
        """
        This will use the fastq_gz_fwd_path and fastq_gz_rev_paths to make a .file file that will be used
        as input to mothurs make.contigs command.
        N.B. Although in theory we can use the fastq_gz_fwd_path and the
        rev path directly as arguments to the mothur.contigs
        there appears to be a bug that doesn't allow this to work. Using a .file file is fine though. The .file file
        is in the format "path_to_file_1 path_to_file_2" i.e the paths only separated by a space.
        """
        # create .file file for the fwd fastq pair and the reverse fastq pair
        self._make_contig_make_and_write_out_dot_file()

        self._make_contig_make_and_write_mothur_batch()

        self._run_mothur_batch_file_command_make_contigs()

        self.fasta_path = self.dot_file_file_path.replace('.file', '.trim.contigs.fasta')
        self.report_path = self.dot_file_file_path.replace('.file', '.contigs.report')

        try:
            num_contigs = int(len(self.thread_safe_general.read_defined_file_to_list(self.fasta_path))/2)
        except FileNotFoundError:
            raise RuntimeError('Make.contigs out fasta not found')

        if num_contigs == 0:
            raise RuntimeError('empty fasta')
        else:
            self.__execute_summary()
            return num_contigs

    def execute_screen_seqs(self):
        """
        This will perform a mothur screen.seqs to get rid of sequences with ambiguous bases.
        """
        self._screen_seqs_make_and_write_mothur_batch_file()
        self._run_mothur_batch_file_command()

        # It may be that no names path is output if there were not seqs removed so we will
        # check the name path that would have been made to see if it exists and only update if it does
        # We will do the same for the fasta path.
        new_fasta_path = self.fasta_path.replace('.fasta', '.good.fasta')
        if self.screening_for == 'ambig':
            new_name_file_path = self.name_file_path.replace('.names', '.good.names')
            if os.path.exists(new_name_file_path):
                self.name_file_path = new_name_file_path

        if os.path.exists(new_fasta_path):
            self.fasta_path = new_fasta_path
            remaining_seqs = int(len(self.thread_safe_general.read_defined_file_to_list(self.fasta_path))/2)
            if remaining_seqs == 0 and self.screening_for == 'overlap':
                raise RuntimeError("no seqs left after overlap/mismatch seq screening")


        self.__execute_summary()

    def execute_pcr(self, do_reverse_pcr_as_well=False):
        """This will perform a mothur pcr.seqs analysis.
        if do_reverse_pcr__as_well is true then we will also reverse complement the fasta a perform the
        """

        self._pcr_make_and_write_oligo_file_if_doesnt_exist()

        self._pcr_make_and_write_mothur_batch_file()

        print(f'{self.name}: starting fwd PCR. This may take some time.')
        # Mothur now automatically removes seqs from the name file
        # This is a problem here, as we were relying on the name file not being
        # To fix this we wil asign a variable to the names file at this point in time
        # And asign this back after the PCRs.

        self._run_mothur_batch_file_command()

        fwd_output_scrapped_fasta_path = self.fasta_path.replace('.fasta', '.scrap.pcr.fasta')
        fwd_output_good_fasta_path = self.fasta_path.replace('.fasta', '.pcr.fasta')
        # Purposefully do not update the name file here.

        self.remove_primer_mismatch_annotations_from_fasta(fwd_output_good_fasta_path)

        # In some uncommon cases, all amplicons gave a PCR match and there is no scrapped fastsa to do a rev PCR on
        # We therefore check to see if the fwd_output_good_scrapped_fasta_path file exists
        # If exists, then we should clean up the output_bad_fasta
        # then reverse complement it
        # then do a pcr on it again using the same oligo set as the first run
        # we should then get the output from that pcr and add it to the previous run
        if do_reverse_pcr_as_well and self._if_scrap_fasta_exists_clean_and_write_out(fwd_output_scrapped_fasta_path):
            self.remove_primer_mismatch_annotations_from_fasta(fwd_output_scrapped_fasta_path)
            self.fasta_path = fwd_output_scrapped_fasta_path
            self._rev_comp_make_and_write_mothur_batch_file()

            self._run_mothur_batch_file_command()
            self.fasta_path = self.fasta_path.replace('.fasta', '.rc.fasta')
            self._pcr_make_and_write_mothur_batch_file()
            print(f'{self.name}: starting rev PCR. This may take some time.')
            self._run_mothur_batch_file_command()
            rev_output_good_fasta_path = self.fasta_path.replace('.fasta', '.pcr.fasta')
            self.remove_primer_mismatch_annotations_from_fasta(rev_output_good_fasta_path)
            self._make_new_fasta_path_for_fwd_rev_combined(rev_output_good_fasta_path)

            # now create a fasta that is the good fasta from both of the pcrs.
            # this will become the new mothuranalysis fasta.
            self.thread_safe_general.combine_two_fasta_files(
                path_one=fwd_output_good_fasta_path,
                path_two=rev_output_good_fasta_path,
                path_for_combined=self.fasta_path
            )
        else:
            self.fasta_path = fwd_output_good_fasta_path
        if len(self.thread_safe_general.read_defined_file_to_list(self.fasta_path)) == 0:
            raise RuntimeError('PCR fasta file is blank')


    def remove_primer_mismatch_annotations_from_fasta(self, fasta_path):
        temp_fasta = []
        fasta_to_clean = self.thread_safe_general.read_defined_file_to_list(fasta_path)
        for i in range(len(fasta_to_clean) - 1):
            if fasta_to_clean[i]:
                if fasta_to_clean[i][0] == '>' and fasta_to_clean[i + 1]:
                    if '|' in fasta_to_clean[i]:
                        temp_fasta.extend([fasta_to_clean[i].split('|')[0], fasta_to_clean[i + 1]])
                    else:
                        temp_fasta.extend([fasta_to_clean[i].split('\t')[0], fasta_to_clean[i + 1]])
        self.thread_safe_general.write_list_to_destination(fasta_path, temp_fasta)
    
    def check_fasta_and_name_valid(self):
        if self.name_file_path is None or self.fasta_path is None:
            raise RuntimeError
        if 'fasta' not in self.fasta_path or 'names' not in self.name_file_path:
            raise RuntimeError

    def execute_unique_seqs(self):

        self._unique_seqs_make_and_write_mothur_batch()

        self._run_mothur_batch_file_command()

        self.name_file_path, self.fasta_path = self._extract_output_path_two_lines()

        self.check_fasta_and_name_valid()

        self.__execute_summary()

    def execute_split_abund(self):

        self._split_abund_make_and_write_mothur_batch()

        self._run_mothur_batch_file_command()

        new_fasta_path = self.fasta_path.replace('.fasta', '.abund.fasta')
        new_name_file_path = self.name_file_path.replace('.names', '.abund.names')

        if os.path.exists(new_fasta_path):
            self.fasta_path = new_fasta_path
        if os.path.exists(new_name_file_path):
            self.name_file_path = new_name_file_path

        self.check_fasta_and_name_valid()

        self.__execute_summary()

    def __execute_summary(self):
        self._summarise_make_and_write_mothur_batch()
        self._run_mothur_batch_file_summary()

    # #####################
    def _if_scrap_fasta_exists_clean_and_write_out(self, fwd_output_scrapped_fasta_path):
        # NB Mothur will not output a scrap fasta file if there are no scrap fasta. Also NB that mothur will output
        # sequence names with no sequence for sequences that have multiple matches for a given primer.
        # we should screen for these and remove them.
        if not os.path.exists(fwd_output_scrapped_fasta_path):
            return False
        else:
            scrapped_fasta_as_list = self.thread_safe_general.read_defined_file_to_list(fwd_output_scrapped_fasta_path)
            if scrapped_fasta_as_list:
                new_scrapped_fasta = self._make_new_fasta_no_multi_match_lines(scrapped_fasta_as_list)
                if new_scrapped_fasta:
                    self.thread_safe_general.write_list_to_destination(fwd_output_scrapped_fasta_path, new_scrapped_fasta)
                    return True
                else:
                    return False
            else:
                return False

    def _make_new_fasta_no_multi_match_lines(self, scrapped_fasta_as_list):
        new_scrapped_fasta = []
        for i in range(0, len(scrapped_fasta_as_list), 2):
            if not 'multipleMatches' in scrapped_fasta_as_list[i] and len(scrapped_fasta_as_list[i + 1]) > 1:
                new_scrapped_fasta.extend([scrapped_fasta_as_list[i], scrapped_fasta_as_list[i + 1]])
        return new_scrapped_fasta


    def _split_abund_extract_output_path_name_and_fasta(self):
        stdout_string_as_list = self.thread_safe_general.decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):

            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i + 2], stdout_string_as_list[i + 4]

    def _split_abund_make_and_write_mothur_batch(self):
        self.mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'split.abund(fasta={self.fasta_path}, name={self.name_file_path}, cutoff=2)'
        ]

        self.thread_safe_general.write_list_to_destination(self.mothur_batch_file_path, self.mothur_batch_file)

    def _unique_seqs_make_and_write_mothur_batch(self):
        if self.name_file_path:
            self.mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'unique.seqs(fasta={self.fasta_path}, name={self.name_file_path})'
            ]
        else:
            self.mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'unique.seqs(fasta={self.fasta_path})'
            ]
        self.thread_safe_general.write_list_to_destination(self.mothur_batch_file_path, self.mothur_batch_file)

    def _summarise_make_and_write_mothur_batch(self):
        if self.name_file_path:
            self.mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'summary.seqs(fasta={self.fasta_path}, name={self.name_file_path}, processors=1)'
            ]
        else:
            self.mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'summary.seqs(fasta={self.fasta_path}, processors=1)'
            ]

        self.thread_safe_general.write_list_to_destination(self.mothur_batch_file_path, self.mothur_batch_file)

    def _make_contig_make_and_write_mothur_batch(self):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'make.contigs(file={self.dot_file_file_path}, mismatch=-4, gapopen=-4, processors=1)'
        ]
        self.thread_safe_general.write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def _make_contig_make_and_write_out_dot_file(self):
        dot_file_file = [f'{self.fastq_gz_fwd_path} {self.fastq_gz_rev_path}']
        self.dot_file_file_path = os.path.join(self.input_dir, f'{self.name}_fastq_pair.file')
        self.thread_safe_general.write_list_to_destination(self.dot_file_file_path, dot_file_file)

    def _make_new_fasta_path_for_fwd_rev_combined(self, rev_output_good_fasta_path):
        self.fasta_path = rev_output_good_fasta_path.replace('.scrap.pcr.rc.pcr', '.pcr.combined')

    def _rev_comp_make_and_write_mothur_batch_file(self):
        self.mothur_batch_file = self._make_rev_complement_mothur_batch_file()
        self.thread_safe_general.write_list_to_destination(self.mothur_batch_file_path, self.mothur_batch_file)

    def _make_rev_complement_mothur_batch_file(self):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'reverse.seqs(fasta={self.fasta_path})'
        ]
        return mothur_batch_file

    def _extract_output_paths_screen_seqs_with_name_file_command(self, output_as_list=None):
        if output_as_list is None:
            stdout_string_as_list = self.thread_safe_general.decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        else:
            stdout_string_as_list = output_as_list
        for i in range(len(stdout_string_as_list)):
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i+1], stdout_string_as_list[i+3]

    def _extract_output_path_first_line_command(self, output_as_list=None):
        if output_as_list is None:
            stdout_string_as_list = self.thread_safe_general.decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        else:
            stdout_string_as_list = output_as_list
        for i in range(len(stdout_string_as_list)):
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i+1]

    def _extract_output_path_second_line_command(self):
        stdout_string_as_list = self.thread_safe_general.decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i+2]

    def _extract_output_path_first_line_summary(self):
        stdout_string_as_list = self.thread_safe_general.decode_utf8_binary_to_list(self.latest_completed_process_summary.stdout)
        for i in range(len(stdout_string_as_list)):
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i+1]

    def _extract_output_path_two_lines(self):
        stdout_string_as_list = self.thread_safe_general.decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i + 1], stdout_string_as_list[i + 2]

    def _pcr_extract_good_and_scrap_output_paths(self):
        stdout_string_as_list = self.thread_safe_general.decode_utf8_binary_to_list(self.latest_completed_process_command.stdout)
        for i in range(len(stdout_string_as_list)):
            if 'Output File Names' in stdout_string_as_list[i]:
                output_good_fasta_path = stdout_string_as_list[i + 1]
                output_scrapped_fasta_path = stdout_string_as_list[i + 3]
                return output_scrapped_fasta_path, output_good_fasta_path

    def _run_mothur_batch_file_command(self):
        if self.stdout_and_sterr_to_pipe:
            self.latest_completed_process_command = subprocess.run(
                ['mothur', self.mothur_batch_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        else:
            self.latest_completed_process_command = subprocess.run(
                ['mothur', self.mothur_batch_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            for line in self.thread_safe_general.decode_utf8_binary_to_list(self.latest_completed_process_command.stdout):
                print(line)

    def _run_mothur_batch_file_command_make_contigs(self):
        """Now that we are running the new version of mothur, we are hoping that this process
        no longer gets stuck in a loop. We will let it run and check the output fasta to see
        if we need to raise a warning based on the number """
        self.stdout_as_list = []

        warning_count = 0
        self.latest_completed_process_command = subprocess.Popen(['mothur', self.mothur_batch_file_path],
                                                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for byte_line in self.latest_completed_process_command.stdout:
            byte_line_as_string = byte_line.decode('ISO-8859-1')
            self.stdout_as_list.append(byte_line_as_string.rstrip())
            if 'WARNING' in byte_line_as_string:
                warning_count += 1
                if warning_count > 10000000:
                    self.latest_completed_process_command.kill()
                    raise RuntimeError('bad fastq, mothur stuck in loop')
            if 'ERROR' in byte_line_as_string:
                self.make_contigs_error = True
        self.latest_completed_process_command.wait()
        if self.latest_completed_process_command.returncode == 1:
            raise RuntimeError('error in make.contigs')
        if not self.stdout_and_sterr_to_pipe:
            for line in self.stdout_as_list:
                print(line)


    def _run_mothur_batch_file_summary(self):
        if self.stdout_and_sterr_to_pipe:
            self.latest_completed_process_summary = subprocess.run(
                ['mothur', self.mothur_batch_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        else:
            self.latest_completed_process_summary = subprocess.run(
                ['mothur', self.mothur_batch_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            for line in self.thread_safe_general.decode_utf8_binary_to_list(self.latest_completed_process_summary.stdout):
                print(line)

    def _pcr_make_and_write_mothur_batch_file(self):
        self._pcr_make_mothur_batch_file()
        self.thread_safe_general.write_list_to_destination(self.mothur_batch_file_path, self.mothur_batch_file)

    def _screen_seqs_make_and_write_mothur_batch_file(self):
        if self.screening_for == 'overlap':
            self._screen_seqs_make_mothur_batch_file_overlap()
        elif self.screening_for == 'ambig':
            self._screen_seqs_make_mothur_batch_file_ambig()
        else:
            raise RuntimeError("Unrecognized 'screening for' value")
        self.thread_safe_general.write_list_to_destination(self.mothur_batch_file_path, self.mothur_batch_file)


    def _screen_seqs_make_mothur_batch_file_ambig(self):
        self.mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'screen.seqs(fasta={self.fasta_path}, name={self.name_file_path}, maxambig=0, processors=1)'
        ]

    def _screen_seqs_make_mothur_batch_file_overlap(self):
        self.mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'screen.seqs(fasta={self.fasta_path}, '
            f'contigsreport={self.report_path}, minoverlap=30, processors=1)'
        ]

    def _pcr_make_mothur_batch_file(self):
        self.mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'pcr.seqs(fasta={self.fasta_path}, name={self.name_file_path}, oligos={self.pcr_oligo_file_path}, '
            f'pdiffs={self.pcr_fwd_primer_mismatch}, rdiffs={self.pcr_rev_primer_mismatch}, '
            f'processors=1)'
        ]

    def _pcr_make_and_write_oligo_file_if_doesnt_exist(self):
        if self.pcr_oligo_file_path is None:
            oligo_file = [
                f'forward\t{self.pcr_fwd_primer}',
                f'reverse\t{self.pcr_rev_primer}'
            ]
            self.pcr_oligo_file_path = os.path.join(self.input_dir, f'{self.name}_oligo_file.oligo')
            self.thread_safe_general.write_list_to_destination(self.pcr_oligo_file_path, oligo_file)


class NucleotideSequence:
    def __init__(self, sequence, name=None, abundance=None):
        self.sequence = sequence
        self.length = len(sequence)
        self.name = name
        self.abundance = abundance


