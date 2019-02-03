from dbApp.models import (SymportalFramework, DataSet, ReferenceSequence,
                          DataSetSampleSequence, DataSetSample, CladeCollection)
import os
import shutil
import subprocess
import glob

class DataLoading:
    # The clades refer to the phylogenetic divisions of the Symbiodiniaceae. Most of them are represented at the genera
    # level. E.g. All members of clade C belong to the genus Cladocopium.
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    def __init__(self, data_set_uid, user_input_path, datasheet_path):
        self.dataset_object = DataSet.objects.get(id=data_set_uid)
        self.user_input_path = user_input_path
        self.output_directory = self.setup_output_directory()
        self.temp_working_directory = self.infer_temp_working_directory_path()
        self.create_temp_wkd()
        # Used for the pre minimum entropy decomposition (MED), quality controlled (QC), sequences dump
        # Within this directory we will have a directory for each sample that will contain clade
        # separated name and fasta pairs
        self.pre_med_sequence_output_directory_path = self.create_pre_med_write_out_directory_path()
        # data can be loaded either as paired fastq or fastq.gz files or as a single compressed file containing
        # the afore mentioned paired files.
        self.is_single_file_or_paired_input = self.determine_if_single_file_or_paired_input()
        self.datasheet_path = datasheet_path

    def execute(self):
        self.copy_and_decompress_input_files_to_temp_wkd()



    def copy_and_decompress_input_files_to_temp_wkd(self):
        if not self.is_single_file_or_paired_input:
            self.copy_fastq_files_from_input_dir_to_temp_wkd()

        else:
            self.extract_single_compressed_file_to_temp_wkd()

    def extract_single_compressed_file_to_temp_wkd(self):
        ext_components = self.user_input_path.split('.')
        if ext_components[-1] == 'zip':  # .zip
            subprocess.run(["unzip", self.user_input_path, '-d', self.temp_working_directory], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        elif ext_components[-2] == 'tar' and ext_components[-1] == 'gz':  # .tar.gz
            subprocess.run(["tar", "-xf", self.user_input_path, "-C", self.temp_working_directory],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        elif ext_components[-1] == 'gz':  # .gz
            subprocess.run(["gunzip", "-c", self.user_input_path, ">", self.temp_working_directory],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)

    def copy_fastq_files_from_input_dir_to_temp_wkd(self):
        file = os.listdir(self.user_input_path)[0]
        os.chdir('{}'.format(self.user_input_path))
        # * asterix are only expanded in the shell and so don't work through subprocess
        # need to use the glob library instead
        # https://stackoverflow.com/questions/13875978/python-subprocess-popen-why-does-ls-txt-not-work
        # files could be compressed (fastq.gz) or uncompressed (fastq). Either way, they should contain fastq.
        if 'fastq' in file:
            subprocess.run(['cp'] + glob.glob('*.fastq*') + [self.temp_working_directory], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        elif 'fq' in file:
            subprocess.run(['cp'] + glob.glob('*.fq*') + [self.temp_working_directory], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)

    def determine_if_single_file_or_paired_input(self):
        for file in os.listdir(self.user_input_path):
            if 'fastq' in file or 'fq' in file:
                # Then there is a fastq.gz or fastq file already uncompressed in this folder
                # In this case we will assume that the seq data is not a single file containing the pairs of files
                # rather the pairs of files themselves.
                return False
    def create_pre_med_write_out_directory_path(self):
        pre_med_write_out_directory_path = self.temp_working_directory.replace('tempData', 'pre_MED_seqs')
        os.makedirs(pre_med_write_out_directory_path, exist_ok=True)
        return pre_med_write_out_directory_path

    def setup_output_directory(self):
        output_directory = os.path.join(os.path.dirname(__file__),
                                        'outputs/data_set_submissions/{}'.format(self.dataset_object.id))
        os.makedirs(output_directory, exist_ok=True)
        return output_directory

    def infer_temp_working_directory_path(self):
        # working directory will be housed in a temp folder within the directory in which the sequencing data
        # is currently housed
        if '.' in self.user_input_path.split('/')[-1]:
            # then this path points to a file rather than a directory and we should pass through the path only
            wkd = os.path.abspath(
                '{}/tempData/{}'.format(os.path.dirname(self.user_input_path), self.dataset_object.id))
        else:
            # then we assume that we are pointing to a directory and we can directly use that to make the wkd
            wkd = os.path.abspath('{}/tempData/{}'.format(self.user_input_path, self.dataset_object.id))
        return wkd

    def create_temp_wkd(self):
        # if the directory already exists remove it and start from scratch
        if os.path.exists(self.temp_working_directory):
            shutil.rmtree(self.temp_working_directory)
        os.makedirs(self.temp_working_directory)