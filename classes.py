from dbApp.models import (SymportalFramework, DataSet, ReferenceSequence,
                          DataSetSampleSequence, DataSetSample, CladeCollection)
import sys
import os
import shutil
import subprocess
import glob
import pandas as pd
import general
import json

class DataLoading:
    # The clades refer to the phylogenetic divisions of the Symbiodiniaceae. Most of them are represented at the genera
    # level. E.g. All members of clade C belong to the genus Cladocopium.
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    def __init__(self, data_set_uid, user_input_path, datasheet_path):
        self.dataset_object = DataSet.objects.get(id=data_set_uid)
        self.user_input_path = user_input_path
        self.output_directory = self.setup_output_directory()
        self.set_temp_working_directory()
        self.create_temp_wkd()
        # Used for the pre minimum entropy decomposition (MED), quality controlled (QC), sequences dump
        # Within this directory we will have a directory for each sample that will contain clade
        # separated name and fasta pairs
        self.pre_med_sequence_output_directory_path = self.create_pre_med_write_out_directory_path()
        # data can be loaded either as paired fastq or fastq.gz files or as a single compressed file containing
        # the afore mentioned paired files.
        self.is_single_file_or_paired_input = self.determine_if_single_file_or_paired_input()
        self.datasheet_path = datasheet_path
        self.sample_meta_info_df = None
        self.list_of_samples_names = None
        self.list_of_fastq_file_names_in_wkd = None
        self.list_of_fastq_files_in_wkd = None
        self.path_to_mothur_batch_file_for_dot_file_creation = None
        self.fastqs_are_gz_compressed = None
        self.path_to_latest_mothur_batch_file = None
        self.fastq_file_to_sample_name_dict = None
        self.num_of_samples = None

    def execute(self):
        self.copy_and_decompress_input_files_to_temp_wkd()


        if self.datasheet_path:
            self.generate_stability_file_and_data_set_sample_objects_with_datasheet()
        else:
            self.generate_stability_file_and_data_set_sample_objects_without_datasheet()

    def generate_stability_file_and_data_set_sample_objects_without_datasheet(self):

        self.list_of_fastq_file_names_in_wkd = [a for a in os.listdir(self.temp_working_directory) if 'fastq' in a]

        self.generate_and_write_mothur_batch_file_for_dotfile_creation()
        general.execute_mothur_batch_file_with_piped_stoud_sterr(self.path_to_latest_mothur_batch_file)

        self.generate_and_write_new_stability_file_without_datasheet()

        self.create_data_set_sample_objects_in_bulk_without_datasheet()

    def create_data_set_sample_objects_in_bulk_without_datasheet(self):
        list_of_sample_objects = []
        sys.stdout.write('\nCreating data_set_sample objects\n')
        for sampleName in self.list_of_samples_names:
            print('\rCreating data_set_sample {}'.format(sampleName))
            # Create the data_set_sample objects in bulk.
            # The cladal_seq_totals property of the data_set_sample object keeps track of the seq totals for the
            # sample divided by clade. This is used in the output to keep track of sequences that are not
            # included in cladeCollections
            clade_zeroes_list = [0 for _ in self.clade_list]
            empty_cladal_seq_totals = json.dumps(clade_zeroes_list)

            dss = DataSetSample(name=sampleName, data_submission_from=self.dataset_object,
                                cladal_seq_totals=empty_cladal_seq_totals)
            list_of_sample_objects.append(dss)
        # http://stackoverflow.com/questions/18383471/django-bulk-create-function-example
        DataSetSample.objects.bulk_create(list_of_sample_objects)

    def generate_and_write_new_stability_file_without_datasheet(self):
        sample_fastq_pairs = self.read_in_mothur_dot_file_creation_output()
        # We will use the stability file in the rest of the mothur qc but also to make the DataSetSamples (below)
        end_index = self.identify_sample_names_without_datasheet()
        new_stability_file = self.generate_new_stability_file_without_datasheet(end_index, sample_fastq_pairs)
        # write out the new stability file
        general.write_list_to_destination(f'{self.temp_working_directory}/stability.files', new_stability_file)

    @staticmethod
    def generate_new_stability_file_without_datasheet(end_index, sample_fastq_pairs):
        new_stability_file = []
        for stability_file_line in sample_fastq_pairs:
            pair_components = stability_file_line.split('\t')
            # I am going to use '[dS]' as a place holder for a dash in the sample names
            # Each line of the stability file is a three column format with the first
            # column being the sample name. The second and third are the full paths of the .fastq files
            # the sample name at the moment is garbage, we will extract the sample name from the
            # first fastq path using the end_index that we determined above

            new_stability_file.append(
                '{}\t{}\t{}'.format(
                    pair_components[1].split('/')[-1][:-end_index].replace('-', '[dS]'),
                    pair_components[1],
                    pair_components[2]))
        return new_stability_file

    def get_num_chars_in_common_with_fastq_names(self):
        i = 1
        while 1:
            list_of_endings = []
            for file in self.list_of_fastq_files_in_wkd:
                list_of_endings.append(file[-i:])
            if len(set(list_of_endings)) > 2:
                break
            else:
                i += 1
                # then this is one i too many and our magic i was i-1
        end_index = i - 1
        return end_index

    def get_sample_names_from_fastq_files_using_index(self, end_index):
        list_of_names_non_unique = []
        for file in self.list_of_fastq_files_in_wkd:
            list_of_names_non_unique.append(file[:-end_index])
        list_of_sample_names = list(set(list_of_names_non_unique))
        if len(list_of_sample_names) != len(self.list_of_fastq_files_in_wkd) / 2:
            warning_str = 'Error in sample name extraction'
            self.exit_and_del_data_set_sample(warning_str)
        self.list_of_samples_names = list_of_sample_names

    def identify_sample_names_without_datasheet(self):
        # I think the simplest way to get sample names is to find what parts are common between all samples
        # well actually 50% of the samples so that we also remove the R1 and R2 parts.
        end_index = self.get_num_chars_in_common_with_fastq_names()
        self.get_sample_names_from_fastq_files_using_index(end_index)

        return end_index
    def generate_stability_file_and_data_set_sample_objects_with_datasheet(self):
        # Create a pandas df from the data_sheet if it was provided
        # allow the data_sheet to be in a .csv format or .xlsx format. This is so that we can store a datasheet
        # in the github repo in a non-binary format
        # The sample_meta_df that is created from the data_sheet should be identical irrespective of whether a .csv
        # or a .xlsx is submitted.
        sample_meta_info_df = self.create_sample_meta_info_dataframe_from_datasheet_path()

        # if we are given a data_sheet then use the sample names given as the DataSetSample object names
        self.list_of_samples_names = sample_meta_info_df.index.values.tolist()

        # we should verify that all of the fastq files listed in the sample_meta_df
        # are indeed found in the directory that we've been given
        self.check_all_fastqs_in_datasheet_exist()

        self.generate_and_write_mothur_batch_file_for_dotfile_creation()

        # noinspection PyPep8
        general.execute_mothur_batch_file_with_piped_stoud_sterr(self.path_to_latest_mothur_batch_file)

        # we will also need to know how to relate the sample names to the fastq files
        # for this we will make a dict of fastq file name to sample
        self.create_fastq_file_to_sample_name_dict()
        # We will use the stability file in the rest of the mothur qc but also to make the DataSetSamples (below)
        self.generate_and_write_new_stability_file_with_data_sheet()

        # TODO we got this far in the refactoring. We are working through this create_data_sub code first
        # we still need to transfer the equivalent into the DataLoading class.

        # Create data_set_sample instances
        self.create_data_set_sample_objects_in_bulk_with_datasheet()


    def create_data_set_sample_objects_in_bulk_with_datasheet(self):
        list_of_data_set_sample_objects = []
        sys.stdout.write('\nCreating data_set_sample objects\n')
        for sampleName in self.list_of_samples_names:
            print('\rCreating data_set_sample {}'.format(sampleName))
            # Create the data_set_sample objects in bulk.
            # The cladal_seq_totals property of the data_set_sample object keeps track of the seq totals for the
            # sample divided by clade. This is used in the output to keep track of sequences that are not
            # included in cladeCollections
            empty_cladal_seq_totals = json.dumps([0 for _ in self.clade_list])

            dss = DataSetSample(name=sampleName, data_submission_from=self.dataset_object,
                                cladal_seq_totals=empty_cladal_seq_totals,
                                sample_type=self.sample_meta_info_df.loc[sampleName, 'sample_type'],
                                host_phylum=self.sample_meta_info_df.loc[sampleName, 'host_phylum'],
                                host_class=self.sample_meta_info_df.loc[sampleName, 'host_class'],
                                host_order=self.sample_meta_info_df.loc[sampleName, 'host_order'],
                                host_family=self.sample_meta_info_df.loc[sampleName, 'host_family'],
                                host_genus=self.sample_meta_info_df.loc[sampleName, 'host_genus'],
                                host_species=self.sample_meta_info_df.loc[sampleName, 'host_species'],
                                collection_latitude=self.sample_meta_info_df.loc[sampleName, 'collection_latitude'],
                                collection_longitude=self.sample_meta_info_df.loc[sampleName, 'collection_longitude'],
                                collection_date=self.sample_meta_info_df.loc[sampleName, 'collection_date'],
                                collection_depth=self.sample_meta_info_df.loc[sampleName, 'collection_depth']
                                )
            list_of_data_set_sample_objects.append(dss)
        # http://stackoverflow.com/questions/18383471/django-bulk-create-function-example
        DataSetSample.objects.bulk_create(list_of_data_set_sample_objects)


    def generate_and_write_new_stability_file_with_data_sheet(self):
        # Convert the group names in the stability.files so that the dashes are converted to '[ds]',
        # So for the mothur we have '[ds]'s. But for all else we convert these '[ds]'s to dashes
        sample_fastq_pairs = self.read_in_mothur_dot_file_creation_output()

        new_stability_file = self.generate_new_stability_file_with_data_sheet(sample_fastq_pairs)

        general.write_list_to_destination(
            r'{0}/stability.files'.format(self.temp_working_directory),
            new_stability_file
        )


    def generate_new_stability_file_with_data_sheet(self, sample_fastq_pairs):
        new_stability_file = []
        for stability_file_line in sample_fastq_pairs:
            pair_components = stability_file_line.split('\t')
            # I am going to use '[dS]' as a place holder for a dash in the sample names
            # Each line of the stability file is a three column format with the first
            # column being the sample name. The second and third are the full paths of the .fastq files
            # the sample name at the moment is garbage, we will identify the sample name from the
            # first fastq path using the fastq_file_to_sample_name_dict
            new_stability_file.append(
                '{}\t{}\t{}'.format(
                    self.fastq_file_to_sample_name_dict[pair_components[1].split('/')[-1]].replace('-', '[dS]'),
                    pair_components[1],
                    pair_components[2]))
        return new_stability_file

    def read_in_mothur_dot_file_creation_output(self):
        return general.read_defined_file_to_list(f'{self.temp_working_directory}/stability.files')

    def create_fastq_file_to_sample_name_dict(self):
        fastq_file_to_sample_name_dict = {}
        for sample_index in self.sample_meta_info_df.index.values.tolist():
            fastq_file_to_sample_name_dict[self.sample_meta_info_df.loc[sample_index, 'fastq_fwd_file_name']] = sample_index
            fastq_file_to_sample_name_dict[self.sample_meta_info_df.loc[sample_index, 'fastq_rev_file_name']] = sample_index
        self.fastq_file_to_sample_name_dict = fastq_file_to_sample_name_dict

    def execute_mothur_batch_file_with_piped_stoud_sterr(path_to_mothur_batch_file):
        subprocess.run(['mothur', path_to_mothur_batch_file],
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)

    def generate_and_write_mothur_batch_file_for_dotfile_creation(self):
        self.check_if_fastqs_are_gz_compressed()
        mothur_batch_file_as_list = self.generate_mothur_batch_file_for_dotfile_creation_as_list()
        self.path_to_latest_mothur_batch_file = f'{self.temp_working_directory}/mothur_batch_file_makeFile'
        general.write_list_to_destination(self.path_to_latest_mothur_batch_file, mothur_batch_file_as_list)

    def generate_mothur_batch_file_for_dotfile_creation_as_list(self):
        if self.fastqs_are_gz_compressed:
            mothur_batch_file = [
                r'set.dir(input={0})'.format(self.temp_working_directory),
                r'set.dir(output={0})'.format(self.temp_working_directory),
                r'make.file(inputdir={0}, type=gz, numcols=3)'.format(self.temp_working_directory)
            ]
        else:
            mothur_batch_file = [
                r'set.dir(input={0})'.format(self.temp_working_directory),
                r'set.dir(output={0})'.format(self.temp_working_directory),
                r'make.file(inputdir={0}, type=fastq, numcols=3)'.format(self.temp_working_directory)
            ]
        return mothur_batch_file

    def check_if_fastqs_are_gz_compressed(self):
        if self.list_of_fastq_files_in_wkd[0].endswith('fastq.gz'):
            self.fastqs_are_gz_compressed = True
        elif self.list_of_fastq_files_in_wkd[0].endswith('fastq'):
            self.fastqs_are_gz_compressed = False
        else:
            warning_str = f'Unrecognised format of sequecing file: {self.list_of_fastq_files_in_wkd[0]}'
            self.exit_and_del_data_set_sample(warning_str)

    def check_all_fastqs_in_datasheet_exist(self):
        self.list_of_fastq_files_in_wkd = [_ for _ in os.listdir(self.temp_working_directory) if 'fastq' in _]
        list_of_meta_gz_files = self.get_list_of_fastq_file_names_that_should_be_in_directory()
        self.if_fastq_files_missing_sys_exit(list_of_meta_gz_files)

    def if_fastq_files_missing_sys_exit(self, list_of_meta_gz_files):
        for fastq in list_of_meta_gz_files:
            if fastq not in self.list_of_fastq_files_in_wkd:
                warning_str = f'{fastq} not found'
                self.exit_and_del_data_set_sample(warning_str)

    def exit_and_del_data_set_sample(self, warning_str):
        self.dataset_object.delete()
        sys.exit(warning_str)

    def get_list_of_fastq_file_names_that_should_be_in_directory(sample_meta_info_df):
        list_of_meta_gz_files = []
        list_of_meta_gz_files.extend(sample_meta_info_df['fastq_fwd_file_name'].values.tolist())
        list_of_meta_gz_files.extend(sample_meta_info_df['fastq_rev_file_name'].values.tolist())
        return list_of_meta_gz_files

    def create_sample_meta_info_dataframe_from_datasheet_path(self):
        if self.datasheet_path.endswith('.xlsx'):
            self.sample_meta_info_df = pd.read_excel(io=self.datasheet_path, header=0, index_col=0, usecols='A:N', skiprows=[0])
        elif self.datasheet_path.endswith('.csv'):
            self.sample_meta_info_df = pd.read_csv(filepath_or_buffer=self.datasheet_path, header=0, index_col=0, skiprows=[0])
        else:
            sys.exit('Data sheet: {} is in an unrecognised format. '
                     'Please ensure that it is either in .xlsx or .csv format.')


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

    def set_temp_working_directory(self):
        # working directory will be housed in a temp folder within the directory in which the sequencing data
        # is currently housed
        if '.' in self.user_input_path.split('/')[-1]:
            # then this path points to a file rather than a directory and we should pass through the path only
            self.temp_working_directory = os.path.abspath(
                '{}/tempData/{}'.format(os.path.dirname(self.user_input_path), self.dataset_object.id))
        else:
            # then we assume that we are pointing to a directory and we can directly use that to make the wkd
            self.temp_working_directory = os.path.abspath('{}/tempData/{}'.format(self.user_input_path, self.dataset_object.id))
        self.dataset_object.working_directory = self.temp_working_directory
        self.dataset_object.save()


    def create_temp_wkd(self):
        # if the directory already exists remove it and start from scratch
        if os.path.exists(self.temp_working_directory):
            shutil.rmtree(self.temp_working_directory)
        os.makedirs(self.temp_working_directory)