#!/usr/bin/env python3.6

import argparse
import os
from general import return_list_of_file_paths_in_directory
from plumbum import local
class DataSetSubSampler:
    def __init__(self):
        self.seqtk = local["seqtk"]
        self.parser = argparse.ArgumentParser(
            description='Intragenomic analysis of the ITS2 region of the nrDNA',
            epilog='For support email: benjamin.hume@uni-konstanz.de')
        self.parser.add_argument(
            '--input_dir', type=str, required=True,
            help='The full path to the directory containing the fastq or fastq.gz files to be subsampled')
        self.parser.add_argument(
            '--output_dir', type=str, required=True,
            help='The output path to the directory where the subsampled fastq files will be output')
        self.parser.add_argument('--depth', help='subsampling depth.', required=True)
        self.parser.add_argument('-s', '--suffix', help='The string that will be appened between the base of '
                                                        'the original file name and its extension. E.g. if original '
                                                        'was seq_file_one.fastq.gz and \'_sub\' is passed to --sufix, '
                                                        'then the output file name would be seq_file_one_sub.fastq.gz. '
                                                        'Default = \'_sub\'', default='_sub')
        self.parser.add_argument('-f', '--force_overwrite', help='Overwrite existing files.', action='store_true')
        # self.parser.add_argument('--gz', action='store_true', default='inferred')
        self.args = self.parser.parse_args()
        os.makedirs(self.args.output_dir, exist_ok=True)
        print(f'output_dir is set to:{self.args.output_dir}')
        print(f'input_dir is set to:{self.args.input_dir}')
        self.output_dir = os.path.abspath(self.args.output_dir)
        self.input_dir = os.path.abspath(self.args.input_dir)


    def sub_sample_data(self):
        for file_path in return_list_of_file_paths_in_directory(self.input_dir):
            extension = self._determine_extension_of_file_name(file_path)
            # nb seqtk outputs non-compressed only.
            out_extension = extension.replace('.gz', '')
            out_name = os.path.basename(file_path).replace(extension, f'{self.args.suffix}{out_extension}')
            out_path = os.path.join(self.output_dir, out_name)
            if not self.args.force_overwrite and os.path.exists(out_path):
                print(f'Output file {out_path} already exists. File will not be overwritten. '
                      f'Pass -f or --force_overwrite to overwrite.')
            else:
                print(f'Subsampling {file_path} to {self.args.depth}')
                (self.seqtk['sample', '-s100', file_path, self.args.depth] > out_path)()
        print('subsampling complete')


    def _determine_extension_of_file_name(self, file_path):
        if file_path.endswith('.fastq.gz'):
            return '.fastq.gz'
        elif file_path.endswith('.fastq'):
            return '.fastq'
        elif file_path.endswith('.fq.gz'):
            return '.fq'
        else:
            return '.fq'

if __name__ == "__main__":
    dsss = DataSetSubSampler()
    dsss.sub_sample_data()
