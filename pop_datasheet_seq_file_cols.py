#!/usr/bin/env python3.6
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(
            description='populate the sequencing fwd and rev fastq.gz files in the symportal datasheet\n'
                        'by running this script providing the directory fo the sequencing files\n'
                        'and the datasheet that already has the sample names populated.\n'
                        'It relies on the datasample names being present in the seq file names.')
parser.add_argument('--datasheet_path', help='Full path to the datasheet', required=True)
parser.add_argument('--seq_dir', help='The directory in which the sequencing files are contained', required=True)
# custom_args_list = [
#         '--datasheet_path', 'XXX',
#         '--seq_dir', 'XXX'
#                         ]
# args = parser.parse_args(custom_args_list)
args = parser.parse_args()

class PopDS:
    def __init__(self, ds_path, seq_dir):
        if not os.path.exists(ds_path):
            raise RuntimeError(f'{ds_path} does not exist')
        self.ds_path = ds_path
        if ds_path.endswith('.xlsx'):
            self.ds_df = pd.read_excel(io=ds_path, header=0, skiprows=[0])
            self.new_ds_out_path = ds_path.replace('.xlsx', '_pop.csv')
        elif ds_path.endswith('.csv'):
            self.ds_df = pd.read_csv(ds_path, skiprows=[0], header=0)
            self.new_ds_out_path = ds_path.replace('.csv', '_pop.csv')
        self.sample_names_ordered_by_len = sorted(self.ds_df.iloc[:, 0].values.tolist(), key=lambda x: len(x), reverse=True)
        self.ds_df.set_index('sample_name', drop=True, inplace=True)
        self.ds_df['fastq_fwd_file_name'] = self.ds_df['fastq_fwd_file_name'].astype(str)
        self.ds_df['fastq_rev_file_name'] = self.ds_df['fastq_rev_file_name'].astype(str)
        self.file_names = self._get_file_names(seq_dir)

    def _get_file_names(self, dir_to_walk):
        return [_ for _ in os.listdir(dir_to_walk) if _.endswith('.fastq') or _.endswith('.fastq.gz')]

    def pop_ds(self):
        """Go in order of the length of the sample names.
        In this way we should avoid the problem of having sample names that are subsets of other samples.
        E.g. Egypt_1 and Eygpyt_13. We will also need to actively remove names from a list that have already had
        their corresponding sequencing files found."""
        for i, sample in enumerate(self.sample_names_ordered_by_len):
            count = 0
            files_to_remove = []
            for file in self.file_names:
                if file.startswith(sample.replace('_', '-').rstrip().lstrip()) or file.startswith(sample.rstrip().lstrip()):
                    if 'R1' in file:
                        self.ds_df.at[sample, 'fastq_fwd_file_name'] = file
                        files_to_remove.append(file)
                        count += 1
                    elif 'R2' in file:
                        self.ds_df.at[sample, 'fastq_rev_file_name'] = file
                        files_to_remove.append(file)
                        count += 1
                    else:
                        print(f'Dodgy file name {file}')
                        raise RuntimeError
            if count != 2:
                print(f'count was {count} for {sample}')
                raise RuntimeError
            else:
                for f_to_d in files_to_remove:
                    self.file_names.remove(f_to_d)
                    
        # now write the populated df back out
        with open(self.new_ds_out_path, 'w') as f:
            self.ds_df.to_csv(f)

pop = PopDS(ds_path = args.datasheet_path, seq_dir= args.seq_dir)
pop.pop_ds()
