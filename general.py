import os
import pickle
import subprocess
import sys
import pandas as pd
from plumbum import local
import numpy as np

def return_list_of_file_names_in_directory(directory_to_list):
    """
    return a list that contains the filenames found in the specified directory
    :param directory_to_list: the directory that the file names should be returned from
    :return: list of strings that are the file names found in the directory_to_list
    """
    list_of_file_names_in_directory = []
    for (dirpath, dirnames, filenames) in os.walk(directory_to_list):
        list_of_file_names_in_directory.extend(filenames)
    return list_of_file_names_in_directory


def return_list_of_file_paths_in_directory(directory_to_list):
    """
    return a list that contains the full paths of each of the files found in the specified directory
    :param directory_to_list: the directory that the file paths should be returned from
    :return: list of strings that are the file paths found in the directory_to_list
    """
    list_of_file_paths_in_directory = []
    for (dirpath, dirnames, filenames) in os.walk(directory_to_list):
        list_of_file_paths_in_directory.extend([os.path.join(directory_to_list, file_name) for file_name in filenames])
        return list_of_file_paths_in_directory


def write_list_to_destination(destination, list_to_write):
    with open(destination, mode='w') as writer:
        for line in list_to_write:
            writer.write(f'{line}\n')


def write_byte_object_to_defined_directory(directory, byte_object):
    f = open(directory, 'wb+')
    pickle.dump(byte_object, f)


def read_defined_file_to_list(filename):
    with open(filename, mode='r') as reader:
        return [line.rstrip() for line in reader]


def read_defined_file_to_generator(filename):
    with open(filename, mode='r') as reader:
        return (line.rstrip() for line in reader)


def read_byte_object_from_defined_directory(directory):
    f = open(directory, 'rb')
    return pickle.load(f)


def convert_interleaved_to_sequencial_fasta_first_line_removal(fasta_in):
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
            index = i % num_seqs
            list_seq_sequences[index] += ''.join(fasta_cropped[i].split()[1:])

    out_fasta = []
    for name, seq in zip(list_seq_names, list_seq_sequences):
        out_fasta.extend(['>{}'.format(name), seq])

    return out_fasta


def convert_interleaved_to_sequencial_fasta(fasta_as_list):
    new_fasta = []
    temp_seq_string_list = []
    for i, fasta_line in enumerate(fasta_as_list):
        if fasta_line.startswith('>'):
            if temp_seq_string_list:
                new_fasta.append(''.join(temp_seq_string_list))
                temp_seq_string_list = []
                new_fasta.append(fasta_line)
            else:
                new_fasta.append(fasta_line)
        else:
            temp_seq_string_list.append(fasta_line)
    new_fasta.append(''.join(temp_seq_string_list))
    return new_fasta


def execute_mothur_batch_file_with_piped_stoud_sterr(path_to_mothur_batch_file, mothur_exec_str='mothur'):
    return subprocess.run(
        [mothur_exec_str, path_to_mothur_batch_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )


def create_no_space_fasta_file(fasta_list):
    temp_list = []
    i = 0
    while i < len(fasta_list):
        temp_list.extend([fasta_list[i].split('\t')[0], fasta_list[i + 1]])
        i += 2
    return temp_list


def create_dict_from_fasta(fasta_list=None, fasta_path=None):
    if fasta_list is None and fasta_path is None:
        sys.exit('Please provide either a fasta_as_list OR a fasta_path as arguments to create_dict_from_fasta')
    elif fasta_list and fasta_path:
        sys.exit('Please provide either a fasta_as_list OR a fasta_path as arguments to create_dict_from_fasta')
    else:
        if fasta_list:
            temporary_dictionary = {}
            i = 0
            while i < len(fasta_list):
                sequence = fasta_list[i][1:]
                temporary_dictionary[sequence] = fasta_list[i + 1]
                i += 2
            return temporary_dictionary
        if fasta_path:
            fasta_file_as_list = read_defined_file_to_list(fasta_path)
            temporary_dictionary = {}
            i = 0
            while i < len(fasta_file_as_list):
                sequence = fasta_file_as_list[i][1:]
                temporary_dictionary[sequence] = fasta_file_as_list[i + 1]
                i += 2
            return temporary_dictionary


def create_seq_name_to_abundance_dict_from_name_file(name_file_list = None, name_file_path = None):
    if name_file_list is None and name_file_path is None:
        sys.exit(
            'Please provide either a name_file_list OR a name_file_path as '
            'arguments to create_seq_name_to_abundance_dict_from_name_file')
    elif name_file_list and name_file_path:
        sys.exit(
            'Please provide either a name_file_list OR a name_file_path '
            'as arguments to create_seq_name_to_abundance_dict_from_name_file')
    else:
        if name_file_list:
            temporary_dictionary = {}
            for i in range(len(name_file_list)):
                temporary_dictionary[name_file_list[i].split('\t')[0]] = len(
                    name_file_list[i].split('\t')[1].split(','))
            return temporary_dictionary
        if name_file_path:
            name_file_as_list = read_defined_file_to_list(name_file_path)
            temporary_dictionary = {}
            for i in range(len(name_file_as_list)):
                temporary_dictionary[name_file_as_list[i].split('\t')[0]] = len(
                    name_file_as_list[i].split('\t')[1].split(','))
            return temporary_dictionary


def make_new_blast_db(
        input_fasta_to_make_db_from, db_title, db_type='nucl',
        makeblastdb_exec_str='makeblastdb', pipe_stdout_sterr=True):
    if pipe_stdout_sterr:
        completed_process = subprocess.run(
            [makeblastdb_exec_str, '-in', input_fasta_to_make_db_from, '-dbtype', db_type, '-title', db_title],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    else:
        completed_process = subprocess.run(
            [makeblastdb_exec_str, '-in', input_fasta_to_make_db_from, '-dbtype', db_type, '-title', db_title]
        )
    return completed_process


def decode_utf8_binary_to_list(bin_to_decode):
    return bin_to_decode.decode('ISO-8859-1').split('\n')


def mafft_align_fasta(input_path, output_path, method='auto', mafft_exec_string='mafft', num_proc=1, iterations=1000):
    # http://plumbum.readthedocs.io/en/latest/local_commands.html#pipelining
    print(f'Aligning {input_path}')
    if method == 'auto':
        mafft = local[f'{mafft_exec_string}']
        (mafft['--auto', '--thread', f'{num_proc}', input_path] > output_path)()
    elif method == 'linsi':
        mafft = local[f'{mafft_exec_string}']
        (mafft['--localpair', '--maxiterate', f'{iterations}', '--thread', f'{num_proc}', input_path] > output_path)()
    print(f'Writing to {output_path}')


def remove_gaps_from_fasta(fasta_as_list):
    gapless_fasta = []
    for fasta_line in fasta_as_list:
        if fasta_line.startswith('>'):
            gapless_fasta.append(fasta_line)
        else:
            gapless_fasta.append(fasta_line.replace('-', ''))
    return gapless_fasta


def fasta_to_pandas_df(fasta_as_list):
    temp_df = pd.DataFrame([list(line) for line in fasta_as_list if not line.startswith('>')])
    seq_names = [line[1:] for line in fasta_as_list if line.startswith('>')]
    temp_df.index=seq_names
    return temp_df


def pandas_df_to_fasta(cropped_fasta_df):
    temp_fasta = []
    for ind in cropped_fasta_df.index.tolist():
        temp_fasta.extend(['>{}'.format(ind), ''.join(list(cropped_fasta_df.loc[ind]))])
    return temp_fasta


def combine_two_fasta_files(path_one, path_two, path_for_combined):
    one_file_one = read_defined_file_to_list(path_one)
    one_file_two = read_defined_file_to_list(path_two)
    one_file_one.extend(one_file_two)
    write_list_to_destination(path_for_combined, one_file_one)


def return_list_of_directory_names_in_directory(directory_to_list):
    """
        return a list that contains the directory names found in the specified directory
        :param directory_to_list: the directory that the directory names should be returned from
        :return: list of strings that are the directory names found in the directory_to_list
        """
    list_of_directory_names_in_directory = []
    for (dirpath, dirnames, filenames) in os.walk(directory_to_list):
        list_of_directory_names_in_directory.extend(dirnames)
        return list_of_directory_names_in_directory


def return_list_of_directory_paths_in_directory(directory_to_list):
    """
        return a list that contains the full paths of each of the directories found in the specified directory
        :param directory_to_list: the directory that the directory paths should be returned from
        :return: list of strings that are the directory paths found in the directory_to_list
        """
    list_of_directory_paths_in_directory = []
    for (dirpath, dirnames, filenames) in os.walk(directory_to_list):
        list_of_directory_paths_in_directory.extend([
            os.path.join(directory_to_list, dir_name) for dir_name in dirnames])
        return list_of_directory_paths_in_directory


def sqrt_transform_abundance_df(df):
    new_df = df.apply(np.sqrt)
    new_df['sum'] = new_df.sum(axis=1)
    new_df = new_df.iloc[:, :-1].div(new_df['sum'], axis=0)
    return new_df


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]