import os
import pickle
# ####### Setup Django DB and Models ########
# Ensure settings are read
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
import subprocess
import sys
# Your application specific imports
from dbApp.models import DataSet, DataAnalysis
# ###########################################


def delete_data_set(uid):
    DataSet.objects.get(id=uid).delete()


def delete_data_analysis(uid):
    DataAnalysis.objects.get(id=uid).delete()


def write_list_to_destination(destination, list_to_write):
    # print('Writing list to ' + destination)
    try:
        os.makedirs(os.path.dirname(destination))
    except FileExistsError:
        pass

    with open(destination, mode='w') as writer:
        i = 0
        while i < len(list_to_write):
            if i != len(list_to_write) - 1:
                writer.write(list_to_write[i] + '\n')
            elif i == len(list_to_write) - 1:
                writer.write(list_to_write[i])
            i += 1


def read_defined_file_to_list(filename):
    with open(filename, mode='r') as reader:
        temp_list = [line.rstrip() for line in reader]
    return temp_list


def convert_interleaved_to_sequencial_fasta(fasta_in):
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

def execute_mothur_batch_file_with_piped_stoud_sterr(path_to_mothur_batch_file, mothur_exec_str='mothur'):
    return subprocess.run(
        [mothur_exec_str, path_to_mothur_batch_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )



def read_byte_object_from_defined_directory(directory):
    f = open(directory, 'rb')
    return pickle.load(f)


def write_byte_object_to_defined_directory(directory, byte_object):
    f = open(directory, 'wb+')
    pickle.dump(byte_object, f)


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

def make_new_blast_db(input_fasta_to_make_db_from, db_title, db_type='nucl',makeblastdb_exec_str='makeblastdb', pipe_stdout_sterr=True):
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
