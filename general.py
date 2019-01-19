import os
import pickle
def writeListToDestination(destination, listToWrite):
    # print('Writing list to ' + destination)
    try:
        os.makedirs(os.path.dirname(destination))
    except FileExistsError:
        pass

    with open(destination, mode='w') as writer:
        i = 0
        while i < len(listToWrite):
            if i != len(listToWrite) - 1:
                writer.write(listToWrite[i] + '\n')
            elif i == len(listToWrite) - 1:
                writer.write(listToWrite[i])
            i += 1


def readDefinedFileToList(filename):
    temp_list = []
    with open(filename, mode='r') as reader:
        temp_list = [line.rstrip() for line in reader]
    return temp_list


def convert_interleaved_to_sequencial_fasta(fasta_in):
    fasta_dict = {}
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


def readByteObjectFromDefinedDirectory(directory):
    f = open(directory, 'rb')
    return pickle.load(f)

def writeByteObjectToDefinedDirectory(directory,object):
    f = open(directory , 'wb+')
    pickle.dump(object, f)

def createNoSpaceFastaFile(fastaList):
    tempList = []
    i = 0
    while i < len(fastaList):
        tempList.extend([fastaList[i].split('\t')[0], fastaList[i+1]])
        i += 2
    return tempList

def createDictFromFasta(fastaList):
    temporary_dictionary = {}
    i = 0
    while i < len(fastaList):
        sequence = fastaList[i][1:]
        temporary_dictionary[sequence] = fastaList[i+1]
        i += 2
    return temporary_dictionary

def createNewFile(pathtofile):
    try:
        os.makedirs(os.path.dirname(pathtofile))
    except FileExistsError:
        pass

    open(pathtofile, mode='w')

def writeLinesToFile(path_to_file, listoflines):
    with open(path_to_file, mode='a') as f:
        for line in listoflines:
            f.write(line + '\n')

def checkIfFileEmpty(filepath):
    count = 0
    with open(filepath, mode='r') as reader:
        for line in reader:
            if count != 0:
                return False  # not empty
            count += 1
    return True