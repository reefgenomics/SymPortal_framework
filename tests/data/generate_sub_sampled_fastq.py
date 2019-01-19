from plumbum import local
import os

def generate_subsampled_fastqs():
    """This script will create fastaq files that are subsampled versions of the original Ed data set.
    I will use these as part of the testing framework for symportal. I will aim to subsample down to 500 sequences.
    I'm hoping that this will be a sufficiently small number of sequences to allow for the files to be stored
    easily as part of the repository. I will use the seqtk package to do the subsampling as I have tested this
    on a small scale and it appears to work well. Note that seqtk must be in your path for this to work.
    """

    # get a list of the files that need to be converted
    list_of_file_paths_to_subsample = return_list_of_file_paths_to_be_convereted(
        input_directory=os.path.abspath(os.path.dirname(__file__))
    )

    # setup seqtk
    seqtk_local_var = local["seqtk"]

    subsample_fastqs(list_of_file_paths_to_subsample = list_of_file_paths_to_subsample,
                     sub_sampling_depth = 500,
                     random_seed_int = 100,
                     plumbum_seqtk_variable = seqtk_local_var,
                     name_split_char = '_',
                     name_split_int = 1)

def return_list_of_file_paths_to_be_convereted(input_directory):
    """ This function will look in our current directory (the directory of this file) and grab a list
    of the .gz files. It is important that only the files for subsampling are in this directory. It will return a list
    that contains the full paths for each of the files that are to be subsampled.
    """

    list_of_files_in_directory = return_list_of_file_names_in_directory(input_directory)

    list_of_file_names_in_directory_that_are_gz = [
        file_name for file_name in list_of_files_in_directory if '.gz' in file_name
    ]

    return [
        os.path.join(input_directory, file_name) for file_name in list_of_file_names_in_directory_that_are_gz
    ]


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


def subsample_fastqs(list_of_file_paths_to_subsample,
                     sub_sampling_depth, random_seed_int,
                     plumbum_seqtk_variable, name_split_char,
                     name_split_int):
    """
    This function will perform the actual subsampling of the fastq files that have been provided in the
    list_of_file_paths_to_subsample using the seqtk application (https://github.com/lh3/seqtk). The seqtk
    application by default takes in either fastq or fastq.gz files but always writes out fastq files. As such,
    we will need to recompress the fastq files once they have been generated. This will return a list of the full path
    names for each of the fastq files that has been created so that the files can be gzipped if required.
    :param list_of_file_paths_to_subsample: the list that contains the full path of the files to be subsampled
    :param sub_sampling_depth: the number of sequences to be pulled out of the fastq
    :param random_seed_int: the seed integer that will be used to generate the subsampling. Note, that the
    same seed must be used for each pair of subsamplings if the extracted sequences are to match and therefore allow
    contig creation.
    :param plumbum_seqtk_variable: the variable that has been assigned to the plumbum local instance of seqtk
    :param name_split_char: when we generate the name for the newly created fastq.gz we will generate this from the
    file name that we are subsampling. We will do this by doing a str.split and using a first n number of items
    returned from this split as the base of the name. The character by which we split is given by the name_split_char
    parameter. The number of items from the split is given by the name_split_int parameter.
    :param name_split_int: see name_split_char.
    :return: a list of the full path names for each of the newly created subsampled fastq files
    """

    list_of_paths_for_newly_created_fastq_files = []

    for fastq_to_subsample_path in list_of_file_paths_to_subsample:
        fastq_out_full_path = derive_and_return_full_path_of_new_fastq_file(
            fastq_to_subsample_path, name_split_char, name_split_int)
        list_of_paths_for_newly_created_fastq_files.append(fastq_out_full_path)

        # now run seqtk
        (plumbum_seqtk_variable[
             'sample', '-s', random_seed_int, fastq_to_subsample_path, sub_sampling_depth
         ] > fastq_out_full_path)()

    return  list_of_paths_for_newly_created_fastq_files


def derive_and_return_full_path_of_new_fastq_file(fastq_to_subsample_path, name_split_char, name_split_int):
    """
    Generates the path of the new fastq file that we are going to create by subsampling the fastq_to_subsample_path
    file
    :param fastq_to_subsample_path: The path of the fastq file that will be sampled
    :param name_split_char: when we generate the name for the newly created fastq.gz we will generate this from the
    file name that we are subsampling. We will do this by doing a str.split and using a first n number of items
    returned from this split as the base of the name. The character by which we split is given by the name_split_char
    parameter. The number of items from the split is given by the name_split_int parameter.
    :param name_split_int: see name_split_char.
    :return: the full path of the fastq to be created through subsampling
    """
    directory_of_input_fastq = os.path.dirname(fastq_to_subsample_path)
    fastq_to_subsample_file_name = fastq_to_subsample_path.split('/')[-1]
    fastq_out_name_base = '_'.join(fastq_to_subsample_file_name.split(name_split_char)[:name_split_int])
    fastq_out_name = '{}_subsampled.fastq'.format(fastq_out_name_base)
    fastq_out_full_path = '{}/{}'.format(directory_of_input_fastq, fastq_out_name)
    return fastq_out_full_path

generate_subsampled_fastqs()