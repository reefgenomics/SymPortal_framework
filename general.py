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
    elif method == 'unifrac':  # These are the alignment settings specifically for doing the unifrac alignments
        mafft = local[f'{mafft_exec_string}']
        (mafft[
             '--thread', f'{num_proc}', '--maxiterate', f'{iterations}',
             '--ep', '0', '--genafpair', input_path] > output_path)()
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


def chunks(l, n=500):
    """Yield successive n-sized chunks from l.
    Modified to explicitly cast to list to cover the case that a set is passed in.
    Also modified to default to 500 which should be suitably below the SQLITE_MAX_VARIABLE_NUMBER
    """
    in_list = list(l)
    for i in range(0, len(in_list), n):
        yield in_list[i:i + n]


def set_seq_colour_dict(ordered_list_of_seqs):
    """Create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
    to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
    If we are only going to have a legend that is cols x rows as shown below, then we should only use
    that many colours in the plotting."""
    temp_colour_dict = {}
    predefined_colour_dict = get_pre_def_colour_dict()

    for seq_name, color_hash in predefined_colour_dict.items():
        if seq_name in ordered_list_of_seqs:
            temp_colour_dict[seq_name] = color_hash

    colour_palette, grey_palette = get_colour_lists()

    remaining_seqs = [seq for seq in ordered_list_of_seqs if seq not in predefined_colour_dict.keys()]

    for i, seq_name in enumerate(remaining_seqs):
        if i < len(colour_palette):
            temp_colour_dict[seq_name] = colour_palette[i]
        else:
            grey_index = i % len(grey_palette)
            temp_colour_dict[seq_name] = grey_palette[grey_index]

    return temp_colour_dict

def make_js_function_to_return_json_file(function_name, json_path=None, json_file_as_str=None):
    temp_js_file_as_list = []
    temp_js_file_as_list.append('function ' + function_name + '(){')
    if json_path:
        temp_js_file_as_list.extend(read_defined_file_to_list(json_path))
    else:
        temp_js_file_as_list.append(json_file_as_str)
    temp_js_file_as_list[1] = 'return ' + temp_js_file_as_list[1]
    temp_js_file_as_list.append('};')
    return temp_js_file_as_list

def make_json_object_array_from_python_dictionary(p_dict):
    json_str = ''
    json_str += '['
    for k, v in p_dict.items():
        json_str += '{'
        json_str += f"\"d_key\":\"{k}\", "
        json_str += f"\"d_value\":\"{v}\""
        json_str += '}, '
    # remove last commar and space
    json_str = json_str[:-2]
    json_str += ']'
    return json_str

def output_js_color_objects_array(output_directory, colour_dict, js_file_name):
    # write out the colour dict to the html file
    json_col_dict_object_array_as_list = make_js_function_to_return_json_file(
        json_file_as_str=make_json_object_array_from_python_dictionary(colour_dict),
        function_name='getColDictObjArr')
    write_list_to_destination(
        destination=os.path.join(output_directory, js_file_name),
        list_to_write=json_col_dict_object_array_as_list)


def get_colour_lists():
    colour_palette = get_colour_list()
    grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
    return colour_palette, grey_palette

def get_pre_def_colour_dict():
    """These are the top 40 most abundnant named sequences. I have hardcoded their color."""
    return {
        'A1'   :"#FFFF00", 'C3'  :"#1CE6FF", 'C15'  :"#FF34FF", 'A1bo':"#FF4A46", 'D1'   :"#008941",
        'C1'   :"#006FA6", 'C27' :"#A30059", 'D4'   :"#FFDBE5", 'C3u' :"#7A4900", 'C42.2':"#0000A6",
        'A1bp' :"#63FFAC", 'C115':"#B79762", 'C1b'  :"#004D43", 'C1d' :"#8FB0FF", 'A1c'  :"#997D87",
        'C66'  :"#5A0007", 'A1j' :"#809693", 'B1'   :"#FEFFE6", 'A1k' :"#1B4400", 'A4'   :"#4FC601",
        'A1h'  :"#3B5DFF", 'C50a':"#4A3B53", 'C39'  :"#FF2F80", 'C3dc':"#61615A", 'D4c'  :"#BA0900",
        'C3z'  :"#6B7900", 'C21' :"#00C2A0", 'C116' :"#FFAA92", 'A1cc':"#FF90C9", 'C72'  :"#B903AA",
        'C15cl':"#D16100", 'C31' :"#DDEFFF", 'C15cw':"#000035", 'A1bv':"#7B4F4B", 'D6'   :"#A1C299",
        'A4m'  :"#300018", 'C42a':"#0AA6D8", 'C15cr':"#013349", 'C50l':"#00846F", 'C42g' :"#372101" }


def get_colour_list():
    colour_list = [
        "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
        "#0CBD66",
        "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
        "#BEC459",
        "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9",
        "#FF913F",
        "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7",
        "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
        "#201625",
        "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
        "#CB7E98",
        "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489",
        "#806C66",
        "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5",
        "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58", "#7A7BFF", "#D68E01",
        "#353339",
        "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A",
        "#001325",
        "#02525F", "#0AA3F7", "#E98176", "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75",
        "#8D8546",
        "#9695C5", "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
        "#00005F",
        "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058",
        "#A45B02",
        "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
        "#F4D749",
        "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE",
        "#C6DC99",
        "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
        "#C6005A",
        "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183", "#AA5199", "#B5D6C3",
        "#A38469",
        "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433", "#789EC9", "#6D80BA", "#953F00", "#5EFF03",
        "#E4FFFC",
        "#1BE177", "#BCB1E5", "#76912F", "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213",
        "#A76F42",
        "#89412E", "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
        "#BDC9D2",
        "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71",
        "#868E7E",
        "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C",
        "#00B57F",
        "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
    return colour_list