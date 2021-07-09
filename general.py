import os
import pickle
import subprocess
import sys
import pandas as pd
from plumbum import local
import json
import numpy as np
import random
import re

class ThreadSafeGeneral:
    def __init__(self):
        return
    
    @staticmethod
    def write_list_to_destination(destination, list_to_write):
        with open(destination, mode='w') as writer:
            for line in list_to_write:
                writer.write(f'{line}\n')

    @staticmethod
    def read_defined_file_to_list(filename):
        with open(filename, mode='r') as reader:
            return [line.rstrip() for line in reader]

    def create_dict_from_fasta(self, fasta_list=None, fasta_path=None):
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
                fasta_file_as_list = self.read_defined_file_to_list(fasta_path)
                temporary_dictionary = {}
                i = 0
                while i < len(fasta_file_as_list):
                    sequence = fasta_file_as_list[i][1:]
                    temporary_dictionary[sequence] = fasta_file_as_list[i + 1]
                    i += 2
                return temporary_dictionary

    @staticmethod
    def decode_utf8_binary_to_list(bin_to_decode):
        return bin_to_decode.decode('ISO-8859-1').split('\n')

    def combine_two_fasta_files(self, path_one, path_two, path_for_combined):
        one_file_one = self.read_defined_file_to_list(path_one)
        one_file_two = self.read_defined_file_to_list(path_two)
        one_file_one.extend(one_file_two)
        self.write_list_to_destination(path_for_combined, one_file_one)

    def create_seq_name_to_abundance_dict_from_name_file(self, name_file_list = None, name_file_path = None):
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
                name_file_as_list = self.read_defined_file_to_list(name_file_path)
                temporary_dictionary = {}
                for i in range(len(name_file_as_list)):
                    temporary_dictionary[name_file_as_list[i].split('\t')[0]] = len(
                        name_file_as_list[i].split('\t')[1].split(','))
                return temporary_dictionary

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def write_byte_object_to_defined_directory(directory, byte_object):
        f = open(directory, 'wb+')
        pickle.dump(byte_object, f)

    @staticmethod
    def read_defined_file_to_generator(filename):
        with open(filename, mode='r') as reader:
            return (line.rstrip() for line in reader)

    @staticmethod
    def read_byte_object_from_defined_directory(directory):
        f = open(directory, 'rb')
        return pickle.load(f)

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def execute_mothur_batch_file_with_piped_stoud_sterr(path_to_mothur_batch_file, mothur_exec_str='mothur'):
        return subprocess.run(
            [mothur_exec_str, path_to_mothur_batch_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

    @staticmethod
    def create_no_space_fasta_file(fasta_list):
        temp_list = []
        i = 0
        while i < len(fasta_list):
            temp_list.extend([fasta_list[i].split('\t')[0], fasta_list[i + 1]])
            i += 2
        return temp_list

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def remove_gaps_from_fasta(fasta_as_list):
        gapless_fasta = []
        for fasta_line in fasta_as_list:
            if fasta_line.startswith('>'):
                gapless_fasta.append(fasta_line)
            else:
                gapless_fasta.append(fasta_line.replace('-', ''))
        return gapless_fasta

    @staticmethod
    def fasta_to_pandas_df(fasta_as_list):
        temp_df = pd.DataFrame([list(line) for line in fasta_as_list if not line.startswith('>')])
        seq_names = [line[1:] for line in fasta_as_list if line.startswith('>')]
        temp_df.index=seq_names
        return temp_df

    @staticmethod
    def pandas_df_to_fasta(cropped_fasta_df):
        temp_fasta = []
        for ind in cropped_fasta_df.index.tolist():
            temp_fasta.extend(['>{}'.format(ind), ''.join(list(cropped_fasta_df.loc[ind]))])
        return temp_fasta

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def sqrt_transform_abundance_df(df):
        new_df = df.apply(np.sqrt)
        new_df['sum'] = new_df.sum(axis=1)
        new_df = new_df.iloc[:, :-1].div(new_df['sum'], axis=0)
        return new_df

    @staticmethod
    def chunks(l, n=500):
        """Yield successive n-sized chunks from l.
        Modified to explicitly cast to list to cover the case that a set is passed in.
        Also modified to default to 500 which should be suitably below the SQLITE_MAX_VARIABLE_NUMBER
        """
        in_list = list(l)
        for i in range(0, len(in_list), n):
            yield in_list[i:i + n]

    def set_seq_colour_dict(self, ordered_list_of_seqs):
        """Create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
        to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
        If we are only going to have a legend that is cols x rows as shown below, then we should only use
        that many colours in the plotting."""
        temp_colour_dict = {}
        predefined_colour_dict = self.get_pre_def_colour_dict()

        for seq_name, color_hash in predefined_colour_dict.items():
            if seq_name in ordered_list_of_seqs:
                temp_colour_dict[seq_name] = color_hash

        colour_palette, grey_palette = self.get_colour_lists()

        remaining_seqs = [seq for seq in ordered_list_of_seqs if seq not in predefined_colour_dict.keys()]

        for i, seq_name in enumerate(remaining_seqs):
            if i < len(colour_palette):
                temp_colour_dict[seq_name] = colour_palette[i]
            else:
                grey_index = i % len(grey_palette)
                temp_colour_dict[seq_name] = grey_palette[grey_index]

        return temp_colour_dict

    def set_seq_colour_dict_w_reference_c_dict(self, ordered_list_of_seqs, existing_dict):
        """This method has the same purpose as set_seq_colour_dict, but in this case takes into
        account another colour dict that has already been created. The aim here is to create a second colour dict
        that is in complete sync with the first colour dict. If there is an item shared in common between
        the two dictionaries, we can guarantee that they are assigned the same colour"""
        temp_colour_dict = {}
        predefined_colour_dict = self.get_pre_def_colour_dict()

        # first assign the colours of the items that are already in the other dict
        for item, colour in existing_dict.items():
            if item in ordered_list_of_seqs:
                temp_colour_dict[item] = colour

        # now check to see if any of the remaining ordered_list_of_seqs items are in the
        # predefined_colour_dict. If they are then assign them this colour
        items_to_check_still = [item for item in ordered_list_of_seqs if item not in temp_colour_dict]
        for item in items_to_check_still:
            if item in predefined_colour_dict:
                temp_colour_dict[item] = predefined_colour_dict[item]

        # now all of the remaining seqs need to be given any colours that are not already being used by the second
        # dictionary, or otherwise they need to be given a grey colour
        colour_palette, grey_palette = self.get_colour_lists()
        items_to_check_still = [item for item in items_to_check_still if item not in predefined_colour_dict]
        colour_palette = [colour for colour in colour_palette if colour not in existing_dict.values()]

        for i, seq_name in enumerate(items_to_check_still):
            if i < len(colour_palette):
                temp_colour_dict[seq_name] = colour_palette[i]
            else:
                grey_index = i % len(grey_palette)
                temp_colour_dict[seq_name] = grey_palette[grey_index]

        return temp_colour_dict


    def make_js_function_to_return_json_file(self, function_name, json_path=None, json_file_as_str=None):
        temp_js_file_as_list = []
        temp_js_file_as_list.append('function ' + function_name + '(){')
        if json_path:
            temp_js_file_as_list.extend(self.read_defined_file_to_list(json_path))
        else:
            temp_js_file_as_list.append(json_file_as_str)
        temp_js_file_as_list[1] = 'return ' + temp_js_file_as_list[1]
        temp_js_file_as_list.append('};')
        return temp_js_file_as_list

    @staticmethod
    def write_out_js_file_to_return_python_objs_as_js_objs(list_of_func_obj_dicts, js_outpath):
        '''This function writes out a javascript file that will the javascript version of one or more
        python objects. The list_of_func_obj_dicts should ab a list of dictionaries, where each dictionary
        has a key of function_name and python_obj, the values of these keys will be used below.'''
        temp_js_file_as_list = []
        for python_obj_dict in list_of_func_obj_dicts:
            obj_as_json = json.dumps(python_obj_dict['python_obj'])

            temp_js_file_as_list.append('function ' + python_obj_dict['function_name'] + '(){')
            temp_js_file_as_list.append("\treturn " + obj_as_json + ';')
            temp_js_file_as_list.append("}")

        if os.path.isfile(js_outpath):
            # study_data.js already exists so append to it
            with open(js_outpath, 'a') as f:
                for line in temp_js_file_as_list:
                    f.write(f'{line}\n')
        else:
            # study_data.js doesn't exist yet so create and write to it
            with open(js_outpath, 'w') as f:
                for line in temp_js_file_as_list:
                    f.write(f'{line}\n')

    @staticmethod
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


    def output_js_color_objects_array(self, output_directory, colour_dict, js_file_name, function_name):
        # write out the colour dict to the html file
        json_col_dict_object_array_as_list = self.make_js_function_to_return_json_file(
            json_file_as_str=self.make_json_object_array_from_python_dictionary(colour_dict),
            function_name=function_name)
        self.write_list_to_destination(
            destination=os.path.join(output_directory, js_file_name),
            list_to_write=json_col_dict_object_array_as_list)

    @staticmethod
    def json_out_df(path_to_json_file, df_to_json, remove_last_row):
        # manually add the sample_uid to the df to be jsoned as the orient='record' option does not output the index
        temp_df = df_to_json.copy()
        temp_df['sample_uid'] = df_to_json.index.values.tolist()
        cols = temp_df.columns.values.tolist()
        cols_reorder = cols[-1:] + cols[:-1]
        temp_df = temp_df[cols_reorder]
        # remove the assension numbers row at the end
        if remove_last_row:
            temp_df.iloc[:-1].to_json(path_or_buf=path_to_json_file, orient='records')
        else:
            temp_df.to_json(path_or_buf=path_to_json_file, orient='records')


    def get_colour_lists(self):
        colour_palette = self.get_colour_list()
        grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
        return colour_palette, grey_palette

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def create_colour_list(
            sq_dist_cutoff=None, mix_col=None, num_cols=50, time_out_iterations=10000, avoid_black_and_white=True):
        new_colours = []
        min_dist = []
        attempt = 0
        while len(new_colours) < num_cols:
            attempt += 1
            # Check to see if we have run out of iteration attempts to find a colour that fits into the colour space
            if attempt > time_out_iterations:
                sys.exit('Colour generation timed out. We have tried {} iterations of colour generation '
                         'and have not been able to find a colour that fits into your defined colour space.\n'
                         'Please lower the number of colours you are trying to find, '
                         'the minimum distance between them, or both.'.format(attempt))
            if mix_col:
                r = int((random.randint(0, 255) + mix_col[0]) / 2)
                g = int((random.randint(0, 255) + mix_col[1]) / 2)
                b = int((random.randint(0, 255) + mix_col[2]) / 2)
            else:
                r = random.randint(0, 255)
                g = random.randint(0, 255)
                b = random.randint(0, 255)

            # now check to see whether the new colour is within a given distance
            # if the avoids are true also
            good_dist = True
            if sq_dist_cutoff:
                dist_list = []
                for i in range(len(new_colours)):
                    distance = (new_colours[i][0] - r)**2 + (new_colours[i][1] - g)**2 + (new_colours[i][2] - b)**2
                    dist_list.append(distance)
                    if distance < sq_dist_cutoff:
                        good_dist = False
                        break
                # now check against black and white
                d_to_black = (r - 0)**2 + (g - 0)**2 + (b - 0)**2
                d_to_white = (r - 255)**2 + (g - 255)**2 + (b - 255)**2
                if avoid_black_and_white:
                    if d_to_black < sq_dist_cutoff or d_to_white < sq_dist_cutoff:
                        good_dist = False
                if dist_list:
                    min_dist.append(min(dist_list))
            if good_dist:
                new_colours.append((r, g, b))
                attempt = 0

        return new_colours

def chunks(l, n=500):
        """Yield successive n-sized chunks from l.
        Modified to explicitly cast to list to cover the case that a set is passed in.
        Also modified to default to 500 which should be suitably below the SQLITE_MAX_VARIABLE_NUMBER
        """
        in_list = list(l)
        for i in range(0, len(in_list), n):
            yield in_list[i:i + n]

# https://stackoverflow.com/a/3431835/5516420
def hash_bytestr_iter(bytesiter, hasher, ashexstr=False):
    for block in bytesiter:
        hasher.update(block)
    return hasher.hexdigest() if ashexstr else hasher.digest()

def file_as_blockiter(afile, blocksize=65536):
    with afile:
        block = afile.read(blocksize)
        while len(block) > 0:
            yield block
            block = afile.read(blocksize)

def check_lat_lon(lat, lon):
    """
    Takes a dirty lat and lon value and either converts to decimial degrees or raises a run time error.
    Should be able to handle:
        decimal degrees that have a N S W or E character added to it
        degrees decimal minute (with N S W or E)
        degrees minutes seconds (with N S W or E)
    """
    if lat == 'nan' or lon == 'nan':
        raise RuntimeError
    try:
        if np.isnan(lat) or np.isnan(lon):
            raise RuntimeError
    except TypeError:
        pass
    try:
        lat_float = float(lat)
        lon_float = float(lon)
    except ValueError as e:
        # Three options.
        # 1 - we have been given decimal degrees with a hemisphere sign
        # 2 - we have been given degree decimal minutes I.e. 27°38.611'N or N27°38.611'
        # 3 - we have been given degree minutes seconds
        try:
            if chr(176) in lat and chr(176) in lon:
                # Then we are working with either 2 or 3
                if "\"" in lat and "\"" in lon:
                    # Then we are working with degree minutes seconds (DMS)
                    # and we should convert using self.dms2dec
                    lat_float = dms2dec(lat)
                    lon_float = dms2dec(lon)
                elif "\"" not in lat and "\"" not in lon:
                    # Then we are working with degree degree minutes (DDM)
                    # To convert to DD we simply need to remove the NWES characters,
                    # divide the decimal part by 60 and add it to the degrees part
                    if "N" in lat:
                        lat = lat.replace("N", "").replace("'", "")
                        (lat_deg, lat_degmin) = lat.split(chr(176))
                        lat_float = int(lat_deg) + (float(lat_degmin) / 60)
                        if lat_float < 0:
                            lat_float = lat_float * -1
                    elif "S" in lat:
                        lat = lat.replace("S", "").replace("'", "")
                        (lat_deg, lat_degmin) = lat.split(chr(176))
                        lat_float = int(lat_deg) + (float(lat_degmin) / 60)
                        if lat_float > 0:
                            lat_float = lat_float * -1
                    else:
                        raise RuntimeError
                    if "E" in lon:
                        lon = lon.replace("E", "").replace("'", "")
                        (lon_deg, lon_degmin) = lon.split(chr(176))
                        lon_float = int(lon_deg) + (float(lon_degmin) / 60)
                        if lon_float < 0:
                            lon_float = lon_float * -1
                    elif "W" in lon:
                        lon = lon.replace("W", "").replace("'", "")
                        (lon_deg, lon_degmin) = lon.split(chr(176))
                        lon_float = int(lon_deg) + (float(lon_degmin) / 60)
                        if lon_float > 0:
                            lon_float = lon_float * -1
                    else:
                        raise RuntimeError
                else:
                    # Then the lat and lon are in different formats
                    raise RuntimeError

            elif chr(176) not in lat and chr(176) not in lon:
                # Then we are working with 1
                if "N" in lat:
                    lat_float = float(lat.replace("N", ""))
                    if lat_float < 0:
                        lat_float = lat_float * -1
                elif "S" in lat:
                    lat_float = float(lat.replace("S", ""))
                    if lat_float > 0:
                        lat_float = lat_float * -1
                else:
                    raise RuntimeError
                if "E" in lon:
                    lon_float = float(lon.replace("E", ""))
                    if lon_float < 0:
                        lon_float = lon_float * -1
                elif "W" in lon:
                    lon_float = float(lon.replace("W", ""))
                    if lon_float > 0:
                        lon_float = lon_float * -1
                else:
                    raise RuntimeError
            elif chr(176) in lat or chr(176) in lon:
                # THen there is a degree sign in only one of them and we should raise an error
                raise RuntimeError
        except Exception:
            raise RuntimeError
    return lat_float, lon_float

def dms2dec(dms_str):
    """Return decimal representation of DMS

        dms2dec(utf8(48°53'10.18"N))
        48.8866111111F

        dms2dec(utf8(2°20'35.09"E))
        2.34330555556F

        dms2dec(utf8(48°53'10.18"S))
        -48.8866111111F

        dms2dec(utf8(2°20'35.09"W))
        -2.34330555556F

        """

    dms_str = re.sub(r'\s', '', dms_str)

    sign = -1 if re.search('[swSW]', dms_str) else 1

    numbers = [*filter(len, re.split('\D+', dms_str, maxsplit=4))]

    degree = numbers[0]
    minute = numbers[1] if len(numbers) >= 2 else '0'
    second = numbers[2] if len(numbers) >= 3 else '0'
    frac_seconds = numbers[3] if len(numbers) >= 4 else '0'

    second += "." + frac_seconds
    return sign * (int(degree) + float(minute) / 60 + float(second) / 3600)
