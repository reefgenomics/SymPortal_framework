from dbApp.models import (DataSet, ReferenceSequence, DataSetSampleSequence, AnalysisType, DataSetSample,
                          DataAnalysis, DataSetSampleSequencePM, CladeCollectionType)
from multiprocessing import Queue as mp_Queue, Process, Manager, Lock as mp_Lock
from queue import Queue as mt_Queue
from threading import Thread, Lock as mt_Lock
import sys
from django import db
import os
import json
from collections import defaultdict
import pandas as pd
import numpy as np
import sp_config
import virtual_objects
from general import ThreadSafeGeneral
from exceptions import NoDataSetSampleSequencePMObjects


class OutputProfileCountTable:
    def __init__(
            self, num_proc, within_clade_cutoff, call_type, output_dir, html_dir, js_output_path_dict, date_time_str,
            force_basal_lineage_separation,
            data_set_uids_to_output=None, data_set_sample_uid_set_to_output=None,
            data_analysis_obj=None, data_analysis_uid=None, virtual_object_manager=None):
        self.force_basal_lineage_separation = force_basal_lineage_separation
        self.thread_safe_general = ThreadSafeGeneral()
        self.data_set_uid_set_to_output, self.data_set_sample_uid_set_to_output = self._init_dss_and_ds_uids(
            data_set_sample_uid_set_to_output, data_set_uids_to_output)

        self.data_analysis_obj = self._init_da_object(data_analysis_obj, data_analysis_uid)

        # Need to pass in the passed attributes rather than the self. attributes so we know which one is None
        self.virtual_object_manager = self._init_virtual_object_manager(
            virtual_object_manager, data_set_uids_to_output, data_set_sample_uid_set_to_output,
            num_proc, within_clade_cutoff)

        self.vcc_uids_to_output = self._set_vcc_uids_to_output()
        self.date_time_str = date_time_str
        self.clades_of_output = set()
        # A sorting of the vats of the output only by the len of the vccs of the output that they associate with
        # i.e. before they are then sorted by clade. This will be used when calculating the order of the samples
        self.overall_sorted_list_of_vats = None
        # The above overall_sorted_list_of_vats is then ordered by clade to produce clade_sorted_list_of_vats_to_output
        self.clade_sorted_list_of_vats_to_output = self._set_clade_sorted_list_of_vats_to_output()
        self.sorted_list_of_vdss_uids_to_output = self._set_sorted_list_of_vdss_to_output()
        self.number_of_samples = None
        self.call_type = call_type
        self.pre_headers = None
        # the number of meta rows that have been added to the dataframe.
        # This number will be used to remove the meta rows when making the abundance only dataframes
        self.number_of_meta_rows_added = 0
        # the number of header rows that have been added to the dataframe.
        # This number will be used to remove the header rows when making the abundance only dataframes
        # set of all of the species found in the vats
        self.number_of_header_rows_added = 0
        self.rel_abund_output_df, self.abs_abund_output_df = self._init_dfs()
        self.additional_info_file_as_list = []
        self.species_set = set()
        self.species_ref_dict = self._set_species_ref_dict()
        self.output_path_list = []
        self.output_dir = output_dir
        self.profiles_output_dir = os.path.join(self.output_dir, 'its2_type_profiles')
        self.html_dir = html_dir
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.profiles_output_dir, exist_ok=True)
        # dict that will hold the output file types and paths for the DataExplorer
        self.js_output_path_dict = js_output_path_dict
        self._init_output_paths()
        self.output_path_list.extend([
            self.path_to_relative_count_table_profiles_abund_and_meta,
            self.path_to_absolute_count_table_profiles_abund_and_meta])


    def _init_output_paths(self):
        self.path_to_absolute_count_table_profiles_abund_and_meta = os.path.join(
            self.profiles_output_dir, f'{self.data_analysis_obj.id}_'
                             f'{self.data_analysis_obj.name}_'
                             f'{self.date_time_str}.profiles.absolute.abund_and_meta.txt')
        self.js_output_path_dict[
            "profile_absolute_abund_meta_count"] = self.path_to_absolute_count_table_profiles_abund_and_meta

        self.path_to_absolute_count_table_profiles_abund_only = os.path.join(
            self.profiles_output_dir, f'{self.data_analysis_obj.id}_'
                             f'{self.data_analysis_obj.name}_'
                             f'{self.date_time_str}.profiles.absolute.abund_only.txt')
        self.js_output_path_dict[
            "profile_absolute_abund_only_count"] = self.path_to_absolute_count_table_profiles_abund_only

        self.path_to_absolute_count_table_profiles_meta_only = os.path.join(
            self.profiles_output_dir, f'{self.data_analysis_obj.id}_'
                             f'{self.data_analysis_obj.name}_'
                             f'{self.date_time_str}.profiles.meta_only.txt')
        self.js_output_path_dict[
            "profile_meta"] = self.path_to_absolute_count_table_profiles_meta_only

        self.path_to_relative_count_table_profiles_abund_and_meta = os.path.join(
            self.profiles_output_dir, f'{self.data_analysis_obj.id}_'
                             f'{self.data_analysis_obj.name}_'
                             f'{self.date_time_str}.profiles.relative.abund_and_meta.txt')
        self.js_output_path_dict[
            "profile_relative_abund_meta_count"] = self.path_to_relative_count_table_profiles_abund_and_meta

        self.path_to_relative_count_table_profiles_abund_only = os.path.join(
            self.profiles_output_dir, f'{self.data_analysis_obj.id}_'
                             f'{self.data_analysis_obj.name}_'
                             f'{self.date_time_str}.profiles.relative.abund_only.txt')
        self.js_output_path_dict[
            "profile_relative_abund_only_count"] = self.path_to_relative_count_table_profiles_abund_only

        self.path_to_additional_info_file = os.path.join(self.profiles_output_dir, 'additional_info.txt')
        self.js_output_path_dict[
            "profile_additional_info_file"] = self.path_to_additional_info_file

    def _set_vcc_uids_to_output(self):
        list_of_sets_of_vcc_uids_in_vdss = [
            self.virtual_object_manager.vdss_manager.vdss_dict[vdss_uid].set_of_cc_uids for vdss_uid in
            self.data_set_sample_uid_set_to_output]
        vcc_uids_to_output = list_of_sets_of_vcc_uids_in_vdss[0].union(*list_of_sets_of_vcc_uids_in_vdss[1:])
        return vcc_uids_to_output

    def output_types(self):
        print('\n\nOutputting ITS2 type profile abundance count tables\n')
        self._populate_main_body_of_dfs()

        self._populate_sample_name_series()

        self._populate_meta_info_of_dfs()

        self._write_out_dfs()

    def _populate_sample_name_series(self):
        dss_name_ordered_list = [self.virtual_object_manager.vdss_manager.vdss_dict[vdss_uid].name for vdss_uid in self.sorted_list_of_vdss_uids_to_output]

        sample_name_series_data = []

        for _ in range(len(self.pre_headers)):
            sample_name_series_data.append(np.nan)

        for name in dss_name_ordered_list:
            sample_name_series_data.append(name)

        # add two nan for the remainder
        for _ in range(len(self.abs_abund_output_df.index.tolist()) - (len(self.pre_headers) + len(dss_name_ordered_list))):
            sample_name_series_data.append(np.nan)

        sample_name_series = pd.Series(
            name='sample_name',
            data=sample_name_series_data,
            index=self.abs_abund_output_df.index.tolist())

        self.abs_abund_output_df.insert(loc=0, column='sample_name', value=sample_name_series)
        self.rel_abund_output_df.insert(loc=0, column='sample_name', value=sample_name_series)

    def _populate_meta_info_of_dfs(self):
        self._add_species_info_to_addition_info_file()

        if self.call_type == 'analysis':
            self._append_meta_info_for_analysis_call_type()
        else:
            # call_type=='stand_alone'
            self._append_meta_info_for_stand_alone_call_type()

    def _write_out_dfs(self):
        self._write_out_abund_and_meta_dfs_profiles()

        abund_row_indices = self._write_out_abund_only_dfs_profiles()

        prof_meta_only = self._write_out_meta_only_dfs_profiles(abund_row_indices)

        self._write_out_js_profiles_data_file(prof_meta_only)

        self._write_out_additional_info_profiles()

    def _write_out_additional_info_profiles(self):
        self.thread_safe_general.write_list_to_destination(destination=self.path_to_additional_info_file,
                                          list_to_write=self.additional_info_file_as_list)
        print(self.path_to_additional_info_file)

    def _write_out_js_profiles_data_file(self, prof_meta_only):
        """Here we produce the javascript objects that return the data required for viewing in the DataExplorer.
        We will produce the profile meta info object that will hold the meta information for each of the
        resolved ITS2 type profiles. The first level of keys will be profile UIDs.
        We will also output the data for drawing the rectangles, including the maximum sequence abundance of
        an individual sample. And we will output a colour dictionary of profile UID to colour.
        Finally we will output analysis meta information. We should also put out the corresponding data submission
        meta information. These will both be displayed in the
        ."""

        prof_colour_dict, sorted_profile_uids_by_local_abund = self._output_profile_meta_information(prof_meta_only)

        self._make_profile_rect_array(prof_colour_dict, prof_meta_only, sorted_profile_uids_by_local_abund)

        self._output_data_analysis_meta_info(prof_meta_only, sorted_profile_uids_by_local_abund)

    def _output_data_analysis_meta_info(self, prof_meta_only, sorted_profile_uids_by_local_abund):
        # Here output the DataAnalysis information for the DataExplorer
        # Details that we would like to output
        # UID
        # NAME
        # TIMESTAMP
        # TOTAL SAMPLES IN ANALYSIS
        # Number of unique profiles local
        # number of profile instances local
        # number of unique profiles total analysis
        # number of profile instances total analysis
        num_samples_in_analysis = len(DataSetSample.objects.filter(data_submission_from__in=[int(_) for _ in self.data_analysis_obj.list_of_data_set_uids.split(',')]))
        total_local_abund = sum(prof_meta_only['ITS2 profile abundance local'].astype(int).values)
        unique_types_in_analysis = list(AnalysisType.objects.filter(data_analysis_from=self.data_analysis_obj))
        num_unique_profiles_in_analysis = len(unique_types_in_analysis)
        num_profile_instances = len(self._chunk_query_cct_objs_from_at_objs(unique_types_in_analysis))
        data_analysis_meta_info_dict = {
            "uid": str(self.data_analysis_obj.id), "name": self.data_analysis_obj.name,
            "time_stamp": self.data_analysis_obj.time_stamp,
            "samples_in_output": len(self.data_set_sample_uid_set_to_output),
            "samples_in_analysis": str(num_samples_in_analysis),
            "unique_profile_local": str(len(sorted_profile_uids_by_local_abund)),
            "instance_profile_local": str(total_local_abund),
            "unique_profile_analysis": str(num_unique_profiles_in_analysis),
            "instances_profile_analysis": str(num_profile_instances)}
        da_meta_js_path = os.path.join(self.html_dir, 'study_data.js')
        self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
            [{'function_name': 'getDataAnalysisMetaInfo', 'python_obj': data_analysis_meta_info_dict}],
            js_outpath=da_meta_js_path)

    def _chunk_query_cct_objs_from_at_objs(self, da_obj_list):
        cct_obj_list = []
        for at_obj in self.thread_safe_general.chunks(da_obj_list):
            cct_obj_list.extend(list(CladeCollectionType.objects.filter(analysis_type_of__in=at_obj)))
        return cct_obj_list

    def _output_profile_meta_information(self, prof_meta_only):
        # js output path for profile meta
        profile_meta_js_path = os.path.join(self.html_dir, 'study_data.js')
        sorted_profile_uids_by_local_abund = prof_meta_only.sort_values(
            'ITS2 profile abundance local', ascending=False).index.values.tolist()
        # make color dict
        prof_colour_dict = self._set_colour_dict(sorted_profile_uids_by_local_abund)
        with open(os.path.join(self.html_dir, 'prof_color_dict.json'), 'w') as f:
            json.dump(fp=f, obj=prof_colour_dict)
        # first the meta information
        genera_annotation_dict = {
            'A': 'Symbiodinium', 'B': 'Breviolum', 'C': 'Cladocopium', 'D': 'Durusdinium',
            'E': 'Effrenium', 'F': 'Clade F', 'G': 'Clade G', 'H': 'Clade H', 'I': 'Clade I'}
        profile_meta_dict = {uid: {} for uid in prof_meta_only.index.values.tolist()}
        for k in profile_meta_dict.keys():
            profile_meta_dict[k]['uid'] = k
            profile_meta_dict[k]['name'] = prof_meta_only.at[k, 'ITS2 type profile']
            profile_meta_dict[k]['genera'] = genera_annotation_dict[prof_meta_only.at[k, 'Clade']]
            profile_meta_dict[k]['maj_its2_seq'] = prof_meta_only.at[k, 'Majority ITS2 sequence']
            profile_meta_dict[k]['assoc_species'] = prof_meta_only.at[k, 'Associated species']
            profile_meta_dict[k]['local_abund'] = str(prof_meta_only.at[k, 'ITS2 profile abundance local'])
            profile_meta_dict[k]['db_abund'] = str(prof_meta_only.at[k, 'ITS2 profile abundance DB'])
            profile_meta_dict[k]['seq_uids'] = prof_meta_only.at[k, 'Sequence accession / SymPortal UID']
            profile_meta_dict[k]['seq_abund_string'] = prof_meta_only.at[
                k, 'Average defining sequence proportions and [stdev]']
            profile_meta_dict[k]['color'] = prof_colour_dict[k]
        self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
            [{'function_name': 'getProfileMetaInfo', 'python_obj': profile_meta_dict}],
            js_outpath=profile_meta_js_path)
        return prof_colour_dict, sorted_profile_uids_by_local_abund

    def _set_colour_dict(self, sorted_profile_uids_by_local_abund,):
        colour_palette_pas = ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                              self.thread_safe_general.create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=1000, num_cols=50,
                                                 time_out_iterations=10000)]

        grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']

        # We will use the col headers of the df as the its2 type profile order for plotting but we
        # we should colour according to the abundance of the its2 type profiles
        # as we don't want to run out of colours by the time we get to profiles that are very abundant.
        # The sorted_type_prof_names_by_local_abund object has the names of the its2 type profile in order of abundance
        # we will use the index order as the order of samples to plot

        # create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
        # to the most abundant seqs first and after that cycle through the grey_pallette assigning colours

        colour_dict = {}
        for i in range(len(sorted_profile_uids_by_local_abund)):
            if i < len(colour_palette_pas):
                colour_dict[sorted_profile_uids_by_local_abund[i]] = colour_palette_pas[i]
            else:
                grey_index = i % len(grey_palette)
                colour_dict[sorted_profile_uids_by_local_abund[i]] = grey_palette[grey_index]
        return  colour_dict

    def _make_profile_rect_array(self, prof_colour_dict, prof_meta_only, sorted_profile_uids_by_local_abund):
        # get rid of the meta information at the top of the df
        rows_to_keep_by_index = []
        abundance_row_indices = list(range(self.number_of_header_rows_added, (
                len(self.abs_abund_output_df.index.values.tolist()) - self.number_of_meta_rows_added)))
        rows_to_keep_by_index += abundance_row_indices
        absolute_df_abund_only = self.abs_abund_output_df.iloc[rows_to_keep_by_index, :]
        relative_df_abund_only = self.rel_abund_output_df.iloc[rows_to_keep_by_index, :]
        # now we need to create the rect array for the profiles
        profile_rect_dict = {sample_uid: [] for sample_uid in absolute_df_abund_only.index.values.tolist()}
        max_cumulative_abs = 0
        for sample_uid in absolute_df_abund_only.index:
            new_rect_list = []
            abs_series = absolute_df_abund_only.loc[sample_uid]
            rel_series = relative_df_abund_only.loc[sample_uid]
            cumulative_count_abs = 0
            cumulative_count_rel = 0
            for profile_uid in sorted_profile_uids_by_local_abund:
                prof_abund_abs = abs_series.at[profile_uid]
                prof_abund_rel = rel_series.at[profile_uid]
                if prof_abund_abs:
                    cumulative_count_abs += int(prof_abund_abs)
                    cumulative_count_rel += float(prof_abund_rel)
                    new_rect_list.append({

                        "profile_name": prof_meta_only.at[profile_uid, 'ITS2 type profile'],
                        "y_abs": cumulative_count_abs,
                        "y_rel": f'{cumulative_count_rel:.3f}',
                        "height_rel": f'{float(prof_abund_rel):.3f}',
                        "height_abs": int(prof_abund_abs),

                    })
            profile_rect_dict[sample_uid] = new_rect_list
            if cumulative_count_abs > max_cumulative_abs:
                max_cumulative_abs = cumulative_count_abs
        # now we have the dictionary that holds the rectangle arrays populated
        # and we have the maximum abundance
        # now write these out as js file and functions to return.
        js_file_path = os.path.join(self.html_dir, 'study_data.js')
        self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
            [{'function_name': 'getRectDataProfileBySample', 'python_obj': profile_rect_dict},
             {'function_name': 'getRectDataProfileBySampleMaxSeq', 'python_obj': max_cumulative_abs},
             {'function_name': 'getProfColor', 'python_obj': prof_colour_dict}],
            js_outpath=js_file_path)

    def _write_out_meta_only_dfs_profiles(self, abundance_row_indices):
        # now output the meta_only dfs
        # delete the abundance info rows and drop the sample_name column
        # then output the transposed matrix so that this is essentially the meta info for each of the
        # ITS2 type profiles
        rows_to_keep_by_index = list(range(self.number_of_header_rows_added))
        rows_to_keep_by_index += list(range(self.number_of_header_rows_added + len(abundance_row_indices),
                                            len(self.abs_abund_output_df.index.values.tolist())))
        profile_meta_only = self.abs_abund_output_df.iloc[rows_to_keep_by_index, :]
        profile_meta_only.drop(columns='sample_name', inplace=True)
        profile_meta_only = profile_meta_only.T
        profile_meta_only.to_csv(self.path_to_absolute_count_table_profiles_meta_only, sep="\t", header=True,
                                 index=False)
        print(self.path_to_absolute_count_table_profiles_meta_only)
        # Make cols numerical
        profile_meta_only['ITS2 profile abundance local'] = profile_meta_only['ITS2 profile abundance local'].astype(
            int)
        profile_meta_only['ITS2 profile abundance DB'] = profile_meta_only['ITS2 profile abundance DB'].astype(
            int)
        return profile_meta_only

    def _write_out_abund_only_dfs_profiles(self):
        # write out the abund_only dfs
        # get rid of the header rows other than the profile uid that is the first row
        # also get rid of all meta rows at the end
        # also get rid of the sample_name column
        rows_to_keep_by_index = [0]
        abundance_row_indices = list(range(self.number_of_header_rows_added, (
                    len(self.abs_abund_output_df.index.values.tolist()) - self.number_of_meta_rows_added)))
        rows_to_keep_by_index += abundance_row_indices
        absolute_df_abund_only = self.abs_abund_output_df.iloc[rows_to_keep_by_index, :]
        relative_df_abund_only = self.rel_abund_output_df.iloc[rows_to_keep_by_index, :]
        absolute_df_abund_only.drop(columns='sample_name', inplace=True)
        relative_df_abund_only.drop(columns='sample_name', inplace=True)
        # replace the top item of the index with sample_uid rather than 'ITS2 type profile UID'
        new_index = absolute_df_abund_only.index.values.tolist()
        new_index[0] = 'sample_uid'
        absolute_df_abund_only.index = new_index
        relative_df_abund_only.index = new_index
        absolute_df_abund_only.to_csv(self.path_to_absolute_count_table_profiles_abund_only, sep="\t", header=False)
        relative_df_abund_only.to_csv(self.path_to_relative_count_table_profiles_abund_only, sep="\t", header=False)
        print(self.path_to_absolute_count_table_profiles_abund_only)
        print(self.path_to_relative_count_table_profiles_abund_only)
        return abundance_row_indices

    def _write_out_abund_and_meta_dfs_profiles(self):
        # write out the abund_and_meta dfs
        print('\n\nITS2 type profile count tables output to:')
        self.abs_abund_output_df.to_csv(
            self.path_to_absolute_count_table_profiles_abund_and_meta, sep="\t", header=False)
        print(f'{self.path_to_absolute_count_table_profiles_abund_and_meta}\n\n')
        self.rel_abund_output_df.to_csv(
            self.path_to_relative_count_table_profiles_abund_and_meta, sep="\t", header=False)
        print(self.path_to_relative_count_table_profiles_abund_and_meta)

    def _append_meta_info_for_stand_alone_call_type(self):
        data_sets_of_analysis = len(self.data_analysis_obj.list_of_data_set_uids.split(','))
        if self.call_type == 'stand_alone_data_sets':
            meta_info_string = self._make_stand_alone_data_set_meta_info_string(data_sets_of_analysis)
        else:
            meta_info_string = self._make_stand_alone_data_set_samples_meta_info_string(data_sets_of_analysis)
        self._append_meta_info_summary_string_to_additional_info_file(meta_info_string)
        self._append_data_set_info_to_additional_info()

    def _append_meta_info_for_analysis_call_type(self):
        meta_info_string_items = self._make_analysis_meta_info_string()
        self._append_meta_info_summary_string_to_additional_info_file(meta_info_string_items)
        self._append_data_set_info_to_additional_info()

    def _add_species_info_to_addition_info_file(self):
        # add a blank row with just the header Species reference
        self.additional_info_file_as_list.append('Species references')
        # now add the references for each of the associated species
        # we put an NaN in the first column so that we can delete the sample_name column when making the
        # meta only output table. If we put the species reference directly in that first column it would be deleted
        # when deleting this column.
        for species in self.species_set:
            if species in self.species_ref_dict.keys():
                self.additional_info_file_as_list.append(self.species_ref_dict[species])

    def _append_data_set_info_to_additional_info(self):
        ds_obj_list = self._chunk_query_ds_objs_from_ds_uids()
        for data_set_object in ds_obj_list:
            data_set_meta_str = f'Data_set ID: {data_set_object.id}; ' \
                                f'Data_set name: {data_set_object.name}; ' \
                                f'submitting_user: {data_set_object.submitting_user}; ' \
                                f'time_stamp: {data_set_object.time_stamp}'

            self.additional_info_file_as_list.append(data_set_meta_str)

    def _chunk_query_ds_objs_from_ds_uids(self):
        ds_obj_list = []
        for uid_list in self.thread_safe_general.chunks(self.data_set_uid_set_to_output):
            ds_obj_list.extend(list(DataSet.objects.filter(id__in=uid_list)))
        return ds_obj_list

    def _make_analysis_meta_info_string(self):
        meta_info_string_items = [
            f'Output as part of data_analysis ID: {self.data_analysis_obj.id}; '
            f'Number of data_set objects as part of analysis = {len(self.data_set_uid_set_to_output)}; '
            f'submitting_user: {self.data_analysis_obj.submitting_user}; '
            f'time_stamp: {self.data_analysis_obj.time_stamp}']
        return meta_info_string_items

    def _append_meta_info_summary_string_to_additional_info_file(self, meta_info_string_items):
        self.additional_info_file_as_list.append(meta_info_string_items)

    def _make_stand_alone_data_set_meta_info_string(self, data_sets_of_analysis):
        meta_info_string = f'Stand_alone_data_sets output by {sp_config.user_name} on {self.date_time_str}; ' \
                           f'data_analysis ID: {self.data_analysis_obj.id}; ' \
                           f'Number of data_set objects as part of output = {len(self.data_set_uid_set_to_output)}; ' \
                           f'Number of data_set objects as part of analysis = {data_sets_of_analysis}'
        return meta_info_string

    def _make_stand_alone_data_set_samples_meta_info_string(self, data_sets_of_analysis):
        # self.call_type == 'stand_alone_data_set_samples'
        meta_info_string = f'Stand_alone_data_set_samples output by {sp_config.user_name} on {self.date_time_str}; ' \
                           f'data_analysis ID: {self.data_analysis_obj.id}; ' \
                           f'Number of data_set objects as part of output = {len(self.data_set_uid_set_to_output)}; ' \
                           f'Number of data_set objects as part of analysis = {data_sets_of_analysis}'
        return meta_info_string

    def _populate_main_body_of_dfs(self):
        print('\nPopulating output dfs:')
        for vat in self.clade_sorted_list_of_vats_to_output:
            sys.stdout.write(f'\r{vat.name}')
            tosp = self.TypeOutputSeriesPopulation(parent_output_type_count_table=self, vat=vat)
            data_relative_list, data_absolute_list = tosp.make_output_series()
            self.rel_abund_output_df[vat.id] = data_relative_list
            self.abs_abund_output_df[vat.id] = data_absolute_list
            if vat.species != 'None':
                self.species_set.update(vat.species.split(','))

    class TypeOutputSeriesPopulation:
        """will create a relative abundance and absolute abundance
        output pandas series for a given VirtualAnalysisType
        """
        def __init__(self, parent_output_type_count_table, vat):
            self.output_type_count_table = parent_output_type_count_table
            self.vat = vat
            self.data_relative_list = []
            self.data_absolute_list = []

        def make_output_series(self):
            """type_uid"""
            self._pop_type_uid()

            self._pop_type_clade()

            self._pop_maj_seq_str()

            self._pop_species()

            self._pop_type_local_and_global_abundances()

            self._pop_type_name()

            self._pop_type_abundances()

            self._pop_vat_accession_name()

            self._pop_av_and_stdev_abund()

            return self.data_relative_list, self.data_absolute_list


        def _pop_av_and_stdev_abund(self):
            average_abund_and_sd_string = ''
            for rs_id in list(self.vat.multi_modal_detection_rel_abund_df):
                if average_abund_and_sd_string == '':
                    average_abund_and_sd_string = self._append_rel_abund_and_sd_str_for_rs(
                        average_abund_and_sd_string, rs_id)
                else:
                    average_abund_and_sd_string = self._append_dash_or_slash_if_maj_seq(average_abund_and_sd_string,
                                                                                        rs_id)
                    average_abund_and_sd_string = self._append_rel_abund_and_sd_str_for_rs(
                        average_abund_and_sd_string, rs_id)
            self.data_relative_list.append(average_abund_and_sd_string)
            self.data_absolute_list.append(average_abund_and_sd_string)

        def _append_dash_or_slash_if_maj_seq(self, average_abund_and_sd_string, rs_id):
            if rs_id in self.vat.majority_reference_sequence_uid_set:
                average_abund_and_sd_string += '/'
            else:
                average_abund_and_sd_string += '-'
            return average_abund_and_sd_string

        def _append_rel_abund_and_sd_str_for_rs(self, average_abund_and_sd_string, rs_id):
            average_abund_str = "{0:.3f}".format(self.vat.multi_modal_detection_rel_abund_df[rs_id].mean())
            std_dev_str = "{0:.3f}".format(self.vat.multi_modal_detection_rel_abund_df[rs_id].std())
            average_abund_and_sd_string += f'{average_abund_str}[{std_dev_str}]'
            return average_abund_and_sd_string

        def _pop_vat_accession_name(self):
            vat_accession_name = self.vat.generate_name(
                at_df=self.vat.multi_modal_detection_rel_abund_df,
                use_rs_ids_rather_than_names=True)
            self.data_relative_list.append(vat_accession_name)
            self.data_absolute_list.append(vat_accession_name)

        def _pop_type_abundances(self):
            # type abundances
            temp_rel_abund_holder_list = []
            temp_abs_abund_holder_list = []
            for vdss_uid in self.output_type_count_table.sorted_list_of_vdss_uids_to_output:
                count = 0
                vdss_obj = self.output_type_count_table.virtual_object_manager.vdss_manager.vdss_dict[vdss_uid]
                for vcc_uid in vdss_obj.set_of_cc_uids:
                    if vcc_uid in self.vat.type_output_rel_abund_series:
                        count += 1
                        temp_rel_abund_holder_list.append(self.vat.type_output_rel_abund_series[vcc_uid])
                        temp_abs_abund_holder_list.append(self.vat.type_output_abs_abund_series[vcc_uid])

                if count == 0:  # type not found in vdss
                    temp_rel_abund_holder_list.append(0)
                    temp_abs_abund_holder_list.append(0)
                if count > 1:  # more than one vcc from vdss associated with type
                    raise RuntimeError('More than one vcc of vdss matched vat in output')
            self.data_relative_list.extend(temp_rel_abund_holder_list)
            self.data_absolute_list.extend(temp_abs_abund_holder_list)

        def _pop_type_name(self):
            # name
            self.data_absolute_list.append(self.vat.name)
            self.data_relative_list.append(self.vat.name)

        def _pop_type_local_and_global_abundances(self):
            # local_output_type_abundance
            # all analysis_type_abundance
            vccs_of_type = self.vat.clade_collection_obj_set_profile_assignment
            vccs_of_type_from_output = [vcc for vcc in vccs_of_type if
                                        vcc.vdss_uid in self.output_type_count_table.sorted_list_of_vdss_uids_to_output]
            abund_db = self.vat.grand_tot_num_instances_of_vat_in_analysis
            self.data_absolute_list.extend([str(len(vccs_of_type_from_output)), str(abund_db)])
            self.data_relative_list.extend([str(len(vccs_of_type_from_output)), str(abund_db)])

        def _pop_species(self):
            # species
            self.data_absolute_list.append(self.vat.species)
            self.data_relative_list.append(self.vat.species)

        def _pop_maj_seq_str(self):
            # majority sequences string e.g. C3/C3b
            ordered_maj_seq_names = []
            for rs_id in list(self.vat.multi_modal_detection_rel_abund_df):
                for rs in self.vat.footprint_as_ref_seq_objs_set:
                    if rs.id == rs_id and rs in self.vat.majority_reference_sequence_obj_set:
                        ordered_maj_seq_names.append(rs.name)
            maj_seq_str = '/'.join(ordered_maj_seq_names)
            self.data_absolute_list.append(maj_seq_str)
            self.data_relative_list.append(maj_seq_str)

        def _pop_type_clade(self):
            # clade
            self.data_absolute_list.append(self.vat.clade)
            self.data_relative_list.append(self.vat.clade)

        def _pop_type_uid(self):
            # Type uid
            self.data_absolute_list.append(self.vat.id)
            self.data_relative_list.append(self.vat.id)


    def _init_da_object(self, data_analysis_obj, data_analysis_uid):
        if data_analysis_uid:
            self.data_analysis_obj = DataAnalysis.objects.get(id=data_analysis_uid)
        else:
            self.data_analysis_obj = data_analysis_obj
        return self.data_analysis_obj

    def _init_dfs(self):
        self.pre_headers = ['ITS2 type profile UID', 'Clade', 'Majority ITS2 sequence',
                       'Associated species', 'ITS2 profile abundance local', 'ITS2 profile abundance DB', 'ITS2 type profile']
        self.number_of_header_rows_added += len(self.pre_headers)
        post_headers = ['Sequence accession / SymPortal UID', 'Average defining sequence proportions and [stdev]']
        self.number_of_meta_rows_added += len(post_headers)
        self.df_index = self.pre_headers + self.sorted_list_of_vdss_uids_to_output + post_headers
        return pd.DataFrame(index=self.df_index), pd.DataFrame(index=self.df_index)

    def _init_dss_and_ds_uids(self, data_set_sample_uid_set_to_output, data_set_uids_to_output):
        if data_set_sample_uid_set_to_output:
            self.data_set_sample_uid_set_to_output = data_set_sample_uid_set_to_output

            temp_data_set_obj_list = self._chunk_query_distinct_dss_objs_from_dss_uids()
            self.data_set_uid_set_to_output = [ds.id for ds in temp_data_set_obj_list]
        else:
            self.data_set_uid_set_to_output = data_set_uids_to_output
            temp_data_set_sample_obj_list = self._chunk_query_dss_objs_from_ds_uids()

            self.data_set_sample_uid_set_to_output = [dss.id for dss in temp_data_set_sample_obj_list]
        return self.data_set_uid_set_to_output, self.data_set_sample_uid_set_to_output

    def _chunk_query_dss_objs_from_ds_uids(self):
        temp_data_set_sample_obj_list = []
        for uid_list in self.thread_safe_general.chunks(self.data_set_uid_set_to_output):
            temp_data_set_sample_obj_list.extend(list(DataSetSample.objects.filter(data_submission_from__in=uid_list)))
        return temp_data_set_sample_obj_list

    def _chunk_query_distinct_dss_objs_from_dss_uids(self):
        temp_data_set_obj_set = set()
        for uid_list in self.thread_safe_general.chunks(self.data_set_sample_uid_set_to_output):
            temp_data_set_obj_set.update(list(DataSet.objects.filter(datasetsample__in=uid_list)))
        temp_data_set_obj_list = list(temp_data_set_obj_set)
        return temp_data_set_obj_list

    def _init_virtual_object_manager(
            self, virtual_object_manager, data_set_uids_to_output, data_set_sample_uid_set_to_output,
            num_proc, within_clade_cutoff):
        if virtual_object_manager:
            return virtual_object_manager
        else:
            if data_set_uids_to_output:
                self.virtual_object_manager = virtual_objects.VirtualObjectManager(
                    num_proc=num_proc, within_clade_cutoff=within_clade_cutoff,
                    list_of_data_set_uids=data_set_uids_to_output,
                    force_basal_lineage_separation=self.force_basal_lineage_separation)
            else:
                self.virtual_object_manager = virtual_objects.VirtualObjectManager(
                    num_proc=num_proc, within_clade_cutoff=within_clade_cutoff,
                    list_of_data_set_sample_uids=data_set_sample_uid_set_to_output,
                    force_basal_lineage_separation=self.force_basal_lineage_separation)

            print('\nInstantiating VirtualAnalysisTypes')

            temp_analysis_type_obj_list = self._chunk_query_distinct_at_obj_from_dss_uids()

            for at in temp_analysis_type_obj_list:
                sys.stdout.write(f'\r{at.name}')
                self.virtual_object_manager.vat_manager.make_vat_post_profile_assignment_from_analysis_type(at)

            self._associate_vat_to_vcc()

        return self.virtual_object_manager

    def _chunk_query_distinct_at_obj_from_dss_uids(self):
        temp_analysis_type_obj_set = set()
        for uid_list in self.thread_safe_general.chunks(self.data_set_sample_uid_set_to_output):
            temp_analysis_type_obj_set.update(list(AnalysisType.objects.filter(
                cladecollectiontype__clade_collection_found_in__data_set_sample_from__in=uid_list,
                data_analysis_from=self.data_analysis_obj)))
        temp_analysis_type_obj_list = list(temp_analysis_type_obj_set)
        return temp_analysis_type_obj_list

    def _associate_vat_to_vcc(self):
        """The CladeCollections held on disc have know info on which AnalysisTypes were found in them except
        for by association with CladeCollectionTypes. As such we have to populate the
        vcc.analysis_type_obj_to_representative_rel_abund_in_cc_dict manually from the vat."""
        print('\nAssociating VirtualCladeCollections to VirtualAnalysisTypes')
        for vat in self.virtual_object_manager.vat_manager.vat_dict.values():
            for vcc in vat.clade_collection_obj_set_profile_assignment:
                vcc_rel_abund_dict = vcc.ref_seq_id_to_rel_abund_dict
                total_seq_rel_abund_for_cc = []
                for rs_uid in vat.ref_seq_uids_set:
                    total_seq_rel_abund_for_cc.append(vcc_rel_abund_dict[rs_uid])
                vcc.analysis_type_obj_to_representative_rel_abund_in_cc_dict[vat] = sum(total_seq_rel_abund_for_cc)

    def _get_data_set_uids_of_data_sets(self):
        vds_uid_set = set()
        for vdss in [vdss for vdss in self.virtual_object_manager.vdss_manager.vdss_dict.values() if
                     vdss.uid in self.data_set_sample_uid_set_to_output]:
            vds_uid_set.add(vdss.data_set_id)
        return vds_uid_set

    def _set_sorted_list_of_vdss_to_output(self):
        """Generate the list of dss uids that will be the order that we will use for the index
        of the output dataframes. The samples will be ordered by the most abundant type profiles first, and for
        each type profile the samples that had this profile as their most abundant type profile should be listed,
        within this, they should additionally be sorted by the relative abundance at which the profiles were found
        within the samples."""
        sorted_vdss_uid_list = []
        # a sanity checking set just to make sure we are not finding more vccs with more than one most abundant type
        set_of_vcc_uids_already_added = set()
        for vat in self.overall_sorted_list_of_vats:
            temp_set_vcc_rel_abund_tups = set()
            # For each vcc that had this profile found in it
            for vcc in vat.clade_collection_obj_set_profile_assignment:
                vdss_uid_of_vcc = vcc.vdss_uid
                vdss_obj = self.virtual_object_manager.vdss_manager.vdss_dict[vdss_uid_of_vcc]
                # if the clade of the vat is the most abundant clade in the vdss
                sorted_clade_abundances_tups = [tup for tup in sorted(vdss_obj.cladal_abundances_dict.items(), key=lambda x:x[1], reverse=True)]
                most_abund_clade, rel_abund_of_clade = sorted_clade_abundances_tups[0]
                if not vat.clade == most_abund_clade:
                    continue
                # this vcc is not part of the output
                if not vdss_uid_of_vcc in self.data_set_sample_uid_set_to_output:
                    continue
                # see if this vat was the most abundant vat of the vcc
                # get sorted list of the vats
                most_abundant_vat, within_clade_rel_abund = sorted(vcc.analysis_type_obj_to_representative_rel_abund_in_cc_dict.items(), key=lambda x: x[1], reverse=True)[0]
                rel_abund_in_sample = within_clade_rel_abund * rel_abund_of_clade
                # if this type was most abundant type of the vcc, add it to the profile
                if most_abundant_vat == vat:
                    if vcc.id not in set_of_vcc_uids_already_added:
                        set_of_vcc_uids_already_added.add(vcc.id)
                    else:
                        raise RuntimeError('The vcc was already added to the list_of_vcc_uids_already_added')

                    if vdss_uid_of_vcc in self.data_set_sample_uid_set_to_output:
                        if vdss_uid_of_vcc not in sorted_vdss_uid_list:
                            temp_set_vcc_rel_abund_tups.add((vdss_uid_of_vcc, rel_abund_in_sample))
                        else:
                            pass
                    else:
                        raise RuntimeError('vdss associated with vcc seems to not be part of the output')
            temp_set_vcc_rel_abund_tups = list(temp_set_vcc_rel_abund_tups)
            temp_set_vcc_rel_abund_tups.sort(key=lambda x: x[0], reverse=True)
            temp_set_vcc_rel_abund_tups.sort(key=lambda x: x[1], reverse=True)
            sorted_vdss_uid_list.extend([tup[0] for tup in temp_set_vcc_rel_abund_tups])
        # add the samples that didn't have a type associated to them in a specific order
        samples_to_add = [dss_uid for dss_uid in self.data_set_sample_uid_set_to_output if dss_uid not in sorted_vdss_uid_list]
        samples_to_add.sort(reverse=True)
        sorted_vdss_uid_list.extend(samples_to_add)
        return sorted_vdss_uid_list

    def _set_clade_sorted_list_of_vats_to_output(self):
        """Get list of analysis type sorted by clade, and then by
        len of the cladecollections associated to them from the output
        """
        list_of_tup_vat_to_vccs_of_output = []
        for vat in self.virtual_object_manager.vat_manager.vat_dict.values():
            self.clades_of_output.add(vat.clade)
            vccs_of_output_of_vat = []
            for vcc in vat.clade_collection_obj_set_profile_assignment:
                if vcc.vdss_uid in self.data_set_sample_uid_set_to_output:
                    vccs_of_output_of_vat.append(vcc)
            list_of_tup_vat_to_vccs_of_output.append((vat, len(vccs_of_output_of_vat)))

        self.overall_sorted_list_of_vats = [vat for vat, num_vcc_of_output in
                              sorted(list_of_tup_vat_to_vccs_of_output, key=lambda x: x[1], reverse=True) if num_vcc_of_output != 0]

        clade_ordered_type_order = []
        for clade in list('ABCDEFGHI'):
            clade_ordered_type_order.extend([vat for vat in self.overall_sorted_list_of_vats if vat.clade == clade])
        return clade_ordered_type_order

    def _set_species_ref_dict(self):
        return {
        'S. microadriaticum': 'Freudenthal, H. D. (1962). Symbiodinium gen. nov. and Symbiodinium microadriaticum '
                              'sp. nov., a Zooxanthella: Taxonomy, Life Cycle, and Morphology. The Journal of '
                              'Protozoology 9(1): 45-52',
        'S. pilosum': 'Trench, R. (2000). Validation of some currently used invalid names of dinoflagellates. '
                      'Journal of Phycology 36(5): 972-972.\tTrench, R. K. and R. J. Blank (1987). '
                      'Symbiodinium microadriaticum Freudenthal, S. goreauii sp. nov., S. kawagutii sp. nov. and '
                      'S. pilosum sp. nov.: Gymnodinioid dinoflagellate symbionts of marine invertebrates. '
                      'Journal of Phycology 23(3): 469-481.',
        'S. natans': 'Hansen, G. and N. Daugbjerg (2009). Symbiodinium natans sp. nob.: A free-living '
                     'dinoflagellate from Tenerife (northeast Atlantic Ocean). Journal of Phycology 45(1): 251-263.',
        'S. tridacnidorum': 'Lee, S. Y., H. J. Jeong, N. S. Kang, T. Y. Jang, S. H. Jang and T. C. Lajeunesse (2015). '
                            'Symbiodinium tridacnidorum sp. nov., a dinoflagellate common to Indo-Pacific giant clams,'
                            ' and a revised morphological description of Symbiodinium microadriaticum Freudenthal, '
                            'emended Trench & Blank. European Journal of Phycology 50(2): 155-172.',
        'S. linucheae': 'Trench, R. K. and L.-v. Thinh (1995). Gymnodinium linucheae sp. nov.: The dinoflagellate '
                        'symbiont of the jellyfish Linuche unguiculata. European Journal of Phycology 30(2): 149-154.',
        'S. minutum': 'Lajeunesse, T. C., J. E. Parkinson and J. D. Reimer (2012). A genetics-based description of '
                      'Symbiodinium minutum sp. nov. and S. psygmophilum sp. nov. (dinophyceae), two dinoflagellates '
                      'symbiotic with cnidaria. Journal of Phycology 48(6): 1380-1391.',
        'S. antillogorgium': 'Parkinson, J. E., M. A. Coffroth and T. C. LaJeunesse (2015). "New species of Clade B '
                             'Symbiodinium (Dinophyceae) from the greater Caribbean belong to different functional '
                             'guilds: S. aenigmaticum sp. nov., S. antillogorgium sp. nov., S. endomadracis sp. nov., '
                             'and S. pseudominutum sp. nov." Journal of phycology 51(5): 850-858.',
        'S. pseudominutum': 'Parkinson, J. E., M. A. Coffroth and T. C. LaJeunesse (2015). "New species of Clade B '
                            'Symbiodinium (Dinophyceae) from the greater Caribbean belong to different functional '
                            'guilds: S. aenigmaticum sp. nov., S. antillogorgium sp. nov., S. endomadracis sp. nov., '
                            'and S. pseudominutum sp. nov." Journal of phycology 51(5): 850-858.',
        'S. psygmophilum': 'Lajeunesse, T. C., J. E. Parkinson and J. D. Reimer (2012). '
                           'A genetics-based description of '
                           'Symbiodinium minutum sp. nov. and S. psygmophilum sp. nov. (dinophyceae), '
                           'two dinoflagellates '
                           'symbiotic with cnidaria. Journal of Phycology 48(6): 1380-1391.',
        'S. muscatinei': 'No reference available',
        'S. endomadracis': 'Parkinson, J. E., M. A. Coffroth and T. C. LaJeunesse (2015). "New species of Clade B '
                           'Symbiodinium (Dinophyceae) from the greater Caribbean belong to different functional '
                           'guilds: S. aenigmaticum sp. nov., S. antillogorgium sp. nov., S. endomadracis sp. nov., '
                           'and S. pseudominutum sp. nov." Journal of phycology 51(5): 850-858.',
        'S. aenigmaticum': 'Parkinson, J. E., M. A. Coffroth and T. C. LaJeunesse (2015). "New species of Clade B '
                           'Symbiodinium (Dinophyceae) from the greater Caribbean belong to different functional '
                           'guilds: S. aenigmaticum sp. nov., S. antillogorgium sp. nov., S. endomadracis sp. nov., '
                           'and S. pseudominutum sp. nov." Journal of phycology 51(5): 850-858.',
        'S. goreaui': 'Trench, R. (2000). Validation of some currently used invalid names of dinoflagellates. '
                      'Journal of Phycology 36(5): 972-972.\tTrench, R. K. and R. J. Blank (1987). '
                      'Symbiodinium microadriaticum Freudenthal, S. goreauii sp. nov., S. kawagutii sp. nov. and '
                      'S. pilosum sp. nov.: Gymnodinioid dinoflagellate symbionts of marine invertebrates. '
                      'Journal of Phycology 23(3): 469-481.',
        'S. thermophilum': 'Hume, B. C. C., C. D`Angelo, E. G. Smith, J. R. Stevens, J. Burt and J. Wiedenmann (2015).'
                           ' Symbiodinium thermophilum sp. nov., a thermotolerant symbiotic alga prevalent in corals '
                           'of the world`s hottest sea, the Persian/Arabian Gulf. Sci. Rep. 5.',
        'S. glynnii': 'LaJeunesse, T. C., D. T. Pettay, E. M. Sampayo, N. Phongsuwan, B. Brown, D. O. Obura, O. '
                      'Hoegh-Guldberg and W. K. Fitt (2010). Long-standing environmental conditions, geographic '
                      'isolation and host-symbiont specificity influence the relative ecological dominance and '
                      'genetic diversification of coral endosymbionts in the genus Symbiodinium. Journal of '
                      'Biogeography 37(5): 785-800.',
        'S. trenchii': 'LaJeunesse, T. C., D. C. Wham, D. T. Pettay, J. E. Parkinson, S. Keshavmurthy and C. A. Chen '
                       '(2014). Ecologically differentiated stress-tolerant endosymbionts in the dinoflagellate genus'
                       ' Symbiodinium (Dinophyceae) Clade D are different species. Phycologia 53(4): 305-319.',
        'S. eurythalpos': 'LaJeunesse, T. C., D. C. Wham, D. T. Pettay, J. E. Parkinson, '
                          'S. Keshavmurthy and C. A. Chen '
                          '(2014). Ecologically differentiated stress-tolerant '
                          'endosymbionts in the dinoflagellate genus'
                          ' Symbiodinium (Dinophyceae) Clade D are different species. Phycologia 53(4): 305-319.',
        'S. boreum': 'LaJeunesse, T. C., D. C. Wham, D. T. Pettay, J. E. Parkinson, S. Keshavmurthy and C. A. Chen '
                     '(2014). "Ecologically differentiated stress-tolerant endosymbionts in the dinoflagellate genus'
                     ' Symbiodinium (Dinophyceae) Clade D are different species." Phycologia 53(4): 305-319.',
        'S. voratum': 'Jeong, H. J., S. Y. Lee, N. S. Kang, Y. D. Yoo, A. S. Lim, M. J. Lee, H. S. Kim, W. Yih, H. '
                      'Yamashita and T. C. LaJeunesse (2014). Genetics and Morphology Characterize the Dinoflagellate'
                      ' Symbiodinium voratum, n. sp., (Dinophyceae) as the Sole Representative of Symbiodinium Clade E'
                      '. Journal of Eukaryotic Microbiology 61(1): 75-94.',
        'S. kawagutii': 'Trench, R. (2000). Validation of some currently used invalid names of dinoflagellates. '
                        'Journal of Phycology 36(5): 972-972.\tTrench, R. K. and R. J. Blank (1987). '
                        'Symbiodinium microadriaticum Freudenthal, S. goreauii sp. nov., S. kawagutii sp. nov. and '
                        'S. pilosum sp. nov.: Gymnodinioid dinoflagellate symbionts of marine invertebrates. '
                        'Journal of Phycology 23(3): 469-481.'
    }


class SequenceCountTableCreator:
    """ This is essentially broken into two parts. The first part goes through all of the DataSetSamples from
    the DataSets of the output and collects abundance information. The second part then puts this abundance
    information into a dataframe for both the absoulte and the relative abundance.
    This seq output can be run in two ways:
    1 - by providing DataSet uid lists
    2 - by providing DataSetSample uid lists
    Either way, after initial init, we will work on a sample by sample basis.
    """
    def __init__(
            self, symportal_root_dir, call_type, num_proc, html_dir, js_output_path_dict, date_time_str,
            no_pre_med_seqs, multiprocess, dss_uids_output_str=None, ds_uids_output_str=None, output_dir=None,
            sorted_sample_uid_list=None, analysis_obj=None):
        self.multiprocess = multiprocess
        self.thread_safe_general = ThreadSafeGeneral()
        self._init_core_vars(
            symportal_root_dir, analysis_obj, call_type, dss_uids_output_str, ds_uids_output_str, num_proc,
            output_dir, sorted_sample_uid_list, date_time_str, html_dir)
        self._init_seq_abundance_collection_objects()
        self._init_vars_for_putting_together_the_dfs()
        # dict that will hold the file type to file path of output files for the DataExplorer
        self.js_output_path_dict = js_output_path_dict
        self._init_output_paths()
        # dataframes without the meta rows at the bottom for making the javascript objects
        self.df_abs_no_meta_rows = None
        self.df_rel_no_meta_rows = None
        # variables to hold sample uid orders for the js output
        self.profile_based_sample_ordered_uids = None
        self.similarity_based_sample_ordered_uids = None
        # False by default. If true do not output premed seq info
        self.no_pre_med_seqs = no_pre_med_seqs
        


    def _init_core_vars(self, symportal_root_dir, analysis_obj, call_type, dss_uids_output_str,
                        ds_uids_output_str, num_proc,
                        output_dir, sorted_sample_uid_list, date_time_str, html_dir):
        self._check_either_dss_or_dsss_uids_provided(dss_uids_output_str, ds_uids_output_str)
        if dss_uids_output_str:
            dss_uids_for_query = [int(a) for a in dss_uids_output_str.split(',')]
            self.list_of_dss_objects = self._chunk_query_set_dss_objs_from_dss_uids(dss_uids_for_query)

            self.ds_objs_to_output = self._chunk_query_set_distinct_ds_objs_from_dss_objs()

        elif ds_uids_output_str:
            uids_of_data_sets_to_output = [int(a) for a in ds_uids_output_str.split(',')]

            self.ds_objs_to_output = self._chunk_query_set_ds_objs_from_ds_uids(uids_of_data_sets_to_output)

            self.list_of_dss_objects = self._chunk_query_set_dss_objs_from_ds_objs()

        self.ref_seqs_in_datasets = self._chunk_query_set_rs_objs_from_dss_objs()

        self.num_proc = num_proc

        self.date_time_str = date_time_str
        self._set_output_dirs(call_type, ds_uids_output_str, output_dir, symportal_root_dir, html_dir)
        self.sorted_sample_uid_list = sorted_sample_uid_list
        self.analysis_obj = analysis_obj
        self.call_type = call_type
        self.output_user = sp_config.user_name
        self.clade_list = list('ABCDEFGHI')
        set_of_clades_found = {ref_seq.clade for ref_seq in self.ref_seqs_in_datasets}
        self.ordered_list_of_clades_found = [clade for clade in self.clade_list if clade in set_of_clades_found]

    def _chunk_query_set_rs_objs_from_dss_objs(self):
        temp_ref_seqs_in_datasets_set = set()
        for uid_list in self.thread_safe_general.chunks(self.list_of_dss_objects):
            temp_ref_seqs_in_datasets_set.update(
                list(ReferenceSequence.objects.filter(datasetsamplesequence__data_set_sample_from__in=uid_list)))
        return list(temp_ref_seqs_in_datasets_set)

    def _chunk_query_set_dss_objs_from_ds_objs(self):
        temp_list_of_dss_objects = []
        for uid_list in self.thread_safe_general.chunks(self.ds_objs_to_output):
            temp_list_of_dss_objects.extend(list(DataSetSample.objects.filter(data_submission_from__in=uid_list)))
        return temp_list_of_dss_objects

    def _chunk_query_set_ds_objs_from_ds_uids(self, uids_of_data_sets_to_output):
        temp_ds_objs_to_output = []
        for uid_list in self.thread_safe_general.chunks(uids_of_data_sets_to_output):
            temp_ds_objs_to_output.extend(list(DataSet.objects.filter(id__in=uid_list)))
        return temp_ds_objs_to_output

    def _chunk_query_set_distinct_ds_objs_from_dss_objs(self):
        temp_ds_objs_to_output_set = set()
        for uid_list in self.thread_safe_general.chunks(self.list_of_dss_objects):
            temp_ds_objs_to_output_set.update(list(DataSet.objects.filter(datasetsample__in=uid_list)))
        return list(temp_ds_objs_to_output_set)

    def _chunk_query_set_dss_objs_from_dss_uids(self, dss_uids_for_query):
        temp_list_of_dss_objects = []
        for uid_list in self.thread_safe_general.chunks(dss_uids_for_query):
            temp_list_of_dss_objects.extend(list(DataSetSample.objects.filter(id__in=uid_list)))
        return temp_list_of_dss_objects

    @staticmethod
    def _check_either_dss_or_dsss_uids_provided(data_set_sample_ids_to_output_string, data_set_uids_to_output_as_comma_sep_string):
        if data_set_sample_ids_to_output_string is not None and data_set_uids_to_output_as_comma_sep_string is not None:
            raise RuntimeError('Provide either dss uids or ds uids for outputing sequence count tables')

    def _set_output_dirs(self, call_type, data_set_uids_to_output_as_comma_sep_string, output_dir, symportal_root_dir, html_dir):
        """Set both the standard output_dir where the count tables will be output and the directory to output
        the resources for the browser based data explorer"""
        if call_type == 'submission':
            self.output_dir = os.path.abspath(os.path.join(
                symportal_root_dir, 'outputs', 'loaded_data_sets', data_set_uids_to_output_as_comma_sep_string, self.date_time_str))
        else:  # call_type == 'analysis or call_type == 'stand_alone'
            self.output_dir = output_dir
        # the directory where all count tables and plots that are of the post-med seqs will be output.
        self.post_med_output_dir = os.path.join(self.output_dir, 'post_med_seqs')
        self.html_dir = html_dir
        os.makedirs(self.post_med_output_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)


    def _init_seq_abundance_collection_objects(self):
        """Output objects from first worker to be used by second worker"""
        self.dss_id_to_list_of_dsss_objects_dict_mp_dict = None
        self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict = None
        self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict = None
        # this is the list that we will use the self.annotated_dss_name_to_cummulative_rel_abund_mp_dict to create
        # it is a list of the ref_seqs_ordered first by clade then by abundance.
        self.clade_abundance_ordered_ref_seq_list = []

    def _init_vars_for_putting_together_the_dfs(self):
        # variables concerned with putting together the dataframes
        self.dss_id_to_pandas_series_results_list_dict = None
        self.output_df_absolute_post_med = None
        self.output_df_relative_post_med = None
        self.output_df_relative_pre_med = None
        self.additional_info_file = []
        self.output_seqs_fasta_as_list = []
        # the number of series (rows) that have been added to the dataframes.
        # this number will be used to delete the appropriate number of rows when producing
        # the abundance only tables.
        self.number_of_meta_rows_added = 0
        # the number of columns that are associated with meta data
        # (including the sample names, qc data, no_name summaries and user supplied data)
        # this number will be used to delete the appropriate number of columns when producing
        # the abundance only tables.
        self.number_of_meta_cols_added = 0

    def _init_output_paths(self):
        self.output_paths_list = []
        if self.analysis_obj:
            self._set_analysis_seq_table_output_paths()
        else:
            self._set_non_analysis_seq_table_output_paths()
        self.pre_med_absolute_df_path = None
        self.pre_med_relative_df_path = None
        self.pre_med_fasta_out_path = None

    def _set_non_analysis_seq_table_output_paths(self):
        self._set_non_analysis_abs_count_tab_output_paths()
        self._set_non_analysis_rel_count_tab_output_paths()
        self.output_fasta_path = os.path.join(self.post_med_output_dir, f'{self.date_time_str}.seqs.fasta')
        self.additional_info_file_path = os.path.join(self.post_med_output_dir, f'{self.date_time_str}.additional_info.txt')

    def _set_analysis_seq_table_output_paths(self):
        self._set_analysis_abs_count_tab_output_paths()
        self._set_analysis_rel_count_tab_output_paths()
        self.output_fasta_path = os.path.join(
            self.post_med_output_dir, f'{self.analysis_obj.id}_{self.analysis_obj.name}_{self.date_time_str}.seqs.fasta')
        self.js_output_path_dict["post_med_fasta"] = self.output_fasta_path
        self.additional_info_file_path = os.path.join(
            self.post_med_output_dir,
            f'{self.analysis_obj.id}_{self.analysis_obj.name}_{self.date_time_str}.additional_info.txt')
        self.js_output_path_dict["post_med_additional_info"] = self.additional_info_file_path

    def _set_analysis_rel_count_tab_output_paths(self):
        self.path_to_seq_output_abund_and_meta_df_relative = os.path.join(
            self.post_med_output_dir,
            f'{self.analysis_obj.id}_{self.analysis_obj.name}_{self.date_time_str}.seqs.relative.abund_and_meta.txt')
        self.js_output_path_dict["post_med_relative_abund_meta_count"] = self.path_to_seq_output_abund_and_meta_df_relative

        self.path_to_seq_output_abund_only_df_relative = os.path.join(
            self.post_med_output_dir,
            f'{self.analysis_obj.id}_{self.analysis_obj.name}_{self.date_time_str}.seqs.relative.abund_only.txt')
        self.js_output_path_dict["post_med_relative_abund_only_count"] = self.path_to_seq_output_abund_only_df_relative


        self.path_to_seq_output_meta_only_df_relative = os.path.join(
            self.post_med_output_dir,
            f'{self.analysis_obj.id}_{self.analysis_obj.name}_{self.date_time_str}.seqs.relative.meta_only.txt')
        self.js_output_path_dict["post_med_relative_meta_only_count"] = self.path_to_seq_output_meta_only_df_relative

    def _set_analysis_abs_count_tab_output_paths(self):
        self.path_to_seq_output_abund_and_meta_df_absolute = os.path.join(
            self.post_med_output_dir,
            f'{self.analysis_obj.id}_{self.analysis_obj.name}_{self.date_time_str}.seqs.absolute.abund_and_meta.txt')
        self.js_output_path_dict["post_med_absolute_abund_meta_count"] = self.path_to_seq_output_abund_and_meta_df_absolute

        self.path_to_seq_output_abund_only_df_absolute = os.path.join(
            self.post_med_output_dir,
            f'{self.analysis_obj.id}_{self.analysis_obj.name}_{self.date_time_str}.seqs.absolute.abund_only.txt')
        self.js_output_path_dict["post_med_absolute_abund_only_count"] = self.path_to_seq_output_abund_only_df_absolute

        self.path_to_seq_output_meta_only_df_absolute = os.path.join(
            self.post_med_output_dir,
            f'{self.analysis_obj.id}_{self.analysis_obj.name}_{self.date_time_str}.seqs.absolute.meta_only.txt')
        self.js_output_path_dict["post_med_absolute_meta_only_count"] = self.path_to_seq_output_meta_only_df_absolute

    def _set_non_analysis_rel_count_tab_output_paths(self):
        self.path_to_seq_output_abund_and_meta_df_relative = os.path.join(
            self.post_med_output_dir, f'{self.date_time_str}.seqs.relative.abund_and_meta.txt')
        self.js_output_path_dict[
            "post_med_relative_abund_meta_count"] = self.path_to_seq_output_abund_and_meta_df_relative

        self.path_to_seq_output_abund_only_df_relative = os.path.join(
            self.post_med_output_dir, f'{self.date_time_str}.seqs.relative.abund_only.txt')
        self.js_output_path_dict[
            "post_med_relative_abund_only_count"] = self.path_to_seq_output_abund_only_df_relative

        self.path_to_seq_output_meta_only_df_relative = os.path.join(
            self.post_med_output_dir, f'{self.date_time_str}.seqs.relative.meta_only.txt')
        self.js_output_path_dict[
            "post_med_relative_meta_only_count"] = self.path_to_seq_output_meta_only_df_relative


    def _set_non_analysis_abs_count_tab_output_paths(self):
        # Path to output table that contains both the meta info and the abundance info
        self.path_to_seq_output_abund_and_meta_df_absolute = os.path.join(
            self.post_med_output_dir, f'{self.date_time_str}.seqs.absolute.abund_and_meta.txt')
        self.js_output_path_dict[
            "post_med_absolute_abund_meta_count"] = self.path_to_seq_output_abund_and_meta_df_absolute

        # Path to output table that contains only the abundance info
        self.path_to_seq_output_abund_only_df_absolute = os.path.join(
            self.post_med_output_dir, f'{self.date_time_str}.seqs.absolute.abund_only.txt')
        self.js_output_path_dict[
            "post_med_absolute_abund_only_count"] = self.path_to_seq_output_abund_only_df_absolute

        # Path to output table that contains only the meta info
        self.path_to_seq_output_meta_only_df_absolute = os.path.join(
            self.post_med_output_dir, f'{self.date_time_str}.seqs.absolute.meta_only.txt')
        self.js_output_path_dict[
            "post_med_absolute_meta_only_count"] = self.path_to_seq_output_meta_only_df_absolute

    def _make_output_tables_pre_med(self):
        # 20221013 we are getting futher SSL connection time out errors on a big dataset
        # I'm going to see if resetting the db connection here will help.
        db.connections.close_all()
        pre_med_output = self.PreMedSeqOutput(parent=self)
        pre_med_output.make_pre_med_count_tables()
        self.output_df_relative_pre_med = pre_med_output.rel_count_df
        db.connections.close_all()
        

    def make_seq_output_tables(self):
        self._make_output_tables_post_med()
        if not self.no_pre_med_seqs:
            try:
                self._make_output_tables_pre_med()
            except NoDataSetSampleSequencePMObjects:
                print('\n There were no DataSetSampleSequencePM objects detected for one of the samples'
                      ' so we will assume that these objects were not generated during data loading. As such,'
                      ' there will be no pre med sequence outputs.')
        else:
            print('\n\nPre med sequence objects not generatred. Skipping output.\n\n')

    def _make_output_tables_post_med(self):
        print('\n\nOutputting sequence abundance count tables\n')
        self._collect_abundances_for_creating_the_output()

        self._generate_sample_output_series()

        self._create_ordered_output_dfs_from_series()

        self._add_uids_for_seqs_to_dfs()

        self._append_meta_info_to_additional_info_file()

        self._write_out_dfs_and_fasta()

    def _write_out_dfs_and_fasta(self):
        self._write_out_abund_and_meta_dfs_seqs()

        self._write_out_js_seq_data_files_post_med()

        self._write_out_abund_only_dfs_seqs()

        self._write_out_meta_only_dfs_seqs()

        self._write_out_seq_fasta_for_loading()

        self._write_out_additional_info_seqs()

        print('\n\nITS2 sequence output files:')
        for path_item in self.output_paths_list:
            print(path_item)

    def _write_out_additional_info_seqs(self):
        self.thread_safe_general.write_list_to_destination(destination=self.additional_info_file_path,
                                          list_to_write=self.additional_info_file)
        self.output_paths_list.append(self.additional_info_file_path)

    def _write_out_meta_only_dfs_seqs(self):
        # now do the meta output table
        # we need to drop the accession row and then we need to drop the seq abund columns
        df_abs_meta_only = self.output_df_absolute_post_med.drop(index='seq_accession', columns=self.clade_abundance_ordered_ref_seq_list)
        df_rel_meta_only = self.output_df_relative_post_med.drop(index='seq_accession', columns=self.clade_abundance_ordered_ref_seq_list)
        df_abs_meta_only.to_csv(self.path_to_seq_output_meta_only_df_absolute, sep="\t", index_label='sample_uid')
        self.output_paths_list.append(self.path_to_seq_output_meta_only_df_absolute)
        df_rel_meta_only.to_csv(self.path_to_seq_output_meta_only_df_relative, sep="\t", index_label='sample_uid')
        self.output_paths_list.append(self.path_to_seq_output_meta_only_df_relative)

    def _write_out_seq_fasta_for_loading(self):
        # we created the fasta above.
        self.thread_safe_general.write_list_to_destination(self.output_fasta_path, self.output_seqs_fasta_as_list)
        self.output_paths_list.append(self.output_fasta_path)

    def _write_out_abund_only_dfs_seqs(self):
        # get rid of the meta info rows and columns
        # we will delete the number of meta row added plus one for the accession row
        cols_to_drop = [col for col in list(self.output_df_absolute_post_med) if col not in self.clade_abundance_ordered_ref_seq_list]
        df_abs_abund_only = self.output_df_absolute_post_med.iloc[:-1 * (self.number_of_meta_rows_added + 1)].drop(columns=cols_to_drop)
        df_rel_abund_only = self.output_df_relative_post_med.iloc[:-1 * (self.number_of_meta_rows_added + 1)].drop(columns=cols_to_drop)
        df_abs_abund_only.to_csv(self.path_to_seq_output_abund_only_df_absolute, sep="\t", index_label='sample_uid')
        self.output_paths_list.append(self.path_to_seq_output_abund_only_df_absolute)
        df_rel_abund_only.to_csv(self.path_to_seq_output_abund_only_df_relative, sep="\t", index_label='sample_uid')
        self.output_paths_list.append(self.path_to_seq_output_abund_only_df_relative)

    def _write_out_abund_and_meta_dfs_seqs(self):
        self.output_df_absolute_post_med.to_csv(self.path_to_seq_output_abund_and_meta_df_absolute, sep="\t", index_label='sample_uid')
        self.output_paths_list.append(self.path_to_seq_output_abund_and_meta_df_absolute)
        self.output_df_relative_post_med.to_csv(self.path_to_seq_output_abund_and_meta_df_relative, sep="\t", index_label='sample_uid')
        self.output_paths_list.append(self.path_to_seq_output_abund_and_meta_df_relative)

    def _write_out_js_seq_data_files_post_med(self):
        '''Here we will want to create four js methods files.
        1 - sample meta information,
        2 - rectangle information
        3 - sorted sample UIDs according to orders of various properties.
        4 - the DataSubmission meta information.

        The sample meta information will be an object that has the sample
        uid as key and for each of these there will be another object that has the various property and value
        pairs.
        Rectangle array will be sample uid as key and then as single array for each sample. In each of these samples we will
        have objects, each one representing a rectange that will represent a sequence found in the sample.
        '''

        sample_meta_dict, index_of_first_seq, sample_clade_proportion_dict = self._populate_sample_meta_info_dict()

        # here we have the dictionary that will become the sample meta information

        # next we want to produce arrays that are the order of the sample uids according to various sortings
        sorted_sample_arrays = self._populate_sorted_sample_uid_arrays()

        js_file_path = os.path.join(self.html_dir, 'study_data.js')

        self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
            [{'function_name': 'getSampleMetaInfo', 'python_obj': sample_meta_dict},
             {'function_name': 'getSampleSortedArrays', 'python_obj': sorted_sample_arrays},
             {'function_name': 'getSampleProportions', 'python_obj': sample_clade_proportion_dict}],
            js_outpath=js_file_path)

        self._make_post_med_rect_array(index_of_first_seq)

        self._output_data_set_meta_info(js_file_path)

    def _output_data_set_meta_info(self, js_file_path):
        #  here make the DataSet meta info object
        # It is of course possible that we will be working with samples from multiple data sets
        # so when populating this information we will need to look to see how many datasets we were working with
        # Items to output
        # Dataset names
        # Dataset UIDs
        # Dataset time stamps
        # Number of samples
        # Average sequencing depth
        # Average Symbiodiniaceae seqs
        # self.ds_objs_to_output
        data_set_meta_info = {
            "num_associated_data_sets": str(len(self.ds_objs_to_output)),
            "ds_names": ';'.join(ds.name for ds in self.ds_objs_to_output),
            "ds_uids": ';'.join(str(ds.id) for ds in self.ds_objs_to_output),
            "ds_time_stamps": ';'.join(ds.time_stamp for ds in self.ds_objs_to_output),
            "num_samples": str(len(self.list_of_dss_objects)),
            "seq_depth_av": str(int(self.df_abs_no_meta_rows['raw_contigs'].mean())),
            'sym_seqs_absolute_av': str(int(self.df_abs_no_meta_rows['post_med_absolute'].mean())),
            "sym_seqs_unique_av": str(int(self.df_abs_no_meta_rows['post_med_unique'].mean()))
        }
        self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
            [{'function_name': 'getDataSetMetaData', 'python_obj': data_set_meta_info}],
            js_outpath=js_file_path)

    def _populate_sorted_sample_uid_arrays(self):
        # we should be able to make these quickly from the output_df_absolute by doing sorts
        # we will create a single dictionary that has the sorting prerequisite as the key and then an array of the
        # sample uids.
        # we will create a copy of the df so that we can mess it up without affecting the output df
        df_to_sort = self.df_abs_no_meta_rows.copy()
        sorted_sample_arrays = {}
        if self.profile_based_sample_ordered_uids:
            sorted_sample_arrays['profile_based'] = self.profile_based_sample_ordered_uids
        sorted_sample_arrays['similarity'] = self.similarity_based_sample_ordered_uids
        sorted_sample_arrays['sample_name'] = df_to_sort.sort_values('sample_name',
                                                                     ascending=True).index.values.tolist()
        sorted_sample_arrays['raw_contigs'] = df_to_sort.sort_values('raw_contigs',
                                                                     ascending=False).index.values.tolist()
        sorted_sample_arrays['post_med_absolute'] = df_to_sort.sort_values('post_med_absolute',
                                                                           ascending=False).index.values.tolist()
        sorted_sample_arrays['post_med_unique'] = df_to_sort.sort_values('post_med_unique',
                                                                         ascending=False).index.values.tolist()
        # Check to see if there is any taxa information else don't put this into the sorted lists
        # Go through each of the taxa keys and check to see if there is any data for them. Store false if not
        taxa_keys = ["host_phylum", "host_class", "host_order", "host_family", "host_genus", "host_species"]
        taxa_to_output_dict = {taxa_k: True for taxa_k in taxa_keys}
        for t_k in taxa_keys:
            values = df_to_sort[t_k].values.tolist()
            if values.count('NoData') == len(values):
                taxa_to_output_dict[t_k] = False
        if list(taxa_to_output_dict.values()).count(False) != len(taxa_keys):
            # Then at least some of the taxa should be output
            # Get a list of lists that will be what we make the string from
            taxa_l_l = [df_to_sort[t_k] for t_k in taxa_keys if taxa_to_output_dict[t_k]]
            df_to_sort['taxa_string'] = [';'.join(i) for i in zip(*taxa_l_l)]
            sorted_sample_arrays['taxa_string'] = df_to_sort.sort_values('taxa_string',
                                                                         ascending=True).index.values.tolist()
        # if all lats have the default value of 999 then do not output lat_lon_str
        lat_count = len([False for lat in df_to_sort['collection_latitude'].values.tolist() if lat > 90])
        if lat_count != len(df_to_sort['collection_latitude'].values.tolist()):
            df_to_sort['lat_lon_str'] = [';'.join(i) for i in zip(
                df_to_sort["collection_latitude"].map(str),
                df_to_sort["collection_longitude"].map(str))]
            sorted_sample_arrays['lat_lon'] = df_to_sort.sort_values('lat_lon_str',
                                                                     ascending=True).index.values.tolist()

        # if all collection_date are the default value of NoData then do not output
        if df_to_sort['collection_date'].values.tolist().count('NoData') != len(
                df_to_sort['collection_date'].values.tolist()):
            sorted_sample_arrays['collection_date'] = df_to_sort.sort_values('collection_date',
                                                                             ascending=True).index.values.tolist()

        # if all collection_depth are the default value of NoData then do not output
        if df_to_sort['collection_date'].values.tolist().count('NoData') != len(
                df_to_sort['collection_date'].values.tolist()):
            sorted_sample_arrays['collection_depth'] = df_to_sort.sort_values('collection_depth',
                                                                                  ascending=True).index.values.tolist()
        return sorted_sample_arrays

    def _make_post_med_rect_array(self, index_of_first_seq):
        # now we need to create the rectangle array
        post_med_rect_dict = {sample_uid: [] for sample_uid in self.df_abs_no_meta_rows.index.values.tolist()}
        # The sequences in the output df are ordered by clade and then abundance. We need to provide the
        # y_abs and y_rel properties of the rect objects in order of the most abundant sequences first.
        # We want this order to be independent of clade so we need to do a sorting.
        abundance_dict = {}
        for col in list(self.df_abs_no_meta_rows)[index_of_first_seq:]:
            abundance_dict[col] = sum(self.df_abs_no_meta_rows[col])
        # get the names of the sequences sorted according to their totalled abundance
        sorted_sorted_seq_names = list(abundance_dict.items())
        sorted_sorted_seq_names.sort(key=lambda x: x[0], reverse=True)
        sorted_sorted_seq_names.sort(key=lambda x: x[1], reverse=True)
        sorted_seq_names = [x[0] for x in sorted_sorted_seq_names]
        # The col_dict output happens before the plotting where the col dict is
        # created so we will have to create the col
        # dict here and then read it back in for making the pre_seq rect array.
        # We will then read these in for the plotting.
        seq_colour_dict = self.thread_safe_general.set_seq_colour_dict(sorted_seq_names)
        with open(os.path.join(self.html_dir, 'color_dict_post_med.json'), 'w') as f:
            json.dump(fp=f, obj=seq_colour_dict)

        max_cumulative_abs = self._populate_post_med_rect_dict(post_med_rect_dict, sorted_seq_names)
        # now we have the dictionary that holds the rectangle arrays populated
        # and we have the maximum abundance
        # now write these out as js file and functions to return.
        js_file_path = os.path.join(self.html_dir, 'study_data.js')
        self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
            [{'function_name': 'getRectDataPostMEDBySample', 'python_obj': post_med_rect_dict},
             {'function_name': 'getRectDataPostMEDBySampleMaxSeq', 'python_obj': max_cumulative_abs},
             {'function_name': 'getSeqColorPostMED', 'python_obj': seq_colour_dict}],
            js_outpath=js_file_path)

    def _populate_post_med_rect_dict(self, post_med_rect_dict, sorted_seq_names):
        max_cumulative_abs = 0
        for sample_uid in self.df_abs_no_meta_rows.index:
            new_rect_list = []
            abs_series = self.df_abs_no_meta_rows.loc[sample_uid]
            rel_series = self.df_rel_no_meta_rows.loc[sample_uid]
            cumulative_count_abs = 0
            cumulative_count_rel = 0

            for seq in sorted_seq_names:
                seq_abund_abs = abs_series.at[seq]
                seq_abund_rel = rel_series.at[seq]
                if seq_abund_abs:
                    cumulative_count_abs += int(seq_abund_abs)
                    cumulative_count_rel += float(seq_abund_rel)
                    new_rect_list.append({
                        "seq_name": seq,
                        "y_abs": cumulative_count_abs,
                        "y_rel": f'{cumulative_count_rel:.3f}',
                        "height_rel": f'{float(seq_abund_rel):.3f}',
                        "height_abs": int(seq_abund_abs),
                    })

            post_med_rect_dict[sample_uid] = new_rect_list
            if cumulative_count_abs > max_cumulative_abs:
                max_cumulative_abs = cumulative_count_abs
        return max_cumulative_abs

    def _populate_sample_meta_info_dict(self):
        # first lets produce the meta information.
        # dictionary of where the sample UID is the key to a second dictionary that contains the other properties
        # We also want to produce a dictionary where genera is key to another dict, where there are
        # key values that are samples, which then have a dict as value
        # that is k, v pairs that are the absolute and relative abundnce of that sample and clade
        # These will be used when annotating the sample information in the betwn sample distance
        # plots for tooltips and in the info section at the bottom.
        self.df_abs_no_meta_rows = self.output_df_absolute_post_med.iloc[:-1 * (self.number_of_meta_rows_added + 1)]
        self.df_rel_no_meta_rows = self.output_df_relative_post_med.iloc[:-1 * (self.number_of_meta_rows_added + 1)]
        sample_clade_proportion_dict = defaultdict(dict)
        sample_meta_dict = {int(uid): {} for uid in self.df_abs_no_meta_rows.index.values.tolist()}
        taxa_fields = ['host_phylum', 'host_class', 'host_order', 'host_family', 'host_genus', 'host_species']
        index_of_first_seq = list(self.df_abs_no_meta_rows).index('collection_depth') + 1
        for k in sample_meta_dict.keys():
            sample_meta_dict[k]['uid'] = k
            sample_meta_dict[k]['name'] = self.df_abs_no_meta_rows.at[k, 'sample_name']
            sample_meta_dict[k]['raw_contigs'] = self.df_abs_no_meta_rows.at[k, 'raw_contigs']
            sample_meta_dict[k]['post_taxa_id_absolute_symbiodiniaceae_seqs'] = self.df_abs_no_meta_rows.at[
                k, 'post_taxa_id_absolute_symbiodiniaceae_seqs']
            sample_meta_dict[k]['post_taxa_id_unique_symbiodiniaceae_seqs'] = self.df_abs_no_meta_rows.at[
                k, 'post_taxa_id_unique_symbiodiniaceae_seqs']
            sample_meta_dict[k]['post_taxa_id_absolute_non_symbiodiniaceae_seqs'] = self.df_abs_no_meta_rows.at[
                k, 'post_taxa_id_absolute_non_symbiodiniaceae_seqs']
            sample_meta_dict[k]['post_taxa_id_unique_non_symbiodiniaceae_seqs'] = self.df_abs_no_meta_rows.at[
                k, 'post_taxa_id_unique_non_symbiodiniaceae_seqs']
            sample_meta_dict[k]['post_med_absolute'] = self.df_abs_no_meta_rows.at[k, 'post_med_absolute']
            sample_meta_dict[k]['post_med_unique'] = self.df_abs_no_meta_rows.at[k, 'post_med_unique']
            sample_meta_dict[k]['sample_type'] = self.df_abs_no_meta_rows.at[k, 'sample_type']
            taxa_string = ';'.join([self.df_abs_no_meta_rows.at[k, taxa_field] for taxa_field in taxa_fields])
            sample_meta_dict[k]['taxa_string'] = taxa_string
            sample_meta_dict[k]['lat'] = str(self.df_abs_no_meta_rows.at[k, 'collection_latitude'])
            sample_meta_dict[k]['lon'] = str(self.df_abs_no_meta_rows.at[k, 'collection_longitude'])
            sample_meta_dict[k]['collection_date'] = self.df_abs_no_meta_rows.at[k, 'collection_date']
            sample_meta_dict[k]['collection_depth'] = self.df_abs_no_meta_rows.at[k, 'collection_depth']

            self._pop_clade_prop_string_and_prop_dict(k, sample_meta_dict, index_of_first_seq,
                                                      sample_clade_proportion_dict)

        return sample_meta_dict, index_of_first_seq, sample_clade_proportion_dict

    def _pop_clade_prop_string_and_prop_dict(self, k, sample_meta_dict, index_of_first_seq, sample_clade_proportion_dict):
        # get the clade breakdown from the abundances of the sequences
        # both absolute and relative
        genera_annotation_dict = {'A': 'Symbiodinium', 'B': 'Breviolum', 'C': 'Cladocopium', 'D': 'Durusdinium',
                                  'E': 'Effrenium', 'F': 'Clade F', 'G': 'Clade G', 'H': 'Clade H', 'I': 'Clade I'}
        prop_counting_dict = {clade: 0 for clade in list('ABCDEFGHI')}
        for seq_index in list(self.output_df_absolute_post_med)[index_of_first_seq:]:
            seq_abund = self.output_df_absolute_post_med.at[k, seq_index]
            if seq_index[0] in list('ABCDEFGHI'):
                prop_counting_dict[seq_index[0]] += seq_abund
            elif seq_index[-1] in list('ABCDEFGHI'):
                prop_counting_dict[seq_index[-1]] += seq_abund
            else:
                raise RuntimeError('Unrecognised clade identifier in seq name')
        # now turn these into proportions and a string semicolon seperated
        tot_seqs = sum(prop_counting_dict.values())
        if tot_seqs:
            props = [prop_counting_dict[clade] / tot_seqs for clade in list('ABCDEFGHI')]
        else:
            props = [0 for clade in list('ABCDEFGHI')]
        sample_meta_dict[k]['clade_prop_string'] = ';'.join([f'{prop:.2f}' for prop in props])
        sample_meta_dict[k]['clade_abs_abund_string'] = ';'.join([str(prop_counting_dict[clade]) for clade in list('ABCDEFGHI')])
        for clade_key, value in prop_counting_dict.items():
            if value:
                # then there were some sequences from this clade
                sample_clade_proportion_dict[genera_annotation_dict[clade_key]][k] = {'absolute':value, 'relative':value/tot_seqs}


    def _append_meta_info_to_additional_info_file(self):
        # Now append the meta infromation for each of the data_sets that make up the output contents
        # this is information like the submitting user, what the uids of the datasets are etc.
        # There are several ways that this can be called.
        # it can be called as part of the submission: call_type = submission
        # part of an analysis output: call_type = analysis
        # or stand alone: call_type = 'stand_alone'
        # we should have an output for each scenario
        if self.call_type == 'submission':
            self._append_meta_info_to_additional_info_file_submission()
        elif self.call_type == 'analysis':
            self._append_meta_info_to_additional_info_file_analysis()
        else:
            # =='stand_alone'call_type
            self._append_meta_info_to_additional_info_file_stand_alone()

    def _append_meta_info_to_additional_info_file_submission(self):
        data_set_object = self.ds_objs_to_output[0]
        # there will only be one data_set object
        meta_info_string = f'Output as part of data_set submission ID: {data_set_object.id}; ' \
                           f'submitting_user: {data_set_object.submitting_user}; ' \
                           f'time_stamp: {data_set_object.time_stamp}'

        self.additional_info_file.append(meta_info_string)

    def _increase_number_meta_series_added(self):
        # keep track of the number of meta series we will need to remove from the dataframe when producing the
        # abundance only output table
        self.number_of_meta_rows_added += 1

    def _append_meta_info_to_additional_info_file_analysis(self):

        num_data_set_objects_as_part_of_analysis = len(self.analysis_obj.list_of_data_set_uids.split(','))
        meta_info_string = f'Output as part of data_analysis ID: {self.analysis_obj.id}; ' \
                           f'Number of data_set objects as part of analysis = {num_data_set_objects_as_part_of_analysis}; ' \
                           f'submitting_user: {self.analysis_obj.submitting_user}; time_stamp: {self.analysis_obj.time_stamp}'
        self.additional_info_file.append(meta_info_string)

        for data_set_object in self.ds_objs_to_output:
            data_set_meta_str = f'Data_set ID: {data_set_object.id}; ' \
                                f'Data_set name: {data_set_object.name}; ' \
                                f'submitting_user: {data_set_object.submitting_user}; ' \
                                f'time_stamp: {data_set_object.time_stamp}'

            self.additional_info_file.append(data_set_meta_str)

    def _append_meta_info_to_additional_info_file_stand_alone(self):
        meta_info_string = f'Stand_alone output by {self.output_user} on {self.date_time_str}; ' \
                           f'Number of data_set objects as part of output = {len(self.ds_objs_to_output)}'
        self.additional_info_file.append(meta_info_string)

        for data_set_object in self.ds_objs_to_output:
            data_set_meta_str = f'Data_set ID: {data_set_object.id}; Data_set name: {data_set_object.name}; ' \
                                f'submitting_user: {data_set_object.submitting_user}; ' \
                                f'time_stamp: {data_set_object.time_stamp}'
            self.additional_info_file.append(data_set_meta_str)

    def _add_uids_for_seqs_to_dfs(self):
        """Now add the UID for each of the sequences"""
        sys.stdout.write('\nGenerating accession and fasta\n')
        reference_sequences_in_data_sets_no_name = self._chunk_query_rs_objs_no_name_from_dss_objs()
        reference_sequences_in_data_sets_has_name = self._chunk_query_rs_objs_with_name_from_dss_objs()

        no_name_dict = {rs.id: rs.sequence for rs in reference_sequences_in_data_sets_no_name}
        has_name_dict = {rs.name: (rs.id, rs.sequence) for rs in reference_sequences_in_data_sets_has_name}
        accession_list = []
        num_cols = len(list(self.output_df_relative_post_med))
        for i, col_name in enumerate(list(self.output_df_relative_post_med)):
            sys.stdout.write('\rAppending accession info and creating fasta {}: {}/{}'.format(col_name, i, num_cols))
            if col_name in self.clade_abundance_ordered_ref_seq_list:
                if col_name[-2] == '_':
                    col_name_id = int(col_name[:-2])
                    accession_list.append(str(col_name_id))
                    self.output_seqs_fasta_as_list.append('>{}'.format(col_name))
                    self.output_seqs_fasta_as_list.append(no_name_dict[col_name_id])
                else:
                    col_name_tup = has_name_dict[col_name]
                    accession_list.append(str(col_name_tup[0]))
                    self.output_seqs_fasta_as_list.append('>{}'.format(col_name))
                    self.output_seqs_fasta_as_list.append(col_name_tup[1])
            else:
                accession_list.append(np.nan)
        temp_series = pd.Series(accession_list, name='seq_accession', index=list(self.output_df_relative_post_med))
        self.output_df_absolute_post_med = self.output_df_absolute_post_med.append(temp_series)
        self.output_df_relative_post_med = self.output_df_relative_post_med.append(temp_series)

    def _chunk_query_rs_objs_with_name_from_dss_objs(self):
        reference_sequences_in_data_sets_no_name_set = set()
        for uid_list in self.thread_safe_general.chunks(self.list_of_dss_objects):
            reference_sequences_in_data_sets_no_name_set.update(list(
                ReferenceSequence.objects.filter(datasetsamplesequence__data_set_sample_from__in=uid_list,
                                                 has_name=True)))
        reference_sequences_in_data_sets_no_name = list(reference_sequences_in_data_sets_no_name_set)
        return reference_sequences_in_data_sets_no_name

    def _chunk_query_rs_objs_no_name_from_dss_objs(self):
        reference_sequences_in_data_sets_no_name_set = set()
        for uid_list in self.thread_safe_general.chunks(self.list_of_dss_objects):
            reference_sequences_in_data_sets_no_name_set.update(list(
                ReferenceSequence.objects.filter(datasetsamplesequence__data_set_sample_from__in=uid_list,
                                                 has_name=False)))
        reference_sequences_in_data_sets_no_name = list(reference_sequences_in_data_sets_no_name_set)
        return reference_sequences_in_data_sets_no_name

    def _create_ordered_output_dfs_from_series(self):
        """Put together the pandas series that hold sequences abundance outputs for each sample in order of the samples
        either according to a predefined ordered list or by an order that will be generated below.
        For the javascript outputs we want to capture these sample orders. Ideally it would be good
        to have both the profile-based order and the similarity order. If a sorted_sample_uid_list
        is provided then this is the sample-based. If this is not provided (because there is not an analysis
        associated with this output) then we can only get the similarity output."""
        if self.sorted_sample_uid_list:
            sys.stdout.write('\nValidating sorted sample list and ordering dataframe accordingly\n')
            self._check_sorted_sample_list_is_valid()
            self.profile_based_sample_ordered_uids = self.sorted_sample_uid_list
            # even though we have the profile order of samples, still calculate the similarity order
            self.similarity_based_sample_ordered_uids = self._generate_ordered_sample_list()
            self._create_ordered_output_dfs_from_series_with_sorted_sample_list()

        else:
            sys.stdout.write('\nGenerating ordered sample list and ordering dataframe accordingly\n')
            self.sorted_sample_uid_list = self._generate_ordered_sample_list()
            self.similarity_based_sample_ordered_uids = self.sorted_sample_uid_list
            self._create_ordered_output_dfs_from_series_with_sorted_sample_list()

    def _generate_ordered_sample_list(self):
        """ Returns a list which is simply the ids of the samples ordered
        This will order the samples according to which sequence is their most abundant.
        I.e. samples found to have the sequence which is most abundant in the largest number of sequences
        will be first. Within each maj sequence, the samples will be sorted by the abundance of that sequence
        in the sample.
        At the moment we are also ordering by clade just so that you see samples with the A's at the top
        of the output so that we minimise the number of 0's in the top left of the output
        honestly I think we could perhaps get rid of this and just use the over all abundance of the sequences
        discounting clade. This is what we do for the clade order when plotting.
        """
        output_df_relative = self._make_raw_relative_abund_df_from_series()
        ordered_sample_list = self._get_sample_order_from_rel_seq_abund_df(output_df_relative)
        return ordered_sample_list

    def _make_raw_relative_abund_df_from_series(self):
        output_df_relative = pd.concat(
            [list_of_series[1] for list_of_series in self.dss_id_to_pandas_series_results_list_dict.values()],
            axis=1)
        output_df_relative = output_df_relative.T
        # now remove the rest of the non abundance columns
        cols_to_drop = [col for col in list(output_df_relative) if col not in self.clade_abundance_ordered_ref_seq_list]
        output_df_relative.drop(columns=cols_to_drop, inplace=True)
        return output_df_relative

    def _get_sample_order_from_rel_seq_abund_df(self, sequence_only_df_relative):

        max_seq_ddict, no_maj_samps, seq_to_samp_ddict = self._generate_most_abundant_sequence_dictionaries(
            sequence_only_df_relative)

        return self._generate_ordered_sample_list_from_most_abund_seq_dicts(max_seq_ddict, no_maj_samps,
                                                                            seq_to_samp_ddict)

    @staticmethod
    def _generate_ordered_sample_list_from_most_abund_seq_dicts(max_seq_ddict, no_maj_samps, seq_to_samp_ddict):
        # then once we have compelted this for all sequences go clade by clade
        # and generate the sample order
        ordered_sample_list_by_uid = []
        sys.stdout.write('\nGoing clade by clade sorting by abundance\n')
        for clade in list('ABCDEFGHI'):
            sys.stdout.write(f'\rGetting clade {clade} seqs')
            tup_list_of_clade = []
            # get the clade specific list of the max_seq_ddict
            for k, v in max_seq_ddict.items():
                sys.stdout.write('\r{}'.format(k))
                if k.startswith(clade) or k[-2:] == '_{}'.format(clade):
                    tup_list_of_clade.append((k, v))

            if not tup_list_of_clade:
                continue
            # now get an ordered list of the sequences for this clade
            sys.stdout.write('\rOrdering clade {} seqs'.format(clade))
            tup_list_of_clade.sort(key=lambda x: x[0], reverse=True)
            tup_list_of_clade.sort(key=lambda x: x[1], reverse=True)
            ordered_sequence_of_clade_list = [x[0] for x in tup_list_of_clade]

            for seq_to_order_samples_by in ordered_sequence_of_clade_list:
                sys.stdout.write('\r{}'.format(seq_to_order_samples_by))
                tup_list_of_samples_that_had_sequence_as_most_abund = seq_to_samp_ddict[seq_to_order_samples_by]
                tup_list_of_samples_that_had_sequence_as_most_abund.sort(key=lambda x: x[0], reverse=True)
                tup_list_of_samples_that_had_sequence_as_most_abund.sort(key=lambda x: x[1], reverse=True)
                ordered_list_of_samples_for_seq_ordered = \
                    [x[0] for x in
                     tup_list_of_samples_that_had_sequence_as_most_abund]
                ordered_sample_list_by_uid.extend(ordered_list_of_samples_for_seq_ordered)
        # finally add in the samples that didn't have a maj sequence
        no_maj_samps.sort(reverse=True)
        ordered_sample_list_by_uid.extend(no_maj_samps)
        return ordered_sample_list_by_uid

    def _generate_most_abundant_sequence_dictionaries(self, sequence_only_df_relative):
        # {sequence_name_found_to_be_most_abund_in_sample: num_samples_it_was_found_to_be_most_abund_in}
        max_seq_ddict = defaultdict(int)
        # {most_abundant_seq_name: [(dss.id, rel_abund_of_most_abund_seq) for samples with that seq as most abund]}
        seq_to_samp_ddict = defaultdict(list)
        # a list to hold the names of samples in which there was no most abundant sequence identified
        no_maj_samps = []
        for sample_to_sort_uid in sequence_only_df_relative.index.values.tolist():
            sys.stdout.write(f'\r{sample_to_sort_uid}: Getting maj seq for sample')
            sample_series_as_float = self._get_sample_seq_abund_info_as_pd_series_float_type(
                sample_to_sort_uid, sequence_only_df_relative)
            max_rel_abund = self._get_rel_abund_of_most_abund_seq(sample_series_as_float)
            if not max_rel_abund > 0:
                no_maj_samps.append(sample_to_sort_uid)
            else:
                max_abund_seq = self._get_name_of_most_abundant_seq(sample_series_as_float)
                # add a tup of sample name and rel abund of seq to the seq_to_samp_dict
                seq_to_samp_ddict[max_abund_seq].append((sample_to_sort_uid, max_rel_abund))
                # add this to the ddict count
                max_seq_ddict[max_abund_seq] += 1
        return max_seq_ddict, no_maj_samps, seq_to_samp_ddict

    @staticmethod
    def _get_sample_seq_abund_info_as_pd_series_float_type(sample_to_sort_uid, sequence_only_df_relative):
        return sequence_only_df_relative.loc[sample_to_sort_uid].astype('float')

    @staticmethod
    def _get_rel_abund_of_most_abund_seq(sample_series_as_float):
        return sample_series_as_float.max()

    @staticmethod
    def _get_name_of_most_abundant_seq(sample_series_as_float):
        max_abund_seq = sample_series_as_float.idxmax()
        return max_abund_seq

    def _create_ordered_output_dfs_from_series_with_sorted_sample_list(self):
        # NB I was originally performing the concat directly on the managedSampleOutputDict (i.e. the mp dict)
        # but this was starting to produce errors. Starting to work on the dss_id_to_pandas_series_results_list_dict
        #  (i.e. normal, not mp, dict) seems to not produce these errors.
        sys.stdout.write('\rPopulating the absolute dataframe with series. This could take a while...')
        output_df_absolute = pd.concat(
            [list_of_series[0] for list_of_series in
             self.dss_id_to_pandas_series_results_list_dict.values()], axis=1)
        sys.stdout.write('\rPopulating the relative dataframe with series. This could take a while...')
        output_df_relative = pd.concat(
            [list_of_series[1] for list_of_series in
             self.dss_id_to_pandas_series_results_list_dict.values()], axis=1)
        # now transpose
        output_df_absolute = output_df_absolute.T
        output_df_relative = output_df_relative.T
        # now make sure that the order is correct.
        self.output_df_absolute_post_med = output_df_absolute.reindex(self.sorted_sample_uid_list)
        self.output_df_relative_post_med = output_df_relative.reindex(self.sorted_sample_uid_list)

    def _check_sorted_sample_list_is_valid(self):
        if len(self.sorted_sample_uid_list) != len(self.list_of_dss_objects):
            raise RuntimeError({'message': 'Number of items in sorted_sample_list do not match those to be outputted!'})
        if self._smpls_in_sorted_smpl_list_not_in_list_of_samples():
            raise RuntimeError(
                {'message': 'Sample list passed in does not match sample list from db query'})

    def _smpls_in_sorted_smpl_list_not_in_list_of_samples(self):
        return list(
            set(self.sorted_sample_uid_list).difference(set([dss.id for dss in self.list_of_dss_objects])))

    def _generate_sample_output_series(self):
        """This generate a pandas series for each of the samples. It uses the ordered ReferenceSequence list created
         in the previous method as well as the other two dictionaries made.
         One df for absolute abundances and one for relative abundances. These series will be put together
         and ordered to construct the output data frames that will be written out for the user.
         These dataframes contain the sequence abundances and all of the sample meta information.
        """
        seq_count_table_output_series_generator_handler = SeqOutputSeriesGeneratorHandler(parent=self)
        seq_count_table_output_series_generator_handler.execute_sequence_count_table_dataframe_contructor_handler()
        self.dss_id_to_pandas_series_results_list_dict = \
            dict(seq_count_table_output_series_generator_handler.dss_id_to_pandas_series_results_list_mp_dict)

    def _collect_abundances_for_creating_the_output(self):
        seq_collection_handler = SequenceCountTableCollectAbundanceHandler(parent_seq_count_tab_creator=self)
        seq_collection_handler.execute_sequence_count_table_ordered_seqs_worker()
        # update the dictionaries that will be used in the second worker from the first worker
        self.update_dicts_for_the_second_worker_from_first_worker(seq_collection_handler)

    def update_dicts_for_the_second_worker_from_first_worker(self, seq_collection_handler):
        self.dss_id_to_list_of_dsss_objects_dict_mp_dict = \
            seq_collection_handler.dss_id_to_list_of_dsss_objects_mp_dict

        self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict = \
            seq_collection_handler.\
                dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict

        self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict = \
            seq_collection_handler.\
                dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict

        self.clade_abundance_ordered_ref_seq_list = \
            seq_collection_handler.clade_abundance_ordered_ref_seq_list

    class PreMedSeqOutput:
        def __init__(self, parent):
            self.parent = parent
            # the directory that is specific to the pre-med outputs
            self.pre_med_dir = os.path.join(self.parent.output_dir, 'pre_med_seqs')
            os.makedirs(self.pre_med_dir, exist_ok=True)
            # the directory one above the pre_med_dir that contins the main data_loading outputs
            self.root_output_dir = self.parent.output_dir
            # the directory that will house the ouputs for html
            self.html_dir = self.parent.html_dir

            # dict that will have dss uid as key and an abundance dictionary as value
            self.master_sample_uid_to_abund_dict = {}
            # dict that will have the rs uid as key and the cummulative relative abundance as value
            self.master_sequence_count_dict = defaultdict(float)

            # The dataframes that will become the output tables
            self.abs_count_df = None
            self.rel_count_df = None

            # The directories for the output tables
            self.parent.pre_med_absolute_df_path = os.path.join(self.pre_med_dir, 'pre_med_absolute_abundance_df.csv')
            self.parent.pre_med_relative_df_path = os.path.join(self.pre_med_dir, 'pre_med_relative_abundance_df.csv')
            self.parent.pre_med_fasta_out_path = os.path.join(self.pre_med_dir, 'pre_med_master_seqs.fasta')
            self.parent.js_output_path_dict["pre_med_absolute_count"] = self.parent.pre_med_absolute_df_path
            self.parent.js_output_path_dict["pre_med_relative_count"] = self.parent.pre_med_relative_df_path
            self.parent.js_output_path_dict["pre_med_fasta"] = self.parent.pre_med_fasta_out_path
            # The dictionary that will be used to create the master fasta file for all pre_med_seqs
            self.master_fasta_dict = {}

            # dictionary to keep track of the id to dss name so we can output name in the dataframe as well.
            self.sample_uid_to_name_dict = {}



        def make_pre_med_count_tables(self):
            # output two tables, one relative abundance and one absolute
            # will also ouput a single .js file that allows us to retrieve the absolute and relative abund info
            # we will use the sample order that we can inherit from the outer sequence output
            # we will calculate our own sequence order purely according to the abundance of reference sequences
            # associated with each of the DataSetSampleSequencePM objects
            # we go through each sample by sample and populate a master dict that holds abundance
            # to reference sequence uid, and then for each sample we can have a dict of that is accessed through the
            # uid of the sample and then holds another dictionary that is refseq uid to abundance in that sample

            self._count_sequences_in_each_sample()

            self._output_pre_med_master_fasta()

            self._populate_dfs_in_parent_sample_order()

            self._write_out_dfs_as_csv()

            self._write_out_js_seq_data_file_pre_med()

        def _write_out_js_seq_data_file_pre_med(self):
            '''Here we will want to put out one js file.
            The js file will be the the rectangle array. This
            is going to be sample uid as key and then as single array for each sample. In each of these samples we will
            have objects, each one representing a rectange that will represent a sequence found in the sample.'''


            self._make_pre_med_rect_array()

        def _make_pre_med_rect_array(self):
            # now we need to create the rectangle array
            pre_med_rect_dict = {sample_uid: [] for sample_uid in
                                  self.abs_count_df.index.values.tolist()}


            # get the names of the sequences sorted according to their totalled abundance
            # NB these are already sorted purely by abundance in the premed df
            sorted_seq_names = list(self.abs_count_df)[1:]

            # get the colour dict
            with open(os.path.join(self.html_dir, 'color_dict_post_med.json'), 'r') as f:
                c_dict_post_med = json.load(f)
            seq_colour_dict = self.parent.thread_safe_general.set_seq_colour_dict_w_reference_c_dict(sorted_seq_names, c_dict_post_med)
            with open(os.path.join(self.html_dir, 'color_dict_pre_med.json'), 'w') as f:
                json.dump(fp=f, obj=seq_colour_dict)
            # merge the two colour dictionaries and add to the html output
            combi_color_dict = {**c_dict_post_med, **seq_colour_dict}

            max_cumulative_abs = self._populate_pre_med_rect_dict(
                pre_med_rect_dict, seq_colour_dict, sorted_seq_names)
            # now we have the dictionary that holds the rectangle arrays populated
            # and we have the maximum abundance
            # now write these out as js file and functions to return.
            js_file_path = os.path.join(self.html_dir, 'study_data.js')
            # For the time being we are going to no longer output the pre_med_rect_dict data
            # as it very slow to render online using d3 and takes up a huge amount of space.
            # However, if we want to reimpleent its output in the future we can simply uncomment the below
            # general.write_out_js_file_to_return_python_objs_as_js_objs(
            #     [{'function_name': 'getRectDataPreMEDBySample', 'python_obj': pre_med_rect_dict},
            #      {'function_name': 'getRectDataPreMEDBySampleMaxSeq', 'python_obj': max_cumulative_abs},
            #      {'function_name': 'getSeqColor', 'python_obj': combi_color_dict}],
            #     js_outpath=js_file_path)
            self.parent.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
                [{'function_name': 'getRectDataPreMEDBySampleMaxSeq', 'python_obj': max_cumulative_abs},
                 {'function_name': 'getSeqColor', 'python_obj': combi_color_dict}],
                js_outpath=js_file_path)

        def _populate_pre_med_rect_dict(self, post_med_rect_dict, seq_colour_dict,
                                         sorted_seq_names):
            max_cumulative_abs = 0
            for sample_uid in self.abs_count_df.index:
                new_rect_list = []
                abs_series = self.abs_count_df.loc[sample_uid]
                rel_series = self.rel_count_df.loc[sample_uid]
                cumulative_count_abs = 0
                cumulative_count_rel = 0

                for seq in sorted_seq_names:
                    seq_abund_abs = abs_series.at[seq]
                    seq_abund_rel = rel_series.at[seq]
                    if seq_abund_abs:
                        cumulative_count_abs += int(seq_abund_abs)
                        cumulative_count_rel += float(seq_abund_rel)
                        new_rect_list.append({
                            "seq_name": seq,
                            "y_abs": cumulative_count_abs,
                            "y_rel": f'{cumulative_count_rel:.3f}',
                            "height_rel": f'{float(seq_abund_rel):.3f}',
                            "height_abs": int(seq_abund_abs),
                        })

                post_med_rect_dict[sample_uid] = new_rect_list
                if cumulative_count_abs > max_cumulative_abs:
                    max_cumulative_abs = cumulative_count_abs
            return max_cumulative_abs

        def _output_pre_med_master_fasta(self):

            fasta_out = []
            for seq_name, seq in self.master_fasta_dict.items():
                fasta_out.append(f'>{seq_name}')
                fasta_out.append(f'{seq}')

            self.parent.thread_safe_general.write_list_to_destination(self.parent.pre_med_fasta_out_path, fasta_out)

        def _write_out_dfs_as_csv(self):
            print('\nWriting out pre-med sequence count tables: \n')
            print(self.parent.pre_med_absolute_df_path)
            print(self.parent.pre_med_relative_df_path)
            self.rel_count_df.to_csv(self.parent.pre_med_relative_df_path)
            self.abs_count_df.to_csv(self.parent.pre_med_absolute_df_path)

        def _populate_dfs_in_parent_sample_order(self):
            print('\nPopulating pre-MED dfs\n')
            ordered_seq_names = list(self.master_sequence_count_dict.items())
            ordered_seq_names.sort(key=lambda x: x[0], reverse=True)
            ordered_seq_names.sort(key=lambda x: x[1], reverse=True)
            ordered_seq_names = [_[0] for _ in ordered_seq_names]

            # fastest way to create the df from the dictionaries is to create a list of dictionaries and then
            # fill in the nan values and rearrange the columns
            list_of_dicts = []
            for sample_uid in self.parent.sorted_sample_uid_list:
                abund_dict = self.master_sample_uid_to_abund_dict[sample_uid]
                temp_dict = {seq_name: abund for seq_name, abund in abund_dict.items()}
                temp_dict['sample_uid'] = sample_uid
                list_of_dicts.append(temp_dict)

            self.abs_count_df = pd.DataFrame(list_of_dicts)
            self.abs_count_df.set_index('sample_uid', inplace=True, drop=True)
            self.abs_count_df = self.abs_count_df.reindex(ordered_seq_names, axis=1).fillna(0).astype(int)

            self.rel_count_df = self.abs_count_df.div(self.abs_count_df.sum(axis=1), axis=0).fillna(0).astype(float)
            sample_names_list = [self.sample_uid_to_name_dict[sample_uid] for sample_uid in
                                 self.abs_count_df.index.values.tolist()]
            self.abs_count_df['sample_name'] = sample_names_list
            self.rel_count_df['sample_name'] = sample_names_list
            # reorder the columns of the df
            self.abs_count_df = self.abs_count_df.reindex(['sample_name'] + list(self.abs_count_df)[:-1], axis=1)
            self.rel_count_df = self.rel_count_df.reindex(['sample_name'] + list(self.rel_count_df)[:-1], axis=1)

            return ordered_seq_names

        def _count_sequences_in_each_sample(self):
            """Populate the master_sequence_count_dict that holds the cummulative relative abundance of every
            reference sequence in the pre-med seq collection across all of the samples.
            Also populate the master_sample_uid_to_abund_dict, that holds an absolute abundance dict for each
            sample where the sample.uid is the key and a temporary dictionary is the value where the key is the
            rs.name and the value is the abosulte abundance of that sequence in the sample.
            Through the process of counting, also populate the master_fasta_dict be rs_name to rs_sequence dictionary.
            This will be used to make the master fasta that will represent every sequence in the pre-med sequence
            collection."""
            print('\nCounting pre-MED sequences in DataSetSamples')
            count = 0
            num_samples = len(self.parent.list_of_dss_objects)
            for dss_obj in self.parent.list_of_dss_objects:
                count += 1
                sys.stdout.write(f'\rcounting pre-MED sequences for {dss_obj.name}: {count} out of {num_samples} samples')
                self.sample_uid_to_name_dict[dss_obj.id] = dss_obj.name
                dsspm_objs_of_sample = DataSetSampleSequencePM.objects.filter(data_set_sample_from=dss_obj)
                dss_objs_of_sample = DataSetSampleSequence.objects.filter(data_set_sample_from=dss_obj)
                # Check to see that there are DataSetSampleSequencePM associated with the sample
                # If there are no sequenecs associated with the sample, then we will assume that
                # DataSetSampleSequencePM objects were not generated during the loading of this dataset
                # and we will raise a NoDataSetSampleSequencePMObjects exception that will allow us to skip generation
                # of this output.
                # Also check that there are data set sample sequences. Else there could be no DataSetSamplePM objects
                # due to the fact that there were no symbiodiniaceae sequences in this sample.
                if len(dsspm_objs_of_sample) < 1 and len(dss_objs_of_sample) > 1:
                    raise NoDataSetSampleSequencePMObjects(f'No DataSetSampleSequence objects found')
                # total_seqs = dss_obj.non_sym_absolute_num_seqs
                sample_temp_abundance_dict = {}
                for dsspm_obj in dsspm_objs_of_sample:
                    if str(dsspm_obj.reference_sequence_of) not in self.master_fasta_dict:
                        self.master_fasta_dict[str(dsspm_obj.reference_sequence_of)] = dsspm_obj.reference_sequence_of.sequence
                    sample_temp_abundance_dict[str(dsspm_obj.reference_sequence_of)] = dsspm_obj.abundance
                    self.master_sequence_count_dict[str(dsspm_obj.reference_sequence_of)] += dsspm_obj.abundance / dss_obj.absolute_num_sym_seqs
                self.master_sample_uid_to_abund_dict[dss_obj.id] = sample_temp_abundance_dict
            print('\nPre-MED sequence counting complete')


class SequenceCountTableCollectAbundanceHandler:
    """The purpose of this handler and the associated worker is to populate three dictionaries that will be used
    in making the count table output.
    1 - dict(ref_seq_name : cumulative relative abundance for each sequence across all samples)
    2 - sample_id : list(
                         dict(ref_seq_of_sample_name:absolute_abundance_of_dsss_in_sample),
                         dict(ref_seq_of_sample_name:relative_abundance_of_dsss_in_sample)
                         )
    3 - sample_id : list(
                         dict(clade:total_abund_of_no_name_seqs_of_clade_in_q_),
                         dict(clade:relative_abund_of_no_name_seqs_of_clade_in_q_)
                         )
    Abbreviations:
    ds = DataSet
    dss = DataSetSample
    dsss = DataSetSampleSequence
    ref_seq = ReferenceSeqeunce
    The end product of this method will be returned to the count table creator. The first dict will be used to create a
    list of the ReferenceSequence objects of this output ordered first by clade and then by cumulative relative
    abundance across all samples in the output.
    """
    def __init__(self, parent_seq_count_tab_creator):
        self.seq_count_table_creator = parent_seq_count_tab_creator
        if self.seq_count_table_creator.multiprocess:
            self.mp_manager = Manager()
            self.input_dss_mp_queue = mp_Queue()
            self._populate_input_dss_mp_queue()
            self.ref_seq_names_clade_annotated = [
                ref_seq.name if ref_seq.has_name else str(ref_seq) for
                ref_seq in self.seq_count_table_creator.ref_seqs_in_datasets]
            self.dss_id_to_list_of_dsss_objects_mp_dict = self.mp_manager.dict()
            self._populate_dss_id_to_list_of_dsss_objects()
            self.annotated_dss_name_to_cummulative_rel_abund_mp_dict = self.mp_manager.dict(
            {refSeq_name: 0 for refSeq_name in self.ref_seq_names_clade_annotated})
            self.lock = mp_Lock()
            self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict = self.mp_manager.dict()
            self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict = self.mp_manager.dict()
        else:
            self.input_dss_mp_queue = mt_Queue()
            self._populate_input_dss_mp_queue()
            self.ref_seq_names_clade_annotated = [
                ref_seq.name if ref_seq.has_name else str(ref_seq) for
                ref_seq in self.seq_count_table_creator.ref_seqs_in_datasets]
            self.dss_id_to_list_of_dsss_objects_mp_dict = {}
            self._populate_dss_id_to_list_of_dsss_objects()
            self.annotated_dss_name_to_cummulative_rel_abund_mp_dict = {refSeq_name: 0 for refSeq_name in self.ref_seq_names_clade_annotated}
            self.lock = mt_Lock()
            self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict = {}
            self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict = {}

        # this is the list that we will use the self.annotated_dss_name_to_cummulative_rel_abund_mp_dict to create
        # it is a list of the ref_seqs_ordered first by clade then by abundance.
        self.clade_abundance_ordered_ref_seq_list = []

    def execute_sequence_count_table_ordered_seqs_worker(self):
        all_processes = []
        
        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        for n in range(self.seq_count_table_creator.num_proc):
            if self.seq_count_table_creator.multiprocess:
                p = Process(target=self._sequence_count_table_ordered_seqs_worker, args=(
                    self.input_dss_mp_queue, 
                    self.dss_id_to_list_of_dsss_objects_mp_dict, 
                    self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict, 
                    self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict,
                    self.annotated_dss_name_to_cummulative_rel_abund_mp_dict, self.ref_seq_names_clade_annotated, self.lock))
            else:
                p = Thread(target=self._sequence_count_table_ordered_seqs_worker, args=(
                    self.input_dss_mp_queue, 
                    self.dss_id_to_list_of_dsss_objects_mp_dict, 
                    self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict, 
                    self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict,
                    self.annotated_dss_name_to_cummulative_rel_abund_mp_dict, self.ref_seq_names_clade_annotated, self.lock))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

        self._generate_clade_abundance_ordered_ref_seq_list_from_seq_name_abund_dict()

    @staticmethod
    def _sequence_count_table_ordered_seqs_worker(
        in_q, 
        dss_id_to_list_of_dsss_objects_mp_dict, 
        dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict, 
        dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict, 
        annotated_dss_name_to_cummulative_rel_abund_mp_dict, ref_seq_names_clade_annotated, lock):
        for dss in iter(in_q.get, 'STOP'):
            sys.stdout.write(f'\r{dss.name}: collecting seq abundances')
            sequence_count_table_ordered_seqs_worker_instance = SequenceCountTableCollectAbundanceWorker(
                dss_id_to_list_of_dsss_objects_mp_dict=dss_id_to_list_of_dsss_objects_mp_dict,
                dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict=dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict,
                dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict=dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict,
                annotated_dss_name_to_cummulative_rel_abund_mp_dict=annotated_dss_name_to_cummulative_rel_abund_mp_dict, 
                ref_seq_names_clade_annotated=ref_seq_names_clade_annotated,
                dss=dss, lock=lock)
            sequence_count_table_ordered_seqs_worker_instance.start_seq_abund_collection()
    
    def _generate_clade_abundance_ordered_ref_seq_list_from_seq_name_abund_dict(self):
        for i in range(len(self.seq_count_table_creator.ordered_list_of_clades_found)):
            temp_within_clade_list_for_sorting = []
            for seq_name, abund_val in self.annotated_dss_name_to_cummulative_rel_abund_mp_dict.items():
                if seq_name.startswith(
                        self.seq_count_table_creator.ordered_list_of_clades_found[i]) or seq_name[-2:] == \
                        f'_{self.seq_count_table_creator.ordered_list_of_clades_found[i]}':
                    # then this is a seq of the clade in Q and we should add to the temp list
                    temp_within_clade_list_for_sorting.append((seq_name, abund_val))
            # now sort the temp_within_clade_list_for_sorting and add to the cladeAbundanceOrderedRefSeqList
            # We want this sort order to be constant. To enusre this we should sort by both cummulative rel abund
            # and then by the seq name
            # https://docs.python.org/3/howto/sorting.html
            # The python8 docs say to sort by secondary parameter first then the primary parameter
            temp_within_clade_list_for_sorting.sort(key=lambda x: x[0], reverse=True)
            temp_within_clade_list_for_sorting.sort(key=lambda x: x[1], reverse=True)
            sorted_within_clade = [a[0] for a in temp_within_clade_list_for_sorting]
            self.clade_abundance_ordered_ref_seq_list.extend(sorted_within_clade)

    def _populate_input_dss_mp_queue(self):
        for dss in self.seq_count_table_creator.list_of_dss_objects:
            self.input_dss_mp_queue.put(dss)

        for N in range(self.seq_count_table_creator.num_proc):
            self.input_dss_mp_queue.put('STOP')

    def _populate_dss_id_to_list_of_dsss_objects(self):
        for dss in self.seq_count_table_creator.list_of_dss_objects:
            sys.stdout.write(f'\r{dss.name}')
            self.dss_id_to_list_of_dsss_objects_mp_dict[dss.id] = list(
                DataSetSampleSequence.objects.filter(data_set_sample_from=dss))


class SequenceCountTableCollectAbundanceWorker:
    def __init__(
        self,
        dss_id_to_list_of_dsss_objects_mp_dict, 
        dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict, 
        dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict, 
        annotated_dss_name_to_cummulative_rel_abund_mp_dict, ref_seq_names_clade_annotated,
        dss, lock):
        self.dss_id_to_list_of_dsss_objects_mp_dict = dss_id_to_list_of_dsss_objects_mp_dict
        self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict = dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict
        self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict = dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict
        self.annotated_dss_name_to_cummulative_rel_abund_mp_dict = annotated_dss_name_to_cummulative_rel_abund_mp_dict
        self.ref_seq_names_clade_annotated = ref_seq_names_clade_annotated
        self.dss = dss
        self.total_abundance_of_sequences_in_sample = sum([int(a) for a in json.loads(self.dss.cladal_seq_totals)])
        self.lock = lock

    def start_seq_abund_collection(self):
        clade_summary_absolute_dict, clade_summary_relative_dict = \
            self._generate_empty_noname_seq_abund_summary_by_clade_dicts()

        smple_seq_count_aboslute_dict, smple_seq_count_relative_dict = self._generate_empty_seq_name_to_abund_dicts()

        dsss_in_sample = self.dss_id_to_list_of_dsss_objects_mp_dict[self.dss.id]

        for dsss in dsss_in_sample:
            # determine what the name of the seq will be in the output
            name_unit = self._determine_output_name_of_dsss_and_pop_noname_clade_dicts(
                clade_summary_absolute_dict, clade_summary_relative_dict, dsss)

            self._populate_abs_and_rel_abundances_for_dsss(dsss, name_unit, smple_seq_count_aboslute_dict,
                                                           smple_seq_count_relative_dict)

        self._associate_sample_abundances_to_mp_dicts(
            clade_summary_absolute_dict,
            clade_summary_relative_dict,
            smple_seq_count_aboslute_dict,
            smple_seq_count_relative_dict)

    def _associate_sample_abundances_to_mp_dicts(self, clade_summary_absolute_dict, clade_summary_relative_dict,
                                                 smple_seq_count_aboslute_dict, smple_seq_count_relative_dict):
        with self.lock:
            self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict[self.dss.id] = [smple_seq_count_aboslute_dict,
                                                                                                            smple_seq_count_relative_dict]
            self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict[self.dss.id] = [
                clade_summary_absolute_dict, clade_summary_relative_dict]

    def _populate_abs_and_rel_abundances_for_dsss(self, dsss, name_unit, smple_seq_count_aboslute_dict,
                                                  smple_seq_count_relative_dict):
        rel_abund_of_dsss = dsss.abundance / self.total_abundance_of_sequences_in_sample
        with self.lock:
            self.annotated_dss_name_to_cummulative_rel_abund_mp_dict[name_unit] += rel_abund_of_dsss
        smple_seq_count_aboslute_dict[name_unit] += dsss.abundance
        smple_seq_count_relative_dict[name_unit] += rel_abund_of_dsss

    def _determine_output_name_of_dsss_and_pop_noname_clade_dicts(
            self, clade_summary_absolute_dict, clade_summary_relative_dict, dsss):
        if not dsss.reference_sequence_of.has_name:
            name_unit = str(dsss.reference_sequence_of.id) + f'_{dsss.reference_sequence_of.clade}'
            # the clade summries are only for the noName seqs
            clade_summary_absolute_dict[dsss.reference_sequence_of.clade] += dsss.abundance
            clade_summary_relative_dict[
                dsss.reference_sequence_of.clade] += dsss.abundance / self.total_abundance_of_sequences_in_sample
        else:
            name_unit = str(dsss.reference_sequence_of)
        return name_unit

    def _generate_empty_seq_name_to_abund_dicts(self):
        smple_seq_count_aboslute_dict = {seq_name: 0 for seq_name in self.ref_seq_names_clade_annotated}
        smple_seq_count_relative_dict = {seq_name: 0 for seq_name in self.ref_seq_names_clade_annotated}
        return smple_seq_count_aboslute_dict, smple_seq_count_relative_dict

    @staticmethod
    def _generate_empty_noname_seq_abund_summary_by_clade_dicts():
        clade_summary_absolute_dict = {clade: 0 for clade in list('ABCDEFGHI')}
        clade_summary_relative_dict = {clade: 0 for clade in list('ABCDEFGHI')}
        return clade_summary_absolute_dict, clade_summary_relative_dict


class SeqOutputSeriesGeneratorHandler:
    def __init__(self, parent):
        self.seq_count_table_creator = parent
        self.output_df_header = self._create_output_df_header()
        if self.seq_count_table_creator.multiprocess:
            self.worker_manager = Manager()
            # dss.id : [pandas_series_for_absolute_abundace, pandas_series_for_absolute_abundace]
            self.dss_id_to_pandas_series_results_list_mp_dict = self.worker_manager.dict()
            self.dss_input_queue = mp_Queue()
            self._populate_dss_input_queue()
            self.lock = mp_Lock()
        else:
            # dss.id : [pandas_series_for_absolute_abundace, pandas_series_for_absolute_abundace]
            self.dss_id_to_pandas_series_results_list_mp_dict = {}
            self.dss_input_queue = mt_Queue()
            self._populate_dss_input_queue()
            self.lock = mt_Lock()

    def execute_sequence_count_table_dataframe_contructor_handler(self):
        all_processes = []
        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        sys.stdout.write('\n\nOutputting seq data\n')
        for n in range(self.seq_count_table_creator.num_proc):
            if self.seq_count_table_creator.multiprocess:
                p = Process(target=self._output_df_contructor_worker, args=(
                    self.dss_input_queue,
                    self.seq_count_table_creator.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict,
                    self.seq_count_table_creator.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict,
                    self.seq_count_table_creator.clade_abundance_ordered_ref_seq_list,
                    self.dss_id_to_pandas_series_results_list_mp_dict,
                    self.output_df_header, self.lock))
            else:
                p = Thread(target=self._output_df_contructor_worker, args=(
                self.dss_input_queue,
                self.seq_count_table_creator.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict,
                self.seq_count_table_creator.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict,
                self.seq_count_table_creator.clade_abundance_ordered_ref_seq_list,
                self.dss_id_to_pandas_series_results_list_mp_dict,
                self.output_df_header, self.lock))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

    @staticmethod
    def _output_df_contructor_worker(
        in_q, 
        dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict,
        dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict,
        clade_abundance_ordered_ref_seq_list,
        dss_id_to_pandas_series_results_list_mp_dict,
        output_df_header, lock):
        for dss in iter(in_q.get, 'STOP'):
            seq_output_series_generator_worker = SeqOutputSeriesGeneratorWorker(
                dss=dss,
                list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs=dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict[dss.id],
                list_of_abs_and_rel_abund_of_contained_dsss_dicts=dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict[dss.id],
                clade_abundance_ordered_ref_seq_list=clade_abundance_ordered_ref_seq_list,
                dss_id_to_pandas_series_results_list_mp_dict=dss_id_to_pandas_series_results_list_mp_dict,
                output_df_header=output_df_header, lock=lock)
            seq_output_series_generator_worker.make_series()

    def _populate_dss_input_queue(self):
        for dss in self.seq_count_table_creator.list_of_dss_objects:
            self.dss_input_queue.put(dss)

        for N in range(self.seq_count_table_creator.num_proc):
            self.dss_input_queue.put('STOP')

    def _create_output_df_header(self):
        header_pre = self.seq_count_table_creator.clade_abundance_ordered_ref_seq_list
        no_name_summary_strings = ['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']]
        qc_stats = [
            'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs', 'post_taxa_id_absolute_symbiodiniaceae_seqs',
            'post_taxa_id_unique_symbiodiniaceae_seqs', 'size_screening_violation_absolute',
            'size_screening_violation_unique',
            'post_taxa_id_absolute_non_symbiodiniaceae_seqs', 'post_taxa_id_unique_non_symbiodiniaceae_seqs',
            'post_med_absolute',
            'post_med_unique']
        user_supplied_stats = [
            'sample_type', 'host_phylum', 'host_class', 'host_order', 'host_family', 'host_genus', 'host_species',
            'collection_latitude', 'collection_longitude', 'collection_date', 'collection_depth']

        # we add the plus one to take account of the 'sample_name' header
        self.seq_count_table_creator.number_of_meta_cols_added = \
            len(qc_stats) + len(no_name_summary_strings) + len(user_supplied_stats) + 1
        return ['sample_name', 'fastq_fwd_file_name', 'fastq_fwd_sha256_file_hash', 'fastq_rev_file_name',
                'fastq_rev_sha256_file_hash', 'data_set_uid', 'data_set_name'] + qc_stats + no_name_summary_strings + user_supplied_stats + header_pre


class SeqOutputSeriesGeneratorWorker:
    def __init__(
            self, dss, list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs,
            list_of_abs_and_rel_abund_of_contained_dsss_dicts,
            dss_id_to_pandas_series_results_list_mp_dict,
            clade_abundance_ordered_ref_seq_list, output_df_header, lock
    ):

        self.dss = dss
        # dss.id : [{dsss:absolute abundance in dss}, {dsss:relative abundance in dss}]
        # dss.id : [{clade:total absolute abundance of no name seqs from that clade},
        #           {clade:total relative abundance of no name seqs from that clade}
        #          ]
        self.list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs = list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs
        self.list_of_abs_and_rel_abund_of_contained_dsss_dicts = list_of_abs_and_rel_abund_of_contained_dsss_dicts
        self.dss_id_to_pandas_series_results_list_mp_dict = dss_id_to_pandas_series_results_list_mp_dict
        self.clade_abundance_ordered_ref_seq_list = clade_abundance_ordered_ref_seq_list
        self.output_df_header = output_df_header
        self.sample_row_data_absolute = []
        self.sample_row_data_relative = []
        self.sample_seq_tot = sum([int(a) for a in json.loads(dss.cladal_seq_totals)])
        self.lock = lock

    def make_series(self):
        sys.stdout.write(f'\r{self.dss.name}: Creating data ouput row')
        ds = self.dss.data_submission_from
        self.sample_row_data_absolute.append(self.dss.name)
        self.sample_row_data_absolute.append(self.dss.fastq_fwd_file_name)
        self.sample_row_data_absolute.append(self.dss.fastq_fwd_file_hash)
        self.sample_row_data_absolute.append(self.dss.fastq_rev_file_name)
        self.sample_row_data_absolute.append(self.dss.fastq_rev_file_hash)
        self.sample_row_data_absolute.append(ds.id)
        self.sample_row_data_absolute.append(ds.name)

        self.sample_row_data_relative.append(self.dss.name)
        self.sample_row_data_relative.append(self.dss.fastq_fwd_file_name)
        self.sample_row_data_relative.append(self.dss.fastq_fwd_file_hash)
        self.sample_row_data_relative.append(self.dss.fastq_rev_file_name)
        self.sample_row_data_relative.append(self.dss.fastq_rev_file_hash)
        self.sample_row_data_relative.append(ds.id)
        self.sample_row_data_relative.append(ds.name)

        if self._dss_had_problem_in_processing():

            self._populate_quality_control_data_of_failed_sample()

            self._output_the_failed_sample_pandas_series()
            return

        self._populate_quality_control_data_of_successful_sample()

        self._output_the_successful_sample_pandas_series()

    def _output_the_successful_sample_pandas_series(self):
        sample_series_absolute = pd.Series(self.sample_row_data_absolute, index=self.output_df_header, name=self.dss.id)
        sample_series_relative = pd.Series(self.sample_row_data_relative, index=self.output_df_header, name=self.dss.id)
        with self.lock:
            self.dss_id_to_pandas_series_results_list_mp_dict[self.dss.id] = [
                sample_series_absolute, sample_series_relative]

    def _populate_quality_control_data_of_successful_sample(self):
        self._populate_qc_meta_successful_sample()

        self._populate_no_name_seq_clade_summaries_successful_sample()

        self._populate_user_supplied_meta()

        self._populate_seq_abunds_successful_sample()

    def _populate_seq_abunds_successful_sample(self):
        # and append these abundances in order of cladeAbundanceOrderedRefSeqList to
        # the sampleRowDataCounts and the sampleRowDataProps
        for seq_name in self.clade_abundance_ordered_ref_seq_list:
            sys.stdout.write('\rOutputting seq data for {}: sequence {}'.format(self.dss.name, seq_name))
            self.sample_row_data_absolute.append(
                self.list_of_abs_and_rel_abund_of_contained_dsss_dicts[0][seq_name])
            self.sample_row_data_relative.append(
                self.list_of_abs_and_rel_abund_of_contained_dsss_dicts[1][seq_name])

    def _populate_no_name_seq_clade_summaries_successful_sample(self):
        # now add the clade divided summaries of the clades
        for clade in list('ABCDEFGHI'):
            
            self.sample_row_data_absolute.append(
                self.list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs[0][clade])
            self.sample_row_data_relative.append(
                self.list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs[1][clade])

    def _populate_qc_meta_successful_sample(self):
        # Here we add in the post qc and post-taxa id counts
        # For the absolute counts we will report the absolute seq number
        # For the relative counts we will report these as proportions of the sampleSeqTot.
        # I.e. we will have numbers larger than 1 for many of the values and the symbiodiniaceae seqs should be 1

        # CONTIGS
        # This is the absolute number of sequences after make.contigs
        contig_num = self.dss.num_contigs
        self.sample_row_data_absolute.append(contig_num)
        self.sample_row_data_relative.append(contig_num / self.sample_seq_tot)
        # POST-QC
        # store the aboslute number of sequences after sequencing QC at this stage
        post_qc_absolute = self.dss.post_qc_absolute_num_seqs
        self.sample_row_data_absolute.append(post_qc_absolute)
        self.sample_row_data_relative.append(post_qc_absolute / self.sample_seq_tot)
        # This is the unique number of sequences after the sequencing QC
        post_qc_unique = self.dss.post_qc_unique_num_seqs
        self.sample_row_data_absolute.append(post_qc_unique)
        self.sample_row_data_relative.append(post_qc_unique / self.sample_seq_tot)
        # POST TAXA-ID
        # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
        tax_id_symbiodiniaceae_absolute = self.dss.absolute_num_sym_seqs
        self.sample_row_data_absolute.append(tax_id_symbiodiniaceae_absolute)
        self.sample_row_data_relative.append(tax_id_symbiodiniaceae_absolute / self.sample_seq_tot)
        # Same as above but the number of unique seqs
        tax_id_symbiodiniaceae_unique = self.dss.unique_num_sym_seqs
        self.sample_row_data_absolute.append(tax_id_symbiodiniaceae_unique)
        self.sample_row_data_relative.append(tax_id_symbiodiniaceae_unique / self.sample_seq_tot)
        # store the absolute number of sequences lost to size cutoff violations
        size_violation_aboslute = self.dss.size_violation_absolute
        self.sample_row_data_absolute.append(size_violation_aboslute)
        self.sample_row_data_relative.append(size_violation_aboslute / self.sample_seq_tot)
        # store the unique size cutoff violations
        size_violation_unique = self.dss.size_violation_unique
        self.sample_row_data_absolute.append(size_violation_unique)
        self.sample_row_data_relative.append(size_violation_unique / self.sample_seq_tot)
        # store the abosolute number of sequenes that were not considered Symbiodinium
        tax_id_non_symbiodinum_abosulte = self.dss.non_sym_absolute_num_seqs
        self.sample_row_data_absolute.append(tax_id_non_symbiodinum_abosulte)
        self.sample_row_data_relative.append(tax_id_non_symbiodinum_abosulte / self.sample_seq_tot)
        # This is the number of unique sequences that were not considered Symbiodinium
        tax_id_non_symbiodiniaceae_unique = self.dss.non_sym_unique_num_seqs
        self.sample_row_data_absolute.append(tax_id_non_symbiodiniaceae_unique)
        self.sample_row_data_relative.append(tax_id_non_symbiodiniaceae_unique / self.sample_seq_tot)
        # Post MED absolute
        post_med_absolute = self.dss.post_med_absolute
        self.sample_row_data_absolute.append(post_med_absolute)
        self.sample_row_data_relative.append(post_med_absolute / self.sample_seq_tot)
        # Post MED unique
        post_med_unique = self.dss.post_med_unique
        self.sample_row_data_absolute.append(post_med_unique)
        self.sample_row_data_relative.append(post_med_unique / self.sample_seq_tot)

    def _output_the_failed_sample_pandas_series(self):
        sample_series_absolute = pd.Series(self.sample_row_data_absolute, index=self.output_df_header, name=self.dss.id)
        sample_series_relative = pd.Series(self.sample_row_data_relative, index=self.output_df_header, name=self.dss.id)
        with self.lock:
            self.dss_id_to_pandas_series_results_list_mp_dict[self.dss.id] = [sample_series_absolute,
                                                                                    sample_series_relative]

    def _populate_quality_control_data_of_failed_sample(self):
        # Add in the qc totals if possible
        # For the proportions we will have to add zeros as we cannot do proportions

        self._populate_qc_meta_failed_sample()

        self._populate_no_name_seq_clade_summaries_failed_sample()

        self._populate_user_supplied_meta()

        self._populate_seq_abunds_failed_sample()

    def _populate_seq_abunds_failed_sample(self):
        # All sequences get 0s
        for _ in self.clade_abundance_ordered_ref_seq_list:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)

    def _populate_user_supplied_meta(self):
        # insert the user supplied meta stats
        # sample_type
        self.sample_row_data_absolute.append(self.dss.sample_type)
        self.sample_row_data_relative.append(self.dss.sample_type)

        # host_phylum
        self.sample_row_data_absolute.append(self.dss.host_phylum)
        self.sample_row_data_relative.append(self.dss.host_phylum)

        # host_class
        self.sample_row_data_absolute.append(self.dss.host_class)
        self.sample_row_data_relative.append(self.dss.host_class)

        # host_order
        self.sample_row_data_absolute.append(self.dss.host_order)
        self.sample_row_data_relative.append(self.dss.host_order)

        # host_family
        self.sample_row_data_absolute.append(self.dss.host_family)
        self.sample_row_data_relative.append(self.dss.host_family)

        # host_genus
        self.sample_row_data_absolute.append(self.dss.host_genus)
        self.sample_row_data_relative.append(self.dss.host_genus)

        # host_species
        self.sample_row_data_absolute.append(self.dss.host_species)
        self.sample_row_data_relative.append(self.dss.host_species)

        # collection_latitude
        self.sample_row_data_absolute.append(self.dss.collection_latitude)
        self.sample_row_data_relative.append(self.dss.collection_latitude)

        # collection_longitude
        self.sample_row_data_absolute.append(self.dss.collection_longitude)
        self.sample_row_data_relative.append(self.dss.collection_longitude)

        # collection_date
        self.sample_row_data_absolute.append(self.dss.collection_date)
        self.sample_row_data_relative.append(self.dss.collection_date)

        # collection_depth
        self.sample_row_data_absolute.append(self.dss.collection_depth)
        self.sample_row_data_relative.append(self.dss.collection_depth)


    def _populate_no_name_seq_clade_summaries_failed_sample(self):
        # no name clade summaries get 0.
        for _ in list('ABCDEFGHI'):
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)

    def _populate_qc_meta_failed_sample(self):
        if self.dss.num_contigs:
            contig_num = self.dss.num_contigs
            self.sample_row_data_absolute.append(contig_num)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # POST-QC
        # store the aboslute number of sequences after sequencing QC at this stage
        if self.dss.post_qc_absolute_num_seqs:
            post_qc_absolute = self.dss.post_qc_absolute_num_seqs
            self.sample_row_data_absolute.append(post_qc_absolute)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # This is the unique number of sequences after the sequencing QC
        if self.dss.post_qc_unique_num_seqs:
            post_qc_unique = self.dss.post_qc_unique_num_seqs
            self.sample_row_data_absolute.append(post_qc_unique)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # POST TAXA-ID
        # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
        if self.dss.absolute_num_sym_seqs:
            tax_id_symbiodiniaceae_absolute = self.dss.absolute_num_sym_seqs
            self.sample_row_data_absolute.append(tax_id_symbiodiniaceae_absolute)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # Same as above but the number of unique seqs
        if self.dss.unique_num_sym_seqs:
            tax_id_symbiodiniaceae_unique = self.dss.unique_num_sym_seqs
            self.sample_row_data_absolute.append(tax_id_symbiodiniaceae_unique)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # size violation absolute
        if self.dss.size_violation_absolute:
            size_viol_ab = self.dss.size_violation_absolute
            self.sample_row_data_absolute.append(size_viol_ab)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # size violation unique
        if self.dss.size_violation_unique:
            size_viol_uni = self.dss.size_violation_unique
            self.sample_row_data_absolute.append(size_viol_uni)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # store the abosolute number of sequenes that were not considered Symbiodinium
        if self.dss.non_sym_absolute_num_seqs:
            tax_id_non_symbiodinum_abosulte = self.dss.non_sym_absolute_num_seqs
            self.sample_row_data_absolute.append(tax_id_non_symbiodinum_abosulte)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # This is the number of unique sequences that were not considered Symbiodinium
        if self.dss.non_sym_unique_num_seqs:
            tax_id_non_symbiodiniaceae_unique = self.dss.non_sym_unique_num_seqs
            self.sample_row_data_absolute.append(tax_id_non_symbiodiniaceae_unique)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # post-med absolute
        if self.dss.post_med_absolute:
            post_med_abs = self.dss.post_med_absolute
            self.sample_row_data_absolute.append(post_med_abs)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # post-med absolute
        if self.dss.post_med_unique:
            post_med_uni = self.dss.post_med_unique
            self.sample_row_data_absolute.append(post_med_uni)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)

    def _dss_had_problem_in_processing(self):
        return self.dss.error_in_processing or self.sample_seq_tot == 0



