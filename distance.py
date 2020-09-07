import itertools
import math
import os
import subprocess
import sys
import logging
import numpy as np
import pandas as pd
from scipy.spatial.distance import braycurtis
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.tree import TreeNode
import django_general
from general import ThreadSafeGeneral
from dbApp.models import (
    ReferenceSequence, DataSetSampleSequence, AnalysisType, DataSetSample,
    CladeCollection, CladeCollectionType)
from exceptions import InsufficientSequencesInAlignment, EigenValsTooSmallError


class BaseUnifracDistPCoACreator:
    """Base class for TypeUnifracDistPCoACreator and SampleUnifracDistPCoACreator.
    These classes are used for generating UniFrac distances between either ITS2 type profiles
    or Samples."""
    def __init__(
            self, num_proc, output_dir, data_set_uid_list, js_output_path_dict,
            html_dir, data_set_sample_uid_list, cct_set_uid_list,
            date_time_str):
        self.thread_safe_general = ThreadSafeGeneral()
        self.num_proc = num_proc
        self.output_dir = output_dir
        self.data_set_sample_uid_list, self.clade_col_uid_list = self._set_data_set_sample_uid_list(
            data_set_sample_uid_list=data_set_sample_uid_list, data_set_uid_list=data_set_uid_list,
            cct_set_uid_list=cct_set_uid_list)
        self.output_path_list = []
        self.date_time_str = date_time_str
        self.clade_output_dir = None
        self.html_dir = html_dir
        # for the js dict
        self.genera_annotation_dict = {
            'A': 'Symbiodinium', 'B': 'Breviolum', 'C': 'Cladocopium', 'D': 'Durusdinium',
            'E': 'Effrenium', 'F': 'Clade F', 'G': 'Clade G', 'H': 'Clade H', 'I': 'Clade I'
        }
        self.pc_coordinates_dict_sqrt = {}
        self.pc_coordinates_dict_no_sqrt = {}
        self.pc_variances_dict_sqrt = {}
        self.pc_variances_dict_no_sqrt = {}
        self.pc_availabaility_dict_sqrt = {}
        self.pc_availabaility_dict_no_sqrt = {}
        self.js_file_path = os.path.join(self.html_dir, 'study_data.js')
        self.js_output_path_dict = js_output_path_dict

    def _set_data_set_sample_uid_list(self, data_set_sample_uid_list, data_set_uid_list, cct_set_uid_list):
        if data_set_sample_uid_list:
            data_set_samples_of_output = self._chunk_query_dss_objs_from_dss_uids(data_set_sample_uid_list)

            clade_col_uids_of_output = self._get_cc_uids_of_output_from_dss_objs(data_set_samples_of_output)

            return data_set_sample_uid_list, clade_col_uids_of_output
        elif cct_set_uid_list:
            clade_col_uids_of_output, clade_cols_of_output = self._get_distinct_cc_uids_of_output_from_cct_uids(
                cct_set_uid_list)

            data_set_sample_uid_of_output = self._get_distinct_dss_objs_uids_of_output_from_cc_objs(
                clade_cols_of_output)

            return data_set_sample_uid_of_output, clade_col_uids_of_output
        else:
            data_set_samples_of_output = self._chunk_query_dss_objs_from_ds_uids(data_set_uid_list)

            clade_col_uids_of_output = self._get_cc_uids_of_output_from_dss_objs(data_set_samples_of_output)

            data_set_sample_uid_of_output = [dss.id for dss in data_set_samples_of_output]

            return data_set_sample_uid_of_output, clade_col_uids_of_output

    def _get_distinct_dss_objs_uids_of_output_from_cc_objs(self, clade_cols_of_output):
        data_set_samples_of_output = self._chunk_query_distinct_dss_objs_from_cc_objs(clade_cols_of_output)
        data_set_sample_uid_of_output = [dss.id for dss in data_set_samples_of_output]
        return data_set_sample_uid_of_output

    def _get_distinct_cc_uids_of_output_from_cct_uids(self, cct_set_uid_list):
        clade_cols_of_output = self._chunk_query_distinct_cc_objs_from_cct_uids(cct_set_uid_list)
        clade_col_uids_of_output = [
            cc.id for cc in clade_cols_of_output]
        return clade_col_uids_of_output, clade_cols_of_output

    def _chunk_query_distinct_dss_objs_from_cc_objs(self, clade_cols_of_output):
        data_set_samples_of_output_set = set()
        for uid_list in self.thread_safe_general.chunks(clade_cols_of_output):
            data_set_samples_of_output_set.update(
                list(DataSetSample.objects.filter(cladecollection__in=uid_list)))
        return list(data_set_samples_of_output_set)

    def _chunk_query_distinct_cc_objs_from_cct_uids(self, cct_set_uid_list):
        clade_cols_of_output_set = set()
        for uid_list in self.thread_safe_general.chunks(cct_set_uid_list):
            clade_cols_of_output_set.update(list(CladeCollection.objects.filter(cladecollectiontype__id__in=uid_list)))
        return list(clade_cols_of_output_set)

    def _chunk_query_distinct_rs_objs_from_rs_uids(self, rs_uid_list):
        rs_obj_of_output_set = set()
        for rs_list in self.thread_safe_general.chunks(rs_uid_list):
            rs_obj_of_output_set.update(list(ReferenceSequence.objects.filter(id__in=rs_list)))
        return list(rs_obj_of_output_set)

    def _get_cc_uids_of_output_from_dss_objs(self, data_set_samples_of_output):
        clade_cols_of_output = self._chunk_query_cc_objs_from_dss_objs(data_set_samples_of_output)
        clade_col_uids_of_output = [
            cc.id for cc in clade_cols_of_output]
        return clade_col_uids_of_output

    def _chunk_query_cc_objs_from_dss_objs(self, data_set_samples_of_output):
        clade_cols_of_output = []
        for uid_list in self.thread_safe_general.chunks(data_set_samples_of_output):
            clade_cols_of_output.extend(list(CladeCollection.objects.filter(data_set_sample_from__in=uid_list)))
        return clade_cols_of_output

    def _chunk_query_dss_objs_from_dss_uids(self, data_set_sample_uid_list):
        data_set_samples_of_output = []
        for uid_list in self.thread_safe_general.chunks(data_set_sample_uid_list):
            data_set_samples_of_output.extend(list(DataSetSample.objects.filter(id__in=uid_list)))
        return data_set_samples_of_output

    def _chunk_query_dss_objs_from_ds_uids(self, data_set_uid_list):
        data_set_samples_of_output = []
        for uid_list in self.thread_safe_general.chunks(data_set_uid_list):
            data_set_samples_of_output.extend(list(DataSetSample.objects.filter(data_submission_from__in=uid_list)))
        return data_set_samples_of_output

    def _write_output_paths_to_stdout(self):
        print('UniFrac and PCoA computation complete. Ouput files:')
        if self.output_path_list:
            print('Output files:')
            for output_path in self.output_path_list:
                print(output_path)
        else:
            print("There are no output files")

    @staticmethod
    def _rescale_array(max_val, min_val):
        # work through the magnitudes of order and see what the bigest scaler we can work with is
        # whilst still remaining below 1
        # Return the scaler by which we should multiply
        query = 0.1
        scaler = 10
        while 1:
            if max_val > query:
                # then we cannot multiply by the scaler
                # revert back and break
                scaler /= 10
                break
            else:
                # then we can safely multiply by the scaler
                # increase by order of magnitude and test again
                # we also need to test the negative if it is negative
                if min_val < 0:
                    if min_val > (-1 * query):
                        scaler *= 10
                        query /= 10
                    else:
                        scaler /= 10
                        break
                else:
                    scaler *= 10
                    query /= 10
        # now scale the df by the scaler unless it is 1
        return scaler

    def _scale_and_compute_pcoa(self, wu):
        dist_array_scaler = self._rescale_array(max_val=wu.data.max(), min_val=wu.data.min())
        wu_data = wu.data * dist_array_scaler
        pcoa_output = pcoa(wu_data)
        # When the pcoa calculation converts very small eigen values to 0
        # In doing this, if there were not large enougher eigen values,
        # the sum of the eigen values will add to 0. This will cause a TrueDivide error.
        if pcoa_output.eigvals.sum() == 0:
            raise EigenValsTooSmallError
        pcoa_scaler = self._rescale_array(max_val=pcoa_output.samples.max().max(),
                                          min_val=pcoa_output.samples.min().min())
        pcoa_output.samples = pcoa_output.samples * pcoa_scaler
        return pcoa_output

class TypeUnifracDistPCoACreator(BaseUnifracDistPCoACreator):
    """Class for calculating the UniFrac distances between ITS2 type profiles, producing PCoA coordinates from
    these distances and plotting the resultant PCoA coordinates. These are all done on a clade by clade basis.
    Done: implement the --local flag. Also implement a --cct_uids flag which will allow users calculate distances
     for specific sets of its2 type profiles where the distances are calculated using specific sets of CladeCollections
     to work out the average relative abundances of the DIVs.
    """
    def __init__(
            self, num_processors, data_analysis_obj, js_output_path_dict, html_dir,
            output_dir, date_time_str=None, data_set_uid_list=None, data_set_sample_uid_list=None,
            cct_set_uid_list=None, local_abunds_only=False):

        super().__init__(
            num_proc=num_processors, output_dir=output_dir,
            data_set_uid_list=data_set_uid_list, data_set_sample_uid_list=data_set_sample_uid_list,
            cct_set_uid_list=cct_set_uid_list,
            date_time_str=date_time_str, js_output_path_dict=js_output_path_dict,
            html_dir=html_dir)

        self.thread_safe_general = ThreadSafeGeneral()
        self.data_analysis_obj = data_analysis_obj
        self.cct_set_uid_list = cct_set_uid_list
        if self.cct_set_uid_list is not None:
            # if working with specific profile/sample sets
            self.clade_col_type_objects = self._chunk_query_set_cct_objs_from_cct_uids()
            self.at_list_for_output = self._chunk_query_set_distinct_at_list_for_output_from_cct_uids(cct_set_uid_list)
        else:
            self.clade_col_type_objects = None
            self.at_list_for_output = self._chunk_query_set_distinct_at_list_for_output_from_dss_uids()
        self.at_id_to_at_name = {at.id: at.name for at in self.at_list_for_output}
        self.clades_for_dist_calcs = list(set([at.clade for at in self.at_list_for_output]))
        self.output_dir = os.path.join(output_dir, 'between_profile_distances')
        # whether to only use the abundances of DIVs in Types that are from the data set samples form this output only
        # i.e. rather than all instances of the type found in all samples (including samples outside of this output)
        self.local_abunds_only = local_abunds_only

    def _chunk_query_set_distinct_at_list_for_output_from_dss_uids(self):
        temp_at_set = set()
        for uid_list in self.thread_safe_general.chunks(self.data_set_sample_uid_list):
            temp_at_set.update(list(AnalysisType.objects.filter(
                data_analysis_from=self.data_analysis_obj,
                cladecollectiontype__clade_collection_found_in__data_set_sample_from__in=uid_list)))
        return list(temp_at_set)

    def _chunk_query_set_distinct_at_list_for_output_from_cct_uids(self, cct_set_uid_list):
        temp_at_set = set()
        for uid_list in self.thread_safe_general.chunks(cct_set_uid_list):
            temp_at_set.update(list(AnalysisType.objects.filter(cladecollectiontype__id__in=uid_list)))
        return list(temp_at_set)

    def _chunk_query_set_cct_objs_from_cct_uids(self):
        temp_clade_col_type_objs_list = []
        for uid_list in self.thread_safe_general.chunks(self.cct_set_uid_list):
            temp_clade_col_type_objs_list.extend(list(CladeCollectionType.objects.filter(id__in=uid_list)))
        return temp_clade_col_type_objs_list

    def compute_unifrac_dists_and_pcoa_coords(self):
        for clade_in_question in self.clades_for_dist_calcs:

            print(f'Calculating UniFrac Distances for clade: {clade_in_question}')

            analysis_type_objs_of_clade = [at for at in self.at_list_for_output if at.clade == clade_in_question]

            if len(analysis_type_objs_of_clade) < 2:
                print(f'There are insufficient samples in clade {clade_in_question} to calculate distances '
                      f'({len(analysis_type_objs_of_clade)}).'
                      f'\nNo distances will be calculated for this clade.')
                continue

            self.clade_output_dir = os.path.join(self.output_dir, clade_in_question)

            os.makedirs(self.clade_output_dir, exist_ok=True)

            # First step is to create an abundance dataframe where the sequences are in the columns
            # and the AnalysisType objects or possibly the DataSetSample objects uids are in the rows
            clade_abund_df_no_sqrt, clade_abund_df_sqrt, set_of_ref_seq_uids = self._create_profile_abundance_df(
                clade_in_question)

            try:
                tree = self._create_tree(clade_in_question, set_of_ref_seq_uids)
            except InsufficientSequencesInAlignment:
                print(f"There are insufficient sequences in clade {clade_in_question}.")
                print("Distance information will not be computed")
                continue

            if str(tree).count(')') == 1:
                print(f'There are no internal nodes on the rooted tree. '
                      f'This is likely caused by a lack of variation in the sequences used to build the tree. '
                      f'UniFrac distances cannot be calculated for clade {clade_in_question}.')
                continue

            wu_no_sqrt = None
            wu_sqrt = None
            try:
                wu_no_sqrt = self._perform_unifrac(clade_abund_df_no_sqrt, tree)
                wu_sqrt = self._perform_unifrac(clade_abund_df_sqrt, tree)
            except ValueError as e:
                if 'must be rooted' in str(e):
                    logging.error('a tree rooting error occured')
                    logging.error(f"Between profile Unifrac distance information will not be computed for clade {clade_in_question}")
                    continue

            clade_dist_file_path_no_sqrt, ordered_at_names_no_sqrt = self._write_out_dist_df(
                clade_abund_df_no_sqrt, wu_no_sqrt, clade_in_question, sqrt=False)
            clade_dist_file_path_sqrt, ordered_at_names_sqrt = self._write_out_dist_df(
                clade_abund_df_sqrt, wu_sqrt, clade_in_question, sqrt=True)

            try:
                pcoa_output_no_sqrt = self._compute_pcoa(wu_no_sqrt)
                pcoa_output_sqrt = self._compute_pcoa(wu_sqrt)
            except EigenValsTooSmallError:
                logging.error(f"The eigenvalues for the clade {clade_in_question} PCoA were too small and were "
                              f"converted to 0s by skbio's implementation of PCoA.")
                logging.error(f" Between profile Unifrac distances cannot be calculated for clade {clade_in_question}")
                continue


            clade_pcoa_file_path_no_sqrt, pcoa_coords_df_no_sqrt = self._write_out_pcoa(
                ordered_at_names_no_sqrt, pcoa_output_no_sqrt, clade_in_question, sqrt=False)
            clade_pcoa_file_path_sqrt, pcoa_coords_df_sqrt = self._write_out_pcoa(
                ordered_at_names_sqrt, pcoa_output_sqrt, clade_in_question, sqrt=True)

            self._populate_js_output_objects(clade_in_question, pcoa_coords_df_no_sqrt, sqrt=False)
            self._populate_js_output_objects(clade_in_question, pcoa_coords_df_sqrt, sqrt=True)

            self.output_path_list.extend([clade_dist_file_path_no_sqrt, clade_pcoa_file_path_no_sqrt,
                                          clade_dist_file_path_sqrt, clade_pcoa_file_path_sqrt])

            self.js_output_path_dict[
                f"btwn_profile_unifrac_{clade_in_question}_dist_no_sqrt"] = clade_dist_file_path_no_sqrt
            self.js_output_path_dict[
                f"btwn_profile_unifrac_{clade_in_question}_pcoa_no_sqrt"] = clade_pcoa_file_path_no_sqrt
            self.js_output_path_dict[
                f"btwn_profile_unifrac_{clade_in_question}_dist_sqrt"] = clade_dist_file_path_sqrt
            self.js_output_path_dict[
                f"btwn_profile_unifrac_{clade_in_question}_pcoa_sqrt"] = clade_pcoa_file_path_sqrt

        self._write_out_js_objects()
        self._write_output_paths_to_stdout()

    def _write_out_js_objects(self):
        self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
            [{'function_name': 'getBtwnProfileDistCoordsUFNoSqrt', 'python_obj': self.pc_coordinates_dict_no_sqrt},
             {'function_name': 'getBtwnProfileDistPCVariancesUFNoSqrt', 'python_obj': self.pc_variances_dict_no_sqrt},
             {'function_name': 'getBtwnProfileDistPCAvailableUFNoSqrt',
              'python_obj': self.pc_availabaility_dict_no_sqrt},
             {'function_name': 'getBtwnProfileDistCoordsUFSqrt', 'python_obj': self.pc_coordinates_dict_sqrt},
             {'function_name': 'getBtwnProfileDistPCVariancesUFSqrt', 'python_obj': self.pc_variances_dict_sqrt},
             {'function_name': 'getBtwnProfileDistPCAvailableUFSqrt', 'python_obj': self.pc_availabaility_dict_sqrt}
             ], js_outpath=self.js_file_path)

    def _populate_js_output_objects(self, clade_in_question, pcoa_coords_df, sqrt):
        # set the variance dict
        # and set the available pcs
        pcoa_coords_df.set_index('analysis_type_uid', drop=True, inplace=True)
        available_pcs = list(pcoa_coords_df)
        if len(available_pcs) > 6:
            available_pcs = available_pcs[:6]
        if sqrt:
            self.pc_availabaility_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = available_pcs
            variances = [pcoa_coords_df.iloc[-1][pc] for pc in available_pcs]
            self.pc_variances_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = variances
        else:
            self.pc_availabaility_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = available_pcs
            variances = [pcoa_coords_df.iloc[-1][pc] for pc in available_pcs]
            self.pc_variances_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = variances
        # set the coordinates data holder dict here
        genera_pc_coords_dict = {}
        for profile_uid in pcoa_coords_df.index.values.tolist()[:-1]:
            profile_pc_coords_dict = {}
            for pc in available_pcs:
                profile_pc_coords_dict[pc] = f'{pcoa_coords_df.at[profile_uid, pc]:.3f}'
            genera_pc_coords_dict[int(profile_uid)] = profile_pc_coords_dict
        if sqrt:
            self.pc_coordinates_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = genera_pc_coords_dict
        else:
            self.pc_coordinates_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = genera_pc_coords_dict

    def _write_out_pcoa(self, ordered_at_names, pcoa_output, clade_in_question, sqrt):
        # rename the pcoa dataframe index as the sample names
        pcoa_output.samples['sample'] = ordered_at_names
        renamed_pcoa_dataframe = pcoa_output.samples.set_index('sample')
        # now add the variance explained as a final row to the renamed_dataframe
        var_explained_list = [0]
        var_explained_list.extend(pcoa_output.proportion_explained.values.tolist())

        ser = pd.Series(var_explained_list, index=list(renamed_pcoa_dataframe), name='proportion_explained')
        renamed_pcoa_dataframe = renamed_pcoa_dataframe.append(ser)
        renamed_pcoa_dataframe['analysis_type_uid'] = renamed_pcoa_dataframe['analysis_type_uid'].astype(int)
        # now output the pcoa
        if sqrt:
            clade_pcoa_file_path = os.path.join(
                self.clade_output_dir,
                f'{self.date_time_str}_unifrac_profiles_PCoA_coords_{clade_in_question}_sqrt.csv')
        else:
            clade_pcoa_file_path = os.path.join(
                self.clade_output_dir,
                f'{self.date_time_str}_unifrac_profiles_PCoA_coords_{clade_in_question}_no_sqrt.csv')
        renamed_pcoa_dataframe.to_csv(header=True, index=True, path_or_buf=clade_pcoa_file_path, sep=',')
        return clade_pcoa_file_path, renamed_pcoa_dataframe

    def _compute_pcoa(self, wu):
        pcoa_output = self._scale_and_compute_pcoa(wu)
        pcoa_output.samples['analysis_type_uid'] = wu.ids
        pcoa_output.samples = pcoa_output.samples[list(pcoa_output.samples)[-1:] + list(pcoa_output.samples)[:-1]]
        return pcoa_output

    @staticmethod
    def _rescale_pcoa(pcoa_output):
        # work through the magnitudes of order and see what the bigest scaler we can work with is
        # whilst still remaining below 1
        query = 0.1
        scaler = 10
        while 1:
            if pcoa_output.samples.max().max() > query:
                # then we cannot multiply by the scaler
                # revert back and break
                scaler /= 10
                break
            else:
                # then we can safely multiply by the scaler
                # increase by order of magnitude and test again
                # we also need to test the negative if it is negative
                min_val = pcoa_output.samples.min().min()
                if min_val < 0:
                    if min_val > (-1 * query):
                        scaler *= 10
                        query /= 10
                    else:
                        scaler /= 10
                        break
                else:
                    scaler *= 10
                    query /= 10
        # now scale the df by the scaler unless it is 1
        if scaler != 1:
            pcoa_output.samples = pcoa_output.samples * scaler

    def _write_out_dist_df(self, clade_abund_df, wu, clade_in_question, sqrt):
        # get the names of the at types to ouput in the df so that the user can relate distances
        ordered_at_names = list(self.at_id_to_at_name[at_id] for at_id in clade_abund_df.index)
        # create df from the numpy 2d array
        dist_df = pd.DataFrame(data=wu.data, columns=ordered_at_names, index=ordered_at_names)
        # add in the uids of the profiles
        dist_df['profile_uid'] = clade_abund_df.index.values.tolist()
        dist_df = dist_df[list(dist_df)[-1:] + list(dist_df)[:-1]]
        # write out the df
        if sqrt:
            clade_dist_file_path = os.path.join(
                self.clade_output_dir,
                f'{self.date_time_str}_unifrac_profile_distances_{clade_in_question}_sqrt.dist')
        else:
            clade_dist_file_path = os.path.join(
                self.clade_output_dir,
                f'{self.date_time_str}_unifrac_profile_distances_{clade_in_question}_no_sqrt.dist')
        dist_df.to_csv(header=False, index=True, path_or_buf=clade_dist_file_path, sep='\t')
        return clade_dist_file_path, ordered_at_names

    @staticmethod
    def _perform_unifrac(clade_abund_df, tree):
        # perform unifrac
        print('Performing unifrac calculations')
        wu = beta_diversity(
            metric='weighted_unifrac', counts=clade_abund_df.to_numpy(),
            ids=[str(_) for _ in list(clade_abund_df.index)],
            tree=tree, otu_ids=[str(_) for _ in list(clade_abund_df.columns)])
        return wu

    def _create_tree(self, clade_in_question, set_of_ref_seq_uids):
        print(
            f'Generating phylogentic tree from {len(set_of_ref_seq_uids)} DIV its2 '
            f'sequences found in its2 type profiles of clade: {clade_in_question}')

        tree_creator = TreeCreatorForUniFrac(
            parent=self, set_of_ref_seq_uids=set_of_ref_seq_uids, clade=clade_in_question)

        tree_creator.make_tree()
        tree = tree_creator.rooted_tree
        return tree

    def _create_profile_abundance_df(self, clade_in_question):
        print(f'Generating abundance dataframe for its2 type profile in clade: {clade_in_question}')
        tadfg = self.TypeAbundanceDFGenerator(parent=self, clade=clade_in_question)
        tadfg.generate_abundance_dataframe_for_clade()
        clade_abund_df_no_sqrt = tadfg.abundance_df_no_sqrt
        clade_abund_df_sqrt = tadfg.abundance_df_sqrt
        set_of_ref_seq_uids = tadfg.reference_seq_uid_set
        return clade_abund_df_no_sqrt, clade_abund_df_sqrt, set_of_ref_seq_uids

    class TypeAbundanceDFGenerator:
        """This class will, for a given clade, produce a dataframe where the sequence uids are in the columns
        and the AnalysisType object uids are in the rows and the values are the normalised abundances
        (normalised to 10000 sequences)."""

        def __init__(self, parent, clade):
            self.parent = parent
            self.clade = clade
            self.thread_safe_general = ThreadSafeGeneral()
            # Added as part of the mcminds enhancement
            # A default dict that will have clade as key and list of reference sequence objects as value
            self.reference_seq_uid_set = set()

            # Type uid is the key and the value is a dictionary of abundances normalised to 10000
            self.seq_abundance_dict_no_sqrt = {}
            self.seq_abundance_dict_sqrt = {}
            # This will be the dataframe that we populate and can be returned from the class
            self.abundance_df_no_sqrt = None
            self.abundance_df_sqrt = None

            self.analysis_type_objs_of_clade = [at for at in self.parent.at_list_for_output if at.clade == self.clade]

        def generate_abundance_dataframe_for_clade(self):
            for at_obj in self.analysis_type_objs_of_clade:
                # the order of rows in this df is the CladeCollections in list_of_clade_collections
                # the columns is order of ordered_footprint_list (ref_seq_ids)
                df_no_sqrt = pd.DataFrame(at_obj.get_ratio_list())
                df_sqrt = pd.DataFrame(at_obj.get_ratio_list())

                ref_seq_uids_of_analysis_type = [int(b) for b in at_obj.ordered_footprint_list.split(',')]

                self.reference_seq_uid_set.update(ref_seq_uids_of_analysis_type)

                df_sqrt = self.thread_safe_general.sqrt_transform_abundance_df(df_sqrt)

                # A dictionary that will hold the abundance of the reference sequences in the given AnalysisType
                # object. Key is RefSeq_obj uid and value is abundance normalised to 10000 reads.
                # Calculate the UniFrac using only the ccts from the specified samples
                if self.parent.local_abunds_only:
                    normalised_abund_of_divs_dict_no_sqrt = self._create_norm_abund_dict_from_local_clade_cols_unifrac(
                        at_obj, df_no_sqrt, ref_seq_uids_of_analysis_type)
                    normalised_abundance_of_divs_dict_sqrt = self._create_norm_abund_dict_from_local_clade_cols_unifrac(
                        at_obj, df_sqrt, ref_seq_uids_of_analysis_type)

                # Caculate the abundance for only the specific cct that have been provided
                elif self.parent.clade_col_type_objects:
                    normalised_abund_of_divs_dict_no_sqrt = self._create_norm_abund_dict_from_cct_set_unifrac(
                        at_obj, df_no_sqrt, ref_seq_uids_of_analysis_type)
                    normalised_abundance_of_divs_dict_sqrt = self._create_norm_abund_dict_from_cct_set_unifrac(
                        at_obj, df_sqrt, ref_seq_uids_of_analysis_type)

                else:  # use all abund info to calculate av div rel abund
                    normalised_abund_of_divs_dict_no_sqrt = self._create_norm_abund_dict_from_all_clade_cols_unifrac(
                        df_no_sqrt, ref_seq_uids_of_analysis_type)
                    normalised_abundance_of_divs_dict_sqrt = self._create_norm_abund_dict_from_all_clade_cols_unifrac(
                        df_sqrt, ref_seq_uids_of_analysis_type)

                self.seq_abundance_dict_no_sqrt[at_obj.id] = normalised_abund_of_divs_dict_no_sqrt
                self.seq_abundance_dict_sqrt[at_obj.id] = normalised_abundance_of_divs_dict_sqrt

            # Here we have the self.seq_abundance_dict populated and we can use this dict of dicts to build
            # the dataframe
            self.abundance_df_no_sqrt = pd.DataFrame.from_dict(self.seq_abundance_dict_no_sqrt, orient='index')
            self.abundance_df_no_sqrt[pd.isna(self.abundance_df_no_sqrt)] = 0

            self.abundance_df_sqrt = pd.DataFrame.from_dict(self.seq_abundance_dict_sqrt, orient='index')
            self.abundance_df_sqrt[pd.isna(self.abundance_df_sqrt)] = 0

        @staticmethod
        def _create_norm_abund_dict_from_all_clade_cols_unifrac(df, ref_seq_uids_of_analysis_type):
            normalised_abundance_of_divs_dict = {
                ref_seq_uids_of_analysis_type[i]: math.ceil(df[i].mean() * 10000) for
                i in range(len(ref_seq_uids_of_analysis_type))}
            return normalised_abundance_of_divs_dict

        def _create_norm_abund_dict_from_cct_set_unifrac(self, at, df, ref_seq_uids_of_analysis_type):
            # Then we are working with a custom set of CladeCollection-AnalysisType associations. i.e. for every
            # analysis type we are working with we will have to check exactly which CladeCollections we should
            # be infering the average DIV abundances from.
            # Because the list of AnalysisTypes that we are looking through are defined by the CladeCollectionTypes
            # that we're working with, we know that there are no types in which we aren't looking at at least
            # one CladeCollection from.
            clade_collection_uid_list_of_type = [
                int(a) for a in at.list_of_clade_collections.split(',')]
            clade_collection_uids_to_output_for_this_type = []
            for cct_obj in self.parent.clade_col_type_objects:
                if cct_obj.analysis_type_of.id == at.id:
                    # then this is a CladeCollectionType of the at
                    clade_collection_uids_to_output_for_this_type.append(cct_obj.clade_collection_found_in.id)
            indices = []
            for i in range(len(clade_collection_uid_list_of_type)):
                if clade_collection_uid_list_of_type[i] in clade_collection_uids_to_output_for_this_type:
                    indices.append(i)
            normalised_abundance_of_divs_dict = {
                ref_seq_uids_of_analysis_type[i]: math.ceil(df.iloc[indices, i].mean() * 10000) for
                i in range(len(ref_seq_uids_of_analysis_type))}
            return normalised_abundance_of_divs_dict

        def _create_norm_abund_dict_from_local_clade_cols_unifrac(self, at, df, ref_seq_uids_of_analysis_type):
            # we want to limit the inference of the average relative abundance of the DIVs to only those
            # div abundance values that were found in samples from the output in question
            # as such we need to first get a list of the clade collections and see which of these
            # clade collections were from self.data_set_sample_uid_list
            clade_collection_uid_list_of_type = [
                int(a) for a in at.list_of_clade_collections.split(',')]
            # the indices of the clade collections that are of this output
            indices = []
            for i in range(len(clade_collection_uid_list_of_type)):
                if clade_collection_uid_list_of_type[i] in self.parent.clade_col_uid_list:
                    indices.append(i)
            normalised_abundance_of_divs_dict = {
                ref_seq_uids_of_analysis_type[i]: math.ceil(df.iloc[indices, i].mean() * 10000) for
                i in range(len(ref_seq_uids_of_analysis_type))}
            return normalised_abundance_of_divs_dict


class SampleUnifracDistPCoACreator(BaseUnifracDistPCoACreator):
    """

    This method will generate a distance matrix between samples using the UniFrac method.
    One for each clade.
    It will also perform a PCoA for each distance matrix. a .dist file and a .csv with the pcoa coords will be output

    I have been giving some thought as to which level we should be generating these distance matrices on.
    If we do them on the all sequence level, i.e. mixture of clades then we are going to end up with the ordination
    separations based primarily on the presence or absence of clades and this will mask the within clade differences
    Biologically it is more important to be able to tell the difference between within clade different types
    than be able to see if they contain different cladal proportions. After all, SP is about the increased resolution
    So... what I propose is that the researcher do a clade level analysis seperately to look at the differnces between
    clades and then they do ordinations clade by clade.

    Now, with the within-clade analyses we have to chose between ITS2 type profile collection based ordinations
    (i.e. only looking at the sequences contained in a clade_collection_type) or working with all of the sequnces that
    are found within a sample of the given clade. With the former, this will require that the sample has gone through
    an analysis, this is not possible for environmental samples. For coral samples this is of course possible, and we
    would then need to choose about whether you would want to have one clade_collection_type / sample pairing per clade
    collection type or whether you would plot a sample with all of its clade_collection_types plotted. I think the
    latter would be better as this would really be telling the difference between the samples rather than the types
    found within the samples. However, I'm now thinking about the cases where samples are missing one or two DIVs and
    maybe SymPortal hasn't done the best job of resolving the true type they contain. In this case the distance
    should still give a true representation of the similarity of this sample to another sample. Then you might start
    to wonder what the point of the ITS2 type profiles are. Well, they would be to tell us WHAT the types are
    whereas the distance matrix would not be able to do this but give us quantification of similarity. So they would
    in fact be nicely complementary.

    So for the time being, given the above logic, we will work on creating cladally separated distance matrices
    for samples of a given data_set or collection of dataSubmissions. For each sample we will use all sequences
    of the given clade found within the sample to calculate the distances. This will output the distance matrix
    in the outputs folder.
    """

    def __init__(
            self, num_processors, html_dir, js_output_path_dict, output_dir, date_time_str,
            data_set_uid_list=None, data_set_sample_uid_list=None):
        super().__init__(
            num_proc=num_processors, output_dir=output_dir,
            data_set_uid_list=data_set_uid_list, data_set_sample_uid_list=data_set_sample_uid_list,
            date_time_str=date_time_str, cct_set_uid_list=None, html_dir=html_dir,
            js_output_path_dict=js_output_path_dict)

        self.clade_collections_from_data_set_samples = self._chunk_query_set_cc_obj_from_dss_uids()
        self.cc_id_to_sample_name_dict = {
            cc_obj.id: cc_obj.data_set_sample_from.name for cc_obj in self.clade_collections_from_data_set_samples}
        self.cc_id_to_sample_id = {
            cc_obj.id: cc_obj.data_set_sample_from.id for cc_obj in self.clade_collections_from_data_set_samples}

        self.clades_for_dist_calcs = list(set([a.clade for a in self.clade_collections_from_data_set_samples]))

        self.output_dir = os.path.join(output_dir, 'between_sample_distances')
        os.makedirs(self.output_dir, exist_ok=True)

    def _chunk_query_set_cc_obj_from_dss_uids(self):
        temp_clade_col_objs = []
        for uid_list in self.thread_safe_general.chunks(self.data_set_sample_uid_list, 100):
            temp_clade_col_objs.extend(list(CladeCollection.objects.filter(data_set_sample_from__in=uid_list)))
        return temp_clade_col_objs

    def compute_unifrac_dists_and_pcoa_coords(self):
        for clade_in_question in self.clades_for_dist_calcs:

            print(f'Calculating UniFrac Distances for clade: {clade_in_question}')

            clade_collections_of_clade = [
                cc for cc in self.clade_collections_from_data_set_samples if cc.clade == clade_in_question]

            if len(clade_collections_of_clade) < 2:
                print(f'There are insufficient samples in clade {clade_in_question} to calculate distances '
                      f'({len(clade_collections_of_clade)}).'
                      f'\nNo distances will be calculated for this clade.')
                continue

            self.clade_output_dir = os.path.join(self.output_dir, clade_in_question)

            os.makedirs(self.clade_output_dir, exist_ok=True)

            # First step is to create an abundance dataframe where the sequences are in the columns
            # and the CladeCollection object uids are in the rows
            clade_abund_df_no_sqrt, clade_abund_df_sqrt, set_of_ref_seq_uids = self._create_sample_abundance_df(
                clade_in_question)

            try:
                tree = self._create_tree(clade_in_question, set_of_ref_seq_uids)
            except InsufficientSequencesInAlignment:
                print(f"There are insufficient sequences in clade {clade_in_question}.")
                print("Distance information will not be computed")
                continue

            if str(tree).count(')') == 1:
                print(f'There are no internal nodes on the rooted tree. '
                      f'This is likely caused by a lack of variation in the sequences used to build the tree. '
                      f'Between sample UniFrac distances cannot be calculated for clade {clade_in_question}.')
                continue

            try:
                wu_no_sqrt = self._perform_unifrac(clade_abund_df_no_sqrt, tree)
                wu_sqrt = self._perform_unifrac(clade_abund_df_sqrt, tree)
            except ValueError as e:
                if 'must be rooted' in str(e):
                    logging.error('a tree rooting error occured')
                    logging.error(f"Distance information will not be computed for clade {clade_in_question}")
                    continue

            clade_dist_file_path_no_sqrt, ordered_sample_names_no_sqrt = self._write_out_dist_df(
                clade_abund_df_no_sqrt, wu_no_sqrt, clade_in_question, sqrt=False)
            clade_dist_file_path_sqrt, ordered_sample_names_sqrt = self._write_out_dist_df(
                clade_abund_df_sqrt, wu_sqrt, clade_in_question, sqrt=True)

            try:
                pcoa_output_no_sqrt = self._compute_pcoa(wu_no_sqrt)
                pcoa_output_sqrt = self._compute_pcoa(wu_sqrt)
            except EigenValsTooSmallError:
                logging.error(f"The eigenvalues for the clade {clade_in_question} PCoA were too small and were "
                              f"converted to 0s by skbio's implementation of PCoA.")
                logging.error(f"Between sample Unifrac Distances cannot be calculated for clade {clade_in_question}")
                continue

            clade_pcoa_file_path_no_sqrt, pcoa_coords_df_no_sqrt = self._write_out_pcoa(
                ordered_sample_names_no_sqrt, pcoa_output_no_sqrt, clade_in_question, sqrt=False)
            clade_pcoa_file_path_sqrt, pcoa_coords_df_sqrt = self._write_out_pcoa(
                ordered_sample_names_sqrt, pcoa_output_sqrt, clade_in_question, sqrt=True)

            self._populate_js_output_objects(clade_in_question, pcoa_coords_df_no_sqrt, sqrt=False)
            self._populate_js_output_objects(clade_in_question, pcoa_coords_df_sqrt, sqrt=True)

            self.output_path_list.extend([clade_dist_file_path_no_sqrt, clade_pcoa_file_path_no_sqrt,
                                          clade_dist_file_path_sqrt, clade_pcoa_file_path_sqrt])

            self.js_output_path_dict[
                f"btwn_sample_unifrac_{clade_in_question}_dist_no_sqrt"] = clade_dist_file_path_no_sqrt
            self.js_output_path_dict[
                f"btwn_sample_unifrac_{clade_in_question}_pcoa_no_sqrt"] = clade_pcoa_file_path_no_sqrt
            self.js_output_path_dict[
                f"btwn_sample_unifrac_{clade_in_question}_dist_sqrt"] = clade_dist_file_path_sqrt
            self.js_output_path_dict[
                f"btwn_sample_unifrac_{clade_in_question}_pcoa_sqrt"] = clade_pcoa_file_path_sqrt

        self._write_out_js_objects()
        self._write_output_paths_to_stdout()

    def _write_out_js_objects(self):
        self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
            [{'function_name': 'getBtwnSampleDistCoordsUFNoSqrt', 'python_obj': self.pc_coordinates_dict_no_sqrt},
             {'function_name': 'getBtwnSampleDistPCVariancesUFNoSqrt', 'python_obj': self.pc_variances_dict_no_sqrt},
             {'function_name': 'getBtwnSampleDistPCAvailableUFNoSqrt',
              'python_obj': self.pc_availabaility_dict_no_sqrt},
             {'function_name': 'getBtwnSampleDistCoordsUFSqrt', 'python_obj': self.pc_coordinates_dict_sqrt},
             {'function_name': 'getBtwnSampleDistPCVariancesUFSqrt', 'python_obj': self.pc_variances_dict_sqrt},
             {'function_name': 'getBtwnSampleDistPCAvailableUFSqrt', 'python_obj': self.pc_availabaility_dict_sqrt}
             ], js_outpath=self.js_file_path)

    def _populate_js_output_objects(self, clade_in_question, pcoa_coords_df, sqrt):
        # set the variance dict
        # and set the available pcs
        pcoa_coords_df.set_index('sample_uid', drop=True, inplace=True)
        available_pcs = list(pcoa_coords_df)
        if len(available_pcs) > 6:
            available_pcs = available_pcs[:6]
        if sqrt:
            self.pc_availabaility_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = available_pcs
            variances = [pcoa_coords_df.iloc[-1][pc] for pc in available_pcs]
            self.pc_variances_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = variances
        else:
            self.pc_availabaility_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = available_pcs
            variances = [pcoa_coords_df.iloc[-1][pc] for pc in available_pcs]
            self.pc_variances_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = variances

        # set the coordinates data holder dict here
        genera_pc_coords_dict = {}
        for sample_uid in pcoa_coords_df.index.values.tolist()[:-1]:
            sample_pc_coords_dict = {}
            for pc in available_pcs:
                sample_pc_coords_dict[pc] = f'{pcoa_coords_df.at[sample_uid, pc]:.3f}'
            genera_pc_coords_dict[int(sample_uid)] = sample_pc_coords_dict

        if sqrt:
            self.pc_coordinates_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = genera_pc_coords_dict
        else:
            self.pc_coordinates_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = genera_pc_coords_dict

    def _write_out_pcoa(self, ordered_sample_names, pcoa_output, clade_in_question, sqrt):
        # rename the pcoa dataframe index as the sample names
        pcoa_output.samples['sample'] = ordered_sample_names
        renamed_pcoa_dataframe = pcoa_output.samples.set_index('sample')
        # now add the variance explained as a final row to the renamed_dataframe
        var_explained_list = [0]
        var_explained_list.extend(pcoa_output.proportion_explained.values.tolist())
        ser = pd.Series(var_explained_list, index=list(renamed_pcoa_dataframe), name='proportion_explained')
        renamed_pcoa_dataframe = renamed_pcoa_dataframe.append(ser)
        renamed_pcoa_dataframe['sample_uid'] = renamed_pcoa_dataframe['sample_uid'].astype(int)
        # now output the pcoa
        if sqrt:
            clade_pcoa_file_path = os.path.join(
                self.clade_output_dir,
                f'{self.date_time_str}_unifrac_sample_PCoA_coords_{clade_in_question}_sqrt.csv')
        else:
            clade_pcoa_file_path = os.path.join(
                self.clade_output_dir,
                f'{self.date_time_str}_unifrac_sample_PCoA_coords_{clade_in_question}_no_sqrt.csv')
        renamed_pcoa_dataframe.to_csv(header=True, index=True, path_or_buf=clade_pcoa_file_path, sep=',')
        return clade_pcoa_file_path, renamed_pcoa_dataframe

    def _compute_pcoa(self, wu):
        pcoa_output = self._scale_and_compute_pcoa(wu)
        pcoa_output.samples['sample_uid'] = [self.cc_id_to_sample_id[int(cc_id)] for cc_id in wu.ids]
        pcoa_output.samples = pcoa_output.samples[list(pcoa_output.samples)[-1:] + list(pcoa_output.samples)[:-1]]
        return pcoa_output

    def _write_out_dist_df(self, clade_abund_df, wu, clade_in_question, sqrt):
        # get the names of the samples that contained the CladeCollections
        # to ouput in the df so that the user can relate distances
        ordered_sample_names = list(self.cc_id_to_sample_name_dict[cc_uid] for cc_uid in clade_abund_df.index)
        ordered_sample_uids = list(self.cc_id_to_sample_id[cc_uid] for cc_uid in clade_abund_df.index)
        # create df from the numpy 2d array
        dist_df = pd.DataFrame(data=wu.data, columns=ordered_sample_names, index=ordered_sample_names)
        dist_df['sample_uid'] = ordered_sample_uids
        dist_df = dist_df[list(dist_df)[-1:] + list(dist_df)[:-1]]
        # write out the df
        if sqrt:
            clade_dist_file_path = os.path.join(
                self.clade_output_dir,
                f'{self.date_time_str}_unifrac_sample_distances_{clade_in_question}_sqrt.dist')
        else:
            clade_dist_file_path = os.path.join(
                self.clade_output_dir,
                f'{self.date_time_str}_unifrac_sample_distances_{clade_in_question}_no_sqrt.dist')
        dist_df.to_csv(header=False, index=True, path_or_buf=clade_dist_file_path, sep='\t')
        return clade_dist_file_path, ordered_sample_names

    @staticmethod
    def _perform_unifrac(clade_abund_df, tree):
        print('Performing unifrac calculations')
        wu = beta_diversity(
            metric='weighted_unifrac', counts=clade_abund_df.to_numpy(),
            ids=[str(_) for _ in list(clade_abund_df.index)],
            tree=tree, otu_ids=[str(_) for _ in list(clade_abund_df.columns)])
        return wu

    def _create_tree(self, clade_in_question, set_of_ref_seq_uids):
        print(
            f'Generating phylogentic tree from {len(set_of_ref_seq_uids)} its2 '
            f'sequences found in CladeCollections of clade: {clade_in_question}')
        tree_creator = TreeCreatorForUniFrac(
            parent=self, set_of_ref_seq_uids=set_of_ref_seq_uids, clade=clade_in_question)
        tree_creator.make_tree()
        tree = tree_creator.rooted_tree
        return tree

    def _create_sample_abundance_df(self, clade_in_question):
        print(f'Generating abundance dataframe for its2 type profile in clade: {clade_in_question}')
        sadfg = self.SampleAbundanceDFGenerator(parent=self, clade=clade_in_question)
        sadfg.generate_abundance_dataframe_for_clade()
        clade_abund_df_no_sqrt = sadfg.abundance_df_no_sqrt
        clade_abund_df_sqrt = sadfg.abundance_df_sqrt
        set_of_ref_seq_uids = sadfg.reference_seq_uid_set
        return clade_abund_df_no_sqrt, clade_abund_df_sqrt, set_of_ref_seq_uids

    class SampleAbundanceDFGenerator:
        """This class will, for a given clade, produce a dataframe where the sequence uids are in the columns
        and the CladeCollection object uids are in the rows and the values are the normalised abundances
        (normalised to 10000 sequences)."""
        def __init__(self, parent, clade):
            self.parent = parent
            self.clade = clade
            # Added as part of the mcminds enhancement
            # A default dict that will have clade as key and list of reference sequence objects as value
            self.reference_seq_uid_set = set()

            # Type uid is the key and the value is a dictionary of abundances normalised to 10000
            self.seq_abundance_dict_no_sqrt = {}
            self.seq_abundance_dict_sqrt = {}
            # This will be the dataframe that we populate and can be returned from the class
            self.abundance_df_no_sqrt = None
            self.abundance_df_sqrt = None

            self.clade_collections_of_clade = [
                cc for cc in self.parent.clade_collections_from_data_set_samples if cc.clade == self.clade]

        def generate_abundance_dataframe_for_clade(self):
            for cc_obj in self.clade_collections_of_clade:
                list_of_dsss_in_cc = list(DataSetSampleSequence.objects.filter(
                    clade_collection_found_in=cc_obj))
                self.reference_seq_uid_set.update([dsss.reference_sequence_of.id for dsss in list_of_dsss_in_cc])
                normalised_abund_dict_sqrt = self._make_norm_abund_dict_sqrt(list_of_dsss_in_cc)
                normalised_abund_dict_no_sqrt = self._make_norm_abund_dict_no_sqrt(list_of_dsss_in_cc)
                self.seq_abundance_dict_no_sqrt[cc_obj.id] = normalised_abund_dict_no_sqrt
                self.seq_abundance_dict_sqrt[cc_obj.id] = normalised_abund_dict_sqrt

            # Here we have the self.seq_abundance_dict populated and we can use this dict of dicts to build
            # the dataframe
            self.abundance_df_no_sqrt = pd.DataFrame.from_dict(self.seq_abundance_dict_no_sqrt, orient='index')
            self.abundance_df_no_sqrt[pd.isna(self.abundance_df_no_sqrt)] = 0
            self.abundance_df_sqrt = pd.DataFrame.from_dict(self.seq_abundance_dict_sqrt, orient='index')
            self.abundance_df_sqrt[pd.isna(self.abundance_df_sqrt)] = 0

        @staticmethod
        def _make_norm_abund_dict_no_sqrt(list_of_dsss_in_cc, normalisation_sequencing_depth=10000):
            total_seqs_of_cc = sum([dss.abundance for dss in list_of_dsss_in_cc])
            normalised_abund_dict = {
                dsss.reference_sequence_of.id: int(
                    (dsss.abundance / total_seqs_of_cc) * normalisation_sequencing_depth) for
                dsss in list_of_dsss_in_cc}
            return normalised_abund_dict

        @staticmethod
        def _make_norm_abund_dict_sqrt(list_of_dsss_in_cc, normalisation_sequencing_depth=10000):
            total_seqs_of_cc = sum([dss.abundance for dss in list_of_dsss_in_cc])
            rel_abund_dict = {dsss.id: (dsss.abundance / total_seqs_of_cc) for
                              dsss in list_of_dsss_in_cc}
            dsss_uid_to_sqrt_rel_abund_dict = {
                dsss_uid: math.sqrt(rel_abund) for dsss_uid, rel_abund in rel_abund_dict.items()}

            sqr_total = sum(dsss_uid_to_sqrt_rel_abund_dict.values())

            normalised_abund_dict = {
                dsss.reference_sequence_of.id: int(
                    (dsss_uid_to_sqrt_rel_abund_dict[dsss.id] / sqr_total) * normalisation_sequencing_depth)
                for dsss in list_of_dsss_in_cc}

            return normalised_abund_dict


class TreeCreatorForUniFrac:
    """Class responsible for generating a tree using iqtree to use in the calculation of weighted unifrac
    distances for both between sample and between its2 type profile sequences."""
    def __init__(self, parent, set_of_ref_seq_uids, clade):
        self.parent = parent
        self.clade = clade
        self.ref_seq_objs = self.parent._chunk_query_distinct_rs_objs_from_rs_uids(rs_uid_list=set_of_ref_seq_uids)
        self.num_seqs = len(self.ref_seq_objs)
        self.fasta_unaligned_path = os.path.join(
            self.parent.clade_output_dir, f'clade_{self.clade}_seqs.unaligned.fasta')
        self.fasta_aligned_path = self.fasta_unaligned_path.replace('unaligned', 'aligned')
        # self.iqtree = local['iqtree']
        self.tree_out_path_unrooted = self.fasta_aligned_path + '.treefile'
        self.tree_out_path_rooted = self.tree_out_path_unrooted.replace('.treefile', '.rooted.treefile')
        self.rooted_tree = None
        self.thread_safe_general = ThreadSafeGeneral()

    def make_tree(self):

        # write out the sequences unaligned
        print(f'Writing out {self.num_seqs} unaligned sequences')
        self._write_out_unaligned_seqs()

        if len(self.thread_safe_general.read_defined_file_to_list(self.fasta_unaligned_path)) < 5:
            raise InsufficientSequencesInAlignment

        # align the sequences
        print(f'Aligning {self.num_seqs} sequences')
        self.thread_safe_general.mafft_align_fasta(
            input_path=self.fasta_unaligned_path, output_path=self.fasta_aligned_path,
            method='unifrac', num_proc=self.parent.num_proc)

        # make the tree
        print('Testing models and making phylogenetic tree')
        print('This could take some time...')
        subprocess.run(
            ['iqtree', '-T', 'AUTO', '--threads-max', '2', '-s', f'{self.fasta_aligned_path}'])

        # root the tree
        print('Tree creation complete')
        print('Rooting the tree at midpoint')
        self.rooted_tree = TreeNode.read(self.tree_out_path_unrooted).root_at_midpoint()
        self.rooted_tree.write(self.tree_out_path_rooted)

    def _write_out_unaligned_seqs(self):
        django_general.write_ref_seq_objects_to_fasta(
            path=self.fasta_unaligned_path, list_of_ref_seq_objs=self.ref_seq_objs, identifier='id')


# BrayCurtis classes
class BaseBrayCurtisDistPCoACreator:
    def __init__(self, date_time_str, profiles_or_samples, js_output_path_dict, html_dir):
        self.date_time_str = date_time_str
        self.output_path_list = []
        self.clade_output_dir = None
        # path to the .csv file that will hold the PCoA coordinates
        self.clade_pcoa_coord_file_path_no_sqrt = None
        self.clade_pcoa_coord_file_path_sqrt = None
        # path to the .dist file that holds the unifrac or braycurtis derived sample clade-separated paired distances
        self.clade_dist_file_path_no_sqrt = None
        self.clade_dist_file_path_sqrt = None
        # The actual contents of the distance file in list form
        self.clade_dist_file_as_list_no_sqrt = None
        self.clade_dist_file_as_list_sqrt = None
        # either 'profiles' or 'samples'
        self.profiles_or_samples = profiles_or_samples

        # clade specific attributes. Will be updated for every clade processe
        self.objs_of_clade = None
        self.clade_rs_uid_to_normalised_abund_clade_dict_sqrt = {}
        self.clade_rs_uid_to_normalised_abund_clade_dict_no_sqrt = {}
        self.clade_within_clade_distances_dict_no_sqrt = {}
        self.clade_within_clade_distances_dict_sqrt = {}
        self.js_output_path_dict = js_output_path_dict
        self.html_dir = html_dir
        self.genera_annotation_dict = {
            'A': 'Symbiodinium', 'B': 'Breviolum', 'C': 'Cladocopium', 'D': 'Durusdinium',
            'E': 'Effrenium', 'F': 'Clade F', 'G': 'Clade G', 'H': 'Clade H', 'I': 'Clade I'
        }
        self.pc_coordinates_dict_no_sqrt = {}
        self.pc_coordinates_dict_sqrt = {}
        self.pc_variances_dict_no_sqrt = {}
        self.pc_variances_dict_sqrt = {}
        self.pc_availabaility_dict_no_sqrt = {}
        self.pc_availabaility_dict_sqrt = {}
        self.js_file_path = os.path.join(self.html_dir, 'study_data.js')
        self.thread_safe_general = ThreadSafeGeneral()

    def _compute_pcoa_coords(self, clade, sqrt):
        # simultaneously grab the sample names in the order of the distance matrix and put the matrix into
        # a twoD list and then convert to a numpy array
        if sqrt:
            self.clade_pcoa_coord_file_path_sqrt = os.path.join(
                self.clade_output_dir,
                f'{self.date_time_str}_braycurtis_{self.profiles_or_samples}_PCoA_coords_{clade}_sqrt.csv')
            raw_dist_file = self.thread_safe_general.read_defined_file_to_list(self.clade_dist_file_path_sqrt)
        else:
            self.clade_pcoa_coord_file_path_no_sqrt = os.path.join(
                self.clade_output_dir,
                f'{self.date_time_str}_braycurtis_{self.profiles_or_samples}_PCoA_coords_{clade}_no_sqrt.csv')
            raw_dist_file = self.thread_safe_general.read_defined_file_to_list(self.clade_dist_file_path_no_sqrt)

        temp_two_d_list = []
        object_names_from_dist_matrix = []
        object_ids_from_dist_matrix = []
        for line in raw_dist_file:
            temp_elements = line.split('\t')
            object_names_from_dist_matrix.append(temp_elements[0].replace(' ', ''))
            object_ids_from_dist_matrix.append(int(temp_elements[1]))
            temp_two_d_list.append([float(a) for a in temp_elements[2:]])

        dist_as_np_array = np.array(temp_two_d_list)

        sys.stdout.write('\rcalculating PCoA coordinates')

        dist_array_scaler = self._rescale_array(max_val=dist_as_np_array.max(), min_val=dist_as_np_array.min())

        dist_as_np_array = dist_as_np_array * dist_array_scaler

        pcoa_output = pcoa(dist_as_np_array)
        # When the pcoa calculation converts very small eigen values to 0
        # In doing this, if there were not large enougher eigen values,
        # the sum of the eigen values will add to 0. This will cause a TrueDivide error.
        if pcoa_output.eigvals.sum() == 0:
            raise EigenValsTooSmallError
        pcoa_scaler = self._rescale_array(max_val=pcoa_output.samples.max().max(),
                                          min_val=pcoa_output.samples.min().min())
        pcoa_output.samples = pcoa_output.samples * pcoa_scaler

        # rename the pcoa dataframe index as the sample names
        pcoa_output.samples['sample'] = object_names_from_dist_matrix
        renamed_pcoa_dataframe = pcoa_output.samples.set_index('sample')

        # now add the variance explained as a final row to the renamed_dataframe
        renamed_pcoa_dataframe = renamed_pcoa_dataframe.append(
            pcoa_output.proportion_explained.rename('proportion_explained'))

        object_ids_from_dist_matrix.append(0)
        if self.profiles_or_samples == 'samples':
            renamed_pcoa_dataframe.insert(loc=0, column='sample_uid', value=object_ids_from_dist_matrix)
            renamed_pcoa_dataframe['sample_uid'] = renamed_pcoa_dataframe['sample_uid'].astype('int')
        else:  # 'profiles'
            renamed_pcoa_dataframe.insert(loc=0, column='analysis_type_uid', value=object_ids_from_dist_matrix)
            renamed_pcoa_dataframe['analysis_type_uid'] = renamed_pcoa_dataframe['analysis_type_uid'].astype('int')

        if sqrt:
            renamed_pcoa_dataframe.to_csv(self.clade_pcoa_coord_file_path_sqrt, index=True, header=True, sep=',')
        else:
            renamed_pcoa_dataframe.to_csv(self.clade_pcoa_coord_file_path_no_sqrt, index=True, header=True, sep=',')
        return renamed_pcoa_dataframe

    @staticmethod
    def _rescale_array(max_val, min_val):
        # work through the magnitudes of order and see what the bigest scaler we can work with is
        # whilst still remaining below 1
        # Return the scaler by which we should multiply
        query = 0.1
        scaler = 10
        while 1:
            if max_val > query:
                # then we cannot multiply by the scaler
                # revert back and break
                scaler /= 10
                break
            else:
                # then we can safely multiply by the scaler
                # increase by order of magnitude and test again
                # we also need to test the negative if it is negative
                if min_val < 0:
                    if min_val > (-1 * query):
                        scaler *= 10
                        query /= 10
                    else:
                        scaler /= 10
                        break
                else:
                    scaler *= 10
                    query /= 10
        # now scale the df by the scaler unless it is 1
        return scaler

    def _set_data_set_sample_uid_list(self, data_set_sample_uid_list, data_set_uid_list, cct_set_uid_list):
        if data_set_sample_uid_list:
            data_set_samples_of_output = self._chunk_query_dss_objs_from_dss_uids(data_set_sample_uid_list)

            clade_col_uids_of_output = self._get_cc_uids_of_output_from_dss_objs(data_set_samples_of_output)

            return data_set_sample_uid_list, clade_col_uids_of_output
        elif cct_set_uid_list:
            clade_col_uids_of_output, clade_cols_of_output = self._get_distinct_cc_uids_of_output_from_cct_uids(
                cct_set_uid_list)

            data_set_sample_uid_of_output = self._get_distinct_dss_objs_uids_of_output_from_cc_objs(
                clade_cols_of_output)

            return data_set_sample_uid_of_output, clade_col_uids_of_output
        else:
            data_set_samples_of_output = self._chunk_query_dss_objs_from_ds_uids(data_set_uid_list)

            clade_col_uids_of_output = self._get_cc_uids_of_output_from_dss_objs(data_set_samples_of_output)

            data_set_sample_uid_of_output = [dss.id for dss in data_set_samples_of_output]

            return data_set_sample_uid_of_output, clade_col_uids_of_output

    def _get_distinct_dss_objs_uids_of_output_from_cc_objs(self, clade_cols_of_output):
        data_set_samples_of_output = self._chunk_query_distinct_dss_objs_from_cc_objs(clade_cols_of_output)
        data_set_sample_uid_of_output = [dss.id for dss in data_set_samples_of_output]
        return data_set_sample_uid_of_output

    def _get_distinct_cc_uids_of_output_from_cct_uids(self, cct_set_uid_list):
        clade_cols_of_output = self._chunk_query_distinct_cc_objs_from_cct_uids(cct_set_uid_list)
        clade_col_uids_of_output = [
            cc.id for cc in clade_cols_of_output]
        return clade_col_uids_of_output, clade_cols_of_output

    def _chunk_query_distinct_dss_objs_from_cc_objs(self, clade_cols_of_output):
        data_set_samples_of_output_set = set()
        for uid_list in self.thread_safe_general.chunks(clade_cols_of_output):
            data_set_samples_of_output_set.update(
                list(DataSetSample.objects.filter(cladecollection__in=uid_list)))
        return list(data_set_samples_of_output_set)

    def _chunk_query_distinct_cc_objs_from_cct_uids(self, cct_set_uid_list):
        clade_cols_of_output_set = set()
        for uid_list in self.thread_safe_general.chunks(cct_set_uid_list):
            clade_cols_of_output_set.update(list(CladeCollection.objects.filter(cladecollectiontype__id__in=uid_list)))
        return list(clade_cols_of_output_set)

    def _get_cc_uids_of_output_from_dss_objs(self, data_set_samples_of_output):
        clade_cols_of_output = self._chunk_query_cc_objs_from_dss_objs(data_set_samples_of_output)
        clade_col_uids_of_output = [
            cc.id for cc in clade_cols_of_output]
        return clade_col_uids_of_output

    def _chunk_query_cc_objs_from_dss_objs(self, data_set_samples_of_output):
        clade_cols_of_output = []
        for uid_list in self.thread_safe_general.chunks(data_set_samples_of_output):
            clade_cols_of_output.extend(list(CladeCollection.objects.filter(data_set_sample_from__in=uid_list)))
        return clade_cols_of_output

    def _chunk_query_dss_objs_from_dss_uids(self, data_set_sample_uid_list):
        data_set_samples_of_output = []
        for uid_list in self.thread_safe_general.chunks(data_set_sample_uid_list):
            data_set_samples_of_output.extend(list(DataSetSample.objects.filter(id__in=uid_list)))
        return data_set_samples_of_output

    def _chunk_query_dss_objs_from_ds_uids(self, data_set_uid_list):
        data_set_samples_of_output = []
        for uid_list in self.thread_safe_general.chunks(data_set_uid_list):
            data_set_samples_of_output.extend(list(DataSetSample.objects.filter(data_submission_from__in=uid_list)))
        return data_set_samples_of_output

    def _compute_braycurtis_btwn_obj_pairs(self, sqrt):
        for obj_one, obj_two in itertools.combinations(list(self.objs_of_clade), 2):
            # let's work to a virtual subsample of 100 000
            if sqrt:
                obj_one_seq_rel_abundance_dict = self.clade_rs_uid_to_normalised_abund_clade_dict_sqrt[
                    obj_one.id]
                obj_two_seq_rel_abundance_dict = self.clade_rs_uid_to_normalised_abund_clade_dict_sqrt[
                    obj_two.id]
            else:
                obj_one_seq_rel_abundance_dict = self.clade_rs_uid_to_normalised_abund_clade_dict_no_sqrt[
                    obj_one.id]
                obj_two_seq_rel_abundance_dict = self.clade_rs_uid_to_normalised_abund_clade_dict_no_sqrt[
                    obj_two.id]

            # for each comparison. Get a set of all of the sequence and convert this to a list.
            set_of_rs_uids = set(list(obj_one_seq_rel_abundance_dict.keys()))
            set_of_rs_uids.update(list(obj_two_seq_rel_abundance_dict.keys()))
            list_of_rs_uids = list(set_of_rs_uids)

            # then iter through the list to get the rel abundances for each of the samples, putting 0 if not found in
            # the sample.
            rs_abundance_list_one = []
            rs_abundance_list_two = []
            for rs_uid in list_of_rs_uids:
                # populate the abundance list cc one
                if rs_uid in obj_one_seq_rel_abundance_dict:
                    rs_abundance_list_one.append(obj_one_seq_rel_abundance_dict[rs_uid])
                else:
                    rs_abundance_list_one.append(0)

                # populate the abundance list cc two
                if rs_uid in obj_two_seq_rel_abundance_dict:
                    rs_abundance_list_two.append(obj_two_seq_rel_abundance_dict[rs_uid])
                else:
                    rs_abundance_list_two.append(0)

            # once you have this we should simply be able to crunch the bray-curtis.
            distance = braycurtis(rs_abundance_list_one, rs_abundance_list_two)

            # these distances can be stored in a dictionary by the 'id1/id2' and 'id2/id1'
            if sqrt:
                self.clade_within_clade_distances_dict_sqrt[frozenset({obj_one.id, obj_two.id})] = distance
            else:
                self.clade_within_clade_distances_dict_no_sqrt[frozenset({obj_one.id, obj_two.id})] = distance

    def _generate_distance_file(self, sqrt):
        if sqrt:
            self.clade_dist_file_as_list_sqrt = []
        else:
            self.clade_dist_file_as_list_no_sqrt = []
        for obj_outer in self.objs_of_clade:
            temp_at_string = [obj_outer.id]

            for obj_inner in self.objs_of_clade:
                if obj_outer == obj_inner:
                    temp_at_string.append(0)
                else:
                    if sqrt:
                        temp_at_string.append(
                            self.clade_within_clade_distances_dict_sqrt[frozenset({obj_outer.id, obj_inner.id})])
                    else:
                        temp_at_string.append(
                            self.clade_within_clade_distances_dict_no_sqrt[frozenset({obj_outer.id, obj_inner.id})])
            if sqrt:
                self.clade_dist_file_as_list_sqrt.append(
                    '\t'.join([str(distance_item) for distance_item in temp_at_string]))
            else:
                self.clade_dist_file_as_list_no_sqrt.append(
                    '\t'.join([str(distance_item) for distance_item in temp_at_string]))

    def _add_obj_uids_to_dist_file_and_write(self, sqrt):
        # for the output version lets also append the sample name to each line so that we can see which sample it is
        # it is important that we otherwise work eith the sample ID as the sample names may not be unique.
        if sqrt:
            dist_with_obj_name = []
            list_of_obj_uids = [int(line.split('\t')[0]) for line in self.clade_dist_file_as_list_sqrt]
        else:
            dist_with_obj_name = []
            list_of_obj_uids = [int(line.split('\t')[0]) for line in self.clade_dist_file_as_list_no_sqrt]

        if self.profiles_or_samples == 'samples':
            objs_of_outputs = self._chunk_query_dss_objs_from_dss_uids(list_of_obj_uids)
            dict_of_obj_id_to_obj_name = {obj.id: obj.name for obj in objs_of_outputs}
        else:  # 'profiles'
            objs_of_outputs = self._chunk_query_at_obj_from_at_uids(list_of_obj_uids)
            dict_of_obj_id_to_obj_name = {obj.id: obj.name for obj in objs_of_outputs}

        if sqrt:
            for line in self.clade_dist_file_as_list_sqrt:
                self._append_obj_name_to_dist_line(dict_of_obj_id_to_obj_name, dist_with_obj_name, line)
            self.clade_dist_file_as_list_sqrt = dist_with_obj_name
            self.thread_safe_general.write_list_to_destination(self.clade_dist_file_path_sqrt, self.clade_dist_file_as_list_sqrt)
        else:
            for line in self.clade_dist_file_as_list_no_sqrt:
                self._append_obj_name_to_dist_line(dict_of_obj_id_to_obj_name, dist_with_obj_name, line)
            self.clade_dist_file_as_list_no_sqrt = dist_with_obj_name
            self.thread_safe_general.write_list_to_destination(self.clade_dist_file_path_no_sqrt, self.clade_dist_file_as_list_no_sqrt)

    @staticmethod
    def _append_obj_name_to_dist_line(dict_of_obj_id_to_obj_name, dist_with_obj_name, line):
        temp_list = []
        obj_id = int(line.split('\t')[0])
        obj_name = dict_of_obj_id_to_obj_name[obj_id]
        temp_list.append(obj_name)
        temp_list.extend(line.split('\t'))
        new_line = '\t'.join(temp_list)
        dist_with_obj_name.append(new_line)

    def _chunk_query_at_obj_from_at_uids(self, list_of_obj_uids):
        objs_of_outputs = []
        for uid_list in self.thread_safe_general.chunks(list_of_obj_uids):
            objs_of_outputs.extend(list(AnalysisType.objects.filter(id__in=uid_list)))
        return objs_of_outputs

    def _chunk_query_cc_objs_from_cc_uids(self, list_of_cc_ids):
        cc_of_outputs = []
        for uid_list in self.thread_safe_general.chunks(list_of_cc_ids):
            cc_of_outputs.extend(list(CladeCollection.objects.filter(id__in=uid_list)))
        return cc_of_outputs

    def _write_output_paths_to_stdout(self):
        print('\n\nBrayCurtis distances and PCoA computation complete. Output files:')
        if self.output_path_list:
            print('Output files:')
            for output_path in self.output_path_list:
                print(output_path)
        else:
            print("There are no output files")

    def _append_output_files_to_output_list(self):
        self.output_path_list.extend(
            [self.clade_dist_file_path_no_sqrt, self.clade_pcoa_coord_file_path_no_sqrt,
             self.clade_dist_file_path_sqrt, self.clade_pcoa_coord_file_path_sqrt])


class SampleBrayCurtisDistPCoACreator(BaseBrayCurtisDistPCoACreator):
    """This will produce braycurtis based clade separated distance matrices and compute PCoAs from these
    matrices. It will allow input of data set ids or as sample uid lists that are comma separated.
    It may be that we can have a bas class here which will allow us to have a bray curtis that is specific
    to types and samples.

    Abbreviations
    -------------
    ds = DataSet
    dss = DataSetSample
    dsss = DataSetSampleSequence
    cc = CladeCollection
    """

    def __init__(
            self, js_output_path_dict, html_dir, output_dir, date_time_str=None,
            data_set_sample_uid_list=None,
            data_set_uid_list=None, cct_set_uid_list=None):
        super().__init__(
            date_time_str=date_time_str,
            profiles_or_samples='samples', js_output_path_dict=js_output_path_dict, html_dir=html_dir)

        self.data_set_sample_uid_list, self.clade_col_uid_list = self._set_data_set_sample_uid_list(
            data_set_sample_uid_list=data_set_sample_uid_list, data_set_uid_list=data_set_uid_list,
            cct_set_uid_list=cct_set_uid_list)

        self.cc_list_for_output = self._chunk_query_set_cc_list_from_dss_uids()

        self.clades_of_ccs = list(set([a.clade for a in self.cc_list_for_output]))
        self.output_dir = os.path.join(output_dir, 'between_sample_distances')
        self.thread_safe_general = ThreadSafeGeneral()
        os.makedirs(self.output_dir, exist_ok=True)

    def _chunk_query_set_cc_list_from_dss_uids(self):
        temp_cc_list_for_output = []
        for uid_list in self.thread_safe_general.chunks(self.data_set_sample_uid_list):
            temp_cc_list_for_output.extend(list(CladeCollection.objects.filter(data_set_sample_from__in=uid_list)))
        return temp_cc_list_for_output

    def compute_braycurtis_dists_and_pcoa_coords(self):
        print('\n\nComputing sample pairwise distances and PCoA coordinates using the BrayCurtis method\n')
        for clade_in_question in self.clades_of_ccs:
            dss_obj_to_cct_obj_dict = {cc_obj.data_set_sample_from: cc_obj for cc_obj in
                                       self.cc_list_for_output if cc_obj.clade == clade_in_question}
            self.objs_of_clade = list(dss_obj_to_cct_obj_dict.keys())
            if len(self.objs_of_clade) < 2:
                continue
            self._init_clade_dirs_and_paths(clade_in_question)
            self._create_rs_uid_to_normalised_abund_dict_for_each_obj_samples(dss_obj_to_cct_obj_dict, sqrt=True)
            self._create_rs_uid_to_normalised_abund_dict_for_each_obj_samples(dss_obj_to_cct_obj_dict, sqrt=False)
            self._compute_braycurtis_btwn_obj_pairs(sqrt=True)
            self._compute_braycurtis_btwn_obj_pairs(sqrt=False)
            self._generate_distance_file(sqrt=True)
            self._generate_distance_file(sqrt=False)
            self._add_obj_uids_to_dist_file_and_write(sqrt=True)
            self._add_obj_uids_to_dist_file_and_write(sqrt=False)
            try:
                pcoa_coords_df_sqrt = self._compute_pcoa_coords(clade=clade_in_question, sqrt=True)
                pcoa_coords_df_no_sqrt = self._compute_pcoa_coords(clade=clade_in_question, sqrt=False)
            except EigenValsTooSmallError:
                logging.error(f"The eigenvalues for the clade {clade_in_question} PCoA were too small and were "
                              f"converted to 0s by skbio's implementation of PCoA.")
                logging.error(f"Between sample Bray-Curtis distances cannot be calculated for clade {clade_in_question}")
                continue

            self._populate_js_output_objects(clade_in_question, pcoa_coords_df_sqrt, sqrt=True)
            self._populate_js_output_objects(clade_in_question, pcoa_coords_df_no_sqrt, sqrt=False)
            self._append_output_files_to_output_list()
            self.js_output_path_dict[
                f"btwn_sample_braycurtis_{clade_in_question}_dist_sqrt"] = self.clade_dist_file_path_sqrt
            self.js_output_path_dict[
                f"btwn_sample_braycurtis_{clade_in_question}_pcoa_sqrt"] = self.clade_pcoa_coord_file_path_sqrt
            self.js_output_path_dict[
                f"btwn_sample_braycurtis_{clade_in_question}_dist_no_sqrt"] = self.clade_dist_file_path_no_sqrt
            self.js_output_path_dict[
                f"btwn_sample_braycurtis_{clade_in_question}_pcoa_no_sqrt"] = self.clade_pcoa_coord_file_path_no_sqrt

        self._write_out_js_objects()
        self._write_output_paths_to_stdout()

    def _write_out_js_objects(self):
        self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
            [{'function_name': 'getBtwnSampleDistCoordsBCSqrt', 'python_obj': self.pc_coordinates_dict_sqrt},
             {'function_name': 'getBtwnSampleDistPCVariancesBCSqrt', 'python_obj': self.pc_variances_dict_sqrt},
             {'function_name': 'getBtwnSampleDistPCAvailableBCSqrt', 'python_obj': self.pc_availabaility_dict_sqrt},
             {'function_name': 'getBtwnSampleDistCoordsBCNoSqrt', 'python_obj': self.pc_coordinates_dict_no_sqrt},
             {'function_name': 'getBtwnSampleDistPCVariancesBCNoSqrt', 'python_obj': self.pc_variances_dict_no_sqrt},
             {'function_name': 'getBtwnSampleDistPCAvailableBCNoSqrt', 'python_obj': self.pc_availabaility_dict_no_sqrt}
             ],
            js_outpath=self.js_file_path)

    def _populate_js_output_objects(self, clade_in_question, pcoa_coords_df, sqrt):
        # set the variance dict
        # and set the available pcs
        pcoa_coords_df.set_index('sample_uid', drop=True, inplace=True)
        available_pcs = list(pcoa_coords_df)
        if len(available_pcs) > 6:
            available_pcs = available_pcs[:6]
        if sqrt:
            self.pc_availabaility_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = available_pcs
            variances = [pcoa_coords_df.iloc[-1][pc] for pc in available_pcs]
            self.pc_variances_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = variances
        else:
            self.pc_availabaility_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = available_pcs
            variances = [pcoa_coords_df.iloc[-1][pc] for pc in available_pcs]
            self.pc_variances_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = variances
        # set the coordinates data holder dict here
        genera_pc_coords_dict = {}
        for sample_uid in pcoa_coords_df.index.values.tolist()[:-1]:
            sample_pc_coords_dict = {}
            for pc in available_pcs:
                sample_pc_coords_dict[pc] = f'{pcoa_coords_df.at[sample_uid, pc]:.3f}'
            genera_pc_coords_dict[int(sample_uid)] = sample_pc_coords_dict
        if sqrt:
            self.pc_coordinates_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = genera_pc_coords_dict
        else:
            self.pc_coordinates_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = genera_pc_coords_dict

    def _init_clade_dirs_and_paths(self, clade_in_question):
        self.clade_output_dir = os.path.join(self.output_dir, clade_in_question)

        self.clade_dist_file_path_sqrt = os.path.join(
            self.clade_output_dir, f'{self.date_time_str}_braycurtis_sample_distances_{clade_in_question}_sqrt.dist')

        self.clade_dist_file_path_no_sqrt = os.path.join(
            self.clade_output_dir,
            f'{self.date_time_str}_braycurtis_sample_distances_{clade_in_question}_no_sqrt.dist')
        os.makedirs(self.clade_output_dir, exist_ok=True)

    def _create_rs_uid_to_normalised_abund_dict_for_each_obj_samples(self, dss_obj_to_cct_obj_dict, sqrt):
        # Go through each of the clade collections and create a dict
        # that has key as ref_seq_uid and relative abundance of that sequence
        # we can then store these dict in a dict where the key is the sample ID.
        for dss_obj, clade_col in dss_obj_to_cct_obj_dict.items():
            temp_dict = {}
            list_of_dsss_in_cc = list(DataSetSampleSequence.objects.filter(
                clade_collection_found_in=clade_col))

            if sqrt:
                total_seqs_ind_clade_col = sum([dsss.abundance for dsss in list_of_dsss_in_cc])
                rel_abund_dict = {dsss.id: dsss.abundance/total_seqs_ind_clade_col for dsss in list_of_dsss_in_cc}
                dsss_uid_to_sqrt_rel_abund_dict = {
                    dsss_uid: math.sqrt(rel_abund) for dsss_uid, rel_abund in rel_abund_dict.items()}
                sqr_total = sum(dsss_uid_to_sqrt_rel_abund_dict.values())

                for dsss in list_of_dsss_in_cc:
                    temp_dict[dsss.reference_sequence_of.id] = (dsss_uid_to_sqrt_rel_abund_dict[dsss.id] /
                                                                sqr_total) * 10000
                self.clade_rs_uid_to_normalised_abund_clade_dict_sqrt[dss_obj.id] = temp_dict
            else:
                total_seqs_ind_clade_col = sum([dsss.abundance for dsss in list_of_dsss_in_cc])
                for dsss in list_of_dsss_in_cc:
                    temp_dict[dsss.reference_sequence_of.id] = (dsss.abundance / total_seqs_ind_clade_col)*10000
                self.clade_rs_uid_to_normalised_abund_clade_dict_no_sqrt[dss_obj.id] = temp_dict

    @staticmethod
    def _infer_is_dataset_of_datasetsample(smpl_id_list_str, data_set_string):
        if smpl_id_list_str is not None and data_set_string is not None:
            raise RuntimeError('Please input smpl_id_list_str OR data_set_string')
        elif smpl_id_list_str and data_set_string:
            raise RuntimeError('Please input smpl_id_list_str OR data_set_string')
        else:
            if smpl_id_list_str:
                return True
            else:
                return False


class TypeBrayCurtisDistPCoACreator(BaseBrayCurtisDistPCoACreator):
    """This will produce braycurtis based clade separated distance matrices and compute PCoAs from these
    matrices. It will allow input of data set ids or as sample uid lists that are comma separated.
    It may be that we can have a bas class here which will allow us to have a bray curtis that is specific
    to types and samples.

    Abbreviations
    -------------
    ds = DataSet
    dss = DataSetSample
    dsss = DataSetSampleSequence
    cc = CladeCollection
    """

    def __init__(
            self, data_analysis_obj, js_output_path_dict, html_dir, output_dir,
            date_time_str, data_set_sample_uid_list=None,
            data_set_uid_list=None, cct_set_uid_list=None, local_abunds_only=False):
        super().__init__(
            date_time_str=date_time_str,
            profiles_or_samples='profiles', js_output_path_dict=js_output_path_dict, html_dir=html_dir)

        self.data_set_sample_uid_list, self.clade_col_uid_list = self._set_data_set_sample_uid_list(
            data_set_sample_uid_list=data_set_sample_uid_list, data_set_uid_list=data_set_uid_list,
            cct_set_uid_list=cct_set_uid_list)
        self.data_analysis_obj = data_analysis_obj

        self.cct_set_uid_list = cct_set_uid_list
        if self.cct_set_uid_list is not None:
            # if working with specific profile/sample sets
            self.clade_col_type_objects = self._chunk_query_set_clade_col_type_objs_from_cct_uids()
            self.at_list_for_output = self._chunk_query_set_distinct_at_list_for_output_from_cct_uids(cct_set_uid_list)
        else:
            self.clade_col_type_objects = None
            self.at_list_for_output = self._chunk_query_set_at_list_for_output_from_dss_uids()
        self.clades_of_ats = list(set([at.clade for at in self.at_list_for_output]))
        self.output_dir = os.path.join(output_dir, 'between_profile_distances')
        self.local = local_abunds_only

    def _chunk_query_set_distinct_at_list_for_output_from_cct_uids(self, cct_set_uid_list):
        temp_at_set = set()
        for uid_list in self.thread_safe_general.chunks(cct_set_uid_list):
            temp_at_set.update(list(AnalysisType.objects.filter(cladecollectiontype__id__in=uid_list)))
        return list(temp_at_set)

    def _chunk_query_set_clade_col_type_objs_from_cct_uids(self):
        temp_clade_col_type_objs_list = []
        for uid_list in self.thread_safe_general.chunks(self.cct_set_uid_list):
            temp_clade_col_type_objs_list.extend(list(CladeCollectionType.objects.filter(id__in=uid_list)))
        return temp_clade_col_type_objs_list

    def _chunk_query_set_at_list_for_output_from_dss_uids(self):
        temp_at_set = set()
        for uid_list in self.thread_safe_general.chunks(self.data_set_sample_uid_list):
            temp_at_set.update(list(AnalysisType.objects.filter(
                data_analysis_from=self.data_analysis_obj,
                cladecollectiontype__clade_collection_found_in__data_set_sample_from__in=uid_list)))
        return list(temp_at_set)

    def compute_braycurtis_dists_and_pcoa_coords(self):
        print('\n\nComputing ITS2 type profile pairwise distances and PCoA coordinates using the BrayCurtis method\n')
        for clade_in_question in self.clades_of_ats:
            self.objs_of_clade = [at for at in self.at_list_for_output if at.clade == clade_in_question]
            if len(self.objs_of_clade) < 2:
                continue
            self._init_clade_dirs_and_paths(clade_in_question)
            self._create_rs_uid_to_normalised_abund_dict_for_each_obj_profiles()
            self._compute_braycurtis_btwn_obj_pairs(sqrt=True)
            self._compute_braycurtis_btwn_obj_pairs(sqrt=False)
            self._generate_distance_file(sqrt=True)
            self._generate_distance_file(sqrt=False)
            self._add_obj_uids_to_dist_file_and_write(sqrt=True)
            self._add_obj_uids_to_dist_file_and_write(sqrt=False)

            try:
                pcoa_coords_df_sqrt = self._compute_pcoa_coords(clade=clade_in_question, sqrt=True)
                pcoa_coords_df_no_sqrt = self._compute_pcoa_coords(clade=clade_in_question, sqrt=False)
            except EigenValsTooSmallError:
                logging.error(f"The eigenvalues for the clade {clade_in_question} PCoA were too small and were "
                              f"converted to 0s by skbio's implementation of PCoA.")
                logging.error(f"Between profile Bray-Curtis distances cannot be calculated for clade {clade_in_question}")
                continue


            self._populate_js_output_objects(clade_in_question, pcoa_coords_df_sqrt, sqrt=True)
            self._populate_js_output_objects(clade_in_question, pcoa_coords_df_no_sqrt, sqrt=False)

            self._append_output_files_to_output_list()
            self.js_output_path_dict[
                f"btwn_profile_braycurtis_{clade_in_question}_dist_no_sqrt"] = self.clade_dist_file_path_no_sqrt
            self.js_output_path_dict[
                f"btwn_profile_braycurtis_{clade_in_question}_pcoa_no_sqrt"] = self.clade_pcoa_coord_file_path_no_sqrt
            self.js_output_path_dict[
                f"btwn_profile_braycurtis_{clade_in_question}_dist_sqrt"] = self.clade_dist_file_path_sqrt
            self.js_output_path_dict[
                f"btwn_profile_braycurtis_{clade_in_question}_pcoa_sqrt"] = self.clade_pcoa_coord_file_path_sqrt

        self._write_out_js_objects()
        self._write_output_paths_to_stdout()

    def _write_out_js_objects(self):
        self.thread_safe_general.write_out_js_file_to_return_python_objs_as_js_objs(
            [{'function_name': 'getBtwnProfileDistCoordsBCNoSqrt', 'python_obj': self.pc_coordinates_dict_no_sqrt},
             {'function_name': 'getBtwnProfileDistPCVariancesBCNoSqrt', 'python_obj': self.pc_variances_dict_no_sqrt},
             {'function_name': 'getBtwnProfileDistPCAvailableBCNoSqrt',
              'python_obj': self.pc_availabaility_dict_no_sqrt},
             {'function_name': 'getBtwnProfileDistCoordsBCSqrt', 'python_obj': self.pc_coordinates_dict_sqrt},
             {'function_name': 'getBtwnProfileDistPCVariancesBCSqrt', 'python_obj': self.pc_variances_dict_sqrt},
             {'function_name': 'getBtwnProfileDistPCAvailableBCSqrt', 'python_obj': self.pc_availabaility_dict_sqrt}
             ], js_outpath=self.js_file_path)

    def _populate_js_output_objects(self, clade_in_question, pcoa_coords_df, sqrt):
        # set the variance dict
        # and set the available pcs
        pcoa_coords_df.set_index('analysis_type_uid', drop=True, inplace=True)
        available_pcs = list(pcoa_coords_df)
        if len(available_pcs) > 6:
            available_pcs = available_pcs[:6]
        if sqrt:
            self.pc_availabaility_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = available_pcs
            variances = [pcoa_coords_df.iloc[-1][pc] for pc in available_pcs]
            self.pc_variances_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = variances
        else:
            self.pc_availabaility_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = available_pcs
            variances = [pcoa_coords_df.iloc[-1][pc] for pc in available_pcs]
            self.pc_variances_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = variances

        # set the coordinates data holder dict here
        genera_pc_coords_dict = {}
        for profile_uid in pcoa_coords_df.index.values.tolist()[:-1]:
            profile_pc_coords_dict = {}
            for pc in available_pcs:
                profile_pc_coords_dict[pc] = f'{pcoa_coords_df.at[profile_uid, pc]:.3f}'
            genera_pc_coords_dict[int(profile_uid)] = profile_pc_coords_dict

        if sqrt:
            self.pc_coordinates_dict_sqrt[self.genera_annotation_dict[clade_in_question]] = genera_pc_coords_dict
        else:
            self.pc_coordinates_dict_no_sqrt[self.genera_annotation_dict[clade_in_question]] = genera_pc_coords_dict

    def _init_clade_dirs_and_paths(self, clade_in_question):
        self.clade_output_dir = os.path.join(self.output_dir, clade_in_question)
        # path for the sqrt transformed distance file
        self.clade_dist_file_path_sqrt = os.path.join(
            self.clade_output_dir,
            f'{self.date_time_str}_braycurtis_profile_distances_{clade_in_question}_sqrt.dist')
        # path for the non transformed distance file
        self.clade_dist_file_path_no_sqrt = os.path.join(
            self.clade_output_dir,
            f'{self.date_time_str}_braycurtis_profile_distances_{clade_in_question}_no_sqrt.dist')
        os.makedirs(self.clade_output_dir, exist_ok=True)

    def _create_rs_uid_to_normalised_abund_dict_for_each_obj_profiles(self):
        # Go through each of the clade collections and create a dict
        # that has key as ref_seq_uid and relative abundance of that sequence
        # we can then store these dict in a dict where the key is the sample ID.
        for at in self.objs_of_clade:
            ref_seq_uids_of_analysis_type = [int(b) for b in at.ordered_footprint_list.split(',')]

            df_no_sqrt = pd.DataFrame(at.get_ratio_list())
            df_sqrt = pd.DataFrame(at.get_ratio_list())

            df_sqrt = self.thread_safe_general.sqrt_transform_abundance_df(df_sqrt)

            if self.local:
                normalised_abundance_of_divs_dict_no_sqrt = self._create_norm_abund_dict_from_local_clade_cols(
                    at, df_no_sqrt, ref_seq_uids_of_analysis_type)
                normalised_abundance_of_divs_dict_sqrt = self._create_norm_abund_dict_from_local_clade_cols(
                    at, df_sqrt, ref_seq_uids_of_analysis_type)

            elif self.cct_set_uid_list:
                normalised_abundance_of_divs_dict_no_sqrt = self._create_norm_abund_dict_from_cct_set(
                    at, df_no_sqrt, ref_seq_uids_of_analysis_type)
                normalised_abundance_of_divs_dict_sqrt = self._create_norm_abund_dict_from_cct_set(
                    at, df_sqrt, ref_seq_uids_of_analysis_type)
            else:
                normalised_abundance_of_divs_dict_no_sqrt = self._create_norm_abund_dict_from_all_clade_cols(
                    df_no_sqrt, ref_seq_uids_of_analysis_type)
                normalised_abundance_of_divs_dict_sqrt = self._create_norm_abund_dict_from_all_clade_cols(
                    df_sqrt, ref_seq_uids_of_analysis_type)

            self.clade_rs_uid_to_normalised_abund_clade_dict_no_sqrt[at.id] = normalised_abundance_of_divs_dict_no_sqrt
            self.clade_rs_uid_to_normalised_abund_clade_dict_sqrt[at.id] = normalised_abundance_of_divs_dict_sqrt

    def _create_norm_abund_dict_from_local_clade_cols(self, at, df, ref_seq_uids_of_analysis_type):
        # we want to limit the inference of the average relative abundance of the DIVs to only those
        # div abundance values that were found in samples from the output in question
        # as such we need to first get a list of the clade collections and see which of these
        # clade collections were from self.data_set_sample_uid_list
        clade_collection_uid_list_of_type = [
            int(a) for a in at.list_of_clade_collections.split(',')]
        # the indices of the clade collections that are of this output
        indices = []
        for i in range(len(clade_collection_uid_list_of_type)):
            if clade_collection_uid_list_of_type[i] in self.clade_col_uid_list:
                indices.append(i)
        normalised_abundance_of_divs_dict = {
            ref_seq_uids_of_analysis_type[i]: math.ceil(df.iloc[indices, i].mean() * 10000) for
            i in range(len(ref_seq_uids_of_analysis_type))}
        return normalised_abundance_of_divs_dict

    def _create_norm_abund_dict_from_cct_set(
            self, at, df, ref_seq_uids_of_analysis_type):
        # Then we are working with a custom set of CladeCollection-AnalysisType associations. i.e. for every
        # analysis type we are working with we will have to check exactly which CladeCollections we should
        # be infering the average DIV abundances from.
        # Because the list of AnalysisTypes that we are looking through are defined by the CladeCollectionTypes
        # that we're working with, we know that there are no types in which we aren't looking at at least
        # one CladeCollection from.
        clade_collection_uid_list_of_type = [
            int(a) for a in at.list_of_clade_collections.split(',')]
        clade_collection_uids_to_output_for_this_type = []
        for cct_obj in self.clade_col_type_objects:
            if cct_obj.analysis_type_of.id == at.id:  # then this is a CladeCollectionType of the at
                clade_collection_uids_to_output_for_this_type.append(cct_obj.clade_collection_found_in.id)
        indices = []
        for i in range(len(clade_collection_uid_list_of_type)):
            if clade_collection_uid_list_of_type[i] in clade_collection_uids_to_output_for_this_type:
                indices.append(i)
        normalised_abundance_of_divs_dict = {
            ref_seq_uids_of_analysis_type[i]: math.ceil(df.iloc[indices, i].mean() * 10000) for
            i in range(len(ref_seq_uids_of_analysis_type))}
        return normalised_abundance_of_divs_dict

    @staticmethod
    def _create_norm_abund_dict_from_all_clade_cols(df, ref_seq_uids_of_analysis_type):
        normalised_abundance_of_divs_dict = {
            ref_seq_uids_of_analysis_type[i]: math.ceil(df[i].mean() * 100000) for
            i in range(len(ref_seq_uids_of_analysis_type))}
        return normalised_abundance_of_divs_dict