import pandas as pd
from multiprocessing import Queue, Manager, Process, current_process
import sys
from dbApp.models import DataSetSampleSequence
import pickle
import os
class VirtualObjectManager():
    """This class will link together an instance of a VirtualCladeCollectionManger and a VirtualAnalaysisTypeManger.
    I will therefore allow VirtualAnalysisTypes to access the information in the VirtualCladeCollections."""
    def __init__(self, parent_sp_data_analysis):
        self.sp_data_analysis = parent_sp_data_analysis
        self.cc_manager = VirtualCladeCollectionManager(obj_manager=self)
        # with open(os.path.join(self.sp_data_analysis.workflow_manager.symportal_root_directory, 'tests', 'objects', 'cc_manager.p'), 'rb') as f:
        #     self.cc_manager = pickle.load(f)
        self.vat_manager = VirtualAnalysisTypeManager(obj_manager=self)


class VirtualCladeCollectionManager():
    """Unlike the VirtualAnalysisType the VirtualCladeCollection will be a proxy for an object that already exists
    in the datbase already. As such we won't need to generate pks."""
    def __init__(self, obj_manager):
        self.obj_manager = obj_manager
        self.clade_collection_instances_dict = {}

        self._populate_virtual_cc_manager_from_db()


    def _populate_virtual_cc_manager_from_db(self):
        """When first instantiated we should grab all of the CladeCollections from the database that are part of
        this DataAnalysis and make VirtualCladeCollectionsFrom them.
        We will need to populate the VirtualAnalysisTypeManager vat_dict before
        we can populate the analysis_type_obj_to_representative_rel_abund_in_cc_dicts for each of the
        VirtualCladeCollections so we will do this in a seperate method.
        """

        self.clade_collection_instances_dict = self._create_cc_info_dict()

    def _create_cc_info_dict(self):
        print('Instantiating VirtualCladeCollectionManager')
        cc_input_mp_queue = Queue()
        mp_manager = Manager()
        cc_to_info_items_mp_dict = mp_manager.dict()

        for cc in self.obj_manager.sp_data_analysis.ccs_of_analysis:
            cc_input_mp_queue.put(cc)

        for n in range(self.obj_manager.sp_data_analysis.workflow_manager.args.num_proc):
            cc_input_mp_queue.put('STOP')

        all_processes = []
        for n in range(self.obj_manager.sp_data_analysis.workflow_manager.args.num_proc):
            p = Process(target=self._vcc_id_to_vcc_obj_worker, args=(cc_input_mp_queue, cc_to_info_items_mp_dict))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

        return dict(cc_to_info_items_mp_dict)

    def _vcc_id_to_vcc_obj_worker(self, cc_input_mp_queue, cc_to_info_items_mp_dict):
        for clade_collection_object in iter(cc_input_mp_queue.get, 'STOP'):
            sys.stdout.write(f'\r{clade_collection_object.data_set_sample_from.name} {current_process().name}')

            dss_objects_of_cc_list = list(DataSetSampleSequence.objects.filter(
                clade_collection_found_in=clade_collection_object))

            sorted_dss_objects_of_cc_list = [dsss for dsss in sorted(dss_objects_of_cc_list, key=lambda x: x.abundance, reverse=True)]

            list_of_ref_seq_uids_in_cc = [
                dsss.reference_sequence_of.id for dsss in dss_objects_of_cc_list]

            above_cutoff_ref_seqs_obj_set = clade_collection_object.cutoff_footprint(
                self.obj_manager.sp_data_analysis.data_analysis_obj.within_clade_cutoff)

            total_sequences_in_cladecollection = sum([dsss.abundance for dsss in dss_objects_of_cc_list])

            list_of_rel_abundances = [dsss.abundance / total_sequences_in_cladecollection for dsss in
                                  dss_objects_of_cc_list]

            ref_seq_frozen_set = frozenset(dsss.reference_sequence_of.id for dsss in dss_objects_of_cc_list)

            ref_seq_id_to_rel_abund_dict = {}
            for i in range(len(dss_objects_of_cc_list)):
                ref_seq_id_to_rel_abund_dict[list_of_ref_seq_uids_in_cc[i]] = list_of_rel_abundances[i]

            cc_to_info_items_mp_dict[clade_collection_object.id] = VirtualCladeCollection(
                clade=clade_collection_object.clade,
                footprint_as_frozen_set_of_ref_seq_uids=ref_seq_frozen_set,
                ref_seq_id_to_rel_abund_dict=ref_seq_id_to_rel_abund_dict,
                total_seq_abundance=total_sequences_in_cladecollection,
                cc_object=clade_collection_object,
                above_cutoff_ref_seqs_obj_set=above_cutoff_ref_seqs_obj_set,
                ordered_dsss_objs=sorted_dss_objects_of_cc_list,
                sample_from_name=str(clade_collection_object))


class VirtualCladeCollection:
    """A RAM stored representation of a CladeCollection object that already exists in the DB"""
    def __init__(
            self, clade, footprint_as_frozen_set_of_ref_seq_uids, ref_seq_id_to_rel_abund_dict,
            total_seq_abundance, cc_object, above_cutoff_ref_seqs_obj_set, ordered_dsss_objs, sample_from_name=None):

        self.clade = clade
        self.cc_object = cc_object
        self.id = self.cc_object.id
        # This is the ref seq uids for all dss found in the cc as oposed to just those above the
        # within_clade_cutoff. The above cutoff equivalents are stored below
        self.footprint_as_frozen_set_of_ref_seq_uids = footprint_as_frozen_set_of_ref_seq_uids
        self.ref_seq_id_to_rel_abund_dict = ref_seq_id_to_rel_abund_dict
        self.total_seq_abundance = total_seq_abundance
        self.sample_from_name = sample_from_name
        self.above_cutoff_ref_seqs_obj_set = above_cutoff_ref_seqs_obj_set
        self.above_cutoff_ref_seqs_id_set = [rs.id for rs in self.above_cutoff_ref_seqs_obj_set]
        self.ordered_dsss_objs = ordered_dsss_objs

        # key = AnalysisType object, value = the relative abundance of the cc that this AnalysisType represents
        self.analysis_type_obj_to_representative_rel_abund_in_cc_dict = {}

    def __str__(self):
        try:
            return self.sample_from_name
        except Exception:
            return f'VirtualCladeCollection uid: {self.id}'



class VirutalAnalysisTypeInit:
    def __init__(self, parent_vat_manager, vat_to_init):
        self.vat = vat_to_init
        self.vat_manager = parent_vat_manager

    def init_vat_post_profile_assignment(self):
        self._make_rel_abund_dfs_post_prof_assignment()

        self._generate_maj_ref_seq_set_and_infer_codom(self.vat.post_prof_assignment_rel_abund_df)

        self._generate_name(self.vat.post_prof_assignment_rel_abund_df)

    def _make_rel_abund_dfs_post_prof_assignment(self):
        at_df = pd.DataFrame(index=[cc.id for cc in self.vat.clade_collection_obj_set_profile_assignment],
                             columns=[rs.id for rs in self.vat.footprint_as_ref_seq_objs_set])
        for cc in self.vat.clade_collection_obj_set_profile_assignment:
            ref_seq_abund_dict_for_cc = self.vat_manager.obj_manager.cc_manager.clade_collection_instances_dict[
                cc.id].ref_seq_id_to_rel_abund_dict
            at_df.loc[cc.id] = pd.Series(
                {rs_uid_key: rs_rel_abund_val for rs_uid_key, rs_rel_abund_val in ref_seq_abund_dict_for_cc.items()
                 if
                 rs_uid_key in list(at_df)})
        at_df["sum"] = at_df.sum(axis=1)
        at_df = at_df.iloc[:, 0:-1].div(at_df["sum"], axis=0)
        self.vat.post_prof_assignment_rel_abund_df = at_df.reindex(
            at_df.sum().sort_values(ascending=False).index, axis=1).astype('float')

    def init_vat_pre_profile_assignment(self):
        self._make_rel_abund_dfs()

        self._generate_maj_ref_seq_set_and_infer_codom(self.vat.relative_seq_abund_profile_assignment_df)

        self._populate_artefact_set()

        self._populate_max_min_profile_assignment_requirement_dict()

        self.vat.non_artefact_ref_seq_uid_set = set([
            rs_id for rs_id in self.vat.ref_seq_uids_set if rs_id not in self.vat.artefact_ref_seq_uid_set])

        self._set_basal_seq()

        self._generate_name(self.vat.relative_seq_abund_profile_assignment_df)

        return self.vat

    def _populate_artefact_set(self):
        # Identify those DIVs of the analysis type that are artefact seqs for the purporses of checking artefacts
        # NB when we are checking for AnalysisTypes that are caused due to artefact DIVs we are working with rel
        # abundances as a proportion of all of the sequences in a CladeCollection. This is in contrast
        # to when we are working with the required abundances during TypeAssignment when we are working with rel
        # abundances of DIVs as a proportion of the CladeCollection sequences that DIVs in the type in question.
        for rs_col in self.vat.relative_seq_abund_profile_discovery_df:
            min = self.vat.relative_seq_abund_profile_discovery_df[rs_col].min()
            if min < 0.06:
                self.vat.artefact_ref_seq_uid_set.add(rs_col)

    def _populate_max_min_profile_assignment_requirement_dict(self):
        # populate the max_min dict that will be used during profile assignment
        for rs_col in self.vat.relative_seq_abund_profile_assignment_df:
            max = self.vat.relative_seq_abund_profile_assignment_df[rs_col].max()
            min = self.vat.relative_seq_abund_profile_assignment_df[rs_col].min()
            if min < 0.06:
                min = 0.0001
            self.vat.prof_assignment_required_rel_abund_dict[rs_col] = self.RefSeqReqAbund(
                max_rel_abund=max, min_rel_abund=min)

    def _generate_maj_ref_seq_set_and_infer_codom(self, vat_df):
        # get the most abund rs for each cc
        majority_reference_sequence_uid_set = set()
        for index, row in vat_df.iterrows():
            majority_reference_sequence_uid_set.add(row.idxmax())
        self.vat.majority_reference_sequence_uid_set = majority_reference_sequence_uid_set
        if len(self.vat.majority_reference_sequence_uid_set) > 1:
            self.vat.co_dominant = True
        else:
            self.vat.co_dominant = False
        self.majority_reference_sequence_obj_set = set(
            [rs for rs in self.vat.footprint_as_ref_seq_objs_set if
             rs.id in self.vat.majority_reference_sequence_uid_set])

    def _make_rel_abund_dfs(self):
        at_df = self._create_rel_seq_abund_profile_disco_df()
        prof_ass_df = self._create_rel_seq_abund_prof_assign_df(at_df)
        self._reorder_dfs(at_df, prof_ass_df)

    def _reorder_dfs(self, at_df, prof_ass_df):
        # We will sort both DataFrames according to the summed abundances in the relative_seq_abund_prof_assign.
        # https://stackoverflow.com/questions/26537878/pandas-sum-across-columns-and-divide-each-cell-from-that-value
        self.vat.relative_seq_abund_profile_discovery_df = at_df.reindex(
            prof_ass_df.sum().sort_values(ascending=False).index, axis=1).astype('float')
        self.vat.relative_seq_abund_profile_assignment_df = prof_ass_df.reindex(
            prof_ass_df.sum().sort_values(ascending=False).index, axis=1).astype('float')

    def _create_rel_seq_abund_prof_assign_df(self, at_df):
        # compute the relative_seq_abund_profile_assignment_df from the relative_seq_abund_profile_discovery_df
        prof_ass_df = at_df.copy()
        prof_ass_df["sum"] = prof_ass_df.sum(axis=1)
        prof_ass_df = prof_ass_df.iloc[:, 0:-1].div(prof_ass_df["sum"], axis=0)
        return prof_ass_df

    def _create_rel_seq_abund_profile_disco_df(self):
        # create and populate the relative_seq_abund_profile_discovery_df
        at_df = pd.DataFrame(index=[cc.id for cc in self.vat.clade_collection_obj_set_profile_discovery],
                             columns=[rs.id for rs in self.vat.footprint_as_ref_seq_objs_set])
        for cc in self.vat.clade_collection_obj_set_profile_discovery:
            ref_seq_abund_dict_for_cc = self.vat_manager.obj_manager.cc_manager.clade_collection_instances_dict[
                cc.id].ref_seq_id_to_rel_abund_dict
            at_df.loc[cc.id] = pd.Series(
                {rs_uid_key: rs_rel_abund_val for rs_uid_key, rs_rel_abund_val in ref_seq_abund_dict_for_cc.items()
                 if
                 rs_uid_key in list(at_df)})
        return at_df

    class RefSeqReqAbund:
        """A very simple object that holds the maximum and mimum relative abundances for a DIV of an AnalysisType """
        def __init__(self, max_rel_abund, min_rel_abund):
            # The maximum allowable relative abundance of the RefSeq in the CC in question
            self.max_abund = max_rel_abund
            # The minimum allowable relative abundance of the RefSeq in the CC in question
            self.min_abund = min_rel_abund


    def _set_basal_seq(self):
        basal_set = set()
        found_c15_a = False
        for rs in self.vat.footprint_as_ref_seq_objs_set:
            if rs.name == 'C3':
                basal_set.add('C3')
            elif rs.name == 'C1':
                basal_set.add('C1')
            elif 'C15' in rs.name and not found_c15_a:
                basal_set.add('C15')
                found_c15_a = True

        if len(basal_set) == 1:
            self.vat.basal_seq = list(basal_set)[0]
        elif len(basal_set) > 1:
            raise RuntimeError(f'basal seq set {basal_set} contains more than one ref seq')
        else:
            self.vat.basal_seq = None

    def _generate_name(self, at_df):
        if self.vat.co_dominant:
            list_of_maj_ref_seq = [rs for rs in self.vat.footprint_as_ref_seq_objs_set if rs.id in self.vat.majority_reference_sequence_uid_set]
            # Start the name with the co_dominant intras in order of abundance.
            # Then append the nonco_dominant intras in order of abundance
            ordered_list_of_co_dom_ref_seq_obj = []
            for ref_seq_id in list(at_df):
                for ref_seq in list_of_maj_ref_seq:
                    if ref_seq.id == ref_seq_id:
                        ordered_list_of_co_dom_ref_seq_obj.append(ref_seq)

            co_dom_name_part = '/'.join(rs.name for rs in ordered_list_of_co_dom_ref_seq_obj)

            list_of_remaining_ref_seq_objs = []
            for ref_seq_id in list(at_df):
                for ref_seq in self.vat.footprint_as_ref_seq_objs_set:
                    if ref_seq not in ordered_list_of_co_dom_ref_seq_obj and ref_seq.id == ref_seq_id:
                        list_of_remaining_ref_seq_objs.append(ref_seq)

            if list_of_remaining_ref_seq_objs:
                co_dom_name_part += '-{}'.format('-'.join([rs.name for rs in list_of_remaining_ref_seq_objs]))
            self.vat.name = co_dom_name_part
        else:
            ordered_list_of_ref_seqs = []
            for ref_seq_id in list(at_df):
                for ref_seq in self.vat.footprint_as_ref_seq_objs_set:
                    if ref_seq.id == ref_seq_id:
                        ordered_list_of_ref_seqs.append(ref_seq)
            self.vat.name = '-'.join(rs.name for rs in ordered_list_of_ref_seqs)


class VirtualAnalysisTypeManager():
    """This is a class that will manage the collection of VirtualAnalysisType instances that exist in memory.
    It will be used to generate a new type (so that it can assign a uid) and it will be used to delete too."""
    def __init__(self, obj_manager):
        self.obj_manager = obj_manager
        self.next_uid = 1
        # key = uid of at, value = VirtualAnalysisType instance
        self.vat_dict = {}

    def make_vat_post_profile_assignment(self, clade_collection_obj_list, ref_seq_obj_list):
        new_vat = self.VirtualAnalysisType(
            clade_collection_obj_list_post_prof_assignment=clade_collection_obj_list,
            ref_seq_obj_list=ref_seq_obj_list, id=self.next_uid)
        vat_init = VirutalAnalysisTypeInit(parent_vat_manager=self, vat_to_init=new_vat)
        vat_init.init_vat_post_profile_assignment()

        self.vat_dict[new_vat.id] = new_vat

        self.next_uid += 1

        return new_vat

    def make_vat_pre_profile_assignment(self, clade_collection_obj_list, ref_seq_obj_list):

        new_vat = self.VirtualAnalysisType(
            clade_collection_obj_list_pre_prof_assignment=clade_collection_obj_list,
            ref_seq_obj_list=ref_seq_obj_list,
            id=self.next_uid)
        vat_init = VirutalAnalysisTypeInit(parent_vat_manager=self, vat_to_init=new_vat)
        vat_init.init_vat_pre_profile_assignment()

        self.vat_dict[new_vat.id] = new_vat

        self.next_uid += 1

        return new_vat

    def reinit_vat_post_profile_assignment(self, vat_to_reinit, new_clade_collection_obj_set):

        vat_to_reinit.clade_collection_obj_set_profile_assignment = set(new_clade_collection_obj_set)

        vat_init = VirutalAnalysisTypeInit(parent_vat_manager=self, vat_to_init=vat_to_reinit)

        vat_init.init_vat_post_profile_assignment()

    def reinit_vat_pre_profile_assignment(self, vat_to_reinit, new_clade_collection_obj_set):
        vat_to_reinit.clade_collection_obj_set_profile_discovery = set(new_clade_collection_obj_set)

        vat_init = VirutalAnalysisTypeInit(parent_vat_manager=self, vat_to_init=vat_to_reinit)

        vat_init.init_vat_pre_profile_assignment()



    def delete_virtual_analysis_type(self, virtual_analysis_type):
        try:
            del self.vat_dict[virtual_analysis_type.id]
        except KeyError:
            raise RuntimeError(
                f'VirtualAnalysisType {virtual_analysis_type} '
                f'not found in the VirtualAnalysisTypeManager\'s collection')

    def add_ccs_and_reinit_virtual_analysis_type(self, vat_to_add_ccs_to, list_of_clade_collection_objs_to_add):

        new_clade_collection_obj_set_profile_discovery = \
            vat_to_add_ccs_to.clade_collection_obj_set_profile_discovery.union(
            set(list_of_clade_collection_objs_to_add))

        self.reinit_vat_pre_profile_assignment(
            vat_to_reinit=vat_to_add_ccs_to,
            new_clade_collection_obj_set=new_clade_collection_obj_set_profile_discovery)

    def remove_cc_and_reinit_vat_pre_profile_assignment(
            self, vat_to_remove_ccs_from, list_of_clade_collection_objs_to_remove):

        new_clade_collection_obj_set_profile_discovery = \
            vat_to_remove_ccs_from.clade_collection_obj_set_profile_discovery - set(
                list_of_clade_collection_objs_to_remove)

        self.reinit_vat_pre_profile_assignment(
            vat_to_reinit=vat_to_remove_ccs_from,
            new_clade_collection_obj_set=new_clade_collection_obj_set_profile_discovery)



    class VirtualAnalysisType():
        """A RAM stored representation of the AnalysisType object. Instances of these objects do not yet
        exist in the database. We will eventually use these instances to make make AnalysisType objects that can be
        stored in the db. I am hoping that by using these virtual objects that we will be able to cut down on some of the
        attributes held in the AnalysisType model fields as many of these are used in the actual ananlysis. The only
        attributes we would need to keep hold of those are those that are used in the outputs.

        When doing init pre-profile assignment we will init from the clade_collection_obj_set_profile_discovery and
        populate the relative_seq_abund_profile_discovery_df and the relative_seq_abund_profile_assignment_df.

        When doing init post-profile assignment we will init from the clade_collection_obj_set_profile_assignment and
        populate the post_prof_assignment_rel_abund_df
        """

        def __init__(
                self, ref_seq_obj_list, id, clade_collection_obj_list_pre_prof_assignment=None,
                clade_collection_obj_list_post_prof_assignment=None):
            self.id = id
            # There will be two different clade_collection_obj_sets. Firstly there is the set of CCs that are associated
            # to this VirtualAnalysisType during ProfileDiscovery. These CCs are used to define the max and min abundances
            # that ref seqs need to be found at.
            # Secondly there will be the list of CladeCollections in which this VirtualAnalysisType is found during
            # ProfileAssignment.
            if clade_collection_obj_list_pre_prof_assignment is not None:
                self.clade_collection_obj_set_profile_discovery = set(clade_collection_obj_list_pre_prof_assignment)
                self.clade_collection_obj_set_profile_assignment = set()
                self.clade = list(clade_collection_obj_list_pre_prof_assignment)[0].clade
            else:
                self.clade_collection_obj_set_profile_discovery = set()
                self.clade_collection_obj_set_profile_assignment = set(clade_collection_obj_list_post_prof_assignment)
                self.clade = list(clade_collection_obj_list_post_prof_assignment)[0].clade

            self.footprint_as_ref_seq_objs_set = ref_seq_obj_list
            self.ref_seq_uids_set = set([rs.id for rs in self.footprint_as_ref_seq_objs_set])
            # NB in the type discovery part the DataAnalysis we will be concerned with relative sequence abundances
            # as a proportion of all of the sequences found within a CladeCollection. But, as we move into ProfileAssignment
            # we will be concerned with the relative abundances of the sequences as a proportion of only those sequences
            # in the CladeCollection that are found within the AnalysisType in Question.
            # E.g. in a CladeCollection that contains C3-0.4, C3b-0.1, C15-0.4, C15b-0.1, when working with an
            # AnalysisType of footprint C3-C3b we will use the relabundances of 0.4 and 0.1, respectively. When
            # working in ProfileAssignment we will use C3-0.8, C3b-0.2.
            # To work out the relative abundances for ProfileAssignment we can simply divide the rel abundances
            # for ProfileDiscovery by their summed rel abundances. We will hold two seperate DataFrame objects representing
            # each of these differnt relative abundances.
            self.relative_seq_abund_profile_discovery_df = None
            self.relative_seq_abund_profile_assignment_df = None

            # To do the multimodal detection we need to create a final relative_seq_abund_df. This will be based on the
            # CladeCollections in which the VirtualAnalysisType was found (clade_collection_obj_set_profile_assignment).
            # It will also be based on the CladeCollection
            # total seqs that are DIVs of the VirtualAnalysisType rather than all of the seqs in the CladeCollection.
            self.post_prof_assignment_rel_abund_df = None

            self.artefact_ref_seq_uid_set = set()
            self.non_artefact_ref_seq_uid_set = set()
            self.co_dominant = None
            self.majority_reference_sequence_uid_set = set()
            self.majority_reference_sequence_obj_set = set()
            self.name = None

            # key = ref seq id, val=RefSeqReqAbund object
            self.prof_assignment_required_rel_abund_dict = {}

            # will be used to hold the species information associated at the end of the analysis
            self.species = None

        def generate_name(self, at_df):
            if self.co_dominant:
                list_of_maj_ref_seq = [rs for rs in self.footprint_as_ref_seq_objs_set if
                                       rs.id in self.majority_reference_sequence_uid_set]
                # Start the name with the co_dominant intras in order of abundance.
                # Then append the nonco_dominant intras in order of abundance
                ordered_list_of_co_dom_ref_seq_obj = []
                for ref_seq_id in list(at_df):
                    for ref_seq in list_of_maj_ref_seq:
                        if ref_seq.id == ref_seq_id:
                            ordered_list_of_co_dom_ref_seq_obj.append(ref_seq)

                co_dom_name_part = '/'.join(rs.name for rs in ordered_list_of_co_dom_ref_seq_obj)

                list_of_remaining_ref_seq_objs = []
                for ref_seq_id in list(at_df):
                    for ref_seq in self.footprint_as_ref_seq_objs_set:
                        if ref_seq not in ordered_list_of_co_dom_ref_seq_obj and ref_seq.id == ref_seq_id:
                            list_of_remaining_ref_seq_objs.append(ref_seq)

                if list_of_remaining_ref_seq_objs:
                    co_dom_name_part += '-{}'.format('-'.join([rs.name for rs in list_of_remaining_ref_seq_objs]))
                self.name = co_dom_name_part
            else:
                ordered_list_of_ref_seqs = []
                for ref_seq_id in list(at_df):
                    for ref_seq in self.vat.footprint_as_ref_seq_objs_set:
                        if ref_seq.id == ref_seq_id:
                            ordered_list_of_ref_seqs.append(ref_seq)
                self.name = '-'.join(rs.name for rs in ordered_list_of_ref_seqs)


        def __str__(self):
            return self.name

