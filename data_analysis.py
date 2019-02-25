import os
import shutil
from multiprocessing import Queue, Manager, Process, current_process
import sys
from dbApp.models import AnalysisType, DataSetSampleSequence
import itertools
from collections import defaultdict

from django import db
class SPDataAnalysis:
    def __init__(self, workflow_manager_parent, data_analysis_obj):
        self.parent = workflow_manager_parent
        self.temp_wkd = os.path.join(self.parent.symportal_root_directory)
        self._setup_temp_wkd()
        self.data_analysis_obj = data_analysis_obj
        # The abundance that a given DIV must be found at when it has been considered 'unlocked'
        # https://github.com/didillysquat/SymPortal_framework/wiki/The-SymPortal-logic#type-profile-assignment---logic
        self.unlocked_abundance = 0.0001
        self.clade_list = list('ABCDEFGH')
        self.ccs_of_analysis = self.data_analysis_obj.get_clade_collections()
        # List that will hold a dictionary for each clade
        # Each dictionary will hold key = footprint (set of sequences)
        # value = [[] []] where []0 = list of cladeCollections containing given footprint
        # and []1 = list of majority sequence for given sample
        self.clade_footp_dicts_list = [{} for _ in self.clade_list]
        self.collapsed_footprint_dict = None
        self.current_clade = None

    def analyse_data(self):
        print('Beginning profile discovery')
        self._populate_clade_fp_dicts_list()

        self._make_analysis_types()

        self._check_for_additional_artefact_types()

    def _check_for_additional_artefact_types(self):
        # Generate a dict that simply holds the total number of seqs per clade_collection_object
        # This will be used when working out relative proportions of seqs in the clade_collection_object

        cc_to_total_seqs_dict = create_total_seqs_dict_for_all_ccs(num_processors)

        # Generate dict per cc for listof reference sequences and their abundances in the clade_collection_object
        cc_to_ref_seq_list_and_abundances = create_ref_seq_rel_abunds_for_all_ccs_dict(
            cc_to_total_seqs_dict, num_processors)

        # Get list of all types
        all_types_from_data_analysis = AnalysisType.objects.filter(data_analysis_from=analysis_object)

        # Get list of clades that are represented by the types
        clade_list = set()
        for at in all_types_from_data_analysis:
            clade_list.add(at.clade)

        # Create dict for each type of refSeqs, artefact and non-artefact DIVs, and footprint
        type_footprint_dict = create_type_artefact_div_info_dict(all_types_from_data_analysis)

        # 08/12/17 this is going to need to be looked at as each clade_collection_object can now have multiple initial types
        # Create a dict that is clade_collection_object: inital type found in
        cc_to_initial_type_dict = create_clade_collection_to_initial_type_dictionary(all_types_from_data_analysis)

        # Do pairwise comparison of types within clades
        # 08/12/17
        # I have gone through all of the remaining code in this function and fixed the cc to initial type dict
        # so that it is compatible with ccts having multiple types.
        # we still need to do the other artefact checking
        for clade in clade_list:
            all_types_from_data_analysis = AnalysisType.objects.filter(data_analysis_from=analysis_object, clade=clade)

            static_list_of_types = list(AnalysisType.objects.filter(data_analysis_from=analysis_object, clade=clade))
            done_list = []
            while 1:
                restart = False

                # Pairwise comparison of types within clade

                # For implementing multiprocessing in this part of the code
                # We may be able to create a list of all of the pairwise comparisons and then add this to a queue
                # In essence this would mean that all of the code below the a,b iter would go into the worker.
                # However, it might be that the time limiting part of all of this is going to be
                # the check_type_pairing_for_artefact_type
                # So, perhaps we should look to see whether this component can be multiprocessed.
                # BUT because we may be deleting in and adding types this might be a problem
                # so it might be best if we can try to speed up the internals of the checkTypePa...
                for a, b, in itertools.combinations(
                        [at for at in all_types_from_data_analysis if at in static_list_of_types], 2):
                    # Check to see if the pairs have already been done as we will be running through this process several
                    # times in all likelihood

                    if (a.id, b.id) not in done_list:
                        list_of_non_artefact_intras_a = type_footprint_dict[a.id][0]
                        list_of_non_artefact_intras_b = type_footprint_dict[b.id][0]

                        # 08/12/17 here we can simply check to see if the pair of types contain incompatible types
                        # if they do then we don't compare them and simply add to the done list
                        # returns False if it is OK to compare, i.e. non incompatabiites
                        # returns True if there are incompatabilities and they cannot be compared
                        if not check_if_type_pairs_contain_incompatible_basal_seqs(a, b):

                            # If all of the nonArtefact intras of A are found in the footprint of B
                            # and if all of the nonArtefact intras of B are found in footprint of A.
                            # i.e. if only artefact intras differentiate the types
                            if set(list_of_non_artefact_intras_a).issubset(type_footprint_dict[b.id][1]) and set(
                                    list_of_non_artefact_intras_b).issubset(type_footprint_dict[a.id][1]):
                                # Check to see whether either a or b are subsets of each other
                                a_subset_b = type_footprint_dict[a.id][1].issubset(type_footprint_dict[b.id][1])
                                b_subset_a = type_footprint_dict[b.id][1].issubset(type_footprint_dict[a.id][1])

                                if not a_subset_b and not b_subset_a:
                                    list_of_artefact_intras_in_both_types = list(a.artefact_intras.split(','))
                                    list_of_artefact_intras_in_both_types.extend(b.artefact_intras.split(','))
                                    if list_of_artefact_intras_in_both_types:
                                        # Here we have finally found a type pairing that needs to be checked
                                        print('\nChecking {} and {} for additional artefactual profiles'.format(a, b))
                                        if check_type_pairing_for_artefact_type(
                                                a, b, type_footprint_dict, clade, cc_to_initial_type_dict,
                                                cc_to_ref_seq_list_and_abundances, num_processors):
                                            # If we did find a new type then we need to start the a,b comparison again
                                            # It should be quick though due to the done_list
                                            restart = True
                                            # We need to recalculate this as we may have delted and created types
                                            # DONE this may be stuck infinitely looping by searching the new types created.
                                            # Maybe only work with the types we started with.
                                            # We could make a list of original types and then just look through the types
                                            # that are found in common between the original types list and an
                                            # updated types list. This way no newly added types will be analysed.
                                            # and we should avoid a run away situation
                                            all_types_from_data_analysis = AnalysisType.objects.filter(
                                                data_analysis_from=analysis_object, clade=clade)
                                            break
                                        else:
                                            done_list.extend([(a.id, b.id), (b.id, a.id)])

                                    else:
                                        done_list.extend([(a.id, b.id), (b.id, a.id)])
                                else:
                                    done_list.extend([(a.id, b.id), (b.id, a.id)])
                            else:
                                # Whenever we fall out of the series of conditionals
                                # we should add both possible combinations
                                # of type a and type b to the done list so that we don't waste time going through them again
                                done_list.extend([(a.id, b.id), (b.id, a.id)])
                        else:
                            # then they have incompatable basal seqs and must not be compared
                            done_list.extend([(a.id, b.id), (b.id, a.id)])
                # If we make it to here we did a full run through the types without making a new type
                # Time to do the same for the next type
                if not restart:
                    break


        return cc_to_total_seqs_dict, cc_to_ref_seq_list_and_abundances, type_footprint_dict, cc_to_initial_type_dict

    def _make_analysis_types(self):
        for clade_fp_dict in self.clade_footp_dicts_list:
            self.current_clade = self.parent.clade_list[self.clade_footp_dicts_list.index(clade_fp_dict)]
            if self._there_are_footprints_of_this_clade(
                    clade_fp_dict):  # If there are some clade collections for the given clade

                sfi = SupportedFootPrintIdentifier(clade_footprint_dict=clade_fp_dict, parent_sp_data_analysis=self)
                self.collapsed_footprint_dict = sfi.identify_supported_footprints()

                analysis_type_creator = AnalysisTypeCreator(parent_sp_data_analysis=self)
                analysis_type_creator.create_analysis_types()

    def _there_are_footprints_of_this_clade(self, clade_fp_dict):
        return clade_fp_dict

    def _populate_clade_fp_dicts_list(self):
        footprint_dict_pop_handler = FootprintDictPopHandler(sp_data_analysis_parent=self)
        footprint_dict_pop_handler.populate_clade_footprint_dicts()

    def _setup_temp_wkd(self):
        if os.path.exists(self.temp_wkd):
            shutil.rmtree(self.temp_wkd)
        os.makedirs(self.temp_wkd, exist_ok=True)

class CladeCollectionInfoHoder:
    def __init__(
            self, total_seq_abundance, footprint_as_frozen_set_of_ref_seq_uids, ref_seq_rel_abund_dict, clade):
        self.total_seq_abundance = total_seq_abundance
        self.footprint_as_frozen_set_of_ref_seq_uids = footprint_as_frozen_set_of_ref_seq_uids
        self.ref_seq_rel_abund_dict = ref_seq_rel_abund_dict
        self.analysis_type_uid = []
        self.clade = clade

class ArtefactAssessor:
    def __init__(self, parent_sp_data_analysis):
        self.parent = parent_sp_data_analysis
        # TODO we already have three cc.id dicts and we are about to add a fourth for the footprint as frozen set
        # it seems to me that we should have these in a single dictionary
        # this will make it easier to keep track of too.
        # key:cc.id, value:int(total_sequences_in_clade_collection)
        self.cc_info_dict = self._create_cc_info_dict()
        self.cc_to_total_seqs_dict = self._create_total_seqs_dict_for_all_ccs()
        # key:cc.id, value:dict(key:ref_seq_obj_of_cc, value:int(abundance_of_the_ref_seq))
        self.cc_to_ref_seq_abunds_dict = self._create_ref_seq_rel_abunds_for_all_ccs_dict()
        self.analysis_types_of_analysis = AnalysisType.objects.filter(data_analysis_from=self.parent.data_analysis_obj)
        self.set_of_clades_from_analysis = self._set_set_of_clades_from_analysis()
        # key:AnalysisType.id, value:AnalysisTypeAretfactInfoHolder
        self.analysis_type_artefact_info_dict = self._create_analysis_type_aretefact_info_dict()
        # set(frozenset(AnalysisTypeArtefactInfoHolder.ref_seq_uids_set) for AnalysisTypeArtefactInfoHolder in
        # self.analysis_teyp_artefact_info_dict
        self.analysis_type_footprint_set = self._init_analysis_type_footprint_set()
        # key:CladeCollection.id, value:list(associated AnalysisType uids)
        self.cc_id_to_analysis_type_id_dict = self._create_clade_collection_to_starting_analysis_type_dictionary()

        # attributes updated on an iterative basis
        self.current_clade = None
        # NB we have the two lists below as we only want to check combinations of the original AnalysiTypes and
        # not the new AnalysisTypes that will be created as part of this process. This is to prevent any infinite
        # loops occuring.
        # A query that will be coninually updated
        self.types_of_clade_dynamic = None
        # A fixed list of the original types that we stated with
        self.types_of_clade_static = None
        # A list that holds a tuple of ids that have already been compared.
        self.already_compared_analysis_type_uid_set = set()
        # Bool whether the pair comparisons need to be restarted.
        # This will be true when we have modified a type in anyway
        self.restart_pair_comparisons = None

    def _create_cc_info_dict(self):
        print('Generating CladeCollection info objects')
        cc_input_mp_queue = Queue()
        mp_manager = Manager()
        cc_to_info_items_mp_dict = mp_manager.dict()

        for cc in self.parent.ccs_of_analysis:
            cc_input_mp_queue.put(cc)

        for n in range(self.parent.parent.args.num_proc):
            cc_input_mp_queue.put('STOP')

        all_processes = []
        for n in range(self.parent.parent.args.num_proc):
            p = Process(target=self._cc_id_to_cc_info_obj_worker, args=(cc_input_mp_queue, cc_to_info_items_mp_dict))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

        return dict(cc_to_info_items_mp_dict)

    def _cc_id_to_cc_info_obj_worker(self, cc_input_mp_queue, cc_to_info_items_mp_dict):
        for clade_collection_object in iter(cc_input_mp_queue.get, 'STOP'):
            sys.stdout.write(f'\r{clade_collection_object} {current_process().name}')

            dss_objects_of_cc_list = list(DataSetSampleSequence.objects.filter(
                clade_collection_found_in=clade_collection_object))

            list_of_ref_seqs_in_cc = [
                dsss.reference_sequence_of for dsss in dss_objects_of_cc_list]

            total_sequences_in_cladecollection = sum([dsss.abundance for dsss in dss_objects_of_cc_list])

            list_of_rel_abundances = [dsss.abundance / total_sequences_in_cladecollection for dsss in
                                  dss_objects_of_cc_list]

            ref_seq_frozen_set = frozenset(dsss.reference_sequence_of.id for dsss in dss_objects_of_cc_list)

            ref_seq_rel_abund_dict = {}
            for i in range(len(dss_objects_of_cc_list)):
                ref_seq_rel_abund_dict[list_of_ref_seqs_in_cc[i]] = list_of_rel_abundances[i]

            cc_to_info_items_mp_dict[clade_collection_object.id] = CladeCollectionInfoHoder(
                clade=clade_collection_object.clade,
                footprint_as_frozen_set_of_ref_seq_uids=ref_seq_frozen_set,
                ref_seq_rel_abund_dict=ref_seq_rel_abund_dict,
                total_seq_abundance=total_sequences_in_cladecollection)


    def _init_analysis_type_footprint_set(self):
        return set(
            [frozenset(anal_info_obj.ref_seq_uids_set) for anal_info_obj in self.analysis_type_artefact_info_dict])

    def _types_should_be_checked(self, artefact_info_a, artefact_info_b):
        """Check
        1 - the non-artefact_ref_seqs match
        2 - neither of the types is a subset of the other
        3 - there are artefact ref seqs in at least one of the types"""
        if artefact_info_a.non_artefact_ref_seq_uids_set == artefact_info_b.non_artefact_ref_seq_uids_set:
            if not set(artefact_info_a.ref_seq_uids_set).issubset(artefact_info_b.ref_seq_uids_set):
                if not set(artefact_info_b.ref_seq_uids_set).issubset(artefact_info_a.ref_seq_uids_set):
                    if artefact_info_a.artefact_ref_seq_uids_set.union(artefact_info_b.artefact_ref_seq_uids_set):
                        return True
        return False

    def assess_within_clade_cutoff_artefacts(self):
        for clade in self.set_of_clades_from_analysis:
            self.current_clade = clade
            self.types_of_clade_dynamic = AnalysisType.objects.filter(
                data_analysis_from=self.parent.data_analysis_obj,
                clade=self.current_clade)
            self.types_of_clade_static = list(AnalysisType.objects.filter(
                data_analysis_from=self.parent.data_analysis_obj,
                clade=self.current_clade))
            while 1:
                self.restart_pair_comparisons = False
                for analysis_type_a, analysis_type_b in itertools.combinations(
                    [at for at in self.types_of_clade_dynamic if at in self.types_of_clade_static], 2):
                    if {analysis_type_a.id, analysis_type_b.id} not in self.already_compared_analysis_type_uid_set:
                        artefact_info_a = self.analysis_type_artefact_info_dict[analysis_type_a.id]
                        artefact_info_b = self.analysis_type_artefact_info_dict[analysis_type_b.id]
                        if self._types_should_be_checked(artefact_info_a, artefact_info_b):

                            print('\nChecking {} and {} for additional artefactual profiles'.format(analysis_type_a, analysis_type_b))
                            ctph = CheckTypePairingHandler(parent_artefact_assessor=self)
                            ctph.check_type_pairing()
                            if check_type_pairing_for_artefact_type(
                                    analysis_type_a, analysis_type_b, type_footprint_dict, clade, cc_to_initial_type_dict,
                                    cc_to_ref_seq_list_and_abundances, num_processors):

                                self.restart_pair_comparisons = True

                                self.types_of_clade_dynamic = AnalysisType.objects.filter(
                                    data_analysis_from=self.parent.data_analysis_obj, clade=clade)
                                break
                            else:
                                self.already_compared_analysis_type_uid_set.add(
                                    frozenset({analysis_type_a.id, analysis_type_b.id}))
                        else:
                            self.already_compared_analysis_type_uid_set.add(
                                frozenset({analysis_type_a.id, analysis_type_b.id}))



                    # If we make it to here we did a full run through the types without making a new type
                    # Time to do the same for the next type
                if not self.restart_pair_comparisons:
                    break


    def _type_non_artefacts_are_subsets_of_each_other(self, analysis_type_a, analysis_type_b):
        if set(self.analysis_type_artefact_info_dict[
                            analysis_type_a.id].ref_seq_uids_set).issubset(
            self.analysis_type_artefact_info_dict[analysis_type_b.id].ref_seq_uids_set):
            if set(self.analysis_type_artefact_info_dict[
                            analysis_type_b.id].ref_seq_uids_set).issubset(
                    self.analysis_type_artefact_info_dict[analysis_type_a.id].ref_seq_uids_set):
                return True
        return False

    def _create_clade_collection_to_starting_analysis_type_dictionary(self):
        """Create a dictionary that links CladeCollection to the AnalysisTypes that they associated with.
        NB it used to be that only one AnalysisType could be associated to one CladeCollection, but since
        the introduction of the basal type theory a single CladeCollection can now associate with more than one
        AnalysisType. As such, we use a defaultdict(list) as the default rather than a direct key to value association
        between CladeCollection and AnlaysisType."""
        cc_to_initial_type_dict = defaultdict(list)

        clade_collection_to_type_tuple_list = []
        for at in self.analysis_types_of_analysis:
            type_uid = at.id
            initial_clade_collections = [int(x) for x in at.list_of_clade_collections_found_in_initially.split(',')]
            for CCID in initial_clade_collections:
                clade_collection_to_type_tuple_list.append((CCID, type_uid))

        for ccid, atid in clade_collection_to_type_tuple_list:
            cc_to_initial_type_dict[ccid].append(atid)

        return dict(cc_to_initial_type_dict)

    def _create_analysis_type_aretefact_info_dict(self):
        """Create a dict for each of the AnalysisTypes that will be kept updated throughout the artefact checking.
        The dict will be key AnalysisType.id, value will be an AnalysisTypeAretefactInfoHolder."""
        analysis_type_artefact_info_dict = {}
        for at in self.analysis_types_of_analysis:

            ref_seqs_uids_of_analysis_type_set = set([ref_seq.id for ref_seq in at.get_ordered_footprint_list()])
            artefact_ref_seq_uids_set = set([int(x) for x in at.artefact_intras.split(',') if x != ''])
            non_artefact_ref_seq_uids_set = set([uid for uid in ref_seqs_uids_of_analysis_type_set if uid not in artefact_ref_seq_uids_set])
            footprint_as_ref_seq_objects = at.get_ordered_footprint_list()
            analysis_type_artefact_info_dict[at.id] = AnalysisTypeArtefactInfoHolder(
                artefact_ref_seq_uids_set=artefact_ref_seq_uids_set,
                non_artefact_ref_seq_uids_set=non_artefact_ref_seq_uids_set,
                ref_seq_uids_set=ref_seqs_uids_of_analysis_type_set,
                footprint_as_ref_seq_objects=footprint_as_ref_seq_objects,
                basal_seqs_set=self._generate_basal_seqs_set(footprint=footprint_as_ref_seq_objects),
                clade=at.clade)

        return analysis_type_artefact_info_dict

    def _generate_basal_seqs_set(self, footprint):
        basal_set = set()
        found_c15_a = False
        for rs in footprint:
            if rs.name == 'C3':
                basal_set.add('C3')
            elif rs.name == 'C1':
                basal_set.add('C1')
            elif 'C15' in rs.name and not found_c15_a:
                basal_set.add('C15')
                found_c15_a = True

        if basal_set:
            return basal_set
        else:
            return None

    def _set_set_of_clades_from_analysis(self):
        self.set_of_clades_from_analysis = set()
        for at in self.analysis_types_of_analysis:
            self.set_of_clades_from_analysis.add(at.clade)
        return self.set_of_clades_from_analysis

    def _create_ref_seq_rel_abunds_for_all_ccs_dict(self):
        cc_input_mp_queue = Queue()
        mp_manager = Manager()
        cc_to_ref_seq_abunds_mp_dict = mp_manager.dict()

        for cc in self.parent.ccs_of_analysis:
            cc_input_mp_queue.put(cc)

        for n in range(self.parent.parent.args.num_proc):
            cc_input_mp_queue.put('STOP')

        all_processes = []
        for n in range(self.parent.parent.args.num_proc):
            p = Process(target=self._cc_to_ref_seq_list_and_abund_worker, args=(cc_input_mp_queue, cc_to_ref_seq_abunds_mp_dict))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

        return dict(cc_to_ref_seq_abunds_mp_dict)

    def _cc_to_ref_seq_list_and_abund_worker(self, cc_input_mp_queue, cc_to_ref_seq_abunds_mp_dict):
        for clade_collection_object in iter(cc_input_mp_queue.get, 'STOP'):

            sys.stdout.write(f'\r{clade_collection_object} {current_process().name}')

            list_of_dataset_sample_sequences_in_clade_collection = [
                dsss for dsss in DataSetSampleSequence.objects.filter(
                    clade_collection_found_in=clade_collection_object)]

            list_of_ref_seqs_in_cladecollection = [
                dsss.reference_sequence_of for dsss in list_of_dataset_sample_sequences_in_clade_collection]

            list_of_abundances = [dsss.abundance / self.cc_to_total_seqs_dict[clade_collection_object.id] for dsss in
                                  list_of_dataset_sample_sequences_in_clade_collection]
            inner_dict = {}
            for i in range(len(list_of_dataset_sample_sequences_in_clade_collection)):
                inner_dict[list_of_ref_seqs_in_cladecollection[i]] = list_of_abundances[i]
            cc_to_ref_seq_abunds_mp_dict[clade_collection_object.id] = inner_dict

    def _create_total_seqs_dict_for_all_ccs(self):
        """Generate a dict that simply holds the total number of seqs per clade_collection_object.
        This will be used when working out relative proportions of seqs in the clade_collection_object
        """
        print('Generating CladeCollection to total sequence dictionary')
        cc_input_mp_queue = Queue()
        mp_manager = Manager()
        cc_to_total_seqs_mp_dict = mp_manager.dict()

        for cc in self.parent.ccs_of_analysis:
            cc_input_mp_queue.put(cc)

        for n in range(self.parent.parent.args.num_proc):
            cc_input_mp_queue.put('STOP')

        all_processes = []
        for n in range(self.parent.parent.args.num_proc):
            p = Process(target=self._cc_to_total_seqs_dict_worker, args=(cc_input_mp_queue, cc_to_total_seqs_mp_dict))
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

        return dict(cc_to_total_seqs_mp_dict)

    def _cc_to_total_seqs_dict_worker(self, cc_input_mp_queue, cc_to_total_seqs_mp_dict):
        for clade_collection_object in iter(cc_input_mp_queue.get, 'STOP'):
            sys.stdout.write(f'\r{clade_collection_object} {current_process().name}')
            total_sequences_in_cladecollection = sum(
                [dsss.abundance for dsss in
                 DataSetSampleSequence.objects.filter(clade_collection_found_in=clade_collection_object)])
            cc_to_total_seqs_mp_dict[clade_collection_object.id] = total_sequences_in_cladecollection

class PotentialNewType:
    def __init__(self, artefact_ref_seq_set, non_artefact_ref_seq_set, ref_seq_uids_set, name):
        self.artefact_ref_seq_set = artefact_ref_seq_set
        self.non_artefact_ref_seq_set = non_artefact_ref_seq_set
        self.ref_seq_uids_set = ref_seq_uids_set
        self.name = name

class CheckTypePairingHandler:
    def __int__(self, parent_artefact_assessor, artefact_info_a, artefact_info_b):
        self.parent = parent_artefact_assessor
        # AnalysisTypeAretefactInfoHolder for types a and b
        self.info_a = artefact_info_a
        self.info_b = artefact_info_b
        # NB that the non_artefact ref seqs being the same is a prerequisite of doing a comparison
        # as such we can set the below pnt non artefact ref seqs to match just one of the types being compared
        self.pnt = self._init_pnt(self.info_a, self.info_b)
        self.list_of_ccs_to_check = None

    def _init_pnt(self, artefact_info_a, artefact_info_b):
        return PotentialNewType(
            ref_seq_uids_set=self.info_a.ref_seq_uids_set.union(self.info_b.artefact_ref_seq_uids_set),
            artefact_ref_seq_set=self.info_a.artefact_ref_seq_uids_set.union(self.info_b.artefact_ref_seq_uids_set),
            non_artefact_ref_seq_set=self.info_a.non_artefact_ref_seq_uids_set,
            name=','.join(
                [ref_seq.name for ref_seq in
                 set(list(
                     artefact_info_a.footprint_as_ref_seq_objects +
                     artefact_info_b.footprint_as_ref_seq_objects
                 ))]))

    def check_type_pairing(self):
        if self._pnt_profile_already_an_existing_analysis_type_profile():
            print(f'Assessing new type:{pnt.name}')
            print('Potential new type already exists')
            return False

        # TODO This is perhpas stupid but I'm going to remove the second conditional from this
        # as we were supposed to be checking all of the ccs
        self.list_of_ccs_to_check = [cc for cc in self.parent.parent.ccs_of_analysis if
                                     cc.clade == self.info_a.clade if
                                     cc.id in self.parent.cc_id_to_analysis_type_id_dict]

        print(f'Assessing support for potential new type:{self.pnt.name}')

    def _pnt_profile_already_an_existing_analysis_type_profile(self):
        return self.pnt.ref_seq_uids_set in self.parent.analysis_type_footprint_set

class AnalysisTypeArtefactInfoHolder:
    """An object to aid in fast access of the following information for each analysis type"""
    def __init__(self, artefact_ref_seq_uids_set, non_artefact_ref_seq_uids_set, ref_seq_uids_set, footprint_as_ref_seq_objects, basal_seqs_set, clade):
        self.artefact_ref_seq_uids_set = artefact_ref_seq_uids_set
        self.non_artefact_ref_seq_uids_set = non_artefact_ref_seq_uids_set
        self.ref_seq_uids_set = ref_seq_uids_set
        self.footprint_as_ref_seq_objects = footprint_as_ref_seq_objects
        self.basal_seqs_set = basal_seqs_set
        self.clade = clade

class AnalysisTypeCreator:
    """Create AnalysisType objects from the supported initial type profiles that have been generated in the
    SupportedFootprintIdentifier"""
    def __init__(self, parent_sp_data_analysis):
        self.parent = parent_sp_data_analysis

    def create_analysis_types(self):
        print(f'\n\nCreating analysis types clade {self.parent.current_clade}')
        for initial_type in self.parent.collapsed_footprint_dict:
            if self._initial_type_is_codom(initial_type):
                self._create_new_analysis_type_co_dom(initial_type)
            else:
                self._create_new_analysis_type_non_co_dom(initial_type)

    def _create_new_analysis_type_non_co_dom(self, initial_type):
        new_analysis_type = AnalysisType(co_dominant=False, data_analysis_from=self.parent.data_analysis_obj,
                                         clade=self.parent.current_clade)
        new_analysis_type.set_maj_ref_seq_set(initial_type.set_of_maj_ref_seqs)
        new_analysis_type.init_type_attributes(initial_type.clade_collection_list, initial_type.profile)
        new_analysis_type.save()
        sys.stdout.write(f'\rCreating analysis type: {new_analysis_type.name}')

    def _create_new_analysis_type_co_dom(self, initial_type):
        new_analysis_type = AnalysisType(
            co_dominant=True,
            data_analysis_from=self.parent.data_analysis_obj,
            clade=self.parent.current_clade)
        new_analysis_type.set_maj_ref_seq_set(initial_type.set_of_maj_ref_seqs)
        new_analysis_type.init_type_attributes(initial_type.clade_collection_list, initial_type.profile)
        new_analysis_type.save()
        print('\rCreating analysis type: {}'.format(new_analysis_type.name), end='')

    def _initial_type_is_codom(self, initial_type):
        return len(initial_type.set_of_maj_ref_seqs) > 1

class FootprintDictPopHandler:
    """Will handle the execusion of the FootprintDictWorker. This worker will populate the
    SPDataAnalysis master_cladal_list_of_footpinrt_dicts"""
    def __init__(self, sp_data_analysis_parent):
        self.parent = sp_data_analysis_parent
        self.cc_mp_queue = Queue()
        self.output_mp_queue = Queue()
        self._populate_queue()
        self.all_procs = []

    def _populate_queue(self):
        for cc in self.parent.ccs_of_analysis:
            self.cc_mp_queue.put(cc)

        for N in range(self.parent.parent.args.num_proc):
            self.cc_mp_queue.put('STOP')

    def populate_clade_footprint_dicts(self):

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        for n in range(self.parent.parent.args.num_proc):
            p = Process(target=self._start_footprint_dict_workers,
                        args=())
            self.all_procs.append(p)
            p.start()

        self._collect_cc_info()

    def _collect_cc_info(self):
        kill_number = 0
        while 1:
            cc_info_holder = self.output_mp_queue.get()
            if cc_info_holder == 'kill':
                kill_number += 1
                if kill_number == self.parent.parent.args.num_procs:
                    break
            else:
                if self._footprint_in_dict_already(cc_info_holder):
                    self._add_cc_and_maj_seq_to_clade_existant_fp_dict(cc_info_holder)
                else:
                    self._add_cc_and_maj_seq_to_new_fp_dict(cc_info_holder)


        for p in self.all_procs:
            p.join()

    def _add_cc_and_maj_seq_to_new_fp_dict(self, cc_info_holder):
        self.parent.clade_footp_dicts_list[
            cc_info_holder.clade_index][
            cc_info_holder.footprint] = FootprintRepresentative(
            cc=cc_info_holder.cc, cc_maj_ref_seq=cc_info_holder.maj_ref_seq)

    def _add_cc_and_maj_seq_to_clade_existant_fp_dict(self, cc_info_holder):
        self.parent.clade_footp_dicts_list[cc_info_holder.clade_index][
            cc_info_holder.footprint].cc_list.append(
            cc_info_holder.cc)
        self.parent.clade_footp_dicts_list[cc_info_holder.clade_index][
            cc_info_holder.footprint].maj_seq_list.append(
            cc_info_holder.maj_ref_seq)

    def _footprint_in_dict_already(self, cc_info_holder):
        return cc_info_holder.footprint in self.parent.clade_footp_dicts_list[
                    cc_info_holder.clade_index]

    def _start_footprint_dict_workers(self):
        for cc in iter(self.cc_mp_queue.get, 'STOP'):

            footprint = cc.cutoff_footprint(self.parent.parent.data_analysis_obj.within_clade_cutoff)

            self.output_mp_queue.put(CCInfoHolder(
                footprint=footprint,
                clade_index=self.parent.parent.clade_list.index(cc.clade), cc=cc,
                maj_ref_seq=cc.maj()))

            sys.stdout.write(f'\rFound footprint {footprint}')

        self.output_mp_queue.put('kill')


class SupportedFootPrintIdentifier:
    """This class is responsible for identifying the footprints that are found in a sufficient number of clade
    collections to warrant becoming AnalysisType objects.
    Operates by working with the longest footprints and trying to collapse those that aren't already supported into
    those of length n-1 etc etc. It also tries to find support amongst the collection of footprints of length n
    e.g. if you have 1-2-3-4, 1-2-3-5, 1-2-3-6, 1-2-3-7, 1-2-3-8, it will pull out the 1-2-3 footprint as supported
    even if the 1-2-3 footprint doesnt already exist as an n=3 footprint.
    We are also taking into account not allowing analysis types to contain the C3 and the C15 sequences. We refer to
    these as the basal sequences. All footprints that contain more than one basal sequencs ill automatically be
    put into the unsupported list at first, irrespective of their support. Then they will attempt to be collapsed into
    other footprints just like all other unsupported footprints. If a suitable footprint is foud that they can be
    collapsed into then they will be but the pulled out sequqences can only contain one of the basal
    sequences. In this way, a single clade collection can go towards the support of two different
    footprints, one for with C3 and one with C15.
    as """
    def __init__(self, clade_footprint_dict, parent_sp_data_analysis):
        self.parent = parent_sp_data_analysis
        self.clade_fp_dict = clade_footprint_dict
        self.supported_list = []
        self.unsupported_list = []
        self.initial_types_list = []
        self._init_initial_types_list()
        # number of clade collections a footprint must be found in to be supported
        self.required_support = 4

        # arguments that are used in the _populate_collapse_dict_for_next_n_mehod
        # nb these arguments are regularly updated
        # Bool that represents whether we should conitnue to iterate through the larger types trying to find
        # shorter types to collapse into
        self.repeat = None
        # list holding the initial types that are length n-1
        self.n_minu_one_list = []
        # key = big initial type, value = small initial type that it should be collapsed into
        self.collapse_dict = {}
        # the number of support that a potetially collapsed type, i.e. support of the large initial type
        # plus the cc support of the small inital type it will be collapsed into
        # this is used to assesss which collapse will happen in the case that there are several viable collapse
        # options. The bigest score will be collapsed.
        self.top_score = 0
        # the size of initial type footprint we are currnetly working with
        self.current_n = None
        self.large_fp_to_collapse_list = None

        # Attributes used in the synthetic footprint generation
        self.list_of_initial_types_of_len_n = None
        self.synthetic_fp_dict = None

        # Attributes used in synthetic footprint collapse
        self.ordered_sig_synth_fps = None

    def _init_initial_types_list(self):
        for footprint_key, footprint_representative in self.clade_fp_dict.items():
            self.initial_types_list.append(
                InitialType(
                    footprint_key, footprint_representative.cc_list, footprint_representative.maj_seq_list))

    def identify_supported_footprints(self):
        # for each length starting at max and dropping by 1 with each increment
        longest_footprint = self._return_len_of_longest_fp()
        for n in range(longest_footprint, 0, -1):
            self.current_n = n
            self._update_supported_unsupported_lists_for_n()

            self._collapse_n_len_initial_types_into_minus_one_initial_types()

            if self.current_n >2:
                # generate insilico intial types and try to collapse to these
                # We only need to attempt the further collapsing if there are unsupported types to collapse
                if len(self.unsupported_list) > 1:
                    # For every initial type (supported and unsupported) of length n,
                    # generate all the permutations of the reference sequnces that are n-1 in length
                    # Then try to collapse into these.
                    self.list_of_initial_types_of_len_n = [
                        initial_t for initial_t in self.initial_types_list if initial_t.profile_length == self.current_n]
                    # Only carry on if we have lengthN footprints to get sequences from
                    if self.list_of_initial_types_of_len_n:
                        self._generate_synth_footprints()
                        self._associate_un_sup_init_types_to_synth_footprints()


                        # Here we have a populated dict.
                        # We are only interseted in n-1 len footprints that were found in the unsupported types
                        # Because each of the synth footprints will be found in the
                        # initial_types they originated from we only need concern ourselves with synthetic
                        # footprints associated with more than 1 cladecollection
                        sig_synth_fps = [
                            kmer for kmer in self.synthetic_fp_dict.keys() if len(self.synthetic_fp_dict[kmer]) > 1]
                        if sig_synth_fps:
                            # parse through the synth fps in order of the number of initial types they associated with
                            self.ordered_sig_synth_fps = sorted(
                                sig_synth_fps, key=lambda x: len(self.synthetic_fp_dict[x]), reverse=True)
                            sfc = SyntheticFootprintCollapser(parent_supported_footprint_identifier=self)
                            sfc.collapse_to_synthetic_footprints()

            else:
                uffc = UnsupportedFootprintFinalCollapser(parent_supported_footprint_identifier=self)
                uffc.collapse_unsup_footprints_to_maj_refs()
        return self.initial_types_list

    def _collapse_type_to_matching_init_type_if_exists(self, k, synth_fp):
        if self._extracted_initial_type_fp_now_matches_existing_initial_type(
                synth_fp=synth_fp, init_type_index=k):
            self._absorb_matching_init_type_and_delete(k, synth_fp)
        self._eval_remove_type_from_unsup_list(k, synth_fp)

    def _del_type_to_collapse(self, k, synth_fp):
        # then the big init_type no longer contains any
        # of its original maj ref seqs and
        # so should be delted.
        self.initial_types_list.remove(self.synthetic_fp_dict[synth_fp][k])
        self.unsupported_list.remove(self.synthetic_fp_dict[synth_fp][k])

    def _eval_remove_type_from_unsup_list(self, k, synth_fp):
        """We now need to decide if the footprint to collapse should be removed from the unsupported_list.
        This will depend on if it is longer than n or not.
        """
        if self.synthetic_fp_dict[synth_fp][
            k].profile_length < self.current_n:
            self.unsupported_list.remove(
                self.synthetic_fp_dict[synth_fp][k])

    def _absorb_matching_init_type_and_delete(self, k, synth_fp):
        """Here we have found an intial type that exactly matches the initial type to
        collapse's new footprint (i.e. after extraction).
        Now absorb the found match initial type in to the initial type to be collapsed.
        We do it this way around the initial type to collapse stays in the corect place
        i.e. in the unsupported list of not.
        After absorption, remove the matched initial type from the initial types list and unsupported list
        """
        self.synthetic_fp_dict[synth_fp][
            k].absorb_large_init_type(
            self.initial_types_list[j])
        if self.initial_types_list[j] in self.unsupported_list:
            self.unsupported_list.remove(self.initial_types_list[j])
        self.initial_types_list.remove(self.initial_types_list[j])

    def _synth_fp_matches_existing_intial_type(self, synth_fp):
        """If an non-synthetic init_type already exists with the same footprint as the synthetic footprint in
        question then add the init_type to be collapsed to it rather than creating a new initial type from the
        synthetic footprint.
        """
        for i in range(len(self.initial_types_list)):
            # for initT_one in initial_types_list:
            if self.initial_types_list[i].profile == synth_fp:
                return True
        return False

    def _extracted_initial_type_fp_now_matches_existing_initial_type(self, synth_fp, init_type_index):
        """If the initial type to collapse still contains maj
        ref sequences, then it is still a profile that we needs to be assessed for collapse.
        Now need to check if it's new profile (i.e. after extraction) matches that of any of the other initial types."""
        for j in range(len(self.initial_types_list)):
            if self.initial_types_list[j].profile == self.synthetic_fp_dict[synth_fp][init_type_index].profile:
                if self.initial_types_list[j] != self.synthetic_fp_dict[synth_fp][init_type_index]:
                    return True
        return False

    def _validate_init_type_is_unsup_and_synth_fp_is_subset(self, k, synth_fp):
        """Then this footprint hasn't been collapsed anywhere yet.
        We also check to make sure that the initial type's profile hasn't been modified (i.e. by extraction) and check
        that the synthetic fp is still a sub set of the initial type's profile.
        """
        if self.synthetic_fp_dict[synth_fp][k] in self.unsupported_list:
            if synth_fp.issubset(self.synthetic_fp_dict[synth_fp][k].profile):
                return True
        return False

    def _generate_synth_footprints(self):
        sfp = SyntheticFootprintPermutator(parent_supported_footprint_identifier=self)
        sfp.permute_synthetic_footprints()
        # NB its necessary to convert the mp dict to standard dict for lists that are the values
        # to behave as expected with the .append() function.
        self.synthetic_fp_dict = dict(sfp.collapse_n_mer_mp_dict)

    def _associate_un_sup_init_types_to_synth_footprints(self):
        # Now go through each of the (n-1) footprints and see if they
        # fit into a footprint in the unsuported list
        print('Checking new set of synthetic types')
        for synth_fp in self.synthetic_fp_dict.keys():

            if self._does_synth_fp_have_multi_basal_seqs(synth_fp):
                continue

            for un_sup_initial_type in self.unsupported_list:
                # For each of the synth footprints see if they fit within the unsupported types.
                # If so then add the initial type into the list associated with that synth footprint
                # in the synthetic_fp_dict.
                # For a match, at least one maj ref seqs need to be in common between the two footprints
                if synth_fp.issubset(un_sup_initial_type.profile):
                    if len(un_sup_initial_type.set_of_maj_ref_seqs & synth_fp) >= 1:
                        # Then associate un_sup_initial_type to the synth fp
                        self.synthetic_fp_dict[synth_fp].append(un_sup_initial_type)

    def _does_synth_fp_have_multi_basal_seqs(self, frozen_set_of_ref_seqs):
        basal_count = 0
        c15_found = False
        for ref_seq in frozen_set_of_ref_seqs:
            if 'C15' in ref_seq.name and not c15_found:
                basal_count += 1
                c15_found = True
                continue
            elif ref_seq.name == 'C3':
                basal_count += 1
                continue
            elif ref_seq.name == 'C1':
                basal_count += 1
                continue
        if basal_count > 1:
            return True
        else:
            return False

    def _collapse_n_len_initial_types_into_minus_one_initial_types(self):
        # Try to collapse length n footprints into size n-1 footprints
        # we will try iterating this as now that we have the potential to find two types in one profile, e.g.
        # a C15 and C3, we may only extract the C3 on the first iteration but there may still be a C15 in initial
        # type.
        repeat = True
        while repeat:
            self._populate_collapse_dict_for_next_n()

            self.large_fp_to_collapse_list = list(self.collapse_dict.keys())

            for q in range(len(self.large_fp_to_collapse_list)):
                large_fp_to_collapse = self.large_fp_to_collapse_list[q]

                if self._remove_fp_to_collapse_from_unsupported_if_now_supported(large_fp_to_collapse):
                    continue

                fp_collapser = FootprintCollapser(
                    footprint_to_collapse_index=q,
                    parent_supported_footprint_identifier=self)
                fp_collapser.collapse_footprint()

    def _remove_fp_to_collapse_from_unsupported_if_now_supported(self, large_fp_to_collapse):
        if large_fp_to_collapse.support >= \
                self.required_support and not large_fp_to_collapse.contains_multiple_basal_sequences:
            # Then this type has already had some other leftovers put into it so that it now has the required
            # support. In this case we can remove the type from the unsupported list and continue to the next
            self.unsupported_list.remove(large_fp_to_collapse)
            return True
        else:
            return False

    def _populate_collapse_dict_for_next_n(self):
        self._set_attributes_for_collapse_dict_population()
        if self.n_minus_one_list:
            for longer_initial_type in self.unsupported_list:
                sys.stdout.write(f'Assessing footprint {longer_initial_type} for supported type\r')
                self.top_score = 0
                for shorter_initial_type in self.n_minus_one_list:
                    collapse_assessor = CollapseAssessor(
                        parent_supported_footprint_identifier=self,
                        longer_intial_type=longer_initial_type,
                        shorter_initial_type=shorter_initial_type)
                    collapse_assessor.assess_collapse()


    def _set_attributes_for_collapse_dict_population(self):
        self.collapse_dict = {}
        self.repeat = False
        self.n_minus_one_list = self._get_n_minus_one_list()

    def _get_n_minus_one_list(self):
        n_minus_one_list = [initial_type for initial_type in self.initial_types_list if
                            initial_type.profile_length == self.current_n - 1]
        return n_minus_one_list

    def _return_len_of_longest_fp(self):
        longest_footprint = max([initial_type.profile_length for initial_type in self.initial_types_list])
        return longest_footprint

    def _update_supported_unsupported_lists_for_n(self):
        # populate supported and unsupported list for the next n
        n_list = [
            initial_type for initial_type in self.initial_types_list if initial_type.profile_length == self.current_n]

        for initial_type in n_list:
            if initial_type.support >= self.required_support and not initial_type.contains_multiple_basal_sequences:
                self.supported_list.append(initial_type)
            else:
                self.unsupported_list.append(initial_type)

class UnsupportedFootprintFinalCollapser:
    def __init__(self, parent_supported_footprint_identifier):
        self.parent = parent_supported_footprint_identifier
        self.fp_to_collapse = None
        self.matching_initial_type = None

    def collapse_unsup_footprints_to_maj_refs(self):
        while self.parent.unsupported_list:
            self.fp_to_collapse = self.parent.unsupported_list[0]
            for i in range(len(self.fp_to_collapse.clade_collection_list)):  # for each cc
                for maj_dss in self.fp_to_collapse.majority_sequence_list[i]:  # for each maj_ref_seq
                    if self._inital_type_exists_with_maj_ref_seq_as_profile(maj_dss):
                        self._add_unsup_type_info_to_smll_match_type(i, maj_dss)
                    else:
                        self._create_new_maj_seq_init_type(i, maj_dss)

            self._del_type_to_collapse()

    def _create_new_maj_seq_init_type(self, i, maj_dss):
        new_initial_type = InitialType(
            reference_sequence_set=frozenset([maj_dss.reference_sequence_of]),
            clade_collection_list=[self.fp_to_collapse.clade_collection_list[i]])
        self.parent.initial_types_list.append(new_initial_type)

    def _add_unsup_type_info_to_smll_match_type(self, i, maj_dss):
        self.matching_initial_type.clade_collection_list.append(self.fp_to_collapse.clade_collection_list[i])
        self.matching_initial_type.majority_sequence_list.append([maj_dss])

    def _inital_type_exists_with_maj_ref_seq_as_profile(self, maj_dss):
        """Check to see if an initial type with profile of that maj_dss refseq already exists"""
        for initT in [init for init in self.parent.initial_types_list if init.profile_length == 1]:
            if maj_dss.reference_sequence_of in initT.profile:
                self.matching_initial_type = initT
                return True
        return False

    def _del_type_to_collapse(self):
        """Once we have associated each of the cc of an initial type to collapse to existing or new initial types,
        delete.
        """
        self.parent.initial_types_list.remove(self.fp_to_collapse)
        if self.fp_to_collapse in self.parent.unsupported_list:
            self.parent.unsupported_list.remove(self.fp_to_collapse)

class CollapseAssessor:
    """Responsible for assessing whether an unsupported large initial type can be collapsed into a given
    small initial type.
    """
    def __init__(self, parent_supported_footprint_identifier, longer_intial_type, shorter_initial_type):
        self.parent = parent_supported_footprint_identifier
        self.longer_intial_type = longer_intial_type
        self.shorter_initial_type = shorter_initial_type

    def assess_collapse(self):
        # see docstring of method for more info
        if self._if_short_initial_type_suitable_for_collapse():
            if self.longer_intial_type.contains_multiple_basal_sequences:
                if self.does_small_footprint_contain_the_required_ref_seqs_of_the_large_footprint():
                    # score = number of samples big was found in plus num samples small was found in
                    self._if_highest_score_so_far_assign_big_fp_to_smll_fp_for_collapse()

            else:
                if self.longer_intial_type.set_of_maj_ref_seqs.issubset(self.shorter_initial_type.profile):
                    self._if_highest_score_so_far_assign_big_fp_to_smll_fp_for_collapse()

    def _if_highest_score_so_far_assign_big_fp_to_smll_fp_for_collapse(self):
        score = self.longer_intial_type.support + self.shorter_initial_type.support
        if score > self.parent.top_score:
            self.parent.top_score = score
            self.parent.repeat = True
            self.parent.collapse_dict[self.longer_intial_type] = self.shorter_initial_type

    def does_small_footprint_contain_the_required_ref_seqs_of_the_large_footprint(self):
        set_of_seqs_to_find = set()
        ref_seqs_in_big_init_type = list(self.longer_intial_type.set_of_maj_ref_seqs)
        if self.shorter_initial_type.basalSequence_list:
            if 'C15' in self.shorter_initial_type.basalSequence_list[0]:
                # then this is a C15x basal type and we will need to find all sequences that are not C1 or C3
                for ref_seq in ref_seqs_in_big_init_type:
                    if ref_seq.name in ['C1', 'C3']:
                        # then this is a squence we don't need to find
                        pass
                    else:
                        set_of_seqs_to_find.add(ref_seq)

            elif self.shorter_initial_type.basalSequence_list[0] == 'C1':
                # then this is a C1 basal type and we need to find all sequence that are not C15x or C3
                for ref_seq in ref_seqs_in_big_init_type:
                    if 'C15' in ref_seq.name or ref_seq.name == 'C3':
                        # then this is a squence we don't need to find
                        pass
                    else:
                        set_of_seqs_to_find.add(ref_seq)

            elif self.shorter_initial_type.basalSequence_list[0] == 'C3':
                # then this is a C3 basal type and we need to find all sequence that are not C15x or C1
                for ref_seq in ref_seqs_in_big_init_type:
                    if 'C15' in ref_seq.name or ref_seq.name == 'C1':
                        # then this is a squence we don't need to find
                        pass
                    else:
                        set_of_seqs_to_find.add(ref_seq)

            # Here we have the list of the ref_seqs that we need to find in the small_init_type.profile
            if set_of_seqs_to_find.issubset(self.shorter_initial_type.profile):
                return True
            else:
                return False
        else:
            # if the small_init_type doesn't contain a basal sequence sequence, then we need to find all of the seqs
            # in the big_intit_type.set_of_maj_ref_seqs that are not C15x, C1 or C3
            for ref_seq in ref_seqs_in_big_init_type:
                if 'C15' in ref_seq.name or ref_seq.name in ['C1', 'C3']:
                    # then this is a squence we don't need to find
                    pass
                else:
                    set_of_seqs_to_find.add(ref_seq)
            # Here we have the list of the ref_seqs that we need to find in the small_init_type.profile
            if set_of_seqs_to_find.issubset(self.shorter_initial_type.profile):
                return True
            else:
                return False


    def _if_short_initial_type_suitable_for_collapse(self):
        """Consider this for collapse only if the majsequences of the smaller are a subset of the maj sequences of
        the larger e.g. we don't want A B C D being collapsed into B C D when A is the maj of the first and B is a
        maj of the second simplest way to check this is to take the setOfMajRefSeqsLarge which is a set of all of the
        ref seqs that are majs in the cc that the footprint is found in and make sure that it is a subset of the
        smaller footprint in question.

        10/01/18 what we actually need is quite complicated. If the big type is not multi basal, then we have no
        problem and we need to find all of the set of maj ref seqs in the small profile but if the large type is
        multi basal then it gets a little more complicated if the large type is multi basal then which of its set of
        maj ref seqs we need to find in the small profile is dependent on what the basal seq of the smallfootprint is.

        If the small has no basal seqs in it, then we need to find every sequence in the large's set of maj ref seqs
        that is NOT a C15x, or the C3 or C1 sequences.

        If small basal = C15x then we need to find every one of the large's seqs that isn't C1 or C3

        If small basal = C1 then we need to find every on of the large's seqs that isn't C3 or C15x

        If small basal = C3 then we need to find every on of the large's seqs that isn't C1 or C15x we should
        put this decision into a new function.
        """

        multi_basal = self.shorter_initial_type.contains_multiple_basal_sequences
        return self.shorter_initial_type.profile.issubset(self.longer_intial_type.profile) and not multi_basal

class SyntheticFootprintCollapser:
    def __init__(self, parent_supported_footprint_identifier):
        self.parent = parent_supported_footprint_identifier
        self.current_synth_fp = None
        self.current_fp_to_collapse = None
        self.matching_existing_init_type_outer = None
        self.matching_existing_init_type_inner = None

    def collapse_to_synthetic_footprints(self):
        for synth_fp in self.parent.ordered_sig_synth_fps:  # for each synth_fp
            self.current_synth_fp = synth_fp
            for k in range(len(self.parent.synthetic_fp_dict[synth_fp])):  # for each assoc. initial type
                self.current_fp_to_collapse = self.parent.synthetic_fp_dict[synth_fp][k]
                if self._validate_init_type_is_unsup_and_synth_fp_is_subset():
                    if self._synth_fp_matches_existing_intial_type():
                        if self._should_extract_rather_than_absorb():
                            self.matching_existing_init_type_outer.extract_support_from_large_initial_type(
                                self.current_fp_to_collapse)
                            if self.current_fp_to_collapse.set_of_maj_ref_seqs:
                                self._collapse_type_to_matching_init_type_if_exists()
                            else:
                                self._del_type_to_collapse()
                        else:
                            self._absorb_type_into_match_and_del()
                    else:
                        # then the synth footprint was not already represented by an existing initial type
                        # Check to see if the big type is contains mutiple basal.
                        # If it does then we should extract as above. This should be exactly the
                        # same code as above but extracting into a new type rather than an existing one
                        if self.current_fp_to_collapse.contains_multiple_basal_sequences:
                            new_blank_initial_type = self._create_new_init_type_and_add_to_init_type_list()

                            # Now remove the above new type's worth of info from the current big footprint
                            self.current_fp_to_collapse.substract_init_type_from_other_init_type(
                                new_blank_initial_type)

                            if self.current_fp_to_collapse.set_of_maj_ref_seqs:
                                self._collapse_type_to_matching_init_type_if_exists()
                            else:
                                self._del_type_to_collapse()
                        else:
                            self._create_new_initial_type_from_synth_type_and_del_type_to_collapse()

    def _create_new_init_type_and_add_to_init_type_list(self):
        new_blank_initial_type = InitialType(
            reference_sequence_set=self.current_synth_fp, clade_collection_list=list(
                self.current_fp_to_collapse.clade_collection_list))
        self.parent.initial_types_list.append(new_blank_initial_type)
        return new_blank_initial_type

    def _create_new_initial_type_from_synth_type_and_del_type_to_collapse(self):
        self._create_new_init_type_and_add_to_init_type_list()
        self._del_type_to_collapse()

    def _absorb_type_into_match_and_del(self):
        self.matching_existing_init_type_outer.absorb_large_init_type(self.current_fp_to_collapse)
        self.parent.initial_types_list.remove(self.current_fp_to_collapse)
        self.parent.unsupported_list.remove(self.current_fp_to_collapse)

    def _del_type_to_collapse(self):
        """Then the initial type to collapse no longer contains any of its original maj ref seqs and should be deleted.
        """
        self.parent.initial_types_list.remove(self.current_fp_to_collapse)
        if self.current_fp_to_collapse in self.parent.unsupported_list:
            self.parent.unsupported_list.remove(self.current_fp_to_collapse)

    def _collapse_type_to_matching_init_type_if_exists(self):
        """ If the type to collapse still has ref seqs:
        Check to see if its new profile (i.e. after extraction) matches an initial type that already exists.
        """
        if self._extracted_initial_type_fp_now_matches_existing_initial_type():
            self._absorb_matching_init_type_and_delete()
        self._eval_remove_type_from_unsup_list()

    def _eval_remove_type_from_unsup_list(self):
        """We now need to decide if the footprint to collapse should be removed from the unsupported_list.
        This will depend on if it is longer than n or not.
        """
        if self.current_fp_to_collapse.profile_length < self.parent.current_n:
            self.parent.unsupported_list.remove(self.current_fp_to_collapse)

    def _absorb_matching_init_type_and_delete(self):
        """Here we have found an intial type that exactly matches the initial type to
        collapse's new footprint (i.e. after extraction).
        Now absorb the found match initial type in to the initial type to be collapsed.
        We do it this way around the initial type to collapse stays in the corect place
        i.e. in the unsupported list of not.
        After absorption, remove the matched initial type from the initial types list and unsupported list
        """
        self.current_fp_to_collapse.absorb_large_init_type(
            self.matching_existing_init_type_inner)
        if self.matching_existing_init_type_inner in self.parent.unsupported_list:
            self.parent.unsupported_list.remove(self.matching_existing_init_type_inner)
        self.parent.initial_types_list.remove(self.matching_existing_init_type_inner)

    def _extracted_initial_type_fp_now_matches_existing_initial_type(self):
        """If the initial type to collapse still contains maj
        ref sequences, then it is still a profile that we needs to be assessed for collapse.
        Now need to check if it's new profile (i.e. after extraction) matches that of any of the other initial types."""
        for j in range(len(self.parent.initial_types_list)):
            if self.parent.initial_types_list[j].profile == self.current_fp_to_collapse.profile:
                if self.parent.initial_types_list[j] != self.current_fp_to_collapse:
                    self.matching_existing_init_type_inner = self.parent.initial_types_list[j]
                    return True
        return False

    def _should_extract_rather_than_absorb(self):
        """We have to check whether the init_type to collapse is a multiple basal seqs and therefore
        whether this is an extraction or a absorption bear in mind that it doesn't matter if the matching
        smaller n-1 initial type we are absorbing or extracting into is a multi basal. We will worry about
        that in the next iteration.
        """
        return self.current_fp_to_collapse.contains_multiple_basal_sequences

    def _synth_fp_matches_existing_intial_type(self):
        """If an non-synthetic init_type already exists with the same footprint as the synthetic footprint in
        question then add the init_type to be collapsed to it rather than creating a new initial type from the
        synthetic footprint.
        """
        for i in range(len(self.parent.initial_types_list)):
            # for initT_one in initial_types_list:
            if self.parent.initial_types_list[i].profile == self.current_synth_fp:
                self.matching_existing_init_type_outer = self.parent.initial_types_list[i]
                return True
        return False


    def _validate_init_type_is_unsup_and_synth_fp_is_subset(self):
        """Then this footprint hasn't been collapsed anywhere yet.
        We also check to make sure that the initial type's profile hasn't been modified (i.e. by extraction) and check
        that the synthetic fp is still a sub set of the initial type's profile.
        """
        if self.current_fp_to_collapse in self.parent.unsupported_list:
            if self.current_synth_fp.issubset(self.current_fp_to_collapse.profile):
                return True
        return False


class SyntheticFootprintPermutator:
    def __init__(self, parent_supported_footprint_identifier):
        self.parent = parent_supported_footprint_identifier
        self.len_n_profile_input_mp_list = Queue()
        self.mp_manager = Manager()
        self.collapse_n_mer_mp_dict = self.mp_manager.dict()
        self._populate_mp_input_queue()

    def _populate_mp_input_queue(self):
        for len_n_initial_type in self.parent.list_of_initial_types_of_len_n:
            self.len_n_profile_input_mp_list.put(len_n_initial_type.profile)
        for n in range(self.parent.parent.parent.args.num_proc):
            self.len_n_profile_input_mp_list.put('STOP')

    def permute_synthetic_footprints(self):
        all_processes = []

        for N in range(self.parent.parent.parent.args.num_proc):
            p = Process(target=self.permute_synthetic_footprints_worker, args=())
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

        sys.stdout.write(f'\rGenerated {len(self.collapse_n_mer_mp_dict)} synthetic footprints')

    def permute_synthetic_footprints_worker(self):
        for footprint_set in iter(self.len_n_profile_input_mp_list.get, 'STOP'):
            temp_dict = {
                frozenset(tup): [] for tup in itertools.combinations(footprint_set, self.parent.current_n - 1)}
            self.collapse_n_mer_mp_dict.update(temp_dict)
            sys.stdout.write(f'\rGenerated iterCombos using {current_process().name}')


class FootprintCollapser:
    """Responsible for collapsing the long initial type into a short initial type."""
    def __init__(self, parent_supported_footprint_identifier, footprint_to_collapse_index):
        self.parent = parent_supported_footprint_identifier
        self.fp_index = footprint_to_collapse_index
        self.long_initial_type = self.parent.large_fp_to_collapse_list[self.fp_index]
        self.short_initial_type = self.parent.collapse_dict[self.long_initial_type]
        # bool on whether we are simply collapsing big into small or whether we need to
        # extract certain sequences of the big footprint to go into the small (i.e. if the long
        # fp contains multiple basal seqs
        self.should_extract_not_delete = self.long_initial_type.contains_multiple_basal_sequences
        # whether an extracted
        self.match = False

    def collapse_footprint(self):
        if not self.should_extract_not_delete:
            self._collapse_long_into_short_no_extraction()
        else:
            self._collapse_long_into_short_by_extraction()
            self.match = False
            for p in range(len(self.parent.initial_types_list)):
                intial_type_being_checked = self.parent.initial_types_list[p]
                if self._extracted_long_initial_type_now_has_same_fp_as_another_initial_type(intial_type_being_checked):
                    self.match = True
                    if self._matching_initial_type_in_initial_types_still_to_be_collapsed(intial_type_being_checked, p):
                        self._absorb_long_initial_type_into_matching_intial_type(intial_type_being_checked)
                        break
                    else:
                        self._absorb_matching_initial_type_into_long_intial_type(intial_type_being_checked)

                        if self.long_initial_type.profile_length < self.parent.current_n:
                            # If the left over type is less than n then we need to now remove it from the un
                            # supported list as it will be collapsed on another iteration than this one.
                            self.parent.unsupported_list.remove(self.long_initial_type)
                        else:
                            if self._long_initial_type_has_sufficient_support():
                                self.parent.unsupported_list.remove(self.long_initial_type)
                            else:
                                # If insufficient support then leave it in the unsupportedlist
                                # and it will go on to be seen if it can be collapsed into one of the insilico
                                # types that are genearted.
                                pass
                        break
            if not self.match:
                if self.long_initial_type.profile_length < self.parent.current_n:
                    self.parent.unsupported_list.remove(self.long_initial_type)

    def _long_initial_type_has_sufficient_support(self):
        if self.long_initial_type.support >= self.parent.required_support:
            if not self.long_initial_type.support.contains_multiple_basal_sequences:
                return True
        return False

    def _absorb_matching_initial_type_into_long_intial_type(self, intial_type_being_checked):
        self.long_initial_type.absorb_large_init_type(intial_type_being_checked)
        if intial_type_being_checked in self.parent.unsupported_list:
            self.parent.unsupported_list.remove(intial_type_being_checked)
        self.parent.initial_types_list.remove(intial_type_being_checked)

    def _absorb_long_initial_type_into_matching_intial_type(self, intial_type_being_checked):
        """Collapse the long initial type into the matching initial type.
        Because this collapsing can cause the matching type that is also in the collapse
        dict to gain support we will also check to if each of the types in the collapse dict
        have gained sufficient support to no longer need collapsing.
        We will do this earlier in the process, not here.
        """
        intial_type_being_checked.absorb_large_init_type(self.long_initial_type)
        self.parent.unsupported_list.remove(self.long_initial_type)
        self.parent.initial_types_list.remove(self.long_initial_type)

    def _matching_initial_type_in_initial_types_still_to_be_collapsed(self, intial_type_being_checked, index):
        return intial_type_being_checked in self.parent.large_fp_to_collapse_list[index + 1:]

    def _extracted_long_initial_type_now_has_same_fp_as_another_initial_type(self, intial_type_being_checked):
        """Check if the new profile created from the original footprintToCollapse that has now had the short
        initial_type extracted from it is already shared with another initial type.
        If so then we need to combine the initial types.
        Else then we don't need to do anything
        If the matching type is also in the collapse dict then we will collapse to that type
        Else if the type is not in the collapse dict then we will absorb that type.
        There should only be a maximum of one initT that has the same footprint.
        """
        if intial_type_being_checked.profile == self.long_initial_type.profile:
            if intial_type_being_checked != self.long_initial_type:
                return True
        return False

    def _collapse_long_into_short_by_extraction(self):
        self.short_initial_type.extract_support_from_large_initial_type(self.long_initial_type)

    def _collapse_long_into_short_no_extraction(self):
        """If long type does not contain multiple basal sequences then we do not need to extract and we
        can collapse the large into the small. We need to simply extend the clade collection list of the short
        type with that of the large. We also need to add the ref_seq lists of the large to the small.
        Finally, remove the large_init_type from the init_type_list and from the unsupported_list
        """
        self.short_initial_type.absorb_large_init_type(self.long_initial_type)
        self.parent.initial_types_list.remove(self.long_initial_type)
        self.parent.unsupported_list.remove(self.long_initial_type)

class InitialType:
    def __init__(self, reference_sequence_set, clade_collection_list, maj_dsss_list=False):
        self.profile = reference_sequence_set
        self.profile_length = len(self.profile)
        self.contains_multiple_basal_sequences, self.basalSequence_list = \
            self.check_if_initial_type_contains_basal_sequences()
        self.clade_collection_list = list(clade_collection_list)
        self.support = len(self.clade_collection_list)
        # We may move away from using the dsss but for the time being we will use it
        if maj_dsss_list:
            self.majority_sequence_list, self.set_of_maj_ref_seqs = self.create_majority_sequence_list_for_initial_type(
                maj_dsss_list)
        else:
            self.majority_sequence_list, self.set_of_maj_ref_seqs = \
                self.create_majority_sequence_list_for_inital_type_from_scratch()

    def __repr__(self):
        return str(self.profile)

    def check_if_initial_type_contains_basal_sequences(self):
        """This function will return two items, firstly a list a bool if there are multiple basal sequences contained
        within the profile_set and secondly it will return a list of the
        I will just check the profile sequence"""
        basal_seq_list = []
        found_c15 = False
        for rs in self.profile:
            if rs.name == 'C3':
                basal_seq_list.append('C3')
            elif rs.name == 'C1':
                basal_seq_list.append('C1')
            elif 'C15' in rs.name and not found_c15:
                basal_seq_list.append('C15')
                found_c15 = True

        if len(basal_seq_list) > 1:
            return True, basal_seq_list
        else:
            return False, basal_seq_list

    def substract_init_type_from_other_init_type(self, other_init_type):
        self.profile = self.profile.difference(other_init_type.profile)
        self.profile_length = len(self.profile)
        self.basalSequence_list = list(set(self.basalSequence_list).difference(set(other_init_type.basalSequence_list)))
        if len(self.basalSequence_list) > 1:
            self.contains_multiple_basal_sequences = True
        else:
            self.contains_multiple_basal_sequences = False
        self.majority_sequence_list, self.set_of_maj_ref_seqs = \
            self.create_majority_sequence_list_for_inital_type_from_scratch()

    def create_majority_sequence_list_for_initial_type(self, maj_dsss_list):
        # I'm trying to remember what form this takes. I think we'll need to be looking
        # This should be a list of lists. There should be a list for each
        # cladeCollection in the self.clade_collection_list
        # Within in each of the lists we should have a list of dataSetSampleSequence objects
        # We should already have a list of the dsss's with one dsss for each
        # of the cladeCollections found in the maj_dsss_list
        # we will look to see if there are multiple basal sequences
        # If there are multiple basal sequences then for each cladeCollection within the intial type we will ad
        # the dsss to the list. If there are not multiple basal sequences then we will simply add the dss to a list
        set_of_majority_reference_sequences = set()
        master_dsss_list = []
        if self.contains_multiple_basal_sequences:
            for clade_collection_obj in self.clade_collection_list:
                temp_dsss_list = []
                data_set_sample_sequence_list = list(
                    DataSetSampleSequence.objects.filter(clade_collection_found_in=clade_collection_obj).order_by(
                        '-abundance'))
                # for each of the basal seqs in the basal seqs list, find the dsss representative
                for basal_seq in self.basalSequence_list:
                    if basal_seq == 'C15':
                        # Then we just need to find the most abundnant dsss that's name contains the C15

                        for dsss in data_set_sample_sequence_list:
                            if 'C15' in dsss.reference_sequence_of.name:
                                temp_dsss_list.append(dsss)
                                # Important to break so that we only add the first and most abundant C15 seq
                                break
                    else:
                        # then we are looking for exact matches
                        for dsss in data_set_sample_sequence_list:
                            if dsss.reference_sequence_of.name == basal_seq:
                                temp_dsss_list.append(dsss)
                                break
                # We should also make sure that the original maj sequence is found in the list
                if data_set_sample_sequence_list[0] not in temp_dsss_list:
                    temp_dsss_list.append(data_set_sample_sequence_list[0])
                # Make sure that each of the refSeqs for each of the basal or majs are in the maj set
                for dss in temp_dsss_list:
                    set_of_majority_reference_sequences.add(dss.reference_sequence_of)
                # Here we should have a list of the dsss instances that represent the basal
                # sequences for the clade_collection_object in Q
                master_dsss_list.append(temp_dsss_list)
        else:
            # Then there is ony one basal sequence in this initial type and so we simply need to surround the maj
            # with a list.
            for i in range(len(self.clade_collection_list)):
                master_dsss_list.append(maj_dsss_list[i])
                set_of_majority_reference_sequences.add(maj_dsss_list[i][0].reference_sequence_of)

        return master_dsss_list, set_of_majority_reference_sequences

    def create_majority_sequence_list_for_inital_type_from_scratch(self):
        # This will be like above but will not start with the maj_dsss_list
        # we will go through each of the cladeCollections of the type and get the maj sequence for the type

        # if the init type has multiple basal sequences then we will have to find the actual maj and the basal seq dsss
        set_of_majority_reference_sequences = set()
        master_dsss_list = []
        if self.contains_multiple_basal_sequences:
            for clade_collection_obj in self.clade_collection_list:
                temp_dsss_list = []
                dsss_in_cc = list(
                    DataSetSampleSequence.objects.filter(clade_collection_found_in=clade_collection_obj).order_by(
                        '-abundance'))
                dsss_in_cc_in_profile = [dsss for dsss in dsss_in_cc if dsss.reference_sequence_of in self.profile]
                # first find the dsss that are the representatives of the basal types
                for basal_seq in self.basalSequence_list:
                    if basal_seq == 'C15':
                        # Then we just need to find the most abundnant dsss that's name contains the C15

                        for dsss in dsss_in_cc_in_profile:
                            if 'C15' in dsss.reference_sequence_of.name:
                                temp_dsss_list.append(dsss)
                                # Important to break so that we only add the first and most abundant C15 seq
                                break
                    else:
                        # then we are looking for exact matches
                        for dsss in dsss_in_cc_in_profile:
                            if dsss.reference_sequence_of.name == basal_seq:
                                temp_dsss_list.append(dsss)
                                break
                # now add the actual maj dsss if not one of the basal seqs
                basal_dsss = dsss_in_cc_in_profile[0]
                if basal_dsss not in temp_dsss_list:
                    # Then the actual maj is not already in the list
                    temp_dsss_list.append(basal_dsss)
                for dsss in temp_dsss_list:
                    set_of_majority_reference_sequences.add(dsss.reference_sequence_of)
                master_dsss_list.append(temp_dsss_list)

        # else we are just going to be looking for the actual maj dsss
        else:
            for clade_collection_obj in self.clade_collection_list:
                temp_dsss_list = []
                dsss_in_cc = list(
                    DataSetSampleSequence.objects.filter(clade_collection_found_in=clade_collection_obj).order_by(
                        '-abundance'))
                dsss_in_cc_in_profile = [dsss for dsss in dsss_in_cc if dsss.reference_sequence_of in self.profile]
                basal_dsss = dsss_in_cc_in_profile[0]
                temp_dsss_list.append(basal_dsss)
                set_of_majority_reference_sequences.add(basal_dsss.reference_sequence_of)
                master_dsss_list.append(temp_dsss_list)

        return master_dsss_list, set_of_majority_reference_sequences

    def absorb_large_init_type(self, large_init_type):
        """The aim of this function is simply to add the infomation of the large init type to that of the small init
        type"""
        self.clade_collection_list.extend(large_init_type.clade_collection_list)
        self.majority_sequence_list.extend(large_init_type.majority_sequence_list)
        self.support = len(self.clade_collection_list)
        self.set_of_maj_ref_seqs.update(large_init_type.set_of_maj_ref_seqs)

    def extract_support_from_large_initial_type(self, large_init_type):
        """The aim of this function differs from above. We are extracting support for this small init_type from
        the large_init type. Once we have extracted the support then we will need to reinitialise the bigtype"""

        # 1 - create the list of maj dss lists that will be added to the small init type from the large init type
        # do this by sending over any dss from the big type that is a refseq of the refseqs that the small and
        # large init types have in common
        large_init_type_ref_seqs = large_init_type.set_of_maj_ref_seqs
        small_init_type_ref_seqs = self.set_of_maj_ref_seqs
        ref_seqs_in_common = large_init_type_ref_seqs & small_init_type_ref_seqs
        temp_majdss_list_list = []
        # Keep track of whether some new maj_ref_seqs have been added to the small init type
        # TODO not sure if we need to do this
        new_maj_seq_set = set()

        for i in range(len(large_init_type.majority_sequence_list)):
            list_of_dsss_to_remove_from_large_init_type = []
            temp_dss_list = []
            for j in range(len(large_init_type.majority_sequence_list[i])):
                if large_init_type.majority_sequence_list[i][j].reference_sequence_of in ref_seqs_in_common:
                    # Then this is one of the dsss that we should remove from the clade_collection_object
                    # in the big init type and
                    # add to the small init type
                    temp_dss_list.append(large_init_type.majority_sequence_list[i][j])
                    # NB important to add the reference_sequence_of before removing the dsss from the large type
                    new_maj_seq_set.add(large_init_type.majority_sequence_list[i][j].reference_sequence_of)
                    list_of_dsss_to_remove_from_large_init_type.append(large_init_type.majority_sequence_list[i][j])
            for dsss in list_of_dsss_to_remove_from_large_init_type:
                large_init_type.majority_sequence_list[i].remove(dsss)

            temp_majdss_list_list.append(temp_dss_list)
        # At this point we should have a list of maj dss lists that we can extend the small init type with
        # we have also removed the dss in question from the large init type

        # 2 Now extract into the small init type
        self.clade_collection_list.extend(large_init_type.clade_collection_list)
        self.majority_sequence_list.extend(temp_majdss_list_list)
        self.set_of_maj_ref_seqs.update(new_maj_seq_set)

        # 3 Now modify the large init_type
        # we have already removed the maj seqs from the large init_type
        # need to change, profile, profile length
        # clade_collection_list should not change
        # put through check_if_initial_type_contains_basal_seqs
        # re-initialise set of set_of_maj_ref_seqs

        # new large profile should simply be the ref seqs of the small and large profiles not found in common
        # essentially we extract the small init_types' profile from the large
        large_init_type.profile = large_init_type.profile.difference(self.profile)
        large_init_type.profile_length = len(large_init_type.profile)
        large_init_type.contains_multiple_basal_sequences, large_init_type.basalSequence_list = \
            large_init_type.check_if_initial_type_contains_basal_sequences()
        large_init_type.majority_sequence_list, large_init_type.set_of_maj_ref_seqs = \
            large_init_type.create_majority_sequence_list_for_inital_type_from_scratch()

    def remove_small_init_type_from_large(self, small_init_type):
        prof_reference_sequence = list(small_init_type.profile)[0]
        self.profile = self.profile.difference(small_init_type.profile)
        self.profile_length = len(self.profile)
        self.contains_multiple_basal_sequences, self.basalSequence_list = \
            self.check_if_initial_type_contains_basal_sequences()
        # remove the ref_seq and dsss from the majority_sequence_list and set_of_maj_ref_seqs
        # remove ref_seq from set_of_maj_ref_seqs
        if prof_reference_sequence in self.set_of_maj_ref_seqs:
            self.set_of_maj_ref_seqs.remove(prof_reference_sequence)
        for i in range(len(small_init_type.clade_collection_list)):  # For each list of dsss
            for j in small_init_type.clade_collection_list[i]:  # for each dsss
                if j.reference_sequence_of == prof_reference_sequence:
                    del small_init_type.clade_collection_list[i][j]

        self.majority_sequence_list, self.set_of_maj_ref_seqs = \
            self.create_majority_sequence_list_for_inital_type_from_scratch()

    def __str__(self):
        return str(self.profile)


class FootprintRepresentative:
    def __init__(self, cc, cc_maj_ref_seq):
        self.cc_list = [cc]
        self.maj_seq_list = [cc_maj_ref_seq]

class CCInfoHolder:
    """An object used in the FootprintDictPopWorker to hold information for each of the clade collections
    that will be used to popualte the clade footprint dicts"""
    def __init__(self, footprint, clade_index, cc, maj_ref_seq):
        self.footprint = footprint
        self.clade_index = clade_index
        self.cc = cc
        self.maj_ref_seq = maj_ref_seq

