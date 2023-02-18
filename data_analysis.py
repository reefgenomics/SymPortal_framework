import os
import shutil
import sys
from typing import FrozenSet
from dbApp.models import AnalysisType, ReferenceSequence, CladeCollectionType, CladeCollection
import itertools
from collections import defaultdict
import virtual_objects
import numpy as np
from scipy.stats import gaussian_kde
import symportal_utils
from general import ThreadSafeGeneral
import string
import re
import sp_config
import json
from django import db

class SPDataAnalysis:
    def __init__(self, workflow_manager_parent, data_analysis_obj, force_basal_lineage_separation):
        self.workflow_manager = workflow_manager_parent
        self.force_basal_lineage_separation = force_basal_lineage_separation
        self.temp_wkd = os.path.join(self.workflow_manager.symportal_root_directory, 'temp')
        self._del_and_remake_temp_wkd()
        self.data_analysis_obj = data_analysis_obj
        # The abundance that a given DIV must be found at when it has been considered 'unlocked'
        # https://github.com/didillysquat/SymPortal_framework/wiki/The-SymPortal-logic#type-profile-assignment---logic
        self.unlocked_abundance = 0.0001
        self.clade_list = list('ABCDEFGHI')
        self.ccs_of_analysis = self.data_analysis_obj.get_clade_collections()
        # List that will hold a dictionary for each clade
        # Each dictionary will hold key = footprint (set of sequences)
        # value = [[] []] where []0 = list of cladeCollections containing given footprint
        # and []1 = list of majority sequence for given sample
        self.clade_footp_dicts_list = [{} for _ in self.clade_list]
        self.list_of_initial_types_after_collapse = None
        self.current_clade = None
        self.list_of_data_set_uids = [
                int(ds_id_str) for ds_id_str in self.data_analysis_obj.list_of_data_set_uids.split(',')]
        self.virtual_object_manager = virtual_objects.VirtualObjectManager(
            within_clade_cutoff=self.workflow_manager.within_clade_cutoff,
            num_proc=self.workflow_manager.args.num_proc,
            list_of_data_set_uids=self.list_of_data_set_uids,
            force_basal_lineage_separation=self.force_basal_lineage_separation)
        self.thread_safe_general = ThreadSafeGeneral()

    def analyse_data(self):
        print('\n\nBeginning profile discovery')
        self._populate_clade_fp_dicts_list()

        self._collapse_footprints_and_make_analysis_types()

        self._associate_vats_to_vccs()

        self._check_for_artefacts()
        print('TYPE DISCOVERY COMPLETE')

        self._reset_vcc_vat_rep_abund_dicts()

        self._profile_assignment()

        self._update_grand_tot_attribute_for_vats()

        self._name_divs()

        self._associate_species_designations()

        self._del_and_remake_temp_wkd()

        print('DATA ANALYSIS COMPLETE')
        self._make_analysis_type_objects_from_vats()

    def _update_grand_tot_attribute_for_vats(self):
        """We need to populate the grand_tot_num_instances_of_vat_in_analysis attribute of the vats after
        # the profile assignment."""
        for vat in self.virtual_object_manager.vat_manager.vat_dict.values():
            self.virtual_object_manager.vat_manager.vat_dict[vat.id].grand_tot_num_instances_of_vat_in_analysis = len(
                vat.clade_collection_obj_set_profile_assignment)

    def _reset_vcc_vat_rep_abund_dicts(self):
        for vcc_uid in self.virtual_object_manager.vcc_manager.vcc_dict.keys():
            self.virtual_object_manager.vcc_manager.vcc_dict[
                vcc_uid].analysis_type_obj_to_representative_rel_abund_in_cc_dict = {}

    def _make_analysis_type_objects_from_vats(self):
        print('\nConverting VirtualAnalysisTypes to database AnalysisTypes')
        for vat in self.virtual_object_manager.vat_manager.vat_dict.values():
            sys.stdout.write(f'\r{vat.name}')
            new_at = self._create_analysis_type_from_vat(vat)

            self._update_uid_of_vat(new_at, vat)
        self._update_keys_of_vat_dict()

        # now create the clade_collection_types
        # we will need a cct for each vat, vcc combination
        clade_collection_type_list_for_bulk_create = []
        for vat in self.virtual_object_manager.vat_manager.vat_dict.values():
            for vcc in vat.clade_collection_obj_set_profile_assignment:
                clade_collection_type_list_for_bulk_create.append(
                    CladeCollectionType(
                        analysis_type_of=AnalysisType.objects.get(id=vat.id),
                        clade_collection_found_in=CladeCollection.objects.get(id=vcc.id)))
        for cct_chunk in self.thread_safe_general.chunks(clade_collection_type_list_for_bulk_create):
            CladeCollectionType.objects.bulk_create(cct_chunk)

    def _update_keys_of_vat_dict(self):
        # now update remake the vat dict so that the correct ids are used
        new_dict = {}
        for vat in self.virtual_object_manager.vat_manager.vat_dict.values():
            new_dict[vat.id] = vat
        self.virtual_object_manager.vat_manager.vat_dict = new_dict

    def _update_uid_of_vat(self, new_at, vat):
        # now update the id of the vat
        self.virtual_object_manager.vat_manager.vat_dict[vat.id].id = new_at.id

    def _create_analysis_type_from_vat(self, vat):
        ordered_footprint_list = ','.join(str(rs_id) for rs_id in list(vat.multi_modal_detection_rel_abund_df))
        majority_reference_sequence_set = ','.join([str(rs_id) for rs_id in vat.majority_reference_sequence_uid_set])
        list_of_clade_collections = ','.join(
            [str(cc.id) for cc in vat.clade_collection_obj_set_profile_assignment])
        footprint_sequence_abundances = json.dumps(vat.abs_abund_of_ref_seqs_in_assigned_vccs_df.values.tolist())
        footprint_sequence_ratios = json.dumps(vat.multi_modal_detection_rel_abund_df.values.tolist())
        artefact_intras = ','.join(str(rs_uid) for rs_uid in vat.artefact_ref_seq_uid_set)
        new_at = AnalysisType(
            data_analysis_from=self.data_analysis_obj,
            ordered_footprint_list=ordered_footprint_list,
            majority_reference_sequence_set=majority_reference_sequence_set,
            list_of_clade_collections=list_of_clade_collections,
            footprint_sequence_abundances=footprint_sequence_abundances,
            footprint_sequence_ratios = footprint_sequence_ratios,
            clade=vat.clade, co_dominant=vat.co_dominant, name=vat.name,
            species=vat.species, artefact_intras=artefact_intras
        )
        new_at.save()
        return new_at

    def _del_and_remake_temp_wkd(self):
        if os.path.exists(self.temp_wkd):
            shutil.rmtree(self.temp_wkd)
        os.makedirs(self.temp_wkd, exist_ok=True)

    def _associate_species_designations(self):
        for vat in self.virtual_object_manager.vat_manager.vat_dict.values():
            # For the time being I am disabling the species association
            # We may reimplement this if we can get a separate species description
            # platform up and running
            # I am disabling it in the _associate_species_info_to_vat method of 
            # the SpeciesAssociation class.
            species_association = self.SpeciesAssociation(vat=vat)
            species_association.assign_species()

    class SpeciesAssociation:
        def __init__(self, vat):
            self.vat = vat
            self.assigned_species = []
            self.maj_seq_names = [rs.name for rs in self.vat.majority_reference_sequence_obj_set]
            self.all_seq_names = [rs.name for rs in self.vat.footprint_as_ref_seq_objs_set]

        def assign_species(self):
            """For each analysis type check and assign the associated species"""

            self._designate_species()

            self._associate_species_info_to_vat()

        def _designate_species(self):
            if self.vat.clade == 'A':
                self._clade_a_associations()
            elif self.vat.clade == 'B':
                self._clade_b_associations()
            elif self.vat.clade == 'C':
                self._clade_c_associations()
            elif self.vat.clade == 'D':
                self._clade_d_associations()
            elif self.vat.clade == 'E':
                self._clade_e_associations()
            elif self.vat.clade == 'F':
                self._clade_f_associations()
            elif self.vat.clade == 'G':
                pass
            elif self.vat.clade == 'H':
                pass
            elif self.vat.clade == 'I':
                pass

        def _associate_species_info_to_vat(self):
            # Disable species association
            self.vat.species = 'None'
            # if not self.assigned_species:  # If no suggested species have been associated
            #     self.vat.species = 'None'
            # else:
            #     self.vat.species = ','.join(self.assigned_species)

        def _clade_f_associations(self):
            if 'F1' in self.maj_seq_names:
                self.assigned_species.append('S. kawagutii')

        def _clade_e_associations(self):
            self.assigned_species.append('S. voratum')

        def _clade_d_associations(self):
            # I have decided that we are not going to take into account the abundance of non-maj intragenomic
            # defining sequences. e.g. D6 when calling associated species.
            # This is because there can be a large difference sample to sample in the abundance of the sequences
            # Rather we will assign both clade D species if the required sequences are present
            # We are also giving the researcher the average abundances and SDs for each output type
            if 'D1' in self.maj_seq_names:
                if 'D4' not in self.all_seq_names:
                    self.assigned_species.append('S. glynnii')
                else:  # There is a significant abundance of D4
                    if 'D6' in self.all_seq_names:
                        # Then there is a significant amount of D6
                        self.assigned_species.extend(['S. glynnii', 'S. trenchii'])
                    else:
                        # THere is D1, D4 but not D6
                        self.assigned_species.append('S. trenchii')
            if 'D8' in self.maj_seq_names or 'D12' in self.maj_seq_names or 'D13' in self.maj_seq_names:
                self.assigned_species.append('S. eurythalpos')
            if 'D15' in self.maj_seq_names:
                self.assigned_species.append('S. boreum')

        def _clade_c_associations(self):
            if 'C1' in self.maj_seq_names:
                self.assigned_species.append('S. goreaui')
            if 'C3' in self.all_seq_names and 'C3gulf' in self.all_seq_names:
                self.assigned_species.append('S. thermophilum')

        def _clade_b_associations(self):
            if 'B1' in self.maj_seq_names:
                self.assigned_species.extend(['S. minutum', 'S. antillogorgium', 'S. pseudominutum'])
            if 'B2' in self.maj_seq_names:
                self.assigned_species.append('S. psygmophilum')
            if 'B4' in self.maj_seq_names:
                self.assigned_species.append('S. muscatinei')
            if 'B7' in self.maj_seq_names or 'B13' in self.maj_seq_names:
                self.assigned_species.append('S. endomadracis')
            if 'B2a' in self.maj_seq_names:
                self.assigned_species.append('S. aenigmaticum')

        def _clade_a_associations(self):
            if 'A1' in self.maj_seq_names:
                self.assigned_species.append('S. microadriaticum')
            if 'A2' in self.maj_seq_names:
                self.assigned_species.append('S. pilosum')
            if 'A3' in self.maj_seq_names:
                self.assigned_species.extend(['S. natans', 'S. tridacnidorum'])
            if 'A4' in self.maj_seq_names:
                self.assigned_species.append('S. linucheae')

    def _name_divs(self):
        if sp_config.system_type == 'remote':
            print('Naming unamed DIVs')
            div_namer = self.DIVNamer(parent_sp_data_analysis=self)
            div_namer.name_unamed_div_seqs()
            print('\nDIV naming complete')
        else:
            print('Automatic sequence name generation is currently disabled for local instances of SymPortal.\n'
                  'This is to prevent naming conlifcts between the remote and the '
                  'local instances of SymPortal from arising\n')

    class DIVNamer:
        def __init__(self, parent_sp_data_analysis):
            self.thread_safe_general = ThreadSafeGeneral()
            self.sp_data_analysis = parent_sp_data_analysis
            self.query_fasta_as_list = []
            self.query_fasta_path = os.path.join(
                self.sp_data_analysis.workflow_manager.symportal_root_directory,
                'symbiodiniaceaeDB', 'unnamedRefSeqs.fasta')
            self.db_fasta_as_list = []
            self.db_fasta_path = os.path.join(
                self.sp_data_analysis.workflow_manager.symportal_root_directory,
                'symbiodiniaceaeDB', 'named_seqs_in_SP_remote_db.fa')
            self.blast_output_path = os.path.join(
                self.sp_data_analysis.workflow_manager.symportal_root_directory,
                'symbiodiniaceaeDB', 'blast.out')

            self.blast_analysis_object = symportal_utils.BlastnAnalysis(
                input_file_path=self.query_fasta_path, output_file_path=self.blast_output_path,
                db_path=self.db_fasta_path, output_format_string='6 qseqid sseqid evalue pident qcovs', num_threads=20)
            self.blast_output_dict = None
            self.list_of_sequence_names_that_already_exist = self._set_exist_seq_names()
            self.unamed_div_uid_to_div_obj = {}

        def _set_exist_seq_names(self):
            # This is giving us the strange SSL EOF error: django.db.utils.OperationalError: SSL SYSCALL error: EOF detected
            # I have managed to duplicate this error by filling up the linode server's RAM using python:
            # https://stackoverflow.com/questions/6317818/eat-memory-using-python
            # (I did this from within the test.py script; it didn't work doing manage.py shell)
            # And then run this same query and we get the SSL EOF error. So at least we know what the problem is.
            # The strange thing is, the query itself doesn't seem to use up too much memory so it must
            # be other parts of the script that are using up memory on the linode machine.
            # In support of this when i run a smaller analysis I also don't get the error being raised
            # despite this query being run.
            # So one thing I will try is to reset all of the db connections and see if this helps this command pass.
            db.connections.close_all()
            list_of_sequence_names_that_already_exist = [ref_seq.name for ref_seq in
                                                         ReferenceSequence.objects.filter(has_name=True)]
            list_of_sequence_names_that_already_exist.append('D1a')
            return list_of_sequence_names_that_already_exist

        def name_unamed_div_seqs(self):
            """ Generate names for the DIV ReferenceSequences that currently have no names. This will only happen
            on 'remote' type systems.
            """

            for vat in self.sp_data_analysis.virtual_object_manager.vat_manager.vat_dict.values():
                for unamed_div in [rs for rs in vat.footprint_as_ref_seq_objs_set if not rs.has_name]:
                    self.unamed_div_uid_to_div_obj[unamed_div.id] = unamed_div

            if self.unamed_div_uid_to_div_obj:

                self._create_and_write_query_fasta()

                self._create_and_write_db_fasta()

                self.blast_analysis_object.make_db(title_for_db='name_ref_seqs')

                self.blast_analysis_object.execute_blastn_analysis()

                self.blast_output_dict = {blast_line.split('\t')[0]:blast_line.split('\t')[1:] for blast_line in
                                          self.blast_analysis_object.return_blast_output_as_list()}

                # It will be possible that some of the queries did not return a match
                # as such, we should assert that all did
                assert(set([int(_) for _ in self.blast_output_dict.keys()]) == set(self.unamed_div_uid_to_div_obj.keys()))
                
                self._generate_and_assign_new_names()

                self._regenerate_vat_names()

        def _regenerate_vat_names(self):
            """
            Some of the DIVs were not names and so their ID was being used in the vat name.
            Now that all DIVs have names, use these names. Not that we will need to refreash
            the ReferenceSequence objects that didn't have names now that we have done the naming."""
            print('Regenerating VirtualAnalysisType names')
            for vat in self.sp_data_analysis.virtual_object_manager.vat_manager.vat_dict.values():
                vat.generate_name(at_df=vat.multi_modal_detection_rel_abund_df)
                sys.stdout.write(f'\r{vat.name}')

        def _generate_and_assign_new_names(self):
            # Now assign names to those that aren't exact matches
            # NB this was causing us issues as although we have updated the db object,
            # and we have updated one of th instances of the ref seq we are hold in memory
            # if there were multiple instances of the refseq objects, these other instances
            # will not have been updated.
            # We will update these in the following method
            for no_name_ref_seq_id, output_items in self.blast_output_dict.items():
                ref_seq_in_question = self.unamed_div_uid_to_div_obj[int(no_name_ref_seq_id)]
                if not ref_seq_in_question.has_name:
                    new_name = self._create_new_reference_sequence_name(output_items[0])
                    ref_seq_in_question.name = new_name
                    ref_seq_in_question.has_name = True
                    ref_seq_in_question.save()
                    self.list_of_sequence_names_that_already_exist.append(new_name)

        def _create_new_reference_sequence_name(self, closest_match):
            match_object = re.match("^[A-I]{1}[0-9]{1,3}", closest_match)
            base_name = match_object.group(0)

            #https://stackoverflow.com/questions/23686398/iterate-a-to-zzz-in-python
            for x in range(1, 4):
                for combo in itertools.product(string.ascii_lowercase, repeat=x):
                    alpha = ''.join(combo)
                    if f'{base_name}{alpha}' not in self.list_of_sequence_names_that_already_exist:
                        return f'{base_name}{alpha}'
            return False

        def _create_and_write_db_fasta(self):
            # create the fasta that will be the database to blast against
            for rs in ReferenceSequence.objects.filter(has_name=True):
                self.db_fasta_as_list.extend(['>{}'.format(rs.name), rs.sequence])
            self.thread_safe_general.write_list_to_destination(destination=self.db_fasta_path,
                                              list_to_write=self.db_fasta_as_list)

        def _create_and_write_query_fasta(self):
            # create the fasta as a file that will be queried in the blast
            for rs in self.unamed_div_uid_to_div_obj.values():
                self.query_fasta_as_list.extend([f'>{rs.id}', rs.sequence])
            self.thread_safe_general.write_list_to_destination(
                destination=self.query_fasta_path, list_to_write=self.query_fasta_as_list)

    class ProfileAssigner:
        """Responsible for searching a given VirtualCladeCollection for VirtualAnalysisTypes and associating
        the found VirtualAnalysisTypes to the VirtualCladeCollection"""
        def __init__(self, virtual_clade_collection, parent_sp_data_analysis):
            self.sp_data_analysis = parent_sp_data_analysis
            self.vcc = virtual_clade_collection
            self.vat_match_object_list = []

            # transient objects updated during vat checks
            self.potential_match_object = None

        def assign_profiles(self):
            print(f'\nAssigning ITS2 type profiles to {self.vcc}:')
            list_of_vats_to_search = self._get_list_of_vats_to_search()

            self._find_vats_in_vcc(list_of_vats_to_search)

            # # TODO here we want to make sure that the most abundant sequence of the
            # # VCC is represented by one of the profiles in the self.vat_match_objects_list
            # # If it is not, then we should create a new 1 DIV VAT that is 
            # # the most abundant sequence and assign this to the sample.
            # if not self._maj_seq_is_represented_in_matched_vats():
            #     self._create_vat_of_maj_seq()
            #     self._add_maj_seq_vat_to_matched_list()

            self._associate_vcc_to_vats()



        def _associate_vcc_to_vats(self):
            for vat_match in self.vat_match_object_list:
                print(f'Assigning {vat_match.at.name}')
                vat_match.at.clade_collection_obj_set_profile_assignment.add(vat_match.cc)
                vat_match.cc.analysis_type_obj_to_representative_rel_abund_in_cc_dict[
                    vat_match.at] = vat_match.rel_abund_of_at_in_cc

        def _find_vats_in_vcc(self, list_of_vats_to_search):
            for vat in list_of_vats_to_search:
                if self._search_for_vat_in_vcc(vat=vat):
                    if self._vat_has_divs_in_common_with_other_vats(vat=vat):
                        self._add_new_vat_to_list_if_highest_rel_abund_representative()
                    else:
                        self.vat_match_object_list.append(self.potential_match_object)

        def _get_list_of_vats_to_search(self):
            return [
                vat for vat in
                self.sp_data_analysis.virtual_object_manager.vat_manager.vat_dict.values() if
                vat.ref_seq_uids_set.issubset(self.vcc.footprint_as_frozen_set_of_ref_seq_uids)]

        def _add_new_vat_to_list_if_highest_rel_abund_representative(self):
            """Get a list of the current matches that have refseqs in common with the potential match.
            Iter through this list and compare the represented abundances of the current matches to potential match.
            If the potential match has a lower abundance than any of them, do not accept. If it has a higher abundance
            than all of then, accept and be sure to remove the current matches that shared a div with it from
            the current matches of the VirtualCladeCollection.
            """
            shared_div_match_list = []
            for match_obj in self.vat_match_object_list:
                if self._vats_have_divs_in_common(vat_one=match_obj.at, vat_two=self.potential_match_object.at):
                    shared_div_match_list.append(match_obj)
            # if any one of the current matches reps a greater proportion then do not accept the potential match
            for match_obj in shared_div_match_list:
                if self.potential_match_object.rel_abund_of_at_in_cc < match_obj.rel_abund_of_at_in_cc:
                    return
            # if we reach here then we should delete all of the matches in the shared_div_match list and
            # add the potential new type in their place.
            for match_obj in shared_div_match_list:
                self.vat_match_object_list.remove(match_obj)
            self.vat_match_object_list.append(self.potential_match_object)

        def _vats_have_divs_in_common(self, vat_one, vat_two):
            return vat_one.ref_seq_uids_set.intersection(vat_two.ref_seq_uids_set)

        def _vat_has_divs_in_common_with_other_vats(self, vat):
            if self.vat_match_object_list:
                for vat_match_obj in self.vat_match_object_list:
                    if vat.ref_seq_uids_set.intersection(vat_match_obj.at.ref_seq_uids_set):
                        return True
                return False
            else:
                return False

        def _search_for_vat_in_vcc(self, vat):
            """This will check whether the DIV relative abundance requirements of
            a vat are met by the vcc in question. If met we will keep track of what proportion of the CladeCollection
            this set of refseqs represents. If the VAT is a single sequence VAT then we will require it to
            be present at an abundance of at least 0.05."""
            vcc_rel_abund_dict = self.vcc.ref_seq_id_to_rel_abund_dict

            if len(vat.ref_seq_uids_set) > 1:
                # NB here, when looking to see if the DIVs are found at the right proportions in the CladeCollection
                # to match the VirtualAnalysisType we need to work with the seq abundances as a proportion of
                # only the sequences in the CC that are found in the VAT. However, when we want to get an
                # idea of whether this VAT match is better than another one, then we need to work with the VAT DIV
                # rel abundances as a proportion of all of the sequences in the CladeCollection
                total_abundance_of_seqs_of_vat = sum([rel_abund for ref_seq_uid, rel_abund in vcc_rel_abund_dict.items() if
                                                      ref_seq_uid in vat.ref_seq_uids_set])

                vcc_rel_abund_dict_for_vat = {ref_seq_uid : rel_abund/total_abundance_of_seqs_of_vat for
                                              ref_seq_uid, rel_abund in vcc_rel_abund_dict.items() if
                                              ref_seq_uid in vat.ref_seq_uids_set}

                total_seq_rel_abund_for_cc = []
                for ref_seq_id, ref_seq_req_abund_obj in vat.prof_assignment_required_rel_abund_dict.items():
                    rel_abund_of_div_in_vat_seqs = vcc_rel_abund_dict_for_vat[ref_seq_id]
                    total_seq_rel_abund_for_cc.append(vcc_rel_abund_dict[ref_seq_id])

                    if ref_seq_req_abund_obj.max_abund <= rel_abund_of_div_in_vat_seqs <= ref_seq_req_abund_obj.min_abund:
                        return False

                self.potential_match_object = CCToATMatchInfoHolder(
                    vat=vat, vcc=self.vcc, rel_abund_of_at_in_cc=sum(total_seq_rel_abund_for_cc))
                return True
            else:
                abund_of_vat_in_vcc = vcc_rel_abund_dict[list(vat.ref_seq_uids_set)[0]]
                if abund_of_vat_in_vcc > 0.05:
                    self.potential_match_object = CCToATMatchInfoHolder(
                        vat=vat, vcc=self.vcc, rel_abund_of_at_in_cc=abund_of_vat_in_vcc)
                    return True
                else:
                    return False

    def _profile_assignment(self):
        print('\n\nBeginning profile assignment')
        for virtual_clade_collection in self.virtual_object_manager.vcc_manager.vcc_dict.values():
            if virtual_clade_collection.id == 111118:
                foo = "bar"
            profile_assigner = self.ProfileAssigner(virtual_clade_collection = virtual_clade_collection,
                parent_sp_data_analysis = self)
            profile_assigner.assign_profiles()

        # Reinit the VirtualAnalysisTypes to populate the post-profile assignment objects
        self.reinit_vats_post_profile_assignment()

        self.multimodal_detection()
        print('Profile Assignment Complete')

    def multimodal_detection(self):
        mmd = self.MultiModalDetection(parent_sp_data_analysis=self)
        mmd.run_multimodal_detection()

    class MultiModalDetection:
        def __init__(self, parent_sp_data_analysis):
            self.sp_data_analysis = parent_sp_data_analysis
            self.vat_uids_checked_set = set()
            self.restart = True
            # attributes that will be updated with each vat checked
            self.current_vat = None
            # The two lists that will hold the VirtualCladeCollection objects belonging to each of the potential
            # new VATs resulting from a splitting occurrence.
            self.list_of_vcc_uids_one = []
            self.list_of_vcc_uids_two = []

        def run_multimodal_detection(self):
            print('\nStarting MultiModalDetection')
            while self.restart:
                self.restart = False
                for vat_uid in self.sp_data_analysis.virtual_object_manager.vat_manager.vat_dict.keys():
                    self.current_vat = self.sp_data_analysis.virtual_object_manager.vat_manager.vat_dict[vat_uid]
                    sys.stdout.write(f'\rChecking {self.current_vat.name}')
                    if self.current_vat.id in self.vat_uids_checked_set:
                        continue
                    if len(self.current_vat.ref_seq_uids_set) == 1:
                        self.vat_uids_checked_set.add(self.current_vat.id)
                        continue
                    if len(self.current_vat.clade_collection_obj_set_profile_assignment) < 8:
                        self.vat_uids_checked_set.add(self.current_vat.id)
                        continue
                    for ref_seq_uid_col in list(self.current_vat.multi_modal_detection_rel_abund_df):
                        if not self.restart:
                            self._assess_if_div_multimodal(ref_seq_uid_col)
                            if not self.restart:
                                self.vat_uids_checked_set.add(self.current_vat.id)
                    if self.restart:
                        break

        def _assess_if_div_multimodal(self, ref_seq_uid_col):
            c, modes, pdf, x_grid = self._find_modes_of_abundances(ref_seq_uid_col)
            if modes == 2:
                x_diff_valid = self._assess_if_modes_sufficiently_separated(c, pdf, x_grid)
                if x_diff_valid:
                    self._assign_vccs_to_modes(pdf, ref_seq_uid_col, x_grid)

                    if self._sufficient_support_of_each_mode():
                        self._split_vat_into_two_new_vats()

        def _assess_if_modes_sufficiently_separated(self, c, pdf, x_grid):
            # Must be sufficient separation between the peaks in x axis
            x_diff_valid = False
            if x_grid[c[1]] - x_grid[c[0]] > 0.7:
                x_diff_valid = True
            # plotHists(pdf, x_grid, listOfRatios, listOfTypesToAnalyse[k].name)
            # Must also be sufficient diff between minima y and small peak y
            # This represents the x spread and overlap of the two peaks
            d = list((np.diff(np.sign(np.diff(pdf))) != 0).nonzero()[0] + 1)  # max and min indices
            if min([pdf[d[0]], pdf[d[2]]]) == 0:
                x_diff_valid = False
            else:
                if pdf[d[1]] / min([pdf[d[0]], pdf[d[2]]]) > 0.85:  # Insufficient separation of peaks
                    x_diff_valid = False
            return x_diff_valid

        def _assign_vccs_to_modes(self, pdf, ref_seq_uid_col, x_grid):
            # Then we have found modes that are sufficiently separated.
            self.list_of_vcc_uids_one = []
            self.list_of_vcc_uids_two = []
            min_x = x_grid[list(((np.diff(np.sign(np.diff(pdf))) != 0).nonzero()[0] + 1))[1]]
            for vcc_uid in self.current_vat.multi_modal_detection_rel_abund_df.index.tolist():
                if self.current_vat.multi_modal_detection_rel_abund_df.at[
                    vcc_uid, ref_seq_uid_col] < min_x:
                    self.list_of_vcc_uids_one.append(vcc_uid)
                else:
                    self.list_of_vcc_uids_two.append(vcc_uid)

        def _sufficient_support_of_each_mode(self):
            return len(self.list_of_vcc_uids_one) >= 4 and len(self.list_of_vcc_uids_two) >= 4

        def _update_vccs_rep_abund_dict_for_split_type(self, list_of_vcc_objs, resultant_vat):
            for vcc in list_of_vcc_objs:
                del vcc.analysis_type_obj_to_representative_rel_abund_in_cc_dict[self.current_vat]
                rep_rel_abund_of_resultant_type = sum([vcc.ref_seq_id_to_rel_abund_dict[ref_seq_uid] for ref_seq_uid in resultant_vat.ref_seq_uids_set])
                vcc.analysis_type_obj_to_representative_rel_abund_in_cc_dict[resultant_vat] = rep_rel_abund_of_resultant_type

        def _split_vat_into_two_new_vats(self):
            print(f'\n\nMultiModalDetection: Splitting {self.current_vat.name}')

            list_of_vcc_objs_one = [
                vcc for vcc in self.current_vat.clade_collection_obj_set_profile_assignment if
                vcc.id in self.list_of_vcc_uids_one]
            resultant_type_one = self.sp_data_analysis.virtual_object_manager.vat_manager. \
                make_vat_post_profile_assignment(
                clade_collection_obj_list=list_of_vcc_objs_one,
                ref_seq_obj_list=self.current_vat.footprint_as_ref_seq_objs_set)
            self._update_vccs_rep_abund_dict_for_split_type(
                list_of_vcc_objs=list_of_vcc_objs_one, resultant_vat=resultant_type_one)
            print(f'Created {resultant_type_one.name}')

            list_of_vcc_objs_two = [
                vcc for vcc in self.current_vat.clade_collection_obj_set_profile_assignment if
                vcc.id in self.list_of_vcc_uids_two]
            resultant_type_two = self.sp_data_analysis.virtual_object_manager.vat_manager. \
                make_vat_post_profile_assignment(
                clade_collection_obj_list=list_of_vcc_objs_two,
                ref_seq_obj_list=self.current_vat.footprint_as_ref_seq_objs_set)
            self._update_vccs_rep_abund_dict_for_split_type(
                list_of_vcc_objs=list_of_vcc_objs_two, resultant_vat=resultant_type_two)
            print(f'Created {resultant_type_two.name}')

            print(f'Destroyed {self.current_vat.name}\n')
            self.sp_data_analysis.virtual_object_manager.vat_manager. \
                delete_virtual_analysis_type(self.current_vat)
            self.restart = True

        def _find_modes_of_abundances(self, ref_seq_uid_col):
            rel_abunds_of_ref_seq = self.current_vat.multi_modal_detection_rel_abund_df.loc[
                                    :, ref_seq_uid_col].values.tolist()
            x_grid = np.linspace(min(rel_abunds_of_ref_seq) - 1, max(rel_abunds_of_ref_seq) + 1, 2000)
            kde = gaussian_kde(rel_abunds_of_ref_seq)
            pdf = kde.evaluate(x_grid)
            c = list((np.diff(np.sign(np.diff(pdf))) < 0).nonzero()[0] + 1)
            modes = len(c)
            return c, modes, pdf, x_grid

    def reinit_vats_post_profile_assignment(self):
        print('\nReinstantiating VirtualAnalysisTypes')
        for vat in self.virtual_object_manager.vat_manager.vat_dict.values():
            sys.stdout.write(f'\r{vat.name}')
            self.virtual_object_manager.vat_manager.reinit_vat_post_profile_assignment(
                vat_to_reinit=vat, new_clade_collection_obj_set=vat.clade_collection_obj_set_profile_assignment)
        print('\nReinstantiation complete')


    def _associate_vats_to_vccs(self):
        """Populate the analysis_type_obj_to_representative_rel_abund_in_cc_dict of the VirtualCladeCollection
        using the VirtualAnalysisTypes."""
        print('Populating starting analysis type info to cc info dict')
        clade_collection_to_type_tuple_list = []
        for vat in self.virtual_object_manager.vat_manager.vat_dict.values():
            initial_clade_collections = vat.clade_collection_obj_set_profile_discovery
            for cc in initial_clade_collections:
                clade_collection_to_type_tuple_list.append((cc, vat))

        for cc, vat in clade_collection_to_type_tuple_list:
            virtual_cc = self.virtual_object_manager.vcc_manager.vcc_dict[cc.id]
            current_type_seq_rel_abund_for_cc = []
            cc_ref_seq_abundance_dict = virtual_cc.ref_seq_id_to_rel_abund_dict
            for ref_seq in vat.footprint_as_ref_seq_objs_set:
                rel_abund = cc_ref_seq_abundance_dict[ref_seq.id]
                current_type_seq_rel_abund_for_cc.append(rel_abund)
            current_type_seq_tot_rel_abund_for_cc = sum(current_type_seq_rel_abund_for_cc)
            virtual_cc.analysis_type_obj_to_representative_rel_abund_in_cc_dict[
                vat] = current_type_seq_tot_rel_abund_for_cc
            sys.stdout.write(f'\rCladeCollection:{cc} AnalysisType:{vat}')

    def _check_for_artefacts(self):
        artefact_assessor = ArtefactAssessor(parent_sp_data_analysis=self)
        artefact_assessor.assess_within_clade_cutoff_artefacts()
        artefact_assessor.reassess_support_of_artefact_div_containing_types()

    def _collapse_footprints_and_make_analysis_types(self):
        for i, clade_fp_dict in enumerate(self.clade_footp_dicts_list):
            self.current_clade = self.clade_list[i]
            if self._there_are_footprints_of_this_clade(clade_fp_dict):
                sfi = SupportedFootPrintIdentifier(clade_footprint_dict=clade_fp_dict, parent_sp_data_analysis=self)
                self.list_of_initial_types_after_collapse = sfi.identify_supported_footprints()
                analysis_type_creator = AnalysisTypeCreator(parent_sp_data_analysis=self)
                analysis_type_creator.create_analysis_types()
        self._verify_all_ccs_associated_to_analysis_type()


    def _verify_all_ccs_associated_to_analysis_type(self):
        print('\nVerifying all CladeCollections have been associated to an AnalysisType...')
        clade_collections_represented_by_types = set()
        for vat in self.virtual_object_manager.vat_manager.vat_dict.values():
            clade_collections_represented_by_types.update([vcc.id for vcc in vat.clade_collection_obj_set_profile_discovery])
        ccs_of_data_analysis_dict = {cc.id : cc for cc in self.ccs_of_analysis}
        if set(ccs_of_data_analysis_dict.keys()).issuperset(clade_collections_represented_by_types):
            set_of_unassociated_cc_uids = set(ccs_of_data_analysis_dict.keys()).difference(clade_collections_represented_by_types)
            diff_list = []
            for diff_uid in set_of_unassociated_cc_uids:
                diff_list.append(ccs_of_data_analysis_dict[diff_uid])
            if len(set_of_unassociated_cc_uids) == 0:
                print('All CladeCollections successfuly associated to at least one AnalysisType')
            else:
                raise RuntimeError(
                    f'{len(set_of_unassociated_cc_uids)} CladeCollections are unassociated from an AnalysisType')


    def _there_are_footprints_of_this_clade(self, clade_fp_dict):
        return clade_fp_dict

    def _populate_clade_fp_dicts_list(self):
        for cc_id, vcc in self.virtual_object_manager.vcc_manager.vcc_dict.items():
            clade_index = self.clade_list.index(vcc.clade)
            if vcc.above_cutoff_ref_seqs_obj_set in self.clade_footp_dicts_list[clade_index]:
                self.clade_footp_dicts_list[clade_index][vcc.above_cutoff_ref_seqs_obj_set].cc_list.append(vcc)
                self.clade_footp_dicts_list[clade_index][vcc.above_cutoff_ref_seqs_obj_set].maj_dss_seq_list.append(vcc.ordered_dsss_objs[0])
            else:
                self.clade_footp_dicts_list[clade_index][vcc.above_cutoff_ref_seqs_obj_set] = FootprintRepresentative(
                cc=vcc, maj_dss_seq_list=vcc.ordered_dsss_objs[0])


class ArtefactAssessor:
    def __init__(self, parent_sp_data_analysis):
        self.sp_data_analysis = parent_sp_data_analysis
        self.list_of_vccs = list(self.sp_data_analysis.virtual_object_manager.\
            vcc_manager.vcc_dict.values())
        # key:VirtualAnalysisType.id, value:VirtualAnalysisType
        self.virtual_analysis_type_dict = self.sp_data_analysis.virtual_object_manager.vat_manager.vat_dict
        self.analysis_types_of_analysis = list(self.virtual_analysis_type_dict.values())
        self.set_of_clades_from_analysis = self._set_set_of_clades_from_analysis()

        # key = set of ref_seq_objects, value = VirtualAnalysisType
        self.ref_seq_fp_set_to_analysis_type_obj_dict = self._init_fp_to_at_dict()
        # Attributes updated on an iterative basis
        self.current_clade = None
        # NB we have the two lists below as we only want to check combinations of the original AnalysiTypes and
        # not the new AnalysisTypes that will be created as part of this process. This is to prevent any infinite
        # loops occuring.
        # A query that will be coninually updated
        self.vat_uids_of_clade_dynamic = None
        # A fixed list of the original types that we stated with
        self.vat_uids_of_clade_static = None
        # A list that holds a tuple of ids that have already been compared.
        self.already_compared_analysis_type_uid_set = set()
        # Bool whether the pair comparisons need to be restarted.
        # This will be true when we have modified a type in anyway
        self.restart_pair_comparisons = True

        # reassess support of artefact DIV containing analysis types attributes
        self.already_checked_vat_uid_list = []

    def _init_fp_to_at_dict(self):
        ref_seq_fp_set_to_analysis_type_obj_dict = {}
        for vat in self.virtual_analysis_type_dict.values():
            ref_seq_fp_set_to_analysis_type_obj_dict[
                frozenset(vat.footprint_as_ref_seq_objs_set)] = vat
        return ref_seq_fp_set_to_analysis_type_obj_dict

    def _types_should_be_checked(self, vat_a, vat_b):
        """Check
        1 - the non-artefact_ref_seqs match
        2 - neither of the types is a subset of the other
        3 - there are artefact ref seqs in at least one of the types"""
        if vat_a.basal_seq == vat_b.basal_seq:
            if vat_a.non_artefact_ref_seq_uid_set and vat_b.non_artefact_ref_seq_uid_set:
                if vat_a.non_artefact_ref_seq_uid_set == vat_b.non_artefact_ref_seq_uid_set:
                    if not set(vat_a.ref_seq_uids_set).issubset(vat_b.ref_seq_uids_set):
                        if not set(vat_b.ref_seq_uids_set).issubset(vat_a.ref_seq_uids_set):
                            if vat_a.artefact_ref_seq_uid_set.union(vat_b.artefact_ref_seq_uid_set):
                                return True
        return False

    def assess_within_clade_cutoff_artefacts(self):
        """Check through all of the types to see if there are super types that have not been identified due to the
        withincladecutoff. Please see:
        https://github.com/didillysquat/SymPortal_framework/wiki/
        The-SymPortal-logic#artefacts-during-the-its2-type-profile-discovery-phase
        For further details"""
        for clade in self.set_of_clades_from_analysis:
            self.current_clade = clade
            self._set_static_and_dynamic_vat_lists()
            self.restart_pair_comparisons = True
            while self.restart_pair_comparisons:
                self.restart_pair_comparisons = False
                self._set_vat_to_compare_pairwise()
                for analysis_type_a_uid, analysis_type_b_uid in itertools.combinations(self.vat_uids_to_check_of_clade, 2):
                    if {analysis_type_a_uid, analysis_type_b_uid} not in self.already_compared_analysis_type_uid_set:
                        vat_a = self.virtual_analysis_type_dict[analysis_type_a_uid]
                        vat_b = self.virtual_analysis_type_dict[analysis_type_b_uid]
                        if self._types_should_be_checked(vat_a, vat_b):
                            print(f'\n\nChecking {vat_a.name} and {vat_b.name} for additional artefactual profiles')

                            ctph = CheckTypePairingHandler(parent_artefact_assessor=self, vat_a=vat_a, vat_b=vat_b)
                            if ctph.check_type_pairing():
                                self.restart_pair_comparisons = True
                                self._reset_dynamic_vat_list()
                                break
                            else:
                                self._log_completed_comparison(analysis_type_a_uid, analysis_type_b_uid)
                        else:
                            self._log_completed_comparison(analysis_type_a_uid, analysis_type_b_uid)

    def reassess_support_of_artefact_div_containing_types(self):
        """Check to see how the association of VirtualCladeCollections to VirtualAnalysisTypes changes when taking into
        account that some sequencs are now 'artefact' sequences and thus 'unlocked'. Please see:
        https://github.com/didillysquat/SymPortal_framework/wiki
        /The-SymPortal-logic#artefacts-during-the-its2-type-profile-assignment-phase
        for further details."""
        for clade in self.set_of_clades_from_analysis:
            self.current_clade = clade
            self.already_checked_vat_uid_list = []
            self.restart_pair_comparisons = True

            while self.restart_pair_comparisons:
                self._set_artefact_div_vat_to_check()
                self.restart_pair_comparisons = False

                for vat_uid in self.vat_uids_to_check_of_clade:
                    if vat_uid not in self.already_checked_vat_uid_list:

                        vat_to_check = self.virtual_analysis_type_dict[vat_uid]
                        print(f'\n\nChecking associations of VirtualCladeCollections to {vat_to_check.name}')
                        cadivvata = CheckArtefactDIVVATAssociations(parent_artefact_assessor=self, vat_to_check=vat_to_check)
                        if cadivvata.check_artefact_div_vat_associations():
                            self.restart_pair_comparisons = True
                            self.already_checked_vat_uid_list.append(vat_to_check.id)
                            break
                        else:
                            self.already_checked_vat_uid_list.append(vat_to_check.id)

    def _set_artefact_div_vat_to_check(self):
        self.vat_uids_to_check_of_clade = [
            vat.id for vat in self.virtual_analysis_type_dict.values() if
            vat.clade == self.current_clade if
            vat.artefact_ref_seq_uid_set]

    def _reset_dynamic_vat_list(self):
        self.vat_uids_of_clade_dynamic = [
            at_id for at_id in self.virtual_analysis_type_dict.keys() if
            self.virtual_analysis_type_dict[at_id].clade == self.current_clade]

    def _log_completed_comparison(self, analysis_type_a_uid, analysis_type_b_uid):
        self.already_compared_analysis_type_uid_set.add(
            frozenset({analysis_type_a_uid, analysis_type_b_uid}))

    def _set_vat_to_compare_pairwise(self):
        self.vat_uids_to_check_of_clade = [at_id for at_id in self.vat_uids_of_clade_dynamic if at_id in
                                     self.vat_uids_of_clade_static]


    def _set_static_and_dynamic_vat_lists(self):
        self.vat_uids_of_clade_dynamic = [
            vat_id for vat_id in self.virtual_analysis_type_dict.keys() if
            self.virtual_analysis_type_dict[vat_id].clade == self.current_clade]
        self.vat_uids_of_clade_static = [
            at_id for at_id in self.virtual_analysis_type_dict.keys() if
            self.virtual_analysis_type_dict[at_id].clade == self.current_clade]

    def _set_set_of_clades_from_analysis(self):
        self.set_of_clades_from_analysis = set()
        for at in self.virtual_analysis_type_dict.values():
            self.set_of_clades_from_analysis.add(at.clade)
        return self.set_of_clades_from_analysis

class PotentialNewType:
    def __init__(
            self, artefact_ref_seq_uid_set, non_artefact_ref_seq_uid_set,
            ref_seq_uids_set, list_of_ref_seq_names, resf_seq_obj_set, force_basal_lineage_separation):
        self.artefact_ref_seq_uid_set = artefact_ref_seq_uid_set
        self.non_artefact_ref_seq_uid_set = non_artefact_ref_seq_uid_set
        self.ref_seq_uids_set = ref_seq_uids_set
        self.name = ','.join(list_of_ref_seq_names)
        if force_basal_lineage_separation:
            self.basal_seq = self._set_basal_seq(list_of_ref_seq_names)
        else:
            self.basal_seq = None
        self.ref_seq_objects_set = resf_seq_obj_set

    def _set_basal_seq(self, list_of_ref_seq_names):
        basal_set = set()
        found_c15_a = False
        for name in list_of_ref_seq_names:
            if name == 'C3':
                basal_set.add('C3')
            elif name == 'C1':
                basal_set.add('C1')
            elif 'C15' in name and not found_c15_a:
                basal_set.add('C15')
                found_c15_a = True

        if len(basal_set) == 1:
            return list(basal_set)[0]
        elif len(basal_set) > 1:
            raise RuntimeError(f'basal seq set {basal_set} contains more than one ref seq')
        else:
            return None


class CheckVCCToVATAssociations:
    """Base class for use by the two ArtefactAssessor instances"""
    def __init__(self, parent_artefact_assessor):
        self.artefact_assessor = parent_artefact_assessor
        self.list_of_vcc_objs_to_check = None
        self.list_of_loss_of_support_info_holder_objs = []
        self.at_obj_to_cc_obj_list_to_be_removed = defaultdict(list)
        self.stranded_ccs = []
        self.ref_seqs_in_common_for_stranded_ccs = set()
        self.at_matching_stranded_ccs = None
        self.new_analysis_type_from_stranded_ccs = None

    def _assess_support_of_pnt_or_vat(self, pnt_or_vat):
        for vcc in self.list_of_vcc_objs_to_check:
            cpntsw = self.CheckPNTorVATSupportWorker(
                virtual_clade_collection_object=vcc, parent_check_type_pairing=self, pnt_or_vat=pnt_or_vat)
            cpntsw.check_pnt_support()

    def _pnt_or_vat_has_support(self):
        return len(self.list_of_loss_of_support_info_holder_objs) >= 4

    def _update_cc_info_for_ccs_that_support_new_type(self, vat_to_add_to_vcc):
        for loss_of_support_info_obj in self.list_of_loss_of_support_info_holder_objs:
            self._remove_no_longer_supported_type_from_cc_info(loss_of_support_info_obj)
            self.add_new_type_to_cc_from_match_obj(loss_of_support_info_obj, vat_to_add_to_vcc)
            self._populate_at_obj_to_cc_obj_to_be_removed_dict(loss_of_support_info_obj)

    def _remove_no_longer_supported_type_from_cc_info(self, loss_of_support_info_obj):
        print(f'removing association of {loss_of_support_info_obj.cc} from {loss_of_support_info_obj.at}')
        del loss_of_support_info_obj.cc.analysis_type_obj_to_representative_rel_abund_in_cc_dict[
            loss_of_support_info_obj.at]

    def add_new_type_to_cc_from_match_obj(self, match_info_obj, vat_to_add_to_vcc):
        print(f'associating {match_info_obj.cc} to {vat_to_add_to_vcc.name}')
        match_info_obj.cc.analysis_type_obj_to_representative_rel_abund_in_cc_dict[
            vat_to_add_to_vcc] = match_info_obj.rel_abund_of_at_in_cc

    def _populate_at_obj_to_cc_obj_to_be_removed_dict(self, loss_of_support_info_obj):
        self.at_obj_to_cc_obj_list_to_be_removed[
            loss_of_support_info_obj.at].append(loss_of_support_info_obj.cc)

    def _reinit_or_del_affected_types_and_create_stranded_cc_list(self):
        for vat, cc_obj_list_val in self.at_obj_to_cc_obj_list_to_be_removed.items():
            new_list_of_ccs_to_associate_to = [cc for cc in vat.clade_collection_obj_set_profile_discovery if cc not in cc_obj_list_val]
            # If the analysis type still has support then simply reinitialize it
            if self._afftected_type_still_has_sufficient_support(new_list_of_ccs_to_associate_to):
                print(f'\nType {vat.name} supported by {len(new_list_of_ccs_to_associate_to)} CCs. Reinitiating.')
                self.artefact_assessor.sp_data_analysis.virtual_object_manager.vat_manager.\
                    reinit_vat_pre_profile_assignment(
                    vat_to_reinit=vat, new_clade_collection_obj_set=new_list_of_ccs_to_associate_to)
            else:
                try:
                    self._del_affected_type_and_populate_stranded_cc_list(
                        vat=vat, new_list_of_ccs_to_associate_to=new_list_of_ccs_to_associate_to)
                except:
                    apples = 'asdf'


    def _afftected_type_still_has_sufficient_support(self, new_list_of_ccs_to_associate_to):
        return len(new_list_of_ccs_to_associate_to) >= 4

    def _del_affected_type_and_populate_stranded_cc_list(self, vat, new_list_of_ccs_to_associate_to):
        print(
            f'Type {vat.name} no longer supported. '
            f'Deleting. {len(new_list_of_ccs_to_associate_to)} CCs stranded.')
        del self.artefact_assessor.ref_seq_fp_set_to_analysis_type_obj_dict[
            frozenset(vat.footprint_as_ref_seq_objs_set)]
        del self.artefact_assessor.virtual_analysis_type_dict[vat.id]
        for vcc in new_list_of_ccs_to_associate_to:
            if self._vcc_supports_only_one_vat(vcc):
                self.stranded_ccs.append(vcc)
            else:
                # vcc was supported by more than one item. Simply remove the at in question from the dict
                # and leave as is.
                del vcc.analysis_type_obj_to_representative_rel_abund_in_cc_dict[vat]

    def _vcc_supports_only_one_vat(self, vcc):
        return len(vcc.analysis_type_obj_to_representative_rel_abund_in_cc_dict.items()) == 1

    def _reassociate_stranded_ccs_if_necessary(self):
        if self._sufficient_stranded_ccs_for_new_analysis_type():
            print(f'{len(self.stranded_ccs)} VirtualCladeCollections are stranded. Rehoming...')
            # Get ref_seqs in common
            self._get_ref_seqs_in_common_btw_stranded_ccs()
            if self.ref_seqs_in_common_for_stranded_ccs:
                if self._analysis_type_already_exists_with_profile_of_seqs_in_common():
                    self._add_stranded_ccs_to_existing_at_and_update_dicts()
                else:
                    if self.artefact_assessor.sp_data_analysis.force_basal_lineage_separation:
                        if not self._ref_seqs_in_common_contain_multiple_basal_seqs():
                            self._add_stranded_ccs_to_new_at_made_from_common_ref_seqs_and_update_dicts()
                        else:
                            self._rehome_cc_individually()
                    else:
                        self._rehome_cc_individually()
            else:
                self._rehome_cc_individually()
        else:

            if self.stranded_ccs:
                print(f'{len(self.stranded_ccs)} VirtualCladeCollections are stranded. Rehoming...')
                self._rehome_cc_individually()
            else:
                print(f'There are no stranded VirtualCladeCollections.')

    def _sufficient_stranded_ccs_for_new_analysis_type(self):
        return len(self.stranded_ccs) >= 4

    def _get_ref_seqs_in_common_btw_stranded_ccs(self):
        print('Finding ReferenceSequences in common between the VirtualCladeCollections')
        list_of_sets_of_ref_seqs_above_cutoff = [cc.above_cutoff_ref_seqs_obj_set for cc in self.stranded_ccs]
        self.ref_seqs_in_common_for_stranded_ccs = list_of_sets_of_ref_seqs_above_cutoff[0].intersection(
            *list_of_sets_of_ref_seqs_above_cutoff[1:])
        print(f'{len(self.ref_seqs_in_common_for_stranded_ccs)} ReferenceSequences found in common.')

    def _analysis_type_already_exists_with_profile_of_seqs_in_common(self):
        try:
            self.at_matching_stranded_ccs = self.artefact_assessor.ref_seq_fp_set_to_analysis_type_obj_dict[
                frozenset(self.ref_seqs_in_common_for_stranded_ccs)]
            return True
        except KeyError:
            return False

    def _add_stranded_ccs_to_existing_at_and_update_dicts(self):
        self._reinit_existing_type_with_additional_ccs()
        self._add_exisiting_type_to_stranded_cc_info_objects(vat=self.at_matching_stranded_ccs)

    def _reinit_existing_type_with_additional_ccs(self):
        self.artefact_assessor.sp_data_analysis.virtual_object_manager.vat_manager.\
            add_ccs_and_reinit_virtual_analysis_type(
            vat_to_add_ccs_to=self.at_matching_stranded_ccs,
            list_of_clade_collection_objs_to_add=self.stranded_ccs)

    def _add_exisiting_type_to_stranded_cc_info_objects(self, vat):
        for cc in self.stranded_ccs:
            cc.analysis_type_obj_to_representative_rel_abund_in_cc_dict = {}
            self.add_a_type_to_cc_without_match_obj(cc=cc, vat=vat)

    def add_a_type_to_cc_without_match_obj(self, cc, vat):
        current_type_seq_rel_abund_for_cc = []
        cc_ref_seq_abundance_dict = cc.ref_seq_id_to_rel_abund_dict
        for ref_seq in vat.footprint_as_ref_seq_objs_set:
            rel_abund = cc_ref_seq_abundance_dict[ref_seq.id]
            current_type_seq_rel_abund_for_cc.append(rel_abund)
        current_type_seq_tot_rel_abund_for_cc = sum(current_type_seq_rel_abund_for_cc)
        cc.analysis_type_obj_to_representative_rel_abund_in_cc_dict[
            vat] = current_type_seq_tot_rel_abund_for_cc

    def _ref_seqs_in_common_contain_multiple_basal_seqs(self):
        """Return False if there is only one basal seq in the profile"""
        basal_seq_list = []
        found_c15 = False
        for rs in self.ref_seqs_in_common_for_stranded_ccs:
            if rs.name == 'C3':
                basal_seq_list.append('C3')
            elif rs.name == 'C1':
                basal_seq_list.append('C1')
            elif 'C15' in rs.name and not found_c15:
                basal_seq_list.append('C15')
                found_c15 = True

        if len(basal_seq_list) > 1:
            return True
        else:
            return False

    def _add_stranded_ccs_to_new_at_made_from_common_ref_seqs_and_update_dicts(self):
        self._make_new_analysis_type_from_stranded_ccs()
        self._add_exisiting_type_to_stranded_cc_info_objects(vat=self.new_analysis_type_from_stranded_ccs)
        self._update_fp_to_at_dict(self.new_analysis_type_from_stranded_ccs)

    def _make_new_analysis_type_from_stranded_ccs(self):
        self.new_analysis_type_from_stranded_ccs = self.artefact_assessor.sp_data_analysis.virtual_object_manager.\
            vat_manager.make_vat_pre_profile_assignment(
            clade_collection_obj_list=self.stranded_ccs, ref_seq_obj_list=self.ref_seqs_in_common_for_stranded_ccs)

    def _update_fp_to_at_dict(self, vat):
        self.artefact_assessor.ref_seq_fp_set_to_analysis_type_obj_dict[
            frozenset(vat.footprint_as_ref_seq_objs_set)] = vat

    def _rehome_cc_individually(self):
        sccrh = self.StrandedCCRehomer(parent_check_type_pairing_handler=self)
        sccrh.rehome_stranded_ccs()

    class CheckPNTorVATSupportWorker:
        def __init__(self, parent_check_type_pairing, virtual_clade_collection_object, pnt_or_vat):
            self.check_type_pairing = parent_check_type_pairing
            self.vcc = virtual_clade_collection_object
            self.pnt_or_vat = pnt_or_vat
            self.pnt_seq_rel_abund_total_for_cc = None
            self.rel_abund_of_current_virtual_analysis_type_of_cc = None
            self.current_virtual_analysis_type_of_cc = None
            self.within_clade_cutoff = self.check_type_pairing.artefact_assessor. \
                sp_data_analysis.workflow_manager.within_clade_cutoff
            self.unlocked_abundance = self.check_type_pairing.artefact_assessor.sp_data_analysis.unlocked_abundance

        def check_pnt_support(self):
            sys.stdout.write(f'\rChecking {self.vcc}')

            if not self._pnt_or_vat_abundances_met():
                return

            if self._vcc_has_analysis_types_associated_to_it():
                self._get_rel_abund_represented_by_current_at_of_cc()

                # if self._another_type_with_divs_in_common_exists():
                #     return

                if self.pnt_seq_rel_abund_total_for_cc > self.rel_abund_of_current_virtual_analysis_type_of_cc:
                    self.check_type_pairing.list_of_loss_of_support_info_holder_objs.append(
                        CCToATMatchInfoHolder(
                            vcc=self.vcc, vat=self.current_virtual_analysis_type_of_cc,
                            rel_abund_of_at_in_cc=self.rel_abund_of_current_virtual_analysis_type_of_cc))
                else:
                    # if cc doesn't support pnt then nothing to do
                    pass

            else:
                # This could be buggy.
                raise RuntimeError('Could not find associated AnalysisType with basal seq that matches that of '
                                   'the PotentialNewType')

        def _get_rel_abund_represented_by_current_at_of_cc(self):
            # if pnt has basal type then make sure that the basal type we are comparing to is either the match
            # or none.
            # else if pnt doesn't have basal type. compare to the type with the highest match.
            if self.pnt_or_vat.basal_seq is not None:
                # first look for exact basal match
                for vat, vat_rel_abund in self.vcc.analysis_type_obj_to_representative_rel_abund_in_cc_dict.items():
                    if vat.basal_seq == self.pnt_or_vat.basal_seq:
                        self.rel_abund_of_current_virtual_analysis_type_of_cc = vat_rel_abund
                        self.current_virtual_analysis_type_of_cc = vat
                        return
                # then look for None basal match
                for vat, vat_rel_abund in self.vcc.analysis_type_obj_to_representative_rel_abund_in_cc_dict.items():
                    if vat.basal_seq is None:
                        self.rel_abund_of_current_virtual_analysis_type_of_cc = vat_rel_abund
                        self.current_virtual_analysis_type_of_cc = vat
                        return
                self._set_most_abund_vat_to_compare_to()

            else:
                # then compare to the AnalysisType that has the highest relative abundance currently.
                self._set_most_abund_vat_to_compare_to()

        def _set_most_abund_vat_to_compare_to(self):
            top = 0
            for vat, vat_rel_abund in self.vcc.analysis_type_obj_to_representative_rel_abund_in_cc_dict.items():
                if vat_rel_abund > top:
                    self.rel_abund_of_current_virtual_analysis_type_of_cc = vat_rel_abund
                    self.current_virtual_analysis_type_of_cc = vat

        def _vcc_has_analysis_types_associated_to_it(self):
            return self.vcc.analysis_type_obj_to_representative_rel_abund_in_cc_dict

        def _pnt_or_vat_abundances_met(self):
            """This will check whether the DIV relative abundance requirements of
            either a pnt or vat are met by the vcc in question."""
            vcc_rel_abund_dict = self.vcc.ref_seq_id_to_rel_abund_dict

            total_seq_rel_abund_for_cc = []

            for ref_seq_id in self.pnt_or_vat.non_artefact_ref_seq_uid_set:
                rel_abund = vcc_rel_abund_dict[ref_seq_id]
                total_seq_rel_abund_for_cc.append(rel_abund)
                if rel_abund < self.within_clade_cutoff:
                    return False

            for ref_seq_id in self.pnt_or_vat.artefact_ref_seq_uid_set:
                rel_abund = vcc_rel_abund_dict[ref_seq_id]
                total_seq_rel_abund_for_cc.append(rel_abund)
                if rel_abund < self.unlocked_abundance:
                    return False

            self.pnt_seq_rel_abund_total_for_cc = sum(total_seq_rel_abund_for_cc)
            return True

    class StrandedCCRehomer:
        """Responsible for find new AnalysisTypes for the stranded CladeCollections to be associated with.
        Not having CCs reasigned to a type causes problems. Due to the fact that strict limits are being used to
        assign CladeCollections to the discovered AnalysisTypes it means that these stranded CCs are getting
        bad associations. A such, we need to reassociate them to the best possible types.

        For each CladeCollection object we are going to do a sort of mini type assignment
        Because a lot of the single DIV AnalysisTypes will have been gotten rid of at this point
        there is a possibility that there will not be an AnalysisType for the CCs to fit into.

        Also it may be that the CCs now fit into a lesser intra single intra type.
        e.g. B5c if original type was B5-B5s-B5c.
        In this case we would obviously rather have the type B5 associated with the clade_collection_object.
        So we should check to see if the abundance of the type that the clade_collection_object has been found
        in is larger than the abundance of the CCs Maj intra.

        If it is not then we should simply create a new type of the maj intra and associate the
        CladeCollection object to that.

        We also we need to be mindful of the fact that a CladeCollection object may not find a
        match at all, e.g. the best_type_uid will = 'None'.

        In this case we will also need to make the Maj intra the type.
        """

        def __init__(self, parent_check_type_pairing_handler):
            self.check_type_pairing_handler = parent_check_type_pairing_handler
            self.cc_to_at_match_info_holder_list = []
            self.analysis_types_just_created = []
            self.cc_to_match_object_dict = {}

        def rehome_stranded_ccs(self):
            """For each CladeCollection search through all of the AnalysisTypes and find the analysis type that
            represents the highest proportion of the CladeCollections sequences. If such a type is found,
            associate the CladeCollection to that AnalysisType. Else, create a new AnalysisType that is just the maj seq
            of the CladeCollection"""
            self._find_best_at_match_for_each_stranded_cc()
            self._convert_best_match_list_to_dict()
            self._associate_stranded_ccs()

        def _associate_stranded_ccs(self):
            for cc_obj in list(self.cc_to_match_object_dict.keys()):
                support_obj = self.cc_to_match_object_dict[cc_obj]
                self._associate_stranded_cc_to_existing_analysis_type(cc_obj, support_obj)

        def _associate_stranded_cc_to_new_maj_seq_analysis_type(self, cc_obj):
            """ Make a new type that is simply the maj intra
            NB if a type that was just the CCs maj type already existed then we would have found a suitable match above.
            I.e. at this point we know that we need to create the AnalysisType that is just the maj seq"""
            most_abund_ref_seq_of_clade_collection = cc_obj.ordered_dsss_objs[0].reference_sequence_of

            maj_seq_vat = self.check_type_pairing_handler.artefact_assessor.sp_data_analysis.virtual_object_manager. \
                vat_manager.make_vat_pre_profile_assignment(
                clade_collection_obj_list=[cc_obj], ref_seq_obj_list=[most_abund_ref_seq_of_clade_collection])

            self.add_a_type_to_cc_without_match_obj(cc=cc_obj, vat=maj_seq_vat)

            self._update_fp_to_at_dict(vat=maj_seq_vat)

        def _update_fp_to_at_dict(self, vat):
            self.check_type_pairing_handler.artefact_assessor.ref_seq_fp_set_to_analysis_type_obj_dict[
                frozenset(vat.footprint_as_ref_seq_objs_set)] = vat

        def add_a_type_to_cc_without_match_obj(self, cc, vat):
            current_type_seq_rel_abund_for_cc = []
            cc_ref_seq_abundance_dict = cc.ref_seq_id_to_rel_abund_dict
            for ref_seq in vat.footprint_as_ref_seq_objs_set:
                rel_abund = cc_ref_seq_abundance_dict[ref_seq.id]
                current_type_seq_rel_abund_for_cc.append(rel_abund)
            current_type_seq_tot_rel_abund_for_cc = sum(current_type_seq_rel_abund_for_cc)
            cc.analysis_type_obj_to_representative_rel_abund_in_cc_dict[
                vat] = current_type_seq_tot_rel_abund_for_cc

        def _associate_stranded_cc_to_existing_analysis_type(self, cc_obj, support_obj):
            self.check_type_pairing_handler.artefact_assessor.sp_data_analysis.virtual_object_manager.vat_manager. \
                add_ccs_and_reinit_virtual_analysis_type(
                vat_to_add_ccs_to=support_obj.at, list_of_clade_collection_objs_to_add=[cc_obj])
            self._remove_single_type_from_vcc_vat_dict(cc_obj)
            self.add_new_type_to_cc_with_match_obj(match_info_obj=support_obj)

        def _remove_single_type_from_vcc_vat_dict(self, cc_obj):
            if not len(cc_obj.analysis_type_obj_to_representative_rel_abund_in_cc_dict.items()) == 1:
                raise RuntimeError('Resetting dictionary of cc_obj that has more than 1 type in it')
            cc_obj.analysis_type_obj_to_representative_rel_abund_in_cc_dict = {}

        def add_new_type_to_cc_with_match_obj(self, match_info_obj):
            match_info_obj.cc.analysis_type_obj_to_representative_rel_abund_in_cc_dict[
                match_info_obj.at] = match_info_obj.rel_abund_of_at_in_cc

        def _get_most_abundnant_ref_seq_info_for_cc(self, cc_obj):
            most_abund_dss_of_clade_collection = cc_obj.ordered_dsss_objs[0]
            rel_abund = most_abund_dss_of_clade_collection.abundance / cc_obj.total_seq_abundance
            ref_seq = most_abund_dss_of_clade_collection.reference_sequence_of
            return rel_abund, ref_seq

        def _convert_best_match_list_to_dict(self):
            for support_obj in self.cc_to_at_match_info_holder_list:
                self.cc_to_match_object_dict[support_obj.cc] = support_obj

        def _find_best_at_match_for_each_stranded_cc(self):
            """Start the worker that will search through all of the AnalysisTypes and find the AnalysisType that represents
            the highest relative abundance in the CladeCollection. The function will add a CCToATMatchInfoHolder for
            each CladeCollection and Analysis best match found. Or, if no match is found will add None."""
            for cc in self.check_type_pairing_handler.stranded_ccs:
                print(f'Searching for best match to {cc}')
                sccats = StrandedCladeCollectionAnalysisTypeSearcher(
                    parent_stranded_cc_rehomer=self, virtual_clade_collection_object=cc)
                if not sccats.search_analysis_types():
                    self._associate_stranded_cc_to_new_maj_seq_analysis_type(cc)



class CheckArtefactDIVVATAssociations(CheckVCCToVATAssociations):
    """Check to see whether there are additional VirtualCladeCollections that will support a given VirtualAnalysisType
    now that artefact or unlocked DIVs are being taken into account. If there are additional VCC that support,
    reassociate these VCCs."""
    def __init__(self, parent_artefact_assessor, vat_to_check):
        super().__init__(parent_artefact_assessor=parent_artefact_assessor)
        self.vat_to_check = vat_to_check

    def check_artefact_div_vat_associations(self):
        self._set_vccs_to_check()
        self._assess_support_of_pnt_or_vat(pnt_or_vat=self.vat_to_check)
        if self.list_of_loss_of_support_info_holder_objs:

            print(f'Found {len(self.list_of_loss_of_support_info_holder_objs)} '
                  f'VirtualCladeCollections that support {self.vat_to_check.name}')
            print('Reassociating VirtualCladeCollections...')

            self._update_cc_info_for_ccs_that_support_new_type(vat_to_add_to_vcc=self.vat_to_check)
            self._reinit_or_del_affected_types_and_create_stranded_cc_list()
            self._reassociate_stranded_ccs_if_necessary()
            return True
        else:
            return False

    def _set_vccs_to_check(self):
        # its important that we don't check the VirtualCladeCollections that are already associated with the
        # VirtualAnalysisType. This way, any VirtualCladeCollections that we do find to be a better match to the
        # VAT in question we can add to the VAT and remove from their old VAT, just like we did when looking
        # for support for the PotentialNewType in the first round of artefact assessment.
        self.list_of_vcc_objs_to_check = [
            vcc for vcc in self.artefact_assessor.list_of_vccs if
            vcc.clade == self.vat_to_check.clade if
            vcc.footprint_as_frozen_set_of_ref_seq_uids.issuperset(
                self.vat_to_check.ref_seq_uids_set) if
            vcc not in self.vat_to_check.clade_collection_obj_set_profile_discovery]



class CheckTypePairingHandler(CheckVCCToVATAssociations):
    """Check to see whether a PotentialNewType has sufficient support to become a new VirtualAnalysisType. If so
    then redistribute the VirtualCladeCollections that noew support it (and no long support the VATs they were
    previously supporting."""
    def __init__(self, parent_artefact_assessor, vat_a, vat_b):
        super().__init__(parent_artefact_assessor=parent_artefact_assessor)
        # AnalysisTypeAretefactInfoHolder for types a and b
        self.vat_a = vat_a
        self.vat_b = vat_b
        # NB that the non_artefact ref seqs being the same is a prerequisite of doing a comparison
        # as such we can set the below pnt non artefact ref seqs to match just one of the types being compared
        self.pnt = self._init_pnt(self.vat_a, self.vat_b)

        # Attribute once support found
        self.new_virtual_analysis_type_from_pnt = None

    def check_type_pairing(self):
        if self._pnt_profile_already_an_existing_analysis_type_profile():
            print(f'Assessing support for potential new type:{self.pnt.name}')
            print('Potential new type already exists')
            return False

        self._set_vccs_to_check()

        self._assess_support_of_pnt_or_vat(pnt_or_vat=self.pnt)

        if self._pnt_or_vat_has_support():
            self._make_new_at_from_pnt_and_update_dicts()
            self._reassociate_stranded_ccs_if_necessary()
        else:
            print(f'Assessing support for potential new type:{self.pnt.name}')
            print('Insufficient support for potential new type')
            return False
        return True

    def _set_vccs_to_check(self):
        self.list_of_vcc_objs_to_check = [vcc for vcc in self.artefact_assessor.list_of_vccs if
                                          vcc.clade == self.vat_a.clade if
                                          vcc.footprint_as_frozen_set_of_ref_seq_uids.issuperset(
                                              self.pnt.ref_seq_uids_set)]

    def _pnt_profile_already_an_existing_analysis_type_profile(self):
        return frozenset(
            self.pnt.ref_seq_objects_set) in self.artefact_assessor.ref_seq_fp_set_to_analysis_type_obj_dict.keys()

    def _make_new_at_from_pnt_and_update_dicts(self):
        self._make_analysis_type_from_pnt()
        self._update_fp_to_at_dict_from_pnt()
        self._update_cc_info_for_ccs_that_support_new_type(vat_to_add_to_vcc=self.new_virtual_analysis_type_from_pnt)
        self._reinit_or_del_affected_types_and_create_stranded_cc_list()

    def _make_analysis_type_from_pnt(self):
        print('\nSupport found. Creating new type.')
        self.new_virtual_analysis_type_from_pnt = self.artefact_assessor.sp_data_analysis.virtual_object_manager.vat_manager.make_vat_pre_profile_assignment(
            clade_collection_obj_list=[support_obj.cc for support_obj in self.list_of_loss_of_support_info_holder_objs],
            ref_seq_obj_list=self.pnt.ref_seq_objects_set)
        print(f'{self.new_virtual_analysis_type_from_pnt.name} created.')

    def _update_fp_to_at_dict_from_pnt(self):
        self.artefact_assessor.ref_seq_fp_set_to_analysis_type_obj_dict[frozenset(
            self.new_virtual_analysis_type_from_pnt.footprint_as_ref_seq_objs_set
        )] = self.new_virtual_analysis_type_from_pnt

    def _init_pnt(self, artefact_info_a, artefact_info_b):
        name_of_ref_seqs_in_pnt = [ref_seq.name for ref_seq in
                 artefact_info_a.footprint_as_ref_seq_objs_set.union(artefact_info_b.footprint_as_ref_seq_objs_set)]

        return PotentialNewType(
            ref_seq_uids_set=self.vat_a.ref_seq_uids_set.union(self.vat_b.artefact_ref_seq_uid_set),
            artefact_ref_seq_uid_set=self.vat_a.artefact_ref_seq_uid_set.union(self.vat_b.artefact_ref_seq_uid_set),
            non_artefact_ref_seq_uid_set=self.vat_a.non_artefact_ref_seq_uid_set,
            list_of_ref_seq_names=name_of_ref_seqs_in_pnt,
            resf_seq_obj_set=artefact_info_a.footprint_as_ref_seq_objs_set.union(
                artefact_info_b.footprint_as_ref_seq_objs_set),
            force_basal_lineage_separation=self.artefact_assessor.sp_data_analysis.force_basal_lineage_separation)

class StrandedCladeCollectionAnalysisTypeSearcher:
    """For a given stranded CladeCollection it will search the AnalysisTypes to identify the the AnalysisType that
    representes the greatest proportion of the CladeCollections sequences."""
    def __init__(self, parent_stranded_cc_rehomer, virtual_clade_collection_object):
        self.stranded_cc_rehomer = parent_stranded_cc_rehomer
        self.cc = virtual_clade_collection_object
        self.best_rel_abund_of_at_in_cc = 0
        self.best_match_vat = None

    def search_analysis_types(self):
        self._find_at_found_in_cc_with_highest_rel_abund()
        self._if_match_put_in_output_else_return_false()
        if self.best_match_vat is None:
            return False
        else:
            return True

    def _if_match_put_in_output_else_return_false(self):
        if self.best_match_vat is not None:
            print(f'Best match to {self.cc} is {self.best_match_vat}')
            self.stranded_cc_rehomer.cc_to_at_match_info_holder_list.append(
                CCToATMatchInfoHolder(
                    vat=self.best_match_vat,
                    vcc=self.cc,
                    rel_abund_of_at_in_cc=self.best_rel_abund_of_at_in_cc))
        else:
            print(f'No VirtualAnalysisMatch found for {self.cc}')

    def _find_at_found_in_cc_with_highest_rel_abund(self):
        for vat in self.stranded_cc_rehomer.check_type_pairing_handler.artefact_assessor.virtual_analysis_type_dict.values():
            atficc = AnalysisTypeFoundInCladeCollection(virtual_analysis_type=vat, virtual_calde_collection=self.cc)
            if atficc.search_for_at_in_cc():
                self.best_match_vat = vat
                self.best_rel_abund_of_at_in_cc = atficc.rel_abund_of_at_in_cc





class AnalysisTypeFoundInCladeCollection:
    """Class sees whether a given AnalysisType is found within a CladeCollection
    and if so calculates the relative abundance of sequences that the AnalysisType represents."""
    def __init__(self, virtual_analysis_type, virtual_calde_collection):
        self.vat = virtual_analysis_type
        self.vcc = virtual_calde_collection
        self.rel_abund_of_at_in_cc = None

    def search_for_at_in_cc(self):
        if self._at_non_artefact_seqs_found_in_cc_above_cutoff_seqs():
            if self._at_artefact_seqs_found_in_cc():
                self._get_rel_abund_of_at_in_cc()
                return self.rel_abund_of_at_in_cc
        return False

    def _get_rel_abund_of_at_in_cc(self):
        self.rel_abund_of_at_in_cc = sum(
            [self.vcc.ref_seq_id_to_rel_abund_dict[rs.id] for rs in self.vat.footprint_as_ref_seq_objs_set])

    def _at_artefact_seqs_found_in_cc(self):
        """Here we check to see if the artefact seqs are found in the CC. Perhaps strictly, strictly, strictly speaking
        we should be looking to see if each of these abundances is above the unlocked relative abundance, but given
        that the unlocked rel abund is no 0.0001 this is so low that I think we can just be happy if the seqs
        exist in the CladeCollection."""
        return self.vat.artefact_ref_seq_uid_set.issubset(
            self.vcc.footprint_as_frozen_set_of_ref_seq_uids)

    def _at_non_artefact_seqs_found_in_cc_above_cutoff_seqs(self):
        return self.vat.non_artefact_ref_seq_uid_set.issubset(self.vcc.above_cutoff_ref_seqs_id_set)


class CCToATMatchInfoHolder:
    """Responsible for holding the information of AnalysisType and CladeCollection and relative abundance
    that the AnalysisType represents in the CladeCollection when doing the CheckPNTSupport and we
    find that the PNT is a better fit thatn the current AnalysisType."""
    def __init__(self, vcc, vat, rel_abund_of_at_in_cc=None):
        self.cc = vcc
        self.at = vat
        self.rel_abund_of_at_in_cc = rel_abund_of_at_in_cc






class AnalysisTypeCreator:
    """Create AnalysisType objects from the supported initial type profiles that have been generated in the
    SupportedFootprintIdentifier"""
    def __init__(self, parent_sp_data_analysis):
        self.sp_data_analysis = parent_sp_data_analysis

    def create_analysis_types(self):
        print(f'\n\nCreating analysis types clade {self.sp_data_analysis.current_clade}')
        for initial_type in self.sp_data_analysis.list_of_initial_types_after_collapse:
            self._create_new_virtual_analysis_type_from_initial_type(initial_type)

    def _create_new_virtual_analysis_type_from_initial_type(self, initial_type):
        new_virtual_analysis_type = self.sp_data_analysis.virtual_object_manager.vat_manager.make_vat_pre_profile_assignment(
            clade_collection_obj_list=initial_type.clade_collection_list,
            ref_seq_obj_list=initial_type.profile)

        print(f'Creating virtual analysis type: {new_virtual_analysis_type.name}')


class SupportedFootPrintIdentifier:
    """This class is responsible for identifying the footprints that are found in a sufficient number of clade
    collections to warrant becoming AnalysisType objects.
    Operates by working with the longest footprints and trying to collapse those that aren't already supported into
    those of length n-1 etc etc. It also tries to find support amongst the collection of footprints of length n
    e.g. if you have 1-2-3-4, 1-2-3-5, 1-2-3-6, 1-2-3-7, 1-2-3-8, it will pull out the 1-2-3 footprint as supported
    even if the 1-2-3 footprint doesnt already exist as an n=3 footprint.
    20200511. We are going to make the exlusion of C3/C1/C15 sequences flagged as option. For the time
    being it will be the default NOT to apply it. If we then do choose to apply it, we should do it properly where
    we classify sequences as C3 basal, C15 basal or C1 basal within the cladocopium. We would then not allow basal
    sequences to mix in the creation of type profiles. The easiest way to do this could be to treat them in the same
    way that we treat genera. But... this will be some way down the line. I think we're going to find that actually
    its perfectly viable for profiles to contain both C1 and C3. C15 might be a different case.
    We are also taking into account not allowing analysis types to contain the C3 and the C15 sequences. We refer to
    these as the basal sequences. All footprints that contain more than one basal sequences will automatically be
    put into the unsupported list at first, irrespective of their support. Then they will attempt to be collapsed into
    other footprints just like all other unsupported footprints. If a suitable footprint is foud that they can be
    collapsed into then they will be but the pulled out sequqences can only contain one of the basal
    sequences. In this way, a single clade collection can go towards the support of two different
    footprints, one for with C3 and one with C15.
    as """

    def __init__(self, clade_footprint_dict, parent_sp_data_analysis):
        self.sp_data_analysis = parent_sp_data_analysis
        self.clade_fp_dict = clade_footprint_dict
        self.supported_list = []
        self.unsupported_list = []
        self.initial_types_list = []
        self._init_initial_types_list(force_basal_lineage_separation=self.sp_data_analysis.force_basal_lineage_separation)
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

    def _init_initial_types_list(self, force_basal_lineage_separation):
        for footprint_key, footprint_representative in self.clade_fp_dict.items():
            self.initial_types_list.append(
                InitialType(
                    reference_sequence_set=footprint_key,
                    clade_collection_list=footprint_representative.cc_list,
                    force_basal_lineage_separation=force_basal_lineage_separation,
                    maj_dsss_list=footprint_representative.maj_dss_seq_list))

    def _verify_that_all_cc_associated_to_an_initial_type(self):
        set_of_ccs_of_clade = set([cc.id for cc in self.sp_data_analysis.ccs_of_analysis if cc.clade == self.sp_data_analysis.current_clade])
        set_of_ccs_found_in_init_types = set()
        for init_type in self.initial_types_list:
            set_of_ccs_found_in_init_types.update([cc.id for cc in init_type.clade_collection_list])

        if set_of_ccs_of_clade.issuperset(set_of_ccs_found_in_init_types):
            if len(set_of_ccs_of_clade.difference(set_of_ccs_found_in_init_types)) == 0:
                print("\n\nAll CladeCollections successfuly associated to an InitialType\n")
            else:
                raise RuntimeError(f'{len(set_of_ccs_of_clade.difference(set_of_ccs_found_in_init_types))} '
                                   f'CladeCollections were not associated to an InitialType')
        else:
            raise RuntimeError(f'{len(set_of_ccs_of_clade.difference(set_of_ccs_found_in_init_types))} '
                                   f'CladeCollections were not associated to an InitialType')

    def identify_supported_footprints(self):
        # for each length starting at max and dropping by 1 with each increment
        longest_footprint = self._return_len_of_longest_fp()
        for n in range(longest_footprint, 0, -1):
            self.current_n = n
            self._update_supported_unsupported_lists_for_n()
            self._collapse_n_len_initial_types_into_minus_one_initial_types()

            if self.current_n > 2:
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
        
        self._make_maj_seq_initial_types()

        self._verify_that_all_cc_associated_to_an_initial_type()
        return self.initial_types_list

    def _make_maj_seq_initial_types(self):
        # At this point we also want to make initial types that
        # are the majority sequence of each of the vccs.
        # We do this because we were ending up with the case where a sample
        # was not being assciated to a profile, because the only single DIV
        # profile that matched did not represent > 5% of the rel seq abund.
        # It makes sense that the Maj seq of a vcc should always be a profile itself
        # if it is not already included in another profile.
        
        # maj_seq_ids_of_ccs = set([int(cc.footprint.split(",")[0]) for cc in self.sp_data_analysis.ccs_of_analysis if cc.clade == self.sp_data_analysis.current_clade])
        # maj_seq_of_vccs = set([vcc.ordered_dsss_objs[0].reference_sequence_of for vcc in self.sp_data_analysis.virtual_object_manager.vcc_manager.vcc_dict.values() if vcc.clade == self.sp_data_analysis.current_clade])
        # maj_ref_seq_id_to_rs_obj_dict = 
        # maj_ref_seq_id_to_rs_obj_dict = {rs.id: rs for rs in  ReferenceSequence.objects.filter(id__in=maj_seq_ids_of_ccs)}
        
        # Because we can't modify the initial_types_list at the same time as parsing through it
        # we will use a default dict to store the references sequence to the list of clade collections that should be associated to it
        rs_to_vcc_maj_seq_dict = defaultdict(list)
        for vcc in [_ for _ in self.sp_data_analysis.virtual_object_manager.vcc_manager.vcc_dict.values() if _.clade == self.sp_data_analysis.current_clade]:
            vcc_maj_rs = vcc.ordered_dsss_objs[0].reference_sequence_of
            rs_to_vcc_maj_seq_dict[vcc_maj_rs].append(vcc)
        
        # Now parse through the rs_to_vcc_maj_seq_dict and check if the single seq inital type already existis
        # If it does, then check whether the ccs already support, if not add support, if so, move on
        # if it doesn't, create it and add all of the support
        profile_to_initial_type_dict = {_.profile: _ for _ in self.initial_types_list}
        for rs, vcc_list in rs_to_vcc_maj_seq_dict.items():
            try:
                initial_type = profile_to_initial_type_dict[frozenset([rs])]
                # Create a new initial type to replace the current initial type
                # Make sure to add any vccs that are in the current intial type to the new initial type
                # if not already in the vcc_list
                # So that we maintain that all vcc are associated with an initial type
                vcc_list_to_add = list(set(vcc_list + initial_type.clade_collection_list))
                # Remove the old one from the self.initial_types_list
                # and add the new one
                new_initial_type = InitialType(reference_sequence_set=frozenset([rs]), clade_collection_list=vcc_list_to_add, force_basal_lineage_separation=self.sp_data_analysis.force_basal_lineage_separation)
                self.initial_types_list.remove(initial_type)
                self.initial_types_list.append(new_initial_type)
            except KeyError:
                # Create a new intial type and add it to the self.list_of_initial_types
                new_initial_type = InitialType(reference_sequence_set=frozenset([rs]), clade_collection_list=vcc_list, force_basal_lineage_separation=self.sp_data_analysis.force_basal_lineage_separation)
                self.initial_types_list.append(new_initial_type)

    def _del_type_to_collapse(self, k, synth_fp):
        # then the big init_type no longer contains any
        # of its original maj ref seqs and
        # so should be delted.
        self.initial_types_list.remove(self.synthetic_fp_dict[synth_fp][k])
        self.unsupported_list.remove(self.synthetic_fp_dict[synth_fp][k])

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
        temp_dict = {}
        for footprint_set in [_.profile for _ in self.list_of_initial_types_of_len_n]:
            temp_dict.update({
                frozenset(tup): [] for tup in itertools.combinations(footprint_set, self.current_n-1)})
        sys.stdout.write(f'\rGenerated {len(temp_dict.items())} synthetic footprints')
        self.synthetic_fp_dict = temp_dict

    def _associate_un_sup_init_types_to_synth_footprints(self):
        # Now go through each of the (n-1) footprints and see if they
        # fit into a footprint in the unsuported list
        sys.stdout.write('\rChecking new set of synthetic types')
        for synth_fp in self.synthetic_fp_dict.keys():
            if self.sp_data_analysis.force_basal_lineage_separation:
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

            if not self.large_fp_to_collapse_list:
                break

            for q in range(len(self.large_fp_to_collapse_list)):
                large_fp_to_collapse = self.large_fp_to_collapse_list[q]

                if self._remove_fp_to_collapse_from_unsupported_if_now_supported(large_fp_to_collapse):
                    continue

                fp_collapser = FootprintCollapser(
                    footprint_to_collapse_index=q,
                    parent_supported_footprint_identifier=self)
                fp_collapser.collapse_footprint()

    def _remove_fp_to_collapse_from_unsupported_if_now_supported(self, large_fp_to_collapse):
        if self.sp_data_analysis.force_basal_lineage_separation:
            if large_fp_to_collapse.support >= \
                    self.required_support and not large_fp_to_collapse.contains_multiple_basal_sequences:
                # Then this type has already had some other leftovers put into it so that it now has the required
                # support. In this case we can remove the type from the unsupported list and continue to the next
                self.unsupported_list.remove(large_fp_to_collapse)
                return True
            else:
                return False
        else:
            if large_fp_to_collapse.support >= self.required_support:
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
            if self.sp_data_analysis.force_basal_lineage_separation:
                if initial_type.support >= self.required_support and not initial_type.contains_multiple_basal_sequences:
                    self.supported_list.append(initial_type)
                else:
                    self.unsupported_list.append(initial_type)
            else:
                if initial_type.support >= self.required_support:
                    self.supported_list.append(initial_type)
                else:
                    self.unsupported_list.append(initial_type)

class UnsupportedFootprintFinalCollapser:
    def __init__(self, parent_supported_footprint_identifier):
        self.supported_footprint_identifier = parent_supported_footprint_identifier
        self.fp_to_collapse = None
        self.matching_initial_type = None
        self.ccs_not_already_in_maj_seq_initial_type = True

    def collapse_unsup_footprints_to_maj_refs(self):
        """Collapse unsupported footprints into an InitialType that is just their majority sequence.
        Bear in mind that when we are at n==1 we will still be visiting this code. At this point we should first check
        to see if the InitialType to which a CladeCollection is currently associated is its majority sequence. If
        it is, then there is nothing to do here as this is the InitialType that we were looking to collapse into
        anyway."""
        while self.supported_footprint_identifier.unsupported_list:
            self.fp_to_collapse = self.supported_footprint_identifier.unsupported_list[0]
            if self.supported_footprint_identifier.current_n == 1:
                    if self._all_ccs_already_associated_to_maj_seq_init_type():
                        if self.fp_to_collapse in self.supported_footprint_identifier.unsupported_list:
                            self.supported_footprint_identifier.unsupported_list.remove(self.fp_to_collapse)
                    elif self._only_some_ccs_already_associated_to_maj_seq_init_type():
                        raise RuntimeError(f'InitialType {self.fp_to_collapse} is supported by some CladeCollections that'
                                           f'have a maj seq of its profile, whilst other CladeCollections'
                                           f'have a different maj seqs.')
                    else:
                        # none of the CladeCollections supporting the InitialType have its profile as their maj seqs
                        self._do_maj_ref_init_type_collapse_and_clean_up()
            else:
                self._do_maj_ref_init_type_collapse_and_clean_up()

    def _only_some_ccs_already_associated_to_maj_seq_init_type(self):
        count = 0
        for i in range(len(self.fp_to_collapse.clade_collection_list)):
            if not self.fp_to_collapse.majority_sequence_list[i].reference_sequence_of in self.fp_to_collapse.profile:
                count += 1
        if 0 < count < len(self.fp_to_collapse.clade_collection_list):
            return True
        return False

    def _do_maj_ref_init_type_collapse_and_clean_up(self):
        for i in range(len(self.fp_to_collapse.clade_collection_list)):  # for each cc
            for maj_dss in self.fp_to_collapse.majority_sequence_list[i]:  # for each maj_ref_seq
                if self._inital_type_exists_with_maj_ref_seq_as_profile(maj_dss):
                    self._add_unsup_type_info_to_smll_match_type(i, maj_dss)
                else:
                    self._create_new_maj_seq_init_type(i, maj_dss)
        self._del_type_to_collapse()

    def _all_ccs_already_associated_to_maj_seq_init_type(self):
        """If all of the CCs are already in an InitialType that is their
        maj ref seq then there is no need to do anything.
        """
        for i in range(len(self.fp_to_collapse.clade_collection_list)):
            if not self.fp_to_collapse.majority_sequence_list[i][0].reference_sequence_of in \
                   self.fp_to_collapse.profile:
                return False
        return True
    def _create_new_maj_seq_init_type(self, i, maj_dss):
        new_initial_type = InitialType(
            reference_sequence_set=frozenset([maj_dss.reference_sequence_of]),
            clade_collection_list=[self.fp_to_collapse.clade_collection_list[i]],
            force_basal_lineage_separation=self.supported_footprint_identifier.sp_data_analysis.force_basal_lineage_separation)
        self.supported_footprint_identifier.initial_types_list.append(new_initial_type)

    def _add_unsup_type_info_to_smll_match_type(self, i, maj_dss):
        self.matching_initial_type.clade_collection_list.append(self.fp_to_collapse.clade_collection_list[i])
        self.matching_initial_type.majority_sequence_list.append([maj_dss])
        self.matching_initial_type.support = len(self.matching_initial_type.clade_collection_list)

    def _inital_type_exists_with_maj_ref_seq_as_profile(self, maj_dss):
        """Check to see if an initial type with profile of that maj_dss refseq already exists"""
        for initT in [init for init in self.supported_footprint_identifier.initial_types_list if init.profile_length == 1]:
            if maj_dss.reference_sequence_of in initT.profile:
                self.matching_initial_type = initT
                return True
        return False

    def _del_type_to_collapse(self):
        """Once we have associated each of the cc of an initial type to collapse to existing or new initial types,
        delete.
        """
        self.supported_footprint_identifier.initial_types_list.remove(self.fp_to_collapse)
        if self.fp_to_collapse in self.supported_footprint_identifier.unsupported_list:
            self.supported_footprint_identifier.unsupported_list.remove(self.fp_to_collapse)

class CollapseAssessor:
    """Responsible for assessing whether an unsupported large initial type can be collapsed into a given
    small initial type.
    """
    def __init__(self, parent_supported_footprint_identifier, longer_intial_type, shorter_initial_type):
        self.supported_footprint_identifier = parent_supported_footprint_identifier
        self.longer_intial_type = longer_intial_type
        self.shorter_initial_type = shorter_initial_type

    def assess_collapse(self):
        if self._if_short_initial_type_suitable_for_collapse():
            if self.supported_footprint_identifier.sp_data_analysis.force_basal_lineage_separation:
                if self.longer_intial_type.contains_multiple_basal_sequences:
                    if self.does_small_footprint_contain_the_required_ref_seqs_of_the_large_footprint():
                        # score = number of samples big was found in plus num samples small was found in
                        self._if_highest_score_so_far_assign_big_fp_to_smll_fp_for_collapse()
                else:
                    if self.longer_intial_type.set_of_maj_ref_seqs.issubset(self.shorter_initial_type.profile):
                        self._if_highest_score_so_far_assign_big_fp_to_smll_fp_for_collapse()
            else:
                if self.longer_intial_type.set_of_maj_ref_seqs.issubset(self.shorter_initial_type.profile):
                    self._if_highest_score_so_far_assign_big_fp_to_smll_fp_for_collapse()

    def _if_highest_score_so_far_assign_big_fp_to_smll_fp_for_collapse(self):
        score = self.longer_intial_type.support + self.shorter_initial_type.support
        if score > self.supported_footprint_identifier.top_score:
            self.supported_footprint_identifier.top_score = score
            self.supported_footprint_identifier.repeat = True
            self.supported_footprint_identifier.collapse_dict[self.longer_intial_type] = self.shorter_initial_type

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
        if self.supported_footprint_identifier.sp_data_analysis.force_basal_lineage_separation:
            multi_basal = self.shorter_initial_type.contains_multiple_basal_sequences
            return self.shorter_initial_type.profile.issubset(self.longer_intial_type.profile) and not multi_basal
        else:
            return self.shorter_initial_type.profile.issubset(self.longer_intial_type.profile)

class SyntheticFootprintCollapser:
    def __init__(self, parent_supported_footprint_identifier):
        self.supported_footprint_identifier = parent_supported_footprint_identifier
        self.current_synth_fp = None
        self.current_fp_to_collapse = None
        self.matching_existing_init_type_outer = None
        self.matching_existing_init_type_inner = None

    def collapse_to_synthetic_footprints(self):
        for synth_fp in self.supported_footprint_identifier.ordered_sig_synth_fps:  # for each synth_fp
            self.current_synth_fp = synth_fp
            for k in range(len(self.supported_footprint_identifier.synthetic_fp_dict[synth_fp])):  # for each assoc. initial type
                self.current_fp_to_collapse = self.supported_footprint_identifier.synthetic_fp_dict[synth_fp][k]
                if self._validate_init_type_is_unsup_and_synth_fp_is_subset():
                    if self._synth_fp_matches_existing_intial_type():
                        if self.supported_footprint_identifier.sp_data_analysis.force_basal_lineage_separation:
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
                            self._absorb_type_into_match_and_del()
                    else:
                        # then the synth footprint was not already represented by an existing initial type
                        # Check to see if the big type contains mutiple basal.
                        # If it does then we should extract as above. This should be exactly the
                        # same code as above but extracting into a new type rather than an existing one
                        if self.supported_footprint_identifier.sp_data_analysis.force_basal_lineage_separation:
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
                        else:
                            self._create_new_initial_type_from_synth_type_and_del_type_to_collapse()

    def _create_new_init_type_and_add_to_init_type_list(self):
        new_blank_initial_type = InitialType(
            reference_sequence_set=self.current_synth_fp, clade_collection_list=list(
                self.current_fp_to_collapse.clade_collection_list),
            force_basal_lineage_separation=self.supported_footprint_identifier.sp_data_analysis.force_basal_lineage_separation)
        self.supported_footprint_identifier.initial_types_list.append(new_blank_initial_type)
        return new_blank_initial_type

    def _create_new_initial_type_from_synth_type_and_del_type_to_collapse(self):
        self._create_new_init_type_and_add_to_init_type_list()
        self._del_type_to_collapse()

    def _absorb_type_into_match_and_del(self):
        self.matching_existing_init_type_outer.absorb_large_init_type(self.current_fp_to_collapse)
        self.supported_footprint_identifier.initial_types_list.remove(self.current_fp_to_collapse)
        self.supported_footprint_identifier.unsupported_list.remove(self.current_fp_to_collapse)

    def _del_type_to_collapse(self):
        """Then the initial type to collapse no longer contains any of its original maj ref seqs and should be deleted.
        """
        self.supported_footprint_identifier.initial_types_list.remove(self.current_fp_to_collapse)
        if self.current_fp_to_collapse in self.supported_footprint_identifier.unsupported_list:
            self.supported_footprint_identifier.unsupported_list.remove(self.current_fp_to_collapse)

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
        if self.current_fp_to_collapse.profile_length < self.supported_footprint_identifier.current_n:
            self.supported_footprint_identifier.unsupported_list.remove(self.current_fp_to_collapse)

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
        if self.matching_existing_init_type_inner in self.supported_footprint_identifier.unsupported_list:
            self.supported_footprint_identifier.unsupported_list.remove(self.matching_existing_init_type_inner)
        self.supported_footprint_identifier.initial_types_list.remove(self.matching_existing_init_type_inner)

    def _extracted_initial_type_fp_now_matches_existing_initial_type(self):
        """If the initial type to collapse still contains maj
        ref sequences, then it is still a profile that we needs to be assessed for collapse.
        Now need to check if it's new profile (i.e. after extraction) matches that of any of the other initial types."""
        for j in range(len(self.supported_footprint_identifier.initial_types_list)):
            if self.supported_footprint_identifier.initial_types_list[j].profile == self.current_fp_to_collapse.profile:
                if self.supported_footprint_identifier.initial_types_list[j] != self.current_fp_to_collapse:
                    self.matching_existing_init_type_inner = self.supported_footprint_identifier.initial_types_list[j]
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
        for i in range(len(self.supported_footprint_identifier.initial_types_list)):
            # for initT_one in initial_types_list:
            if self.supported_footprint_identifier.initial_types_list[i].profile == self.current_synth_fp:
                self.matching_existing_init_type_outer = self.supported_footprint_identifier.initial_types_list[i]
                return True
        return False


    def _validate_init_type_is_unsup_and_synth_fp_is_subset(self):
        """Then this footprint hasn't been collapsed anywhere yet.
        We also check to make sure that the initial type's profile hasn't been modified (i.e. by extraction) and check
        that the synthetic fp is still a sub set of the initial type's profile.
        """
        if self.current_fp_to_collapse in self.supported_footprint_identifier.unsupported_list:
            if self.current_synth_fp.issubset(self.current_fp_to_collapse.profile):
                return True
        return False

class FootprintCollapser:
    """Responsible for collapsing the long initial type into a short initial type."""
    def __init__(self, parent_supported_footprint_identifier, footprint_to_collapse_index):
        self.supported_footprint_identifier = parent_supported_footprint_identifier
        self.fp_index = footprint_to_collapse_index
        self.long_initial_type = self.supported_footprint_identifier.large_fp_to_collapse_list[self.fp_index]
        self.short_initial_type = self.supported_footprint_identifier.collapse_dict[self.long_initial_type]
        # bool on whether we are simply collapsing big into small or whether we need to
        # extract certain sequences of the big footprint to go into the small (i.e. if the long
        # fp contains multiple basal seqs
        self.should_extract_not_delete = self.long_initial_type.contains_multiple_basal_sequences
        # whether an extracted
        self.match = False

    def collapse_footprint(self):
        if self.supported_footprint_identifier.sp_data_analysis.force_basal_lineage_separation:
            if not self.should_extract_not_delete:
                self._collapse_long_into_short_no_extraction()
            else:
                self._collapse_long_into_short_by_extraction()
                self.match = False
                for p in range(len(self.supported_footprint_identifier.initial_types_list)):
                    intial_type_being_checked = self.supported_footprint_identifier.initial_types_list[p]
                    if self._extracted_long_initial_type_now_has_same_fp_as_another_initial_type(intial_type_being_checked):
                        self.match = True
                        self._absorb_long_initial_type_into_matching_intial_type(intial_type_being_checked)
                        break
                if not self.match:
                    if self.long_initial_type.profile_length < self.supported_footprint_identifier.current_n:
                        self.supported_footprint_identifier.unsupported_list.remove(self.long_initial_type)
        else:
            self._collapse_long_into_short_no_extraction()

    def _long_initial_type_has_sufficient_support(self):
        if self.long_initial_type.support >= self.supported_footprint_identifier.required_support:
            if not self.long_initial_type.contains_multiple_basal_sequences:
                return True
        return False

    def _absorb_matching_initial_type_into_long_intial_type(self, intial_type_being_checked):
        self.long_initial_type.absorb_large_init_type(intial_type_being_checked)
        if intial_type_being_checked in self.supported_footprint_identifier.unsupported_list:
            self.supported_footprint_identifier.unsupported_list.remove(intial_type_being_checked)
        self.supported_footprint_identifier.initial_types_list.remove(intial_type_being_checked)

    def _absorb_long_initial_type_into_matching_intial_type(self, intial_type_being_checked):
        """Collapse the long initial type into the matching initial type.
        Because this collapsing can cause the matching type that is also in the collapse
        dict to gain support we will also check to if each of the types in the collapse dict
        have gained sufficient support to no longer need collapsing.
        We will do this earlier in the process, not here.
        """
        intial_type_being_checked.absorb_large_init_type(self.long_initial_type)
        self.supported_footprint_identifier.unsupported_list.remove(self.long_initial_type)
        self.supported_footprint_identifier.initial_types_list.remove(self.long_initial_type)

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
        self.supported_footprint_identifier.initial_types_list.remove(self.long_initial_type)
        self.supported_footprint_identifier.unsupported_list.remove(self.long_initial_type)

class InitialType:
    def __init__(self, reference_sequence_set, clade_collection_list, force_basal_lineage_separation, maj_dsss_list=False):
        self.profile = reference_sequence_set
        self.profile_length = len(self.profile)
        if force_basal_lineage_separation:
            self.contains_multiple_basal_sequences, self.basalSequence_list = \
                self.check_if_initial_type_contains_basal_sequences()
        else:
            self.contains_multiple_basal_sequences, self.basalSequence_list = False, []
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
            for vcc in self.clade_collection_list:
                temp_dsss_list = []
                data_set_sample_sequence_list = vcc.ordered_dsss_objs
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
                master_dsss_list.append([maj_dsss_list[i]])
                set_of_majority_reference_sequences.add(maj_dsss_list[i].reference_sequence_of)

        return master_dsss_list, set_of_majority_reference_sequences

    def create_majority_sequence_list_for_inital_type_from_scratch(self):
        # This will be like above but will not start with the maj_dsss_list
        # we will go through each of the cladeCollections of the type and get the maj sequence for the type

        # if the init type has multiple basal sequences then we will have to find the actual maj and the basal seq dsss
        set_of_majority_reference_sequences = set()
        master_dsss_list = []
        if self.contains_multiple_basal_sequences:
            for vcc in self.clade_collection_list:
                temp_dsss_list = []
                dsss_in_cc_in_profile = [dsss for dsss in vcc.ordered_dsss_objs if dsss.reference_sequence_of in self.profile]
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
            for vcc in self.clade_collection_list:
                temp_dsss_list = []
                dsss_in_cc_in_profile = [dsss for dsss in vcc.ordered_dsss_objs if
                                         dsss.reference_sequence_of in self.profile]
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
        self.support = len(self.clade_collection_list)
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
    def __init__(self, cc, maj_dss_seq_list):
        self.cc_list = [cc]
        self.maj_dss_seq_list = [maj_dss_seq_list]

class FootprintDictGenerationCCInfoHolder:
    """An object used in the FootprintDictPopWorker to hold information for each of the clade collections
    that will be used to popualte the clade footprint dicts"""
    def __init__(self, footprint, clade_index, cc, maj_ref_seq):
        self.footprint = footprint
        self.clade_index = clade_index
        self.cc = cc
        self.maj_ref_seq = maj_ref_seq

