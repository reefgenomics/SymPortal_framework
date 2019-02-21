import os
import shutil
from multiprocessing import Queue, Manager, Process
import sys
from dbApp.models import AnalysisType

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



    def analyse_data(self):
        print('Beginning profile discovery')
        self._populate_clade_fp_dicts_list()

        for clade_fp_dict in self.clade_footp_dicts_list:
            if self._there_are_footprints_of_this_clade(clade_fp_dict):  # If there are some clade collections for the given clade


                # FIND WHICH FOOTPRINTS ARE SUPPORTED AND WHICH CCs SUPPORT THEM #######

                collapsed_footprint_dict = self.collapse_potential_profiles_initial_type_objects(
                    footprint_list=clade_fp_dict,
                    reqsupport=4,
                    nprocessors=self.parent.args.num_proc)

                # CREATE ANALYSIS TYPES BASED ON DISCOVRED FOOTPRINTS
                # 08/12/17 we need to be careful here when we initiate the types as the types we were previously
                ''' generating would represent essentially the majoirty of the ccts sequences. but now some of the types
                will be smaller proportions of the ccts so we should check to see how the initial abundnace of the types
                are calculated. e.g. are they calculated as the proportion of the total seqs in the cct or are we 
                already working on as proportions of the seqs of the type in question. Hopefully it is the latter
                and we were just going with the types that represented the largest number of sequences for the cct.
                '''
                # 08/12/17 I have carefully looked through the type initTypeAtributes method
                ''' Firsly it always works in the context of the sequences found in the type. It produces absoulte
                 counts per sequence in the type for each cladeCollection that the type was supported by.
                 It also produces a count that is relative proportions of each sequence of the type for 
                 each clade_collection_object.
                  Hopefully this is what we are using when we do the second artefact check. I.e. we are looking
                  for the types in the CCs again.
                  For each type, the abslute counts per type sequence per clade_collection_object are stored in 
                  type.footprint_sequence_abundances
                  in the order of initalCCs and orderedfootprint list. THere is also the relative version wihich is 
                  stored as type.footprint_sequence_ratios'''
                # for every footprint that will become an analysis_type


                print(f'\n\nCreating analysis types clade {self.parent.clade_list[self.clade_footp_dicts_list.index(clade_fp_dict)]}')
                for initialType in collapsed_footprint_dict:

                    # Work out the corresponding reference_sequence for each Maj
                    # of the samples with that corresponding type
                    # Then do len(set()) and see if it is a co_dominant, i.e. different Maj seqs within the type


                    if len(initialType.set_of_maj_ref_seqs) > 1:  # Then this is a co_dominant

                        # the Counter class (from collections import Counter) may be useful
                        # http://stackoverflow.com/questions/2600191/how-can-i-count-the-occurrences-of-a-list-item-in-python
                        new_analysis_type = AnalysisType(
                            co_dominant=True,
                            data_analysis_from=analysis_object,
                            clade=clade_list[master_cladal_list_of_footprint_dicts.index(clade_fp_dict)])

                        new_analysis_type.set_maj_ref_seq_set(initialType.set_of_maj_ref_seqs)
                        new_analysis_type.init_type_attributes(initialType.clade_collection_list, initialType.profile)

                        new_analysis_type.save()
                        print('\rCreating analysis type: {}'.format(new_analysis_type.name), end='')
                    else:

                        new_analysis_type = AnalysisType(co_dominant=False, data_analysis_from=analysis_object,
                                                         clade=clade_list[
                                                              master_cladal_list_of_footprint_dicts.index(
                                                                  clade_fp_dict)])
                        new_analysis_type.set_maj_ref_seq_set(initialType.set_of_maj_ref_seqs)
                        new_analysis_type.init_type_attributes(initialType.clade_collection_list, initialType.profile)
                        new_analysis_type.save()
                        sys.stdout.write(f'\rCreating analysis type: {new_analysis_type.name}')




    def collapse_potential_profiles_initial_type_objects(self, footprint_list, reqsupport, nprocessors):
        # 10/12/17 I have just thought of something else we will have to consider. When you have a profile that
        # has a genuine codom type so there are e.g. MAJs of C3 and C3a, what will happen in this cricumstances.
        # For the time being we will not think too much about this and in the initialType profile classes that we are
        # creating, if we find a multiple basal type situation we will be sure to add the maj type to each of the
        # dsss lists.
        # 10/12/17, I also think we can take the opportunity to move away from working with dataSetSampleSequences
        # and move towards working with referenceSequences, as it seems that we are always converting them to
        # reference sequences any way. Maybe for the time being we will work with dsss and get things up and running
        # before converting to working with referenceSequences.
        # 09/12/17 I want to make sure that the subsets collapsing is still working and taking into account the
        # basal types compatabilites. Actually the subset collapsing is not taking into account basal sequences
        # but this is not a problem. We are essentially just consolidating the larger types into the smaller ones
        # the smaller ones still then have to found as supported using the basal sequences information.

        # 09/12/17 This has become unbarably slow. We have an awful lot of redundancy in here that I would like to fix
        # especially when it comes to checking for basal types. So I am going to create some clases which will only be
        # held in RAM. ie. not in the database.
        # class will be initialType. will contain, profile, profile length, contains_multiple_basal_sequences,
        # basal_sequence_list, clade_collection_list, maj_seq_list, maj_seq_set
        """This will be our latest battle field 07/12/17. We have a problem where we are getting multiple types
        So the object that we're playing around with here is a dictionary with key that is the set of referencesequences
        objects with a value of a 2d list with containing two lists. the first list is the clade collection object that
        this profile was found in. THe second list is the majority sequence of this set in this clade collection
        so if we want to be able to split up these objects to extract e.g. c15 and c3 sets then we'll need to create another
        dictionary entry. So we were probably relying on their only ever being one cladecollection to every set of sequences
        but if we start to do the extractions then were going to end up with multiple clade collections in the dict
        we will also need to be able to calcualate a new maj sequence for each for each of the extracted set of sequences
        we'll then have to see how this might reek havoc further down stream.
        so we end up creating the analysis types from the information that we have created here so
        I don't think that its going to be a problem that we ahve serveral types for a cgiven clade collection as it simply
        means that one clade collection is going to be found listed in multiple types. So lets give it a go!"""

        # We were having problems with the original.
        # The problem was when you have lots of long footprints that have intras in common i.e. D1 and D1a,
        # these long footprints are collapsed all the way down to D1 or D1a, depending on the maj becuase
        # there aren't any footprints that are JUST D1/D1a. So I have a better approach hopefully that will solve this
        # This approach builds on the original collapseFootprints
        # It starts with an n of the longest type and it tries to collapse each of the types into the n-1 types.
        # However, unlike the original, if it fails to collapse into the n-1 it does not then try to collapse into
        # n-2 types etc. Once we have tested all of the n length types and collapsed those that can be we move those that
        # couldn't into a list of unsupported types.
        # We then get all of the sequences found in the footprints of length n
        # and we do itertools to get all n-1 combinations
        # we then try to fit each of these into the unsupported types. We create
        # new footprints in order of the most supported
        # footprints and assign the associted CCs and Maj seqs. Any n length seqs that still don't have a collapse are
        # kept in the unsupported list. This is iterated dropping through the ns. I think this will be great!
        # 1 - check which supported which not, move supported to supportedlist, move unsupported to unsupported lit
        # 2 - for all unsupported check to see if can be collapsed into fooprints of n-1
        #   if can then collapse into and remove type from not supported list
        # 3 - get all n-1 combos of all seqs found in the n-1 unsupported sequences
        # All the seqs found in a clade_collection_object can now be accessed from the [2] list (third list) of
        # the foorprintList dict value
        # actually I don't think we'll need to caryy around the [2] list as we can simply use the footprint value
        # to get the list of refseqs - duh!
        # see how many footprints have these types, sort by order
        # create new footprints fot these footprints
        # if supported move to supportedlist
        # if not supported move to the non-supporetd list
        # next n...

        # Lists that will hold the footprints that have and do not have support throughout the collapse
        supported_list = []
        unsupported_list = []

        # convert the footprint_list to initalTypes
        initial_types_list = []
        for fpkey, fpvalue in footprint_list.items():
            initial_types_list.append(InitialType(fpkey, fpvalue[0], fpvalue[1]))

        # start with n at the largest footprint in the footprintlist
        # the lowest n that we will work with is 2
        # when we get to the footprints that are two in length we will do the first type of collapsing
        # i.e. we will collapse into existant 1s if feasible but we will not do the combinations
        # section as this will be moot as each will only be one intra big anyway.
        # once we have done the inital collapse of the 2s then we will sort the 1s into
        # supported or unsupported.
        # At this point all of the footprints that were in the footprint_list to begin with will be either in
        # the supported or unsupported list
        # we will then convert these types to the maj types

        # for each length starting at max and dropping by 1 with each increment
        for n in range(max([initT.profile_length for initT in initial_types_list]), 0, -1):
            # This dict will hold the top collapses where bigfootprint = key and smallFootprint = value

            # populate supported and unsupported list for the next n
            n_list = [initT for initT in initial_types_list if initT.profile_length == n]

            for initT in n_list:

                if initT.support >= reqsupport and not initT.contains_multiple_basal_sequences:  # supported
                    supported_list.append(initT)
                else:  # above support threshold
                    unsupported_list.append(initT)

            # TRY TO COLLAPSE SIZE n FOOTPRINTS INTO SIZE n-1 FOOTPRINTS
            # we will try iterating this as now that we have the potential to find two types in one profile, e.g.
            # a C15 and C3, we may only extract the C3 on the first one but there may still be a C15 in there.
            # whichever had the highest score will have been extracted. Eitherway, along the way, the unsupported_list
            # will have been updated
            repeat = True
            while repeat:
                collapse_dict = {}
                repeat = False
                n_minus_one_list = [initT for initT in initial_types_list if initT.profile_length == n - 1]
                if n_minus_one_list:

                    for bigFootprint in unsupported_list:  # For each big check if small will fit in

                        print('Assessing discovered footprint {} for supported type'.format(
                            '-'.join(str(refseq) for refseq in bigFootprint.profile)), end='\r')

                        top_score = 0
                        for smallerFootprint in n_minus_one_list:
                            # If the small foot print is a subset of the big footprint consider for collapse
                            # Only collapse into the small footprint if it doesn't contain multiple basal types.
                            # Only collapse if this is the best option i.e. if it give the largest number of support
                            multi_basal = smallerFootprint.contains_multiple_basal_sequences
                            if smallerFootprint.profile.issubset(bigFootprint.profile) and not multi_basal:
                                # Consider this for collapse only if the majsequences of the
                                # smaller are a subset of the maj sequences of the larger
                                # e.g. we don't want A B C D being collapsed into B C D when
                                # A is the maj of the first and B is a maj of the second

                                # simplest way to check this is to take the setOfMajRefSeqsLarge
                                # which is a set of all of the
                                # ref seqs that are majs in the cc that the footprint is found in
                                # and make sure that it is a subset of the smaller footprint in question

                                # 10/01/18 what we actually need is quite complicated. If the big type is not multi basal,
                                # then we have no problem and we need to find all of the set of
                                # maj ref seqs in the small profile
                                # but if the large type is multi basal then it gets a little more complicated
                                # if the large type is multi basal then which of its set of maj ref
                                # seqs we need to find in the small profile
                                # is dependent on what the basal seq of the smallfootprint is.
                                # if the small has no basal seqs in it, then we need to find every
                                # sequence in the large's set of maj ref seqs
                                # that is NOT a C15x, or the C3 or C1 sequences.
                                # if small basal = C15x then we need to find every one of the
                                # large's seqs that isn't C1 or C3
                                # if small basal = C1 then we need to find every on of the
                                # large's seqs that isn't C3 or C15x
                                # if small basal = C3 then we need to find every on of the
                                # large's seqs that isn't C1 or C15x
                                # we should put this decision into a new function

                                if bigFootprint.contains_multiple_basal_sequences:
                                    if does_small_footprint_contain_the_required_ref_seqs_of_the_large_footprint(
                                            bigFootprint, smallerFootprint):
                                        # score = number of samples big was found in plus num samples small was found in
                                        score = bigFootprint.support + smallerFootprint.support
                                        if score > top_score:
                                            top_score = score
                                            repeat = True
                                            collapse_dict[bigFootprint] = smallerFootprint

                                else:
                                    if bigFootprint.set_of_maj_ref_seqs.issubset(smallerFootprint.profile):
                                        # score = number of samples big was found in plus num samples small was found in
                                        score = bigFootprint.support + smallerFootprint.support
                                        if score > top_score:
                                            top_score = score
                                            repeat = True
                                            collapse_dict[bigFootprint] = smallerFootprint

                # Once here we have tried to collapse all of the large seq into all of the smaller footprinst
                # now if there is a collapse do the collapse and then remove it from the unsupported list
                # (we will see if the smaller footprint we added it to is supported next round of n)
                # Or each key of the collapse_dict add it to the value.
                # then remove the key from the unsupported_list
                # 08/12/17 here is where we will need to start to implement extraction
                # rather than just deletion for the potentially multiple
                # We will need to be careful that we don't start extracting profiles for non basal mixes, e.g. if we have
                # C3, C3a, C3b, C3d in the big and the small is C3, C3a, C3b we don't want to extract this, we want to use
                # the preivous method of deletion. I guess the best way to tell whether we want to do extraction vs deletion
                # is if there is another basal maj type in the large footprint
                # I'm starting to think that it might be worth having an initial type object where we can store useful
                # information like list of cc, list of dsss and list of basal maj
                # Yep, we now have this and I am doing the initial write now.

                collapse_dict_keys_list = list(collapse_dict.keys())
                # 070218 - Go through to make sure that all of the keys are found in the unsupported typelist-they should be

                count = 0
                for t in range(len(collapse_dict_keys_list)):
                    if collapse_dict_keys_list[t] not in unsupported_list:
                        count += 1

                # The collapse_dict footprint is now key=large initial type and value = small initial type

                for q in range(len(collapse_dict_keys_list)):
                    if n == 3:
                        count = 0
                        for t in range(len(collapse_dict_keys_list)):
                            if collapse_dict_keys_list[t] not in unsupported_list:
                                count += 1

                    # returns bool representing whether to extract. If false, delete rather than extract

                    if collapse_dict_keys_list[q].support >= \
                            reqsupport and not collapse_dict_keys_list[q].contains_multiple_basal_sequences:
                        # Then this type has already had some other leftovers put into it so that it now has the required
                        # support. In this case we can remove the type from the unsupported list and continue to the next
                        unsupported_list.remove(collapse_dict_keys_list[q])
                        continue

                    extraction_deletion = collapse_dict_keys_list[q].contains_multiple_basal_sequences
                    if not extraction_deletion:
                        # If large type does not contain multiple basal sequences then we do not need to extract and we
                        # can collapse the large into the small.
                        # we need to simply extend the clade collection list of the large type with that of the small
                        # we also need to add the dss lists of the large to the small as well

                        small_init_type = collapse_dict[collapse_dict_keys_list[q]]
                        small_init_type.absorb_large_init_type(collapse_dict_keys_list[q])

                        # # below = smaller footprint data = [[smallerCCs += bigger CCs],
                        # [smallerMajSeqsForCCs += largerMajSeqsForCCs]]
                        # footprint_list[collapse_dict[footprintToCollapse]] =
                        # [footprint_list[collapse_dict[footprintToCollapse]][0] +
                        # footprint_list[footprintToCollapse][0], footprint_list[collapse_dict[footprintToCollapse]][1] +
                        # footprint_list[footprintToCollapse][1]]

                        # remove the large_init_type from the init_type_list and from the unsupported_list
                        initial_types_list.remove(collapse_dict_keys_list[q])
                        unsupported_list.remove(collapse_dict_keys_list[q])
                    else:
                        # Here we need to implement the extraction. It is important that we extract all of the Majs that
                        # should go into the type we are collapsing into
                        # It is not just enough to send over all of the Majs that the big and
                        # small majref seqs have in common.
                        # Actually I think it is OK to just send over the Majs that the small
                        # and big have in commmon. We need to
                        # remember that we are simply listing ccts for support and the Majs
                        # found in those CCts. We can remove
                        # the subsetted sequences from the footprint of the big and put it back in the dict

                        # 171217 consolidating thoughts
                        # so below we want to extend the clade collection list of the small with that of the large
                        # this is not a problem to do and very straight forwards
                        # what will cause us more of an issue is the dsss, that we extend with
                        # for the dssss we want to send over the maj dss list of each of the large CCs if the dss is of
                        # a refseq that is found in the small init types profile or
                        # should all the ref seq basal types of the large be found in the small.
                        # My thoughts below seem to think that's its ok to just add the maj dsss that are of a refseq
                        # that are found in common between the two types

                        # 1 - Collapse big to small
                        # This will all be done within the extract_support_from_large_initial_type method
                        small_init_type = collapse_dict[collapse_dict_keys_list[q]]
                        small_init_type.extract_support_from_large_initial_type(collapse_dict_keys_list[q])

                        # the large init type should remain in the init_type_list

                        # need to check if the new profile created from the original footprintToCollapse
                        #  that has now had the small)init_type extracted from it is already shared with another init_type
                        # If so then we need to combine the init_types.
                        # else then we don't need to do anything
                        # If the similar type is also in the collapse dict then we will collapse to that type
                        # else if the type is not in the collapse dict then we will absorb that type.
                        # There should only be a maximum of one initT that has the same
                        match = False
                        for p in range(len(initial_types_list)):
                            if initial_types_list[p].profile == \
                                    collapse_dict_keys_list[q].profile and initial_types_list[p] != \
                                    collapse_dict_keys_list[q]:

                                # Then we have found an initT that has an exact match for the large initT in Q
                                # In this case we need to create a single initType that is a combination of the two
                                # 070218 we are going to change this and we are going to extract
                                # Let's add the existing initT to the footPinrtToCollapse
                                # This way we can make sure that the footprintToCoolapse is still in the correct place
                                # i.e. in the unsupported list or not.

                                # If we do find a match then we need to make sure to get rid of the initT that we have
                                # absorbed into the footprintToCollapse

                                # Should this just be the same as when a small initT absorbs a large initT?
                                # I think so but lets check
                                # check to see that this is appropriate
                                # we need to check specifically if the initial_types[p] is found in the types that
                                # still need to be collapsed. so we need to slice here.
                                match = True
                                if initial_types_list[p] in collapse_dict_keys_list[q + 1:]:
                                    # then we need to collapse into the found type
                                    # because this collapsing can cause the matching type that is also in the collapse
                                    # dict to gain support we will also check to if each of the types in the collapse dict
                                    # have gained sufficient support to no longer need collapsing. We will do this earlier
                                    # in the process, not here.
                                    initial_types_list[p].absorb_large_init_type(collapse_dict_keys_list[q])
                                    unsupported_list.remove(collapse_dict_keys_list[q])
                                    initial_types_list.remove(collapse_dict_keys_list[q])
                                    break
                                else:
                                    collapse_dict_keys_list[q].absorb_large_init_type(initial_types_list[p])

                                    # The initT will also need removing from the unsupported_list if it is in there
                                    # This is causing a problem
                                    # Some time the initial_types_list[p] that matches the leftover
                                    # collapse_dict_keys_list[q]
                                    # can also be in the collapse_dict_keys_list.
                                    # In this case it will be removed from the unsupported_list at
                                    # this point and won't be removable
                                    # later on.
                                    # Shouldn't we be checking to see if initial_types_list[p] is smaller than n.
                                    if initial_types_list[p] in unsupported_list:
                                        unsupported_list.remove(initial_types_list[p])

                                    # Delete the initT as this has now been absorbed into the footprint to collapse
                                    initial_types_list.remove(initial_types_list[p])
                                    # If the left over type is less than n then we need to now remove it from the un
                                    # supported list as it will be collapsed on another iteration than this one.
                                    if collapse_dict_keys_list[q].profile_length < n:

                                        unsupported_list.remove(collapse_dict_keys_list[q])

                                    else:
                                        # now we need to check to see if the
                                        # collapse_dict_keys_list[q] type has support bigger then
                                        # the required. If it does, then it should also be removed from the unsupported list

                                        if collapse_dict_keys_list[q].support >= \
                                                reqsupport \
                                                and not collapse_dict_keys_list[q].contains_multiple_basal_sequences:
                                            unsupported_list.remove(collapse_dict_keys_list[q])
                                        else:
                                            # if it doesn't have support then we simply leave it in the unsupportedlist
                                            # and it will go on to be seen if it can be collapsed into one of the insilico
                                            # types that are genearted.
                                            pass
                                    break
                        if not match:
                            if collapse_dict_keys_list[q].profile_length < n:
                                unsupported_list.remove(collapse_dict_keys_list[q])

                        # the large init_type does not need removing from the initial type list.
                        # but we still need to do the check to see if the profile length is SMALLER than n. If it is smaller
                        # then n then we should remove it from the unsupported list else we should leave it in the
                        # unsupported list.

            if n > 2:

                # We only need to attempt the further collapsing if there are unsupported types to collapse
                # else move onto the next n
                if len(unsupported_list) > 1:
                    # Now we are left with the footprints that are still in the unsupported list

                    # we are going to change this.
                    # Currently we are generating an enormous generator in the itertools.combinations
                    # This is because we are using all of the sequences found in all of the length n types
                    # However, I see non reason why we should be using
                    # combinatations of sequences found in different profiles.
                    # This is not in keeping with the biology.
                    # Instead the synthetic types should be made from
                    # combinattions of sequences found in each profile
                    # This will better represent the biology and also greatly cutdown on computational cost
                    # And it will mean that we can multiprocess this

                    list_of_types_of_length_n = [initT for initT in initial_types_list if initT.profile_length == n]
                    # Only carry on if we have lengthN footprints to get sequences from
                    if list_of_types_of_length_n:

                        # here is where we should start the new approach.
                        # add the nLength Types into a MP list
                        # then add the stops
                        # then in the worker create a dict that is same as the one below,
                        # only the using the nLengthtype instead
                        # Create the queues that will hold the sample information
                        print('\rstarting to generate Synthetic footprints for collapse', end='')

                        task_queue = Queue()
                        collapse_n_mer_dictionary_manager = Manager()
                        collapse_n_mer_dictionary = collapse_n_mer_dictionary_manager.dict()
                        # output_queue = Queue()

                        for nLengthType in list_of_types_of_length_n:
                            task_queue.put(nLengthType.profile)

                        for N in range(nprocessors):
                            task_queue.put('STOP')

                        all_processes = []

                        for N in range(nprocessors):
                            p = Process(target=worker_discovery_one, args=(task_queue, collapse_n_mer_dictionary, n))
                            all_processes.append(p)
                            p.start()

                        for p in all_processes:
                            p.join()
                        # The manager(dict) object was not behaving correctly when I was trying to append items to the list
                        # values. However, when converting to a dict, it does.
                        collapse_n_mer_dictionary = dict(collapse_n_mer_dictionary)

                        # Now go through each of the (n-1)mer footprints and see if they
                        # fit into a footprint in the unsuported list
                        # This dict value = the synthetic footprint of length n-1,
                        # value, a list of unsupported types that are mastersets of the synthetic footprint

                        print('\rGenerated {} Synthetic footprints'.format(len(collapse_n_mer_dictionary)), end='')
                        # We create a dict that is key, type and value is the set of the majrefseqs found in the type
                        # this was being calculated once for every synthetic type.
                        # But the number of synthetic types are getting very large
                        # so by having this dict we should speed things up considerably
                        # This info is used to make sure that when collapsing types the maj sequences are shared between the
                        # larger and smaller types

                        # Items for printing out our progress
                        total = len(collapse_n_mer_dictionary) * len(unsupported_list)
                        count = 0
                        print_value = 0.01
                        print('\rChecking new set of synthetic types', end='')

                        for frTupKey in collapse_n_mer_dictionary.keys():  # for each of the synthetic footprints
                            # 10/01/18 we need to check that each of the frTupKeys contains multiple basal seqs
                            # If it does then we can't collapse into it.
                            if does_set_of_ref_seqs_contain_multiple_basal_types(frTupKey):
                                continue
                            for nLengthType in unsupported_list:  # for each of the unsupported init_types
                                # Items for printing out progress
                                count += 1
                                percent_complete = count / total
                                if percent_complete > print_value:
                                    print("\r%.2f" % percent_complete, end='')
                                    print_value = max([print_value + 0.01, percent_complete + 0.01])

                                # For each of the N-1mers see if they fit within the unsupported types
                                # if so then add the footprint into the n-1mers list
                                # we should check for the maj rule again. I.e. we should want all of the Majs of the
                                # footprint in Q to be in the kmer footprint
                                # maybe recycle code for this
                                if frTupKey.issubset(nLengthType.profile):
                                    # Now check to see that at least one of the set_of_maj_ref_seqs seqs is found in the
                                    # frTupKey.
                                    if len(nLengthType.set_of_maj_ref_seqs & frTupKey) > 0:
                                        # Then this nLengthType init_type can be collapsed into the frTupKey
                                        collapse_n_mer_dictionary[frTupKey].append(nLengthType)
                        # Here we have a populated dict.
                        # Order dict
                        # We are only interseted in kmers that were found in the unsupported types
                        # Because alot of the kmers will be found in the sequences they originated from
                        # require the kmer to be associated with more than 1 cladecollection
                        # I am not able to express my logic here perfectly but I am sure that this makes sense

                        list_of_n_kmers_with_collapsable_footprints = [kmer for kmer in collapse_n_mer_dictionary.keys()
                                                                       if
                                                                       len(collapse_n_mer_dictionary[kmer]) > 1]
                        # No need to continue if there are no footprints that match the nKmers
                        if list_of_n_kmers_with_collapsable_footprints:
                            # We now go through the kmers according to the number of footprints they were found in
                            ordered_list_of_populated_kmers = sorted(list_of_n_kmers_with_collapsable_footprints,
                                                                     key=lambda x: len(collapse_n_mer_dictionary[x]),
                                                                     reverse=True)

                            for kmer in ordered_list_of_populated_kmers:
                                # for each initType in the kmer lists

                                for k in range(len(collapse_n_mer_dictionary[kmer])):
                                    # chek to see if the footprint is still in the unsupportedlist
                                    # !! 11/01/18 also need to check that the bigFootprint hasn't already been collapsed
                                    # if it isn't then it has already been associated to a new footprint
                                    #   pass it over
                                    if collapse_n_mer_dictionary[kmer][k] in unsupported_list and kmer.issubset(
                                            collapse_n_mer_dictionary[kmer][
                                                k].profile):  # Then this footprint hasn't been collapsed anywhere yet.
                                        # Here we also check to make sure that the profile hasn't been changed
                                        # if an init_type already exists with the profile (kmer)
                                        # in question then add the big init_type to it
                                        # 231217 we will have to check whether the big init_type is a multiple basal seqs
                                        # and therefore whether this is an extraction or a absorption
                                        exists = False
                                        for i in range(len(initial_types_list)):
                                            # for initT_one in initial_types_list:
                                            if initial_types_list[i].profile == kmer:
                                                exists = True
                                                # 231217 we will have to check whether the big
                                                # init_type is a multiple basal seqs
                                                # and therefore whether this is an extraction or a absorption
                                                # bear in mind that it doesn't matter if the initT we are absorbing or
                                                # extracting into is a multi basal. We will wory about that in the next
                                                # iteration
                                                if collapse_n_mer_dictionary[kmer][k].contains_multiple_basal_sequences:
                                                    # 231317 this still needs writing
                                                    # Then we need to extract
                                                    initial_types_list[i].extract_support_from_large_initial_type(
                                                        collapse_n_mer_dictionary[kmer][k])

                                                    # Once we have extracted into the smaller type
                                                    # check to see if the big init Type still contains set_of_maj_ref_seqs
                                                    if collapse_n_mer_dictionary[kmer][k].set_of_maj_ref_seqs:
                                                        # 10/01/18 we first need to check if the new profile already
                                                        # exists and if it does we need to do as we do above
                                                        # and add it to the similar one
                                                        # if the big init type still exists with maj containing profiles
                                                        # then remove from unsupported list
                                                        for j in range(len(initial_types_list)):
                                                            # for initT_two in initial_types_list:

                                                            if initial_types_list[j].profile == \
                                                                    collapse_n_mer_dictionary[kmer][
                                                                        k].profile and initial_types_list[j] != \
                                                                    collapse_n_mer_dictionary[kmer][k]:

                                                                # Then we have found an initT that has an exact match for
                                                                # the large initT in Q
                                                                # In this case we need to create a single initType that
                                                                # is a combination of the two
                                                                # Let's add the existing initT to the footPinrtToCollapse
                                                                # This way we can make sure that the footprintToCoolapse
                                                                # is still in the correct place
                                                                # i.e. in the unsupported list or not.

                                                                # If we do find a match then we need to make sure to
                                                                # get rid of the initT that we have
                                                                # absorbed into the footprintToCollapse

                                                                # Should this just be the same as when a small initT
                                                                #  absorbs a large initT?
                                                                # I think so but lets check
                                                                # check to see that this is appropriate
                                                                collapse_n_mer_dictionary[kmer][
                                                                    k].absorb_large_init_type(
                                                                    initial_types_list[j])

                                                                # the initT will need removing from the inital_types_list
                                                                # and the unsupported_list as it no longer exists.
                                                                if initial_types_list[j] in unsupported_list:
                                                                    unsupported_list.remove(initial_types_list[j])
                                                                initial_types_list.remove(initial_types_list[j])
                                                                break

                                                        # We now need to decide if the footprint to collapse should be
                                                        # removed from the unsupported_list.
                                                        # This will depend on if it is longer then n or not.
                                                        if collapse_n_mer_dictionary[kmer][k].profile_length < n:
                                                            unsupported_list.remove(collapse_n_mer_dictionary[kmer][k])

                                                    else:
                                                        # then the big init_type no longer contains any
                                                        # of its original maj ref seqs and
                                                        # so should be delted.
                                                        initial_types_list.remove(collapse_n_mer_dictionary[kmer][k])
                                                        unsupported_list.remove(collapse_n_mer_dictionary[kmer][k])
                                                else:
                                                    # 231217 then we need can absorb
                                                    initial_types_list[i].absorb_large_init_type(
                                                        collapse_n_mer_dictionary[kmer][k])
                                                    # make sure that we then get rid of the bigFootprintToCollapse
                                                    initial_types_list.remove(collapse_n_mer_dictionary[kmer][k])
                                                    unsupported_list.remove(collapse_n_mer_dictionary[kmer][k])
                                                break
                                        if not exists:  # then the kmer was not already represented by an existing initT
                                            # 270118 We need to check to see if the big type is contains mutiple
                                            # basal if it does then we should extract as above. This should be exactly the
                                            # same code as above but extracting into a new type rather than an existing one
                                            # Try to create a blank type
                                            # NEW
                                            if collapse_n_mer_dictionary[kmer][k].contains_multiple_basal_sequences:
                                                new_blank_initial_type = InitialType(reference_sequence_set=kmer,
                                                                                     clade_collection_list=list(
                                                                                         collapse_n_mer_dictionary[
                                                                                             kmer][
                                                                                             k].clade_collection_list))
                                                initial_types_list.append(new_blank_initial_type)

                                                # now remove the above new type's worth of info
                                                # from the current big footprint
                                                collapse_n_mer_dictionary[kmer][
                                                    k].substract_init_type_from_other_init_type(
                                                    new_blank_initial_type)

                                                if collapse_n_mer_dictionary[kmer][k].set_of_maj_ref_seqs:
                                                    # we first need to check if the new profile already exists and if
                                                    # it does we need to do as we do above and add it to the similar one
                                                    # if the big init type is still exists with maj containing profiles
                                                    # then remove from unsupported list
                                                    for initT in initial_types_list:
                                                        if initT.profile == collapse_n_mer_dictionary[kmer][
                                                            k].profile and initT != \
                                                                collapse_n_mer_dictionary[kmer][k]:
                                                            # Then we have found an initT that has an exact
                                                            # match for the large initT in Q
                                                            # In this case we need to create a single initType
                                                            # that is a combination of the two
                                                            # Let's add the existing initT to the footPinrtToCollapse
                                                            # This way we can make sure that the footprintToCoolapse
                                                            # is still in the correct place
                                                            # i.e. in the unsupported list or not.

                                                            # If we do find a match then we need to make sure to get
                                                            # rid of the initT that we have
                                                            # absorbed into the footprintToCollapse

                                                            # Should this just be the same as when a small initT
                                                            # absorbs a large initT?
                                                            # I think so but lets check
                                                            # check to see that this is appropriate
                                                            collapse_n_mer_dictionary[kmer][k].absorb_large_init_type(
                                                                initT)

                                                            # the initT will need removing from the inital_types_list and
                                                            # the unsupported_list as it no longer exists.
                                                            initial_types_list.remove(initT)
                                                            if initT in unsupported_list:
                                                                unsupported_list.remove(initT)
                                                            break

                                                    # We now need to decide if the footprint to collapse should be
                                                    # removed from the unsupported_list.
                                                    # This will depend on if it is longer then n or not.
                                                    if collapse_n_mer_dictionary[kmer][k].profile_length < n:
                                                        unsupported_list.remove(collapse_n_mer_dictionary[kmer][k])
                                                        # unsupported_list.remove(bigFootprintToCollapse)

                                                else:
                                                    # then the big init_type no longer contains any of
                                                    # its original maj ref seqs and
                                                    # so should be delted.
                                                    initial_types_list.remove(collapse_n_mer_dictionary[kmer][k])
                                                    unsupported_list.remove(collapse_n_mer_dictionary[kmer][k])
                                            # NEW
                                            # 231217 in this case we should create a new inital_type
                                            # I am going to have to think about how we will do this
                                            # this still needs to be coded

                                            # the profile will be the kmer
                                            # the cladeCollectionList will be the same as the large init type
                                            # the basalSequence_list can come from the
                                            # self.check_if_initial_type_contains_basal_sequences()

                                            # NB i have made it so that if InitialType() doesn't get a
                                            # majdsss list then it will create one from scratch
                                            # ALSO CHANGED
                                            else:
                                                new_initial_type = InitialType(
                                                    reference_sequence_set=kmer,
                                                    clade_collection_list=collapse_n_mer_dictionary[kmer]
                                                    [k].clade_collection_list)
                                                initial_types_list.append(new_initial_type)

                                                # now delete the big init_type from the intial
                                                # types list and from the unsupported_list
                                                initial_types_list.remove(collapse_n_mer_dictionary[kmer][k])
                                                if collapse_n_mer_dictionary[kmer][k] in unsupported_list:
                                                    unsupported_list.remove(collapse_n_mer_dictionary[kmer][k])

                                # each unsupportedfootprint in the nmer has been collapsed
                            # each of the kmers have been gone through and all unsupportedfootprints
                            # that could fit into the kmers have been collapsed
                            # Time for the next n here!
                            # Check to see what is happening once we get to the smaller ns
                            # It may be that we can use the code from below to simply put the maj as
                            # the type for footprints that are still
                            # unassigned
            else:
                # At this point we may have some 2s or larger in the unsupported list still
                # we will collapse these to their maj ref seq type
                # Once we have pushed these into the 1s or made new 1s if the maj ref seq type did not already exist
                # then we should collect all of the 1 length footprints and put them in the supported list
                # Or actually in theory so long as we have deleted collapsed footprints from the footprint_list we should
                # be able to return the footprint list
                while unsupported_list:
                    unsupported_footprint = unsupported_list[0]

                    # For each clade collection
                    # 241217 we need to do each cladeCollection individually here bear in mind.
                    for i in range(len(unsupported_footprint.clade_collection_list)):
                        # check each maj dsss
                        for maj_dss in unsupported_footprint.majority_sequence_list[i]:
                            # Check to see if an initial type with profile of that maj_dss refseq already exists
                            # 241217 it may be worth having a dictionary of the profiles and init types to speed this up.
                            # Let's see how slow it is.
                            found = False
                            for initT in [init for init in initial_types_list if init.profile_length == 1]:
                                if maj_dss.reference_sequence_of in initT.profile:
                                    # Then we have found a profile that has the refseq as a profile
                                    # If the big init type is multi basal then we need to extract
                                    # if the clade maj dss contains multiple sequences
                                    # if unsupported_footprint.contains_multiple_basal_sequences:
                                    #     #241217 check that the large type has had the basal seqs readjusted
                                    #     # so that when we do the next maj_dss it may be non multiple basal
                                    #     initT.extract_support_from_large_initial_type(unsupported_footprint)
                                    # # else we need to absorb
                                    # else:
                                    #     initT.absorb_large_init_type(unsupported_footprint)

                                    # If found then we simply need to add a clade collection to the initT support
                                    # and the appropriate maj dsss
                                    initT.clade_collection_list.append(unsupported_footprint.clade_collection_list[i])
                                    initT.majority_sequence_list.append([maj_dss])

                                    found = True
                                    break
                            if not found:
                                # If the type init_type doesn't already exist we must create one.
                                new_initial_type = InitialType(
                                    reference_sequence_set=frozenset([maj_dss.reference_sequence_of]),
                                    clade_collection_list=[unsupported_footprint.clade_collection_list[i]])
                                initial_types_list.append(new_initial_type)
                                # We must then alter the big_init_type
                                # Actually I don't think we do need to alter the big type
                                # We just put it into the small for each of the cladeCcols and maj dsss
                                # once this is done then we delte the big type

                                # unsupported_footprint.remove_small_init_type_from_large(new_initial_type)

                                # If the large init type no longer contains any maj ref seqs then we can delete it
                    # Here we have completed collapsing one big_init_type and we can now get rid of the type
                    unsupported_list.remove(unsupported_footprint)
                    initial_types_list.remove(unsupported_footprint)

                # All unsupported footprints have been associated to their maj and deleted from the footprint_list
                return initial_types_list

        return False

    def _there_are_footprints_of_this_clade(self, clade_fp_dict):
        return clade_fp_dict

    def _populate_clade_fp_dicts_list(self):
        footprint_dict_pop_handler = FootprintDictPopHandler(sp_data_analysis_parent=self)
        footprint_dict_pop_handler.populate_clade_footprint_dicts()

    def _setup_temp_wkd(self):
        if os.path.exists(self.temp_wkd):
            shutil.rmtree(self.temp_wkd)
        os.makedirs(self.temp_wkd, exist_ok=True)

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
    def __init__(self, clade_footprint_dict):
        self.clade_fp_dict = clade_footprint_dict
        self.supported_list = []
        self.unsupported_list = []
        self.initial_types_list = []
        self._init_initial_types_list()
        # number of clade collections a footprint must be found in to be supported
        self.required_support = 4

    def _init_initial_types_list(self):
        for footprint_key, footprint_representative in self.clade_fp_dict.items():
            self.initial_types_list.append(
                InitialType(
                    footprint_key, footprint_representative.cc_list, footprint_representative.maj_seq_list))

    def identify_supported_footprints(self):
        # for each length starting at max and dropping by 1 with each increment
        longest_footprint = self._return_len_of_longest_fp()
        for n in range(longest_footprint, 0, -1):
            self._update_supported_unsupported_lists_for_n(n)

            # Try to collapse length n footprints into size n-1 footprints
            # we will try iterating this as now that we have the potential to find two types in one profile, e.g.
            # a C15 and C3, we may only extract the C3 on the first iteration but there may still be a C15 in initial
            # type.
            repeat = True
            while repeat:
                collapse_dict, repeat = self._populate_collapse_dict_for_next_n(n, repeat)
                # TODO you got here
                collapse_dict_keys_list = list(collapse_dict.keys())

    def _populate_collapse_dict_for_next_n(self, n, repeat):
        collapse_dict = {}
        repeat = False
        n_minus_one_list = self._get_n_minus_one_list(n)
        if n_minus_one_list:
            for longer_initial_type in self.unsupported_list:
                sys.stdout.write(f'Assessing footprint {longer_initial_type} for supported type\r')
                top_score = 0
                for shorter_initial_type in n_minus_one_list:
                    # If the small foot print is a subset of the big footprint consider for collapse
                    # Only collapse into the small footprint if it doesn't contain multiple basal types.
                    # Only collapse if this is the best option i.e. if it give the largest number of support
                    multi_basal = shorter_initial_type.contains_multiple_basal_sequences
                    if shorter_initial_type.profile.issubset(longer_initial_type.profile) and not multi_basal:
                        # Consider this for collapse only if the majsequences of the
                        # smaller are a subset of the maj sequences of the larger
                        # e.g. we don't want A B C D being collapsed into B C D when
                        # A is the maj of the first and B is a maj of the second

                        # simplest way to check this is to take the setOfMajRefSeqsLarge
                        # which is a set of all of the
                        # ref seqs that are majs in the cc that the footprint is found in
                        # and make sure that it is a subset of the smaller footprint in question

                        # 10/01/18 what we actually need is quite complicated. If the big type is not multi basal,
                        # then we have no problem and we need to find all of the set of
                        # maj ref seqs in the small profile
                        # but if the large type is multi basal then it gets a little more complicated
                        # if the large type is multi basal then which of its set of maj ref
                        # seqs we need to find in the small profile
                        # is dependent on what the basal seq of the smallfootprint is.
                        # if the small has no basal seqs in it, then we need to find every
                        # sequence in the large's set of maj ref seqs
                        # that is NOT a C15x, or the C3 or C1 sequences.
                        # if small basal = C15x then we need to find every one of the
                        # large's seqs that isn't C1 or C3
                        # if small basal = C1 then we need to find every on of the
                        # large's seqs that isn't C3 or C15x
                        # if small basal = C3 then we need to find every on of the
                        # large's seqs that isn't C1 or C15x
                        # we should put this decision into a new function

                        if longer_initial_type.contains_multiple_basal_sequences:
                            if self.does_small_footprint_contain_the_required_ref_seqs_of_the_large_footprint(
                                    longer_initial_type, shorter_initial_type):
                                # score = number of samples big was found in plus num samples small was found in
                                score = longer_initial_type.support + shorter_initial_type.support
                                if score > top_score:
                                    top_score = score
                                    repeat = True
                                    collapse_dict[longer_initial_type] = shorter_initial_type

                        else:
                            if longer_initial_type.set_of_maj_ref_seqs.issubset(shorter_initial_type.profile):
                                # score = number of samples big was found in plus num samples small was found in
                                score = longer_initial_type.support + shorter_initial_type.support
                                if score > top_score:
                                    top_score = score
                                    repeat = True
                                    collapse_dict[longer_initial_type] = shorter_initial_type
        return collapse_dict, repeat

    def does_small_footprint_contain_the_required_ref_seqs_of_the_large_footprint(self, longer_initial_type, shorter_initial_type):
        if shorter_initial_type.basalSequence_list:
            if 'C15' in shorter_initial_type.basalSequence_list[0]:
                # then this is a C15x basal type and we will need to find all sequences that are not C1 or C3
                set_of_seqs_to_find = set()
                ref_seqs_in_big_init_type = list(longer_initial_type.set_of_maj_ref_seqs)

                for ref_seq in ref_seqs_in_big_init_type:
                    if ref_seq.name in ['C1', 'C3']:
                        # then this is a squence we don't need to find
                        pass
                    else:
                        set_of_seqs_to_find.add(ref_seq)
                # Here we have the list of the ref_seqs that we need to find in the small_init_type.profile
                if set_of_seqs_to_find.issubset(shorter_initial_type.profile):
                    return True
                else:
                    return False
            elif shorter_initial_type.basalSequence_list[0] == 'C1':
                # then this is a C1 basal type and we need to find all sequence that are not C15x or C3
                set_of_seqs_to_find = set()
                ref_seqs_in_big_init_type = list(longer_initial_type.set_of_maj_ref_seqs)

                for ref_seq in ref_seqs_in_big_init_type:
                    if 'C15' in ref_seq.name or ref_seq.name == 'C3':
                        # then this is a squence we don't need to find
                        pass
                    else:
                        set_of_seqs_to_find.add(ref_seq)
                # Here we have the list of the ref_seqs that we need to find in the small_init_type.profile
                if set_of_seqs_to_find.issubset(shorter_initial_type.profile):
                    return True
                else:
                    return False
            elif shorter_initial_type.basalSequence_list[0] == 'C3':
                # then this is a C3 basal type and we need to find all sequence that are not C15x or C1
                set_of_seqs_to_find = set()
                ref_seqs_in_big_init_type = list(longer_initial_type.set_of_maj_ref_seqs)

                for ref_seq in ref_seqs_in_big_init_type:
                    if 'C15' in ref_seq.name or ref_seq.name == 'C1':
                        # then this is a squence we don't need to find
                        pass
                    else:
                        set_of_seqs_to_find.add(ref_seq)
                # Here we have the list of the ref_seqs that we need to find in the small_init_type.profile
                if set_of_seqs_to_find.issubset(shorter_initial_type.profile):
                    return True
                else:
                    return False
        else:
            # if the small_init_type doesn't contain a basal sequence sequence, then we need to find all of the seqs
            # in the big_intit_type.set_of_maj_ref_seqs that are not C15x, C1 or C3
            set_of_seqs_to_find = set()
            ref_seqs_in_big_init_type = list(longer_initial_type.set_of_maj_ref_seqs)

            for ref_seq in ref_seqs_in_big_init_type:
                if 'C15' in ref_seq.name or ref_seq.name in ['C1', 'C3']:
                    # then this is a squence we don't need to find
                    pass
                else:
                    set_of_seqs_to_find.add(ref_seq)
            # Here we have the list of the ref_seqs that we need to find in the small_init_type.profile
            if set_of_seqs_to_find.issubset(shorter_initial_type.profile):
                return True
            else:
                return False

    def _get_n_minus_one_list(self, n):
        n_minus_one_list = [initial_type for initial_type in self.initial_types_list if
                            initial_type.profile_length == n - 1]
        return n_minus_one_list

    def _return_len_of_longest_fp(self):
        longest_footprint = max([initial_type.profile_length for initial_type in self.initial_types_list])
        return longest_footprint

    def _update_supported_unsupported_lists_for_n(self, n):
        # populate supported and unsupported list for the next n
        n_list = [initial_type for initial_type in self.initial_types_list if initial_type.profile_length == n]

        for initial_type in n_list:
            if initial_type.support >= self.required_support and not initial_type.contains_multiple_basal_sequences:
                self.supported_list.append(initial_type)
            else:
                self.unsupported_list.append(initial_type)

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