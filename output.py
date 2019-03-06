from dbApp.models import (DataSet, ReferenceSequence, DataSetSampleSequence, AnalysisType, DataSetSample,
                          DataAnalysis, CladeCollection, CladeCollectionType)
from multiprocessing import Queue, Process, Manager
import sys
from django import db
from datetime import datetime
import os
import json
from general import write_list_to_destination
from collections import defaultdict
import pandas as pd
from plotting import generate_stacked_bar_data_analysis_type_profiles, generate_stacked_bar_data_loading
import pickle
from collections import Counter
import numpy as np
import sp_config
import plotting


def output_type_count_tables(
        analysisobj, datasubstooutput, call_type,
        num_processors=1, no_figures=False, output_user=None, time_date_str=None):
    analysis_object = analysisobj
    # This is one of the last things to do before we can use our first dataset
    # The table will have types as columns and rows as samples
    # its rows will give various data about each type as well as how much they were present in each sample

    # It will produce four output files, for DIV abundances and proportions and Type abundances and proportions.
    # found in the given clade collection.
    # Types will be listed first by clade and then by the number of clade collections the types were found in

    # Each type's ID will also be given as a UID.
    # The date the database was accessed and the version should also be noted
    # The formal species descriptions which correspond to the found ITS2 type will be noted for each type
    # Finally the AccessionNumber of each of the defining reference species will also be noted
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    # List of the paths to the files that have been output
    output_files_list = []

    data_submissions_to_output = [int(a) for a in str(datasubstooutput).split(',')]

    # Get collection of types that are specific for the dataSubmissions we are looking at
    query_set_of_data_sets = DataSet.objects.filter(id__in=data_submissions_to_output)

    clade_collections_from_data_sets = CladeCollection.objects.filter(
        data_set_sample_from__data_submission_from__in=query_set_of_data_sets)
    clade_collection_types_from_this_data_analysis_and_data_set = CladeCollectionType.objects.filter(
        clade_collection_found_in__in=clade_collections_from_data_sets,
        analysis_type_of__data_analysis_from=analysis_object).distinct()

    at = set()
    for cct in clade_collection_types_from_this_data_analysis_and_data_set:
        sys.stdout.write('\rCollecting analysis_type from clade collection type {}'.format(cct))
        at.add(cct.analysis_type_of)
    at = list(at)

    # Need to get a list of Samples that are part of the dataAnalysis
    list_of_data_set_samples = list(
        DataSetSample.objects.filter(
            data_submission_from__in=DataSet.objects.filter(id__in=data_submissions_to_output)))

    # Now go through the types, firstly by clade and then by the number of cladeCollections they were found in
    # Populate a 2D list of types with a list per clade
    types_cladal_list = [[] for _ in clade_list]
    for att in at:
        try:
            if len(att.list_of_clade_collections) > 0:
                types_cladal_list[clade_list.index(att.clade)].append(att)
        except:
            pass

    # Items for for creating the new samples sorted output
    across_clade_type_sample_abund_dict = dict()
    across_clade_sorted_type_order = []
    # Go clade by clade
    for i in range(len(types_cladal_list)):
        if types_cladal_list[i]:
            clade_in_question = clade_list[i]
            # ##### CALCULATE REL ABUND AND SD OF DEF INTRAS FOR THIS TYPE ###############
            #     # For each type we want to calculate the average proportions of the defining seqs in that type
            #     # We will also calculate an SD.
            #     # We will do both of these calculations using the footprint_sequence_abundances, list_of_clade_collections
            #     # and orderedFoorprintList attributes of each type

            # ######### MAKE GROUP COUNTER AND DICT ###########
            # These are now going to be managed items for use in the MP
            # We want to name the groups that types are found in sequencially
            # To get the next number to name a group we will use the groupCount
            # To look up what number has been assigned to groups that have
            # already been printed we will use the groupDict
            # groupCount = 0
            # groupDict = {}
            ###################################################

            # sort types by the number of samples they were found in for this output (not across whole analysis)
            # returns list of tuples, [0] = analysis_type object, [1] number of ccs found in for this output
            sorted_list_of_types = sort_list_of_types_by_clade_collections_in_current_output(
                types_cladal_list[i], clade_collection_types_from_this_data_analysis_and_data_set)

            # Here we will MP at the type level
            # i.e. each type will be processed on a differnt core. In order for this to work we should have a managed
            # dictionary where the key can be the types ID and the value can be the datalist
            # In order to get the types output in the correct order we should use the sortedListOfTypes to resort the
            # data once the MPing has been done.

            # This dict will hold all of the output rows that the MPing has created.
            worker_manager = Manager()
            type_output_managed_dict = worker_manager.dict({an_type: None for an_type in sorted_list_of_types})

            # NB using shared items was considerably slowing us down so now I will just use copies of items
            # we don't actually need to use managed items as they don't really need to be shared (i.e. input
            # from different processes.
            # a list that is the uids of each of the samples in the analyis
            sample_uids_list = [smp.id for smp in list_of_data_set_samples]
            # listOfDataSetSampleIDsManagedList = worker_manager.list(sample_uids_list)

            # A corresponding dictionary that is the list of the clade collection uids for each of the samples
            # that are in the ID list above.
            sample_uid_to_clade_collection_uids_of_clade = {}
            print('\nCreating sample_ID_to_cc_IDs dictionary clade {}'.format(clade_in_question))
            for smp in list_of_data_set_samples:
                sys.stdout.write('\rCollecting clade collections for sample {}'.format(smp))
                try:
                    sample_uid_to_clade_collection_uids_of_clade[smp.id] = CladeCollection.objects.get(
                        data_set_sample_from=smp, clade=clade_in_question).id
                except CladeCollection.DoesNotExist:
                    sample_uid_to_clade_collection_uids_of_clade[smp.id] = None
                except Exception as ex:  # Just incase there is some weird stuff going on
                    print(ex)

            type_input_queue = Queue()

            for an_type in sorted_list_of_types:
                type_input_queue.put(an_type)

            for N in range(num_processors):
                type_input_queue.put('STOP')

            all_processes = []

            # close all connections to the db so that they are automatically recreated for each process
            # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
            db.connections.close_all()

            sys.stdout.write('\nCalculating ITS2 type profile abundances clade {}\n'.format(clade_list[i]))
            for N in range(num_processors):
                p = Process(target=output_worker_one,
                            args=(type_input_queue, sample_uids_list, type_output_managed_dict,
                                  sample_uid_to_clade_collection_uids_of_clade))
                all_processes.append(p)
                p.start()

            for p in all_processes:
                p.join()

            # So we have the current list of samples stored in a mangaed list that is the listOfDataSetSamplesManaged
            # the dictionary that we have below is key = analysis_type and the value is two lists
            # the first list is the raw counts and the second is the prop counts
            # In each of these lists the string is a string of the row values for each of the types
            # we can extract the seq abundances from these values.
            # we are currently doing this by clade so we will have to extract this info for each clade
            # I think we can simply add the items for each of the within clade dicts to each other and end up with
            # an across clades dict then we should be able to work with this to get the order we are looking for

            across_clade_type_sample_abund_dict.update(type_output_managed_dict)
            across_clade_sorted_type_order.extend(sorted_list_of_types)
            # ####################

            # for antype in sortedListOfTypes:
            #     outputTableOne.append(typeOutputManagedDict[antype][0])
            #     outputTableTwo.append(typeOutputManagedDict[antype][1])

    # At this stage we have the across_clade_type_sample_abund_dict that contains the type to row values
    # Below is the pseudo code for getting the sorted list of samples
    # get a list of the types for the output
    # Ideally we want to have this list of types sorted according to the most abundant in the output first
    # we currently don't have this except separated by clade
    # We can get this from the dict though
    type_to_abund_list = []
    for analysis_type_obj in across_clade_type_sample_abund_dict.keys():
        type_to_abund_list.append(
            (analysis_type_obj, int(across_clade_type_sample_abund_dict[analysis_type_obj][0].split('\t')[4])))

    # now we can sort this list according to the local abundance and this will give us the order of types that we want
    sorted_analysis_type_abundance_list = [a[0] for a in sorted(type_to_abund_list, key=lambda x: x[1], reverse=True)]

    # for the purposes of doing the sample sorting we will work with the relative df so that we can
    # compare the abundance of types across samples
    list_for_df_absolute = []
    list_for_df_relative = []

    # get list for the aboslute df
    # need to make sure that this is ordered according to across_clade_sorted_type_order
    for an_type in across_clade_sorted_type_order:
        list_for_df_absolute.append(across_clade_type_sample_abund_dict[an_type][0].split('\t'))

    # get list for the relative df
    # need to make sure that this is ordered according to across_clade_sorted_type_order
    for an_type in across_clade_sorted_type_order:
        list_for_df_relative.append(across_clade_type_sample_abund_dict[an_type][1].split('\t'))

    # headers can be same for each
    pre_headers = ['ITS2 type profile UID', 'Clade', 'Majority ITS2 sequence',
                   'Associated species', 'ITS2 type abundance local', 'ITS2 type abundance DB', 'ITS2 type profile']
    sample_headers = [dataSamp.id for dataSamp in list_of_data_set_samples]
    post_headers = ['Sequence accession / SymPortal UID', 'Average defining sequence proportions and [stdev]']
    columns_for_df = pre_headers + sample_headers + post_headers

    # make absolute
    df_absolute = pd.DataFrame(list_for_df_absolute, columns=columns_for_df)
    df_absolute.set_index('ITS2 type profile', drop=False, inplace=True)

    df_relative = pd.DataFrame(list_for_df_relative, columns=columns_for_df)
    df_relative.set_index('ITS2 type profile', drop=False, inplace=True)

    # add a series that gives us the uids of the samples incase we have samples that have the same names
    # this rather complex comprehension puts an nan into the list for the pre_header headers
    # (i.e not sample headers) and then puts an ID value for the samples
    # this seires works because we rely on the fact that it will just automatically
    # put nan values for all of the headers after the samples
    data_list_for_sample_name_series = [np.nan if i < len(pre_headers)
                                        else list_of_data_set_samples[i - len(pre_headers)].name
                                        for i in range(len(pre_headers) + len(sample_headers))]
    sample_name_series = pd.Series(
                        name='sample_name',
                        data=data_list_for_sample_name_series,
                        index=list(df_absolute)[:len(data_list_for_sample_name_series)])

    # it turns out that you cannot have duplicate header values (which makes sense). So we're going to have to
    # work with the sample uids as the header values and put the sample_name in as the secondary series

    # now add the series to the df and then re order the df
    df_absolute = df_absolute.append(sample_name_series)
    df_relative = df_relative.append(sample_name_series)

    # now reorder the index so that the sample_id_series is on top
    index_list = df_absolute.index.values.tolist()
    re_index_index = [index_list[-1]] + index_list[:-1]
    df_absolute = df_absolute.reindex(re_index_index)
    df_relative = df_relative.reindex(re_index_index)

    # at this point we have both of the dfs. We will use the relative df for getting the ordered smpl list
    # now go sample by sample find the samples max type and add to the dictionary where key is types, and value
    # is list of tups, one for each sample which is sample name and rel_abund of the given type
    # We should also produce a dict that holds the ID to sample_name for reordering purposes later on.
    # sample_id_to_sample_name_dict = {}
    type_to_sample_abund_dict = defaultdict(list)
    typeless_samples_list_by_uid = []
    for i in range(len(pre_headers), len(pre_headers) + len(list_of_data_set_samples)):
        sys.stdout.write('\rGetting type abundance information for {}'.format(
            list_of_data_set_samples[i - len(pre_headers)]))
        sample_series = df_relative.iloc[:, i]
        sample_abundances_series = sample_series[1:].astype('float')
        max_type_label = sample_abundances_series.idxmax()
        rel_abund_of_max_type = sample_abundances_series[max_type_label]
        if not rel_abund_of_max_type > 0:
            # append the ID of the sample to the list
            smpl_id = sample_series.name
            typeless_samples_list_by_uid.append(smpl_id)
            # sample_id_to_sample_name_dict[smpl_id] = sample_series['sample_name']
        else:
            # append a tuple that is (sample_id, rel abundance of the max type)
            smpl_id = sample_series.name
            type_to_sample_abund_dict[max_type_label].append((smpl_id, rel_abund_of_max_type))
            # sample_id_to_sample_name_dict[smpl_id] = sample_series.name

    # here we have the dictionary populated. We can now go type by type according
    # to the sorted_analysis_type_abundance_list and put the samples that had the given type as their most abundant
    # type, into the sorted sample list, addtionaly sorted by how abund the type was in each of the samples
    samples_by_uid_that_have_been_sorted = []
    # we are only concerned with the types that had samples that had them as most abundant
    for an_type_name in [at.name for at in sorted_analysis_type_abundance_list
                         if at.name in type_to_sample_abund_dict.keys()]:
        samples_by_uid_that_have_been_sorted.extend(
            [a[0] for a in sorted(type_to_sample_abund_dict[an_type_name], key=lambda x: x[1], reverse=True)])

    # here we should have a list of samples that have been sorted according to the types they were found to
    # have as their most abundant
    # now we just need to add the samples that didn't have a type in them to be associated to. Negative etc.
    samples_by_uid_that_have_been_sorted.extend(typeless_samples_list_by_uid)

    # now pickle out the samples_that_have_been_sorted list if we are running on the remote system
    output_directory = os.path.abspath(
        os.path.join(os.path.dirname(__file__), 'outputs/analyses/{}'.format(analysis_object.id)))
    os.makedirs(output_directory, exist_ok=True)
    if sp_config.system_type == 'remote':
        pickle.dump(samples_by_uid_that_have_been_sorted,
                    open("{}/samples_that_have_been_sorted.pickle".format(output_directory), "wb"))

    # rearange the sample columns so that they are in the new order
    new_sample_headers = samples_by_uid_that_have_been_sorted
    new_cols = pre_headers + new_sample_headers + post_headers

    df_absolute = df_absolute[new_cols]
    df_relative = df_relative[new_cols]

    # transpose
    df_absolute = df_absolute.T
    df_relative = df_relative.T

    os.chdir(output_directory)

    # Finally append the species references to the tables
    species_ref_dict = {
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

    species_set = set()
    for analysis_type_obj in across_clade_sorted_type_order:
        if analysis_type_obj.species != 'None':
            species_set.update(analysis_type_obj.species.split(','))

    # put in the species reference title
    temp_series = pd.Series()
    temp_series.name = 'Species references'
    df_absolute = df_absolute.append(temp_series)
    df_relative = df_relative.append(temp_series)

    # now add the references for each of the associated species
    for species in species_set:
        if species in species_ref_dict.keys():
            temp_series = pd.Series([species_ref_dict[species]], index=[list(df_relative)[0]])
            temp_series.name = species
            df_absolute = df_absolute.append(temp_series)
            df_relative = df_relative.append(temp_series)

    # Now append the meta infromation for the output. i.e. the user running the analysis or the standalone
    # and the data_set submissions used in the analysis.
    # two scenarios in which this can be called
    # from the analysis: call_type = 'analysis'
    # or as a stand alone: call_type = 'stand_alone'

    if call_type == 'analysis':
        meta_info_string_items = [
            'Output as part of data_analysis ID: {}; '
            'Number of data_set objects as part of analysis = {}; '
            'submitting_user: {}; time_stamp: {}'.format(
                analysisobj.id, len(query_set_of_data_sets), analysisobj.submitting_user, analysisobj.time_stamp)]
        temp_series = pd.Series(meta_info_string_items, index=[list(df_relative)[0]], name='meta_info_summary')
        df_absolute = df_absolute.append(temp_series)
        df_relative = df_relative.append(temp_series)

        for data_set_object in query_set_of_data_sets:
            data_set_meta_list = [
                'Data_set ID: {}; Data_set name: {}; submitting_user: {}; time_stamp: {}'.format(
                    data_set_object.id, data_set_object.name,
                    data_set_object.submitting_user, data_set_object.time_stamp)]

            temp_series = pd.Series(data_set_meta_list, index=[list(df_relative)[0]], name='data_set_info')
            df_absolute = df_absolute.append(temp_series)
            df_relative = df_relative.append(temp_series)
    else:
        # call_type=='stand_alone'
        meta_info_string_items = [
            'Stand_alone output by {} on {}; '
            'data_analysis ID: {}; '
            'Number of data_set objects as part of output = {}; '
            'Number of data_set objects as part of analysis = {}'.format(
                output_user, str(datetime.now()).replace(' ', '_').replace(':', '-'), analysisobj.id,
                len(query_set_of_data_sets), len(analysisobj.list_of_data_set_uids.split(',')))]

        temp_series = pd.Series(meta_info_string_items, index=[list(df_relative)[0]], name='meta_info_summary')
        df_absolute = df_absolute.append(temp_series)
        df_relative = df_relative.append(temp_series)
        for data_set_object in query_set_of_data_sets:
            data_set_meta_list = [
                'Data_set ID: {}; Data_set name: {}; submitting_user: {}; time_stamp: {}'.format(
                    data_set_object.id, data_set_object.name,
                    data_set_object.submitting_user, data_set_object.time_stamp)]

            temp_series = pd.Series(data_set_meta_list, index=[list(df_relative)[0]], name='data_set_info')
            df_absolute = df_absolute.append(temp_series)
            df_relative = df_relative.append(temp_series)

    if time_date_str:
        date_time_string = time_date_str
    else:
        date_time_string = str(datetime.now()).replace(' ', '_').replace(':', '-')
    os.makedirs(output_directory, exist_ok=True)

    path_to_profiles_absolute = '{}/{}_{}_{}.profiles.absolute.txt'.format(
        output_directory, analysis_object.id, analysis_object.name, date_time_string)

    df_absolute.to_csv(path_to_profiles_absolute, sep="\t", header=False)
    output_files_list.append(path_to_profiles_absolute)

    del df_absolute

    path_to_profiles_rel = '{}/{}_{}_{}.profiles.relative.txt'.format(
        output_directory, analysis_object.id, analysis_object.name, date_time_string)

    df_relative.to_csv(path_to_profiles_rel, sep="\t", header=False)
    # write_list_to_destination(path_to_profiles_rel, outputTableTwo)
    output_files_list.append(path_to_profiles_rel)

    del df_relative

    # ########################## ITS2 INTRA ABUND COUNT TABLE ################################
    div_output_file_list, date_time_string, number_of_samples = output_sequence_count_tables(
        datasubstooutput=datasubstooutput, num_processors=num_processors, output_dir=output_directory,
        sorted_sample_uid_list=samples_by_uid_that_have_been_sorted, analysis_obj_id=analysisobj.id,
        call_type='analysis', time_date_str=date_time_string)

    print('ITS2 type profile output files:')
    for output_file in output_files_list:
        print(output_file)
        if 'relative' in output_file:
            output_to_plot = output_file
            break

    output_files_list.extend(div_output_file_list)
    # Finally lets produce output plots for the dataoutput. For the time being this should just be a
    # plot for the ITS2 type profiles and one for the sequences
    # as with the data_submission let's pass in the path to the outputfiles that we can use to make the plot with
    output_dir = os.path.dirname(output_to_plot)
    if not no_figures:
        if number_of_samples > 1000:
            print('Too many samples ({}) to generate plots'.format(number_of_samples))
        else:
            svg_path, png_path, sorted_sample_id_list = generate_stacked_bar_data_analysis_type_profiles(
                path_to_tab_delim_count=output_to_plot, output_directory=output_dir,
                analysis_obj_id=analysisobj.id, time_date_str=date_time_string)

            print('Figure output files:')
            print(svg_path)
            print(png_path)
            output_files_list.extend([svg_path, png_path])
            for file in div_output_file_list:
                if 'relative' in file:
                    path_to_plot = file
                    break

            svg_path, png_path = generate_stacked_bar_data_loading(
                path_to_tab_delim_count=path_to_plot, output_directory=output_dir,
                time_date_str=date_time_string, sample_id_order_list=sorted_sample_id_list)

            print('Figure output files:')
            print(svg_path)
            print(png_path)
            output_files_list.extend([svg_path, png_path])

    return output_dir, date_time_string, output_files_list

def output_type_count_tables_data_set_sample_id_input(
        analysisobj, data_set_sample_ids_to_output_string,
        num_processors=1, no_figures=False, output_user=None, time_date_str=None):
    analysis_object = analysisobj
    # This is one of the last things to do before we can use our first dataset
    # The table will have types as columns and rows as samples
    # its rows will give various data about each type as well as how much they were present in each sample

    # It will produce four output files, for DIV abundances and proportions and Type abundances and proportions.
    # found in the given clade collection.
    # Types will be listed first by clade and then by the number of clade collections the types were found in

    # Each type's ID will also be given as a UID.
    # The date the database was accessed and the version should also be noted
    # The formal species descriptions which correspond to the found ITS2 type will be noted for each type
    # Finally the AccessionNumber of each of the defining reference species will also be noted
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    # List of the paths to the files that have been output
    output_files_list = []

    data_set_sample_ids_to_output = [int(a) for a in str(data_set_sample_ids_to_output_string).split(',')]



    # Get collection of types that are specific for the dataSubmissions we are looking at

    query_set_of_data_set_samples = DataSetSample.objects.filter(id__in=data_set_sample_ids_to_output)

    clade_collections_from_data_set_samples = CladeCollection.objects.filter(
        data_set_sample_from__in=query_set_of_data_set_samples)


    clade_collection_types_from_this_data_analysis_and_data_set = CladeCollectionType.objects.filter(
        clade_collection_found_in__in=clade_collections_from_data_set_samples,
        analysis_type_of__data_analysis_from=analysis_object).distinct()

    at = set()
    for cct in clade_collection_types_from_this_data_analysis_and_data_set:
        sys.stdout.write('\rCollecting analysis_type from clade collection type {}'.format(cct))
        at.add(cct.analysis_type_of)
    at = list(at)

    # Need to get a list of Samples that are part of the dataAnalysis
    list_of_data_set_samples = list(query_set_of_data_set_samples)

    # Now go through the types, firstly by clade and then by the number of cladeCollections they were found in
    # Populate a 2D list of types with a list per clade
    types_cladal_list = [[] for _ in clade_list]
    for att in at:
        try:
            if len(att.list_of_clade_collections) > 0:
                types_cladal_list[clade_list.index(att.clade)].append(att)
        except:
            pass

    # Items for for creating the new samples sorted output
    across_clade_type_sample_abund_dict = dict()
    across_clade_sorted_type_order = []
    # Go clade by clade
    for i in range(len(types_cladal_list)):
        if types_cladal_list[i]:
            clade_in_question = clade_list[i]
            # ##### CALCULATE REL ABUND AND SD OF DEF INTRAS FOR THIS TYPE ###############
            #     # For each type we want to calculate the average proportions of the defining seqs in that type
            #     # We will also calculate an SD.
            #     # We will do both of these calculations using the footprint_sequence_abundances, list_of_clade_collections
            #     # and orderedFoorprintList attributes of each type

            # ######### MAKE GROUP COUNTER AND DICT ###########
            # These are now going to be managed items for use in the MP
            # We want to name the groups that types are found in sequencially
            # To get the next number to name a group we will use the groupCount
            # To look up what number has been assigned to groups that have
            # already been printed we will use the groupDict
            # groupCount = 0
            # groupDict = {}
            ###################################################

            # sort types by the number of samples they were found in for this output (not across whole analysis)
            # returns list of tuples, [0] = analysis_type object, [1] number of ccs found in for this output
            sorted_list_of_types = sort_list_of_types_by_clade_collections_in_current_output(
                types_cladal_list[i], clade_collection_types_from_this_data_analysis_and_data_set)

            # Here we will MP at the type level
            # i.e. each type will be processed on a differnt core. In order for this to work we should have a managed
            # dictionary where the key can be the types ID and the value can be the datalist
            # In order to get the types output in the correct order we should use the sortedListOfTypes to resort the
            # data once the MPing has been done.

            # This dict will hold all of the output rows that the MPing has created.
            worker_manager = Manager()
            type_output_managed_dict = worker_manager.dict({an_type: None for an_type in sorted_list_of_types})

            # NB using shared items was considerably slowing us down so now I will just use copies of items
            # we don't actually need to use managed items as they don't really need to be shared (i.e. input
            # from different processes.
            # a list that is the uids of each of the samples in the analyis
            sample_uids_list = [smp.id for smp in list_of_data_set_samples]
            # listOfDataSetSampleIDsManagedList = worker_manager.list(sample_uids_list)

            # A corresponding dictionary that is the list of the clade collection uids for each of the samples
            # that are in the ID list above.
            sample_uid_to_clade_collection_uids_of_clade = {}
            print('\nCreating sample_ID_to_cc_IDs dictionary clade {}'.format(clade_in_question))
            for smp in list_of_data_set_samples:
                sys.stdout.write('\rCollecting clade collections for sample {}'.format(smp))
                try:
                    sample_uid_to_clade_collection_uids_of_clade[smp.id] = CladeCollection.objects.get(
                        data_set_sample_from=smp, clade=clade_in_question).id
                except CladeCollection.DoesNotExist:
                    sample_uid_to_clade_collection_uids_of_clade[smp.id] = None
                except Exception as ex:  # Just incase there is some weird stuff going on
                    print(ex)

            type_input_queue = Queue()

            for an_type in sorted_list_of_types:
                type_input_queue.put(an_type)

            for N in range(num_processors):
                type_input_queue.put('STOP')

            all_processes = []

            # close all connections to the db so that they are automatically recreated for each process
            # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
            db.connections.close_all()

            sys.stdout.write('\nCalculating ITS2 type profile abundances clade {}\n'.format(clade_list[i]))
            for N in range(num_processors):
                p = Process(target=output_worker_one,
                            args=(type_input_queue, sample_uids_list, type_output_managed_dict,
                                  sample_uid_to_clade_collection_uids_of_clade))
                all_processes.append(p)
                p.start()

            for p in all_processes:
                p.join()

            # So we have the current list of samples stored in a mangaed list that is the listOfDataSetSamplesManaged
            # the dictionary that we have below is key = analysis_type and the value is two lists
            # the first list is the raw counts and the second is the prop counts
            # In each of these lists the string is a string of the row values for each of the types
            # we can extract the seq abundances from these values.
            # we are currently doing this by clade so we will have to extract this info for each clade
            # I think we can simply add the items for each of the within clade dicts to each other and end up with
            # an across clades dict then we should be able to work with this to get the order we are looking for

            across_clade_type_sample_abund_dict.update(type_output_managed_dict)
            across_clade_sorted_type_order.extend(sorted_list_of_types)
            # ####################

            # for antype in sortedListOfTypes:
            #     outputTableOne.append(typeOutputManagedDict[antype][0])
            #     outputTableTwo.append(typeOutputManagedDict[antype][1])

    # At this stage we have the across_clade_type_sample_abund_dict that contains the type to row values
    # Below is the pseudo code for getting the sorted list of samples
    # get a list of the types for the output
    # Ideally we want to have this list of types sorted according to the most abundant in the output first
    # we currently don't have this except separated by clade
    # We can get this from the dict though
    type_to_abund_list = []
    for analysis_type_obj in across_clade_type_sample_abund_dict.keys():
        type_to_abund_list.append(
            (analysis_type_obj, int(across_clade_type_sample_abund_dict[analysis_type_obj][0].split('\t')[4])))

    # now we can sort this list according to the local abundance and this will give us the order of types that we want
    sorted_analysis_type_abundance_list = [a[0] for a in sorted(type_to_abund_list, key=lambda x: x[1], reverse=True)]

    # for the purposes of doing the sample sorting we will work with the relative df so that we can
    # compare the abundance of types across samples
    list_for_df_absolute = []
    list_for_df_relative = []

    # get list for the aboslute df
    # need to make sure that this is ordered according to across_clade_sorted_type_order
    for an_type in across_clade_sorted_type_order:
        list_for_df_absolute.append(across_clade_type_sample_abund_dict[an_type][0].split('\t'))

    # get list for the relative df
    # need to make sure that this is ordered according to across_clade_sorted_type_order
    for an_type in across_clade_sorted_type_order:
        list_for_df_relative.append(across_clade_type_sample_abund_dict[an_type][1].split('\t'))

    # headers can be same for each
    pre_headers = ['ITS2 type profile UID', 'Clade', 'Majority ITS2 sequence',
                   'Associated species', 'ITS2 type abundance local', 'ITS2 type abundance DB', 'ITS2 type profile']
    sample_headers = [dataSamp.id for dataSamp in list_of_data_set_samples]
    post_headers = ['Sequence accession / SymPortal UID', 'Average defining sequence proportions and [stdev]']
    columns_for_df = pre_headers + sample_headers + post_headers

    # make absolute
    df_absolute = pd.DataFrame(list_for_df_absolute, columns=columns_for_df)
    df_absolute.set_index('ITS2 type profile', drop=False, inplace=True)

    df_relative = pd.DataFrame(list_for_df_relative, columns=columns_for_df)
    df_relative.set_index('ITS2 type profile', drop=False, inplace=True)

    # add a series that gives us the uids of the samples incase we have samples that have the same names
    # this rather complex comprehension puts an nan into the list for the pre_header headers
    # (i.e not sample headers) and then puts an ID value for the samples
    # this seires works because we rely on the fact that it will just automatically
    # put nan values for all of the headers after the samples
    data_list_for_sample_name_series = [np.nan if i < len(pre_headers)
                                        else list_of_data_set_samples[i - len(pre_headers)].name
                                        for i in range(len(pre_headers) + len(sample_headers))]
    sample_name_series = pd.Series(
                        name='sample_name',
                        data=data_list_for_sample_name_series,
                        index=list(df_absolute)[:len(data_list_for_sample_name_series)])

    # it turns out that you cannot have duplicate header values (which makes sense). So we're going to have to
    # work with the sample uids as the header values and put the sample_name in as the secondary series

    # now add the series to the df and then re order the df
    df_absolute = df_absolute.append(sample_name_series)
    df_relative = df_relative.append(sample_name_series)

    # now reorder the index so that the sample_id_series is on top
    index_list = df_absolute.index.values.tolist()
    re_index_index = [index_list[-1]] + index_list[:-1]
    df_absolute = df_absolute.reindex(re_index_index)
    df_relative = df_relative.reindex(re_index_index)

    # at this point we have both of the dfs. We will use the relative df for getting the ordered smpl list
    # now go sample by sample find the samples max type and add to the dictionary where key is types, and value
    # is list of tups, one for each sample which is sample name and rel_abund of the given type
    # We should also produce a dict that holds the ID to sample_name for reordering purposes later on.
    # sample_id_to_sample_name_dict = {}
    type_to_sample_abund_dict = defaultdict(list)
    typeless_samples_list_by_uid = []
    for i in range(len(pre_headers), len(pre_headers) + len(list_of_data_set_samples)):
        sys.stdout.write('\rGetting type abundance information for {}'.format(
            list_of_data_set_samples[i - len(pre_headers)]))
        sample_series = df_relative.iloc[:, i]
        sample_abundances_series = sample_series[1:].astype('float')
        max_type_label = sample_abundances_series.idxmax()
        rel_abund_of_max_type = sample_abundances_series[max_type_label]
        if not rel_abund_of_max_type > 0:
            # append the ID of the sample to the list
            smpl_id = sample_series.name
            typeless_samples_list_by_uid.append(smpl_id)
            # sample_id_to_sample_name_dict[smpl_id] = sample_series['sample_name']
        else:
            # append a tuple that is (sample_id, rel abundance of the max type)
            smpl_id = sample_series.name
            type_to_sample_abund_dict[max_type_label].append((smpl_id, rel_abund_of_max_type))
            # sample_id_to_sample_name_dict[smpl_id] = sample_series.name

    # here we have the dictionary populated. We can now go type by type according
    # to the sorted_analysis_type_abundance_list and put the samples that had the given type as their most abundant
    # type, into the sorted sample list, addtionaly sorted by how abund the type was in each of the samples
    samples_by_uid_that_have_been_sorted = []
    # we are only concerned with the types that had samples that had them as most abundant
    for an_type_name in [at.name for at in sorted_analysis_type_abundance_list
                         if at.name in type_to_sample_abund_dict.keys()]:
        samples_by_uid_that_have_been_sorted.extend(
            [a[0] for a in sorted(type_to_sample_abund_dict[an_type_name], key=lambda x: x[1], reverse=True)])

    # here we should have a list of samples that have been sorted according to the types they were found to
    # have as their most abundant
    # now we just need to add the samples that didn't have a type in them to be associated to. Negative etc.
    samples_by_uid_that_have_been_sorted.extend(typeless_samples_list_by_uid)

    # now pickle out the samples_that_have_been_sorted list if we are running on the remote system
    output_directory = os.path.abspath(
        os.path.join(os.path.dirname(__file__), 'outputs/analyses/{}'.format(analysis_object.id)))
    os.makedirs(output_directory, exist_ok=True)
    local_or_remote = sp_config.system_type
    if local_or_remote == 'remote':
        pickle.dump(samples_by_uid_that_have_been_sorted,
                    open("{}/samples_that_have_been_sorted.pickle".format(output_directory), "wb"))

    # rearange the sample columns so that they are in the new order
    new_sample_headers = samples_by_uid_that_have_been_sorted
    new_cols = pre_headers + new_sample_headers + post_headers

    df_absolute = df_absolute[new_cols]
    df_relative = df_relative[new_cols]

    # transpose
    df_absolute = df_absolute.T
    df_relative = df_relative.T

    os.chdir(output_directory)

    # Finally append the species references to the tables
    species_ref_dict = {
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

    species_set = set()
    for analysis_type_obj in across_clade_sorted_type_order:
        if analysis_type_obj.species != 'None':
            species_set.update(analysis_type_obj.species.split(','))

    # put in the species reference title
    temp_series = pd.Series()
    temp_series.name = 'Species references'
    df_absolute = df_absolute.append(temp_series)
    df_relative = df_relative.append(temp_series)

    # now add the references for each of the associated species
    for species in species_set:
        if species in species_ref_dict.keys():
            temp_series = pd.Series([species_ref_dict[species]], index=[list(df_relative)[0]])
            temp_series.name = species
            df_absolute = df_absolute.append(temp_series)
            df_relative = df_relative.append(temp_series)

    # Now append the meta infromation for the output. i.e. the user running the analysis or the standalone
    # and the data_set submissions used in the analysis.
    # two scenarios in which this can be called
    # from the analysis: call_type = 'analysis'
    # or as a stand alone: call_type = 'stand_alone'


    # call_type=='stand_alone' call type will always be standalone when doing a data_set_sample id output
    data_sets_of_the_data_set_samples_string = DataSet.objects.filter(
        datasetsample__in=query_set_of_data_set_samples).distinct()
    meta_info_string_items = [
        'Stand_alone output by {} on {}; '
        'data_analysis ID: {}; '
        'Number of data_set objects as part of output = {}; '
        'Number of data_set objects as part of analysis = {}'.format(
            output_user, str(datetime.now()).replace(' ', '_').replace(':', '-'), analysisobj.id,
            len(data_sets_of_the_data_set_samples_string), len(analysisobj.list_of_data_set_uids.split(',')))]

    temp_series = pd.Series(meta_info_string_items, index=[list(df_relative)[0]], name='meta_info_summary')
    df_absolute = df_absolute.append(temp_series)
    df_relative = df_relative.append(temp_series)
    for data_set_object in data_sets_of_the_data_set_samples_string:
        data_set_meta_list = [
            'Data_set ID: {}; Data_set name: {}; submitting_user: {}; time_stamp: {}'.format(
                data_set_object.id, data_set_object.name,
                data_set_object.submitting_user, data_set_object.time_stamp)]

        temp_series = pd.Series(data_set_meta_list, index=[list(df_relative)[0]], name='data_set_info')
        df_absolute = df_absolute.append(temp_series)
        df_relative = df_relative.append(temp_series)

    if time_date_str:
        date_time_string = time_date_str
    else:
        date_time_string = str(datetime.now()).replace(' ', '_').replace(':', '-')
    os.makedirs(output_directory, exist_ok=True)

    path_to_profiles_absolute = '{}/{}_{}_{}.profiles.absolute.txt'.format(
        output_directory, analysis_object.id, analysis_object.name, date_time_string)

    df_absolute.to_csv(path_to_profiles_absolute, sep="\t", header=False)
    output_files_list.append(path_to_profiles_absolute)

    del df_absolute

    path_to_profiles_rel = '{}/{}_{}_{}.profiles.relative.txt'.format(
        output_directory, analysis_object.id, analysis_object.name, date_time_string)

    df_relative.to_csv(path_to_profiles_rel, sep="\t", header=False)
    # write_list_to_destination(path_to_profiles_rel, outputTableTwo)
    output_files_list.append(path_to_profiles_rel)

    del df_relative

    # ########################## ITS2 INTRA ABUND COUNT TABLE ################################
    div_output_file_list, date_time_string, number_of_samples = div_output_pre_analysis_new_meta_and_new_dss_structure_data_set_sample_id_input(
        data_set_sample_ids_to_output_string=data_set_sample_ids_to_output_string, num_processors=num_processors,
        output_dir=output_directory, sorted_sample_uid_list=samples_by_uid_that_have_been_sorted,
        analysis_obj_id=analysisobj.id, time_date_str=date_time_string)

    output_to_plot = None
    print('ITS2 type profile output files:')
    for output_file in output_files_list:
        print(output_file)
        if 'relative' in output_file:
            output_to_plot = output_file
            break

    output_files_list.extend(div_output_file_list)
    # Finally lets produce output plots for the dataoutput. For the time being this should just be a
    # plot for the ITS2 type profiles and one for the sequences
    # as with the data_submission let's pass in the path to the outputfiles that we can use to make the plot with
    output_dir = os.path.dirname(output_to_plot)
    if not no_figures:
        if number_of_samples > 1000:
            print('Too many samples ({}) to generate plots'.format(number_of_samples))
        else:
            svg_path, png_path, sorted_sample_id_list = generate_stacked_bar_data_analysis_type_profiles(
                path_to_tab_delim_count=output_to_plot, output_directory=output_dir,
                analysis_obj_id=analysisobj.id, time_date_str=date_time_string)

            print('Figure output files:')
            print(svg_path)
            print(png_path)
            output_files_list.extend([svg_path, png_path])
            for file in div_output_file_list:
                if 'relative' in file:
                    path_to_plot = file
                    break

            svg_path, png_path = plotting.generate_stacked_bar_data_loading(
                path_to_tab_delim_count=path_to_plot, output_directory=output_dir,
                time_date_str=date_time_string, sample_id_order_list=sorted_sample_id_list)

            print('Figure output files:')
            print(svg_path)
            print(png_path)
            output_files_list.extend([svg_path, png_path])

    return output_dir, date_time_string, output_files_list

def sort_list_of_types_by_clade_collections_in_current_output(
        list_of_analysis_types, clade_collection_types_in_current_output):
    # create a list of tupples that is the type ID and the number of this output's CCs that it was found in.
    sys.stdout.write('\nSorting types by abundance of clade collections\n')
    tuple_list = []
    clade_collections_uid_in_current_output = [
        clade_collection_type_obj.clade_collection_found_in.id
        for clade_collection_type_obj in clade_collection_types_in_current_output]

    for at in list_of_analysis_types:
        sys.stdout.write('\rCounting for type {}'.format(at))
        list_of_clade_collections_found_in = [int(x) for x in at.list_of_clade_collections.split(',')]
        num_clade_collections_of_output = list(
            set(clade_collections_uid_in_current_output).intersection(list_of_clade_collections_found_in))
        tuple_list.append((at.id, len(num_clade_collections_of_output)))

    type_uids_sorted = sorted(tuple_list, key=lambda x: x[1], reverse=True)

    return [AnalysisType.objects.get(id=x[0]) for x in type_uids_sorted]


def output_worker_one(
        input_queue, list_of_data_set_samples_uids, output_dictionary, sample_uid_to_clade_collection_of_clade_uid):
    num_samples = len(list_of_data_set_samples_uids)
    for an_type in iter(input_queue.get, 'STOP'):  # Within each type go through each of the samples

        sys.stdout.write('\rProcessing ITS2 type profile: {}'.format(an_type))

        # ##### CALCULATE REL ABUND AND SD OF DEF INTRAS FOR THIS TYPE ###############
        # For each type we want to calculate the average proportions of the defining seqs in that type
        # We will also calculate an SD.
        # We will do both of these calculations using the footprint_sequence_abundances, list_of_clade_collections
        # and orderedFoorprintList attributes of each type
        footprint_abundances = json.loads(an_type.footprint_sequence_abundances)

        # We have to make a decision as to whether this average should represent all the findings of this type
        # or whether we should only represent the averages of samples found in this dataset.
        # I think it needs to be a global average. Because the type is defined based on all samples in the
        # SymPortal db.

        # Calculate the average proportion of each DIV as a proportion of the absolute abundances of the divs
        # of the type within the samples the type is found in
        div_abundance_df = pd.DataFrame(footprint_abundances)
        # https://stackoverflow.com/questions/26537878/pandas-sum-across-columns-and-divide-each-cell-from-that-value
        # convert each cell to a proportion as a function of the sum of the row
        div_abundance_df_proportion = div_abundance_df.div(div_abundance_df.sum(axis=1),
                                                           axis=0)
        div_abundance_df_proportion_transposed = div_abundance_df_proportion.T

        total_list = list(div_abundance_df_proportion_transposed.mean(axis=1))
        standard_deviation_list = list(div_abundance_df_proportion_transposed.std(axis=1))

        # The total_list now contains the proportions of each def seq
        # The standard_deviation_list now contains the SDs for each of the def seqs proportions
        ###########################################################################

        # Counter that will increment with every sample type is found in
        # This is counter will end up being the type abund local value
        abundance_count = 0

        type_in_question = an_type
        clade = type_in_question.clade

        # For each type create a holder that will hold 0s for each sample until populated below
        data_row_raw = [0 for _ in range(num_samples)]
        data_row_proportion = [0 for _ in range(num_samples)]
        type_clade_collection_uids = [int(a) for a in type_in_question.list_of_clade_collections.split(',')]

        # We also need the type abund db value. We can get this from the type cladeCollections
        global_count = len(type_clade_collection_uids)

        # Within each type go through each of the samples
        # Do we really have to go through every sample? I don't think so.
        # Because we have only one cc ID per sample we can simply identify the
        # sample ID (keys) in the sample_uid_to_clade_collection_of_clade_uid dict where
        # the cc ID (value) is found in the type_clade_collection_uids.
        uids_of_samples_that_had_type = [
            smp_id for smp_id in list_of_data_set_samples_uids if
            sample_uid_to_clade_collection_of_clade_uid[smp_id] in type_clade_collection_uids]
        for ID in uids_of_samples_that_had_type:
            abundance_count += 1
            # Need to work out how many seqs were found from the sample for this type
            # Also need to work this out as a proportion of all of the Symbiodinium seqs found in the sample
            clade_collection_in_question_uid = sample_uid_to_clade_collection_of_clade_uid[ID]

            # The number of sequences that were returned for the sample in question
            total_number_of_sequences = DataSetSample.objects.get(id=ID).absolute_num_sym_seqs

            # The number of sequences that make up the type in q
            # The index of the clade_collection_object in the type's list_of_clade_collections
            clade_collection_index_in_type = type_clade_collection_uids.index(clade_collection_in_question_uid)
            # List containing the abundances of each of the ref seqs that
            # make up the type in the given clade_collection_object
            sequence_abundance_info_for_clade_collection_and_type_in_question = json.loads(
                type_in_question.footprint_sequence_abundances)[clade_collection_index_in_type]
            # total of the above list, i.e. seqs in type
            sum_of_defining_reference_sequence_abundances_in_type = sum(
                sequence_abundance_info_for_clade_collection_and_type_in_question)
            # type abundance as proportion of the total seqs found in the sample
            type_proportion = sum_of_defining_reference_sequence_abundances_in_type / total_number_of_sequences
            # Now populate the dataRow with the sum_of_defining_reference_sequence_abundances_in_type
            # and the type_proportion
            index_of_sample_uid_in_list_of_data_set_sample_uids = list_of_data_set_samples_uids.index(ID)
            data_row_raw[index_of_sample_uid_in_list_of_data_set_sample_uids] = \
                sum_of_defining_reference_sequence_abundances_in_type
            data_row_proportion[index_of_sample_uid_in_list_of_data_set_sample_uids] = type_proportion

        # Type Profile
        type_uid = an_type.id

        type_profile_name = type_in_question.name
        species = type_in_question.species
        type_abundance = abundance_count
        # This is currently putting out the Majs in an odd order due to the fact that the majority_reference_sequence_set is
        # a set. Instead we should get the Majs from the name.
        # majority_its2 = '/'.join([str(refseq) for refseq in sortedListOfTypes[j].getmajority_reference_sequence_set()])
        majority_its2, majority_list = get_maj_list(an_type)
        type_abundance_and_standard_deviation_string = get_abundance_string(
            total_list, standard_deviation_list, majority_list)

        sequence_accession = type_in_question.generate_name(accession=True)
        row_one = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            type_uid, clade, majority_its2, species, str(type_abundance), str(global_count), type_profile_name,
            '\t'.join([str(a) for a in data_row_raw]), sequence_accession,
            type_abundance_and_standard_deviation_string)

        row_two = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            type_uid, clade, majority_its2, species, str(type_abundance), str(global_count), type_profile_name,
            '\t'.join(["{:.3f}".format(prop) for prop in data_row_proportion]), sequence_accession,
            type_abundance_and_standard_deviation_string)

        # Here we will use the output_dictionary instead of adding it to the outputTableOne
        # We will then get these elements in the dict into the right order with the types we will then
        # write out each of the type elements in the order of the orderedTypeLists

        output_dictionary[an_type] = [row_one, row_two]


def get_maj_list(atype):

    # This is a little tricky. Have to take into account that the Maj seqs are not always more abundant
    # than the non-Maj seqs.
    # e.g. type B5/B5e-B5a-B5b, B5e is acutally the least abundant of the seqs
    # Therefore we cannot simply grab the orderedFooprintList seqs in order assuming they are the Majs

    name = atype.name
    count = name.count('/')
    majority_list = []
    # list of the seqs in order of abundance across the type's samples
    uids_of_sequences_in_order_of_abundance = atype.ordered_footprint_list.split(',')
    # list of the maj seqs in the type
    majority_sequeces_uids = atype.majority_reference_sequence_set.split(',')
    for index in range(count + 1):
        for item in range(len(uids_of_sequences_in_order_of_abundance)):
            if uids_of_sequences_in_order_of_abundance[item] in majority_sequeces_uids:
                maj_seq_obj = ReferenceSequence.objects.get(id=int(uids_of_sequences_in_order_of_abundance[item]))
                if maj_seq_obj.has_name:
                    majority_list.append(maj_seq_obj.name)
                else:
                    majority_list.append(str(maj_seq_obj.id))
                del uids_of_sequences_in_order_of_abundance[item]
                break
    majority_string_output = '/'.join(majority_list)
    return majority_string_output, majority_list


def get_abundance_string(totlist, sdlist, majlist):
    maj_comp = []
    less_comp = []
    for i in range(len(totlist)):
        total_string = "{0:.3f}".format(totlist[i])
        standard_deviation_string = "{0:.3f}".format(sdlist[i])
        if i in range(len(majlist)):
            maj_comp.append('{}[{}]'.format(total_string, standard_deviation_string))
        else:
            less_comp.append('{}[{}]'.format(total_string, standard_deviation_string))
    maj_comp_str = '/'.join(maj_comp)
    less_comp_str = '-'.join(less_comp)
    if less_comp_str:
        abund_output_str = '-'.join([maj_comp_str, less_comp_str])
    else:
        abund_output_str = maj_comp_str
    return abund_output_str


def output_sequence_count_tables(
        datasubstooutput, num_processors, output_dir, call_type,
        sorted_sample_uid_list=None, analysis_obj_id=None, output_user=None, time_date_str=None):

    # ######################### ITS2 INTRA ABUND COUNT TABLE ################################
    # This is where we're going to have to work with the sequences that aren't part of a type.
    # Esentially we want to end up with the noName sequences divieded up cladally.
    # So at the moment this will mean that we will divide up the current no names but also
    # add information about the other cladal sequences that didn't make it into a cladeCollection

    # list to hold the paths of the outputted files
    output_path_list = []

    # ############### GET ORDERED LIST OF INTRAS BY CLADE THEN ABUNDANCE #################
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    data_sets_to_output = [int(a) for a in datasubstooutput.split(',')]

    # Get collection of data_sets that are specific for the dataSubmissions we are looking at
    query_set_of_data_sets = DataSet.objects.filter(id__in=data_sets_to_output)

    reference_sequences_in_data_sets = ReferenceSequence.objects.filter(
        datasetsamplesequence__data_set_sample_from__data_submission_from__in=query_set_of_data_sets).distinct()

    # Get list of clades found
    clade_set = set()
    for ref_seq in reference_sequences_in_data_sets:
        clade_set.add(ref_seq.clade)

    # Order the list of clades alphabetically
    # the only purpuse of this worker 2 is to end up with a set of sequences ordered by clade and then
    # by abundance across all samples.
    # I was doing this on a clade by clade basis on a sequence by sequence basis.
    # but now we should get rid of the clade sorting and just do that afterwards.
    # we should instad go by a sample by sample basis. This will allow us also to undertake the work that
    # is happening in worker three.
    # two birds one stone.
    # this is actually also better because now we can count the abundance of the sequences in proportions
    # rather than simply absolute values. Which is far more accurate.
    sub_clade_list = [a for a in clade_list if a in clade_set]

    sys.stdout.write('\n')

    # The manager that we will build all of the shared dictionaries from below.
    worker_manager = Manager()

    sample_list = DataSetSample.objects.filter(data_submission_from__in=query_set_of_data_sets)

    # Dictionary that will hold the list of data_set_sample_sequences for each sample
    sample_to_dsss_list_shared_dict = worker_manager.dict()

    print('Creating sample to data_set_sample_sequence dict:')

    for dss in sample_list:
        sys.stdout.write('\r{}'.format(dss.name))
        sample_to_dsss_list_shared_dict[dss.id] = list(
            DataSetSampleSequence.objects.filter(data_set_sample_from=dss))

    # Queue that will hold the data set samples for the MP
    data_set_sample_queue = Queue()

    # 1 - Seqname to cumulative relative abundance for each sequence across all sampples (for getting the over lying order of ref seqs)
    # 2 - sample_id : list(dict(ref_seq_of_sample_name:absolute_abundance_of_dsss_in_sample), dict(ref_seq_of_sample_name:relative_abundance_of_dsss_in_sample))
    # 3 - sample_id : list(dict(clade:total_abund_of_no_name_seqs_of_clade_in_q_), dict(clade:relative_abund_of_no_name_seqs_of_clade_in_q_)

    reference_sequence_names_clade_annotated = [
        ref_seq.name if ref_seq.has_name
        else str(ref_seq.id) + f'_{ref_seq.clade}' for ref_seq in reference_sequences_in_data_sets]

    generic_seq_to_abund_dict = {refSeq_name: 0 for refSeq_name in reference_sequence_names_clade_annotated}

    list_of_dicts_for_processors = []
    for n in range(num_processors):
        list_of_dicts_for_processors.append(
            (worker_manager.dict(generic_seq_to_abund_dict), worker_manager.dict(), worker_manager.dict()))

    for dss in sample_list:
        data_set_sample_queue.put(dss)

    for N in range(num_processors):
        data_set_sample_queue.put('STOP')

    all_processes = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()

    for n in range(num_processors):
        p = Process(target=output_worker_two, args=(
            data_set_sample_queue, list_of_dicts_for_processors[n][0], list_of_dicts_for_processors[n][1],
            list_of_dicts_for_processors[n][2], reference_sequence_names_clade_annotated, sample_to_dsss_list_shared_dict))
        all_processes.append(p)
        p.start()

    for p in all_processes:
        p.join()

    print('\nCollecting results of data_set_sample_counting across {} dictionaries'.format(num_processors))
    master_seq_abundance_counter = Counter()
    master_smple_seq_dict = dict()
    master_smple_no_name_clade_summary = dict()

    # now we need to do different actions for each of the three dictionary sets. One set for each num_proc
    # for the seqName counter we simply need to add to the master counter as we were doing before
    # for both of the sample-centric dictionaries we simply need to update a master dictionary
    for n in range(len(list_of_dicts_for_processors)):
        sys.stdout.write('\rdictionary {}(0)/{}'.format(n, num_processors))
        master_seq_abundance_counter += Counter(dict(list_of_dicts_for_processors[n][0]))
        sys.stdout.write('\rdictionary {}(1)/{}'.format(n, num_processors))
        master_smple_seq_dict.update(dict(list_of_dicts_for_processors[n][1]))
        sys.stdout.write('\rdictionary {}(2)/{}'.format(n, num_processors))
        master_smple_no_name_clade_summary.update(dict(list_of_dicts_for_processors[n][2]))

    print('Collection complete.')

    # we now need to separate by clade and sort within the clade
    clade_abundance_ordered_ref_seq_list = []
    for i in range(len(sub_clade_list)):
        temp_within_clade_list_for_sorting = []
        for seq_name, abund_val in master_seq_abundance_counter.items():
            if seq_name.startswith(sub_clade_list[i]) or seq_name[-2:] == f'_{sub_clade_list[i]}':
                # then this is a seq of the clade in Q and we should add to the temp list
                temp_within_clade_list_for_sorting.append((seq_name, abund_val))
        # now sort the temp_within_clade_list_for_sorting and add to the cladeAbundanceOrderedRefSeqList
        sorted_within_clade = [
            a[0] for a in sorted(temp_within_clade_list_for_sorting, key=lambda x: x[1], reverse=True)]

        clade_abundance_ordered_ref_seq_list.extend(sorted_within_clade)

    # now delete the master_seq_abundance_counter as we are done with it
    del master_seq_abundance_counter

    # ##### WORKER THREE DOMAIN

    # we will eventually have the outputs stored in pandas dataframes.
    # in the worker below we will create a set of pandas.Series for each of the samples which will hold the abundances
    # one for the absoulte abundance and one for the relative abundances.

    # we will put together the headers piece by piece
    header_pre = clade_abundance_ordered_ref_seq_list
    no_name_summary_strings = ['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']]
    qc_stats = [
        'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs', 'post_taxa_id_absolute_symbiodinium_seqs',
        'post_taxa_id_unique_symbiodinium_seqs', 'size_screening_violation_absolute', 'size_screening_violation_unique',
        'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs', 'post_med_absolute',
        'post_med_unique']

    # append the noName sequences as individual sequence abundances
    output_header = ['sample_name'] + qc_stats + no_name_summary_strings + header_pre

    # #####################################################################################

    # ########### POPULATE TABLE WITH clade_collection_object DATA #############

    # In order to MP this we will have to pay attention to the order. As we can't keep the order as we work with the
    # MPs we will do as we did above for the profies outputs and have an output dict that will have the sample as key
    # and its corresponding data row as the value. Then when we have finished the MPing we will go through the
    # sampleList order and output the data in this order.
    # So we will need a managed output dict.

    worker_manager = Manager()
    managed_sample_output_dict = worker_manager.dict()

    data_set_sample_queue = Queue()

    for dss in sample_list:
        data_set_sample_queue.put(dss)

    for N in range(num_processors):
        data_set_sample_queue.put('STOP')

    all_processes = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()

    sys.stdout.write('\n\nOutputting seq data\n')
    for N in range(num_processors):
        p = Process(target=output_worker_three, args=(
            data_set_sample_queue, managed_sample_output_dict, clade_abundance_ordered_ref_seq_list, output_header,
            master_smple_seq_dict, master_smple_no_name_clade_summary))
        all_processes.append(p)
        p.start()

    for p in all_processes:
        p.join()

    print('\nseq output complete\n')

    managed_sample_output_dict_dict = dict(managed_sample_output_dict)

    # now we need to populate the output dataframe with the sample series in order of the sorted_sample_list
    # if there is one.
    # If there is a sorted sample list, make sure that it matches the samples that we are outputting

    # We were having an issue with data_sets having the same names. To fix this, we should do our ordering
    # accoring to the uids of the samples
    # TODO this is where we're at with the refactoring and class abstraction
    if sorted_sample_uid_list:
        sys.stdout.write('\nValidating sorted sample list and ordering dataframe accordingly\n')
        if len(sorted_sample_uid_list) != len(sample_list):
            sys.exit('Number of items in sorted_sample_list do not match those to be outputted!')
        if list(set(sorted_sample_uid_list).difference(set([s.id for s in sample_list]))):
            # then there is a sample that doesn't match up from the sorted_sample_uid_list that
            # has been passed in and the unordered sample list that we are working with in the code
            sys.exit('Sample list passed in does not match sample list from db query')

        # if we got to here then the sorted_sample_list looks good
        # I was originally performing the concat directly on the managedSampleOutputDict but this was starting
        # to produce errors. Starting to work on the managedSampleOutputDict_dict seems to not produce these
        # errors.
        # it may be a good idea to break this down to series by series instead of a one liner so that we can
        # print out progress
        # we can use the
        sys.stdout.write('\rPopulating the absolute dataframe with series. This could take a while...')
        output_df_absolute = pd.concat(
            [list_of_series[0] for list_of_series in managed_sample_output_dict_dict.values()], axis=1)
        sys.stdout.write('\rPopulating the relative dataframe with series. This could take a while...')
        output_df_relative = pd.concat(
            [list_of_series[1] for list_of_series in managed_sample_output_dict_dict.values()], axis=1)

        # now transpose
        output_df_absolute = output_df_absolute.T
        output_df_relative = output_df_relative.T

        # now make sure that the order is correct.
        output_df_absolute = output_df_absolute.reindex(sorted_sample_uid_list)
        output_df_relative = output_df_relative.reindex(sorted_sample_uid_list)

    else:

        # this returns a list which is simply the names of the samples
        # This will order the samples according to which sequence is their most abundant.
        # I.e. samples found to have the sequence which is most abundant in the largest number of sequences
        # will be first. Within each maj sequence, the samples will be sorted by the abundance of that sequence
        # in the sample.
        # At the moment we are also ordering by clade just so that you see samples with the A's at the top
        # of the output so that we minimise the number of 0's in the top left of the output
        # honestly I think we could perhaps get rid of this and just use the over all abundance of the sequences
        # discounting clade. THis is what we do for the clade order when plotting.
        sys.stdout.write('\nGenerating ordered sample list and ordering dataframe accordingly\n')
        ordered_sample_list_by_uid = generate_ordered_sample_list(managed_sample_output_dict_dict)

        # if we got to here then the sorted_sample_list looks good
        # I was originally performing the concat directly on the managedSampleOutputDict but this was starting
        # to produce errors. Starting to work on the managedSampleOutputDict_dict seems to not produce these
        # errors.
        sys.stdout.write('\rPopulating the absolute dataframe with series. This could take a while...')
        output_df_absolute = pd.concat(
            [list_of_series[0] for list_of_series in managed_sample_output_dict_dict.values()], axis=1)
        sys.stdout.write('\rPopulating the relative dataframe with series. This could take a while...')
        output_df_relative = pd.concat(
            [list_of_series[1] for list_of_series in managed_sample_output_dict_dict.values()], axis=1)

        # now transpose
        sys.stdout.write('\rTransposing...')
        output_df_absolute = output_df_absolute.T
        output_df_relative = output_df_relative.T

        # now make sure that the order is correct.
        sys.stdout.write('\rReordering index...')
        output_df_absolute = output_df_absolute.reindex(ordered_sample_list_by_uid)
        output_df_relative = output_df_relative.reindex(ordered_sample_list_by_uid)

    # when adding the accession numbers below, we have to go through every sequence and look up its object
    # We also have to do this when we are outputting the fasta.
    # to prevent us having to make every look up twice, we should also make the fasta at the same time
    # Output a .fasta for of all of the sequences found in the analysis
    # we will write out the fasta right at the end.
    fasta_output_list = []

    # Now add the accesion number / UID for each of the DIVs
    sys.stdout.write('\nGenerating accession and fasta\n')

    # go column name by column name and if the col name is in seq_annotated_name
    # then get the accession and add to the accession_list
    # else do nothing and a blank should be automatically added for us.
    # This was painfully slow because we were doing individual calls to the dictionary
    # I think this will be much faster if do two queries of the db to get the named and
    # non named refseqs and then make two dicts for each of these and use these to populate the below
    reference_sequences_in_data_sets_no_name = ReferenceSequence.objects.filter(
        datasetsamplesequence__data_set_sample_from__data_submission_from__in=query_set_of_data_sets,
        has_name=False).distinct()
    reference_sequences_in_data_sets_has_name = ReferenceSequence.objects.filter(
        datasetsamplesequence__data_set_sample_from__data_submission_from__in=query_set_of_data_sets,
        has_name=True).distinct()
    # no name dict should be a dict of id to sequence
    no_name_dict = {rs.id: rs.sequence for rs in reference_sequences_in_data_sets_no_name}
    # has name dict should be a dict of name to sequence
    has_name_dict = {rs.name: (rs.id, rs.sequence) for rs in reference_sequences_in_data_sets_has_name}

    # for the time being we are going to ignore whether a refseq has an assession as we have not put this
    # into use yet.
    accession_list = []
    num_cols = len(list(output_df_relative))
    for i, col_name in enumerate(list(output_df_relative)):
        sys.stdout.write('\rAppending accession info and creating fasta {}: {}/{}'.format(col_name, i, num_cols))
        if col_name in clade_abundance_ordered_ref_seq_list:
            if col_name[-2] == '_':
                col_name_id = int(col_name[:-2])
                accession_list.append(str(col_name_id))
                fasta_output_list.append('>{}'.format(col_name))
                fasta_output_list.append(no_name_dict[col_name_id])
            else:
                col_name_tup = has_name_dict[col_name]
                accession_list.append(str(col_name_tup[0]))
                fasta_output_list.append('>{}'.format(col_name))
                fasta_output_list.append(col_name_tup[1])
        else:
            accession_list.append(np.nan)

    temp_series = pd.Series(accession_list, name='seq_accession', index=list(output_df_relative))
    output_df_absolute = output_df_absolute.append(temp_series)
    output_df_relative = output_df_relative.append(temp_series)

    # Now append the meta infromation for each of the data_sets that make up the output contents
    # this is information like the submitting user, what the uids of the datasets are etc.
    # There are several ways that this can be called.
    # it can be called as part of the submission: call_type = submission
    # part of an analysis output: call_type = analysis
    # or stand alone: call_type = 'stand_alone'
    # we should have an output for each scenario

    if call_type == 'submission':
        data_set_object = query_set_of_data_sets[0]
        # there will only be one data_set object
        meta_info_string_items = [
            'Output as part of data_set submission ID: {}; submitting_user: {}; time_stamp: {}'.format(
                data_set_object.id, data_set_object.submitting_user, data_set_object.time_stamp)]

        temp_series = pd.Series(meta_info_string_items, index=[list(output_df_absolute)[0]], name='meta_info_summary')
        output_df_absolute = output_df_absolute.append(temp_series)
        output_df_relative = output_df_relative.append(temp_series)
    elif call_type == 'analysis':
        data_analysis_obj = DataAnalysis.objects.get(id=analysis_obj_id)
        meta_info_string_items = [
            'Output as part of data_analysis ID: {}; Number of data_set objects as part of analysis = {}; '
            'submitting_user: {}; time_stamp: {}'.format(
                data_analysis_obj.id, len(data_analysis_obj.list_of_data_set_uids.split(',')),
                data_analysis_obj.submitting_user, data_analysis_obj.time_stamp)]

        temp_series = pd.Series(meta_info_string_items, index=[list(output_df_absolute)[0]], name='meta_info_summary')
        output_df_absolute = output_df_absolute.append(temp_series)
        output_df_relative = output_df_relative.append(temp_series)
        for data_set_object in query_set_of_data_sets:
            data_set_meta_list = [
                'Data_set ID: {}; Data_set name: {}; submitting_user: {}; time_stamp: {}'.format(
                    data_set_object.id, data_set_object.name,
                    data_set_object.submitting_user, data_set_object.time_stamp)]

            temp_series = pd.Series(data_set_meta_list, index=[list(output_df_absolute)[0]], name='data_set_info')
            output_df_absolute = output_df_absolute.append(temp_series)
            output_df_relative = output_df_relative.append(temp_series)
    else:
        # call_type=='stand_alone'
        meta_info_string_items = [
            'Stand_alone output by {} on {}; Number of data_set objects as part of output = {}'.format(
                output_user, str(datetime.now()).replace(' ', '_').replace(':', '-'), len(query_set_of_data_sets))]

        temp_series = pd.Series(meta_info_string_items, index=[list(output_df_absolute)[0]], name='meta_info_summary')
        output_df_absolute = output_df_absolute.append(temp_series)
        output_df_relative = output_df_relative.append(temp_series)
        for data_set_object in query_set_of_data_sets:
            data_set_meta_list = [
                'Data_set ID: {}; Data_set name: {}; submitting_user: {}; time_stamp: {}'.format(
                    data_set_object.id, data_set_object.name, data_set_object.submitting_user,
                    data_set_object.time_stamp)]

            temp_series = pd.Series(data_set_meta_list, index=[list(output_df_absolute)[0]], name='data_set_info')
            output_df_absolute = output_df_absolute.append(temp_series)
            output_df_relative = output_df_relative.append(temp_series)

    # Here we have the tables populated and ready to output
    if not time_date_str:
        date_time_string = str(datetime.now()).replace(' ', '_').replace(':', '-')
    else:
        date_time_string = time_date_str
    if analysis_obj_id:
        data_analysis_obj = DataAnalysis.objects.get(id=analysis_obj_id)
        path_to_div_absolute = '{}/{}_{}_{}.seqs.absolute.txt'.format(output_dir, analysis_obj_id,
                                                                      data_analysis_obj.name, date_time_string)
        path_to_div_relative = '{}/{}_{}_{}.seqs.relative.txt'.format(output_dir, analysis_obj_id,
                                                                      data_analysis_obj.name, date_time_string)
        fasta_path = '{}/{}_{}_{}.seqs.fasta'.format(output_dir, analysis_obj_id,
                                                     data_analysis_obj.name, date_time_string)

    else:
        path_to_div_absolute = '{}/{}.seqs.absolute.txt'.format(output_dir, date_time_string)
        path_to_div_relative = '{}/{}.seqs.relative.txt'.format(output_dir, date_time_string)
        fasta_path = '{}/{}.seqs.fasta'.format(output_dir, date_time_string)

    os.makedirs(output_dir, exist_ok=True)
    output_df_absolute.to_csv(path_to_div_absolute, sep="\t")
    output_path_list.append(path_to_div_absolute)

    output_df_relative.to_csv(path_to_div_relative, sep="\t")
    output_path_list.append(path_to_div_relative)

    # we created the fasta above.
    write_list_to_destination(fasta_path, fasta_output_list)
    output_path_list.append(fasta_path)

    print('\nITS2 sequence output files:')
    for path_item in output_path_list:
        print(path_item)

    return output_path_list, date_time_string, len(sample_list)

def div_output_pre_analysis_new_meta_and_new_dss_structure_data_set_sample_id_input(
        data_set_sample_ids_to_output_string, num_processors, output_dir,
        sorted_sample_uid_list=None, analysis_obj_id=None, output_user=None, time_date_str=None):

    # ######################### ITS2 INTRA ABUND COUNT TABLE ################################
    # This is where we're going to have to work with the sequences that aren't part of a type.
    # Esentially we want to end up with the noName sequences divieded up cladally.
    # So at the moment this will mean that we will divide up the current no names but also
    # add information about the other cladal sequences that didn't make it into a cladeCollection

    # list to hold the paths of the outputted files
    output_path_list = []

    # ############### GET ORDERED LIST OF INTRAS BY CLADE THEN ABUNDANCE #################
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']


    data_set_samples_ids_to_output = [int(a) for a in data_set_sample_ids_to_output_string.split(',')]


    # Get collection of data_sets that are specific for the dataSubmissions we are looking at

    query_set_of_data_set_samples = DataSetSample.objects.filter(id__in=data_set_samples_ids_to_output)

    reference_sequences_in_data_set_samples = ReferenceSequence.objects.filter(
        datasetsamplesequence__data_set_sample_from__in=query_set_of_data_set_samples).distinct()



    # Get list of clades found
    clade_set = set()
    for ref_seq in     reference_sequences_in_data_set_samples:
        clade_set.add(ref_seq.clade)

    # Order the list of clades alphabetically
    # the only purpuse of this worker 2 is to end up with a set of sequences ordered by clade and then
    # by abundance across all samples.
    # I was doing this on a clade by clade basis on a sequence by sequence basis.
    # but now we should get rid of the clade sorting and just do that afterwards.
    # we should instad go by a sample by sample basis. This will allow us also to undertake the work that
    # is happening in worker three.
    # two birds one stone.
    # this is actually also better because now we can count the abundance of the sequences in proportions
    # rather than simply absolute values. Which is far more accurate.
    sub_clade_list = [a for a in clade_list if a in clade_set]

    sys.stdout.write('\n')

    # The manager that we will build all of the shared dictionaries from below.
    worker_manager = Manager()

    sample_list = DataSetSample.objects.filter(id__in=data_set_samples_ids_to_output)

    # Dictionary that will hold the list of data_set_sample_sequences for each sample
    sample_to_dsss_list_shared_dict = worker_manager.dict()
    print('Creating sample to data_set_sample_sequence dict:')
    for dss in sample_list:
        sys.stdout.write('\r{}'.format(dss.name))
        sample_to_dsss_list_shared_dict[dss.id] = list(
            DataSetSampleSequence.objects.filter(data_set_sample_from=dss))

    # Queue that will hold the data set samples for the MP
    data_set_sample_queue = Queue()

    # I will have a set of three dictionaries to pass into worker 2
    # 1 - Seqname to cumulative abundance of relative abundances (for getting the over lying order of ref seqs)
    # 2 - sample_id : list(dict(seq:abund), dict(seq:rel_abund))
    # 3 - sample_id : list(dict(noNameClade:abund), dict(noNameClade:rel_abund)

    reference_sequence_names_annotated = [
        ref_seq.name if ref_seq.has_name
        else str(ref_seq.id) + '_{}'.format(ref_seq.clade) for ref_seq in     reference_sequences_in_data_set_samples]

    generic_seq_to_abund_dict = {refSeq_name: 0 for refSeq_name in reference_sequence_names_annotated}

    list_of_dicts_for_processors = []
    for n in range(num_processors):
        list_of_dicts_for_processors.append(
            (worker_manager.dict(generic_seq_to_abund_dict), worker_manager.dict(), worker_manager.dict()))

    for dss in sample_list:
        data_set_sample_queue.put(dss)

    for N in range(num_processors):
        data_set_sample_queue.put('STOP')

    all_processes = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()

    for n in range(num_processors):
        p = Process(target=output_worker_two, args=(
            data_set_sample_queue, list_of_dicts_for_processors[n][0], list_of_dicts_for_processors[n][1],
            list_of_dicts_for_processors[n][2], reference_sequence_names_annotated, sample_to_dsss_list_shared_dict))
        all_processes.append(p)
        p.start()

    for p in all_processes:
        p.join()

    print('\nCollecting results of data_set_sample_counting across {} dictionaries'.format(num_processors))
    master_seq_abundance_counter = Counter()
    master_smple_seq_dict = dict()
    master_smple_no_name_clade_summary = dict()

    # now we need to do different actions for each of the three dictionary sets. One set for each num_proc
    # for the seqName counter we simply need to add to the master counter as we were doing before
    # for both of the sample-centric dictionaries we simply need to update a master dictionary
    for n in range(len(list_of_dicts_for_processors)):
        sys.stdout.write('\rdictionary {}(0)/{}'.format(n, num_processors))
        master_seq_abundance_counter += Counter(dict(list_of_dicts_for_processors[n][0]))
        sys.stdout.write('\rdictionary {}(1)/{}'.format(n, num_processors))
        master_smple_seq_dict.update(dict(list_of_dicts_for_processors[n][1]))
        sys.stdout.write('\rdictionary {}(2)/{}'.format(n, num_processors))
        master_smple_no_name_clade_summary.update(dict(list_of_dicts_for_processors[n][2]))

    print('Collection complete.')

    # we now need to separate by clade and sort within the clade
    clade_abundance_ordered_ref_seq_list = []
    for i in range(len(sub_clade_list)):
        temp_within_clade_list_for_sorting = []
        for seq_name, abund_val in master_seq_abundance_counter.items():
            if seq_name.startswith(sub_clade_list[i]) or seq_name[-2:] == '_{}'.format(sub_clade_list[i]):
                # then this is a seq of the clade in Q and we should add to the temp list
                temp_within_clade_list_for_sorting.append((seq_name, abund_val))
        # now sort the temp_within_clade_list_for_sorting and add to the cladeAbundanceOrderedRefSeqList
        sorted_within_clade = [
            a[0] for a in sorted(temp_within_clade_list_for_sorting, key=lambda x: x[1], reverse=True)]

        clade_abundance_ordered_ref_seq_list.extend(sorted_within_clade)

    # now delete the master_seq_abundance_counter as we are done with it
    del master_seq_abundance_counter

    # ##### WORKER THREE DOMAIN

    # we will eventually have the outputs stored in pandas dataframes.
    # in the worker below we will create a set of pandas.Series for each of the samples which will hold the abundances
    # one for the absoulte abundance and one for the relative abundances.

    # we will put together the headers piece by piece
    header_pre = clade_abundance_ordered_ref_seq_list
    no_name_summary_strings = ['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']]
    qc_stats = [
        'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs', 'post_taxa_id_absolute_symbiodinium_seqs',
        'post_taxa_id_unique_symbiodinium_seqs', 'size_screening_violation_absolute', 'size_screening_violation_unique',
        'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs', 'post_med_absolute',
        'post_med_unique']

    # append the noName sequences as individual sequence abundances
    output_header = ['sample_name'] + qc_stats + no_name_summary_strings + header_pre

    # #####################################################################################

    # ########### POPULATE TABLE WITH clade_collection_object DATA #############

    # In order to MP this we will have to pay attention to the order. As we can't keep the order as we work with the
    # MPs we will do as we did above for the profies outputs and have an output dict that will have the sample as key
    # and its corresponding data row as the value. Then when we have finished the MPing we will go through the
    # sampleList order and output the data in this order.
    # So we will need a managed output dict.

    worker_manager = Manager()
    managed_sample_output_dict = worker_manager.dict()

    data_set_sample_queue = Queue()

    for dss in sample_list:
        data_set_sample_queue.put(dss)

    for N in range(num_processors):
        data_set_sample_queue.put('STOP')

    all_processes = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()

    sys.stdout.write('\n\nOutputting seq data\n')
    for N in range(num_processors):
        p = Process(target=output_worker_three, args=(
            data_set_sample_queue, managed_sample_output_dict, clade_abundance_ordered_ref_seq_list, output_header,
            master_smple_seq_dict, master_smple_no_name_clade_summary))
        all_processes.append(p)
        p.start()

    for p in all_processes:
        p.join()

    print('\nseq output complete\n')

    managed_sample_output_dict_dict = dict(managed_sample_output_dict)

    # now we need to populate the output dataframe with the sample series in order of the sorted_sample_list
    # if there is one.
    # If there is a sorted sample list, make sure that it matches the samples that we are outputting

    # We were having an issue with data_sets having the same names. To fix this, we should do our ordering
    # accoring to the uids of the samples

    if sorted_sample_uid_list:
        sys.stdout.write('\nValidating sorted sample list and ordering dataframe accordingly\n')
        if len(sorted_sample_uid_list) != len(sample_list):
            sys.exit('Number of items in sorted_sample_list do not match those to be outputted!')
        if list(set(sorted_sample_uid_list).difference(set([s.id for s in sample_list]))):
            # then there is a sample that doesn't match up from the sorted_sample_uid_list that
            # has been passed in and the unordered sample list that we are working with in the code
            sys.exit('Sample list passed in does not match sample list from db query')

        # if we got to here then the sorted_sample_list looks good
        # I was originally performing the concat directly on the managedSampleOutputDict but this was starting
        # to produce errors. Starting to work on the managedSampleOutputDict_dict seems to not produce these
        # errors.
        # it may be a good idea to break this down to series by series instead of a one liner so that we can
        # print out progress
        # we can use the
        sys.stdout.write('\rPopulating the absolute dataframe with series. This could take a while...')
        output_df_absolute = pd.concat(
            [list_of_series[0] for list_of_series in managed_sample_output_dict_dict.values()], axis=1)
        sys.stdout.write('\rPopulating the relative dataframe with series. This could take a while...')
        output_df_relative = pd.concat(
            [list_of_series[1] for list_of_series in managed_sample_output_dict_dict.values()], axis=1)

        # now transpose
        output_df_absolute = output_df_absolute.T
        output_df_relative = output_df_relative.T

        # now make sure that the order is correct.
        output_df_absolute = output_df_absolute.reindex(sorted_sample_uid_list)
        output_df_relative = output_df_relative.reindex(sorted_sample_uid_list)

    else:

        # this returns a list which is simply the names of the samples
        # This will order the samples according to which sequence is their most abundant.
        # I.e. samples found to have the sequence which is most abundant in the largest number of sequences
        # will be first. Within each maj sequence, the samples will be sorted by the abundance of that sequence
        # in the sample.
        # At the moment we are also ordering by clade just so that you see samples with the A's at the top
        # of the output so that we minimise the number of 0's in the top left of the output
        # honestly I think we could perhaps get rid of this and just use the over all abundance of the sequences
        # discounting clade. THis is what we do for the clade order when plotting.
        sys.stdout.write('\nGenerating ordered sample list and ordering dataframe accordingly\n')
        ordered_sample_list_by_uid = generate_ordered_sample_list(managed_sample_output_dict_dict)

        # if we got to here then the sorted_sample_list looks good
        # I was originally performing the concat directly on the managedSampleOutputDict but this was starting
        # to produce errors. Starting to work on the managedSampleOutputDict_dict seems to not produce these
        # errors.
        sys.stdout.write('\rPopulating the absolute dataframe with series. This could take a while...')
        output_df_absolute = pd.concat(
            [list_of_series[0] for list_of_series in managed_sample_output_dict_dict.values()], axis=1)
        sys.stdout.write('\rPopulating the relative dataframe with series. This could take a while...')
        output_df_relative = pd.concat(
            [list_of_series[1] for list_of_series in managed_sample_output_dict_dict.values()], axis=1)

        # now transpose
        sys.stdout.write('\rTransposing...')
        output_df_absolute = output_df_absolute.T
        output_df_relative = output_df_relative.T

        # now make sure that the order is correct.
        sys.stdout.write('\rReordering index...')
        output_df_absolute = output_df_absolute.reindex(ordered_sample_list_by_uid)
        output_df_relative = output_df_relative.reindex(ordered_sample_list_by_uid)

    # when adding the accession numbers below, we have to go through every sequence and look up its object
    # We also have to do this when we are outputting the fasta.
    # to prevent us having to make every look up twice, we should also make the fasta at the same time
    # Output a .fasta for of all of the sequences found in the analysis
    # we will write out the fasta right at the end.
    fasta_output_list = []

    # Now add the accesion number / UID for each of the DIVs
    sys.stdout.write('\nGenerating accession and fasta\n')

    # go column name by column name and if the col name is in seq_annotated_name
    # then get the accession and add to the accession_list
    # else do nothing and a blank should be automatically added for us.
    # This was painfully slow because we were doing individual calls to the dictionary
    # I think this will be much faster if do two queries of the db to get the named and
    # non named refseqs and then make two dicts for each of these and use these to populate the below
    reference_sequences_in_data_sets_no_name = ReferenceSequence.objects.filter(
        datasetsamplesequence__data_set_sample_from__in=query_set_of_data_set_samples,
        has_name=False).distinct()
    reference_sequences_in_data_sets_has_name = ReferenceSequence.objects.filter(
        datasetsamplesequence__data_set_sample_from__in=query_set_of_data_set_samples,
        has_name=True).distinct()
    # no name dict should be a dict of id to sequence
    no_name_dict = {rs.id: rs.sequence for rs in reference_sequences_in_data_sets_no_name}
    # has name dict should be a dict of name to sequence
    has_name_dict = {rs.name: (rs.id, rs.sequence) for rs in reference_sequences_in_data_sets_has_name}

    # for the time being we are going to ignore whether a refseq has an assession as we have not put this
    # into use yet.
    accession_list = []
    num_cols = len(list(output_df_relative))
    for i, col_name in enumerate(list(output_df_relative)):
        sys.stdout.write('\rAppending accession info and creating fasta {}: {}/{}'.format(col_name, i, num_cols))
        if col_name in clade_abundance_ordered_ref_seq_list:
            if col_name[-2] == '_':
                col_name_id = int(col_name[:-2])
                accession_list.append(str(col_name_id))
                fasta_output_list.append('>{}'.format(col_name))
                fasta_output_list.append(no_name_dict[col_name_id])
            else:
                col_name_tup = has_name_dict[col_name]
                accession_list.append(str(col_name_tup[0]))
                fasta_output_list.append('>{}'.format(col_name))
                fasta_output_list.append(col_name_tup[1])
        else:
            accession_list.append(np.nan)

    temp_series = pd.Series(accession_list, name='seq_accession', index=list(output_df_relative))
    output_df_absolute = output_df_absolute.append(temp_series)
    output_df_relative = output_df_relative.append(temp_series)

    # Now append the meta infromation for each of the data_sets that make up the output contents
    # this is information like the submitting user, what the uids of the datasets are etc.
    # There are several ways that this can be called.
    # it can be called as part of the submission: call_type = submission
    # part of an analysis output: call_type = analysis
    # or stand alone: call_type = 'stand_alone'
    # we should have an output for each scenario


    # call_type=='stand_alone' call type will always be standalone when we are doing a specific data_set_sample output
    meta_info_string_items = [
        'Stand alone output by {} on {}; Number of DataSetSample objects as part of output = {}'.format(
            output_user, str(datetime.now()).replace(' ', '_').replace(':', '-'), len(query_set_of_data_set_samples))]

    temp_series = pd.Series(meta_info_string_items, index=[list(output_df_absolute)[0]], name='meta_info_summary')
    output_df_absolute = output_df_absolute.append(temp_series)
    output_df_relative = output_df_relative.append(temp_series)

    data_sets_of_the_data_set_samples_string = DataSet.objects.filter(datasetsample__in=query_set_of_data_set_samples).distinct()

    for data_set_object in data_sets_of_the_data_set_samples_string:
        # query of the DataSet objects the DataSetSamples are from
        data_set_meta_list = [
            'Data_set ID: {}; Data_set name: {}; submitting_user: {}; time_stamp: {}'.format(
                data_set_object.id, data_set_object.name, data_set_object.submitting_user,
                data_set_object.time_stamp)]

        temp_series = pd.Series(data_set_meta_list, index=[list(output_df_absolute)[0]], name='data_set_info')
        output_df_absolute = output_df_absolute.append(temp_series)
        output_df_relative = output_df_relative.append(temp_series)

    # Here we have the tables populated and ready to output
    if not time_date_str:
        date_time_string = str(datetime.now()).replace(' ', '_').replace(':', '-')
    else:
        date_time_string = time_date_str
    if analysis_obj_id:
        data_analysis_obj = DataAnalysis.objects.get(id=analysis_obj_id)
        path_to_div_absolute = '{}/{}_{}_{}.seqs.absolute.txt'.format(output_dir, analysis_obj_id,
                                                                      data_analysis_obj.name, date_time_string)
        path_to_div_relative = '{}/{}_{}_{}.seqs.relative.txt'.format(output_dir, analysis_obj_id,
                                                                      data_analysis_obj.name, date_time_string)
        fasta_path = '{}/{}_{}_{}.seqs.fasta'.format(output_dir, analysis_obj_id,
                                                     data_analysis_obj.name, date_time_string)

    else:
        path_to_div_absolute = '{}/{}.seqs.absolute.txt'.format(output_dir, date_time_string)
        path_to_div_relative = '{}/{}.seqs.relative.txt'.format(output_dir, date_time_string)
        fasta_path = '{}/{}.seqs.fasta'.format(output_dir, date_time_string)

    os.makedirs(output_dir, exist_ok=True)
    output_df_absolute.to_csv(path_to_div_absolute, sep="\t")
    output_path_list.append(path_to_div_absolute)

    output_df_relative.to_csv(path_to_div_relative, sep="\t")
    output_path_list.append(path_to_div_relative)

    # we created the fasta above.
    write_list_to_destination(fasta_path, fasta_output_list)
    output_path_list.append(fasta_path)

    print('\nITS2 sequence output files:')
    for path_item in output_path_list:
        print(path_item)

    return output_path_list, date_time_string, len(sample_list)


def output_worker_two(input_queue, seq_rel_abund_dict, smpl_seq_dict, sample_no_name_clade_summary_dict,
                      reference_sequence_names_annotated, sample_to_dsss_list_shared_dict):
    # 1 - Seqname to cumulative relative abundance for each sequence across all sampples (for getting the over lying order of ref seqs)
    # 2 - sample_id : list(dict(ref_seq_of_sample_name:absolute_abundance_of_dsss_in_sample), dict(ref_seq_of_sample_name:relative_abundance_of_dsss_in_sample))
    # 3 - sample_id : list(dict(clade:total_abund_of_no_name_seqs_of_clade_in_q_), dict(clade:relative_abund_of_no_name_seqs_of_clade_in_q_)
    for dss in iter(input_queue.get, 'STOP'):
        sys.stdout.write('\rCounting seqs for {}'.format(dss))

        cladal_abundances = [int(a) for a in json.loads(dss.cladal_seq_totals)]

        sample_seq_tot = sum(cladal_abundances)

        # the first dict will hold the absolute abundances, whilst the second will hold the relative abundances
        clade_summary_absolute_dict = {clade: 0 for clade in list('ABCDEFGHI')}
        clade_summary_relative_dict = {clade: 0 for clade in list('ABCDEFGHI')}
        smple_seq_count_aboslute_dict = {seq_name: 0 for seq_name in reference_sequence_names_annotated}
        smple_seq_count_relative_dict = {seq_name: 0 for seq_name in reference_sequence_names_annotated}

        dsss_in_sample = sample_to_dsss_list_shared_dict[dss.id]

        for dsss in dsss_in_sample:
            # determine what the name of the seq will be in the output
            if not dsss.reference_sequence_of.has_name:
                name_unit = str(dsss.reference_sequence_of.id) + '_{}'.format(dsss.reference_sequence_of.clade)
                # the clade summries are only for the noName seqs
                clade_summary_absolute_dict[dsss.reference_sequence_of.clade] += dsss.abundance
                clade_summary_relative_dict[dsss.reference_sequence_of.clade] += dsss.abundance / sample_seq_tot
            else:
                name_unit = dsss.reference_sequence_of.name

            seq_rel_abund_dict[name_unit] += dsss.abundance / sample_seq_tot
            smple_seq_count_aboslute_dict[name_unit] += dsss.abundance
            smple_seq_count_relative_dict[name_unit] += dsss.abundance / sample_seq_tot

        sample_no_name_clade_summary_dict[dss.id] = [clade_summary_absolute_dict, clade_summary_relative_dict]
        smpl_seq_dict[dss.id] = [smple_seq_count_aboslute_dict, smple_seq_count_relative_dict]


def generate_ordered_sample_list(managed_sample_output_dict):
    # create a df from the managedSampleOutputDict. We will use the relative values here

    output_df_relative = pd.concat([list_of_series[1] for list_of_series in managed_sample_output_dict.values()],
                                   axis=1)
    output_df_relative = output_df_relative.T

    # now remove the rest of the non abundance columns
    non_seq_columns = [
        'sample_name', 'raw_contigs', 'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_qc_absolute_seqs',
        'post_qc_unique_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs', 'post_taxa_id_absolute_symbiodinium_seqs',
        'post_taxa_id_unique_symbiodinium_seqs', 'post_med_absolute', 'post_med_unique',
        'size_screening_violation_absolute', 'size_screening_violation_unique']

    no_name_seq_columns = ['noName Clade {}'.format(clade) for clade in list('ABCDEFGHI')]
    cols_to_drop = non_seq_columns + no_name_seq_columns

    output_df_relative.drop(columns=cols_to_drop, inplace=True)
    ordered_sample_list = get_sample_order_from_rel_seq_abund_df(output_df_relative)
    return ordered_sample_list


def get_sample_order_from_rel_seq_abund_df(sequence_only_df_relative):
    max_seq_ddict = defaultdict(int)
    seq_to_samp_dict = defaultdict(list)

    # for each sample get the columns name of the max value of a div
    no_maj_samps = []
    for sample_to_sort_ID in sequence_only_df_relative.index.values.tolist():
        sys.stdout.write('\rGetting maj seq for sample {}'.format(sample_to_sort_ID))
        series_as_float = sequence_only_df_relative.loc[sample_to_sort_ID].astype('float')
        max_rel_abund = series_as_float.max()

        if not max_rel_abund > 0:
            no_maj_samps.append(sample_to_sort_ID)
        else:
            max_abund_seq = series_as_float.idxmax()
            # add a tup of sample name and rel abund of seq to the seq_to_samp_dict
            seq_to_samp_dict[max_abund_seq].append((sample_to_sort_ID, max_rel_abund))
            # add this to the ddict count
            max_seq_ddict[max_abund_seq] += 1

    # then once we have compelted this for all sequences go clade by clade
    # and generate the sample order
    ordered_sample_list_by_uid = []
    sys.stdout.write('\nGoing clade by clade sorting by abundance\n')
    for clade in list('ABCDEFGHI'):
        sys.stdout.write('\rGetting clade {} seqs'.format(clade))
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

        ordered_sequence_of_clade_list = [x[0] for x in sorted(tup_list_of_clade, key=lambda x: x[1], reverse=True)]

        for seq_to_order_samples_by in ordered_sequence_of_clade_list:
            sys.stdout.write('\r{}'.format(seq_to_order_samples_by))
            tup_list_of_samples_that_had_sequence_as_most_abund = seq_to_samp_dict[seq_to_order_samples_by]
            ordered_list_of_samples_for_seq_ordered = \
                [x[0] for x in
                 sorted(tup_list_of_samples_that_had_sequence_as_most_abund, key=lambda x: x[1], reverse=True)]
            ordered_sample_list_by_uid.extend(ordered_list_of_samples_for_seq_ordered)
    # finally add in the samples that didn't have a maj sequence
    ordered_sample_list_by_uid.extend(no_maj_samps)
    return ordered_sample_list_by_uid


def output_worker_three(
        input_queue, out_dict, clade_abundance_ordered_ref_seq_list, output_header,
        smpl_abund_dicts_dict, smpl_clade_summary_dicts_dict):

    clade_list = list('ABCDEFGHI')
    for dss in iter(input_queue.get, 'STOP'):

        sys.stdout.write('\rOutputting seq data for {}'.format(dss.name))
        # List that will hold the row
        sample_row_data_counts = []
        sample_row_data_props = []
        cladal_abundances = [int(a) for a in json.loads(dss.cladal_seq_totals)]
        sample_seq_tot = sum(cladal_abundances)

        if dss.error_in_processing or sample_seq_tot == 0:
            # Then this sample had a problem in the sequencing and we need to just output 0s across the board
            # QC

            # Append the name of the dss for when we have samples of the same name
            sample_row_data_counts.append(dss.name)
            sample_row_data_props.append(dss.name)

            populate_quality_control_data_of_failed_sample(dss, sample_row_data_counts, sample_row_data_props)

            # no name clade summaries get 0.
            for _ in clade_list:
                sample_row_data_counts.append(0)
                sample_row_data_props.append(0)

            # All sequences get 0s
            for _ in clade_abundance_ordered_ref_seq_list:
                sample_row_data_counts.append(0)
                sample_row_data_props.append(0)

            # Here we need to add the string to the output_dictionary rather than the intraAbund table objects
            sample_series_absolute = pd.Series(sample_row_data_counts, index=output_header, name=dss.id)
            sample_series_relative = pd.Series(sample_row_data_counts, index=output_header, name=dss.id)

            out_dict[dss.id] = [sample_series_absolute, sample_series_relative]
            continue

        # get the list of sample specific dicts that contain the clade summaries and the seq abundances
        smpl_seq_abund_absolute_dict = smpl_abund_dicts_dict[dss.id][0]
        smpl_seq_abund_relative_dict = smpl_abund_dicts_dict[dss.id][1]
        smpl_clade_summary_absolute_dict = smpl_clade_summary_dicts_dict[dss.id][0]
        smpl_clade_summary_relative_dict = smpl_clade_summary_dicts_dict[dss.id][1]

        # Here we add in the post qc and post-taxa id counts
        # For the absolute counts we will report the absolute seq number
        # For the relative counts we will report these as proportions of the sampleSeqTot.
        # I.e. we will have numbers larger than 1 for many of the values and the symbiodinium seqs should be 1

        # Append the name of the dss for when we have samples of the same name
        sample_row_data_counts.append(dss.name)
        sample_row_data_props.append(dss.name)

        populate_quality_control_data_of_successful_sample(
            dss, sample_row_data_counts, sample_row_data_props, sample_seq_tot)

        # now add the clade divided summaries of the clades
        for clade in clade_list:
            sys.stdout.write('\rOutputting seq data for {}: clade summary {}'.format(dss.name, clade))
            sample_row_data_counts.append(smpl_clade_summary_absolute_dict[clade])
            sample_row_data_props.append(smpl_clade_summary_relative_dict[clade])

        # and append these abundances in order of cladeAbundanceOrderedRefSeqList to
        # the sampleRowDataCounts and the sampleRowDataProps
        for seq_name in clade_abundance_ordered_ref_seq_list:
            sys.stdout.write('\rOutputting seq data for {}: sequence {}'.format(dss.name, seq_name))
            sample_row_data_counts.append(smpl_seq_abund_absolute_dict[seq_name])
            sample_row_data_props.append(smpl_seq_abund_relative_dict[seq_name])

        # Here we need to add the string to the output_dictionary rather than the intraAbund table objects
        sample_series_absolute = pd.Series(sample_row_data_counts, index=output_header, name=dss.id)
        sample_series_relative = pd.Series(sample_row_data_props, index=output_header, name=dss.id)

        out_dict[dss.id] = [sample_series_absolute, sample_series_relative]


def populate_quality_control_data_of_successful_sample(
        dss, sample_row_data_counts, sample_row_data_props, sample_seq_tot):
    # CONTIGS
    # This is the absolute number of sequences after make.contigs
    contig_num = dss.num_contigs
    sample_row_data_counts.append(contig_num)
    sample_row_data_props.append(contig_num / sample_seq_tot)
    # POST-QC
    # store the aboslute number of sequences after sequencing QC at this stage
    post_qc_absolute = dss.post_qc_absolute_num_seqs
    sample_row_data_counts.append(post_qc_absolute)
    sample_row_data_props.append(post_qc_absolute / sample_seq_tot)
    # This is the unique number of sequences after the sequencing QC
    post_qc_unique = dss.post_qc_unique_num_seqs
    sample_row_data_counts.append(post_qc_unique)
    sample_row_data_props.append(post_qc_unique / sample_seq_tot)
    # POST TAXA-ID
    # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
    tax_id_symbiodinium_absolute = dss.absolute_num_sym_seqs
    sample_row_data_counts.append(tax_id_symbiodinium_absolute)
    sample_row_data_props.append(tax_id_symbiodinium_absolute / sample_seq_tot)
    # Same as above but the number of unique seqs
    tax_id_symbiodinium_unique = dss.unique_num_sym_seqs
    sample_row_data_counts.append(tax_id_symbiodinium_unique)
    sample_row_data_props.append(tax_id_symbiodinium_unique / sample_seq_tot)
    # store the absolute number of sequences lost to size cutoff violations
    size_violation_aboslute = dss.size_violation_absolute
    sample_row_data_counts.append(size_violation_aboslute)
    sample_row_data_props.append(size_violation_aboslute / sample_seq_tot)
    # store the unique size cutoff violations
    size_violation_unique = dss.size_violation_unique
    sample_row_data_counts.append(size_violation_unique)
    sample_row_data_props.append(size_violation_unique / sample_seq_tot)
    # store the abosolute number of sequenes that were not considered Symbiodinium
    tax_id_non_symbiodinum_abosulte = dss.non_sym_absolute_num_seqs
    sample_row_data_counts.append(tax_id_non_symbiodinum_abosulte)
    sample_row_data_props.append(tax_id_non_symbiodinum_abosulte / sample_seq_tot)
    # This is the number of unique sequences that were not considered Symbiodinium
    tax_id_non_symbiodinium_unique = dss.non_sym_unique_num_seqs
    sample_row_data_counts.append(tax_id_non_symbiodinium_unique)
    sample_row_data_props.append(tax_id_non_symbiodinium_unique / sample_seq_tot)
    # Post MED absolute
    post_med_absolute = dss.post_med_absolute
    sample_row_data_counts.append(post_med_absolute)
    sample_row_data_props.append(post_med_absolute / sample_seq_tot)
    # Post MED unique
    post_med_unique = dss.post_med_unique
    sample_row_data_counts.append(post_med_unique)
    sample_row_data_props.append(post_med_unique / sample_seq_tot)


def populate_quality_control_data_of_failed_sample(dss, sample_row_data_counts, sample_row_data_props):
    # Add in the qc totals if possible
    # For the proportions we will have to add zeros as we cannot do proportions
    # CONTIGS
    # This is the absolute number of sequences after make.contigs

    if dss.num_contigs:
        contig_num = dss.num_contigs
        sample_row_data_counts.append(contig_num)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # POST-QC
    # store the aboslute number of sequences after sequencing QC at this stage
    if dss.post_qc_absolute_num_seqs:
        post_qc_absolute = dss.post_qc_absolute_num_seqs
        sample_row_data_counts.append(post_qc_absolute)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # This is the unique number of sequences after the sequencing QC
    if dss.post_qc_unique_num_seqs:
        post_qc_unique = dss.post_qc_unique_num_seqs
        sample_row_data_counts.append(post_qc_unique)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # POST TAXA-ID
    # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
    if dss.absolute_num_sym_seqs:
        tax_id_symbiodinium_absolute = dss.absolute_num_sym_seqs
        sample_row_data_counts.append(tax_id_symbiodinium_absolute)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # Same as above but the number of unique seqs
    if dss.unique_num_sym_seqs:
        tax_id_symbiodinium_unique = dss.unique_num_sym_seqs
        sample_row_data_counts.append(tax_id_symbiodinium_unique)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # size violation absolute
    if dss.size_violation_absolute:
        size_viol_ab = dss.size_violation_absolute
        sample_row_data_counts.append(size_viol_ab)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # size violation unique
    if dss.size_violation_unique:
        size_viol_uni = dss.size_violation_unique
        sample_row_data_counts.append(size_viol_uni)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # store the abosolute number of sequenes that were not considered Symbiodinium
    if dss.non_sym_absolute_num_seqs:
        tax_id_non_symbiodinum_abosulte = dss.non_sym_absolute_num_seqs
        sample_row_data_counts.append(tax_id_non_symbiodinum_abosulte)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # This is the number of unique sequences that were not considered Symbiodinium
    if dss.non_sym_unique_num_seqs:
        tax_id_non_symbiodinium_unique = dss.non_sym_unique_num_seqs
        sample_row_data_counts.append(tax_id_non_symbiodinium_unique)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # post-med absolute
    if dss.post_med_absolute:
        post_med_abs = dss.post_med_absolute
        sample_row_data_counts.append(post_med_abs)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # post-med absolute
    if dss.post_med_unique:
        post_med_uni = dss.post_med_unique
        sample_row_data_counts.append(post_med_uni)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)

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
            self, symportal_root_dir, call_type, num_proc, dss_uids_output_str=None, ds_uids_output_str=None, output_dir=None,
            sorted_sample_uid_list=None, analysis_obj_id=None, time_date_str=None, output_user=None):
        self._init_core_vars(
            symportal_root_dir, analysis_obj_id, call_type, dss_uids_output_str, ds_uids_output_str, num_proc,
            output_dir, output_user, sorted_sample_uid_list, time_date_str)
        self._init_seq_abundance_collection_objects()
        self._init_vars_for_putting_together_the_dfs()
        self._init_output_paths()

    def _init_core_vars(self, symportal_root_dir, analysis_obj_id, call_type, dss_uids_output_str, ds_uids_output_str, num_proc,
                        output_dir, output_user, sorted_sample_uid_list, time_date_str):
        self._check_either_dss_or_dsss_uids_provided(dss_uids_output_str, ds_uids_output_str)
        if dss_uids_output_str:
            self.list_of_dss_objects = DataSetSample.objects.filter(id__in=[int(a) for a in dss_uids_output_str.split(',')])
            self.ds_objs_to_output = DataSet.objects.filter(datasetsample__in=self.list_of_dss_objects).distinct()
        elif ds_uids_output_str:
            uids_of_data_sets_to_output = [int(a) for a in ds_uids_output_str.split(',')]
            self.ds_objs_to_output = DataSet.objects.filter(id__in=uids_of_data_sets_to_output)
            self.list_of_dss_objects = DataSetSample.objects.filter(data_submission_from__in=self.ds_objs_to_output)

        self.ref_seqs_in_datasets = ReferenceSequence.objects.filter(
            datasetsamplesequence__data_set_sample_from__in=
            self.list_of_dss_objects).distinct()

        self.num_proc = num_proc
        self._set_output_dir(call_type, ds_uids_output_str, output_dir, symportal_root_dir)
        self.sorted_sample_uid_list = sorted_sample_uid_list
        self.analysis_obj_id = analysis_obj_id
        if time_date_str:
            self.time_date_str = time_date_str
        else:
            self.time_date_str = str(datetime.now()).replace(' ', '_').replace(':', '-')
        self.call_type = call_type
        self.output_user = output_user
        self.clade_list = list('ABCDEFGHI')


        set_of_clades_found = {ref_seq.clade for ref_seq in self.ref_seqs_in_datasets}
        self.ordered_list_of_clades_found = [clade for clade in self.clade_list if clade in set_of_clades_found]


    @staticmethod
    def _check_either_dss_or_dsss_uids_provided(data_set_sample_ids_to_output_string, data_set_uids_to_output_as_comma_sep_string):
        if data_set_sample_ids_to_output_string is not None and data_set_uids_to_output_as_comma_sep_string is not None:
            raise RuntimeError('Provide either dss uids or ds uids for outputing sequence count tables')

    def _set_output_dir(self, call_type, data_set_uids_to_output_as_comma_sep_string, output_dir, symportal_root_dir):
        if call_type == 'submission':
            self.output_dir = os.path.abspath(os.path.join(
                symportal_root_dir, 'outputs', 'data_set_submissions', data_set_uids_to_output_as_comma_sep_string))
        elif call_type == 'stand_alone':
            self.output_dir = os.path.abspath(os.path.join(symportal_root_dir, 'outputs', 'non_analysis'))
        else:  # call_type == 'analysis
            self.output_dir = output_dir
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
        self.output_df_absolute = None
        self.output_df_relative = None
        self.output_seqs_fasta_as_list = []

    def _init_output_paths(self):
        self.output_paths_list = []
        if self.analysis_obj_id:
            data_analysis_obj = DataAnalysis.objects.get(id=self.analysis_obj_id)
            self.path_to_seq_output_df_absolute = os.path.join(
                self.output_dir,
                f'{self.analysis_obj_id}_{data_analysis_obj.name}_{self.time_date_str}.seqs.absolute.txt')
            self.path_to_seq_output_df_relative = os.path.join(
                self.output_dir,
                f'{self.analysis_obj_id}_{data_analysis_obj.name}_{self.time_date_str}.seqs.relative.txt')

            self.output_fasta_path = os.path.join(
                self.output_dir, f'{self.analysis_obj_id}_{data_analysis_obj.name}_{self.time_date_str}.seqs.fasta')

        else:
            self.path_to_seq_output_df_absolute = os.path.join(self.output_dir,
                                                               f'{self.time_date_str}.seqs.absolute.txt')
            self.path_to_seq_output_df_relative = os.path.join(self.output_dir,
                                                               f'{self.time_date_str}.seqs.relative.txt')
            self.output_fasta_path = os.path.join(self.output_dir, f'{self.time_date_str}.seqs.fasta')

    def make_output_tables(self):
        self._collect_abundances_for_creating_the_output()

        self._generate_sample_output_series()

        self._create_ordered_output_dfs_from_series()

        self._add_uids_for_seqs_to_dfs()

        self._append_meta_info_to_df()

        self._write_out_dfs_and_fasta()

    def _write_out_dfs_and_fasta(self):
        self.output_df_absolute.to_csv(self.path_to_seq_output_df_absolute, sep="\t")
        self.output_paths_list.append(self.path_to_seq_output_df_absolute)
        self.output_df_relative.to_csv(self.path_to_seq_output_df_relative, sep="\t")
        self.output_paths_list.append(self.path_to_seq_output_df_relative)
        # we created the fasta above.
        write_list_to_destination(self.output_fasta_path, self.output_seqs_fasta_as_list)
        self.output_paths_list.append(self.output_fasta_path)
        print('\n\nITS2 sequence output files:')
        for path_item in self.output_paths_list:
            print(path_item)

    def _append_meta_info_to_df(self):
        # Now append the meta infromation for each of the data_sets that make up the output contents
        # this is information like the submitting user, what the uids of the datasets are etc.
        # There are several ways that this can be called.
        # it can be called as part of the submission: call_type = submission
        # part of an analysis output: call_type = analysis
        # or stand alone: call_type = 'stand_alone'
        # we should have an output for each scenario
        if self.call_type == 'submission':
            self._append_meta_info_to_df_submission()
        elif self.call_type == 'analysis':
            self._append_meta_info_to_df_analysis()
        else:
            # call_type=='stand_alone'
            self._append_meta_info_to_df_stand_alone()

    def _append_meta_info_to_df_submission(self):
        data_set_object = self.ds_objs_to_output[0]
        # there will only be one data_set object
        meta_info_string_items = [
            f'Output as part of data_set submission ID: {data_set_object.id}; '
            f'submitting_user: {data_set_object.submitting_user}; '
            f'time_stamp: {data_set_object.time_stamp}']
        temp_series = pd.Series(meta_info_string_items, index=[list(self.output_df_absolute)[0]],
                                name='meta_info_summary')
        self.output_df_absolute = self.output_df_absolute.append(temp_series)
        self.output_df_relative = self.output_df_relative.append(temp_series)

    def _append_meta_info_to_df_analysis(self):
        data_analysis_obj = DataAnalysis.objects.get(id=self.analysis_obj_id)
        num_data_set_objects_as_part_of_analysis = len(data_analysis_obj.list_of_data_set_uids.split(','))
        meta_info_string_items = [
            f'Output as part of data_analysis ID: {data_analysis_obj.id}; '
            f'Number of data_set objects as part of analysis = {num_data_set_objects_as_part_of_analysis}; '
            f'submitting_user: {data_analysis_obj.submitting_user}; time_stamp: {data_analysis_obj.time_stamp}']
        temp_series = pd.Series(meta_info_string_items, index=[list(self.output_df_absolute)[0]],
                                name='meta_info_summary')
        self.output_df_absolute = self.output_df_absolute.append(temp_series)
        self.output_df_relative = self.output_df_relative.append(temp_series)
        for data_set_object in self.ds_objs_to_output:
            data_set_meta_list = [
                f'Data_set ID: {data_set_object.id}; '
                f'Data_set name: {data_set_object.name}; '
                f'submitting_user: {data_set_object.submitting_user}; '
                f'time_stamp: {data_set_object.time_stamp}']

            temp_series = pd.Series(data_set_meta_list, index=[list(self.output_df_absolute)[0]], name='data_set_info')
            self.output_df_absolute = self.output_df_absolute.append(temp_series)
            self.output_df_relative = self.output_df_relative.append(temp_series)

    def _append_meta_info_to_df_stand_alone(self):
        meta_info_string_items = [
            f'Stand_alone output by {self.output_user} on {self.time_date_str}; '
            f'Number of data_set objects as part of output = {len(self.ds_objs_to_output)}']
        temp_series = pd.Series(meta_info_string_items, index=[list(self.output_df_absolute)[0]],
                                name='meta_info_summary')
        self.output_df_absolute = self.output_df_absolute.append(temp_series)
        self.output_df_relative = self.output_df_relative.append(temp_series)
        for data_set_object in self.ds_objs_to_output:
            data_set_meta_list = [
                f'Data_set ID: {data_set_object.id}; '
                f'Data_set name: {data_set_object.name}; '
                f'submitting_user: {data_set_object.submitting_user}; '
                f'time_stamp: {data_set_object.time_stamp}']
            temp_series = pd.Series(data_set_meta_list, index=[list(self.output_df_absolute)[0]], name='data_set_info')
            self.output_df_absolute = self.output_df_absolute.append(temp_series)
            self.output_df_relative = self.output_df_relative.append(temp_series)

    def _add_uids_for_seqs_to_dfs(self):
        """Now add the UID for each of the sequences"""
        sys.stdout.write('\nGenerating accession and fasta\n')
        reference_sequences_in_data_sets_no_name = ReferenceSequence.objects.filter(
            datasetsamplesequence__data_set_sample_from__in=self.list_of_dss_objects,
            has_name=False).distinct()
        reference_sequences_in_data_sets_has_name = ReferenceSequence.objects.filter(
            datasetsamplesequence__data_set_sample_from__in=self.list_of_dss_objects,
            has_name=True).distinct()
        no_name_dict = {rs.id: rs.sequence for rs in reference_sequences_in_data_sets_no_name}
        has_name_dict = {rs.name: (rs.id, rs.sequence) for rs in reference_sequences_in_data_sets_has_name}
        accession_list = []
        num_cols = len(list(self.output_df_relative))
        for i, col_name in enumerate(list(self.output_df_relative)):
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
        temp_series = pd.Series(accession_list, name='seq_accession', index=list(self.output_df_relative))
        self.output_df_absolute = self.output_df_absolute.append(temp_series)
        self.output_df_relative = self.output_df_relative.append(temp_series)

    def _create_ordered_output_dfs_from_series(self):
        """Put together the pandas series that hold sequences abundance outputs for each sample in order of the samples
        either according to a predefined ordered list or by an order that will be generated below."""
        if self.sorted_sample_uid_list:
            sys.stdout.write('\nValidating sorted sample list and ordering dataframe accordingly\n')
            self._check_sorted_sample_list_is_valid()

            self._create_ordered_output_dfs_from_series_with_sorted_sample_list()

        else:
            sys.stdout.write('\nGenerating ordered sample list and ordering dataframe accordingly\n')
            self.sorted_sample_uid_list = self._generate_ordered_sample_list()

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
        non_seq_columns = [
            'sample_name', 'raw_contigs', 'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_qc_absolute_seqs',
            'post_qc_unique_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
            'post_taxa_id_absolute_symbiodinium_seqs',
            'post_taxa_id_unique_symbiodinium_seqs', 'post_med_absolute', 'post_med_unique',
            'size_screening_violation_absolute', 'size_screening_violation_unique']
        no_name_seq_columns = ['noName Clade {}'.format(clade) for clade in list('ABCDEFGHI')]
        cols_to_drop = non_seq_columns + no_name_seq_columns
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

            ordered_sequence_of_clade_list = [x[0] for x in sorted(tup_list_of_clade, key=lambda x: x[1], reverse=True)]

            for seq_to_order_samples_by in ordered_sequence_of_clade_list:
                sys.stdout.write('\r{}'.format(seq_to_order_samples_by))
                tup_list_of_samples_that_had_sequence_as_most_abund = seq_to_samp_ddict[seq_to_order_samples_by]
                ordered_list_of_samples_for_seq_ordered = \
                    [x[0] for x in
                     sorted(tup_list_of_samples_that_had_sequence_as_most_abund, key=lambda x: x[1], reverse=True)]
                ordered_sample_list_by_uid.extend(ordered_list_of_samples_for_seq_ordered)
        # finally add in the samples that didn't have a maj sequence
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
        self.output_df_absolute = output_df_absolute.reindex(self.sorted_sample_uid_list)
        self.output_df_relative = output_df_relative.reindex(self.sorted_sample_uid_list)

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

        self.parent = parent_seq_count_tab_creator
        self.mp_manager = Manager()
        self.input_dss_mp_queue = Queue()
        self._populate_input_dss_mp_queue()
        self.ref_seq_names_clade_annotated = [
        ref_seq.name if ref_seq.has_name else str(ref_seq.id) + '_{}'.format(ref_seq.clade) for
            ref_seq in self.parent.ref_seqs_in_datasets]

        # TODO we were previously creating an MP dictionary for every proc used. We were then collecting them afterwards
        # I'm not sure if there was a good reason for doing this, but I don't see any comments to the contrary.
        # it should not be necessary to have a dict for every proc. Instead we can just have on mp dict.
        # we should check that this is still working as expected.
        # self.list_of_dictionaries_for_processes = self._generate_list_of_dicts_for_processes()
        self.dss_id_to_list_of_dsss_objects_mp_dict = self.mp_manager.dict()
        self._populate_dss_id_to_list_of_dsss_objects()
        self.annotated_dss_name_to_cummulative_rel_abund_mp_dict = self.mp_manager.dict(
            {refSeq_name: 0 for refSeq_name in self.ref_seq_names_clade_annotated})
        self.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict = self.mp_manager.dict()
        self.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict = self.mp_manager.dict()

        # this is the list that we will use the self.annotated_dss_name_to_cummulative_rel_abund_mp_dict to create
        # it is a list of the ref_seqs_ordered first by clade then by abundance.
        self.clade_abundance_ordered_ref_seq_list = []

    def execute_sequence_count_table_ordered_seqs_worker(self):
        all_processes = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        for n in range(self.parent.num_proc):
            p = Process(target=self._sequence_count_table_ordered_seqs_worker, args=())
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

        self._generate_clade_abundance_ordered_ref_seq_list_from_seq_name_abund_dict()

    def _generate_clade_abundance_ordered_ref_seq_list_from_seq_name_abund_dict(self):
        for i in range(len(self.parent.ordered_list_of_clades_found)):
            temp_within_clade_list_for_sorting = []
            for seq_name, abund_val in self.annotated_dss_name_to_cummulative_rel_abund_mp_dict.items():
                if seq_name.startswith(
                        self.parent.ordered_list_of_clades_found[i]) or seq_name[-2:] == \
                        f'_{self.parent.ordered_list_of_clades_found[i]}':
                    # then this is a seq of the clade in Q and we should add to the temp list
                    temp_within_clade_list_for_sorting.append((seq_name, abund_val))
            # now sort the temp_within_clade_list_for_sorting and add to the cladeAbundanceOrderedRefSeqList
            sorted_within_clade = [
                a[0] for a in sorted(temp_within_clade_list_for_sorting, key=lambda x: x[1], reverse=True)]

            self.clade_abundance_ordered_ref_seq_list.extend(sorted_within_clade)

    def _sequence_count_table_ordered_seqs_worker(self):

        for dss in iter(self.input_dss_mp_queue.get, 'STOP'):
            sys.stdout.write(f'\r{dss.name}: collecting seq abundances')
            sequence_count_table_ordered_seqs_worker_instance = SequenceCountTableCollectAbundanceWorker(
                parent_handler=self, dss=dss)
            sequence_count_table_ordered_seqs_worker_instance.start_seq_abund_collection()

    def _populate_input_dss_mp_queue(self):
        for dss in self.parent.list_of_dss_objects:
            self.input_dss_mp_queue.put(dss)

        for N in range(self.parent.num_proc):
            self.input_dss_mp_queue.put('STOP')

    def _populate_dss_id_to_list_of_dsss_objects(self):
        for dss in self.parent.list_of_dss_objects:
            sys.stdout.write(f'\r{dss.name}')
            self.dss_id_to_list_of_dsss_objects_mp_dict[dss.id] = list(
                DataSetSampleSequence.objects.filter(data_set_sample_from=dss))


class SequenceCountTableCollectAbundanceWorker:
    def __init__(self, parent_handler, dss):
        self.parent = parent_handler
        self.dss = dss
        self.total_abundance_of_sequences_in_sample = sum([int(a) for a in json.loads(self.dss.cladal_seq_totals)])

    def start_seq_abund_collection(self):
        clade_summary_absolute_dict, clade_summary_relative_dict = \
            self._generate_empty_noname_seq_abund_summary_by_clade_dicts()

        smple_seq_count_aboslute_dict, smple_seq_count_relative_dict = self._generate_empty_seq_name_to_abund_dicts()

        dsss_in_sample = self.parent.dss_id_to_list_of_dsss_objects_mp_dict[self.dss.id]

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
        self.parent.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict[self.dss.id] = [smple_seq_count_aboslute_dict,
                                                                                         smple_seq_count_relative_dict]
        self.parent.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict[self.dss.id] = [
            clade_summary_absolute_dict, clade_summary_relative_dict]

    def _populate_abs_and_rel_abundances_for_dsss(self, dsss, name_unit, smple_seq_count_aboslute_dict,
                                                  smple_seq_count_relative_dict):
        rel_abund_of_dsss = dsss.abundance / self.total_abundance_of_sequences_in_sample
        self.parent.annotated_dss_name_to_cummulative_rel_abund_mp_dict[name_unit] += rel_abund_of_dsss
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
            name_unit = dsss.reference_sequence_of.name
        return name_unit

    def _generate_empty_seq_name_to_abund_dicts(self):
        smple_seq_count_aboslute_dict = {seq_name: 0 for seq_name in self.parent.ref_seq_names_clade_annotated}
        smple_seq_count_relative_dict = {seq_name: 0 for seq_name in self.parent.ref_seq_names_clade_annotated}
        return smple_seq_count_aboslute_dict, smple_seq_count_relative_dict

    @staticmethod
    def _generate_empty_noname_seq_abund_summary_by_clade_dicts():
        clade_summary_absolute_dict = {clade: 0 for clade in list('ABCDEFGHI')}
        clade_summary_relative_dict = {clade: 0 for clade in list('ABCDEFGHI')}
        return clade_summary_absolute_dict, clade_summary_relative_dict


class SeqOutputSeriesGeneratorHandler:
    def __init__(self, parent):
        self.parent = parent
        self.output_df_header = self._create_output_df_header()
        self.worker_manager = Manager()
        # dss.id : [pandas_series_for_absolute_abundace, pandas_series_for_absolute_abundace]
        self.dss_id_to_pandas_series_results_list_mp_dict = self.worker_manager.dict()
        self.dss_input_queue = Queue()
        self._populate_dss_input_queue()

    def execute_sequence_count_table_dataframe_contructor_handler(self):
        all_processes = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        sys.stdout.write('\n\nOutputting seq data\n')
        for N in range(self.parent.num_proc):
            p = Process(target=self._output_df_contructor_worker, args=())
            all_processes.append(p)
            p.start()

        for p in all_processes:
            p.join()

    def _output_df_contructor_worker(self):
        for dss in iter(self.dss_input_queue.get, 'STOP'):
            seq_output_series_generator_worker = SeqOutputSeriesGeneratorWorker(parent=self, dss=dss)
            seq_output_series_generator_worker.make_series()

    def _populate_dss_input_queue(self):
        for dss in self.parent.list_of_dss_objects:
            self.dss_input_queue.put(dss)

        for N in range(self.parent.num_proc):
            self.dss_input_queue.put('STOP')

    def _create_output_df_header(self):
        header_pre = self.parent.clade_abundance_ordered_ref_seq_list
        no_name_summary_strings = ['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']]
        qc_stats = [
            'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs', 'post_taxa_id_absolute_symbiodinium_seqs',
            'post_taxa_id_unique_symbiodinium_seqs', 'size_screening_violation_absolute',
            'size_screening_violation_unique',
            'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_taxa_id_unique_non_symbiodinium_seqs',
            'post_med_absolute',
            'post_med_unique']

        # append the noName sequences as individual sequence abundances
        return ['sample_name'] + qc_stats + no_name_summary_strings + header_pre


class SeqOutputSeriesGeneratorWorker:
    def __init__(self, parent, dss):
        self.parent = parent
        self.dss = dss
        # dss.id : [{dsss:absolute abundance in dss}, {dsss:relative abundance in dss}]
        # dss.id : [{clade:total absolute abundance of no name seqs from that clade},
        #           {clade:total relative abundance of no name seqs from that clade}
        #          ]
        self.sample_row_data_absolute = []
        self.sample_row_data_relative = []
        self.sample_seq_tot = sum([int(a) for a in json.loads(dss.cladal_seq_totals)])

    def make_series(self):
        sys.stdout.write(f'\r{self.dss.name}: Creating data ouput row')
        if self._dss_had_problem_in_processing():
            self.sample_row_data_absolute.append(self.dss.name)
            self.sample_row_data_relative.append(self.dss.name)

            self._populate_quality_control_data_of_failed_sample()

            self._output_the_failed_sample_pandas_series()
            return

        self._populate_quality_control_data_of_successful_sample()

        self._output_the_successful_sample_pandas_series()

    def _output_the_successful_sample_pandas_series(self):
        sample_series_absolute = pd.Series(self.sample_row_data_absolute, index=self.parent.output_df_header, name=self.dss.id)
        sample_series_relative = pd.Series(self.sample_row_data_relative, index=self.parent.output_df_header, name=self.dss.id)
        self.parent.dss_id_to_pandas_series_results_list_mp_dict[self.dss.id] = [
            sample_series_absolute, sample_series_relative]

    def _populate_quality_control_data_of_successful_sample(self):
        # Here we add in the post qc and post-taxa id counts
        # For the absolute counts we will report the absolute seq number
        # For the relative counts we will report these as proportions of the sampleSeqTot.
        # I.e. we will have numbers larger than 1 for many of the values and the symbiodinium seqs should be 1
        self.sample_row_data_absolute.append(self.dss.name)
        self.sample_row_data_relative.append(self.dss.name)

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
        tax_id_symbiodinium_absolute = self.dss.absolute_num_sym_seqs
        self.sample_row_data_absolute.append(tax_id_symbiodinium_absolute)
        self.sample_row_data_relative.append(tax_id_symbiodinium_absolute / self.sample_seq_tot)
        # Same as above but the number of unique seqs
        tax_id_symbiodinium_unique = self.dss.unique_num_sym_seqs
        self.sample_row_data_absolute.append(tax_id_symbiodinium_unique)
        self.sample_row_data_relative.append(tax_id_symbiodinium_unique / self.sample_seq_tot)
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
        tax_id_non_symbiodinium_unique = self.dss.non_sym_unique_num_seqs
        self.sample_row_data_absolute.append(tax_id_non_symbiodinium_unique)
        self.sample_row_data_relative.append(tax_id_non_symbiodinium_unique / self.sample_seq_tot)
        # Post MED absolute
        post_med_absolute = self.dss.post_med_absolute
        self.sample_row_data_absolute.append(post_med_absolute)
        self.sample_row_data_relative.append(post_med_absolute / self.sample_seq_tot)
        # Post MED unique
        post_med_unique = self.dss.post_med_unique
        self.sample_row_data_absolute.append(post_med_unique)
        self.sample_row_data_relative.append(post_med_unique / self.sample_seq_tot)

        # now add the clade divided summaries of the clades
        for clade in list('ABCDEFGHI'):
            self.sample_row_data_absolute.append(
                self.parent.parent.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict[
                    self.dss.id][0][clade])
            self.sample_row_data_relative.append(
                self.parent.parent.dss_id_to_list_of_abs_and_rel_abund_clade_summaries_of_noname_seqs_mp_dict[
                    self.dss.id][1][clade])

        # and append these abundances in order of cladeAbundanceOrderedRefSeqList to
        # the sampleRowDataCounts and the sampleRowDataProps
        for seq_name in self.parent.parent.clade_abundance_ordered_ref_seq_list:
            sys.stdout.write('\rOutputting seq data for {}: sequence {}'.format(self.dss.name, seq_name))
            self.sample_row_data_absolute.append(
                self.parent.parent.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict[
                    self.dss.id][0][seq_name])
            self.sample_row_data_relative.append(
                self.parent.parent.dss_id_to_list_of_abs_and_rel_abund_of_contained_dsss_dicts_mp_dict[
                    self.dss.id][1][seq_name])

    def _output_the_failed_sample_pandas_series(self):
        sample_series_absolute = pd.Series(self.sample_row_data_absolute, index=self.parent.output_df_header, name=self.dss.id)
        sample_series_relative = pd.Series(self.sample_row_data_relative, index=self.parent.output_df_header, name=self.dss.id)
        self.parent.dss_id_to_pandas_series_results_list_mp_dict[self.dss.id] = [sample_series_absolute,
                                                                          sample_series_relative]

    def _populate_quality_control_data_of_failed_sample(self):
        # Add in the qc totals if possible
        # For the proportions we will have to add zeros as we cannot do proportions
        # CONTIGS
        # This is the absolute number of sequences after make.contigs

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
            tax_id_symbiodinium_absolute = self.dss.absolute_num_sym_seqs
            self.sample_row_data_absolute.append(tax_id_symbiodinium_absolute)
            self.sample_row_data_relative.append(0)
        else:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)
        # Same as above but the number of unique seqs
        if self.dss.unique_num_sym_seqs:
            tax_id_symbiodinium_unique = self.dss.unique_num_sym_seqs
            self.sample_row_data_absolute.append(tax_id_symbiodinium_unique)
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
            tax_id_non_symbiodinium_unique = self.dss.non_sym_unique_num_seqs
            self.sample_row_data_absolute.append(tax_id_non_symbiodinium_unique)
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

        # no name clade summaries get 0.
        for _ in list('ABCDEFGHI'):
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)

        # All sequences get 0s
        for _ in self.parent.parent.clade_abundance_ordered_ref_seq_list:
            self.sample_row_data_absolute.append(0)
            self.sample_row_data_relative.append(0)

    def _dss_had_problem_in_processing(self):
        return self.dss.error_in_processing or self.sample_seq_tot == 0
