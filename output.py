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
from plotting import generate_stacked_bar_data_analysis_type_profiles, generate_stacked_bar_data_submission
import pickle
from collections import Counter
import numpy as np


def output_type_count_tables(
        analysisobj, datasubstooutput, call_type, num_samples,
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
        dataSetSampleFrom__dataSubmissionFrom__in=query_set_of_data_sets)
    clade_collection_types_from_this_data_analysis_and_data_set = CladeCollectionType.objects.filter(
        cladeCollectionFoundIn__in=clade_collections_from_data_sets,
        analysisTypeOf__dataAnalysisFrom=analysis_object).distinct()

    at = set()
    for cct in clade_collection_types_from_this_data_analysis_and_data_set:
        sys.stdout.write('\rCollecting analysis_type from clade collection type {}'.format(cct))
        at.add(cct.analysisTypeOf)
    at = list(at)

    # Need to get a list of Samples that are part of the dataAnalysis
    list_of_data_set_samples = list(
        DataSetSample.objects.filter(
            dataSubmissionFrom__in=DataSet.objects.filter(id__in=data_submissions_to_output)))

    # Now go through the types, firstly by clade and then by the number of cladeCollections they were found in
    # Populate a 2D list of types with a list per clade
    types_cladal_list = [[] for _ in clade_list]
    for att in at:
        try:
            if len(att.listOfCladeCollections) > 0:
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
            #     # We will do both of these calculations using the footprintSeqAbundances, listOfCladeCollections
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
                        dataSetSampleFrom=smp, clade=clade_in_question).id
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
    with open('{}/sp_config'.format(os.path.dirname(__file__))) as f:
        config_dict = json.load(f)
    local_or_remote = config_dict['system_type']
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

    if call_type == 'analysis':
        meta_info_string_items = [
            'Output as part of data_analysis ID: {}; '
            'Number of data_set objects as part of analysis = {}; '
            'submitting_user: {}; time_stamp: {}'.format(
                analysisobj.id, len(query_set_of_data_sets), analysisobj.submittingUser, analysisobj.timeStamp)]
        temp_series = pd.Series(meta_info_string_items, index=[list(df_relative)[0]], name='meta_info_summary')
        df_absolute = df_absolute.append(temp_series)
        df_relative = df_relative.append(temp_series)

        for data_set_object in query_set_of_data_sets:
            data_set_meta_list = [
                'Data_set ID: {}; Data_set name: {}; submitting_user: {}; time_stamp: {}'.format(
                    data_set_object.id, data_set_object.name,
                    data_set_object.submittingUser, data_set_object.timeStamp)]

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
                len(query_set_of_data_sets), len(analysisobj.listOfDataSubmissions.split(',')))]

        temp_series = pd.Series(meta_info_string_items, index=[list(df_relative)[0]], name='meta_info_summary')
        df_absolute = df_absolute.append(temp_series)
        df_relative = df_relative.append(temp_series)
        for data_set_object in query_set_of_data_sets:
            data_set_meta_list = [
                'Data_set ID: {}; Data_set name: {}; submitting_user: {}; time_stamp: {}'.format(
                    data_set_object.id, data_set_object.name,
                    data_set_object.submittingUser, data_set_object.timeStamp)]

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
    div_output_file_list, date_time_string, number_of_samples = div_output_pre_analysis_new_meta_and_new_dss_structure(
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
        if num_samples > 1000:
            print('Too many samples ({}) to generate plots'.format(num_samples))
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

            svg_path, png_path = generate_stacked_bar_data_submission(
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
        clade_collection_type_obj.cladeCollectionFoundIn.id
        for clade_collection_type_obj in clade_collection_types_in_current_output]

    for at in list_of_analysis_types:
        sys.stdout.write('\rCounting for type {}'.format(at))
        list_of_clade_collections_found_in = [int(x) for x in at.listOfCladeCollections.split(',')]
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
        # We will do both of these calculations using the footprintSeqAbundances, listOfCladeCollections
        # and orderedFoorprintList attributes of each type
        footprint_abundances = json.loads(an_type.footprintSeqAbundances)

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
        type_clade_collection_uids = [int(a) for a in type_in_question.listOfCladeCollections.split(',')]

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
            total_number_of_sequences = DataSetSample.objects.get(id=ID).finalTotSeqNum

            # The number of sequences that make up the type in q
            # The index of the clade_collection_object in the type's listOfCladeCollections
            clade_collection_index_in_type = type_clade_collection_uids.index(clade_collection_in_question_uid)
            # List containing the abundances of each of the ref seqs that
            # make up the type in the given clade_collection_object
            sequence_abundance_info_for_clade_collection_and_type_in_question = json.loads(
                type_in_question.footprintSeqAbundances)[clade_collection_index_in_type]
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
        # This is currently putting out the Majs in an odd order due to the fact that the MajRefSeqSet is
        # a set. Instead we should get the Majs from the name.
        # majority_its2 = '/'.join([str(refseq) for refseq in sortedListOfTypes[j].getMajRefSeqSet()])
        majority_its2, majority_list = get_maj_list(an_type)
        type_abundance_and_standard_deviation_string = get_abundance_string(
            total_list, standard_deviation_list, majority_list)

        sequence_accession = type_in_question.generateName(accession=True)
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
    uids_of_sequences_in_order_of_abundance = atype.orderedFootprintList.split(',')
    # list of the maj seqs in the type
    majority_sequeces_uids = atype.MajRefSeqSet.split(',')
    for index in range(count + 1):
        for item in range(len(uids_of_sequences_in_order_of_abundance)):
            if uids_of_sequences_in_order_of_abundance[item] in majority_sequeces_uids:
                maj_seq_obj = ReferenceSequence.objects.get(id=int(uids_of_sequences_in_order_of_abundance[item]))
                if maj_seq_obj.hasName:
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


def div_output_pre_analysis_new_meta_and_new_dss_structure(
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
        data_set_sample_sequence__data_set_sample_from__dataSubmissionFrom__in=query_set_of_data_sets).distinct()

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

    sample_list = DataSetSample.objects.filter(dataSubmissionFrom__in=query_set_of_data_sets)

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
        ref_seq.name if ref_seq.hasName
        else str(ref_seq.id) + '_{}'.format(ref_seq.clade) for ref_seq in reference_sequences_in_data_sets]

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

    sys.stdout.write('\n\nOutputting DIV data\n')
    for N in range(num_processors):
        p = Process(target=output_worker_three, args=(
            data_set_sample_queue, managed_sample_output_dict, clade_abundance_ordered_ref_seq_list, output_header,
            master_smple_seq_dict, master_smple_no_name_clade_summary))
        all_processes.append(p)
        p.start()

    for p in all_processes:
        p.join()

    print('\nDIV output complete\n')

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
        data_set_sample_sequence__data_set_sample_from__dataSubmissionFrom__in=query_set_of_data_sets,
        hasName=False).distinct()
    reference_sequences_in_data_sets_has_name = ReferenceSequence.objects.filter(
        data_set_sample_sequence__data_set_sample_from__dataSubmissionFrom__in=query_set_of_data_sets,
        hasName=True).distinct()
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

    temp_series = pd.Series(accession_list, name='DIV_accession', index=list(output_df_relative))
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
                data_set_object.id, data_set_object.submittingUser, data_set_object.timeStamp)]

        temp_series = pd.Series(meta_info_string_items, index=[list(output_df_absolute)[0]], name='meta_info_summary')
        output_df_absolute = output_df_absolute.append(temp_series)
        output_df_relative = output_df_relative.append(temp_series)
    elif call_type == 'analysis':
        data_analysis_obj = DataAnalysis.objects.get(id=analysis_obj_id)
        meta_info_string_items = [
            'Output as part of data_analysis ID: {}; Number of data_set objects as part of analysis = {}; '
            'submitting_user: {}; time_stamp: {}'.format(
                data_analysis_obj.id, len(data_analysis_obj.listOfDataSubmissions.split(',')),
                data_analysis_obj.submittingUser, data_analysis_obj.timeStamp)]

        temp_series = pd.Series(meta_info_string_items, index=[list(output_df_absolute)[0]], name='meta_info_summary')
        output_df_absolute = output_df_absolute.append(temp_series)
        output_df_relative = output_df_relative.append(temp_series)
        for data_set_object in query_set_of_data_sets:
            data_set_meta_list = [
                'Data_set ID: {}; Data_set name: {}; submitting_user: {}; time_stamp: {}'.format(
                    data_set_object.id, data_set_object.name,
                    data_set_object.submittingUser, data_set_object.timeStamp)]

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
                    data_set_object.id, data_set_object.name, data_set_object.submittingUser,
                    data_set_object.timeStamp)]

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
        path_to_div_absolute = '{}/{}_{}_{}.DIVs.absolute.txt'.format(output_dir, analysis_obj_id,
                                                                      data_analysis_obj.name, date_time_string)
        path_to_div_relative = '{}/{}_{}_{}.DIVs.relative.txt'.format(output_dir, analysis_obj_id,
                                                                      data_analysis_obj.name, date_time_string)
        fasta_path = '{}/{}_{}_{}.DIVs.fasta'.format(output_dir, analysis_obj_id,
                                                     data_analysis_obj.name, date_time_string)

    else:
        path_to_div_absolute = '{}/{}.DIVs.absolute.txt'.format(output_dir, date_time_string)
        path_to_div_relative = '{}/{}.DIVs.relative.txt'.format(output_dir, date_time_string)
        fasta_path = '{}/{}.DIVs.fasta'.format(output_dir, date_time_string)

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
    # 1 - Seqname to cumulative abundance of relative abundances (for getting the over lying order of ref seqs)
    # 2 - sample_id : list(dict(seq:abund), dict(seq:rel_abund))
    # 3 - sample_id : list(dict(noNameClade:abund), dict(noNameClade:rel_abund)
    for dss in iter(input_queue.get, 'STOP'):
        sys.stdout.write('\rCounting seqs for {}'.format(dss))

        cladal_abundances = [int(a) for a in json.loads(dss.cladalSeqTotals)]

        sample_seq_tot = sum(cladal_abundances)

        # the first dict will hold the absolute abundances, whilst the second will hold the relative abundances
        clade_summary_absolute_dict = {clade: 0 for clade in list('ABCDEFGHI')}
        clade_summary_relative_dict = {clade: 0 for clade in list('ABCDEFGHI')}
        smple_seq_count_aboslute_dict = {seq_name: 0 for seq_name in reference_sequence_names_annotated}
        smple_seq_count_relative_dict = {seq_name: 0 for seq_name in reference_sequence_names_annotated}

        dsss_in_sample = sample_to_dsss_list_shared_dict[dss.id]

        for dsss in dsss_in_sample:
            # determine what the name of the seq will be in the output
            if not dsss.referenceSequenceOf.hasName:
                name_unit = str(dsss.referenceSequenceOf.id) + '_{}'.format(dsss.referenceSequenceOf.clade)
                # the clade summries are only for the noName seqs
                clade_summary_absolute_dict[dsss.referenceSequenceOf.clade] += dsss.abundance
                clade_summary_relative_dict[dsss.referenceSequenceOf.clade] += dsss.abundance / sample_seq_tot
            else:
                name_unit = dsss.referenceSequenceOf.name

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
        max_abund_seq = series_as_float.idxmax()
        max_rel_abund = series_as_float.max()
        if not max_rel_abund > 0:
            no_maj_samps.append(sample_to_sort_ID)
        else:
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

        sys.stdout.write('\rOutputting DIV data for {}'.format(dss.name))
        # List that will hold the row
        sample_row_data_counts = []
        sample_row_data_props = []
        cladal_abundances = [int(a) for a in json.loads(dss.cladalSeqTotals)]
        sample_seq_tot = sum(cladal_abundances)

        if dss.errorInProcessing or sample_seq_tot == 0:
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
            sys.stdout.write('\rOutputting DIV data for {}: clade summary {}'.format(dss.name, clade))
            sample_row_data_counts.append(smpl_clade_summary_absolute_dict[clade])
            sample_row_data_props.append(smpl_clade_summary_relative_dict[clade])

        # and append these abundances in order of cladeAbundanceOrderedRefSeqList to
        # the sampleRowDataCounts and the sampleRowDataProps
        for seq_name in clade_abundance_ordered_ref_seq_list:
            sys.stdout.write('\rOutputting DIV data for {}: sequence {}'.format(dss.name, seq_name))
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
    contig_num = dss.initialTotSeqNum
    sample_row_data_counts.append(contig_num)
    sample_row_data_props.append(contig_num / sample_seq_tot)
    # POST-QC
    # store the aboslute number of sequences after sequencing QC at this stage
    post_qc_absolute = dss.post_seq_qc_absolute_num_seqs
    sample_row_data_counts.append(post_qc_absolute)
    sample_row_data_props.append(post_qc_absolute / sample_seq_tot)
    # This is the unique number of sequences after the sequencing QC
    post_qc_unique = dss.initialUniqueSeqNum
    sample_row_data_counts.append(post_qc_unique)
    sample_row_data_props.append(post_qc_unique / sample_seq_tot)
    # POST TAXA-ID
    # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
    tax_id_symbiodinium_absolute = dss.finalTotSeqNum
    sample_row_data_counts.append(tax_id_symbiodinium_absolute)
    sample_row_data_props.append(tax_id_symbiodinium_absolute / sample_seq_tot)
    # Same as above but the number of unique seqs
    tax_id_symbiodinium_unique = dss.finalUniqueSeqNum
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
    tax_id_non_symbiodinium_unique = dss.nonSymSeqsNum
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

    if dss.initialTotSeqNum:
        contig_num = dss.initialTotSeqNum
        sample_row_data_counts.append(contig_num)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # POST-QC
    # store the aboslute number of sequences after sequencing QC at this stage
    if dss.post_seq_qc_absolute_num_seqs:
        post_qc_absolute = dss.post_seq_qc_absolute_num_seqs
        sample_row_data_counts.append(post_qc_absolute)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # This is the unique number of sequences after the sequencing QC
    if dss.initialUniqueSeqNum:
        post_qc_unique = dss.initialUniqueSeqNum
        sample_row_data_counts.append(post_qc_unique)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # POST TAXA-ID
    # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
    if dss.finalTotSeqNum:
        tax_id_symbiodinium_absolute = dss.finalTotSeqNum
        sample_row_data_counts.append(tax_id_symbiodinium_absolute)
        sample_row_data_props.append(0)
    else:
        sample_row_data_counts.append(0)
        sample_row_data_props.append(0)
    # Same as above but the number of unique seqs
    if dss.finalUniqueSeqNum:
        tax_id_symbiodinium_unique = dss.finalUniqueSeqNum
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
    if dss.nonSymSeqsNum:
        tax_id_non_symbiodinium_unique = dss.nonSymSeqsNum
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
