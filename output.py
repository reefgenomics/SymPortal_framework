from dbApp.models import symportal_framework, data_set, reference_sequence, data_set_sample_sequence, analysis_type, analysis_group, data_set_sample, data_analysis, clade_collection, clade_collection_type
from multiprocessing import Queue, Process, Manager
import sys
from django import db
from datetime import datetime
import os
import json
import statistics
import operator
from general import writeListToDestination
from collections import defaultdict
import pandas as pd
from plotting import generate_stacked_bar_data_analysis_type_profiles
import pickle
from collections import Counter
import numpy as np

def formatOutput_ord(analysisobj, datasubstooutput, call_type, numProcessors=1, noFig=False, output_user=None):
    analysisObj = analysisobj
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
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    # List of the paths to the files that have been output
    output_files_list = []

    dataSubmissionsToOutput = [int(a) for a in str(datasubstooutput).split(',')]

    # Get collection of types that are specific for the dataSubmissions we are looking at
    querySetOfDataSubmissions = data_set.objects.filter(id__in=dataSubmissionsToOutput)

    cladeCollectionsFromDSs = clade_collection.objects.filter(
        dataSetSampleFrom__dataSubmissionFrom__in=querySetOfDataSubmissions)
    cladeCollectionTypesFromThisAnalysisObjAndDataSubmission = clade_collection_type.objects.filter(
        cladeCollectionFoundIn__in=cladeCollectionsFromDSs, analysisTypeOf__dataAnalysisFrom=analysisObj).distinct()

    at = set()
    for cct in cladeCollectionTypesFromThisAnalysisObjAndDataSubmission:
        sys.stdout.write('\rCollecting analysis_type from clade collection type {}'.format(cct))
        at.add(cct.analysisTypeOf)
    at = list(at)


    # Need to get a list of Samples that are part of the dataAnalysis
    # listOfDataSetSamples = list(data_set_sample.objects.filter(dataSubmissionFrom__in=data_set.objects.filter(id__in=[int(a) for a in analysisObj.listOfDataSubmissions.split(',')])))
    listOfDataSetSamples = list(
        data_set_sample.objects.filter(dataSubmissionFrom__in=data_set.objects.filter(id__in=dataSubmissionsToOutput)))

    # Now go through the types, firstly by clade and then by the number of cladeCollections they were found in
    # Populate a 2D list of types with a list per clade
    typesCladalList = [[] for clade in cladeList]
    for att in at:
        try:
            if len(att.listOfCladeCollections) > 0:
                typesCladalList[cladeList.index(att.clade)].append(att)
        except:
            pass

    # Items for for creating the new samples sorted output
    across_clade_type_sample_abund_dict = dict()
    across_clade_sorted_type_order = []
    # Go clade by clade
    for i in range(len(typesCladalList)):
        if typesCladalList[i]:
            cladeInQ = cladeList[i]
            ###### CALCULATE REL ABUND AND SD OF DEF INTRAS FOR THIS TYPE ###############
            #     # For each type we want to calculate the average proportions of the defining seqs in that type
            #     # We will also calculate an SD.
            #     # We will do both of these calculations using the footprintSeqAbundances, listOfCladeCollections
            #     # and orderedFoorprintList attributes of each type

            ########## MAKE GROUP COUNTER AND DICT ###########
            # These are now going to be managed items for use in the MP
            # We want to name the groups that types are found in sequencially
            # To get the next number to name a group we will use the groupCount
            # To look up what number has been assigned to groups that have already been printed we will use the groupDict
            # groupCount = 0
            # groupDict = {}
            ###################################################

            # sort types by the number of samples they were found in for this output (not across whole analysis)
            # returns list of tuples, [0] = analysis_type object, [1] number of ccs found in for this output
            sortedListOfTypes = sort_list_of_types_by_clade_collections_in_current_output(typesCladalList[i],
                                                                                          cladeCollectionTypesFromThisAnalysisObjAndDataSubmission)

            # Here we will MP at the type level
            # i.e. each type will be processed on a differnt core. In order for this to work we should have a managed
            # dictionary where the key can be the types ID and the value can be the datalist
            # In order to get the types output in the correct order we should use the sortedListOfTypes to resort the
            # data once the MPing has been done.

            # This dict will hold all of the output rows that the MPing has created.
            worker_manager = Manager()
            typeOutputManagedDict = worker_manager.dict({type: None for type in sortedListOfTypes})


            #NB using shared items was considerably slowing us down so now I will just use copies of items
            # we don't actually need to use managed items as they don't really need to be shared (i.e. input
            # from different processes.
            # a list that is the IDs of each of the samples in the analyis
            sample_IDs_list = [smp.id for smp in listOfDataSetSamples]
            # listOfDataSetSampleIDsManagedList = worker_manager.list(sample_IDs_list)

            # A corresponding dictionary that is the list of the clade collection IDs for each of the samples
            # that are in the ID list above.
            sample_ID_to_cc_IDs_of_clade = {}
            print('\nCreating sample_ID_to_cc_IDs dictionary clade {}'.format(cladeInQ))
            for smp in listOfDataSetSamples:
                sys.stdout.write('\rCollecting clade collections for sample {}'.format(smp))
                try:
                    sample_ID_to_cc_IDs_of_clade[smp.id] = clade_collection.objects.get(dataSetSampleFrom=smp,
                                                                                        clade=cladeInQ).id
                except clade_collection.DoesNotExist:
                    sample_ID_to_cc_IDs_of_clade[smp.id] = None
                except Exception as ex:  # Just incase there is some weird stuff going on
                    print(ex)
            # sample_ID_to_cc_ID_of_clade_shared = worker_manager.dict(sample_ID_to_cc_IDs_of_clade)


            typeInputQueue = Queue()

            for anType in sortedListOfTypes:
                typeInputQueue.put(anType)

            for N in range(numProcessors):
                typeInputQueue.put('STOP')

            allProcesses = []

            # close all connections to the db so that they are automatically recreated for each process
            # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
            db.connections.close_all()

            sys.stdout.write('\nCalculating ITS2 type profile abundances clade {}\n'.format(cladeList[i]))
            for N in range(numProcessors):
                p = Process(target=outputWorkerOne,
                            args=(typeInputQueue, sample_IDs_list, typeOutputManagedDict, sample_ID_to_cc_IDs_of_clade))
                allProcesses.append(p)
                p.start()

            for p in allProcesses:
                p.join()

            # So we have the current list of samples stored in a mangaed list that is the listOfDataSetSamplesManaged
            # the dictionary that we have below is key = analysis_type and the value is two lists
            # the first list is the raw counts and the second is the prop counts
            # In each of these lists the string is a string of the row values for each of the types
            # we can extract the seq abundances from these values.
            # we are currently doing this by clade so we will have to extract this info for each clade
            # I think we can simply add the items for each of the within clade dicts to each other and end up with
            # an across clades dict then we should be able to work with this to get the order we are looking for

            across_clade_type_sample_abund_dict.update(typeOutputManagedDict)
            across_clade_sorted_type_order.extend(sortedListOfTypes)
            typeOutputManagedDict_dict = dict(typeOutputManagedDict)
            #####################

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
    pre_headers = ['ITS2 type profile UID', 'Clade', 'Majority ITS2 sequence', 'Associated species', 'ITS2 type abundance local', 'ITS2 type abundance DB', 'ITS2 type profile']
    sample_headers = [dataSamp.id for dataSamp in listOfDataSetSamples]
    post_headers = ['Sequence accession / SymPortal UID', 'Average defining sequence proportions and [stdev]']
    columns_for_df = pre_headers + sample_headers + post_headers

    # make absolute
    df_absolute = pd.DataFrame(list_for_df_absolute, columns=columns_for_df)
    df_absolute.set_index('ITS2 type profile', drop=False, inplace=True)

    df_relative = pd.DataFrame(list_for_df_relative, columns=columns_for_df)
    df_relative.set_index('ITS2 type profile', drop=False, inplace=True)

    # add a series that gives us the IDs of the samples incase we have samples that have the same names
    # this rather complex comprehension puts an nan into the list for the pre_header headers
    # (i.e not sample headers) and then puts an ID value for the samples
    # this seires works because we rely on the fact that it will just automatically
    # put nan values for all of the headers after the samples
    data_list_for_sample_name_series = [np.nan if i < len(pre_headers) else listOfDataSetSamples[i - len(pre_headers)].name for i in range(len(pre_headers) + len(sample_headers))]
    sample_name_series = pd.Series(
                        name='sample_name',
                        data=data_list_for_sample_name_series,
                        index=list(df_absolute)[:len(data_list_for_sample_name_series)])


    # it turns out that you cannot have duplicate header values (which makes sense). So we're going to have to
    # work with the sample IDs as the header values and put the sample_name in as the secondary series

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
    typeless_samples_list_by_ID = []
    for i in range(len(pre_headers), len(pre_headers) + len(listOfDataSetSamples)):
        sys.stdout.write('\rGetting type abundance information for {}'.format(listOfDataSetSamples[i - len(pre_headers)]))
        sample_series = df_relative.iloc[:,i]
        sample_abundances_series = sample_series[1:].astype('float')
        max_type_label = sample_abundances_series.idxmax()
        rel_abund_of_max_type = sample_abundances_series[max_type_label]
        if not rel_abund_of_max_type > 0:
            # append the ID of the sample to the list
            smpl_id = sample_series.name
            typeless_samples_list_by_ID.append(smpl_id)
            # sample_id_to_sample_name_dict[smpl_id] = sample_series['sample_name']
        else:
            # append a tuple that is (sample_id, rel abundance of the max type)
            smpl_id = sample_series.name
            type_to_sample_abund_dict[max_type_label].append((smpl_id, rel_abund_of_max_type))
            # sample_id_to_sample_name_dict[smpl_id] = sample_series.name

    # here we have the dictionary populated. We can now go type by type according
    # to the sorted_analysis_type_abundance_list and put the samples that had the given type as their most abundant
    # type, into the sorted sample list, addtionaly sorted by how abund the type was in each of the samples
    samples_by_ID_that_have_been_sorted = []
    # we are only concerned with the types that had samples that had them as most abundant
    for an_type_name in [at.name for at in sorted_analysis_type_abundance_list if at.name in type_to_sample_abund_dict.keys()]:
        samples_by_ID_that_have_been_sorted.extend([a[0] for a in sorted(type_to_sample_abund_dict[an_type_name], key=lambda x: x[1], reverse=True)])

    # here we should have a list of samples that have been sorted according to the types they were found to
    # have as their most abundant
    # now we just need to add the samples that didn't have a type in them to be associated to. Negative etc.
    samples_by_ID_that_have_been_sorted.extend(typeless_samples_list_by_ID)

    # now pickle out the samples_that_have_been_sorted list if we are running on the remote system
    outputDir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'outputs/analyses/{}'.format(analysisObj.id)))
    with open('{}/sp_config'.format(os.path.dirname(__file__))) as f:
        config_dict = json.load(f)
    local_or_remote = config_dict['system_type']
    if local_or_remote == 'remote':
        pickle.dump(samples_by_ID_that_have_been_sorted,
                    open("{}/samples_that_have_been_sorted.pickle".format(outputDir), "wb"))


    # rearange the sample columns so that they are in the new order
    new_sample_headers = samples_by_ID_that_have_been_sorted
    new_cols = pre_headers + new_sample_headers + post_headers

    df_absolute = df_absolute[new_cols]
    df_relative = df_relative[new_cols]

    # transpose
    df_absolute = df_absolute.T
    df_relative = df_relative.T


    outputDir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'outputs/analyses/{}'.format(analysisObj.id)))
    os.makedirs(outputDir, exist_ok=True)
    os.chdir(outputDir)

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
        'S. psygmophilum': 'Lajeunesse, T. C., J. E. Parkinson and J. D. Reimer (2012). A genetics-based description of '
                           'Symbiodinium minutum sp. nov. and S. psygmophilum sp. nov. (dinophyceae), two dinoflagellates '
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
        'S. eurythalpos': 'LaJeunesse, T. C., D. C. Wham, D. T. Pettay, J. E. Parkinson, S. Keshavmurthy and C. A. Chen '
                          '(2014). Ecologically differentiated stress-tolerant endosymbionts in the dinoflagellate genus'
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
        meta_info_string_items = ['Output as part of data_analysis ID: {}; '
                                  'Number of data_set objects as part of analysis = {}; '
                                  'submitting_user: {}; '
                                  'time_stamp: {}'
                                      .format(analysisobj.id,
                                              len(querySetOfDataSubmissions),
                                              analysisobj.submittingUser,
                                              analysisobj.timeStamp)]
        temp_series = pd.Series(meta_info_string_items, index=[list(df_relative)[0]], name='meta_info_summary')
        df_absolute = df_absolute.append(temp_series)
        df_relative = df_relative.append(temp_series)

        for data_set_object in querySetOfDataSubmissions:
            data_set_meta_list = ['Data_set ID: {}; Data_set name: {}; '
                                  'submitting_user: {}; '
                                  'time_stamp: {}'
                                      .format(data_set_object.id, data_set_object.name,
                                              data_set_object.submittingUser,
                                              data_set_object.timeStamp)]
            temp_series = pd.Series(data_set_meta_list, index=[list(df_relative)[0]], name='data_set_info')
            df_absolute = df_absolute.append(temp_series)
            df_relative = df_relative.append(temp_series)
    else:
        # call_type=='stand_alone'
        meta_info_string_items = [
            'Stand_alone output by {} on {}; '
            'data_analysis ID: {}; '
            'Number of data_set objects as part of output = {}; Number of data_set objects as part of analysis = {}'
                .format(output_user, str(datetime.now()).replace(' ', '_').replace(':', '-'), analysisobj.id,
                        len(querySetOfDataSubmissions), len(analysisobj.listOfDataSubmissions.split(',')))]
        temp_series = pd.Series(meta_info_string_items, index=[list(df_relative)[0]], name='meta_info_summary')
        df_absolute = df_absolute.append(temp_series)
        df_relative = df_relative.append(temp_series)
        for data_set_object in querySetOfDataSubmissions:
            data_set_meta_list = ['Data_set ID: {}; Data_set name: {}; '
                                  'submitting_user: {}; '
                                  'time_stamp: {}'
                                      .format(data_set_object.id, data_set_object.name,
                                              data_set_object.submittingUser,
                                              data_set_object.timeStamp)]
            temp_series = pd.Series(data_set_meta_list, index=[list(df_relative)[0]], name='data_set_info')
            df_absolute = df_absolute.append(temp_series)
            df_relative = df_relative.append(temp_series)

    date_time_string = str(datetime.now()).replace(' ', '_').replace(':', '-')
    os.makedirs(outputDir, exist_ok=True)

    path_to_profiles_absolute = '{}/{}_{}_{}.profiles.absolute.txt'.format(outputDir, analysisObj.id, analysisObj.name, date_time_string)
    df_absolute.to_csv(path_to_profiles_absolute, sep="\t", header=False)
    output_files_list.append(path_to_profiles_absolute)

    del df_absolute

    path_to_profiles_rel = '{}/{}_{}_{}.profiles.relative.txt'.format(outputDir, analysisObj.id, analysisObj.name, date_time_string)
    df_relative.to_csv(path_to_profiles_rel, sep="\t", header=False)
    # writeListToDestination(path_to_profiles_rel, outputTableTwo)
    output_files_list.append(path_to_profiles_rel)

    del df_relative

    # ########################## ITS2 INTRA ABUND COUNT TABLE ################################
    div_output_pre_analysis_new_meta_and_new_dss_structure(datasubstooutput=datasubstooutput,
                                                           numProcessors=numProcessors, output_dir=outputDir,
                                                           sorted_sample_ID_list=samples_by_ID_that_have_been_sorted,
                                                           analysis_obj_id=analysisobj.id, call_type='analysis')

    print('ITS2 type profile output files:')
    for output_file in output_files_list:
        print(output_file)
        if 'relative' in output_file:
            output_to_plot = output_file

    # Finally lets produce output plots for the dataoutput. For the time being this should just be a
    # plot for the ITS2 type profiles and one for the sequences
    # as with the data_submission let's pass in the path to the outputfiles that we can use to make the plot with
    output_dir = os.path.dirname(output_to_plot)
    if not noFig:
        svg_path, png_path = generate_stacked_bar_data_analysis_type_profiles(path_to_tab_delim_count=output_to_plot,
                                                         output_directory=output_dir,
                                                         data_set_id_str=datasubstooutput,
                                                         analysis_obj_id=analysisobj.id)
    print('Figure output files:')
    print(svg_path)
    print(png_path)

    return output_dir

def sort_list_of_types_by_clade_collections_in_current_output(list_of_analysis_types,
                                                              clade_collection_types_in_current_output):
    # create a list of tupples that is the type ID and the number of this output's CCs that it was found in.
    sys.stdout.write('\nSorting types by abundance of clade collections\n')
    tuple_list = []
    clade_collections_ID_in_current_output = [clade_collection_type_obj.cladeCollectionFoundIn.id for
                                           clade_collection_type_obj in clade_collection_types_in_current_output]
    for at in list_of_analysis_types:
        sys.stdout.write('\rCounting for type {}'.format(at))
        list_of_clade_collections_found_in = [int(x) for x in at.listOfCladeCollections.split(',')]
        num_clade_collections_of_output = list(
            set(clade_collections_ID_in_current_output).intersection(list_of_clade_collections_found_in))
        tuple_list.append((at.id, len(num_clade_collections_of_output)))

    type_IDs_sorted = sorted(tuple_list, key=lambda x: x[1], reverse=True)

    return [analysis_type.objects.get(id=x[0]) for x in type_IDs_sorted]

def outputWorkerOne(input, listOfDataSetSample_IDs, outputDict, sample_ID_to_cc_of_clade_ID):
    num_samples = len(listOfDataSetSample_IDs)
    for anType in iter(input.get, 'STOP'):  # Within each type go through each of the samples

        sys.stdout.write('\rProcessing ITS2 type profile: {}'.format(anType))

        ###### CALCULATE REL ABUND AND SD OF DEF INTRAS FOR THIS TYPE ###############
        # For each type we want to calculate the average proportions of the defining seqs in that type
        # We will also calculate an SD.
        # We will do both of these calculations using the footprintSeqAbundances, listOfCladeCollections
        # and orderedFoorprintList attributes of each type
        footprintAbundances = json.loads(anType.footprintSeqAbundances)
        orderedFootPrintlist = anType.orderedFootprintList.split(',')

        # We have to make a decision as to whether this average should represent all the findings of this type
        # or whether we should only represent the averages of samples found in this dataset.
        # I think it needs to be a global average. Because the type is defined based on all samples in the
        # SymPortal db.

        # Calculate the average proportion of each DIV as a proportion of the absolute abundances of the divs
        # of the type within the samples the type is found in
        div_abundance_df = pd.DataFrame(footprintAbundances)
        # https://stackoverflow.com/questions/26537878/pandas-sum-across-columns-and-divide-each-cell-from-that-value
        # convert each cell to a proportion as a function of the sum of the row
        div_abundance_df_proportion = div_abundance_df.div(div_abundance_df.sum(axis=1),
                                                           axis=0)
        div_abundance_df_proportion_T = div_abundance_df_proportion.T

        TotList = list(div_abundance_df_proportion_T.mean(axis=1))
        SDList = list(div_abundance_df_proportion_T.std(axis=1))

        # The TotList now contains the proportions of each def seq
        # The SDList now contains the SDs for each of the def seqs proportions
        ###########################################################################

        # Counter that will increment with every sample type is found in
        # This is counter will end up being the type abund local value
        abundanceCount = 0

        typeInQ = anType
        clade = typeInQ.clade

        # For each type create a holder that will hold 0s for each sample until populated below
        dataRowRaw = [0 for smpl in range(num_samples)]
        dataRowProp = [0 for smpl in range(num_samples)]
        typeCCIDs = [int(a) for a in typeInQ.listOfCladeCollections.split(',')]

        # We also need the type abund db value. We can get this from the type cladeCollections
        globalCount = len(typeCCIDs)

        # Within each type go through each of the samples
        # Do we really have to go through every sample? I don't think so.
        # Because we have only one cc ID per sample we can simply identify the
        # sample ID (keys) in the sample_ID_to_cc_of_clade_ID dict where the cc ID (value) is found in the typeCCIDs.
        IDs_of_samples_that_had_type = [smp_id for smp_id in listOfDataSetSample_IDs if
                                        sample_ID_to_cc_of_clade_ID[smp_id] in typeCCIDs]
        for ID in IDs_of_samples_that_had_type:
            abundanceCount += 1
            # Need to work out how many seqs were found from the sample for this type
            # Also need to work this out as a proportion of all of the Symbiodinium seqs found in the sample
            try:
                ccInQ_ID = sample_ID_to_cc_of_clade_ID[ID]
            except:
                apples = 'asdf'

            # The number of sequences that were returned for the sample in question
            try:
                totSeqNum = data_set_sample.objects.get(id=ID).finalTotSeqNum
            except:
                apples = 'asdf'

            # The number of sequences that make up the type in q
            # The index of the CC in the type's listOfCladeCollections
            CCindexInType = typeCCIDs.index(ccInQ_ID)
            # List containing the abundances of each of the ref seqs that make up the type in the given CC
            seqAbundanceInfoForCCandTypeInQ = json.loads(typeInQ.footprintSeqAbundances)[CCindexInType]
            # total of the above list, i.e. seqs in type
            sumOfDefRefSeqAbundInType = sum(seqAbundanceInfoForCCandTypeInQ)
            # type abundance as proportion of the total seqs found in the sample
            typeProp = sumOfDefRefSeqAbundInType / totSeqNum
            # Now populate the dataRow with the sumOfDefRefSeqAbundInType and the typeProp
            index_of_sample_ID_in_listOfDataSetSample_IDs = listOfDataSetSample_IDs.index(ID)
            dataRowRaw[index_of_sample_ID_in_listOfDataSetSample_IDs] = sumOfDefRefSeqAbundInType
            dataRowProp[index_of_sample_ID_in_listOfDataSetSample_IDs] = typeProp

        # Type Profile
        typeID = anType.id

        typeProfile = typeInQ.name
        species = typeInQ.species
        typeAbundance = abundanceCount
        # This is currently putting out the Majs in an odd order due to the fact that the MajRefSeqSet is
        # a set. Instead we should get the Majs from the name.
        # majITS2 = '/'.join([str(refseq) for refseq in sortedListOfTypes[j].getMajRefSeqSet()])
        majITS2, majList = getMajList(anType)
        typeAbundAndSDString = getAbundStr(TotList, SDList, majList)

        sequenceAccession = typeInQ.generateName(accession=True)
        rowOne = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(typeID, clade, majITS2,
                                                             species, str(typeAbundance), str(globalCount), typeProfile,
                                                             '\t'.join([str(a) for a in dataRowRaw]),
                                                             sequenceAccession, typeAbundAndSDString)
        rowTwo = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(typeID, clade, majITS2,
                                                             species, str(typeAbundance), str(globalCount), typeProfile,
                                                             '\t'.join(["{:.3f}".format(prop) for prop in
                                                                        dataRowProp]), sequenceAccession,
                                                             typeAbundAndSDString)
        # Here we will use the outputDict instead of adding it to the outputTableOne
        # We will then get these elements in the dict into the right order with the types we will then
        # write out each of the type elements in the order of the orderedTypeLists

        outputDict[anType] = [rowOne, rowTwo]

def getMajList(atype):

    # This is a little tricky. Have to take into account that the Maj seqs are not always more abundant
    # than the non-Maj seqs.
    # e.g. type B5/B5e-B5a-B5b, B5e is acutally the least abundant of the seqs
    # Therefore we cannot simply grab the orderedFooprintList seqs in order assuming they are the Majs

    name = atype.name
    count = name.count('/')
    majList = []
    # list of the seqs in order of abundance across the type's samples
    seqsInOrderOfAbunIDs = atype.orderedFootprintList.split(',')
    # list of the maj seqs in the type
    majSeqsIDs = atype.MajRefSeqSet.split(',')
    for index in range(count + 1):
        for item in range(len(seqsInOrderOfAbunIDs)):
            if seqsInOrderOfAbunIDs[item] in majSeqsIDs:
                maj_seq_obj = reference_sequence.objects.get(id=int(seqsInOrderOfAbunIDs[item]))
                if maj_seq_obj.hasName:
                    majList.append(maj_seq_obj.name)
                else:
                    majList.append(str(maj_seq_obj.id))
                del seqsInOrderOfAbunIDs[item]
                break
    majStringOutput = '/'.join(majList)
    return majStringOutput, majList

def getAbundStr(totlist, sdlist, majlist):
    majComp = []
    lessComp = []
    for i in range(len(totlist)):
        totString = "{0:.3f}".format(totlist[i])
        sdString = "{0:.3f}".format(sdlist[i])
        if i in range(len(majlist)):
            majComp.append('{}[{}]'.format(totString, sdString))
        else:
            lessComp.append('{}[{}]'.format(totString, sdString))
    majCompStr = '/'.join(majComp)
    lessCompStr = '-'.join(lessComp)
    if lessCompStr:
        abundOutputStr = '-'.join([majCompStr, lessCompStr])
    else:
        abundOutputStr = majCompStr
    return abundOutputStr

def div_output_pre_analysis_new_meta_and_new_dss_structure(datasubstooutput, numProcessors, output_dir, call_type,
                                                           sorted_sample_ID_list=None, analysis_obj_id=None, output_user=None ):


    ########################## ITS2 INTRA ABUND COUNT TABLE ################################
    # This is where we're going to have to work with the sequences that aren't part of a type.
    # Esentially we want to end up with the noName sequences divieded up cladally.
    # So at the moment this will mean that we will divide up the current no names but also
    # add information about the other cladal sequences that didn't make it into a cladeCollection

    # list to hold the paths of the outputted files
    output_path_list = []

    ################ GET ORDERED LIST OF INTRAS BY CLADE THEN ABUNDANCE #################
    clade_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    dataSubmissionsToOutput = [int(a) for a in datasubstooutput.split(',')]

    # Get collection of data_sets that are specific for the dataSubmissions we are looking at
    querySetOfDataSubmissions = data_set.objects.filter(id__in=dataSubmissionsToOutput)



    cladeAbundanceOrderedRefSeqList = []
    refSeqsInDSs = reference_sequence.objects.filter(
        data_set_sample_sequence__data_set_sample_from__dataSubmissionFrom__in=querySetOfDataSubmissions).distinct()

    # Get list of clades found
    cladeSet = set()
    for refSeq in refSeqsInDSs:
        cladeSet.add(refSeq.clade)

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
    sub_clade_list = [a for a in clade_list if a in cladeSet]

    sys.stdout.write('\n')

    # The manager that we will build all of the shared dictionaries from below.
    worker_manager = Manager()

    sampleList = data_set_sample.objects.filter(dataSubmissionFrom__in=querySetOfDataSubmissions)

    # Dictionary that will hold the list of data_set_sample_sequences for each sample
    sample_to_dsss_list_shared_dict = worker_manager.dict()
    print('Creating sample to data_set_sample_sequence dict:')
    for dss in sampleList:
        sys.stdout.write('\r{}'.format(dss.name))
        sample_to_dsss_list_shared_dict[dss.id] = list(
            data_set_sample_sequence.objects.filter(data_set_sample_from=dss))

    # Queue that will hold the data set samples for the MP
    dssQueue = Queue()

    # I will have a set of three dictionaries to pass into worker 2
    # 1 - Seqname to cumulative abundance of relative abundances (for getting the over lying order of ref seqs)
    # 2 - sample_id : list(dict(seq:abund), dict(seq:rel_abund))
    # 3 - sample_id : list(dict(noNameClade:abund), dict(noNameClade:rel_abund)

    refSeq_names_annotated = [refSeq.name if refSeq.hasName else str(refSeq.id) + '_{}'.format(refSeq.clade) for refSeq in refSeqsInDSs]
    generic_seq_to_abund_dict = {refSeq_name:0 for refSeq_name in refSeq_names_annotated}

    list_of_dicts_for_processors = []
    for n in range(numProcessors):
        list_of_dicts_for_processors.append((worker_manager.dict(generic_seq_to_abund_dict), worker_manager.dict(), worker_manager.dict()))

    for dss in sampleList:
        dssQueue.put(dss)

    for N in range(numProcessors):
        dssQueue.put('STOP')

    allProcesses = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()


    for n in range(numProcessors):
        p = Process(target=outputWorkerTwo, args=(dssQueue, list_of_dicts_for_processors[n][0],
                                                  list_of_dicts_for_processors[n][1],
                                                  list_of_dicts_for_processors[n][2],
                                                  refSeq_names_annotated, sample_to_dsss_list_shared_dict))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

    print('\nCollecting results of data_set_sample_counting across {} dictionaries'.format(numProcessors))
    master_seq_abundance_counter = Counter()
    master_smple_seq_dict = dict()
    master_smple_noName_clade_summary = dict()

    # now we need to do different actions for each of the three dictionary sets. One set for each numProc
    # for the seqName counter we simply need to add to the master counter as we were doing before
    # for both of the sample-centric dictionaries we simply need to update a master dictionary
    for n in range(len(list_of_dicts_for_processors)):
        sys.stdout.write('\rdictionary {}(0)/{}'.format(n, numProcessors))
        master_seq_abundance_counter += Counter(dict(list_of_dicts_for_processors[n][0]))
        sys.stdout.write('\rdictionary {}(1)/{}'.format(n, numProcessors))
        master_smple_seq_dict.update(dict(list_of_dicts_for_processors[n][1]))
        sys.stdout.write('\rdictionary {}(2)/{}'.format(n, numProcessors))
        master_smple_noName_clade_summary.update(dict(list_of_dicts_for_processors[n][2]))

    print('Collection complete.')


    # we now need to separate by clade and sort within the clade
    cladeAbundanceOrderedRefSeqList = []
    for i in range(len(sub_clade_list)):
        temp_within_clade_list_for_sorting = []
        for seq_name, abund_val in master_seq_abundance_counter.items():
            if seq_name.startswith(sub_clade_list[i]) or seq_name[-2:] == '_{}'.format(sub_clade_list[i]):
                # then this is a seq of the clade in Q and we should add to the temp list
                temp_within_clade_list_for_sorting.append((seq_name, abund_val))
        # now sort the temp_within_clade_list_for_sorting and add to the cladeAbundanceOrderedRefSeqList
        sorted_within_clade = [a[0] for a in sorted(temp_within_clade_list_for_sorting, key=lambda x: x[1], reverse=True)]
        cladeAbundanceOrderedRefSeqList.extend(sorted_within_clade)

    # now delete the master_seq_abundance_counter as we are done with it
    del master_seq_abundance_counter

    ###### WORKER THREE DOMAIN

    # we will eventually have the outputs stored in pandas dataframes.
    # in the worker below we will create a set of pandas.Series for each of the samples which will hold the abundances
    # one for the absoulte abundance and one for the relative abundances.

    # we will put together the headers piece by piece
    headerPre = cladeAbundanceOrderedRefSeqList
    no_name_summary_strings = ['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']]
    qc_stats = ['raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                          'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                          'size_screening_violation_absolute', 'size_screening_violation_unique',
                          'post_taxa_id_absolute_non_symbiodinium_seqs',
                          'post_taxa_id_unique_non_symbiodinium_seqs', 'post_med_absolute', 'post_med_unique']

    # append the noName sequences as individual sequence abundances
    output_header = ['sample_name'] + qc_stats + no_name_summary_strings + headerPre


    ######################################################################################

    ############ POPULATE TABLE WITH CC DATA #############

    # In order to MP this we will have to pay attention to the order. As we can't keep the order as we work with the
    # MPs we will do as we did above for the profies outputs and have an output dict that will have the sample as key
    # and its corresponding data row as the value. Then when we have finished the MPing we will go through the
    # sampleList order and output the data in this order.
    # So we will need a managed output dict.

    worker_manager = Manager()
    managedSampleOutputDict = worker_manager.dict()

    dssQueue = Queue()

    for dss in sampleList:
        dssQueue.put(dss)

    for N in range(numProcessors):
        dssQueue.put('STOP')

    allProcesses = []

    # close all connections to the db so that they are automatically recreated for each process
    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    db.connections.close_all()

    sys.stdout.write('\n\nOutputting DIV data\n')
    for N in range(numProcessors):
        p = Process(target=outputWorkerThree_pre_analysis_new_dss_structure, args=(
            dssQueue, managedSampleOutputDict, cladeAbundanceOrderedRefSeqList, output_header,
            master_smple_seq_dict, master_smple_noName_clade_summary))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

    print('\nDIV output complete\n')

    #### DEVELOPMENT ####
    # So that we can inspect the managedSampleOutputDict
    managedSampleOutputDict_dict = dict(managedSampleOutputDict)

    #####################

    # now we need to populate the output dataframe with the sample series in order of the sorted_sample_list
    # if there is one.
    # If there is a sorted sample list, make sure that it matches the samples that we are outputting

    # TODO we are having an issue with data_sets having the same names. To fix this, we should do our ordering
    # accoring to the IDs of the samples

    if sorted_sample_ID_list:
        sys.stdout.write('\nValidating sorted sample list and ordering dataframe accordingly\n')
        if len(sorted_sample_ID_list) != len(sampleList):
            sys.exit('Number of items in sorted_sample_list do not match those to be outputted!')
        if list(set(sorted_sample_ID_list).difference(set([s.id for s in sampleList]))):
            # then there is a sample that doesn't match up from the sorted_sample_ID_list that
            # has been passed in and the unordered sample list that we are working with in the code
            sys.exit('Sample list passed in does not match sample list from db query')

        # if we got to here then the sorted_sample_list looks good
        sys.stdout.write('\rPopulating the absolute dataframe with series. This could take a second...')
        output_df_absolute = pd.concat([list_of_series[0] for list_of_series in managedSampleOutputDict.values()], axis=1)
        sys.stdout.write('\rPopulating the relative dataframe with series. This could take a second...')
        output_df_relative = pd.concat([list_of_series[1] for list_of_series in managedSampleOutputDict.values()], axis=1)

        # now transpose
        output_df_absolute = output_df_absolute.T
        output_df_relative = output_df_relative.T



        # now make sure that the order is correct.
        output_df_absolute = output_df_absolute.reindex(sorted_sample_ID_list)
        output_df_relative = output_df_relative.reindex(sorted_sample_ID_list)

    else:
        # TODO we should aim to work with IDs here.
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
        ordered_sample_list_by_ID = generate_ordered_sample_list(managedSampleOutputDict)

        # if we got to here then the sorted_sample_list looks good
        sys.stdout.write('\rPopulating the absolute dataframe with series. This could take a second...')
        output_df_absolute = pd.concat([list_of_series[0] for list_of_series in managedSampleOutputDict.values()],
                                       axis=1)
        sys.stdout.write('\rPopulating the relative dataframe with series. This could take a second...')
        output_df_relative = pd.concat([list_of_series[1] for list_of_series in managedSampleOutputDict.values()],
                                       axis=1)

        # now transpose
        sys.stdout.write('\rTransposing...')
        output_df_absolute = output_df_absolute.T
        output_df_relative = output_df_relative.T



        # now make sure that the order is correct.
        sys.stdout.write('\rReordering index...')
        output_df_absolute = output_df_absolute.reindex(ordered_sample_list_by_ID)
        output_df_relative = output_df_relative.reindex(ordered_sample_list_by_ID)

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
    #TODO this is painfully slow because we are doing individual calls to the dictionary
    # I think this will be much faster if do two queries of the db to get the named and
    # non named refseqs and then make two dicts for each of these and use these to populate the below
    refSeqsInDSs_noName = reference_sequence.objects.filter(
        data_set_sample_sequence__data_set_sample_from__dataSubmissionFrom__in=querySetOfDataSubmissions,
        hasName=False).distinct()
    refSeqsInDSs_hasName = reference_sequence.objects.filter(
        data_set_sample_sequence__data_set_sample_from__dataSubmissionFrom__in=querySetOfDataSubmissions,
        hasName=True).distinct()
    # no name dict should be a dict of id to sequence
    no_name_dict = {rs.id: rs.sequence for rs in refSeqsInDSs_noName}
    # has name dict should be a dict of name to sequence
    has_name_dict = {rs.name: (rs.id, rs.sequence) for rs in refSeqsInDSs_hasName}

    # for the time being we are going to ignore whether a refseq has an assession as we have not put this
    # into use yet.
    accession_list = []
    num_cols = len(list(output_df_relative))
    for i, col_name in enumerate(list(output_df_relative)):
        sys.stdout.write('\rAppending accession info and creating fasta {}: {}/{}'.format(col_name, i, num_cols))
        if col_name in cladeAbundanceOrderedRefSeqList:
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
    # this is information like the submitting user, what the IDs of the datasets are etc.
    # There are several ways that this can be called.
    # it can be called as part of the submission: call_type = submission
    # part of an analysis output: call_type = analysis
    # or stand alone: call_type = 'stand_alone'
    # we should have an output for each scenario

    if call_type=='submission':
        data_set_object = querySetOfDataSubmissions[0]
        # there will only be one data_set object
        meta_info_string_items = ['Output as part of data_set submission ID: {}; '
                                  'submitting_user: {}; '
                                  'time_stamp: {}'
                                      .format(data_set_object.id,
                                              data_set_object.submittingUser,
                                              data_set_object.timeStamp)
                                  ]
        temp_series = pd.Series(meta_info_string_items, index=[list(output_df_absolute)[0]], name='meta_info_summary')
        output_df_absolute = output_df_absolute.append(temp_series)
        output_df_relative = output_df_relative.append(temp_series)
    elif call_type=='analysis':
        data_analysis_obj = data_analysis.objects.get(id=analysis_obj_id)
        meta_info_string_items = [
            'Output as part of data_analysis ID: {}; '
            'Number of data_set objects as part of analysis = {}; '
            'submitting_user: {}; time_stamp: {}'.format(
                data_analysis_obj.id, len(data_analysis_obj.listOfDataSubmissions.split(',')), data_analysis_obj.submittingUser, data_analysis_obj.timeStamp)]
        temp_series = pd.Series(meta_info_string_items, index=[list(output_df_absolute)[0]], name='meta_info_summary')
        output_df_absolute = output_df_absolute.append(temp_series)
        output_df_relative = output_df_relative.append(temp_series)
        for data_set_object in querySetOfDataSubmissions:
            data_set_meta_list = ['Data_set ID: {}; Data_set name: {}; '
                                  'submitting_user: {}; '
                                  'time_stamp: {}'
                                      .format(data_set_object.id, data_set_object.name,
                                              data_set_object.submittingUser,
                                              data_set_object.timeStamp)]
            temp_series = pd.Series(data_set_meta_list, index=[list(output_df_absolute)[0]], name='data_set_info')
            output_df_absolute = output_df_absolute.append(temp_series)
            output_df_relative = output_df_relative.append(temp_series)
    else:
        # call_type=='stand_alone'
        meta_info_string_items = [
            'Stand_alone output by {} on {}; '
            'Number of data_set objects as part of output = {}'
                .format(output_user, str(datetime.now()).replace(' ', '_').replace(':', '-'), len(querySetOfDataSubmissions))]
        temp_series = pd.Series(meta_info_string_items, index=[list(output_df_absolute)[0]], name='meta_info_summary')
        output_df_absolute = output_df_absolute.append(temp_series)
        output_df_relative = output_df_relative.append(temp_series)
        for data_set_object in querySetOfDataSubmissions:
            data_set_meta_list = ['Data_set ID: {}; Data_set name: {}; '
                                          'submitting_user: {}; '
                                          'time_stamp: {}'
                                          .format(data_set_object.id, data_set_object.name,
                                                  data_set_object.submittingUser,
                                                  data_set_object.timeStamp)]
            temp_series = pd.Series(data_set_meta_list, index=[list(output_df_absolute)[0]], name='data_set_info')
            output_df_absolute = output_df_absolute.append(temp_series)
            output_df_relative = output_df_relative.append(temp_series)


    # Here we have the tables populated and ready to output
    date_time_string = str(datetime.now()).replace(' ', '_').replace(':', '-')
    if analysis_obj_id:
        data_analysis_obj = data_analysis.objects.get(id=analysis_obj_id)
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
    writeListToDestination(fasta_path, fasta_output_list)
    output_path_list.append(fasta_path)

    print('\nITS2 sequence output files:')
    for path_item in output_path_list:
        print(path_item)

    return output_path_list

def outputWorkerTwo(input, seq_rel_abund_dict, smpl_seq_dict, smpl_noName_clade_summary_dict,
                    refSeq_names_annotated, sample_to_dsss_list_shared_dict):
    # 1 - Seqname to cumulative abundance of relative abundances (for getting the over lying order of ref seqs)
    # 2 - sample_id : list(dict(seq:abund), dict(seq:rel_abund))
    # 3 - sample_id : list(dict(noNameClade:abund), dict(noNameClade:rel_abund)
    for dss in iter(input.get, 'STOP'):
        sys.stdout.write('\rCounting seqs for {}'.format(dss))

        cladalAbundances = [int(a) for a in json.loads(dss.cladalSeqTotals)]

        sampleSeqTot = sum(cladalAbundances)

        # if dss.errorInProcessing or sampleSeqTot == 0:
        #     continue


        # the first dict will hold the absolute abundances, whilst the second will hold the relative abundances
        clade_summary_absolute_dict = {clade:0 for clade in list('ABCDEFGHI')}
        clade_summary_relative_dict = {clade:0 for clade in list('ABCDEFGHI')}
        smple_seq_count_aboslute_dict = {seq_name:0 for seq_name in refSeq_names_annotated}
        smple_seq_count_relative_dict = {seq_name: 0 for seq_name in refSeq_names_annotated}

        dsssInSample = sample_to_dsss_list_shared_dict[dss.id]

        for dsss in dsssInSample:
            # determine what the name of the seq will be in the output
            if not dsss.referenceSequenceOf.hasName:
                name_unit = str(dsss.referenceSequenceOf.id) + '_{}'.format(dsss.referenceSequenceOf.clade)
                # the clade summries are only for the noName seqs
                clade_summary_absolute_dict[dsss.referenceSequenceOf.clade] += dsss.abundance
                clade_summary_relative_dict[dsss.referenceSequenceOf.clade] += dsss.abundance / sampleSeqTot
            else:
                name_unit = dsss.referenceSequenceOf.name

            seq_rel_abund_dict[name_unit] += dsss.abundance / sampleSeqTot
            smple_seq_count_aboslute_dict[name_unit] += dsss.abundance
            smple_seq_count_relative_dict[name_unit] += dsss.abundance / sampleSeqTot

        smpl_noName_clade_summary_dict[dss.id] = [clade_summary_absolute_dict, clade_summary_relative_dict]
        smpl_seq_dict[dss.id] = [smple_seq_count_aboslute_dict, smple_seq_count_relative_dict]



def generate_ordered_sample_list(managedSampleOutputDict):
    # create a df from the managedSampleOutputDict. We will use the relative values here

    output_df_relative = pd.concat([list_of_series[1] for list_of_series in managedSampleOutputDict.values()],
                                   axis=1)
    output_df_relative = output_df_relative.T

    # now remove the rest of the non abundance columns
    non_seq_columns = ['sample_name', 'raw_contigs', 'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_qc_absolute_seqs',
                       'post_qc_unique_seqs',
                       'post_taxa_id_unique_non_symbiodinium_seqs', 'post_taxa_id_absolute_symbiodinium_seqs',
                       'post_taxa_id_unique_symbiodinium_seqs', 'post_med_absolute', 'post_med_unique',
                       'size_screening_violation_absolute', 'size_screening_violation_unique']
    noName_seq_columns = ['noName Clade {}'.format(clade) for clade in list('ABCDEFGHI')]
    cols_to_drop = non_seq_columns + noName_seq_columns

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
    ordered_sample_list_by_ID = []
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
            ordered_sample_list_by_ID.extend(ordered_list_of_samples_for_seq_ordered)
    # finally add in the samples that didn't have a maj sequence
    ordered_sample_list_by_ID.extend(no_maj_samps)
    return ordered_sample_list_by_ID

def outputWorkerThree_pre_analysis_new_dss_structure(input, outDict, cladeAbundanceOrderedRefSeqList, output_header,
                                                     smpl_abund_dicts_dict,
                                                     smpl_clade_summary_dicts_dict):
    cladeList = list('ABCDEFGHI')
    for dss in iter(input.get, 'STOP'):

        sys.stdout.write('\rOutputting DIV data for {}'.format(dss.name))
        # List that will hold the row
        sampleRowDataCounts = []
        sampleRowDataProps = []
        cladalAbundances = [int(a) for a in json.loads(dss.cladalSeqTotals)]
        sampleSeqTot = sum(cladalAbundances)

        if dss.errorInProcessing or sampleSeqTot == 0:
            #Then this sample had a problem in the sequencing and we need to just output 0s across the board
            # QC

            # Append the name of the dss for when we have samples of the same name
            sampleRowDataCounts.append(dss.name)
            sampleRowDataProps.append(dss.name)

            populate_QC_data_of_failed_sample(dss, sampleRowDataCounts, sampleRowDataProps)

            # no name clade summaries get 0.
            for cl in cladeList:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)

            # All sequences get 0s
            for name in cladeAbundanceOrderedRefSeqList:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)

            # Here we need to add the string to the outputDict rather than the intraAbund table objects
            sample_series_absolute = pd.Series(sampleRowDataCounts, index=output_header, name=dss.id)
            sample_series_relative = pd.Series(sampleRowDataCounts, index=output_header, name=dss.id)

            outDict[dss.id] = [sample_series_absolute, sample_series_relative]
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
        sampleRowDataCounts.append(dss.name)
        sampleRowDataProps.append(dss.name)

        populate_QC_data_of_successful_sample(dss, sampleRowDataCounts, sampleRowDataProps, sampleSeqTot)

        # now add the clade divided summaries of the clades
        for clade in cladeList:
            sys.stdout.write('\rOutputting DIV data for {}: clade summary {}'.format(dss.name, clade))
            sampleRowDataCounts.append(smpl_clade_summary_absolute_dict[clade])
            sampleRowDataProps.append(smpl_clade_summary_relative_dict[clade])


        # and append these abundances in order of cladeAbundanceOrderedRefSeqList to
        # the sampleRowDataCounts and the sampleRowDataProps
        for seq_name in cladeAbundanceOrderedRefSeqList:
            sys.stdout.write('\rOutputting DIV data for {}: sequence {}'.format(dss.name, seq_name))
            sampleRowDataCounts.append(smpl_seq_abund_absolute_dict[seq_name])
            sampleRowDataProps.append(smpl_seq_abund_relative_dict[seq_name])


        # Here we need to add the string to the outputDict rather than the intraAbund table objects
        sample_series_absolute = pd.Series(sampleRowDataCounts, index=output_header, name=dss.id)
        sample_series_relative = pd.Series(sampleRowDataProps, index=output_header, name=dss.id)

        outDict[dss.id] = [sample_series_absolute, sample_series_relative]



def populate_QC_data_of_successful_sample(dss, sampleRowDataCounts, sampleRowDataProps, sampleSeqTot):
    # CONTIGS
    # This is the absolute number of sequences after make.contigs
    contig_num = dss.initialTotSeqNum
    sampleRowDataCounts.append(contig_num)
    sampleRowDataProps.append(contig_num / sampleSeqTot)
    # POST-QC
    # store the aboslute number of sequences after sequencing QC at this stage
    post_qc_absolute = dss.post_seq_qc_absolute_num_seqs
    sampleRowDataCounts.append(post_qc_absolute)
    sampleRowDataProps.append(post_qc_absolute / sampleSeqTot)
    # This is the unique number of sequences after the sequencing QC
    post_qc_unique = dss.initialUniqueSeqNum
    sampleRowDataCounts.append(post_qc_unique)
    sampleRowDataProps.append(post_qc_unique / sampleSeqTot)
    # POST TAXA-ID
    # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
    tax_id_symbiodinium_absolute = dss.finalTotSeqNum
    sampleRowDataCounts.append(tax_id_symbiodinium_absolute)
    sampleRowDataProps.append(tax_id_symbiodinium_absolute / sampleSeqTot)
    # Same as above but the number of unique seqs
    tax_id_symbiodinium_unique = dss.finalUniqueSeqNum
    sampleRowDataCounts.append(tax_id_symbiodinium_unique)
    sampleRowDataProps.append(tax_id_symbiodinium_unique / sampleSeqTot)
    # store the absolute number of sequences lost to size cutoff violations
    size_violation_aboslute = dss.size_violation_absolute
    sampleRowDataCounts.append(size_violation_aboslute)
    sampleRowDataProps.append(size_violation_aboslute / sampleSeqTot)
    # store the unique size cutoff violations
    size_violation_unique = dss.size_violation_unique
    sampleRowDataCounts.append(size_violation_unique)
    sampleRowDataProps.append(size_violation_unique / sampleSeqTot)
    # store the abosolute number of sequenes that were not considered Symbiodinium
    tax_id_non_symbiodinum_abosulte = dss.non_sym_absolute_num_seqs
    sampleRowDataCounts.append(tax_id_non_symbiodinum_abosulte)
    sampleRowDataProps.append(tax_id_non_symbiodinum_abosulte / sampleSeqTot)
    # This is the number of unique sequences that were not considered Symbiodinium
    tax_id_non_symbiodinium_unique = dss.nonSymSeqsNum
    sampleRowDataCounts.append(tax_id_non_symbiodinium_unique)
    sampleRowDataProps.append(tax_id_non_symbiodinium_unique / sampleSeqTot)
    # Post MED absolute
    post_med_absolute = dss.post_med_absolute
    sampleRowDataCounts.append(post_med_absolute)
    sampleRowDataProps.append(post_med_absolute / sampleSeqTot)
    # Post MED unique
    post_med_unique = dss.post_med_unique
    sampleRowDataCounts.append(post_med_unique)
    sampleRowDataProps.append(post_med_unique / sampleSeqTot)


def populate_QC_data_of_failed_sample(dss, sampleRowDataCounts, sampleRowDataProps):
    # Add in the qc totals if possible
    # For the proportions we will have to add zeros as we cannot do proportions
    # CONTIGS
    # This is the absolute number of sequences after make.contigs

    if dss.initialTotSeqNum:
        contig_num = dss.initialTotSeqNum
        sampleRowDataCounts.append(contig_num)
        sampleRowDataProps.append(0)
    else:
        sampleRowDataCounts.append(0)
        sampleRowDataProps.append(0)
    # POST-QC
    # store the aboslute number of sequences after sequencing QC at this stage
    if dss.post_seq_qc_absolute_num_seqs:
        post_qc_absolute = dss.post_seq_qc_absolute_num_seqs
        sampleRowDataCounts.append(post_qc_absolute)
        sampleRowDataProps.append(0)
    else:
        sampleRowDataCounts.append(0)
        sampleRowDataProps.append(0)
    # This is the unique number of sequences after the sequencing QC
    if dss.initialUniqueSeqNum:
        post_qc_unique = dss.initialUniqueSeqNum
        sampleRowDataCounts.append(post_qc_unique)
        sampleRowDataProps.append(0)
    else:
        sampleRowDataCounts.append(0)
        sampleRowDataProps.append(0)
    # POST TAXA-ID
    # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
    if dss.finalTotSeqNum:
        tax_id_symbiodinium_absolute = dss.finalTotSeqNum
        sampleRowDataCounts.append(tax_id_symbiodinium_absolute)
        sampleRowDataProps.append(0)
    else:
        sampleRowDataCounts.append(0)
        sampleRowDataProps.append(0)
    # Same as above but the number of unique seqs
    if dss.finalUniqueSeqNum:
        tax_id_symbiodinium_unique = dss.finalUniqueSeqNum
        sampleRowDataCounts.append(tax_id_symbiodinium_unique)
        sampleRowDataProps.append(0)
    else:
        sampleRowDataCounts.append(0)
        sampleRowDataProps.append(0)
    # size violation absolute
    if dss.size_violation_absolute:
        size_viol_ab = dss.size_violation_absolute
        sampleRowDataCounts.append(size_viol_ab)
        sampleRowDataProps.append(0)
    else:
        sampleRowDataCounts.append(0)
        sampleRowDataProps.append(0)
    # size violation unique
    if dss.size_violation_unique:
        size_viol_uni = dss.size_violation_unique
        sampleRowDataCounts.append(size_viol_uni)
        sampleRowDataProps.append(0)
    else:
        sampleRowDataCounts.append(0)
        sampleRowDataProps.append(0)
    # store the abosolute number of sequenes that were not considered Symbiodinium
    if dss.non_sym_absolute_num_seqs:
        tax_id_non_symbiodinum_abosulte = dss.non_sym_absolute_num_seqs
        sampleRowDataCounts.append(tax_id_non_symbiodinum_abosulte)
        sampleRowDataProps.append(0)
    else:
        sampleRowDataCounts.append(0)
        sampleRowDataProps.append(0)
    # This is the number of unique sequences that were not considered Symbiodinium
    if dss.nonSymSeqsNum:
        tax_id_non_symbiodinium_unique = dss.nonSymSeqsNum
        sampleRowDataCounts.append(tax_id_non_symbiodinium_unique)
        sampleRowDataProps.append(0)
    else:
        sampleRowDataCounts.append(0)
        sampleRowDataProps.append(0)
    # post-med absolute
    if dss.post_med_absolute:
        post_med_abs = dss.post_med_absolute
        sampleRowDataCounts.append(post_med_abs)
        sampleRowDataProps.append(0)
    else:
        sampleRowDataCounts.append(0)
        sampleRowDataProps.append(0)
    # post-med absolute
    if dss.post_med_unique:
        post_med_uni = dss.post_med_unique
        sampleRowDataCounts.append(post_med_uni)
        sampleRowDataProps.append(0)
    else:
        sampleRowDataCounts.append(0)
        sampleRowDataProps.append(0)