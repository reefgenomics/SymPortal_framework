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

    # at = analysis_type.objects.filter(dataAnalysisFrom=analysisObj)
    IDs = [att.id for att in at]
    outputTableOne = []
    outputTableTwo = []
    # Need to get a list of Samples that are part of the dataAnalysis
    # listOfDataSetSamples = list(data_set_sample.objects.filter(dataSubmissionFrom__in=data_set.objects.filter(id__in=[int(a) for a in analysisObj.listOfDataSubmissions.split(',')])))
    listOfDataSetSamples = list(
        data_set_sample.objects.filter(dataSubmissionFrom__in=data_set.objects.filter(id__in=dataSubmissionsToOutput)))
    # The header row for the table
    # outputTableOne.append('ITS2 type profile UID\tClade\tSpeciesGrp\tType_Profile\tMajority ITS2 Seq\tAssociated Species\tType abund local\tType abund DB\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))
    # outputTableTwo.append('ITS2 type profile UID\tClade\tSpeciesGrp\tType_Profile\tMajority ITS2 Seq\tAssociated Species\tType abund local\tType abund DB\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))

    # I am making some changes to the output format for the SP manuscript
    # outputTableOne.append('ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))
    # outputTableTwo.append('ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format('\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples])))

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

    # now we can sort this list according to the abundance and this will give us the order of types that we want
    sorted_analysis_type_abundance_list = [a[0] for a in sorted(type_to_abund_list, key=lambda x: x[1], reverse=True)]

    # TODO I think we should try to improve this as it is currently horrifically slow.
    # currently for each sample we collect the abundances of the types for each type. Not sure why
    # then for each analysis type we go through these samples and pull out the samples
    # that had it as it maj type. We then sort these samples based on the abundance of that type. crazy.
    # I think we can pandas this.
    # we already have the

    # need to work with proportions here so that we can compare type abundances betweeen samples
    list_for_df_absolute = []
    list_for_df_relative = []

    # get list for the aboslute df
    for an_type_key in across_clade_type_sample_abund_dict.keys():
        list_for_df_absolute.append(across_clade_type_sample_abund_dict[an_type_key][0].split('\t'))

    # get list for the relative df
    for an_type_key in across_clade_type_sample_abund_dict.keys():
        list_for_df_relative.append(across_clade_type_sample_abund_dict[an_type_key][1].split('\t'))

    # headers can be same for each
    columns_for_df =('ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format(
            '\t'.join([dataSamp.name for dataSamp in listOfDataSetSamples]))).split('\t')

    # make absolute
    df_absolute = pd.DataFrame(list_for_df_absolute, columns=columns_for_df)
    df_absolute.set_index('ITS2 type profile', drop=False, inplace=True)

    df_relative = pd.DataFrame(list_for_df_relative, columns=columns_for_df)
    df_relative.set_index('ITS2 type profile', drop=False, inplace=True)

    # at this point we have both of the dfs. We will use the relative df for getting the ordered smpl list
    # now go sample by sample find the samples max type and add to the dictionary where key is types, and value
    # is list of tups, one for each sample which is sample name and rel_abund of the given type
    type_to_sample_abund_dict = defaultdict.list()
    typeless_samples_list = []
    for i in range(7, 7 + len(listOfDataSetSamples)):
        sample_series = df_relative.iloc[:,i].astype('float')
        max_type_pos = sample_series.values.argmax()
        max_type_label = sample_series.idxmax()
        rel_abund_of_max_type = sample_series[max_type_label]
        if not rel_abund_of_max_type > 0:
            typeless_samples_list.append(sample_series.name)
        type_to_sample_abund_dict[max_type_label].append((sample_series.name, rel_abund_of_max_type))

    # here we have the dictionary populated. We can now go type by type according
    # to the sorted_analysis_type_abundance_list and put the samples that had the given type as their most abundant
    # type, into the sorted sample list, addtionaly sorted by how abund the type was in each of the samples
    samples_that_have_been_sorted = []
    # we are only concerned with the types that had samples that had them as most abundant
    for an_type_name in [at.name for at in sorted_analysis_type_abundance_list if at.name in type_to_sample_abund_dict.keys()]:
        samples_that_have_been_sorted.extend([a[0] for a in sorted(type_to_sample_abund_dict[an_type_name], key=lambda x: x[1], reverse=True)])

    # here we should have a list of samples that have been sorted according to the types they were found to
    # have as their most abundant
    # now we just need to add the samples that didn't have a type in them to be associated to. Negative etc.
    samples_that_have_been_sorted.extend(typeless_samples_list)

    # # now that we have the sorted list of types
    # # we want to extract some information
    # # For each sample, we want to have the abundance of each of the types in that sample in the order of the types
    # print('Collecting type abundances for samples:')
    # dataSetSample_to_type_abundance_dict = dict()
    # for i in range(len(listOfDataSetSamples)):
    #     sys.stdout.write('\r{}'.format(listOfDataSetSamples[i]))
    #     temp_abundance_list = []
    #     for analysis_type_obj in sorted_analysis_type_abundance_list:
    #         temp_abundance = int(across_clade_type_sample_abund_dict[analysis_type_obj][0].split('\t')[7:-1][i])
    #         temp_abundance_list.append(temp_abundance)
    #     dataSetSample_to_type_abundance_dict[listOfDataSetSamples[i]] = temp_abundance_list
    #
    #
    # # list of samples that have been added to the sorted list
    # samples_that_have_been_sorted = []
    #
    # # for each analysis type in order of the most abundant types first
    # sys.stdout.write('\n\nSorting samples according to ITS2 type profile abundances\n')
    # for j in range(len(sorted_analysis_type_abundance_list)):
    #     sys.stdout.write('\rSorting {}'.format(sorted_analysis_type_abundance_list[j]))
    #     # go through all of the samples and if a sample has the type as its most abundant type
    #     # then add it and its relative abundance of the type to to the below list
    #     list_of_samples_containing_analysis_type = []
    #     for data_set_sample_obj in dataSetSample_to_type_abundance_dict.keys():
    #         if data_set_sample_obj not in samples_that_have_been_sorted:
    #             # see if the type in question is the most abundant type in the sample in question
    #             list_of_abundances = dataSetSample_to_type_abundance_dict[data_set_sample_obj]
    #             abundance_of_type_in_sample = list_of_abundances[j]
    #             if abundance_of_type_in_sample != 0:
    #                 # then we need to see if this type is the most abundant in the sample
    #                 if all(i <= abundance_of_type_in_sample for i in list_of_abundances):
    #                     # then this type is the most abundant type in the sample and it should be added to the
    #                     # list_of_samples_containing_analysis_type
    #                     list_of_samples_containing_analysis_type.append(
    #                         (data_set_sample_obj, abundance_of_type_in_sample / sum(list_of_abundances)))
    #     # now that we've been through all of the samples for the given type
    #     # we need to sort the list_of_samples_containing_analysis_type by the second value of the tuple
    #     # we then need to add the samples in this order to the sample_that_have_been_sorted and remove them
    #     # from the samples_still_to_be_sorted
    #     sorted_list_of_samples_containing_analysis_type = sorted(list_of_samples_containing_analysis_type,
    #                                                              key=lambda x: x[1], reverse=True)
    #     sorted_list_of_samples_for_this_type = [a[0] for a in sorted_list_of_samples_containing_analysis_type]
    #
    #     for data_set_sample_obj in sorted_list_of_samples_for_this_type:
    #         samples_that_have_been_sorted.append(data_set_sample_obj)
    #
    # # here the samples_that_have_been_sorted should be populated in the sorted order that we want
    #
    # # Here the samples we should be left with should be samples that don't have a type associated to them
    # # negatives etc.
    # samples_not_ordered = list(set(listOfDataSetSamples) - set(samples_that_have_been_sorted))
    #
    # # add these samples to the list
    # samples_that_have_been_sorted.extend(samples_not_ordered)

    # # these are the tables that the ITS2 type data will be stored in.
    # outputTableOne.append(
    #     'ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format(
    #         '\t'.join([dataSamp.name for dataSamp in samples_that_have_been_sorted])))
    # outputTableTwo.append(
    #     'ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format(
    #         '\t'.join([dataSamp.name for dataSamp in samples_that_have_been_sorted])))

    # rearange the sample columns so that they are in the new order
    new_cols = ('ITS2 type profile UID\tClade\tMajority ITS2 sequence\tAssociated species\tITS2 type abundance local\tITS2 type abundance DB\tITS2 type profile\t{0}\tSequence accession / SymPortal UID\tAverage defining sequence proportions and [stdev]'.format(
            '\t'.join([dataSamp.name for dataSamp in samples_that_have_been_sorted]))).split('\t')
    df_absolute = df_absolute[new_cols]
    df_relative = df_relative[new_cols]

    # transpose
    df_absolute = df_absolute.T
    df_relative = df_relative.T

    # # at this point we need to rearrange the data that is held in the current row strings
    # for analysis_type_obj in across_clade_sorted_type_order:
    #     current_abundance_list_counts = across_clade_type_sample_abund_dict[analysis_type_obj][0].split('\t')[7:-1]
    #     current_abundance_list_props = across_clade_type_sample_abund_dict[analysis_type_obj][1].split('\t')[7:-1]
    #     new_abundance_list_counts = []
    #     new_abundance_list_props = []
    #     for new_sample in samples_that_have_been_sorted:
    #         new_abundance_list_counts.append(current_abundance_list_counts[listOfDataSetSamples.index(new_sample)])
    #         new_abundance_list_props.append(current_abundance_list_props[listOfDataSetSamples.index(new_sample)])
    #
    #     new_string_counts = '{}\t{}\t{}'.format(
    #         '\t'.join(across_clade_type_sample_abund_dict[analysis_type_obj][0].split('\t')[:7]),
    #         '\t'.join(new_abundance_list_counts),
    #         across_clade_type_sample_abund_dict[analysis_type_obj][0].split('\t')[-1])
    #     new_sting_props = '{}\t{}\t{}'.format(
    #         '\t'.join(across_clade_type_sample_abund_dict[analysis_type_obj][1].split('\t')[:7]),
    #         '\t'.join(new_abundance_list_props),
    #         across_clade_type_sample_abund_dict[analysis_type_obj][1].split('\t')[-1])
    #
    #     outputTableOne.append(new_string_counts)
    #     outputTableTwo.append(new_sting_props)

    # At this point we should have gone through each of the clades in order and for each clade we should have
    # created an analysis_type row in order of the sortedListOfTypes.

    # The current out put has samples on the Y and types on the X.
    # We want to reverse this as there will often be more samples than types
    # https://stackoverflow.com/questions/6473679/transpose-list-of-lists

    # # Put into list of lists format to transpose
    # newOutOne = [a.split('\t') for a in outputTableOne]
    # newOutTwo = [a.split('\t') for a in outputTableTwo]
    #
    # # Transpose
    # outputTableOne = list(map(list, zip(*newOutOne)))
    # outputTableTwo = list(map(list, zip(*newOutTwo)))
    #
    # # Put back into tab delim format
    # outputTableOne = ['\t'.join(a) for a in outputTableOne]
    # outputTableTwo = ['\t'.join(a) for a in outputTableTwo]

    outputDir = os.path.join(os.path.dirname(__file__), 'outputs/analyses/{}'.format(analysisObj.id))
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
    df_absolute.append('Species references')
    df_relative.append('Species references')

    # now add the references for each of the associated species
    for species in species_set:
        if species in species_ref_dict.keys():
            df_absolute.append([species, species_ref_dict[species]])
            df_relative.append([species, species_ref_dict[species]])

    # Now append the meta infromation for the output. i.e. the user running the analysis or the standalone
    # and the data_set submissions used in the analysis.
    # two scenarios in which this can be called
    # from the analysis: call_type = 'analysis'
    # or as a stand alone: call_type = 'stand_alone'

    if call_type == 'analysis':
        meta_info_string_items = ['Output as part of data_analysis ID: {}'.format(analysisobj.id),
                                  'Number of data_set objects as part of analysis = {}'.format(len(querySetOfDataSubmissions)),
                                  'submitting_user: {}'.format(analysisobj.submittingUser),
                                  'time_stamp {}'.format(analysisobj.timeStamp)]
        df_absolute.append(meta_info_string_items)
        df_relative.append(meta_info_string_items)
        for data_set_object in querySetOfDataSubmissions:
            data_set_meta_list = ['Data_set ID: {}'.format(data_set_object.id),
                                          'submitting_user: {}'.format(data_set_object.submittingUser),
                                          'time_stamp: {}'.format(data_set_object.timeStamp)]
            df_absolute.append(data_set_meta_list)
            df_relative.append(data_set_meta_list)
    else:
        # call_type=='stand_alone'
        meta_info_string_items = [
            'Stand_alone output by {} on {}'.format(output_user, str(datetime.now())),
            'data_analysis ID: {}'.format(analysisobj.id),
            'Number of data_set objects as part of output = {}'.format(len(querySetOfDataSubmissions))]

        for data_set_object in querySetOfDataSubmissions:
            meta_info_string_items=[
                'Data_set ID: {}'.format(data_set_object.id),
                'submitting_user: {}'.format(data_set_object.submittingUser),
                'time_stamp: {}'.format(data_set_object.timeStamp)]
            df_absolute.append(meta_info_string_items)
            df_relative.append(meta_info_string_items)

    date_time_string = str(datetime.now())
    path_to_profiles_absolute = '{}/{}_{}_{}.profiles.absolute.txt'.format(outputDir, analysisObj.id, analysisObj.name, date_time_string)
    df_absolute.to_csv(path_to_profiles_absolute, sep="\t")
    # writeListToDestination(path_to_profiles_absolute, outputTableOne)
    output_files_list.append(path_to_profiles_absolute)

    path_to_profiles_rel = '{}/{}_{}_{}.profiles.relative.txt'.format(outputDir, analysisObj.id, analysisObj.name, date_time_string)
    df_relative.to_csv(path_to_profiles_rel, sep="\t")
    # writeListToDestination(path_to_profiles_rel, outputTableTwo)
    output_files_list.append(path_to_profiles_rel)

    # ########################## ITS2 INTRA ABUND COUNT TABLE ################################
    div_output_pre_analysis_new_meta_and_new_dss_structure(datasubstooutput=datasubstooutput,
                                                           numProcessors=numProcessors, output_dir=outputDir,
                                                           sorted_sample_list=samples_that_have_been_sorted,
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

        sys.stdout.write('\rprocessing ITS2 type profile: {}'.format(anType))

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
        # TODO do we really have to go through every sample? I don't think so.
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
                                                           sorted_sample_list=None, analysis_obj_id=None, output_user=None ):
    # print out the user and time stamp and dataID if from a submission
    # This function will take in a list of dataSubmissions from which
    # to output the ITS2 sequence information of the samples

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
    sub_clade_list = [a for a in clade_list if a in cladeSet]

    sys.stdout.write('\n')
    for i in range(len(sub_clade_list)):

        # refSeqs of clade
        refSeqsOfClade = refSeqsInDSs.filter(clade=sub_clade_list[i])

        # dataSetSampleSequencesOfClade
        dataSetSampleSequencesOfClade = data_set_sample_sequence.objects.filter(
            data_set_sample_from__dataSubmissionFrom__in=querySetOfDataSubmissions,
            referenceSequenceOf__in=refSeqsOfClade)

        # This counter dict will have key = sequence name (either the sequence name, if named, or 'id_clade' if not)
        # value = the abundnace of the reference sequence across all of the samples in the output
        intraCounterManager = Manager()
        intraCounter = intraCounterManager.dict()

        dsssQueue = Queue()

        for dsss in dataSetSampleSequencesOfClade:
            dsssQueue.put(dsss)

        for N in range(numProcessors):
            dsssQueue.put('STOP')

        allProcesses = []

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        print('\nCreating refSeq counts for clade {}'.format(sub_clade_list[i]))
        for N in range(numProcessors):
            p = Process(target=outputWorkerTwo, args=(dsssQueue, intraCounter))
            allProcesses.append(p)
            p.start()

        for p in allProcesses:
            p.join()

        # Sort counter dict by abundance
        sortedIntraCounter = sorted(intraCounter.items(), key=operator.itemgetter(1), reverse=True)

        # Add the refSeqNames to the master cladeAbundanceOrderedRefSeqList
        cladeAbundanceOrderedRefSeqList.extend([a[0] for a in sortedIntraCounter])

    # Create the two output tables as lists
    intraAbundCountTable = []
    intraAbundPropTable = []
    # Create the header line one for prop one for count
    # add the cladal noName categories directly to the header
    # do not add them to the cladeAbundanceOrderedRefSeqList so that we don't append 0s to them
    # we will add the abundances in the outputWorkerThree method below

    # we will put together the headers piece by piece

    # 1 - the sample header and the named sequence headers
    headerPre = 'Samples\t{}'.format('\t'.join([a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]))
    noNameStrings = '\t'.join(['noName Clade {}'.format(cl) for cl in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']])
    qc_stats = '\t'.join(['raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                          'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                          'size_screening_violation_absolute', 'size_screening_violation_unique',
                          'post_taxa_id_absolute_non_symbiodinium_seqs',
                          'post_taxa_id_unique_non_symbiodinium_seqs', 'post_med_absolute', 'post_med_unique'])
    # append the noName sequences as individual sequence abundances
    break_down_strings = '\t'.join([a for a in cladeAbundanceOrderedRefSeqList if '_' in a])
    output_header = '\t'.join([headerPre, noNameStrings, qc_stats, break_down_strings])

    intraAbundCountTable.append(output_header)
    intraAbundPropTable.append(output_header)
    ######################################################################################

    ############ POPULATE TABLE WITH CC DATA #############
    # For every sample
    sampleList = data_set_sample.objects.filter(dataSubmissionFrom__in=querySetOfDataSubmissions)

    # In order to MP this we will have to pay attention to the order. As we can't keep the order as we work with the
    # MPs we will do as we did above for the profies outputs and have an output dict that will have the sample as key
    # and its corresponding data row as the value. Then when we have finished the MPing we will go through the
    # sampleList order and output the data in this order.
    # So we will need a managed output dict.

    sampleOutputDictManager = Manager()
    managedSampleOutputDict = sampleOutputDictManager.dict()

    cladeAbundanceOrderedRefSeqListManager = Manager()
    cladeAbundanceOrderedRefSeqList = cladeAbundanceOrderedRefSeqListManager.list(cladeAbundanceOrderedRefSeqList)

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
            dssQueue, managedSampleOutputDict, cladeAbundanceOrderedRefSeqList, querySetOfDataSubmissions))
        allProcesses.append(p)
        p.start()

    for p in allProcesses:
        p.join()

    # Now before writing out we need to add to the intraAbund objects from the managedSampleOutputDict in order of
    # sampleList
    #### DEVELOPMENT ####
    # So that we can inspect the managedSampleOutputDict
    managedSampleOutputDict_dict = dict(managedSampleOutputDict)

    #####################

    # If there is a sorted sample list, make sure that it matches the samples that we are outputting

    if sorted_sample_list:
        sorted_sample_list_names = [a.name for a in sorted_sample_list]

        if len(sorted_sample_list_names) != len(sampleList):
            sys.exit('Name in sorted_sample_list do not match those to be outputted!')

        provided_name_to_samp_name_dict = {}
        for provided_name in sorted_sample_list_names:
            match_list = []
            for smp in sampleList:
                if provided_name == smp.name:
                    match_list.append(smp.name)
                    provided_name_to_samp_name_dict[provided_name] = smp
            if len(match_list) > 1:
                sys.exit('Sample name {} matches more than one output sample ({}).'.format(provided_name,
                                                                                           '\t'.join(match_list)))
            if len(match_list) == 0:
                sys.exit('Sample name {} does not match any output sample.'.format(provided_name))

        # if we got to here then the sorted_sample_list looks good

        for provided_name in sorted_sample_list_names:
            dss = provided_name_to_samp_name_dict[provided_name]
            intraAbundCountTable.append(managedSampleOutputDict[dss][0])
            intraAbundPropTable.append(managedSampleOutputDict[dss][1])
    else:
        # this returns a list which is simply the names of the samples
        # for more info on how this is ordered look at the comment in the method
        ordered_sample_list = generate_ordered_sample_list(managedSampleOutputDict, output_header)

        # managedSampleOutputDict currently has data_set_sample objects as the key. This is excessive and the
        # above list returns strings that are the dss object names
        # so quickly convert these keys to just the names
        managedSampleOutputDict = {k.name: v for k, v in managedSampleOutputDict.items()}

        for dss in ordered_sample_list:
            intraAbundCountTable.append(managedSampleOutputDict[dss][0])
            intraAbundPropTable.append(managedSampleOutputDict[dss][1])

    # Now add the accesion number / UID for each of the DIVs
    accession_list = []
    accession_list.append('DIV_accession')
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
        ref_seq = reference_sequence.objects.get(name=ref_seq_name)
        if ref_seq.accession and ref_seq.accession != 'None':
            ref_seq_accession = ref_seq.accession
        else:
            ref_seq_accession = str(ref_seq.id)
        accession_list.append(ref_seq_accession)

    accession_str = '\t'.join(accession_list)

    # now add the blanks for the noName cladal columns
    accession_str += ''.join(['\t' for cld in clade_list])
    # add blanks for the qc stats
    accession_str += ''.join(['\t' for i in range(11)])

    # now add the accession numbers for the noName break down seqs
    accession_list = []
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
        accession_list.append(ref_seq_name.split('_')[0])

    accession_str = '\t'.join([accession_str, '\t'.join(accession_list)])

    intraAbundCountTable.append(accession_str)
    intraAbundPropTable.append(accession_str)

    # Now append the meta infromation for each of the data_sets that make up the output contents
    # this is information like the submitting user, what the IDs of the datasets are etc.
    # There are several ways that this can be called.
    # it can be called as part of the submission: call_type = submission
    # part of an analysis output: call_type = analysis
    # or stand alone: call_type = 'stand_alone'
    # we should have an output for each scenario
    if call_type=='submission':
        meta_info_string_items = []
        data_set_object = querySetOfDataSubmissions[0]
        # there will only be one data_set object
        meta_info_string = 'Output as part of data_set submission ID: {}\tsubmitting_user: {}\ttime_stamp: {}\t'.format(data_set_object.id, data_set_object.submittingUser, data_set_object.timeStamp)
        meta_info_string_items.append(meta_info_string)
    elif call_type=='analysis':
        data_analysis_obj = data_analysis.objects.get(id=analysis_obj_id)
        meta_info_string_items = [
            'Output as part of data_analysis ID: {}\t'
            'Number of data_set objects as part of analysis = {}\t'
            'submitting_user: {}\ttime_stamp {}'.format(
                data_analysis_obj.id, len(querySetOfDataSubmissions), data_analysis_obj.submittingUser, data_analysis_obj.timeStamp)]
        for data_set_object in querySetOfDataSubmissions:
            meta_info_string_items.append('Data_set ID: {}\tsubmitting_user: {}\ttime_stamp: {}'.format(data_set_object.id, data_set_object.submittingUser, data_set_object.timeStamp))
    else:
        # call_type=='stand_alone'
        meta_info_string_items = [
            'Stand_alone output by {} on {}\tNumber of data_set objects as part of output = {}'.format(output_user, str(datetime.now()), len(querySetOfDataSubmissions))]
        for data_set_object in querySetOfDataSubmissions:
            meta_info_string_items.append('Data_set ID: {}\tsubmitting_user: {}\ttime_stamp: {}'.format(data_set_object.id, data_set_object.submittingUser, data_set_object.timeStamp))

    # now add the meta information onto the end of the output table lists
    intraAbundCountTable.extend(meta_info_string_items)
    intraAbundPropTable.extend(meta_info_string_items)



    # Here we have the tables populated and ready to output
    if analysis_obj_id:
        path_to_div_absolute = '{}/{}_{}.DIVs.absolute.txt'.format(output_dir, analysis_obj_id,
                                                                   '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        path_to_div_absolute = '{}/{}.DIVs.absolute.txt'.format(output_dir,
                                                                '_'.join([str(a) for a in dataSubmissionsToOutput]))
    writeListToDestination(path_to_div_absolute, intraAbundCountTable)
    output_path_list.append(path_to_div_absolute)

    if analysis_obj_id:
        path_to_div_relative = '{}/{}_{}.DIVs.relative.txt'.format(output_dir, analysis_obj_id,
                                                                   '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        path_to_div_relative = '{}/{}.DIVs.relative.txt'.format(output_dir,
                                                                '_'.join([str(a) for a in dataSubmissionsToOutput]))

    writeListToDestination(path_to_div_relative, intraAbundPropTable)
    output_path_list.append(path_to_div_relative)

    # Output a .fasta for of all of the sequences found in the analysis
    fasta_output_list = []
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
        ref_seq = reference_sequence.objects.get(name=ref_seq_name)
        fasta_output_list.append('>{}'.format(ref_seq_name))
        fasta_output_list.append(ref_seq.sequence)
    for ref_seq_name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
        ref_seq = reference_sequence.objects.get(id=int(ref_seq_name.split('_')[0]))
        fasta_output_list.append('>{}'.format(ref_seq_name))
        fasta_output_list.append(ref_seq.sequence)

    if analysis_obj_id:
        fasta_path = '{}/{}_{}.DIVs.fasta'.format(output_dir, analysis_obj_id,
                                                  '_'.join([str(a) for a in dataSubmissionsToOutput]))
    else:
        fasta_path = '{}/{}.DIVs.fasta'.format(output_dir, '_'.join([str(a) for a in dataSubmissionsToOutput]))
    writeListToDestination(fasta_path, fasta_output_list)
    output_path_list.append(fasta_path)

    print('\nITS2 sequence output files:')
    for path_item in output_path_list:
        print(path_item)

    return output_path_list
def outputWorkerTwo(input, outDict):
    for dsss in iter(input.get, 'STOP'):
        if not dsss.referenceSequenceOf.hasName:
            name_unit = str(dsss.referenceSequenceOf.id) + '_{}'.format(dsss.referenceSequenceOf.clade)
        else:
            name_unit = dsss.referenceSequenceOf.name
        if name_unit in outDict.keys():
            outDict[name_unit] = outDict[name_unit] + dsss.abundance
        else:
            outDict[name_unit] = dsss.abundance
        sys.stdout.write('\rCounting seqs for {}'.format(name_unit))

def generate_ordered_sample_list(managedSampleOutputDict, output_header):
    # create a df from the managedSampleOutputDict
    two_d_list_for_df = [list_element[1].split('\t') for list_element in managedSampleOutputDict.values()]
    output_df_relative = pd.DataFrame(two_d_list_for_df, columns=output_header.split('\t'))
    # convert the string elements to floats
    # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.astype.html
    # put the sampe names as the index, drop this column and then convert the reamining cols to float
    # n.b the convert to float fails if there are any non convertable elements in the df. Even if error warning are
    # ignored
    output_df_relative = output_df_relative.set_index(keys='Samples', drop=True).astype('float')
    non_seq_columns = ['raw_contigs', 'post_taxa_id_absolute_non_symbiodinium_seqs', 'post_qc_absolute_seqs',
                       'post_qc_unique_seqs',
                       'post_taxa_id_unique_non_symbiodinium_seqs', 'post_taxa_id_absolute_symbiodinium_seqs',
                       'post_taxa_id_unique_symbiodinium_seqs', 'post_med_absolute', 'post_med_unique',
                       'size_screening_violation_absolute', 'size_screening_violation_unique']
    noName_seq_columns = ['noName Clade {}'.format(clade) for clade in list('ABCDEFGHI')]
    cols_to_drop = non_seq_columns + noName_seq_columns
    sequence_only_df_relative = output_df_relative.drop(columns=cols_to_drop)
    ordered_sample_list = get_sample_order_from_rel_seq_abund_df(sequence_only_df_relative)
    return ordered_sample_list

def get_sample_order_from_rel_seq_abund_df(sequence_only_df_relative):
    max_seq_ddict = defaultdict(int)
    seq_to_samp_dict = defaultdict(list)
    # for each sample get the columns name of the max value of a div not including the columns in the following:
    for sample_to_sort in sequence_only_df_relative.index.values.tolist():
        max_abund_seq = sequence_only_df_relative.loc[sample_to_sort].idxmax()
        max_rel_abund = sequence_only_df_relative.loc[sample_to_sort].max()
        # add a tup of sample name and rel abund of seq to the seq_to_samp_dict
        seq_to_samp_dict[max_abund_seq].append((sample_to_sort, max_rel_abund))
        # add this to the ddict count
        max_seq_ddict[max_abund_seq] += 1
    # then once we have compelted this for all sequences go clade by clade
    # and generate the sample order
    ordered_sample_list = []
    for clade in list('ABCDEFGHI'):
        tup_list_of_clade = []
        # get the clade specific list of the max_seq_ddict
        for k, v in max_seq_ddict.items():
            if k.startswith(clade) or k[-2:] == '_{}'.format(clade):
                tup_list_of_clade.append((k, v))
        if not tup_list_of_clade:
            continue
        # now get an ordered list of the sequences for this clade
        ordered_sequence_of_clade_list = [x[0] for x in sorted(tup_list_of_clade, key=lambda x: x[1], reverse=True)]

        for seq_to_order_samples_by in ordered_sequence_of_clade_list:
            tup_list_of_samples_that_had_sequence_as_most_abund = seq_to_samp_dict[seq_to_order_samples_by]
            ordered_list_of_samples_for_seq_ordered = \
                [x[0] for x in
                 sorted(tup_list_of_samples_that_had_sequence_as_most_abund, key=lambda x: x[1], reverse=True)]
            ordered_sample_list.extend(ordered_list_of_samples_for_seq_ordered)
    return ordered_sample_list

def outputWorkerThree_pre_analysis_new_dss_structure(input, outDict, cladeAbundanceOrderedRefSeqList,
                                                     querySetOfDataSubmissions):
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    for dss in iter(input.get, 'STOP'):

        sys.stdout.write('\rOutputting DIV data for {}'.format(dss.name))
        # List that will hold the row
        sampleRowDataCounts = []
        sampleRowDataProps = []
        cladalAbundances = [int(a) for a in json.loads(dss.cladalSeqTotals)]
        sampleSeqTot = sum(cladalAbundances)

        if dss.errorInProcessing or sampleSeqTot == 0:
            # Then this sample had a problem in the sequencing and we need to just output 0s across the board
            # Named sequences get 0s
            for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # Clade totals get 0s
            for cl in cladeList:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
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

            # Breakdown details get 0s
            for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)
            # Here we need to add the string to the outputDict rather than the intraAbund table objects
            countData = '{}\t{}'.format(dss.name, '\t'.join(str(x) for x in sampleRowDataCounts))
            propData = '{}\t{}'.format(dss.name, '\t'.join("{:.3f}".format(x) for x in sampleRowDataProps))
            outDict[dss] = [countData, propData]
            continue

        # Get number of seqs found in Sample
        dsssInSample = data_set_sample_sequence.objects.filter(data_set_sample_from=dss)

        dataSetSampleSeqDict = defaultdict(int)
        # this dict will hold noName seqs abundances separated by clade.
        # when we come to populate the info on the sequences that were not in clade collections then we will
        # check to see if there is info for the clade in question in this dict first.
        # if the clade in question isn't in this dict then we will use the information from the
        # samples cladalSeqTotals parameter
        noNameCladalHolder = defaultdict(int)
        no_name_break_down_dict = defaultdict(int)
        for dsss in dsssInSample:
            ref_seq_obj = dsss.referenceSequenceOf
            ref_seq_obj_name = ref_seq_obj.name
            ref_seq_obj_clade = ref_seq_obj.clade
            if ref_seq_obj.hasName:
                dataSetSampleSeqDict[ref_seq_obj_name] += dsss.abundance
            else:
                noNameCladalHolder[ref_seq_obj_clade] += dsss.abundance
                no_name_break_down_dict[str(ref_seq_obj.id) + '_{}'.format(
                    ref_seq_obj_clade)] += dsss.abundance  # populate row in order of the master sorted intras
        # Do not populate any of the noName categories

        for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' not in a]:
            if name in dataSetSampleSeqDict.keys():
                counts = dataSetSampleSeqDict[name]
                sampleRowDataCounts.append(counts)
                sampleRowDataProps.append(counts / sampleSeqTot)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)

        # for each clade in cladeList, just look to see if there is data.
        for clade in cladeList:
            if clade in noNameCladalHolder.keys():
                sampleRowDataCounts.append(noNameCladalHolder[clade])
                sampleRowDataProps.append(noNameCladalHolder[clade] / sampleSeqTot)
            else:
                # then there were no no_names of this clade
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)

        # Here we add in the post qc and post-taxa id counts
        # For the absolute counts we will report the absolute seq number
        # For the relative counts we will report these as proportions of the sampleSeqTot.
        # I.e. we will have numbers larger than 1 for many of the values and the symbiodinium seqs should be 1

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

        # now add the noName clade breakdown sequence counts
        for name in [a for a in cladeAbundanceOrderedRefSeqList if '_' in a]:
            if name in no_name_break_down_dict.keys():
                counts = no_name_break_down_dict[name]
                sampleRowDataCounts.append(counts)
                sampleRowDataProps.append(counts / sampleSeqTot)
            else:
                sampleRowDataCounts.append(0)
                sampleRowDataProps.append(0)

        # Here we need to add the string to the outputDict rather than the intraAbund table objects
        countData = '{}\t{}'.format(dss.name, '\t'.join(str(x) for x in sampleRowDataCounts))
        propData = '{}\t{}'.format(dss.name, '\t'.join("{:.3f}".format(x) for x in sampleRowDataProps))
        outDict[dss] = [countData, propData]