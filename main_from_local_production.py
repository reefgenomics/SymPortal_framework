#!/usr/bin/env /home/humebc/miniconda3/bin/python3.6

import argparse
import subprocess
import itertools
from datetime import datetime
# Django specific settings
import os
from datetime import datetime
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.conf import settings
######## Setup Django DB and Models ########
# Ensure settings are read
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()

# Your application specific imports
from dbApp.models import symportal_framework, data_set, reference_sequence, data_set_sample_sequence, analysis_type, analysis_group, data_set_sample, data_analysis, clade_collection, clade_collection_type
############################################

import dataSubCollectionRun_from_local_production


import createDataSubmission_from_local_production
from collections import defaultdict
import json


###### Generic functions ######
def readDefinedFileToList(filename):
    tempList = []
    with open(filename, mode='r') as reader:
        tempList = [line.rstrip() for line in reader]
    return tempList

def writeListToDestination(destination, listToWrite):
    #print('Writing list to ' + destination)
    try:
        os.makedirs(os.path.dirname(destination))
    except FileExistsError:
        pass

    with open(destination, mode='w') as writer:
        i = 0
        while i < len(listToWrite):
            if i != len(listToWrite)-1:
                writer.write(listToWrite[i] + '\n')
            elif i == len(listToWrite)-1:
                writer.write(listToWrite[i])
            i += 1
###############################



def main():
    '''
    Data submissions on the server db as of 22.04.18:
    3:20170506_roberto
    15:20171114_mellisa
    4:20171009_roberto
    5:20171017_roberto
    16:20171008_jitErn
    6:20171114_roberto
    7:20170713_alejandro
    17:20171017_tullia
    8:20170719_alejandro
    9:20171114_alejandro
    18:20180117_burkepile
    10:UTS_ONE
    11:UTS_TWO
    19:20171123_caroline
    12:Emily
    13:Voolstra
    20:ryan_first
    14:ed_24_07_17_singapore
    22:20180120_forcioli
    35:120418_buitrago
    23:20171203_ed
    25:EdData
    50:180418_env_res_non_coral
    26:20180316_jit_ern_class
    36:20171025_burkepile
    34:110418_gong
    52:20180418_jit_ern
    56:180418_env_res_coral
    53:240418_buitrago
    54:20180417_18_ryan



'3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,22,23,25,26,35,36,49,52'




    DataSubmissions on the local:
    90:20170506_roberto_SbyS
    91:20171009_roberto_SbyS
    92:20171017_roberto_SbyS
    93:20171114_roberto_SbyS
    94:20170713_alejandro_SbyS
    95:20170719_alejandro_SbyS
    96:20171114_alejandro_SbyS

    97:UTS_ONE_SbyS
    98:UTS_TWO_SbyS
    99:Emily_SbyS
    100:Voolstra_SbyS
    101:ed_24_07_17_singapore_SbyS
    102:20171114_mellisa_SbyS
    103:20171008_jitErn_SbyS

    104:20171017_tullia_SbyS
    105:20171025_burkepile_SbyS
    106:20171123_caroline_SbyS
    107:ryan_first_SbyS
    109:20171203_ed
    114:20180117_burkepile

    115:EdData (DBV1.11)
    20180316_jit_ern_class:116 (DBV1.12)
'''

    # env_res_code.figure_code_two()

    # wthCladeCutOff = 0.03
    # cores = 20
    # # #
    #
    dataSubCollectionRun_from_local_production.formatOutput_ord(data_analysis.objects.get(id=28), numProcessors=20, datasubstooutput='111')
    # # apples = 'juice'
    # #
    # #
    # customDataSetIDs = '3,4,5,6,7,8,9,10,11,13,14,16,17,18,19,20,23,25,26,35,36,52,53,54,87,90,110,111'
    # # try to recreate the DBV from the local DB on this postgres db for the MS
    # # customDataSetIDs = '3,4,5,6,7,8,9,10,11,12,13,14,15,16,36,19,20,25'
    # # DBV_1.16 = '3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,22,23,25,26,35,36,52,53,54,87,90'
    # # 4 x roberto, 3 x alejandro
    # newAnalysisObj = data_analysis(listOfDataSubmissions=str(customDataSetIDs),
    #                                  withinCladeCutOff=float(wthCladeCutOff), name='Emily_new_dss_struct')
    # newAnalysisObj.description = 'same as DBV1.16 only with didier and Emily data with the new dss_structure'
    # newAnalysisObj.save()
    # # # # # # # #
    # # # # # # # # # already_started_analysis = data_analysis.objects.get(id=306)
    # # # # # # # #
    # dataSubCollectionRun_from_local_production.main(newAnalysisObj, cores)
    #
    # # # # # # #
    # print('return code: 0\nAnalysis complete')



    #
    # # # # # ###### DEBUG TO RUN A DBV CREATION without the commandline #####
    # listOfCores = [6, 10, 20]
    # holder = {core:0 for core in listOfCores}
    #
    #
    # for core in listOfCores:
    #     for i in range(20):
    #         wthCladeCutOff = 0.03
    #         cores = core
    #         customDataSetIDs = '105'
    #
    #         newAnalysisObj = data_analysis(listOfDataSubmissions=str(customDataSetIDs),
    #                                          withinCladeCutOff=float(wthCladeCutOff), name='BurkepileTest{}.{}'.format(cores, i))
    #         newAnalysisObj.description = 'This is debugging to see why Burkepile 42g didnt group'
    #         newAnalysisObj.save()
    #         dataSubCollectionRun_from_local_production.main(newAnalysisObj, cores)
    #
    #         print('return code: 0\nAnalysis complete')
    #
    #         # HERE get the summary stats
    #         for at in analysis_type.objects.filter(dataAnalysisFrom=newAnalysisObj):
    #             # here compare the number of clade collections got from the method
    #             lenMethod = len(at.getCladeCollections())
    #             lenCCTs = len(clade_collection_type.objects.filter(analysisTypeOf=at))
    #             if lenMethod != lenCCTs:
    #                 holder[core] += 1
    #
    # for key in holder.keys():
    #     print('core {}: mismatch {}'.format(key, holder[key]))
    # apples = 'pears'



    # # # #
    # dataSubCollectionRun_from_local_production.test()
    # sample_list = readDefinedFileToList('/home/humebc/phylogeneticSoftware/SymPortal_interim_250318/outputs/25/25_90_sorted_sample_list.txt')
    # sample_list = [data_set_sample.objects.get(id=int(a)).name for a in sample_list]
    # outputdir = os.path.join(os.path.dirname(__file__), 'outputs/non_analysis')
    # # # # #
    # # # # dataSubCollectionRun_from_local_production.div_output_pre_analysis(datasubstooutput='91', numProcessors=16, output_dir=outputdir, sorted_sample_list=env_res_code.return_ordered_sample_list())
    # # # # no sorted sample list
    # dataSubCollectionRun_from_local_production.div_output_pre_analysis_new_meta_and_new_dss_structure(datasubstooutput='109', numProcessors=20, output_dir=outputdir)
    # dataSubCollectionRun_from_local_production.div_output_pre_analysis('50,56', 20, env_res_code.return_ordered_sample_list())

    # # # # #################################################################
    #
    # ####### DEBUG PRINT OUTPUT #####
    # ################################
    #
    # # # ######## DEBUG do a submission ####
    # list_of_dirs_to_process = ['/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20170506_roberto',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171009_roberto',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171017_roberto',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171114_roberto',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20170713_alejandro',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20170719_alejandro',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171114_alejandro',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/UTS_ONE',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/UTS_TWO',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/Emily',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/Voolstra',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/ed_24_07_17_singapore',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171114_mellisa',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171008_jitErn',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171017_tullia',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20180117_burkepile',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171123_caroline',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/ryan_first'
    #
    #                        ]



    # list_of_dirs_to_process = ['/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/Emily']
    #
    # for dir in list_of_dirs_to_process:
    #     symportal_framework_object = symportal_framework.objects.get(id=1)
    #     latest_reference_fasta = symportal_framework_object.latest_reference_fasta
    #     iteration_id = symportal_framework_object.next_reference_fasta_iteration
    #     symbiodinium_cutoff_value = symportal_framework_object.required_sub_e_value_seq_support_blast_symbiodinium
    #     # nameForDataSubmission = dir.split('/')[-1]
    #     # nameForDataSubmission = '{}_{}'.format(dir.split('/')[-1], '{}_{}'.format(iteration_id, symbiodinium_cutoff_value))
    #     nameForDataSubmission = '{}_{}'.format(dir.split('/')[-1], '2_2') + '_new_dss_structure'
    #     # path_to_name_file = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/ryan_first/ryan_first_name_list.txt'
    #
    #     newDS = data_set(name=nameForDataSubmission, timeStamp=str(datetime.now()), reference_fasta_database_used='symClade_2_2.fa')
    #     newDS.save()
    #     createDataSubmission_from_local_production.main(dir, newDS.id, 20)
    #     # createDataSubmission_from_local_production.main(dir, 114, 4)
    # bearss = 'berries'

    # dataSubCollectionRun_from_local_production.generate_within_clade_UniFrac_distances_ITS2_type_profiles(data_submission_id_str='25', data_analysis_id=23, num_processors=4, method='mothur', bootstrap_value=10)

    # N.B the fastree tree implementation of this unifrac calculation does not seem to work well
    # The method options are 'mothur', 'fasttree' and 'phylip'.

    # ###### Distance Matrices ######
    # PCoA_paths_list = dataSubCollectionRun_from_local_production.generate_within_clade_UniFrac_distances_ITS2_type_profiles(data_submission_id_str='25', data_analysis_id=26, num_processors=16, method='mothur', bootstrap_value=100)
    # PCoA_paths_list = dataSubCollectionRun_from_local_production.generate_within_clade_UniFrac_distances_samples(dataSubmission_str='110', num_processors=20, method='mothur', bootstrap_value=100)
    # print('distance calculation complete')
    # dataSubCollectionRun_from_local_production.plotPCoA_types(PCoA_paths_list[0])

    ##### DATABASE Iterations #####
    # createDataSubmission_from_local_production.screen_sub_e_value_sequences(data_sub_data_dir='/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/180418_env_res_non_coral', iteration_id=2, seq_sample_support_cut_off=2, previous_reference_fasta_name='symClade_1_2.fa')

    # ###### ENV datasubmitions profile searching ######
    # DBV_id = 25
    # env_data_sub_to_run = '86'
    # num_processors = 16
    # # profile_id_to_cc_id_list = env_res_code.screen_environmental_samples_for_ITS2_type_profiles(analysis_id=DBV_id, env_data_submission_id=env_data_sub_to_run, num_proc=num_processors)
    # with open('/home/humebc/phylogeneticSoftware/SymPortal_interim_250318/outputs/environmental_type_screening/86/screened_list_86.json', 'r') as f:
    #     profile_id_to_cc_id_list = json.load(f)
    # env_res_code.formatOutput_ord_environmental_dataSubmission(analysis_obj_id=DBV_id, num_proc=num_processors, data_env_sub_ids=env_data_sub_to_run, type_to_cc_list=profile_id_to_cc_id_list)
    # ##################################################
    # ###################################
    # ###### args ############
    # # THIS BELOW SECTION DOES WORK AND JUST NEEDS UNCOMMENTING TO IMPLEMENT PROPER COMMAND LINE FUNCTIONALITY
    # # For the time being however I would still like to debug through the creation of the dataSubmissions
    #
    # #These will later be set at the command line but for the time being we will set them here
    # parser = argparse.ArgumentParser(description='Intragenomic analysis of the ITS2 region of the nrDNA', epilog='For support email: support@symportal.org')
    # # #
    # group = parser.add_mutually_exclusive_group(required=True)
    # # choices = ['runSubCol', 'submitDataSet', 'viewDataSubs']
    # group.add_argument('--submit', metavar='path_to_dir', help='Run this to submit data to the SymPortal database. A '
    #                                                            'single argument giving the directory containing the '
    #                                                            'paired sequencing reads in .fastq.gz '
    #                                                            'format. Alternatively, the path can point directly to a'
    #                                                            ' single compressed file containing the '
    #                                                            'same paired fastq.gz files. A name can be associated with'
    #                                                            'the data_set using the --name flag. The number of'
    #                                                            'processes to use can also be specified using the '
    #                                                            '--cores flag')
    # group.add_argument('--r', metavar='--runAgainstDBV', help='Run this to run a data_set '
    #                                                                                 'against the latest database '
    #                                                                                 'version')
    # group.add_argument('--ds', action='store_true', help=' Print the details of data_set currently in the SymPortal database')
    # group.add_argument('--mdb', metavar='dbID', help='Analyse one or more data_set objects '
    #                                                                    'together to create a '
    #                                                                     'new database version. Enter comma separated '
    #                                                                    'IDs of the data_set IDs you which to '
    #                                                                    'analyse. e.g.: 43,44,45. If you wish to use all available dataSubmissions,'
    #                                                                    'you may pass \'all\' as an argument' )
    # #This --printOutput will need to be work on to separate doing an output that was part of --makeDBV vs --runAgainstDBV
    # # for the time being we will just to the --makeDBV output.
    # group.add_argument('--p', metavar='dsID', help='Use this function to output the '
    #                                                                                    'DIV and analysis_type tables for '
    #                                                                                    'a given set of dataSubmissions '
    #                                                                                    'that have been run against a given '
    #                                                                                    'databaseVersion. Give the dataSubmissions '
    #                                                                                    'that you wish to make outputs for as '
    #                                                                                    'arguments to the --printOutput flag. '
    #                                                                                    'Give the ID of the DBV you wish to '
    #                                                                                    'output these from using the --DBVersion '
    #                                                                                    'flag. To output for multiple dataSubmissions, '
    #                                                                                    'comma separate the IDs of the dataSubmissions, '
    #                                                                                    'e.g. 44,45,46.')
    # group.add_argument('--dv', action='store_true', help=' Print the details of current '
    #                                                                          'database versions')
    # parser.add_argument('--cores', type=int, help='Number of cores to use', default=3)
    # parser.add_argument('--name', help='A name for your input or analysis', default='noName')
    # parser.add_argument('--description', help='An optional description', default='No description')
    # parser.add_argument('--DBVersion', type=int, help='The ID of the DBV you wish to output from')
    # group.add_argument('--vacuumDatabase', action='store_true', help='Vacuuming the database will free up memory from '
    #                                                                  'objects that have been deleted recently')
    # args = parser.parse_args()
    #
    #
    #
    # if args.submit:
    #     if args.name == 'noName':
    #         print('Please provide a name using the --name flag. e.g. --name dataSubOne')
    #         return
    #
    #
    #     # This is how we read data into the SP database
    #     # This should be given a string on the command line which is a directory that contains the sequencing data.
    #     # We will need to create a data_set object for each time this command is run.
    #     numProcessors = args.cores
    #     #make it so that we can have the paired .fastq.gz files of a single .zip file
    #     pathToInputFile = args.submit
    #
    #     # pathToInputFile = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/UTS_TWO/UTS_TWO.zip'
    #     nameForDataSubmission = args.name
    #     newDS = data_set(name = nameForDataSubmission, timeStamp=str(datetime.now()))
    #     newDS.save()
    #     createDataSubmission_from_local_production.main(pathToInputFile, newDS.id, numProcessors)
    # elif args.mdb:
    #     if args.name == 'noName':
    #         print('Please provide a name using the --name flag. e.g. --name dataSubOne')
    #         return
    #
    #     wthCladeCutOff = 0.03
    #
    #     cores = args.cores
    #     customDataSetIDs = args.mdb
    #     if args.mdb == 'all':
    #         tempList = []
    #         ds = data_set.objects.all()
    #         for ds in data_set.objects.all():
    #             tempList.append(str(ds.id))
    #         stringList = ','.join(tempList)
    #         customDataSetIDs = stringList
    #     newAnalysisObj = data_analysis(listOfDataSubmissions=str(customDataSetIDs),
    #                                      withinCladeCutOff=float(wthCladeCutOff), name=args.name)
    #     newAnalysisObj.description = args.description
    #     newAnalysisObj.save()
    #     dataSubCollectionRun_from_local_production.main(newAnalysisObj, cores)
    #     print('return code: 0\nAnalysis complete')
    # elif args.p:
    #     # then we will rerun a previously run analysis_type that is a DBV.
    #     # TODO implement checks to make sure that the dataAnalysis ID given is a DBV rather than a run against DBV analysis
    #     # for time being skip this
    #     dataSubCollectionRun_from_local_production.formatOutput(data_analysis.objects.get(id=args.DBVersion), numProcessors=args.cores, datasubstooutput=args.p)
    # elif args.ds:
    #     # Then print out all of the dataSubmissions with names and IDs in the db
    #     for ds in data_set.objects.all():
    #         print('{}: {}'.format(ds.id, ds.name))
    # elif args.dv:
    # # Then print out all of the dataAnalysisTwos with names and IDs in the db
    #     for da in data_analysis.objects.all():
    #         print('{}: {}'.format(da.id, da.name))
    # elif args.vacuumDatabase:
    #     print('Vacuuming database')
    #     vacuum_db()
    #     print('Vacuuming complete')
    # ########################





    # Here we are going to bulk redo submissions by having a list that is a set of directories to be processed
    # listOfDirsToProcess = ['/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20170506_roberto',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171009_roberto',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171017_roberto',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171114_roberto',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20170713_alejandro',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20170719_alejandro',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171114_alejandro',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/UTS_ONE',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/UTS_TWO',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/Emily',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/Voolstra',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/ed_24_07_17_singapore',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171114_mellisa',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171008_jitErn',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171017_tullia',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171025_burkepile',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20171123_caroline',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/ryan_first'
    # ]

    # listOfDirsToProcess = ['/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/ed_24_07_17_singapore',
    #                        '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/mellisa']
    #
    # for dir in listOfDirsToProcess:
    #     # # # if creating a new data_set
    #     numProcessors = 24
    #     #make it so that we can have the paired .fastq.gz files of a single .zip file
    #     # pathToInputFile = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'SymPortal_Data/rawData/20171016_alejandro/'))
    #
    #     # pathToInputFile = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/UTS_TWO/UTS_TWO.zip'
    #     nameForDataSubmission = dir.split('/')[-1] + '_absoluteSize'
    #     newDS = data_set(name = nameForDataSubmission, timeStamp=str(datetime.now()))
    #     newDS.save()
    #     createDataSubmission_from_local_production.main(dir, newDS.id, numProcessors)
    #     # # createDataSubmission_from_local_production.main(pathToInputFile, 51, numProcessors)







    ############################################


    ######## Run collection of data_set ##########
    # This will be one of the major analysis types
    # It will run an analysis on a collection of datasubmissions, doing profile discovery from scratch from the CCs in
    # the data Submissions included in the analysis

    '''
    47: EdData
    65: 20171008_jitErn_absoluteSize
    66: 20170506_roberto_absoluteSize
    67: 20170719_alejandro_absoluteSize
    68: 20170713_alejandro_absoluteSize
    69: 20171009_roberto_absoluteSize
    71: 20171017_roberto_absoluteSize
    72: 20171017_tullia_absoluteSize
    73: ryan_first_absoluteSize
    74: UTS_ONE_absoluteSize
    75: UTS_TWO_absoluteSize
    76: Emily_absoluteSize
    77: Voolstra_absoluteSize
    79: ed_24_07_17_singapore_absoluteSize
    85: 20171025_burkepile
    '''


    # create new analysis object that will be run with dataSubCollectionRun
    # within Clade Cutoff
    # wthCladeCutOff = 0.03
    # #
    # cores = 20
    # customDataSetIDs = '36'
    # newAnalysisObj = data_analysis(listOfDataSubmissions=str(customDataSetIDs), withinCladeCutOff=float(wthCladeCutOff), name='20180120_forcioli_standalone')
    # newAnalysisObj.description = '20180120_forcioli_standalone'
    #
    # newAnalysisObj.save()
    # dataSubCollectionRun_from_local_production.main(newAnalysisObj, cores)

    # DB version 1.0 = analysisTwo id 121
    # DB version 1.3 = data_analysis id 174
    # DB version 1.4 = data_analysis id 175
    # DB version 1.5 = data_analysis id 176
    #Tester to see that local running works id 176+
    # dataSubCollectionRun_from_local_production.main(data_analysis.objects.get(id=175), cores=cores)
    # dataSubCollectionRun_from_local_production.test()
    # dataSubCollectionRun_from_local_production.compareIntraAbundOfTypes([5297,5326],130)
    ######################################################

    # dataSubCollectionRun_from_local_production.formatOutput(data_analysis.objects.get(id=176), numProcessors=10,
    #                                  datasubstooutput='57')
    ########## Query a sequence against the SymPortal DBVersion. ###########
    # What do we want this to acheive.
    # Let's input a fasta file of sequences
    # Let's make a blast dictionary of the reference sequences and blast our sequences against that.
    # BUT for the time being let's just look to see if there is an exact match for the sequences in the Aiptasia





    ########### Run against database version ############
    # This will run a given data_set object against a database version. A database version is simply a previous
    # data_analysis. As such the input will be an ID of a dataAnalysis and the ID of the data_set to be run.

    # The test run will be running Robertos data against a database version that contains, Voolstra, Emily, Ryan and Ed data.
    # databaseVersionID = 344
    # dSID = 115
    # #
    # nameForAnalysis = 'ed_Against_DBV1.11'
    # descriptionForAnalysis = ' This is an against database test to get the logic sorted out for the MS.'
    # newAnalysisObj = data_analysis(dbVersionAnalysis = False, dbVersionToRunAgainstID=databaseVersionID, dSIDToAnalyse=dSID, name=nameForAnalysis, description=descriptionForAnalysis)
    # newAnalysisObj.save()
    # againstDBVersionAnalysis_060318_production.main(newAnalysisObj, cores)

    # againstDBVersionAnalysis.main(data_analysis.objects.get(id=213), cores)

    #####################################################
def vacuum_db():
    from django.db import connection
    cursor = connection.cursor()
    cursor.execute("VACUUM")
    connection.close()


if __name__ == "__main__":
    fn = os.path.join(os.path.dirname(__file__), 'my_file')

    main()

    #vacuum_db()
