#!/usr/bin/env python3.6

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

import data_sub_collection_run
import create_data_submission
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

    ###### args ######

    parser = argparse.ArgumentParser(description='Intragenomic analysis of the ITS2 region of the nrDNA', epilog='For support email: symportal@gmail.com')

    group = parser.add_mutually_exclusive_group(required=True)

    # Mutually exclusive arguments
    group.add_argument('--submit', metavar='path_to_dir', help='Run this to submit data to the framework\'s database. '
                                                               'This command takes a '
                                                               'single argument giving the directory containing the '
                                                               'paired sequencing reads in .fastq.gz '
                                                               'format. Alternatively, the path can point directly to a'
                                                               ' single compressed file containing the '
                                                               'same paired fastq.gz files. A name must be associated with '
                                                               'the data_set using the --name flag. The number of '
                                                               'processes to use can also be specified using the '
                                                               '--num_proc flag')


    group.add_argument('--display_data_sets', action='store_true', help='Display data_sets currently in the framework\'s database')


    group.add_argument('--analyse', metavar='data_set IDs', help='Analyse one or more data_set objects '
                                                                       'together. Enter comma separated '
                                                                       'IDs of the data_set IDs you which to '
                                                                       'analyse. e.g.: 43,44,45. '
                                                           'If you wish to use all available dataSubmissions,'
                                                                       'you may pass \'all\' as an argument. '
                                                           'To display all data_sets currently submitted to the '
                                                           'framework\'s database, including their ids, use the \'show_data_sets\' command' )


    group.add_argument('--print_output', metavar='data_set IDs, analysis ID', help='Use this function to output the '
                                                                                       'ITS2 sequence and ITS2 type profile count tables for '
                                                                                       'a given set of data_sets '
                                                                                       'that have been run in a given analysis. Give the data_set IDs '
                                                                                       'that you wish to make outputs for as '
                                                                                       'arguments to the --print_output flag. To output for multiple data_set objects, '
                                                                                       'comma separate the IDs of the data_set objects, '
                                                                                       'e.g. 44,45,46.'
                                                                                       'Give the ID of the analysis you wish to '
                                                                                       'output these from using the --data_analysis_id '
                                                                                       'flag.')


    group.add_argument('--display_analyses', action='store_true', help=' Display data_analysis objects currently '
                                                                       'stored in the framework\'s database')

    # Additional arguments
    parser.add_argument('--num_proc', type=int, help='Number of processors to use', default=1)
    parser.add_argument('--name', help='A name for your input or analysis', default='noName')
    parser.add_argument('--description', help='An optional description', default='No description')
    parser.add_argument('--data_analysis_id', type=int, help='The ID of the data_analysis you wish to output from')
    group.add_argument('--vacuum_database', action='store_true', help='Vacuuming the database will free up memory from '
                                                                     'objects that have been deleted recently')
    args = parser.parse_args()


    # Code to run the main functionality of SymPortal
    if args.submit:
        if args.name == 'noName':
            print('Please provide a name using the --name flag. e.g. --name splendid_dataset')
            return


        # This is how we read data into the SP database
        # This should be given a string on the command line which is a directory that contains the sequencing data.
        # We will need to create a data_set object for each time this command is run.
        num_proc = args.num_proc

        # input directory should contain either paired .fastq.gz files of a single .zip file
        input_dir = args.submit

        name_for_data_set = args.name

        # If working on the remote server a difference reference_fasta_database_used can be used.
        new_data_set = data_set(name = name_for_data_set, timeStamp=str(datetime.now()), reference_fasta_database_used='symClade.fa')
        new_data_set.save()

        # only perform sub_evalue_screening when working on the remote system
        with open('sp_config') as f:
            config_dict = json.load(f)
        if config_dict['system_type'] == 'remote':
            screen_sub_evalue_bool = True
        else:
            screen_sub_evalue_bool = False

        create_data_submission.main(input_dir, new_data_set.id, num_proc, screen_sub_evalue=screen_sub_evalue_bool)


    elif args.analyse:
        if args.name == 'noName':
            print('Please provide a name using the --name flag. e.g. --name wonderful_analysis')
            return

        within_clade_cutoff = 0.03

        num_proc = args.num_proc
        custom_data_set_ids = args.analyse
        if args.analyse == 'all':
            tempList = []
            for ds in data_set.objects.all():
                tempList.append(str(ds.id))
            stringList = ','.join(tempList)
            custom_data_set_ids = stringList
        new_analysis_object = data_analysis(listOfDataSubmissions=str(custom_data_set_ids),
                                         withinCladeCutOff=float(within_clade_cutoff), name=args.name,
                                            timeStamp=str(datetime.now()))
        new_analysis_object.description = args.description
        new_analysis_object.save()
        data_sub_collection_run.main(dataanalysistwoobject=new_analysis_object, cores=num_proc)
        print('return code: 0\nAnalysis complete')


    elif args.print_output:
        if args.data_analysis_id:
            data_sub_collection_run.formatOutput_ord(data_analysis.objects.get(id=args.data_analysis_id), numProcessors=args.num_proc, datasubstooutput=args.print_output)
        else:
            print('Please provide a data_analysis to ouput from by providing a data_analysis ID to the --data_analysis_id '
                  'argument. To see a list of data_analysis objects in the framework\'s database, use the --display_analyses argument.')


    elif args.display_data_sets:
        # Then print out all of the dataSubmissions with names and IDs in the db
        for ds in data_set.objects.all():
            print('{}: {}\t{}'.format(ds.id, ds.name, ds.timeStamp))


    elif args.display_analyses:
        # Then print out all of the dataAnalysisTwos with names and IDs in the db
        for da in data_analysis.objects.all():
            print('{}: {}\t{}'.format(da.id, da.name, da.timeStamp))


    elif args.vacuumDatabase:
        print('Vacuuming database')
        vacuum_db()
        print('Vacuuming complete')



def vacuum_db():
    from django.db import connection
    cursor = connection.cursor()
    cursor.execute("VACUUM")
    connection.close()


if __name__ == "__main__":
    fn = os.path.join(os.path.dirname(__file__), 'my_file')

    main()

