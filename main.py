#!/usr/bin/env python3.6
''' SymPortal: a novel analytical framework and platform for coral algal
    symbiont next-generation sequencing ITS2 profiling
    Copyright (C) 2018  Benjamin C C Hume

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see
    https://github.com/SymPortal/SymPortal_framework/tree/master/LICENSE.txt.'''
import argparse
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
import output
import plotting
import json
import sys
import distance

def main():

    ###### args ######

    parser = argparse.ArgumentParser(description='Intragenomic analysis of the ITS2 region of the nrDNA', epilog='For support email: symportal@gmail.com')

    group = parser.add_mutually_exclusive_group(required=True)

    # Mutually exclusive arguments
    group.add_argument('--submit', metavar='path_to_dir', help='Run this to submit data to the framework\'s database. '
                                                               'The first argument to this command must be an absolute '
                                                               'path to a directory containing  the '
                                                               'paired sequencing reads in .fastq.gz '
                                                               'format. Alternatively, this path can point directly to a'
                                                               ' single compressed file containing the '
                                                               'same paired fastq.gz files. '
                                                               '\nA name must be associated with '
                                                               'the data_set using the --name flag. '
                                                               '\nThe number of processes to use can also be specified '
                                                               'using the --num_proc flag.'
                                                               '\nA datasheet can also be uploaded using the '
                                                               '--data_sheet flag and the full path to the .xlxs '
                                                               'data_sheet file (RECOMMENDED). \n'
                                                               'To skip the generation of figures pass the '
                                                               '--noFig flag.\n'
                                                               'To skip the generation of ordination files '
                                                               '(pairwise distances and PCoA coordinates) '
                                                               'pass the --noOrd flag')

    group.add_argument('--analyse', metavar='data_set IDs', help='Analyse one or more data_set objects '
                                                                 'together. Enter comma separated '
                                                                 'IDs of the data_set IDs you which to '
                                                                 'analyse. e.g.: 43,44,45. '
                                                                 'If you wish to use all available dataSubmissions,'
                                                                 'you may pass \'all\' as an argument. '
                                                                 'To display all data_sets currently submitted to the '
                                                                 'framework\'s database, including their ids, use the \'show_data_sets\' command\n'
                                                                 'To skip the generation of figures pass the '
                                                                 '--noFig flag.\n'
                                                                 'To skip the generation of ordination files '
                                                                 '(pairwise distances and PCoA coordinates) '
                                                                 'pass the --noOrd flag')


    group.add_argument('--display_data_sets', action='store_true', help='Display data_sets currently in the framework\'s database')

    group.add_argument('--display_analyses', action='store_true', help=' Display data_analysis objects currently '
                                                                       'stored in the framework\'s database')

    group.add_argument('--print_output_seqs', metavar='data_set IDs',
                       help='Use this function to output ITS2 sequence count tables for given data_set instances')

    group.add_argument('--print_output_types', metavar='data_set IDs, analysis ID', help='Use this function to output the '
                                                                                       'ITS2 sequence and ITS2 type profile count tables for '
                                                                                       'a given set of data_sets '
                                                                                       'that have been run in a given analysis. Give the data_set IDs '
                                                                                       'that you wish to make outputs for as '
                                                                                       'arguments to the --print_output_types flag. To output for multiple data_set objects, '
                                                                                       'comma separate the IDs of the data_set objects, '
                                                                                       'e.g. 44,45,46.'
                                                                                       'Give the ID of the analysis you wish to '
                                                                                       'output these from using the --data_analysis_id '
                                                                                       'flag.\nTo skip the generation of figures pass the '
                                                                                       '--noFig flag.')



    group.add_argument('--between_type_distances',
                       metavar='data_set IDs, analysis ID',
                       help='Use this function to output UniFrac pairwise distances '
                            'between ITS2 type profiles clade separated')

    group.add_argument('--between_sample_distances',
                       metavar='data_set IDs',
                       help='Use this function to output UniFrac pairwise distances '
                            'between samples clade separated from a given collection of data_set objects')
    group.add_argument('--between_sample_distances_sample_set',
                       metavar='data_set_sample IDs',
                       help='Use this function to output UniFrac pairwise distances '
                            'between samples clade separated from a given collection of data_set_sample objects')


    # Additional arguments
    parser.add_argument('--num_proc', type=int, help='Number of processors to use', default=1)
    parser.add_argument('--name', help='A name for your input or analysis', default='noName')
    parser.add_argument('--description', help='An optional description', default='No description')
    parser.add_argument('--data_analysis_id', type=int, help='The ID of the data_analysis you wish to output from')
    group.add_argument('--vacuum_database', action='store_true', help='Vacuuming the database will free up memory from '
                                                                     'objects that have been deleted recently')
    parser.add_argument('--bootstrap', type=int, help='Number of bootstrap iterations to perform', default=100)
    parser.add_argument('--data_sheet', help='An absolute path to the .xlxs file containing the meta-data information for the data_set\'s samples')
    parser.add_argument('--noFig', action='store_true', help='Skip figure production')
    parser.add_argument('--noOrd', action='store_true', help='Skip ordination analysis')
    parser.add_argument('--debug', action='store_true', help='Present additional stdout output', default=False)
    parser.add_argument('--noOutput', action='store_true', help='Do no output: count tables, figures, ordinations', default=False)
    parser.add_argument('--distance_method', help='Either \'unifrac\' or \'braycurtis\', default=braycurtis. '
                                                  'The method to use when calculating distances between its2 '
                                                  'type profiles or samples.', default='braycurtis')
    # when run as remote
    parser.add_argument('--submitting_user_name', help='Only for use when running as remote'
                                                       '\nallows the association of a different user_name '
                                                       'to the data_set than the one listed in sp_config',
                        default='not supplied')
    parser.add_argument('--submitting_user_email', help='Only for use when running as remote'
                                                        '\nallows the association of a different user_email '
                                                        'to the data_set than the one listed in sp_config',
                        default='not supplied')
    args = parser.parse_args()


    # Code to run the main functionality of SymPortal
    if args.submit:
        if args.name == 'noName':
            sys.exit('Please provide a name using the --name flag. e.g. --name splendid_dataset')



        # This is how we read data into the SP database
        # This should be given a string on the command line which is a directory that contains the sequencing data.
        # We will need to create a data_set object for each time this command is run.
        num_proc = args.num_proc

        # input directory should contain either paired .fastq.gz files of a single .zip file
        input_dir = args.submit

        name_for_data_set = args.name


        with open('{}/sp_config'.format(os.path.dirname(__file__))) as f:
            config_dict = json.load(f)
        local_or_remote = config_dict['system_type']
        if local_or_remote == 'remote':
            screen_sub_evalue_bool = True
            if args.submitting_user_name and args.submitting_user_email:
                new_data_set_submitting_user = args.submitting_user_name
                new_data_set_user_email = args.submitting_user_email
            else:
                print('Please supply --submitting_user_name and --submitting_user_email.')
                sys.exit(1)
        else:
            new_data_set_submitting_user = config_dict['user_name']
            new_data_set_user_email = config_dict['user_email']
            screen_sub_evalue_bool = False


        # If working on the remote server a difference reference_fasta_database_used can be used.
        new_data_set = data_set(name = name_for_data_set, timeStamp=str(datetime.now()).replace(' ', '_').replace(':', '-'),
                                reference_fasta_database_used='symClade.fa',
                                submittingUser=new_data_set_submitting_user,
                                submitting_user_email=new_data_set_user_email)
        new_data_set.save()

        # only perform sub_evalue_screening when working on the remote system



        if args.data_sheet:
            if os.path.isfile(args.data_sheet):
                create_data_submission.main(input_dir, new_data_set.id, num_proc,
                                            screen_sub_evalue=screen_sub_evalue_bool,
                                            data_sheet_path=args.data_sheet, noFig=args.noFig, noOrd=args.noOrd,
                                            distance_method=args.distance_method, debug=args.debug)
            else:
                sys.exit('{} not found'.format(args.data_sheet))
        else:
            create_data_submission.main(input_dir, new_data_set.id, num_proc,
                                        screen_sub_evalue=screen_sub_evalue_bool, noFig=args.noFig, noOrd=args.noOrd,
                                        distance_method=args.distance_method, debug=args.debug)

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

        with open('{}/sp_config'.format(os.path.dirname(__file__))) as f:
            config_dict = json.load(f)
        new_data_set_submitting_user = config_dict['user_name']
        new_data_set_user_email = config_dict['user_email']

        new_analysis_object = data_analysis(listOfDataSubmissions=str(custom_data_set_ids),
                                         withinCladeCutOff=float(within_clade_cutoff), name=args.name,
                                            timeStamp=str(datetime.now()).replace(' ', '_').replace(':', '-'), submittingUser=new_data_set_submitting_user,
                                            submitting_user_email=new_data_set_user_email)
        new_analysis_object.description = args.description
        new_analysis_object.save()
        data_sub_collection_run.main(dataanalysistwoobject=new_analysis_object, cores=num_proc, noFig=args.noFig,
                                     noOrd=args.noOrd, distance_method=args.distance_method, noOutput=args.noOutput,
                                     debug=args.debug)
        print('return code: 0\nAnalysis complete')

    elif args.print_output_types:
        with open('{}/sp_config'.format(os.path.dirname(__file__))) as f:
            config_dict = json.load(f)
        new_data_set_submitting_user = config_dict['user_name']
        if args.data_analysis_id:
            analysis_object = data_analysis.objects.get(id=args.data_analysis_id)
            dataSubmissionsToOutput = [int(a) for a in args.print_output_types.split(',')]
            querySetOfDataSubmissions = data_set.objects.filter(id__in=dataSubmissionsToOutput)
            num_samples = len(data_set_sample.objects.filter(dataSubmissionFrom__in=querySetOfDataSubmissions))

            data_sub_collection_run.formatOutput_ord(analysisobj=analysis_object,
                                                     numProcessors=args.num_proc, call_type='stand_alone', num_samples=num_samples,
                                                     datasubstooutput=args.print_output_types, noFig=args.noFig,
                                                     output_user=new_data_set_submitting_user)
        else:
            print('Please provide a data_analysis to ouput from by providing a data_analysis ID to the --data_analysis_id '
                  'flag. To see a list of data_analysis objects in the framework\'s database, use the --display_analyses flag.')

    elif args.display_data_sets:
        # Then print out all of the dataSubmissions with names and IDs in the db
        # Let's sort this by id when outputting
        # When just printing the raw query the order is somewhat random
        data_set_id_to_obj_dict = {ds.id: ds for ds in list(data_set.objects.all())}
        sorted_list_of_ids = sorted(list(data_set_id_to_obj_dict.keys()))

        for ds_id in sorted_list_of_ids:
            ds_in_q = data_set_id_to_obj_dict[ds_id]
            print('{}: {}\t{}'.format(ds_in_q.id, ds_in_q.name, ds_in_q.timeStamp))

    elif args.display_analyses:
        # Then print out all of the dataAnalysisTwos with names and IDs in the db
        # Let's sort this by id when outputting
        # When just printing the raw query the order is somewhat random
        data_analysis_id_to_obj_dict = {da.id: da for da in list(data_analysis.objects.all())}
        sorted_list_of_ids = sorted(list(data_analysis_id_to_obj_dict.keys()))

        for da_id in sorted_list_of_ids:
            da_in_q = data_analysis_id_to_obj_dict[da_id]
            print('{}: {}\t{}'.format(da_in_q.id, da_in_q.name, da_in_q.timeStamp))

    elif args.between_type_distances:
        dts = str(datetime.now()).replace(' ', '_').replace(':', '-')
        if args.data_analysis_id:
            if args.distance_method == 'unifrac':
                pcoa_path_list = data_sub_collection_run.generate_within_clade_UniFrac_distances_ITS2_type_profiles(
                    data_submission_id_str=args.between_type_distances, num_processors=args.num_proc,
                    data_analysis_id=args.data_analysis_id, method='mothur', call_type = 'stand_alone', date_time_string=dts, bootstrap_value=args.bootstrap)
            elif args.distance_method == 'braycurtis':
                pcoa_path_list = distance.generate_within_clade_BrayCurtis_distances_ITS2_type_profiles(
                    data_submission_id_str=args.between_type_distances,
                    data_analysis_id=args.data_analysis_id, call_type = 'stand_alone', date_time_string=dts)
        else:
            print('Please provide a data_analysis to ouput from by providing a data_analysis ID to the --data_analysis_id '
                  'argument. To see a list of data_analysis objects in the framework\'s database, use the --display_analyses flag.')

        for pcoa_path in pcoa_path_list:
            if 'PCoA_coords' in pcoa_path:
                sys.stdout.write('\nPlotting between its2 type profile distances clade {}\n'.format(
                    os.path.dirname(pcoa_path).split('/')[-1]))
                # then this is a pcoa csv that we should plot
                plotting.plot_between_its2_type_prof_dist_scatter(pcoa_path, date_time_str=dts)



    elif args.between_sample_distances:
        # we are swaping to bray curtis for the time being
        dts = str(datetime.now()).replace(' ', '_').replace(':', '-')
        if args.distance_method == 'unifrac':
            PCoA_paths_list = distance.generate_within_clade_UniFrac_distances_samples(
                dataSubmission_str=args.between_sample_distances, num_processors=args.num_proc,
                method='mothur', call_type='stand_alone', date_time_string=dts, bootstrap_value=args.bootstrap)
        elif args.distance_method == 'braycurtis':
            PCoA_paths_list = distance.generate_within_clade_BrayCurtis_distances_samples(dataSubmission_str=args.between_sample_distances, call_type='stand_alone', date_time_string=dts)

        for pcoa_path in PCoA_paths_list:
            if 'PCoA_coords' in pcoa_path:
                # then this is a full path to one of the .csv files that contains the coordinates that we can plot
                # we will get the output directory from the passed in pcoa_path
                sys.stdout.write('\n\nGenerating between sample distance plot clade {}\n'.format(
                    os.path.dirname(pcoa_path).split('/')[-1]))
                plotting.plot_between_sample_distance_scatter(pcoa_path, dts)

    elif args.between_sample_distances_sample_set:
        # this is a variation of the between_sample_distances where a set of sample IDs are input rather
        # than data_set IDs. Currently it is only written into unifrac. Once we have this running
        # we can then apply it to BrayCurtis
        dts = str(datetime.now()).replace(' ', '_').replace(':', '-')
        if args.distance_method == 'unifrac':
            PCoA_paths_list = distance.generate_within_clade_UniFrac_distances_samples_sample_list_input(
                smpl_id_list_str=args.between_sample_distances_sample_set, num_processors=args.num_proc,
                method='mothur', bootstrap_value=args.bootstrap,
                date_time_string=dts)
        elif args.distance_method == 'braycurtis':
            PCoA_paths_list = distance.generate_within_clade_BrayCurtis_distances_samples_sample_list_input(
                smpl_id_list_str=args.between_sample_distances_sample_set, date_time_string=dts)
        for pcoa_path in PCoA_paths_list:
            if 'PCoA_coords' in pcoa_path:
                # then this is a full path to one of the .csv files that contains the coordinates that we can plot
                # we will get the output directory from the passed in pcoa_path
                sys.stdout.write('\n\nGenerating between sample distance plot clade {}\n'.format(
                    os.path.dirname(pcoa_path).split('/')[-1]))
                plotting.plot_between_sample_distance_scatter(pcoa_path, dts)

    elif args.print_output_seqs:
        # this is a stand_alone and output and we should grab the user who is requresting it from the config file
        with open('{}/sp_config'.format(os.path.dirname(__file__))) as f:
            config_dict = json.load(f)
        new_data_set_submitting_user = config_dict['user_name']

        outputDir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'outputs/non_analysis'))
        output_file_list, date_time_str, num_samples = output.div_output_pre_analysis_new_meta_and_new_dss_structure(
            datasubstooutput=args.print_output_seqs, numProcessors=args.num_proc, output_dir=outputDir,
            call_type='stand_alone', output_user=new_data_set_submitting_user)
        if num_samples > 1000:
            print('Too many samples ({}) to generate plots'.format(num_samples))
        else:
            for item in output_file_list:
                if 'relative' in item:
                    svg_path, png_path = plotting.\
                        generate_stacked_bar_data_submission(path_to_tab_delim_count=item,
                                                           output_directory=outputDir, time_date_str=date_time_str)
                    print('Output figs:')
                    print(svg_path)
                    print(png_path)
                    break



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

