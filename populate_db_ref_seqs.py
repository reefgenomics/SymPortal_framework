#!/usr/bin/env python3.6

# Django specific settings
import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.conf import settings
######## Setup Django DB and Models ########
# Ensure settings are read
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()

# Your application specific imports
from dbApp.models import symportal_framework, data_set, reference_sequence, data_set_sample_sequence, analysis_type, analysis_group, data_set_sample, data_analysis, clade_collection, clade_collection_type
############################################

###### Generic functions ######
def readDefinedFileToList(filename):
    temp_list = []
    with open(filename, mode='r') as reader:
        temp_list = [line.rstrip() for line in reader]
    return temp_list

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



def populate_db_with_ref_seqs():
    fasta_to_populate_from = readDefinedFileToList(os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB', 'refSeqDB.fa')))
    for i in range(len(fasta_to_populate_from)):
        if fasta_to_populate_from[i][0] == '>':
            try:
                existing_seq = reference_sequence.objects.get(name=fasta_to_populate_from[i][1:])
                print('Sequence {} already in db'.format(fasta_to_populate_from[i][1:]))
            except:
                new_seq = reference_sequence(name=fasta_to_populate_from[i][1:], clade=fasta_to_populate_from[i][1], sequence=fasta_to_populate_from[i+1], hasName=True)
                new_seq.save()
                print('Sequence {} added to db'.format(fasta_to_populate_from[i][1:]))
    return

populate_db_with_ref_seqs()