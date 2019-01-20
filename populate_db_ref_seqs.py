#!/usr/bin/env python3.6

# Django specific settings
import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.conf import settings
# ####### Setup Django DB and Models ########
# Ensure settings are read
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()

# Your application specific imports
from dbApp.models import ReferenceSequence
# ###########################################


def read_defined_file_to_list(filename):
    with open(filename, mode='r') as reader:
        temp_list = [line.rstrip() for line in reader]
    return temp_list


def write_list_to_destination(destination, list_to_write):
    try:
        os.makedirs(os.path.dirname(destination))
    except FileExistsError:
        pass

    with open(destination, mode='w') as writer:
        i = 0
        while i < len(list_to_write):
            if i != len(list_to_write)-1:
                writer.write(list_to_write[i] + '\n')
            elif i == len(list_to_write)-1:
                writer.write(list_to_write[i])
            i += 1


def populate_db_with_ref_seqs():
    fasta_to_populate_from = read_defined_file_to_list(
        os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB', 'refSeqDB.fa')))
    for i in range(len(fasta_to_populate_from)):
        if fasta_to_populate_from[i][0] == '>':
            try:
                print('Sequence {} already in db'.format(fasta_to_populate_from[i][1:]))
            except:
                new_seq = ReferenceSequence(
                    name=fasta_to_populate_from[i][1:], clade=fasta_to_populate_from[i][1],
                    sequence=fasta_to_populate_from[i+1], hasName=True)
                new_seq.save()
                print('Sequence {} added to db'.format(fasta_to_populate_from[i][1:]))
    return


populate_db_with_ref_seqs()
