#!/usr/bin/env python3

# Django specific settings
import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.conf import settings
from general import ThreadSafeGeneral
# ####### Setup Django DB and Models ########
# Ensure settings are read
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
from django.core.exceptions import ObjectDoesNotExist
# Your application specific imports
from dbApp.models import ReferenceSequence
# ###########################################

def populate_db_with_ref_seqs():
    thread_safe_general = ThreadSafeGeneral()
    fasta_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniaceaeDB', 'refSeqDB.fa'))
    fasta_dict = thread_safe_general.create_dict_from_fasta(fasta_path=fasta_path)
    current_ref_seq_names = [rs.name for rs in ReferenceSequence.objects.all()]
    bulk_new_rs_list = []
    created_name_list = []
    for new_name, new_seq in fasta_dict.items():
        if new_name not in current_ref_seq_names:
            bulk_new_rs_list.append(
                ReferenceSequence(name=new_name, clade=new_name[0], sequence=new_seq, has_name=True))
            created_name_list.append(new_name)
    ReferenceSequence.objects.bulk_create(bulk_new_rs_list)
    if created_name_list:
        for name in created_name_list:
            print(f'Sequence {name} added to db')
    else:
        print('All sequences already present in database.')


populate_db_with_ref_seqs()
