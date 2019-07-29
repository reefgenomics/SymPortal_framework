"""A method that relies on Django feature. E.g. model classes."""
import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.conf import settings
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
from dbApp.models import DataSet, DataAnalysis


def delete_data_set(uid):
    DataSet.objects.get(id=uid).delete()


def delete_data_analysis(uid):
    DataAnalysis.objects.get(id=uid).delete()

def write_ref_seq_objects_to_fasta(path, list_of_ref_seq_objs, identifier='name'):
    with open(path, 'w') as f:
        for ref_seq_obj in list_of_ref_seq_objs:
            if identifier == 'name':
                f.write(f'>{ref_seq_obj.name}\n')
            elif identifier == 'id':
                f.write(f'>{ref_seq_obj.id}\n')
            f.write(f'{ref_seq_obj.sequence}\n')

