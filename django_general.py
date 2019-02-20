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

