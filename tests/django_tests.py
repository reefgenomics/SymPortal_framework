#!/usr/bin/env python3.6
from django.test import TestCase
import os
import main
import general
from pathlib import Path
import shutil
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from dbApp.models import DataAnalysis, DataSet, ReferenceSequence

class SPIntegrativeTesting(TestCase):
    """This fixture is a dump of a database with contents specifically setup for testing. It contains
    three loaded DataSet objects with all of the corresponding features and a single analysis that contains
    all three of the DataSets.

    The DataSet IDs are:
    XXX
    XXX
    XXX

    The DataAnalysis ID is:
    XXX
    """
    fixtures = ['testing_fixture.json']

    def test_sequence_count_table_creator(self):
        for ds in DataSet.objects.all():
            print(ds.id)
        print('done')
        # custom_args_list = ['--print_output_seqs', self.test_data_dir_path, '--num_proc', str(self.num_proc)]
        # test_spwfm = main.SymPortalWorkFlowManager(custom_args_list)
        # self.completed_data_loading_object = test_spwfm.start_work_flow()

