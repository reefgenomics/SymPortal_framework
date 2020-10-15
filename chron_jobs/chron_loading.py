#!/usr/bin/env python3
"""This script will be run as a chron job
It will look for submissions that have a status of transfer_to_framework_server_complete.
The files and datasheets associated with the submission will be run loaded into the SymPortal database.
As part of this process, a dataset and study object will be associated to the Submission object.
Once the loading is complete the Submission will have a status of database_loading_complete.
If the Submission object has for_analysis as True it will be picked up by a chron job and submitted as part
of an analysis. Else the results will be picked up by a different chron job and transferred back over to the web
server. Whether or not result file will be output as part of the loadings will depend on the for_analysis of the
Submission object. If results are being output as part of the loading, the output directory will be stored in the
framework_results_dir_path attribute of the Submission object.

As with the other chron jobs we will ensure that only one instance of the script is running at any given time.
"""

import subprocess
import platform
import os
import sys
sys.path.append("..")
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
from dbApp.models import Submission
import main

class ChronLoading:
    def __init__(self):
        self._check_no_other_instance_running()

        self.submissions_to_transfer = list(
            Submission.objects.filter(progress_status="transfer_to_framework_server_complete", error_has_occured=False)
        )

        # Dynamics
        self.submission_to_load = None
        self.work_flow_manager = None

    def load(self):
        """
        Perform a SymPortal loading for each submission in self.submissions_to_transfer.
        Await a successful exit code before starting the next loading.
        """
        # Load each submission
        for sub_to_load in self.submissions_to_transfer:
            self.submission_to_load = sub_to_load
            self._load_submission()

    def _load_submission(self):
        """
        Do a loading by creating a custom args list that can be processed by a SymPortalWorkFlowManager
        """
        # Get the name of the datasheet
        datasheet_path = [obj for obj in os.listdir(self.submission_to_load.framework_local_dir_path) if obj.endswith('.xlsx')]
        # TODO implement error logic
        try:
            assert(len(datasheet_path) == 1)
        except AssertionError as e:
            # TODO make the name of the datasheet standardised
            raise NotImplementedError('More than one .xlsx found in the dataset')
        datasheet_path = os.path.join(self.submission_to_load.framework_local_dir_path, datasheet_path[0])
        # number of proc will be the minimum between 30 and the number of samples in the submission
        # However when running this on the mac for debug we will pull this down to 4
        if platform.system() == 'Linux':
            num_proc = min(30, int(self.submission_to_load.number_samples))
        else:
            num_proc = 4
        # We will only output results if this loading is not proceeding on to an analysis
        # TODO work out how to assign a Study object to the Submission without a command line input
        if self.submission_to_load.for_analysis:
            # No outputs
            custom_args_list = [
                '--load', self.submission_to_load.framework_local_dir_path,
                '--data_sheet', datasheet_path, '--num_proc', str(num_proc), '--no_output',
                '--name', self.submission_to_load.name
            ]
        else:
            # With output
            custom_args_list = [
                '--load', self.submission_to_load.framework_local_dir_path,
                '--data_sheet', datasheet_path, '--num_proc', str(num_proc),
                '--name', self.submission_to_load.name
            ]
        try:
            self.work_flow_manager = main.SymPortalWorkFlowManager(custom_args_list)
            self.work_flow_manager.start_work_flow()
        except Exception as e:
            raise NotImplementedError('An error has occured while trying to load the Submission data.')

        # Once here we will have finished the loading.
        # Update the status of the Submission object and assign results path if not for_analysis


    def _check_no_other_instance_running(self):
        if sys.argv[1] == 'debug':  # For development only
            pass
        else:
            captured_output = subprocess.run(['pgrep', '-f', 'chron_loading.py'], capture_output=True)
            if captured_output.returncode == 0:  # PIDs were returned
                procs = captured_output.stdout.decode('UTF-8').rstrip().split('\n')
                if platform.system() == 'Linux':
                    # Then we expect there to be one PID for the current process
                    if len(procs) > 1:
                        sys.exit()
                else:
                    # Then we are likely on mac and we expect no PIDs
                    sys.exit()
            else:
                # No PIDs returned
                pass