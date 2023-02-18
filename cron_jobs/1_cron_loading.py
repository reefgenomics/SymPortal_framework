#!/usr/bin/env python3
"""This script will be run as a cron job
It will look for submissions that have a status of transfer_to_framework_server_complete.
The files and datasheets associated with the submission will be run loaded into the SymPortal database.
As part of this process, a dataset and study object will be associated to the Submission object.
Once the loading is complete the Submission will have a status of database_loading_complete.
If the Submission object has for_analysis as True it will be picked up by a cron job and submitted as part
of an analysis. Else the results will be picked up by a different cron job and transferred back over to the web
server. Whether or not result file will be output as part of the loadings will depend on the for_analysis of the
Submission object. If results are being output as part of the loading, the output directory will be stored in the
framework_results_dir_path attribute of the Submission object.

As with the other cron jobs we will ensure that only one instance of the script is running at any given time.
"""

import subprocess
import platform
import os
import sys
from pathlib import Path
# We have to add the Symportal_framework path so that the settings.py module
# can be found.
os.chdir(str(Path(__file__).resolve().parent.parent))
sys.path.append(str(Path(__file__).resolve().parent.parent))
sys.path.append(str(Path(__file__).resolve().parent))

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
from dbApp.models import Submission
from django.db.models import Q
import main
from datetime import datetime


class CronLoading:
    def __init__(self):
        self._check_no_other_instance_running()

        self.submissions_to_load = list(
            Submission.objects.filter(progress_status="transfer_to_framework_server_complete", error_has_occured=False, loading_started_date_time=None)
        )

        # Also check for the presence of a submission that has started loading but has not completed
        time_now = str(datetime.utcnow()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')
        submissions_in_progress = Submission.objects.filter(progress_status="transfer_to_framework_server_complete").filter(~Q(loading_started_date_time=None))
        if submissions_in_progress:
            print(f"{time_now}: Incomplete loading detected:")
            for sub in submissions_in_progress:
                print(f"\t{sub.id}: {sub.name}")
            sys.exit("Incomplete loading detected:")

        
        if self.submissions_to_load:
            print(f"{time_now}: The following submission have been found to load.")
            for sub in self.submissions_to_load:
                print(sub.name)
        else:
            print(f"{time_now}: No submission found for loading.")

        # Dynamics
        self.submission_to_load = None
        self.work_flow_manager = None

    def load(self):
        """
        Perform a SymPortal loading for each submission in self.submissions_to_transfer.
        Await a successful exit code before starting the next loading.
        """
        # Load each submission
        for sub_to_load in self.submissions_to_load:
            self.submission_to_load = sub_to_load
            self._load_submission()

    def _load_submission(self):
        """
        Do a loading by creating a custom args list that can be processed by a SymPortalWorkFlowManager
        """
        # Get the name of the datasheet
        datasheet_path = os.path.join(self.submission_to_load.framework_local_dir_path, f"{self.submission_to_load.name}_datasheet.xlsx")
        if not os.path.exists(datasheet_path):
            datasheet_path = os.path.join(self.submission_to_load.framework_local_dir_path, f"{self.submission_to_load.name}_datasheet.csv")
            if not os.path.exists(datasheet_path):
                raise FileNotFoundError(f"Couldn't find {datasheet_path}")
        # TODO implement error logic
        # number of proc will be the minimum between 30 and the number of samples in the submission
        # However when running this on the mac for debug we will pull this down to 4
        if platform.system() == 'Linux':
            num_proc = min(30, int(self.submission_to_load.number_samples))
        else:
            num_proc = 4
        # We will only output results if this loading is not proceeding on to an analysis
        # We will provide the --study_user_string and --study_name arguments for use in creating the Study object
        # We will also pass --is_cron_loading to let SP know that this is being initiated by a cron job
        # Currently we only have a single user associated to each Submission object so we will only be able
        # to associate a single user to the Study object
        study_user_string = self.submission_to_load.submitting_user.name
        if self.submission_to_load.for_analysis:
            # No outputs
            custom_args_list = [
                '--load', self.submission_to_load.framework_local_dir_path,
                '--data_sheet', datasheet_path, '--num_proc', str(num_proc), '--no_output',
                '--name', self.submission_to_load.name, '--is_cron_loading',
                '--study_user_string', study_user_string,
                '--study_name', self.submission_to_load.name
            ]
        else:
            # With output
            custom_args_list = [
                '--load', self.submission_to_load.framework_local_dir_path,
                '--data_sheet', datasheet_path, '--num_proc', str(num_proc),
                '--name', self.submission_to_load.name, '--is_cron_loading',
                '--study_user_string', study_user_string,
                '--study_name', self.submission_to_load.name
            ]
        try:
            self.work_flow_manager = main.SymPortalWorkFlowManager(custom_args_list)
            self.submission_to_load.loading_started_date_time = self.work_flow_manager.date_time_str
            self.submission_to_load.save()
            self.work_flow_manager.start_work_flow()
        except Exception as e:
            # TODO handle errors for Submission objects and cron jobs
            self.submission_to_load.error_has_occured = True
            self.submission_to_load.save()
            raise NotImplementedError('An error has occured while trying to load the Submission data.')

        # Once here we will have finished the loading.
        # Update the status of the Submission object and assign results path if not for_analysis
        if not self.submission_to_load.for_analysis:
            self.submission_to_load.framework_results_dir_path = self.work_flow_manager.data_loading_object.output_directory

        # Assign the associated DataSet and Study objects
        self.submission_to_load.associated_dataset = self.work_flow_manager.data_loading_object.dataset_object
        self.submission_to_load.associated_study = self.work_flow_manager.data_loading_object.study

        # Log the loading complete time
        self.submission_to_load.loading_complete_date_time = str(
                datetime.utcnow()
            ).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')
        if self.submission_to_load.for_analysis:
            self.submission_to_load.progress_status = "framework_loading_complete"
        else:
            # Then there will not be an analysis for this submission
            # We have already complted the output
            self.submission_to_load.progress_status = "framework_output_complete"

        # At this point the loading is complete for a single submission object. Now save and do next.
        self.submission_to_load.save()

    def _check_no_other_instance_running(self):
        try:
            if sys.argv[1] == 'debug':  # For development only
                pass
            else:
                raise RuntimeError('Unknown arg at sys.argv[1]')
        except IndexError:
            captured_output = subprocess.run(['pgrep', '-fa', 'cron_loading.py'], capture_output=True)
            if captured_output.returncode == 0:  # PIDs were returned
                procs = captured_output.stdout.decode('UTF-8').rstrip().split('\n')
                if platform.system() == 'Linux':
                    print("Linux system detected")
                    # Then we expect there to be one PID for the current process
                    # And one for the cron job
                    if len(procs) > 2:
                        print("The following procs were returned:")
                        for p in procs:
                            print(p)
                        time_now = str(datetime.utcnow()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')
                        raise RuntimeError(f"\n{time_now}:More than one instance of cron_loading detected. Killing process.")
                else:
                    # Then we are likely on mac and we expect no PIDs
                    sys.exit()
            else:
                # No PIDs returned
                pass

cron_loading = CronLoading()
cron_loading.load()