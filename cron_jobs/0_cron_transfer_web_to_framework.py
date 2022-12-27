#!/usr/bin/env python3
"""
This script will be managed by cron jobs and will be run once every hour
It will be responsible for transfering successfuly submitted user files from the web sever to the framework server
It will check for Submission objects that have a progress_status of submitted
At the start of the script it will do a pgrep to check that the script is not currently running to prevent
stacking up of the same script.
The script will use the md5sum that is generated at the time of user upload to verify the integrity of the transfer.
After transfer is complete the status of the Submission objects that have been transfered will be updated to
tranfer_to_framework_server_complete. The transfer_to_framework_server_date_time will be logged.
A seperate cron job will handle loading the transferred submissions
"""
# TODO have a log file where all the cron jobs can log their outputs
import subprocess
import platform
import os
import shutil
import sys
from pathlib import Path
# We have to add the Symportal_framework path so that the settings.py module
# can be found.
# os.chdir(str(Path(__file__).resolve().parent.parent))
sys.path.append(str(Path(__file__).resolve().parent.parent))
sys.path.append(str(Path(__file__).resolve().parent))
print(sys.path)
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
from dbApp.models import Submission
import sp_config
import paramiko
from datetime import datetime
from pathlib import Path

class TransferWebToFramework:
    def __init__(self):
        # The very first thing to do is to perform a pgrep to see if another instance of this script is being run
        # However, the pgrep runs differently on mac vs linux.
        # On mac, the self process is not included so it will return nothing if only the one process is running
        # On linux, the self process is included so it will return one PID for the current process
        # When debugging in an IDE, on mac, nothing; on linux, it will return multiple processes (probably 2)
        self._check_no_other_instance_running()

        # Get a list of the Submissions that need to be transferred
        self.submissions_to_transfer = list(
            Submission.objects.filter(progress_status="submitted", error_has_occured=False)
        )
        self.symportal_data_dir = sp_config.symportal_data_dir

        # User paramiko to set up an sftp that we can use to transfer
        self.ssh_client = paramiko.SSHClient()
        self.ssh_client.load_system_host_keys()
        if sp_config.authentication_type == 'pass':
            self.ssh_client.connect(hostname=sp_config.web_ip, username=sp_config.web_user,
                                    password=sp_config.web_pass)
        elif sp_config.authentication_type == 'key':
            self.ssh_client.connect(
                hostname=sp_config.web_ip, username=sp_config.web_user, key_filename=sp_config.key_file
                )
        else:
            raise RuntimeError('Unknown authentication_type from sp_config.')

        # Open sftp client
        self.sftp_client = self.ssh_client.open_sftp()

        # Number of attempts already undertaken to download the data from the webserver and validate md5sum
        self.attempts = 0

        # Dynamics for convenience that will be updated with each Submission object
        self.submission_to_transfer = None
        self.web_source_dir = None
        self.framework_dest_dir = None

    def _check_no_other_instance_running(self):
        try:
            if sys.argv[1] == 'debug':  # For development only
                pass
            else:
                raise RuntimeError('Unknown arg at sys.argv[1]')
        except IndexError:
            captured_output = subprocess.run(['pgrep', '-fa', 'cron_transfer_web_to_framework.py'], capture_output=True)
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
                        raise RuntimeError('\nMore than one instance of cron_transfer_web_to_framework detected. Killing process.')
                else:
                    # Then we are likely on mac and we expect no PIDs
                    raise RuntimeError('More than one instance of cron_transfer_web_to_framework detected. Killing process.')
            else:
                # No PIDs returned
                pass

    def transfer(self):
        """
        Transfer over the files in the Submissions web_local_dir_path from the web server
        Validate the md5sum
        Delete the files from the web server
        Update the status and time log of each Submission object
        """

        # Transfer each submission
        for sub_to_trans in self.submissions_to_transfer:
            self.submission_to_transfer = sub_to_trans
            self._process_submission()

    def _process_submission(self, ):
        self.web_source_dir = self.submission_to_transfer.web_local_dir_path
        self.framework_dest_dir = os.path.join(self.symportal_data_dir, os.path.basename(self.web_source_dir))
        self.submission_to_transfer.framework_local_dir_path = self.framework_dest_dir
        self.submission_to_transfer.save()
        if not os.path.exists(self.framework_dest_dir):
            os.makedirs(self.framework_dest_dir)
        try:
            self._get_files_from_web_dir()
            md5sum_source_dict = self._make_md5sum_web_server_dict()

            self._validate_md5sum(md5sum_source_dict)

            # If we get here then the md5sum has been validated and the transfer is complete

            # Delete the remote files and the remote directory
            for get_file in self.sftp_client.listdir(self.web_source_dir):
                self.sftp_client.remove(os.path.join(self.web_source_dir, get_file))
            self.sftp_client.rmdir(self.web_source_dir)
            
            # we want to delete the parent directory if it is empty
            if len(list(self.sftp_client.listdir(str(Path(self.web_source_dir).parent.absolute())))) == 0:
                self.sftp_client.rmdir(str(Path(self.web_source_dir).parent.absolute()))

            # Update the time of the Submission object transfer_to_framework_server_date_time
            self.submission_to_transfer.transfer_to_framework_server_date_time = str(
                datetime.utcnow()
            ).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')

            # Update the status of the Submission object
            self.submission_to_transfer.progress_status = "transfer_to_framework_server_complete"

            self.submission_to_transfer.save()
        except FileNotFoundError as e:
            # It happens that for some reason a submission has been created
            # But the files ae not available on the linode server
            # In this case, we will delte the framework_dest_dir we just created
            # and delete the submission object from the database
            # and then move on to processing the next submission.
            print(f"A FileNotFoundError was caught for submission:")
            print(f"\tSubmission name: {self.submission_to_transfer.name}")
            print(f"\tSubmission name: {self.submission_to_transfer.id}")
            print(f"This submission will now be deleted from the database and we will move on to the next submission.")
            shutil.rmtree(self.framework_dest_dir)
            self.submission_to_transfer.delete()

        

    def _validate_md5sum(self, md5sum_source_dict):
        # For each key in the md5sum_source_dict we should have the hash
        for web_hash, file_path in md5sum_source_dict.items():
            local_file_path = os.path.join(self.framework_dest_dir, os.path.basename(file_path))
            try:
                framework_md5sum_hash = subprocess.run(['md5sum', local_file_path],
                                                       capture_output=True).stdout.decode("utf-8").split()[0]
            except FileNotFoundError:
                framework_md5sum_hash = subprocess.run(['md5', '-r', local_file_path],
                                                       capture_output=True).stdout.decode("utf-8").split()[0]
            if framework_md5sum_hash != web_hash:
                self.attempts += 1
                if self.attempts < 4:
                    # Delete all files in the framework dir
                    for del_filename in os.listdir(self.framework_dest_dir):
                        os.remove(os.path.join(self.framework_dest_dir, del_filename))
                    self._get_files_from_web_dir()
                else:
                    raise RuntimeError('Max number of download attempts exceeded. Could not verify md5sum')

    def _make_md5sum_web_server_dict(self):
        # There should be one md5sum in the pulled down files.
        # Make a new md5sum with the transfered files and verify
        md5sum_source_path = [
            os.path.join(self.framework_dest_dir, fn) for
            fn in os.listdir(self.framework_dest_dir) if fn.endswith('.md5sum')
        ]
        assert (len(md5sum_source_path) == 1)
        md5sum_source_path = md5sum_source_path[0]
        # Make a dict of the file
        with open(md5sum_source_path, 'r') as f:
            md5sum_source_dict = {line.split()[0]: line.split()[1] for line in f}
        return md5sum_source_dict

    def _get_files_from_web_dir(self):
        # Get every file that is in the web directory
        for get_file in self.sftp_client.listdir(self.web_source_dir):
            if '.' in get_file:
                if os.path.exists(os.path.join(self.framework_dest_dir, get_file)):
                    print(f'{os.path.join(self.framework_dest_dir, get_file)} already exists')
                else:
                    print(f'getting: {os.path.join(self.web_source_dir, get_file)}')
                    self.sftp_client.get(os.path.join(self.web_source_dir, get_file),
                                         os.path.join(self.framework_dest_dir, get_file))


twtf = TransferWebToFramework()
twtf.transfer()





