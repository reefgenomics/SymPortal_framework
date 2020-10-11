#!/usr/bin/env python3
"""
This script will be managed by chron jobs and will be run once every hour
It will be responsible for transfering successfuly submitted user files from the web sever to the framework server
It will check for Submission objects that have a status of progress_status of submitted
At the start of the script it will do a pgrep to check that the script is not currently running to prevent
stacking up of the same script.
The script will use the md5sum that is generated at the time of user upload to verify the integrity of the transfer.
After transfer is complete the status of the Submission objects that have been transfered will be updated to
tranfer_to_framework_server_complete. The transfer_to_framework_server_date_time will be logged.
A seperate chron job will handle loading the transferred submissions
"""

import subprocess
import sys
import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from dbApp.models import Submission
import sp_config
import paramiko
from datetime import datetime

class TransferWebToFramework:
    def __init__(self):
        # The very first thing to do is to perform a pgrep to see if another instance of this script is being run
        captured_output = subprocess.run(['pgrep', 'transfer_web_framework.py'], check=True, capture_output=True)
        if captured_output:
            # Then the script is already running (i.e. there is a transfer in progress)
            sys.exit()
        # Get a list of the Submissions that need to be transfered
        self.submissions_to_transfer = Submission.objects.filter(progress_status="submitted", error_has_occured=False)
        self.symportal_data_dir = sp_config.symportal_data_dir

        # User paramiko to set up an sftp that we can use to transfer
        self.ssh_client = paramiko.SSHClient()
        self.ssh_client.load_system_host_keys()
        self.ssh_client.connect(hostname=sp_config.web_ip, username=sp_config.web_user,
                                password=sp_config.web_pass)
        # Open sftp client
        self.sftp_client = self.ssh_client.open_sftp()

        # Transfer each submission
        for sub_to_trans in self.submissions_to_transfer:
            web_source_dir = sub_to_trans.web_local_dir_path
            framework_dest_dir = os.path.join(self.symportal_data_dir, os.path.basename(web_source_dir))

            # Get every file that is in the web directory
            for get_file in self.sftp_client.listdir(web_source_dir):
                if '.' in get_file:
                    if os.path.exists(os.path.join(framework_dest_dir, get_file)):
                        print(f'{os.path.join(framework_dest_dir, get_file)} already exists')
                    else:
                        print(f'getting: {os.path.join(web_source_dir, get_file)}')
                        self.sftp_client.get(os.path.join(web_source_dir, get_file), os.path.join(framework_dest_dir, get_file))

            # There should be one md5sum in the pulled down files.
            # Make a new md5sum with the transfered files and verify
            md5sum_source_path = [fn for fn in os.listdir(framework_dest_dir) if fn.endswith('.md5sum')]
            assert(len(md5sum_source_path) == 1)
            md5sum_source_path = md5sum_source_path[0]
            # Make a dict of the file
            with open(md5sum_source_path, 'r') as f:
                md5sum_source_dict = {line.split()[0]: line[1] for line in f}

            # For each key in the md5sum_source_dict we should have the hash
            for file_path, web_hash in md5sum_source_dict.items():
                local_file_path = os.path.join(framework_dest_dir, os.path.basename(file_path))
                try:
                    framework_md5sum_hash = subprocess.run(['md5sum', local_file_path],
                                           capture_output=True).stdout.decode("utf-8").split()[1]
                except FileNotFoundError:
                    framework_md5sum_hash = subprocess.run(['md5', '-r', local_file_path],
                                           capture_output=True).stdout.decode("utf-8").split()[1]
                if framework_md5sum_hash != web_hash:
                    # TODO reattmept the transfer
                    raise RuntimeError('Hashes do not match')

            # If we get here then the md5sum has been validated
            # The transfer is complete
            # Delete the remote files and the remote directory
            self.sftp_client.rmdir(web_source_dir)

            # Update the time of the Submission object transfer_to_framework_server_date_time
            sub_to_trans.transfer_to_framework_server_date_time = str(
                datetime.utcnow()
            ).split('.')[0].replace('-','').replace(' ','T').replace(':','')

            # Update the status of the Submission object
            sub_to_trans.progress_status = "transfer_to_framework_server_complete"







