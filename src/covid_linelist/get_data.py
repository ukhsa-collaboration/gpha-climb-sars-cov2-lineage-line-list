import pandas as pd
import sys
import argparse
import os
import paramiko
import fnmatch
import pathlib
import logging

## Data from the SARS-CoV-2 pipeline will be sent to an sFTP
# sftpcol04.unix.phe.gov.uk port 443
# You can access with your email address and normal password but you need to use @phe.gov.uk instead of @ukhsa.gov.uk
# If you do not already have access then you need to put in a service deskrequest

## arguments

def cli():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, add_help=True)
    
    parser.add_argument(
        "--username",
        "-u",
        dest="username",
        type=str,
        required=True,
        help="user account username with phe as address, format = user.name@phe.gov.uk"
    )
    parser.add_argument(
        "--password",
        "-p",
        dest="password",
        type=str,
        required=True,
        help="user account password"
    )
    parser.add_argument(
        "--outdir",
        "-o",
        dest="outdir",
        type=str,
        required=False,
        default=f"{pathlib.Path().resolve()}/results",
        help="path to output directory"
    )
    args = parser.parse_args()
    return args


## connect to sFTP server

class sFTP():
    hostname = "sftpcol04.unix.phe.gov.uk"
    remote_data_folder = "/GPHA/covid_sequencing/"

    def __init__(self, username: str, password: str, outdir: str):
        self.username = username
        self.password = password
        self.outdir = outdir
        self.sftp = self.__connect_to_sFTP()
        self.__check_connection()
        self.get_sFTP_data()
        
    def __connect_to_sFTP(self):      
        host, port = self.hostname, 443
        transport = paramiko.Transport((host, port))
        transport.connect(None, self.username, self.password)  
        sftp = paramiko.SFTPClient.from_transport(transport)
        return sftp
    
    def __check_connection(self):
        host, port = self.hostname, 443
        logging.info("testing transport client")
        try:
            transport = paramiko.Transport((hostname, port))
            logging.info("successful transport client connection")
        except:
            status = "error creating transport client"
            logging.error(status)
            
        logging.info("testing ssh connection")
        try:
            transport.connect(None, self.username, self.password) 
            client = paramiko.SFTPClient.from_transport(transport)
            logging.info(f"successful ssh connection")
        except:
            status = "error with ssh connection"
            logging.error(status)
    
    
    def get_sFTP_data(self):
        # return list of sub-folders within the remote data folder path
        folders = self.sftp.listdir(path=self.remote_data_folder)
        logging.info(f"list of sub-folders found in path {self.remote_data_folder}: {", ".join(folders)}")
        pattern = '*.csv'
        # iterate over sub-folders for csv files and download to a local path copy        
        for i, folder in enumerate(folders):
            # if folder doesn't exist locally create it
            local_folder = f"{self.outdir}/{folder}"
            os.makedirs(local_folder, exist_ok=True)
            # return et csv files within remote folder
            remote_folder = f"{self.remote_data_folder}/{folder}"
            file_list = self.sftp.listdir(remote_folder)
            # return list for .csv files within sub-folder
            matched_files = [f for f in file_list if fnmatch.fnmatch(f, pattern)]
            logging.info(f"{len(matched_files)} csv files found in {remote_folder}: {", ".join(matched_files)}")
            for file_name in matched_files:
                # download files from remote folder -> local folder
                logging.info(f"processing: {remote_folder}/{file_name}")
                self.sftp.get(f'{remote_folder}/{file_name}', f'{local_folder}/{file_name}')
                logging.info(f"copied: {remote_folder}/{file_name} --> {local_folder}/{file_name}")
        # Need to handle error: OSError: Not a directory
        

## run process        
        
def main():
    # return command line inputs
    cmds = cli()
    
    # initiate connection to remote server and download of files. 
    connection = sFTP(
        username=cmds.username,
        password=cmds.password,
        outdir=cmds.outdir
    )
    connection.get_sFTP_data()

if __name__ == "__main__":
    sys.exit(main())