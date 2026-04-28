import pandas as pd
import sys
import argparse
import os
import paramiko
import fnmatch
import pathlib
import logging
import glob
import datetime
import re

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
        # self.get_sFTP_data()
        
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
    
    
    def get_sFTP_data(self, pattern: str):
        # return list of sub-folders within the remote data folder path
        folders = self.sftp.listdir(path=self.remote_data_folder)
        logging.info(f"list of sub-folders found in path {self.remote_data_folder}: {", ".join(folders)}")
        #pattern = '*.csv'
        # iterate over sub-folders for csv files and download to a local path copy   
        for i, folder in enumerate(folders):
            # check formate is as expected, i.e. check folder name contains a digit (date)
            if bool(re.search(r'\d', folder)) is True:
                # if folder doesn't exist locally create it
                local_folder = f"{self.outdir}/{folder}"
                os.makedirs(local_folder, exist_ok=True)
                # return et csv files within remote folder
                remote_folder = f"{self.remote_data_folder}/{folder}"
                file_list = self.sftp.listdir(remote_folder)
                # return list for .csv files within sub-folder
                matched_files = [f for f in file_list if fnmatch.fnmatch(f, pattern)]
                logging.info(f"{len(matched_files)} csv files found in {remote_folder}: {", ".join(matched_files)}")#
                if len(matched_files) == 0:
                    continue
                else:
                    for file_name in matched_files:
                        # download files from remote folder -> local folder
                        logging.info(f"processing: {remote_folder}/{file_name}")
                        self.sftp.get(f'{remote_folder}/{file_name}', f'{local_folder}/{file_name}')
                        logging.info(f"copied: {remote_folder}/{file_name} --> {local_folder}/{file_name}")
            else:
                continue
        # Need to handle error: OSError: Not a directory


def identify_ont_folders(parent_folder: str):
    sub_folders = glob.glob(parent_folder + "/*/")
    #print(sub_folders)
    ont_folders = []
    illumina_folders = []
    for x, folder in enumerate(sub_folders):
        date_identifier = folder.split("/")[-2]
        # print(date_identifier)
        date_identifier = date_identifier.split("_")[0]
        # print(date_identifier)
        date_identifier_length = len(date_identifier)
        date_identifier_is_no = date_identifier.isnumeric()
        # print(date_identifier_is_no)
        if date_identifier_length == 8 and date_identifier_is_no is True:
            ont_folders.append(folder)
        elif date_identifier_length != 8 and date_identifier_is_no is True:
            illumina_folders.append(folder)
    logging.info(f"ont folders discovered: {", ".join(ont_folders)}")
    logging.info(f"illumina folders discovered: {", ".join(illumina_folders)}")
    return ont_folders, illumina_folders


## ONT processing


def process_ont_results_df(ont_results_df, sample_sheet) -> pd.DataFrame:
    # specify useful columns
    cols = ["taxon", "lineage", "scorpio_call", "version", "pangolin_version", "scorpio_version", "qc_status"]
    # import the datafame skipping the 2nd row
    df = pd.read_csv(ont_results_df, usecols=cols, skiprows=[1])
    
    df = match_ont_csv_with_samplesheet(ont_results_filename=ont_results_df,
                                        sample_sheet=sample_sheet,
                                        ont_results_df=df)

    # get date via string split and taking the last 8 digits
    df = df.assign(collection_date=str(df["taxon"]).split("_")[0][-8:])
    #df = df.assign(central_sample_id=str(df["taxon"]).split("_")[3:4])
    df = df.assign(central_sample_id=df["molis_id"])#.split("\n")[0])
    return df
    

def match_ont_csv_with_samplesheet(ont_results_filename:str, sample_sheet:str, ont_results_df:pd.DataFrame) -> pd.DataFrame:
    barcode = ont_results_filename.split("barcode")[1].split(".")[0]
    #print(barcode)
    col_names = ["barcode", "molis_id"]
    df_samplesheet = pd.read_csv(sample_sheet, names=col_names, header=None, converters={'barcode': str})
    #print(df_samplesheet)
    df = ont_results_df.copy()
    df = df.assign(barcode=str(df["taxon"]).split("barcode")[1].split(".")[0][:2])
    #print(df.barcode)
    df = df.merge(df_samplesheet, on="barcode")
    #print(df)
    return df


def process_ont_results_folder(ont_results_folder):
    processed_dfs = []
    illumina_results_dfs = glob.glob(ont_results_folder + "/*report.csv")
    sample_sheet = glob.glob(ont_results_folder + "/*samplesheet.csv")[0]
    #print(sample_sheet)
    for x, file in enumerate(illumina_results_dfs):
        df = process_ont_results_df(ont_results_df=file,
                                    sample_sheet=sample_sheet)
        processed_dfs.append(df)
    df = pd.concat(processed_dfs)
    return df
    

def remove_ont_controls(df: pd.DataFrame) -> pd.DataFrame:
    identifier = "positive|negative|water|control|pos|neg|ctr"
    df = df[~df["molis_id"].str.contains(identifier, case=False)]
    return df
    
    
def process_ont_results(list_ont_results_folders: list) -> pd.DataFrame:
    concat_ont_results_df = []
    for x, folder in enumerate(list_ont_results_folders):
        df = process_ont_results_folder(ont_results_folder=folder)
        concat_ont_results_df.append(df)
    df = pd.concat(concat_ont_results_df)
    #print(df)
    df = remove_ont_controls(df)
    df = df.drop("barcode", axis =1)
    #print(df)
    return df
    
## Illumina processing

def return_date_from_illumina_folder(folder_name: str) -> val:
    depth = folder_name.count("/") -1
    #print(depth)
    date  = folder_name.split("/")[depth].split("_")[0]
    #print(folder_name)
    date = "20" + str(date)
    logging.info(f"retrieved date from {folder_name} as {date}")
    return date


def process_illumina_results_df(illumina_results_df, date) -> pd.DataFrame:
    cols = ["taxon", "lineage", "scorpio_call", "version", "pangolin_version", "scorpio_version", "qc_status"]
    df = pd.read_csv(illumina_results_df, usecols=cols)
    df = df.assign(collection_date=date)

    # central sample id is first part of taxon column
    # central_sample_id = str(df["taxon"]).split("_")[1].split("-")[0]
    central_sample_id = str(df["taxon"]).split("-")[0].split(" ")[0]
    central_sample_id = df["taxon"].str.split("-", expand=True)[0]
    df = df.assign(central_sample_id=central_sample_id)
    
    # molis id is 2nd part of taxon column and is 10 digits long
    # molis_id = str(df["taxon"]).split("_")[1][:10]
    molis_id = str(df["taxon"]).split("_")[1].split(".")[0].split("\nN")[0]
    df = df.assign(molis_id=molis_id)
    # step code is the suffix of the taxon column
    # step_code = str(df["taxon"]).split("_")[1][10:].split("\nN")[0].strip("-")
    # df = df.assign(step_code=step_code)
    return df


def process_illumina_results_folder(illumina_results_folder):
    date = return_date_from_illumina_folder(folder_name=illumina_results_folder)
    processed_dfs = []
    illumina_results_dfs = glob.glob(illumina_results_folder + "/*.csv")
    for x, file in enumerate(illumina_results_dfs):
        df = process_illumina_results_df(illumina_results_df=file,
                                        date=date)
        processed_dfs.append(df)
    df = pd.concat(processed_dfs)
    return df
    

def remove_illumina_controls(df: pd.DataFrame) -> pd.DataFrame:
    identifier = "positive|negative|water|control|pos|neg|ctr"
    df = df[~df["taxon"].str.contains(identifier, case=False)]
    return df


def process_illumina_results(list_of_illumina_folders: list) -> pd.DataFrame:
    concat_illumina_results_df = []
    for x, folder in enumerate(list_of_illumina_folders):
        df = process_illumina_results_folder(illumina_results_folder=folder)
        concat_illumina_results_df.append(df)
    df = pd.concat(concat_illumina_results_df)
    df = remove_illumina_controls(df=df)
    return df
    

def process_results(local_dir: str) -> pd.DataFrame:
    ont_folders, illumina_folders = identify_ont_folders(parent_folder=local_dir)
    df_illumina = process_illumina_results(illumina_folders)
    df_ont = process_ont_results(ont_folders)
    df_results = pd.concat([df_illumina, df_ont], ignore_index=True)
    df_results = df_results.reindex(columns=[*df_results.columns.tolist(), 'Specimen_Number', 'cdr_specimen_request_sk', 'cdr_opie_id'])
    df_results["molis_id"] = df_results["molis_id"].str[0:10]
    # reorder columns
    # ...
    return df_results


# concat fasta files

def concat_fasta(local_dir: str, date:str):
    sub_folders = glob.glob(local_dir + "/*/")
    sub_folders = [x for x in sub_folders if not x.startswith(f'{local_dir}/test.txt')]
    identifier = "positive|negative|water|control|pos|neg|ctr"
    fastas = []
    for x, folder in enumerate(sub_folders):
        print(folder)
        fasta_files = glob.glob("*.fa*")
        fastas.append(fasta_files)
        print(fastas)
    
    output_fas = open(f'{local_dir}/{date}_covid_ll.fasta', 'w')
    file_count = len(fastas)
    print(file_count)
    for f in fastas:
        print(f)
        file_count += 1
        fh = open(f)
        for line in fh:
            output_fas.write(line)
        fh.close()
    output_fas.close()



## run process        
        
def main():
    # return command line inputs
    cmds = cli()     
    
    # initiate connection to remote server and download of files. 
    #connection = sFTP(
    #    username=cmds.username,
    #    password=cmds.password,
    #    outdir=cmds.outdir
    #)
    #connection.get_sFTP_data_fa()
    #connection.get_sFTP_data(pattern='*.csv')
    #connection.get_sFTP_data(pattern='*.fa*')
    
    df = process_results(local_dir=cmds.outdir)
    date = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    df.to_csv(f"{cmds.outdir}/{date}_covid_ll.csv")
    
    #concat_fasta(local_dir=cmds.outdir, date=date)
    
if __name__ == "__main__":
    sys.exit(main())