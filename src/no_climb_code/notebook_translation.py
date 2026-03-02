import pandas as pd
import datetime
import os
import subprocess
from collections import Counter
import re
from dateutil.relativedelta import relativedelta
from datetime import date
from pathlib import Path
import configparser


CONFIG_PATH = Path(__file__).resolve().parent.parent / "config"
print(CONFIG_PATH)


def load_config():
    config_parser = configparser.ConfigParser()
    config_parser.read(CONFIG_PATH)
    return config_parser

# set datestamp to today
datestamp = datetime.datetime.now().strftime('%Y%m%d')
# print(f"Today's date: {datestamp}")
datefound = False
datefound2 = False
day = 1

# USER INPUT: set path to folder containing __year_full_lineage_metadata.csv
def prevelance_data_path() -> str:

    config_parser = load_config()

    dir_prevalence_data = config_parser.get('file-paths', 'prevalence')
    return dir_prevalence_data


def is_local() -> bool:
    while True:
        local_file = input("local copy of genomics_cell_merged.csv? (yes/no): ").lower()
        if local_file not in ["yes", "no"]:
            print("Please type yes or no")
            continue
        else:
            break
    if local_file == ("no"):
        local_file = False
        print("loading genomics cell merged csv from shared-drive")
    elif local_file == ("yes"):
        local_file = True
        print("loading local copy of genomics cell merged csv")
    return local_file

# local_genomics_file = is_local()
# print(local_genomics_file)

def is_windows() -> bool:
    while True:
        system = input("operating system? (windows/linux): ").lower()
        if system not in ["windows" , "linux"]:
            print("Please type Windows or Linux!")
            continue
        else:
            break
    if system == ("windows"):
        system = True
        print("operating system is windows")
    elif system == ("linux"):
        system = False
        print("operating system is linux")
    return system

# op_system = is_windows()
# print(op_system)

def enter_genomics_folder_path() -> str:
    dir_genomics_merged = input("enter path to genomics cell merged csv parent directory (path): ")
    print(f"will load genomics cell merged csv from {dir_genomics_merged}")
    check = input(f"is path {dir_genomics_merged}? correct (yes/no): ")   #### PUT IN LOOP
    return dir_genomics_merged



def load_genomics_cell_merged(folder, datestamp, datefound, day):
    # set useful genomics merged columns
    columns = ['cog_uk_id','Specimen_Date_SK','Specimen_Number','finalid','cdr_specimen_request_sk','cdr_opie_id']
    # get genomics_cell_merged.csv
    while not datefound:
        try:
            # print(f"Genomics_cell_merged_{yesterday}.csv")
            linkage=pd.read_csv(os.path.join(folder, f"Genomics_cell_merged_{datestamp}.csv"), usecols=columns)
            datefound = True
        except:
            datestamp=(datetime.datetime.now() - datetime.timedelta(day)).strftime('%Y%m%d')
            day+=1
    print(f"linkage file loaded from {folder} with date {datestamp}")
    return linkage

def load_linelist(folder, datestamp, datefound, day):
    while not datefound:
        try:
            ll = pd.read_csv(os.path.join(folder,f"{datestamp}__year_full_lineage_metadata.csv"))
            datefound = True
        except:
            datestamp=(datetime.datetime.now() - datetime.timedelta(day)).strftime('%Y%m%d')
            day+=1
    print(f"linelist loaded from {folder} with date {datestamp}")
    return ll

def load_groupings(folder):
    proc = subprocess.Popen(f'ls -1 {os.path.join(folder,"groups*csv")} | tail -n 1', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = proc.communicate()
    json = stdout.decode('utf-8')
    json = json.rstrip()
    groups = pd.read_csv(json)
    print(f"grouping file loaded from {folder}")
    print(f"grouping file name {json}")
    return groups


dict_totidy = {'BA.1.17.2':'BA.1',
                'BA.1.1':'BA.1',
                'BA.5.1':'BA.5',
                'BA.5.2':'BA.5',
                'BA.5.2.1':'BA.5',
                'BQ.1.1':'BQ.1'}


def merged_linelist_with_linkage(linelist, groups, linkage, dict_totidy) -> pd.DataFrame:
    llgroups = linelist.merge(groups,left_on="lineage",right_on="lineage",how="left")
    llgroups.lineage_group = llgroups.lineage_group.fillna('Unassigned')
    lllinked = llgroups.merge(linkage,left_on="cog_id",right_on="cog_uk_id")
    print(Counter(lllinked.lineage_group))
    lllinked['lineage_group'] = lllinked['lineage_group'].replace(dict_totidy)
    print(Counter(lllinked.lineage_group))
    return lllinked

def merged_linelist_with_linkage_no_genomics(linelist, groups, dict_totidy) -> pd.DataFrame:
    llgroups = linelist.merge(groups,left_on="lineage",right_on="lineage",how="left")
    llgroups.lineage_group = llgroups.lineage_group.fillna('Unassigned')
    lllinked = llgroups
    print(Counter(lllinked.lineage_group))
    lllinked['lineage_group'] = lllinked['lineage_group'].replace(dict_totidy)
    print(Counter(lllinked.lineage_group))
    return lllinked
# merged = merged_linelist_with_linkage(linelist, groups, linkage, dict_totidy)
# linelist = load_linelist(dir_prevalence_data, datestamp, datefound, day)
# linkage = load_genomics_cell_merged(dir_genomics_merged, datestamp, datefound, day)
# groups = load_groupings(dir_prevalence_data)
# dir_genomics_merged = load_conditions()

def save_merge(lllinked, folder, datestamp):
    # lincols = ['cog_id','lineage','unaliased_lineage_x','lineage_group','lineages_version','Specimen_Date_SK','adm1','finalid','Specimen_Number','cdr_specimen_request_sk','cdr_opie_id']
    redll = lllinked
    redll = redll.rename(columns={"cog_id": "central_sample_id", "unaliased_lineage_x": "unaliased_lineage"})
    # redll['lineages_version']=linelist['lineages_version'].iloc[0] ## PUSHER-v1.36
    ###alternative lineages_version if required by COVE
    ##redll['lineages_version']='PANGO-v' +linelist['lineages_version'].iloc[0].split('-v')[-1] ## PANGO v-1.36
    redll['published_date']=''
    print(redll.shape)
    redll.to_csv(os.path.join(folder,f"{datestamp}_lineage_epi_line_list.csv"), index = False)
    print(f"file saved to path: {folder}")

dir_prevalence_data = prevelance_data_path()
# dir_genomics_merged = load_conditions()
linelist = load_linelist(dir_prevalence_data, datestamp, datefound, day)
# linkage = load_genomics_cell_merged(dir_genomics_merged, datestamp, datefound, day)
groups = load_groupings(dir_prevalence_data)
# merged = merged_linelist_with_linkage(linelist, groups, linkage, dict_totidy)
merged = merged_linelist_with_linkage_no_genomics(linelist, groups, dict_totidy)
save_merge(merged, dir_prevalence_data, datestamp)

print("now upload to COVE share-drive folder: smb://filecol19.phe.gov.uk/colindale_data/NISICC/WNCoV%20Epi%20Cell/Daily%20COVID19%20Epicell%20Line%20Lists/Variant%20data/genomics_line_lists)")
print("!!!remember to delete any local copies of genomics_cell_merged.csv due to PII!!!")