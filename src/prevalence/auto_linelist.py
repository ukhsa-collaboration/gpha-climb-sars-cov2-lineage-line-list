#!/usr/bin/env python3
from datetime import datetime, timedelta
import subprocess
import os
import glob
import shutil

from src.lineages.lineage_groups import generate_lineage_groups
from src.prevalence.lineage_prevalence import generate_lineage_prevalence


def make_output_folder():
    date_today = datetime.today()
    date_today = date_today.strftime("%Y%m%d")
    out_dir = f'{date_today}-covid-ll'
    print(f"{datetime.now()} Making output directory: {out_dir}")
    os.mkdir(out_dir)
    return out_dir

def move_files(out_dir):
    list_csv = glob.glob('*csv')
    list_json = glob.glob('*.json')
    str_csv = ", ".join(list_csv)
    str_json = ", ".join(list_json)
    print(f"{datetime.now()} Moving generated files to output directory: {str_csv} + {str_json} -> {out_dir}")
    for csv in list_csv:
        shutil.move(csv, out_dir)
    for jsons in list_json:
        shutil.move(jsons, out_dir)


def get_json():
    proc = subprocess.Popen(
        "ls -Art *json | tail -n 1",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
    )
    stdout, stderr = proc.communicate()
    json = stdout.decode("utf-8")#.replace("\\n", "")
    return json


def run_line_list(n_wks=6):
    run_start = datetime.today() - timedelta(days=n_wks * 7)
    run_start = datetime.date(run_start)
    run_start = run_start.strftime("%Y%m%d")
    print(f"{datetime.now()} Setting start date: {run_start}")
    generate_lineage_groups(end=run_start)

    #json=get_json()
    #generate_lineage_groups(start=run_start, threshold=0.05, filename=json)
    print(f"{datetime.now()} Lineage groups generated")


def run_lineage_prevalence(save_loc:str, path_to_alignment:str, path_to_metadata:str):
    print(f"{datetime.now()} Generating lineage prevalence data")
    generate_lineage_prevalence(file_path=path_to_alignment,
                                file_path2=path_to_metadata,
                                save_path=save_loc)


def run_commands():
    pass


def run_scan():
    climb_latest_alignments = "/cephfs/covid/bham/results/msa/latest/alignments/"
    climb_latest_general = "/cephfs/covid/artifacts/elan/latest/"
    out_dir = make_output_folder()
    run_line_list()
    run_lineage_prevalence(save_loc=out_dir,
                           path_to_alignment=climb_latest_alignments,
                           path_to_metadata=climb_latest_general)
    move_files(out_dir)
    print(f"{datetime.now()} Lineage prevalence and line list files should be ready")

# if __name__ == "__main__":
#     n = 6
#     start = datetime.today() - timedelta(days=n * 7)
#     start = datetime.date(start)
#     start = start.strftime("%Y%m%d")
#
#     print(start)
#
#     ## pre-flight checks
#     ## do we have all the scripts we need?
#     assert os.path.isfile("lineage_groups/generate_lineage_groups.py")
#     ##assert os.path.isfile(
#     ##    "sars_cov2_lineage_prevalence_climb/lineage_prevalence_CLIMB4.py"
#     ##)
#     assert os.path.isfile("linprev_12months.py") or os.path.isfile(
#         "linprev_12months.py"
#     )
#
#     ## create the lineage groups, needs to run with two different date ranges / thresholds
#     cmd1 = f"python lineage_groups/generate_lineage_groups.py --end {start} --threshold 0.25"
#     subprocess.run(
#         cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True
#     )
#
#     proc = subprocess.Popen(
#         "ls -Art *json | tail -n 1",
#         stdout=subprocess.PIPE,
#         stderr=subprocess.PIPE,
#         shell=True,
#     )
#     stdout, stderr = proc.communicate()
#     json = stdout.decode("utf-8").replace("\\n", "")
#
#     cmd2 = f"python lineage_groups/generate_lineage_groups.py --start {start} --threshold 0.05 {json}"
#     print("Lineage groups generated...")
#
#     ## run the prevalence tables for the lineage line list
#     if os.path.isfile("linprev_12months.py"):
#         ## if the script is in this directory
#         cmd3 = "python linprev_12months.py"
#     else:
#         ## else assume it is in a subfolder
#         cmd3 = "python climb-lineage-line-list/linprev_12months.py"
#     print("Lineage prevalence for line list generated...")
#
#     ## run the prevalence tables for the 12-week prevalence plot
#     ## cmd4 = "python sars_cov2_lineage_prevalence_climb/lineage_prevalence_CLIMB4.py"
#
#     ## actually run the rest of the commands
#     subprocess.run(
#         cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True
#     )
#     subprocess.run(
#         cmd3, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True
#     )
#     ##subprocess.run(
#     ##    cmd4, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True
#     ##)
#     print("Lineage 12 week prevalence tables generated...")
#     print("Lineage prevalence and line list files should be ready")
#
#     today = datetime.today()
#     today = today.strftime("%Y%m%d")
#     folder = f'{today}-covid-ll'
#     os.mkdir(folder)
#     csv_list = glob.glob('*csv')
#     json_list = glob.glob('*.json')
#
#     for file in csv_list:
#         shutil.move(file, folder)
#     for file in json_list:
#         shutil.move(file, folder)
