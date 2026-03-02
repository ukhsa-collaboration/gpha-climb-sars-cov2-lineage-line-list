#!/usr/bin/env python3
from datetime import datetime, timedelta
import subprocess
import os
import glob
import shutil
import glob
import configparser
import importlib.resources
from pathlib import Path
import src.no_climb_code.lineage_prevalence_class_no_climb
import src.prevalence.auto_linelist as autoll
from src.lineages.lineage_groups import generate_lineage_groups

from src.no_climb_code.lineage_prevalence_no_climb import generate_lineage_prevalence

CONFIG_PATH = Path(__file__).resolve().parent.parent / "config"
print(CONFIG_PATH)





def load_config():
    config_parser = configparser.ConfigParser()
    config_parser.read(CONFIG_PATH)
    return config_parser
def run_scan():
    config_parser = load_config()

    latest_alignments = config_parser.get('file-paths', 'no_climb')
    print(latest_alignments)

    autoll.run_line_list(metadata_path=latest_alignments)

    generate_lineage_prevalence(file_path=latest_alignments)

    print(f"{datetime.now()} Lineage prevalence and line list files should be ready")


if __name__ == '__main__':
    run_scan()