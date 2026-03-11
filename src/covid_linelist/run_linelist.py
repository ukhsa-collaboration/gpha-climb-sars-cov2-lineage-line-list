#!/usr/bin/env python3

"""
Main script to generate covid linelist groupings from covid data received from
core bioinformatics.
"""

# Imports
import argparse
import datetime as dt
import logging
import sys
from pathlib import Path

import linelist_utilities as llu


# Arg parse setup
def get_args():
    """Get command line arguments. Arguments can be added or removed as
    required. It is however recommended to keep the arguments below as
    a minimum for development purposes."""
    parser = argparse.ArgumentParser(
        prog="covid linelist",
        description="""Program for generating covid line list including groupings.
        """,
    )
    parser.add_argument(
        "--input", "-i", type=str, required=True, help="Input csv of pangolin results"
    )
    parser.add_argument("--output", "-o", type=str, required=True, help="Folder to save results to")

    return parser.parse_args()


# Logger set up
def set_up_logger(stdout_file):
    """Example logger set up which can be amended as required. In this example,
    all logging messages go to a stdout log file, and error messages also go to
    stderr log. If the component runs correctly, stderr is empty. The logger is
    set to append mode so logs from older runs are not overwritten.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")

    out_handler = logging.FileHandler(stdout_file, mode="a")
    out_handler.setFormatter(formatter)
    logger.addHandler(out_handler)

    return logger


# Main function
def main():
    "Main function description here"

    # Retrieve command line arguments:
    args = get_args()  # noqa: F841

    # Get today's date
    today_date = dt.date.today().strftime("%Y-%m-%d")

    # Make output dir
    out_dir = Path(args.output) / f"{today_date}_covid_linelist"
    Path(out_dir).mkdir(parents=False, exist_ok=True)

    # Set up log file:
    log_file = Path(out_dir) / f"{today_date}_covid_linelist_logfile.txt"
    set_up_logger(log_file)

    # Read in required columns from input csv with lineage designations
    try:
        lineage_df = llu.read_lineage_reports_csv(args.input)
        logging.info("Lineage reports file read in: %s", args.input)
    except FileNotFoundError:
        logging.error("Lineage reports file %s not found, exiting program", args.input)
        exitcode = 1
        return exitcode
    # Check and fill collection date column.
    # NOTE: - fills empty dates with today's date - decide if correct behaviour
    lineage_df = llu.check_and_update_specimen_dates(
        lineage_df=lineage_df, specimen_date=today_date, fill_blanks=True
    )
    #  Add column to df with week begin date
    lineage_df = llu.add_week_begin_column(lineage_df=lineage_df)

    #  Filter rows to required timeframe.
    filtered_df = llu.filter_lineage_df_to_date_cutoff(
        lineage_df=lineage_df, filter_start_date=today_date, previous_weeks_to_include=52
    )  # TODO - add weeks to include to config if reimplement?

    ## Linelist generation
    # Group data to get lineage counts per week
    counts_by_week_df = llu.get_lineage_counts_per_week(filtered_df)
    # Add timeframe to protect lineages in as column
    counts_by_week_df = llu.add_timeframe_to_protect_column(
        counts_by_week_df=counts_by_week_df, timeframe_start=today_date, timeframe_length=6
    )
    # Add % column to lineage counts table
    counts_by_week_df = llu.add_percentages_column(counts_by_week_df=counts_by_week_df)
    # Retrieve alias key json from COG-UK website
    pango_aliases_dict = llu.get_pango_aliases(do_filter=True)
    # Identify lineages to protect
    lineages_to_protect_list = llu.get_lineages_to_protect(counts_by_week_df=counts_by_week_df,
                                                           timeframe_length=6,
                                                           percent_threshold=5,
                                                           pango_dict=pango_aliases_dict
                                                           )

    # Write to logs if component finished successfully (or not):
    logging.info("Linelist file successfully generated")

    return


# Run
if __name__ == "__main__":
    sys.exit(main())
