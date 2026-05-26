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

from . import linelist_utilities as llu
from .lineage_collapser import LineageCollapser


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
    today_date = dt.date.today()
    # Make output dir
    out_dir = Path(args.output) / f"{today_date}_covid_linelist"
    Path(out_dir).mkdir(parents=False, exist_ok=True)

    # Set up log file:
    log_file = Path(out_dir) / f"{today_date}_covid_linelist_logfile.txt"
    set_up_logger(log_file)

    ## Prepare sample data for prevalence calculations
    # Read in required columns from input csv with lineage designations
    try:
        lineage_df = llu.read_lineage_reports_csv(args.input)
        logging.info("Lineage reports file read in: %s", args.input)
    except FileNotFoundError:
        logging.error("Lineage reports file %s not found, exiting program", args.input)
        exitcode = 1
        return exitcode
    # Get Monday of current week as end of range for any time periods calculated.
    range_end = today_date - dt.timedelta(days=today_date.weekday())
    # Check and fill specimen date column.
    # NOTE: - fills empty dates with Monday of current week - decide if correct behaviour
    lineage_df = llu.check_and_update_specimen_dates(
        lineage_df=lineage_df,
        specimen_date=range_end,
        fill_blanks=True
    )
    #  Filter rows to required timeframe.
    lineage_df = llu.filter_lineage_df_to_date_cutoff(
        lineage_df=lineage_df,
        filter_end_date=range_end,
        previous_weeks_to_include=52 # TODO - add weeks to include to config if reimplement?
    )
    # Add column to df with reporting periods 
    lineage_df = llu.add_reporting_period_column(
        lineage_df=lineage_df,
        end_date=range_end,
        previous_weeks_to_include=52,
        period_size_in_days=14
    )
    ## Calculate lineage counts and percentages for each reporting period
    # Group data to get lineage counts per week
    counts_by_reporting_period_df = llu.get_lineage_counts_per_period(lineage_df)
    # Add % column to lineage counts table
    counts_by_reporting_period_df = llu.add_percentages_column(counts_by_period_df=counts_by_reporting_period_df)
    ## Identify lineages to protect
    # Retrieve alias key json from COG-UK website
    pango_aliases_dict = llu.get_pango_aliases(do_filter=True)
    # Get reporting periods to protect within
    reporting_periods = llu.get_periods_to_protect(
        reporting_periods=counts_by_reporting_period_df['reporting_period'],
        end_date=range_end,
        timeframe_length=6
    )
    # Identify lineages to protect in last x weeks
    lineages_to_protect_list = llu.get_lineages_to_protect(
        counts_by_period_df=counts_by_reporting_period_df,
        periods_to_protect=reporting_periods,
        percent_threshold=5,
        pango_dict=pango_aliases_dict
    )
    # Combine list of lineages to protect from collapse from each of the processes
    combined_protect = set(lineages_to_protect_list + ['Unassigned'])
    # # Identify lineages to protect in last 52 weeks up until 6 weeks ago
    # Based on overall prevalence in this period to get number of lineages to protect
    # up to x value
    # Collapse down lineages not in lineages_to_protect_list
    lc = LineageCollapser(dataframe=counts_by_reporting_period_df,
                          lineages_col='lineage',
                          totals_col='seq_count',
                          min_level=1,
                          protect_lineages=lineages_to_protect_list,
                          collapsed_col='collapsed_alias',
                          pango_aliases=pango_aliases_dict
                          )
    # Get collapsed alias at correct levels for assigning groups
    lc.collapse_recursively_to_at_least_n(n=10)
    # Mask anything not in combined_protect list as 'Other'
    collapsed_masked_counts_df = llu.mask_less_prevalent_values(counts_df=lc.collapsed,
                                                                     lineages_to_leave_unmasked={
                                                                        "collapsed_alias" : combined_protect
                                                                        }
                                                                    )
    ## Add additional lineage columns to sample info df   
    # Add unaliased lineage column
    lineage_df['unaliased_lineage'] = lc.alias_to_lineage(lineage_df.lineage)
    # Add lineage groups
    lineage_df = llu.add_lineage_group_to_metadata(lineage_df=lineage_df, counts_df=collapsed_masked_counts_df)
    # Write result files to csv
    # Masked and collapsed lineage week counts
    llu.write_to_csv(result_df = collapsed_masked_counts_df,
                      outdir=out_dir,
                      filename=f"{today_date}_year_lineage_masking.csv")
    # Sample records with unaliased lineage and collapsed alias
    llu.write_to_csv(result_df = lineage_df,
                      outdir=out_dir,
                      filename=f"{today_date}_year_full_lineage_metadata.csv")

    # Write to logs if component finished successfully (or not):
    logging.info("Linelist file successfully generated")

    return

# Run
if __name__ == "__main__":
    sys.exit(main())
