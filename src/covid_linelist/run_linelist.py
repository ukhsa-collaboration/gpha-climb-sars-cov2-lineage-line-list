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
from importlib import resources
from pathlib import Path

import covid_linelist.linelist_utilities as llu
from covid_linelist.lineage_collapser import LineageCollapser


# Arg parse setup
def get_args():
    """Get command line arguments. Arguments can be added or removed as
    required. It is however recommended to keep the arguments below as
    a minimum for development purposes."""
    parser = argparse.ArgumentParser(
        prog="run_linelist",
        description="""Program for generating the covid linelist including
        logic to collapse lineages into lineage groups if required.
        """,
    )
    parser.add_argument(
        "--input", "-i", type=str, required=True, help="Input csv of pangolin results"
    )
    parser.add_argument("--output", "-o", type=str, required=True, help="Folder to save results to")
    # takes a YYYY-MM-DD date as a string, parsing it to a datetime.date
    # then converts to give the corresponding monday date from that week
    # if not specified, uses todays date then converts that to give the
    # corresponding monday of this week
    parser.add_argument(
        "--reporting-date",
        "-r",
        type=lambda x: date_to_monday_date(dt.datetime.strptime(x, "%Y-%m-%d").date()),
        default=date_to_monday_date(dt.date.today()),
        help="Optional override for run date of the report. YYYY-MM-DD"
    )

    return parser.parse_args()


# Logger set up
def set_up_logger(log_file):
    """Sets up logging for the covid linelist code. Logging level is set
    to "INFO" and all logging messages are sent to a single log file. The
    logger is set to append mode so logs from older runs are not
    overwritten.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")

    out_handler = logging.FileHandler(log_file, mode="a")
    out_handler.setFormatter(formatter)
    logger.addHandler(out_handler)

    return logger


# Helper functions
def date_to_monday_date(input_date:dt.date) -> dt.date:
    """
    Takes in a datetime.date
    
    Finds the date corresponding to the Monday of the
    same week.

    Returns that date as a datetime.date
    """
    return input_date - dt.timedelta(days=input_date.weekday())


# Main function
def main():
    "Entry point for main run_linelist script"

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

    # Read in params from config file
    config_file = resources.files("covid_linelist.lib").joinpath("linelist_parameters.yaml")
    config_dict = llu.read_config_file(config_file)
    ### Set variable names from config
    # Size of period to split data into in days e.g. period_size 14 = 2 weeks. Needs to be
    # specified in days for add_reporting_period_column() to function correctly
    period_size: int = config_dict["global_params"]["period_size"]
    ## Params for % prevalence in recent timeframe
    # Number of weeks to include in recent reporting period timeframe e.g. 6 = 6 weeks
    timeframe_recent: int =  config_dict["recent_reporting_window"]["weeks_to_include"]
    # Max number of lineages to return in recent reporting period e.g. 8 = max of 8 lineages
    # can be returned by get_lineages_to_protect()
    max_lineages_recent: int = config_dict["recent_reporting_window"]["max_lineages"]
    # Min number of lineages to return in recent reporting period e.g. 1 = min of 1 lineages
    # has to be returned by get_lineages_to_protect()
    min_lineages_recent: int = config_dict["recent_reporting_window"]["min_lineages"]
    # Percent prevalence for lineage to be protected e.g. 5 = 5% prevalence in given period
    percent_prevalence: int = config_dict["recent_reporting_window"]["percent_prevalence"]
    ## Params for full timeframe
    # Number of weeks to include in full reporting timeframe e.g. 52 = 52 weeks
    timeframe_full: int = config_dict["full_reporting_window"]["weeks_to_include"]
    # Max number of lineages to return in the full reporting window. E.g. 14 = max of
    # 14 lineages to be returned across the whole reporting period.
    max_lineages_full: int = config_dict["full_reporting_window"]["max_lineages"]
    # Other variants to group as lineages of interest e.g. previously defined variants
    ## Lineage protection
    # List of lineages to collapse no further than e.g. collapse linegaes no further than BA.3
    collapse_limit_list: list[str] = config_dict["lineages_to_protect"]["lineage_collapse_limits"]
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
    # Can be overriden by --reporting-date argument
    range_end = args.reporting_date
    # Check and fill collection date column.
    # NOTE: - fills empty dates with Monday of current week - decide if correct behaviour
    lineage_df = llu.check_and_update_collection_dates(
        lineage_df=lineage_df,
        collection_date=range_end,
        fill_blanks=True
    )
    # Filter rows to required timeframe.
    lineage_df = llu.filter_lineage_df_to_date_cutoff(
        lineage_df=lineage_df,
        filter_end_date=range_end,
        previous_weeks_to_include=timeframe_full
    )
    # Add column to df with reporting periods
    lineage_df = llu.add_reporting_period_column(
        lineage_df=lineage_df,
        end_date=range_end,
        previous_weeks_to_include=timeframe_full,
        period_size=period_size
    )
    ## Calculate lineage counts and percentages for each reporting period
    # Group data to get lineage counts per week
    counts_by_reporting_period_df = llu.get_lineage_counts_per_period(lineage_df)
    # Add % column to lineage counts table
    counts_by_reporting_period_df = llu.add_percentages_column(counts_by_period_df=counts_by_reporting_period_df)
    ## Identify lineages to protect
    # Retrieve alias key json from COG-UK website
    pango_aliases_dict = llu.get_pango_aliases(do_filter=True)
    # Get reporting periods to protect within recent reporting window
    reporting_periods = llu.get_periods_to_protect(
        reporting_periods=counts_by_reporting_period_df['reporting_period'],
        end_date=range_end,
        timeframe_length=timeframe_recent
    )
    # Identify lineages to protect in the recent reporting window/last x weeks
    lineages_to_protect_list = llu.get_lineages_to_protect(
        counts_by_period_df=counts_by_reporting_period_df,
        periods_to_protect=reporting_periods,
        percent_threshold=percent_prevalence,
        max_lineages=max_lineages_recent,
        min_lineages=min_lineages_recent,
        pango_dict=pango_aliases_dict,
        lineage_collapse_limits=collapse_limit_list
    )
    # Add Unassigned to list as don't want to include in lineage to protect total here
    lineages_to_protect_list = list(set(lineages_to_protect_list + ["Unassigned"]))
    # Identify lineages to protect in the full reporting window, minus recent reporting window
    # covered above by get_lineages_to_protect.
    # E.g. samples from the last 52 weeks up until 6 weeks ago.
    # +1 to account for Unassigned as don't want this including in the additional lineages here
    num_additional_lineages = max_lineages_full - len(lineages_to_protect_list) + 1
    # Get lineages to protect for full reporting period, excluding recent reporting period
    total_df = llu.get_lineage_counts_for_full_window(counts_by_reporting_period_df, range_end, timeframe_recent)

    lineages_to_protect_list += llu.get_top_lineages_in_full_window(
        total_counts_df=total_df,
        already_protected=lineages_to_protect_list,
        lineage_collapse_limits=collapse_limit_list,
        pango_aliases=pango_aliases_dict,
        additional_lineages=num_additional_lineages,
        )
    ## Lineage collapsing after identifying lineages to protect
    # TODO: Add to own function
    # Collapse down any lineages not in lineages_to_protect_list
    lc = LineageCollapser(dataframe=counts_by_reporting_period_df,
                          lineages_col='lineage',
                          totals_col='seq_count',
                          min_level=1,
                          protect_lineages=lineages_to_protect_list,
                          collapsed_col='collapsed_alias',
                          pango_aliases=pango_aliases_dict
                          )
    # Get collapsed alias at correct levels for assigning groups
    collapse_threshold = len(lineage_df) + 1
    lc.collapse_based_on_threshold(threshold=collapse_threshold)

    # Mask anything not in lineages to protect list as 'Other'
    collapsed_masked_counts_df = llu.mask_less_prevalent_values(counts_df=lc.collapsed,
                                                                     lineages_to_leave_unmasked={
                                                                        "collapsed_alias" : lineages_to_protect_list
                                                                        }
                                                                    )
    ## Add additional lineage columns to sample info df
    # TODO: Add column with variant groups/non-prevalence based protection
    # Add unaliased lineage column
    lineage_df['unaliased_lineage'] = lc.alias_to_lineage(lineage_df.lineage)
    # Add lineage groups
    lineage_df = llu.add_lineage_group_to_metadata(lineage_df=lineage_df, counts_df=collapsed_masked_counts_df)
    # Write result files to csv
    # Masked and collapsed lineage week counts
    llu.write_to_csv(result_df = collapsed_masked_counts_df,
                      outdir=out_dir,
                      filename=f"{today_date}_year_lineage_masking.csv")
    # Sample records with unaliased lineage and collapsed alias groups
    llu.write_to_csv(result_df = lineage_df,
                     outdir=out_dir,
                     filename=f"{today_date}_year_full_lineage_metadata.csv")

    # Write to logs if component finished successfully (or not):
    logging.info("Linelist file successfully generated")

    return

# Run
if __name__ == "__main__":
    sys.exit(main())
