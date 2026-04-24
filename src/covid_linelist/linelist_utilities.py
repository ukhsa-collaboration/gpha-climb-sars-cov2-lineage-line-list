#!/usr/bin/env python3

"""
Module containing functions used in this repo
"""

# Imports - ordered (can use ruff to do this automatically)
import datetime as dt
import json
import logging
import os
import sys
from pathlib import Path

import pandas as pd
import requests
import yaml

from covid_linelist.lineage_collapser import LineageCollapser


# Functions
def read_config_file(config_file: Path) -> dict:
    """Reads config file containing parameters for linelist code.
    Arguments:
        config_file -- yaml file containing parameters
    Returns:
        dictionary of linelist parameters
    """
    with Path(config_file).open("r") as file:
        linelist_params = yaml.safe_load(file)
    return linelist_params

def read_lineage_reports_csv(input_file: os.path) -> pd.DataFrame:
    """Read input covid lineage designation file and return dataframe
    containing only required columns and records.
    Arguments:
        input_file -- Csv input file containing pangolin lineage designations
                      and sample information.
    Outputs:
        input_df -- Dataframe containing information required for linelist generation
    """
    # Specify columns to keep from input csv
    cols_to_keep = [
        "taxon", "sample_id", "central_sample_id", "molis_id", "collection_date",
        "lineage", "scorpio_call", "version", "pangolin_version", "scorpio_version",
        "qc_status","Specimen_Number","cdr_specimen_request_sk","cdr_opie_id"
        ]
    with Path(input_file).open("r") as file:
        input_df = pd.read_csv(file, 
                               usecols=lambda x: x in cols_to_keep,
                               dtype={'Specimen_Number': str,
                                      'cdr_specimen_request_sk': str,
                                      'cdr_opie_id': str}
                               )
    return input_df

def check_and_update_collection_dates(lineage_df: pd.DataFrame, collection_date: str, fill_blanks: bool) -> pd.DataFrame:
    """Check collection date column exists in input dataframe. If no date
    column present in the file, uses the date provided to generate the date
    column. Also optionally fills in any blank dates with today's date.
    Arguments:
        lineage_df -- Dataframe containing pangolin lineage designations
                      and sample information.
    Outputs:
        lineage_df -- Updated dataframe with checked and filled specimen date column
    """
    if 'collection_date' not in lineage_df.columns:
        lineage_df['collection_date'] = collection_date
        logging.info("""No specimen date column in input dataframe. Adding
                     collection date of %s""", collection_date)
    elif fill_blanks:
        lineage_df['collection_date'].fillna(collection_date)
        logging.info("Filling in blank values in specimen date column with %s", collection_date)

    lineage_df['collection_date'] = pd.to_datetime(lineage_df['collection_date'], format='%Y%m%d')

    return lineage_df

def filter_lineage_df_to_date_cutoff(lineage_df: pd.DataFrame, filter_end_date: str, previous_weeks_to_include: int) -> pd.DataFrame:
    """Filter the lineage_df to the required date range.
    Arguments:
        lineage_df -- Dataframe containing pangolin lineage designations
                      and sample information.
        filter_end_date -- Date to calculate cut off range from.
        previous_weeks_to_include -- Number of weeks to include in the
                                     filtered dataframe. Rows with specimen
                                     dates above this will be filtered out.
    Outputs:
        filtered_df -- lineage_df filtered to date cut off specified
    """
    date_x_weeks_ago = (
        filter_end_date - dt.timedelta(
            days=filter_end_date.weekday(),
            weeks=previous_weeks_to_include
            )).strftime("%Y-%m-%d")
    filtered_df = lineage_df[lineage_df['collection_date'] >= date_x_weeks_ago]
    removed_rows = (len(lineage_df) - len(filtered_df))
    logging.info("""Dataframe filtered to include samples with collection date
                 later than %s. %d rows of data removed from the dataframe,
                 %d remaining.
                 """, date_x_weeks_ago, removed_rows, len(filtered_df))
    return filtered_df

def add_reporting_period_column(
    lineage_df:pd.DataFrame,
    end_date: str,
    previous_weeks_to_include: int,
    period_size: int
    ) -> pd.DataFrame:
    """Add reporting_period column to the lineage dataframe. This replaces the previous
       week_begin column and makes it generic so period can be of varying lengths as required.
        Arguments:
        lineage_df -- Dataframe containing pangolin lineage designations
                      and sample information.
        end_date -- Date to calculate cut off range from.
        previous_weeks_to_include -- Number of weeks to include in the
                                     date range
        period_size -- Size of periods in days to split the date range into
    Outputs:
        filtered_df -- lineage_df with reporting_period column
    """
    # Get reporting periods as a series
    start_date = end_date - dt.timedelta(weeks=previous_weeks_to_include)
    date_range = pd.date_range(start=start_date,
                               end=end_date,
                               freq=f"{period_size}D"
                               ).to_series(name="reporting_period")
    date_range = date_range.astype('datetime64[us]')
    # Sort lineage_df so merge_asof works
    lineage_df.sort_values(by="collection_date", inplace=True)
    # Merge lineage_df with the reporting period series - this allocates a reporting period
    # by going backwards from the collection date until a reporting period start date is
    # encountered and adds this as a column 'reporting_period'.
    lineage_df = pd.merge_asof(lineage_df,
                               date_range,
                               left_on="collection_date",
                               right_on="reporting_period"
                               )
    return lineage_df

def get_lineage_counts_per_period(lineage_df: pd.DataFrame) -> pd.DataFrame:
    """Create a grouped dataframe containing counts per lineage per week
    Arguments:
        lineage_df -- Dataframe containing pangolin lineage designations
                      and sample information.
    Outputs:
        counts_by_period_df -- Dataframe containing lineage counts per week
    """
    counts_by_period_df = lineage_df.groupby(
        ['reporting_period', 'lineage']).size().to_frame('seq_count').reset_index()
    return counts_by_period_df

def add_percentages_column(counts_by_period_df: pd.DataFrame) -> pd.DataFrame:
    """Add boolean column specifying if week falls in date range specified
    Arguments:
        counts_by_week_df -- Dataframe containing lineage counts per week
        timeframe_end -- Most recent date to include in timeframe
        timeframe_length -- Number of weeks to include in the timeframe range
    Outputs:
        counts_by_period_df -- Updated dataframe containing lineage counts per week
    """
    period_total = counts_by_period_df.groupby('reporting_period')["seq_count"].transform('sum')
    counts_by_period_df["pct_of_reporting_period"] = (counts_by_period_df["seq_count"] / period_total).mul(100)
    return counts_by_period_df

def get_pango_aliases(do_filter: bool=True) -> dict:
    """
    Return dict of Pango lineage alias: full lineage name from COG-UK website.
    Arguments:
        do_filter -- Bool of whether to filter alias dict
    Outputs:
        pango_aliases_dict -- Dictionary of pango aliases
    """
    # TODO: take url to config?
    try:
        pango_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json"
        with requests.get(pango_url) as url:
            pango_aliases_dict = json.loads(url.text)
            if do_filter:
                pango_aliases_dict = dict(
                    [t for t in pango_aliases_dict.items() if not t[0].startswith('X') and t[0] not in ['A', 'B']])
                return pango_aliases_dict
            else:
                return pango_aliases_dict
    except Exception:
        logging.error("Could not generate alias dict from remote json file %s", pango_url)
        sys.exit(1)

def get_periods_to_protect(reporting_periods: pd.Series, end_date: dt.date, timeframe_length: int) -> pd.DatetimeIndex:
    """Identify periods in the dataframe to protect based on timeframe.
    Arguments:
        counts_by_week_df -- Dataframe containing lineage counts per week
        timeframe_end -- Most recent date to include in timeframe
        timeframe_length -- Number of weeks to include in the timeframe range
    Outputs:
        periods_to_protect -- DatetimeIndex containing periods to protect
    """
    start_date = end_date - dt.timedelta(weeks=timeframe_length)
    reporting_periods = pd.to_datetime(
        reporting_periods[
            (reporting_periods >= pd.to_datetime(start_date)) &
            (reporting_periods < pd.to_datetime(end_date))
    ].unique(),
    format='%Y-%m-%d'
    )
    return reporting_periods

def get_lineages_above_threshold(
    counts_by_period_df: pd.DataFrame,
    period_to_protect: int,
    percent_threshold: int,
    lineage_column: str
    ) -> list[str]:
    """Check which lineages should be protected based on reaching a percent threshold in a given time period.
    Used by get_lineages_to_protect to re-check lineages reaching protection during collapsing process when trying to reach
    max or min number of lineages required.
    Arguments:
        counts_by_period_df -- Dataframe containing lineage counts per reporting period
        period_to_protect -- Date period to protect
        percent_threshold -- Percent prevalence threshold for lineage to be protected
        lineage_column -- Name of column with lineage info
    Returns:
        to_protect -- List of lineages to protect based on prevalence being above or equal to the percent threshold
    """
    # Check for lineages that qualify for protection
    to_protect = []
    period_df = counts_by_period_df[counts_by_period_df['reporting_period'] == period_to_protect]
    over_threshold_list = period_df[period_df['pct_of_reporting_period'] >= percent_threshold][lineage_column].to_list()
    to_protect += over_threshold_list
    return to_protect

def get_lineages_to_protect(
    counts_by_period_df: pd.DataFrame,
    periods_to_protect: int,
    percent_threshold: int,
    max_lineages: int,
    min_lineages: int,
    pango_dict: dict
    ) -> list[str]:
    """Identify lineages to protect from collapsing in recent reporting period
    based on prevalence.
    Arguments:
        counts_df -- Dataframe containing per week counts of lineages
        timeframe_length -- Timeframe to protect. Used to determine column name in_last_{x}_weeks
        percent_threshold - % of samples of a given lineage needed to protect a lineage
        max_lineages -- Maximum number of lineages to return in lineage list
        min_lineages - Minimum number of lineages to return in lineage list
        pango_dict -- Dict containing pango aliases from COG-UK
    Outputs:
        lineage_list -- List of lineages that should be protected
    """
    lineage_list = []
    # Get initial lineage list based on lineages above the % threshold in each time period before any collapsing
    for period in periods_to_protect:
        lineage_list += get_lineages_above_threshold(counts_by_period_df, period, percent_threshold, "lineage")
    lineage_list = list(set(lineage_list))
    # If Unassigned in initial list, +1 to the max and min lineage values as don't want
    # to include Unassigned in lineage total.
    if "Unassigned" in lineage_list:
        logging.info("""Unassigned in recent reporting window lineage list, adding additional lineage to min
                     and max lineage number to account.""")
        max_lineages += 1
        min_lineages += 1
    # If less than the min number of lineages are over the % threshold, collapse until at least  the minimum
    # number of lineages reaches 5% prevalence:
    if len(lineage_list) < min_lineages:
        # Starting number for threshold is 1 - this will get incremented on first iteration of the while loop
        # so that in the first round of collapsing a lineage has to have at least two samples in it to not
        # be collapsed. The threshold will continue to be incremented by 1 until the min number of lineages is
        # reached.
        threshold = 1
        # Loop through this code, collapsing lineages to minimum threshold size of group and recalculating
        # % prevalence in each time period until the min number of lineages is reached
        while len(lineage_list) < min_lineages:
            # Increment threshold size of lineage group each iteration
            threshold += 1
            collapsed_list = []
            # Collapse lineages to new threshold size
            lc = LineageCollapser(
                dataframe=counts_by_period_df,
                lineages_col='lineage',
                totals_col='seq_count',
                min_level=1,
                collapsed_col='collapsed_alias',
                pango_aliases=pango_dict
                )
            lc.collapse_based_on_threshold(threshold=threshold)
            # Recalculate % prevalences in each time period with new collapsed groups
            collapsed_period_df = pd.DataFrame(lc.collapsed)[['reporting_period','collapsed_alias', 'seq_count']]
            counts_by_period_df = collapsed_period_df.groupby(['reporting_period','collapsed_alias']).sum(numeric_only=True).reset_index()
            counts_by_period_df = add_percentages_column(counts_by_period_df)
            counts_by_period_df.rename(columns={'collapsed_alias': 'lineage'}, inplace=True)
            for period in periods_to_protect:
                collapsed_list += get_lineages_above_threshold(counts_by_period_df, period, percent_threshold, "lineage")
            lineage_list = list(set(collapsed_list))
    # If more than max_lineages >5%, collapse down until have max number
    elif len(lineage_list) > max_lineages:
        # Starting number for threshold is 1 - this will get incremented on first iteration of the while loop
        # so that in the first round of collapsing a lineage has to have at least two samples in it to not
        # be collapsed. The threshold will continue to be incremented by 1 until the max number of lineages is
        # reached.
        threshold = 1
        while len(lineage_list) > max_lineages:
            # Increment threshold size of lineage group each iteration
            threshold += 1
            collapsed_list = []
            lc = LineageCollapser(
                dataframe=counts_by_period_df,
                lineages_col='lineage',
                totals_col='seq_count',
                min_level=1,
                collapsed_col='collapsed_alias',
                pango_aliases=pango_dict
                )
            lc.collapse_based_on_threshold(threshold=threshold)
            collapsed_period_df = pd.DataFrame(lc.collapsed)[['reporting_period','collapsed_alias', 'seq_count']]
            counts_by_period_df = collapsed_period_df.groupby(['reporting_period','collapsed_alias']).sum(numeric_only=True).reset_index()
            counts_by_period_df = add_percentages_column(counts_by_period_df)
            counts_by_period_df.rename(columns={'collapsed_alias': 'lineage'}, inplace=True)
            for period in periods_to_protect:
                collapsed_list += get_lineages_above_threshold(counts_by_period_df, period, percent_threshold, "lineage")
            lineage_list = list(set(collapsed_list))
    # Return lineages if max_lineage number reached without any collapsing
    else:
        lineage_list = list(set(lineage_list))
    logging.info("%s lineages identified to protect: %s", len(lineage_list), lineage_list)
    return lineage_list


def get_top_lineages_in_full_window(
    counts_by_period_df: pd.DataFrame,
    end_date:str,
    weeks_to_exclude: int,
    additional_lineages: int,
    already_protected: list
    ) -> list[str]:
    """Takes counts for reporting periods over the last year and aggregates
    them to return the top lineages by prevalence in the full reporting window,
    excluding the most recent x weeks.
    Arguments:
        counts_by_period_df -- Dataframe of per reporting period counts of lineages
        end_date -- End of the full reporting period
        weeks_to_exclude -- Number of weeks to exclude from calculations.
                            Calculated as end_date - weeks_to_exclude
        additional_lineages -- Number of additional lineages to return,
                               equivalent to top x lineages by prevalence.
        already_protected -- List of lineages already protected due to prevalence
                             in recent reporting window
    Returns:
        to_protect_list -- List of lineages to protect based on prevalence
                           across full reporting period.
    """
    # Filter out last x weeks
    date_x_weeks_ago = (end_date - dt.timedelta(weeks=weeks_to_exclude)).strftime("%Y-%m-%d")
    counts_by_period_df = counts_by_period_df[counts_by_period_df['reporting_period'] <= date_x_weeks_ago]
    logging.info(
        """Dataframe filtered to remove most recent %s weeks of data. Identifying
           the %s most prevalent lineages in this period to retain""",
        weeks_to_exclude,
        additional_lineages,
        )
    # Drop the percentage prevalence column
    counts_by_period_df = counts_by_period_df[['lineage', 'seq_count']]
    # Calculate lineage totals for full period
    totals_df = counts_by_period_df.groupby('lineage').sum(numeric_only=True).reset_index()
    period_total = totals_df["seq_count"].sum()
    # Add percentages column
    totals_df["pct_of_reporting_period"] = (totals_df["seq_count"] / period_total).mul(100)
    # Filter out lineages already protected
    totals_df = totals_df[~ totals_df['lineage'].isin(already_protected)]
    # Return additional lineages to protect by selecting top x in df
    to_protect_list = (totals_df
                       .nlargest(additional_lineages,
                                 columns="pct_of_reporting_period")
                       ['lineage']
                       .to_list()
    )
    logging.info("Identifed %s additional lineages to protect: %s",
                 len(to_protect_list),
                 to_protect_list)
    return to_protect_list

def mask_less_prevalent_values(
    counts_df: pd.DataFrame, lineages_to_leave_unmasked: dict, mask_value='Other'
) -> pd.DataFrame:
    """
    Given data frame and a dict of column names to top values (or values to not
    mask in that column) returns a copy of the data frame with non-top values
    masked with the mask value.
    Arguments:
        collapsed_counts_df -- Dataframe of lineages counts per week
        lineages_to_leave_unmasked -- Dict where key = column name in data, values
                                      are lineages that shouldn't be masked.
        mask_value -- Str to replace masked lineage values with
    Outputs:
        masked_df -- Dataframe containing masked lineages
    """
    for col in lineages_to_leave_unmasked:
        if col not in counts_df.columns:
            logging.error("Column %s to be masked not present in dataframe", col)
            sys.exit(1)
    masked_df = counts_df.copy()
    for col, lineages_no_mask in lineages_to_leave_unmasked.items():
        if isinstance(masked_df[col].dtype, pd.CategoricalDtype) and (
            mask_value not in masked_df[col].cat.categories
        ):
            masked_df[col] = masked_df[col].cat.set_categories(
                lineages_no_mask + [mask_value]
            )
        masked_df[col] = masked_df[col].mask(
            ~ masked_df[col].isin(lineages_no_mask),
            'Other'
        )
    return masked_df

def collapse_lineages_to_protected_levels(
    counts_df: pd.DataFrame,
    lineage_to_protect: list,
    pango_dict: dict
    ) -> pd.DataFrame:
    """
    """


def add_lineage_group_to_metadata(lineage_df: pd.DataFrame, counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds lineage group column to dataframe containing sample metadata. Lineage
    groups are taken from a dataframe containing counts of samples in each lineage
    per reporting period and the lineage group for each.
    Arguments:
        lineage_df -- Dataframe containing pangolin lineage info for each sample
        counts_df -- Dataframe containing counts of lineages and the lineage group
                     each lineage falls into.
    Outputs:
        lineage_df -- Updated lineage df with lineage group column addded
    """
    # Subset required columns - just need the lineage and its group. Drop
    # duplicate entries in the sub-setted table. origintating from where a
    # lineage is present in more than one reporting period in the input counts df.
    counts_df = counts_df[["lineage", "collapsed_alias"]].drop_duplicates()
    # Merge dataframes and rename collapsed alias column to lineage group
    lineage_df = (
        lineage_df
        .merge(counts_df, on="lineage")
        .rename(columns={"collapsed_alias": "lineage_group"})
    )
    return lineage_df

def write_to_csv(result_df: pd.DataFrame, outdir: Path, filename: str) -> Path:
    """
    Writes a given result dataframe to file
    Arguments:
        result_df -- Dataframe to write to file
        outdir -- Directory to write file to
        suffix -- Suffix to use in file name
    Outputs:
        result_file -- Path to result file
    """
    result_file = Path(outdir) / filename
    with Path(result_file).open("w") as file:
        result_df.to_csv(file)
    logging.info("Result dataframe written to file: %s", result_file)
    return result_file
