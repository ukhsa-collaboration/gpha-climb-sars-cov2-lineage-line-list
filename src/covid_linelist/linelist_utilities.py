#!/usr/bin/env python3

"""
Module containing functions used in this repo
"""

# Imports - ordered (can use ruff to do this automatically)
import datetime as dt
import json
import logging
import os
from typing import List

import numpy as np
import pandas as pd
import requests
from lineage_collapser import LineageCollapser


# Functions
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
    cols_to_keep = ["sample_id", "lineage", "pangolin_version", "specimen_date"]
    input_df = pd.read_csv(input_file, 
                           usecols=lambda x: x in cols_to_keep
                           )
    return input_df
    
def check_and_update_specimen_dates(lineage_df: pd.DataFrame, specimen_date: str, fill_blanks: bool) -> pd.DataFrame:
    """Check specimen date column exists in input dataframe. If no date
    column present in the file, uses the date provided to generate the date
    column. Also optionally fills in any blank dates with today's date.
    Arguments:
        lineage_df -- Dataframe containing pangolin lineage designations
                      and sample information.
    Outputs:
        lineage_df -- Updated dataframe with checked and filled specimen date column
    """
    if 'specimen_date' not in lineage_df.columns:
        lineage_df['specimen_date'] = specimen_date
        logging.info("""No specimen date column in input dataframe. Adding
                     specimen date of %s""", specimen_date)
    elif fill_blanks:
        lineage_df['specimen_date'].fillna(specimen_date)
        logging.info("Filling in blank values in specimen date column with %s", specimen_date)

    lineage_df['specimen_date'] = pd.to_datetime(lineage_df['specimen_date'], format='%Y-%m-%d')

    return lineage_df
    

def add_week_begin_column(lineage_df: pd.DataFrame) -> pd.DataFrame:
    """Add week_begin column to lineage_df. This is required for filtering
    and grouping of data later on.
    Arguments:
        lineage_df -- Dataframe containing pangolin lineage designations
                      and sample information.
    Outputs:
        lineage_df -- Updated lineage df with week_begin column 
    """
    lineage_df['week_begin'] = lineage_df.apply(
        lambda x: x.specimen_date - pd.Timedelta(days=x.specimen_date.dayofweek)
        if pd.notnull(x.specimen_date) else np.nan, axis=1)

    return lineage_df

def filter_lineage_df_to_date_cutoff(lineage_df: pd.DataFrame, filter_start_date: str, previous_weeks_to_include: int) -> pd.DataFrame:
    """Filter the lineage_df to the required date range.
    Arguments:
        lineage_df -- Dataframe containing pangolin lineage designations
                      and sample information.
        filter_start_date -- Date to calculate cut off range from. 
        previous_weeks_to_include -- Number of weeks to include in the
                                     filtered dataframe. Rows with specimen
                                     dates above this will be filtered out.
    Outputs:
        filtered_df -- lineage_df filtered to date cut off specified
    """
    filter_start_date = dt.datetime.strptime(filter_start_date, '%Y-%m-%d')
    filter_start_date_monday = filter_start_date - dt.timedelta(days=filter_start_date.weekday())

    date_x_weeks_ago = (
        filter_start_date_monday - dt.timedelta(
            days=filter_start_date_monday.weekday(),
            weeks=previous_weeks_to_include
            )).strftime("%Y-%m-%d")
    filtered_df = lineage_df[lineage_df['specimen_date'] >= date_x_weeks_ago]
    removed_rows = (len(lineage_df) - len(filtered_df))
    logging.info("""Dataframe filtered to include samples with specimen date
                 later than %s. %d rows of data removed from the dataframe,
                 %d remaining. 
                 """, date_x_weeks_ago, removed_rows, len(filtered_df))
    return filtered_df    

def get_lineage_counts_per_week(lineage_df: pd.DataFrame) -> pd.DataFrame:
    """Create a grouped dataframe containing counts per lineage per week
    Arguments:
        lineage_df -- Dataframe containing pangolin lineage designations
                      and sample information.
    Outputs:
        counts_by_week_df -- Dataframe containing lineage counts per week
    """
    counts_by_week_df = lineage_df.groupby(
        ['week_begin', 'lineage']).size().to_frame('seq_count').reset_index()

    return counts_by_week_df

def add_timeframe_to_protect_column(counts_by_week_df: pd.DataFrame, timeframe_start: str, timeframe_length: int) -> pd.DataFrame:
    """Add boolean column specifying if week falls in date range specified
    Arguments:
        counts_by_week_df -- Dataframe containing lineage counts per week
        timeframe_end -- Most recent date to include in timeframe
        timeframe_length -- Number of weeks to include in the timeframe range
    Outputs:
        counts_by_week_df -- Updated dataframe containing lineage counts per week
    """
    # Get start and end of timeframes for boolean column
    timeframe_start = dt.datetime.strptime(timeframe_start, '%Y-%m-%d')
    timeframe_start_monday = timeframe_start - dt.timedelta(days=timeframe_start.weekday())
    timeframe_end = timeframe_start_monday - dt.timedelta(weeks=6)
    # Check if week_begin is within timeframe to protect
    counts_by_week_df[f"in_last_{timeframe_length}_weeks"] = (
        (counts_by_week_df["week_begin"] < timeframe_start) &
        (counts_by_week_df["week_begin"] > timeframe_end)
    )
    return counts_by_week_df

def add_percentages_column(counts_by_week_df: pd.DataFrame) -> pd.DataFrame:
    """Add boolean column specifying if week falls in date range specified
    Arguments:
        counts_by_week_df -- Dataframe containing lineage counts per week
        timeframe_end -- Most recent date to include in timeframe
        timeframe_length -- Number of weeks to include in the timeframe range
    Outputs:
        counts_by_week_df -- Updated dataframe containing lineage counts per week
    """
    week_total = counts_by_week_df.groupby('week_begin')["seq_count"].transform('sum')
    counts_by_week_df["pct_of_week"] = (counts_by_week_df["seq_count"] / week_total).mul(100)
    return counts_by_week_df

def get_pango_aliases(do_filter: bool=True) -> dict:
    """
    Return dict of Pango lineage alias: full lineage name from COG-UK website.
    Arguments:
        do_filter -- Bool of whether to filter alias dict
    Returns:
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

def get_lineages_to_protect(counts_by_week_df: pd.DataFrame, timeframe_length: int, percent_threshold: int, pango_dict: dict) -> List[str]:
    """Identify lineages to protect from collapsing in dataframe
    Arguments:
        counts_df -- Dataframe containing per week counts of lineages
        timeframe_length -- Timeframe to protect. Used to determine column name in_last_{x}_weeks
        percent_threshold - % of samples of a given lineage needed to protect a lineage
        pango_dict -- Dict containing pango aliases from COG-UK
    Outputs:
        to_protect_list -- List of lineages that should be protected
    """
    to_protect_collapsed = []
    weeks_to_protect = pd.to_datetime(
        counts_by_week_df[(counts_by_week_df[f"in_last_{timeframe_length}_weeks"])].week_begin.unique().tolist()
        )
    for week in weeks_to_protect:
        lc_week = LineageCollapser(counts_by_week_df[(counts_by_week_df.week_begin == week)],
                                   lineages_col='lineage',
                                   totals_col='seq_count',
                                   min_level = 2,
                                   pango_aliases=pango_dict)
        over_threshold = lc_week.collapse_based_on_pct(percent_threshold)
        to_protect_week = over_threshold.groupby('collapsed').sum('pct_of_week').reset_index().query('pct_of_week >= 1').collapsed.unique().tolist()
        to_protect_collapsed += to_protect_week

    to_protect_collapsed = set(to_protect_collapsed)
    logging.info("Identified %d lineages to protect in last %d weeks", len(to_protect_collapsed), timeframe_length)
    return list(to_protect_collapsed)

