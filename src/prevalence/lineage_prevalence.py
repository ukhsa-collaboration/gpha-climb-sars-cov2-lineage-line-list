import pandas as pd
import numpy as np
import datetime
from typing import Optional, Union
import logging
from dateutil import relativedelta as rd
import glob
import os
from src.prevalence.lineage_prevalence_class import LineageCollapser

# %% copied genomicslib.utilities functions

def in_last_time_period(
        data: pd.DataFrame, date_col: Union[str, int],
        time_period: Optional[dict] = None, back_from: Optional[str] = None
) -> pd.Series:
    """
    Given data frame with a date column, returns a pandas series with True when
    indices correspond to rows that have a date_col value in the time period
    and False otherwise. time_period is dict with a time period keyword allowed
    by the relativedelta function
    """
    time_period = time_period if time_period else dict(months=6)
    logging.debug(f'data cols:\n{data.columns}\ndate_col:{date_col}')
    if isinstance(back_from, str):
        back_from = datetime.datetime.strptime(back_from, '%Y-%m-%d').date()
    elif back_from is None:
        back_from = datetime.date.today()

    return (
        data[date_col].between(
            pd.to_datetime(back_from - rd.relativedelta(**time_period)),
            pd.to_datetime(back_from))
    )


def calculate_percents(
        data: pd.DataFrame, counts_col: str,
        groupby_cols: Optional[Union[str, list]] = None
) -> pd.Series:
    """
    Calculates the percentage of each line's counts relative to the groupings
    obtained using the groupby columns. If the groupby argument is None or an
    empty iterable, then the percentage of each line relative to the total is
    returned.

    Parameters
    ----------
    data : pd.DataFrame
        Data frame calculating a summary of a dataset where each row contains
        a column with the number of records (for example) matching the values
        for the other columns (Like 'BA.1, UK, 10').
    counts_col: str
        Name of the data frame column containing counts.
    groupby_cols : str or list, default None
        A subset of data frame column names by which the counts will be
        grouped.

    Returns
    -------
    pd.Series
        Percentage of each counts_total value relative to the total records in
        the grouping.
    """
    if groupby_cols:
        pct_series = (
            data[counts_col]
            .div(data
                 .groupby(groupby_cols)
                 [counts_col]
                 .transform('sum')
                 )
            .mul(100)
        )
    else:
        pct_series = (
            data[counts_col]
            .div(data[counts_col].sum())
            .mul(100)
        )
    return pct_series


def mask_less_prevalent_values(
        data: pd.DataFrame, col_to_top_values: dict, mask_value='Other'
) -> pd.DataFrame:
    """
    Given data frame and a dict of column names to top values (or values to not
    mask in that column) returns a copy of the data frame with non-top values
    masked with the mask value.
    """
    assert isinstance(data, pd.DataFrame), (
        'First argument to mask_less_prevalent_values should be a data frame.'
    )
    for col in col_to_top_values:
        assert col in data.columns, (
            f'Column {col} given in the col_to_top_values dict is not in the '
            f'data frame.'
        )
    masked = data.copy()
    for col, col_top_values in col_to_top_values.items():
        if isinstance(masked[col].dtype, pd.CategoricalDtype) and (
                mask_value not in masked[col].cat.categories
        ):
            masked[col] = masked[col].cat.set_categories(
                col_top_values + [mask_value]
            )
        masked[col].mask(
            ~ masked[col].isin(col_top_values),
            'Other',
            inplace=True
        )
    return masked



def return_latest_file(input_path, identifier: str) -> str:
    file = sorted(glob.glob(input_path + f"/*{identifier}"))[-1]
    # todo verify this way of getting the latest works, possibly change to using max(files, key=os.path.getmtime)
    print(f"{datetime.datetime.now()} latest {identifier}: {file}")
    return file


def get_start_of_month_date_x_weeks_ago(input_df, int_weeks: int) -> datetime:
    input_df["collection_date"] = pd.to_datetime(input_df['collection_date'], format='%Y-%m-%d')
    val_today = datetime.datetime.now()
    val_six_months_ago = (val_today - datetime.timedelta(weeks=int_weeks)).replace(day=1)
    val_six_months_ago = val_six_months_ago.strftime("%Y-%m-%d")
    print(f"{datetime.datetime.now()} start of month {int_weeks} ago: {val_six_months_ago}")
    return val_six_months_ago


def get_start_date_of_monday_x_weeks_ago(input_df, int_weeks: int):
    input_df["collection_date"] = pd.to_datetime(input_df['collection_date'], format='%Y-%m-%d')
    val_today = datetime.date.today()
    val_xweeks_ago = (val_today - datetime.timedelta(days=val_today.weekday(), weeks=int_weeks))
    val_xweeks_ago = val_xweeks_ago.strftime("%Y-%m-%d")
    print(f"{datetime.datetime.now()} start of month {int_weeks} ago: {val_xweeks_ago}")
    return val_xweeks_ago


def filter_df_only_last_12_weeks(input_df, date_column_name: str, input_date):
    temp_df = input_df[input_df[date_column_name] >= input_date]
    return temp_df


def get_date_and_lineage_info_v2(input_path, int_weeks: int, identifier: str) -> pd.DataFrame:
    file = return_latest_file(input_path, identifier)
    df = pd.read_csv(file,
                     usecols=["central_sample_id", "collection_date", "adm1", "usher_lineage", "lineage", "mutations",
                              "lineages_version", "usher_lineage", "usher_lineages_version", "lineage_conflict",
                              "lineage_ambiguity_score", "scorpio_call", "scorpio_support", "scorpio_conflict"])
    # replace Scorpio results with un-assigned
    df.loc[df['lineages_version'].str.contains("SCORPIO", na=True), 'lineage'] = "Unassigned"
    df.loc[df['lineages_version'].str.contains("SCORPIO", na=True), 'usher_lineage'] = "Unassigned"
    val_cut_off_date = get_start_of_month_date_x_weeks_ago(df, 53)
    # todo, check that this is correct. 53 weeks ago is an input
    #  here but then below the result is the input for 12 weeks prior, the function input of 13 weeks is never used
    # should the cutoff be hard coded, possibly for config file and passed down.
    df = filter_df_only_last_12_weeks(df, "collection_date", val_cut_off_date)
    df = df[df["adm1"].str.contains("UK", na=False)]
    return df


def return_week_begin_column(input_df: pd.DataFrame, date_column_name: str) -> pd.DataFrame:
    '''
    This lambda returns the date of the monday for the week input
    '''
    input_df['week_begin'] = input_df.apply(
        lambda x: x[date_column_name] - pd.Timedelta(days=x[date_column_name].dayofweek)
        if pd.notnull(x[date_column_name]) else np.nan, axis=1)
    input_df["week_begin"] = input_df["week_begin"].dt.strftime("%Y-%m-%d")
    return input_df


def generate_week_counts(input_df, int_weeks: int, date_column_name: str) -> [int, int]:
    current_week_begin = datetime.datetime.today().date() - pd.Timedelta(days=datetime.datetime.today().weekday())
    # get the date of the Monday 24 weeks ago (i.e 6 months)
    # todo check, comment here says 24 weeks, input is 53 weeks. Comment or input need edited
    # convert to pandas datetime so that time is set to 00:00:00 instead of current time of day
    # so it will pick up stuff where the date matches
    time_since = pd.to_datetime(current_week_begin - datetime.timedelta(weeks=int_weeks))
    # print(time_since)
    input_df[date_column_name] = pd.to_datetime(input_df[date_column_name])
    results_weeks = sorted(input_df[(input_df[date_column_name] >= time_since)].week_begin.unique())

    week_counts = input_df[(input_df[date_column_name].isin(results_weeks))].groupby(
        ['week_begin', 'usher_lineage']).size().to_frame('seq_count').reset_index()
    return week_counts, results_weeks


def return_modelling_mutations():
    mutations = ['R346T', 'N460K', 'K444T', 'G446S', 'L455F', 'F456L', 'F486S', 'R346I', 'K444M', 'N450D', 'V445A',
                 'K444R', 'F490S', 'F486P']
    positions = set([x[1:-1] for x in mutations])
    # positions = set([x[1:-1] for x in mutations])
    df_pos = pd.DataFrame([(int(x[1:-1]), f'S_{x[:-1]}') for x in mutations], columns=['Position', 'column_name'])
    return mutations, positions, df_pos


def generate_table_for_modelling(input_df):
    mutations, positions, df_pos = return_modelling_mutations()
    test = input_df.assign(**dict.fromkeys(df_pos["column_name"], 0))
    for pos in mutations:
        col = f'S_{pos[:-1]}'
        test[col] = test.loc[test["mutations"].str.contains(pos[:-1], na=False), col]
        test[col] = test[col].replace(0, "True")
    return test


def get_lineages_to_protect(counts, pct_threshold=1):
    """
    This method adds two columns to the table:

    - in_last_6: Boolean identifying rows in the dataframe that fall within the last 6 weeks
    - pct_of_week: float giving the percentage of the sequences within week_begin that are repesented by lineage

    The lineages are then collapsed within a single week to produce a list of all lineages that are >=1% after combining with child lineages.
    These are combined over the six week period to produce a single list of lineages to protect from the numeric collapsing.
    """
    counts['in_last_6'] = in_last_time_period(counts, 'week_begin', dict(weeks=6),
                                              datetime.date.today() - datetime.timedelta(
                                                  days=datetime.date.today().weekday()))
    counts['pct_of_week'] = calculate_percents(counts, 'seq_count', 'week_begin')

    to_protect_collapsed = list()
    for week in pd.to_datetime(counts[(counts.in_last_6 == True)].week_begin.unique().tolist()):
        lc_week = LineageCollapser(counts[(counts.week_begin == week)], 'usher_lineage', 'seq_count', min_level=2)
        over_1pct = lc_week.collapse_based_on_pct(pct_threshold)
        to_protect_week = over_1pct.groupby('collapsed').sum('pct_of_week').reset_index().query(
            'pct_of_week >= 1').collapsed.unique().tolist()
        to_protect_collapsed += to_protect_week

    to_protect_collapsed = set(to_protect_collapsed)
    return list(to_protect_collapsed)


def generate_display(row):
    if row['collapsed_alias'] in ['Unassigned', 'Other'] or row['collapsed_alias'].startswith('BA') or row[
        'lineage_clean'] == row['collapsed_alias']:
        return row['collapsed_alias']
    else:
        return f'{row["collapsed_alias"]} ({row["lineage_clean"].replace("B.1.1.529", "BA")})'


def link_climb_ids(df, majora_meta_data_path, identifier: str) -> pd.DataFrame:
    '''
    merge files on "central_sample_id" and on "anonymous_sample_id"
    '''
    file = return_latest_file(majora_meta_data_path, identifier)
    df_link = pd.read_csv(file, sep='\t', on_bad_lines='skip', usecols=["central_sample_id", "anonymous_sample_id"])
    dfdf = df.merge(df_link, left_on="central_sample_id", right_on="anonymous_sample_id")
    dfdf = dfdf.drop_duplicates(subset=['anonymous_sample_id'])
    dfdf = dfdf.rename(columns={"central_sample_id_x": "central_sample_id", "central_sample_id_y": "cog_id"})
    return dfdf


def generate_counts(results_df, weeks_to_include, prefix, local_directory, mount_point, mount_folder, min_level=2,
                    threshold=5000) -> None:
    # results_df["collection_date"] = results_df["collection_date"].dt.date
    # results_df["week_begin"] = results_df["week_begin"].dt.date
    # group by week/lineage to get counts per week
    week_counts = results_df[(results_df["week_begin"].isin(weeks_to_include))].groupby(
        ['week_begin', 'usher_lineage']).size().to_frame('seq_count').reset_index()
    # any lineages that are >=1% of the total for that week in the last 6 weeks will not be collapsed regardless of
    # total number of sequences. This will stop rapidly emerging lineages ending up in 'Other' for ages
    to_protect = get_lineages_to_protect(week_counts)

    # Use the LineageCollapser class in utilities to collapse down lineage that are not in to_protect and are <5000
    # total sequences to a maximum level of 2 - this is a change from previous because the recombinant lineage XBB is
    # not an alias and therefore doesn't have as many levels as the previous BA lineages. Those lineages are at a small
    # enough level now that they do not need to be protected.
    # Returns an updated dataframe with the collapsed lineage column and a list of lineages that are >= 5000 sequences
    # after collapsing, so they can be added to the list of lineages to protect.
    print(f"{datetime.datetime.now()} lineages to protect: {to_protect}")
    lc = LineageCollapser(week_counts, 'usher_lineage', 'seq_count', min_level=min_level, protect_lineages=to_protect,
                          collapsed_col='lineage_clean')

    week_counts_collapsed_masked = get_weeks_counts_collapsed_masked(lc, threshold, to_protect)
    # combine the alias and uncollapsed version of the lineage to make the label used in the figure legend
    week_counts_collapsed_masked['display_label'] = week_counts_collapsed_masked.apply(lambda x: generate_display(x),
                                                                                       axis=1)
    # merge the counts where the collapsed lineage and week_begin are the same to generate total counts for each
    # unique final lineage per week for the plots
    combined_counts = week_counts_collapsed_masked[['week_begin', 'display_label',
                                                    'seq_count', 'pct_of_week']].groupby(['week_begin',
                                                                                          'display_label']).sum().reset_index()

    # add unaliased lineage column to results dataframe before saving
    results_df['unaliased_lineage'] = lc.alias_to_lineage(results_df.usher_lineage)

    results_merged = get_results_merged(results_df, week_counts_collapsed_masked)

    create_output_csvs(combined_counts,
                       local_directory,
                       prefix,
                       results_df,
                       results_merged,
                       week_counts_collapsed_masked)


def get_results_merged(results_df, week_counts_collapsed_masked) -> pd.DataFrame:
    # modelling data has an older time limit to help predictions so use full results_df
    mutations = generate_table_for_modelling(results_df)
    results_merged = results_df.merge(mutations, on='central_sample_id', how='outer', suffixes=('', '_DROP')).filter(
        regex='^(?!.*_DROP)').drop_duplicates()
    collapsing = week_counts_collapsed_masked[
        ['usher_lineage', 'lineage_clean', 'collapsed_alias']].drop_duplicates().rename(
        {'lineage_clean': 'collapsed_lineage_full'})
    results_merged = results_merged.merge(collapsing, on='usher_lineage', how='left')
    return results_merged


def create_output_csvs(combined_counts, local_directory, prefix, results_df, results_merged,
                       week_counts_collapsed_masked):
    date = datetime.date.today().strftime("%Y%m%d")
    # os.chdir(local_directory)
    results_merged.to_csv(f"{local_directory}/{date}_year_full_lineage_metadata_with_mutations.csv")
    results_df.to_csv(f"{local_directory}/{date}_{prefix}_year_full_lineage_metadata.csv")
    week_counts_collapsed_masked.to_csv(f"{local_directory}/{date}_{prefix}_year_lineage_masking.csv")
    combined_counts.to_csv(f"{local_directory}/{date}_{prefix}_year_combined_counts_for_plot.csv")


def get_weeks_counts_collapsed_masked(lc, threshold, to_protect) -> pd.DataFrame:
    week_counts_collapsed = lc.collapse_based_on_threshold(threshold, convert_back=True, new_alias_col=True)
    over_threshold = week_counts_collapsed[['collapsed_alias', 'seq_count']].groupby(
        'collapsed_alias').sum().reset_index().query(f'seq_count >= {threshold}').collapsed_alias.unique().tolist()
    # combine list of lineages to protect from collapse from each of the processes
    combined_protect = set(to_protect + over_threshold + ['Unassigned'])
    # mask anything that is not in either the to_protect or over_threshold lists as Other
    week_counts_collapsed_masked = mask_less_prevalent_values(week_counts_collapsed,
                                                              dict(collapsed_alias=combined_protect))
    return week_counts_collapsed_masked


def generate_lineage_prevalence(file_path: str, file_path2: str, save_path=os.getcwd()) -> None:
    # todo the file names can probably go in the config
    df = get_date_and_lineage_info_v2(file_path, 13, "all_metadata.csv")
    df_results = return_week_begin_column(df, "collection_date")
    df_week_counts, df_results_weeks = generate_week_counts(df_results, 53, "week_begin")

    # todo can the next two lines be removed, they arent used
    df_modelling = generate_table_for_modelling(df_results)
    to_protect = get_lineages_to_protect(df_week_counts)

    last_week = (datetime.datetime.now() - datetime.timedelta(days=7))
    if df_results_weeks[-1] >= np.datetime64(last_week) or \
            df_results[(df_results.week_begin == df_results_weeks[-1])].shape[0] < 100:
        df_results_weeks = df_results_weeks[:-1]

    df_results = link_climb_ids(df_results, file_path2, "majora.metadata.tsv")
    generate_counts(results_df=df_results,
                    weeks_to_include=df_results_weeks,
                    prefix="",
                    local_directory=save_path,
                    mount_point="",
                    mount_folder="")
