# Databricks notebook source
# MAGIC %md
# MAGIC # Lineage Prevalence
# MAGIC 
# MAGIC This notebook generates the data tables used for the lineage prevalence plot in the variant technical briefing and WHO COVID update slides
# MAGIC The output table needs to be processed with the R script to produce the plot because sorting out colours is a nightmare 

# COMMAND ----------

# MAGIC %pip install lxml

# COMMAND ----------

# MAGIC %md
# MAGIC ## Setup

# COMMAND ----------

import pandas as pd
import numpy as np
import datetime

from os import mkdir, chdir
from os.path import sep, exists, join

from genomicslib.utilities import in_last_time_period, calculate_percents, mask_less_prevalent_values, LineageCollapser, get_latest_line_list_filename, process_nested_variants
from genomicslib.database import get_genomics_creds, get_pheip_creds, spark_read_sql
from genomicslib.storage import mount_genomics, read_from_storage, save_to_storage

secret_scope = 'genomics-creds'
local_working_directory = join(sep, 'lineages')

genomics_mount_point = mount_genomics()

if not exists(local_working_directory):
    mkdir(local_working_directory)

chdir(local_working_directory)

genomics_creds = get_genomics_creds()
# folder to write to on storage explorer
output_folder = 'lineage-reporting'


# COMMAND ----------

# MAGIC %md
# MAGIC ## Methods

# COMMAND ----------

def get_date_and_lineage_info(creds):
    """
    Method queries genomics mart for sample date and pangolin lineage for each sequence
    for sample date, it pulls date from linkage table if available and genome table (climb) if not
    if lineage info is not available it filters the sequence out
    """
    # get lineage  and date info from genomics mart
    sql = """SELECT COALESCE(genomics.genomics_genome_table.[COG-ID], cog_uk_id) AS [COG-ID], 
    CASE WHEN genomics.genomics_linkage.Specimen_Date_SK IS NULL THEN genomics.genomics_genome_table.Sample_date ELSE genomics.genomics_linkage.Specimen_Date_SK END AS specimen_date, 
    genomics_linkage.Specimen_Number,
    genomics_linkage.cdr_specimen_request_sk,
    Adm1,
    lin.Cluster_value AS lineage
    FROM genomics.genomics_genome_table 
    FULL OUTER JOIN genomics.genomics_linkage on genomics.genomics_linkage.cog_uk_id = genomics.genomics_genome_table.[COG-ID]
    INNER JOIN (SELECT * FROM genomics.genomics_cluster_table WHERE genomics_cluster_table.Cluster_key = 'usher_lineage') lin ON lin.[COG-ID] = genomics_genome_table.[COG-ID]
    WHERE Cluster_value IS NOT NULL
    """

    results = spark_read_sql(creds, sql)
    results['specimen_date'] = pd.to_datetime(results.specimen_date, format='%Y-%m-%d')
    results['week_begin'] = results.apply(
        lambda x: x.specimen_date - pd.Timedelta(days=x.specimen_date.dayofweek)
        if pd.notnull(x.specimen_date) else np.nan, axis=1)

    return results


def generate_table_for_modelling():
    """
    Modelling need a specific list of spike mutations so that they can be modelled separately from lineages.
    They also requested the "processed" variant calls as well.
    """
    # add amino acid calls to the results table for modelling team
    mutations = ['R346T', 'N460K', 'K444T', 'G446S', 'L455F', 'F456L', 'F486S', 'R346I', 'K444M', 'N450D', 'V445A',
                 'K444R', 'F490S', 'F486P']
    positions = set([x[1:-1] for x in mutations])
    sql = f"SELECT * from phegi.genomics_proteins WHERE Position in ({', '.join(positions)}) and region_id = 16"
    mut_long = spark_read_sql(get_pheip_creds(), sql)

    pos_df = pd.DataFrame([(int(x[1:-1]), f'S_{x[:-1]}') for x in mutations], columns=['Position', 'column_name'])
    pos_df = pos_df.sort_values('Position').drop_duplicates()

    mut_long = mut_long.merge(pos_df, on='Position', how='left')
    mut_wide = mut_long[['COG-ID', 'column_name', 'Alt']].pivot(index='COG-ID', columns='column_name',
                                                                values='Alt').reset_index()

    latest_line_list = get_latest_line_list_filename(0)
    read_from_storage(genomics_mount_point, join('line-list', latest_line_list), local_working_directory)

    linelist = pd.read_csv(join(local_working_directory, latest_line_list), low_memory=False).drop(
        ["published_date", "collection_pillar", "adm1", "E484K", "K417N", "Q493R", "Sanger_Provisional"], axis=1)
    # process linelist to get "final" variant call
    linelist_long = process_nested_variants(linelist).drop_duplicates()
    # remove variant calls for anything with duplicate calls
    linelist_long = linelist_long[(~linelist_long.central_sample_id.duplicated(keep=False))]
    modelling_table = pd.merge(mut_wide, linelist_long, right_on='central_sample_id',
                               left_on='COG-ID', how='outer').drop('central_sample_id', axis=1)

    return modelling_table


def get_lineages_to_protect(counts, pct_threshold=1):
    """
    This method adds two columns to the table:

    - in_last_6: Boolean identifying rows in the dataframe that fall within the last 6 weeks
    - pct_of_week: float giving the percentage of the sequences within week_begin that are repesented by lineage

    The lineages are then collapsed within a single week to produce a list of all lineages that are >=1% after combining with child lineages.
    These are combined over the six week period to produce a single list of lineages to protect from the numeric collapsing.
    """
    counts['in_last_6'] = in_last_time_period(counts, 'week_begin', dict(weeks=6),
                                               datetime.date.today() - datetime.timedelta(days=datetime.date.today().weekday()))
    counts['pct_of_week'] = calculate_percents(counts, 'seq_count', 'week_begin')

    to_protect_collapsed = list()
    for week in pd.to_datetime(counts[(counts.in_last_6 == True)].week_begin.unique().tolist()):
        lc_week = LineageCollapser(counts[(counts.week_begin == week)], 'lineage', 'seq_count', min_level = 2)
        over_1pct = lc_week.collapse_based_on_pct(pct_threshold)
        to_protect_week = over_1pct.groupby('collapsed').sum('pct_of_week').reset_index().query('pct_of_week >= 1').collapsed.unique().tolist()
        to_protect_collapsed += to_protect_week

    to_protect_collapsed = set(to_protect_collapsed)
    return list(to_protect_collapsed)


def generate_display(row):
    if row['collapsed_alias'] in ['Unassigned', 'Other'] or row['collapsed_alias'].startswith('BA') or row[
        'lineage_clean'] == row['collapsed_alias']:
        return row['collapsed_alias']
    else:
        return f'{row["collapsed_alias"]} ({row["lineage_clean"].replace("B.1.1.529", "BA")})'


def generate_counts(results_df, weeks_to_include, prefix, local_directory, mount_point, mount_folder, min_level=2,
                    threshold=5000):
    # group by week/lineage to get counts per week
    week_counts = results_df[(results_df.week_begin.isin(weeks_to_include))].groupby(
        ['week_begin', 'lineage']).size().to_frame('seq_count').reset_index()
    # any lineages that are >=1% of the total for that week in the last 6 weeks will not be collapsed regardless of
    # total number of sequences. This will stop rapidly emerging lineages ending up in 'Other' for ages
    to_protect = get_lineages_to_protect(week_counts)

    # Use the LineageCollapser class in utilities to collapse down lineage that are not in to_protect and are <5000
    # total sequences to a maximum level of 2 - this is a change from previous because the recombinant lineage XBB is
    # not an alias and therefore doesn't have as many levels as the previous BA lineages. Those lineages are at a small
    # enough level now that they do not need to be protected.
    # Returns an updated dataframe with the collapsed lineage column and a list of lineages that are >= 5000 sequences
    # after collapsing, so they can be added to the list of lineages to protect.

    lc = LineageCollapser(week_counts, 'lineage', 'seq_count', min_level=min_level, protect_lineages=to_protect,
                          collapsed_col='lineage_clean')
    week_counts_collapsed = lc.collapse_based_on_threshold(threshold, convert_back=True, new_alias_col=True)
    over_threshold = week_counts_collapsed[['collapsed_alias', 'seq_count']].groupby(
        'collapsed_alias').sum().reset_index().query(f'seq_count >= {threshold}').collapsed_alias.unique().tolist()

    # combine list of lineages to protect from collapse from each of the processes
    combined_protect = set(to_protect + over_threshold + ['Unassigned'])

    # mask anything that is not in either the to_protect or over_threshold lists as Other
    week_counts_collapsed_masked = mask_less_prevalent_values(week_counts_collapsed,
                                                              dict(collapsed_alias=combined_protect))
    # combine the alias and uncollapsed version of the lineage to make the label used in the figure legend
    week_counts_collapsed_masked['display_label'] = week_counts_collapsed_masked.apply(lambda x: generate_display(x),
                                                                                       axis=1)
    # merge the counts where the collapsed lineage and week_begin are the same to generate total counts for each
    # unique final lineage per week for the plots
    combined_counts = week_counts_collapsed_masked[['week_begin', 'display_label',
                                                    'seq_count', 'pct_of_week']].groupby(['week_begin',
                                                                                          'display_label']).sum().reset_index()

    # add unaliased lineage column to results dataframe before saving
    results_df['unaliased_lineage'] = lc.alias_to_lineage(results_df.lineage)

    # modelling data has an older time limit to help predictions so use full results_df
    mutations = generate_table_for_modelling()
    results_merged = results_df[(results_df.specimen_date >= '2021-11-01')].merge(mutations, on='COG-ID', how='left')
    collapsing = week_counts_collapsed_masked[['lineage', 'lineage_clean', 'collapsed_alias']].drop_duplicates().rename(
        {'lineage_clean': 'collapsed_lineage_full'})
    results_merged = results_merged.merge(collapsing, on='lineage', how='left')

    # save csv files to blob store
    save_to_storage(results_merged, local_directory,
                    f'{datetime.date.today().strftime("%Y%m%d")}_full_lineage_metadata_with_mutations.csv',
                    mount_point, mount_folder)

    save_to_storage(results_df, local_directory,
                    f'{datetime.date.today().strftime("%Y%m%d")}_{prefix}full_lineage_metadata.csv', mount_point,
                    mount_folder)
    save_to_storage(week_counts_collapsed_masked, local_directory,
                    f'{datetime.date.today().strftime("%Y%m%d")}_{prefix}lineage_masking.csv', mount_point,
                    mount_folder)
    save_to_storage(combined_counts, local_directory,
                    f'{datetime.date.today().strftime("%Y%m%d")}_{prefix}combined_counts_for_plot.csv', mount_point,
                    mount_folder)


# COMMAND ----------

# MAGIC %md
# MAGIC ## Do the stuff

# COMMAND ----------

# get raw data for analysis
results = get_date_and_lineage_info(genomics_creds)

# filter all the old stuff because it makes the plot too big - this will be adjusted over time
# get the date of the Monday in the current week
current_week_begin = datetime.datetime.today().date() - pd.Timedelta(days=datetime.datetime.today().weekday())
# get the date of the Monday 24 weeks ago (i.e 6 months)
# convert to pandas datetime so that time is set to 00:00:00 instead of current time of day
# so it will pick up stuff where the date matches
time_since = pd.to_datetime(current_week_begin - datetime.timedelta(weeks=24))
# get all the week_beginning values that are after (and including) the time since date
results_weeks = sorted(results[(results.week_begin >= time_since)].week_begin.unique())

# remove the last week of results if it is within 7 days as the numbers are too small and
# increase the number of lineages in the plot because of the 1% threshold
last_week = (datetime.datetime.now() - datetime.timedelta(days=7))
if results_weeks[-1] >= np.datetime64(last_week) or results[(results.week_begin == results_weeks[-1])].shape[0] < 100:
    results_weeks = results_weeks[:-1]
# generate counts and save files
generate_counts(results, results_weeks, '', local_working_directory, genomics_mount_point, output_folder)

# COMMAND ----------
