import json
import requests
import base64
import datetime
from dateutil import relativedelta as rd
import IPython
import logging
import pandas as pd
from typing import Optional, Union
import numpy as np
from os.path import join
from genomicslib.storage import read_from_storage
from genomicslib.database import spark_read_sql, get_genomics_creds
from genomicslib.linelist import get_variant_yaml_git, load_yaml_paths_glob_object_into_dict
import uuid

def get_latest_line_list_filename(
    days_delta=0, *, ll_location="/mnt/phe/genomics/line-list/", ll_suffix="epi_lne_list.csv", filename_len=32):
    dbutils = IPython.get_ipython().user_ns["dbutils"]
    delta_date = (datetime.datetime.today() - datetime.timedelta(days=days_delta)).strftime("%Y%m%d")
    delta_date_filename = list()
    for fileobject in dbutils.fs.ls(ll_location):
        if fileobject.name.startswith(delta_date) and fileobject.name.endswith(ll_suffix) and len(
                fileobject.name) == filename_len:
            delta_date_filename.append(fileobject.name)
    if len(delta_date_filename) > 0:
        return sorted(delta_date_filename)[-1]
    else:
        print("No linelist found with date", delta_date, "trying 1 day earlier")
        days_delta += 1
        return get_latest_line_list_filename(
            days_delta=days_delta, ll_location=ll_location, ll_suffix=ll_suffix, filename_len=filename_len)


def get_latest_lowqc_review_filename(days_delta=0):
    dbutils = IPython.get_ipython().user_ns["dbutils"]
    delta_date = (datetime.datetime.today() - datetime.timedelta(days=days_delta)).strftime("%Y%m%d")
    delta_date_filename = list()
    for fileobject in dbutils.fs.ls("/mnt/phe/genomics/line-list/"):
        if (delta_date in fileobject.name) and fileobject.name.endswith(".csv") and len(
                fileobject.name) == 32 and fileobject.name.startswith("lowqc_review_"):
            delta_date_filename.append(fileobject.name)
    if len(delta_date_filename) > 0:
        return sorted(delta_date_filename)[-1]
    else:
        print("No lowqc_review found with date", delta_date, "trying 1 day earlier")
        days_delta += 1
        return get_latest_lowqc_review_filename(days_delta)


def get_lineage_level(lineage_series: pd.Series) -> pd.Series:
    """
    Given pandas series with lineages, returns series of integers matching
    depth of lineage name.

    Examples
    --------
    B.1                  2
    B.1.1.7     >>       4
    AY.1                 2
    """
    return (
        lineage_series
        .str.split('.')
        .map(lambda x: len(x if isinstance(x, list) else [x]))  # deal w/ NaNs
    )


class PangoAliasError(Exception):
    pass


def get_pango_aliases(do_filter=True) -> dict:
    """
    Return dict of Pango lineage alias: full lineage name from COG-UK website.
    """
    try:
        with requests.get("https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json") as url:
            data = json.loads(url.text)
            if do_filter:
                pango_aliases_filter = dict([t for t in data.items() if not t[0].startswith('X') and t[0] not in ['A', 'B']])
                return pango_aliases_filter
            else:
                return data
    except:
        raise PangoAliasError('Cannot generate alias dict from json file')


def reverse_alias_dict(alias_dict):
    """
    When changing lineages back to alias, dict needs to be sorted by lineage length to cope with nested alias
    e.g. B.1.1.529.5.2.1.5 == BA.5.2.1.5 == BF.5
    Therefore need to replace B.1.1.529.5.2.1 in list before replacing B.1.1.529
    """
    alias_swap = {v: k for k, v in alias_dict.items()}
    sorted_swap = dict()
    for k in sorted(alias_swap, key=len, reverse=True):
        sorted_swap[k] = alias_swap[k]
    return sorted_swap


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


def unalias_lineage(lineage_converter, lineage_series):
    replaced_series = pd.Series(lineage_series.copy())
    for to_repl, repl_with in lineage_converter.items():
        # match lineage exact match and sub-lineages
        lin_match = (
            replaced_series.str.startswith(f'{to_repl}.') &
            (replaced_series != to_repl)
        )
        # replace matched sub-string by replacement sub-string (avoiding
        # slower regex substitution)
        replaced_series.loc[lin_match] = (
            repl_with + replaced_series.loc[lin_match].str[len(to_repl):]
        )
    return replaced_series


class LineageCollapser:
    # this is set to None so that importing the lib is not prevented if there is a problem with generating the alias
    # dict (e.g. web page availability) and will only produce an error when the dict is actually required
    pango_aliases = None  # alias: lineage

    def __init__(
        self, dataframe: pd.DataFrame, lineages_col: Union[str, int],
        totals_col: Union[str, int], min_level: int = 2,
        cols_to_aggregate: Union[str, list, None] = None,
        protect_lineages: Union[list, tuple, set, pd.Series, None
        ] = None, collapsed_col: str = 'collapsed'
    ) -> None:
        assert isinstance(dataframe, pd.DataFrame), (
            'dataframe should be a pandas data frame.')
        assert isinstance(lineages_col, (str, int)), (
            f'lineage_column_name should be a str or int. {type(lineages_col)}'
            f' given.'
        )
        dataframe_cols = ', '.join(dataframe.columns.to_list())
        assert lineages_col in dataframe.columns, (
            f'{lineages_col} is not a column in the dataframe provided.'
            f' columns in dataframe: {dataframe_cols}.'
        )
        assert isinstance(totals_col, (str, int)), (
            f'totals_column_name should be a str or int. '
            f'{type(totals_col)} given')
        assert totals_col in dataframe.columns, (
            f'{totals_col} is not a column in the dataframe provided.'
            f' columns in dataframe: {dataframe_cols}.'
        )
        assert isinstance(min_level, int), (
            f'min_level should be an integer. {type(min_level)} given.')
        if isinstance(cols_to_aggregate, str):
            cols_to_aggregate = [cols_to_aggregate]
        elif cols_to_aggregate is None:
            cols_to_aggregate = []
        else:
            cols_to_aggregate = cols_to_aggregate
        assert isinstance(cols_to_aggregate, list), (
            'Extra columns should be None, str or list.'
        )
        for col in cols_to_aggregate:
            assert isinstance(col, (str, int)), (
                f'{col} should be a str or int. {type(collapsed_col)} '
                f'given.'
            )
            assert col in dataframe.columns, (
                f'{col} is not a column in the dataframe provided.'
                f' columns in dataframe: {dataframe_cols}.'
            )
        if protect_lineages is not None:
            protect_lineages = pd.Series(protect_lineages, dtype='object')
        assert (
            protect_lineages is None
            or (isinstance(protect_lineages, pd.Series)
                and not protect_lineages.empty)), (
            f'protect_lineages should be None, list, tuple, set or pandas '
            f'series. {type(protect_lineages)} given: {protect_lineages}'
        )
        assert isinstance(collapsed_col, (str, int)), (
            f'collapsed_col should be a str or int. {type(collapsed_col)} '
            f'given.'
        )
        assert collapsed_col not in dataframe.columns, (
            f'{collapsed_col} should not be a column in the dataframe provided'
            f' as it will be created to store the lineage each record will be '
            f'collapsed into.'
        )

        # create a copy of the dataframe and setup a copy of the column
        # containing the lineage info. Lineage aliases are replaced with the
        # corresponding pangolin lineages so collapser works on real levels
        self.data = dataframe.copy()
        self.lineage_col = lineages_col
        self.data[collapsed_col] = self.alias_to_lineage(
            self.data[self.lineage_col].copy(), reverse=False)
        self.totals_col = totals_col
        self.min_level = min_level
        if cols_to_aggregate:
            self.cols_to_aggregate = [
                col for col in cols_to_aggregate
                if col not in (lineages_col, collapsed_col)
            ]
        else:
            self.cols_to_aggregate = []
        if protect_lineages is not None:
            protect_lineages = self.alias_to_lineage(
                protect_lineages, reverse=False).to_list()
            self.protect_lineages = protect_lineages
        else:
            self.protect_lineages = []
        self.collapsed_col = collapsed_col
        self.collapsed = self.data.copy()

    def alias_to_lineage(self, lineage_series, reverse: bool = False) -> pd.Series:
        """
        Given pandas series of string objects and dict of sub-strings to replace
        and replacement sub-strings, returns series where rows where strings
        starting with the dict key sub-string are replaced by the string starting
        with the `start_strs_to_replace` dict value.

        Examples
        --------
        B.1.1               B.1.1           (B.1 not an alias)
        AY.4.1      >>      B.1.617.2.4.1   (AY.4 alias of B.1.617.2.4)
        BA.1                B.1.1.529.1     (BA.1 alias of B.1.1.529.1)
        """
        if not type(self).pango_aliases:
            type(self).pango_aliases = get_pango_aliases()

        lineage_converter = (
            type(self).pango_aliases
            if not reverse
            else reverse_alias_dict(type(self).pango_aliases)
        )

        replaced_series = unalias_lineage(lineage_converter, lineage_series)

        return replaced_series

    def get_lineage_level(self) -> pd.Series:
        """
        Returns series of integers matching depth of lineage name for the
        collapsed_col if collapsed is True, otherwise for the lineages_col.

        Examples
        --------
        B.1                  2
        B.1.1.7     >>       4
        AY.1                 2
        """
        return get_lineage_level(self.collapsed[self.collapsed_col])

    def __get_lineages_to_collapse(
        self, level: int, threshold: Union[int, pd.Series]
    ) -> pd.Series:
        """
        Return boolean pandas series to filter data by based on lineage level
        and a count threshold.
        """
        # only collapse lineages with fewer sequences than the threshold
        # get totals for cols to group
        below_threshold = threshold > (
            self.collapsed
            .groupby(self.cols_to_aggregate + [self.collapsed_col])
            [self.totals_col]
            .transform('sum')
        )
        logging.debug('below', below_threshold, '\tthr', threshold)
        # select records with lineage level being collapsed
        in_level = self.collapsed['level'] == level

        if self.protect_lineages:
            protected = (
                self.collapsed[self.collapsed_col].isin(self.protect_lineages)
            )
        else:
            protected = pd.Series(
                [False] * self.collapsed.shape[0], index=self.collapsed.index
            )
        return below_threshold & in_level & ~ protected

    def __remove_highest_lineage_level(self, indices) -> pd.Series:
        """
        Get lineage parent for each lineage.

        Examples
        --------
        'B.1.1'     >>      'B.1'
        """
        return (
            self.collapsed
            .loc[indices, self.collapsed_col]
            .str.rsplit('.', n=1)
            .str[0]
        )

    def __collapse_lineages(
        self, level: int, threshold: Union[int, pd.Series]
    ) -> pd.DataFrame:
        to_collapse = self.__get_lineages_to_collapse(
            level, threshold=threshold
        )
        self.collapsed.loc[
            to_collapse, self.collapsed_col
        ] = self.__remove_highest_lineage_level(to_collapse)
        self.collapsed['level'] = self.get_lineage_level()
        return self.collapsed

    def __convert_back_to_alias_and_drop_level(self, new_col: bool = False) -> pd.DataFrame:
        # convert back to aliases
        col_name = (self.collapsed_col
                    if not new_col
                    else 'collapsed_alias')
        self.collapsed[col_name] = self.alias_to_lineage(
            self.collapsed[self.collapsed_col], reverse=True
        )
        # no further need for level column
        self.collapsed.drop(columns='level', inplace=True)
        return self.collapsed

    def collapse_based_on_threshold(
        self, threshold: Union[int, pd.Series], convert_back: bool = True, new_alias_col: bool = False
    ) -> pd.DataFrame:
        # process each level starting from highest to avoid collapsing a lower
        # level lineage that would accumulate enough sequences from higher levels:
        # e.g.                          should result in         and not
        # B.1       10                   B.1       10          B.1       11
        # B.1.1     1           >>       B.1.1     11          B.1.1     10
        # B.1.1.4   4
        # B.1.1.5   6                   (threshold of 10)
        self.collapsed = (
            self.data.copy()
            .assign(
                level=get_lineage_level(self.collapsed[self.collapsed_col]))
        )
        for level in range(self.collapsed['level'].max(), self.min_level, -1):
            self.collapsed = self.__collapse_lineages(level, threshold)
        if convert_back:
            self.__convert_back_to_alias_and_drop_level(new_col=new_alias_col)
        return self.collapsed

    def __get_thresholds_based_on_pct(
        self, records_pct: Union[int, float]
    ) -> pd.Series:
        if self.cols_to_aggregate:
            thresholds = (
                self.collapsed
                .groupby(self.cols_to_aggregate)
                [self.totals_col]
                .transform('sum')
                .mul(records_pct / 100)
            ).round()
        else:
            thresholds = round(
                self.collapsed[self.totals_col].sum() * (records_pct / 100)
            )
        return thresholds

    def collapse_based_on_pct(
        self, records_pct: Union[int, float], convert_back: bool = True
    ) -> pd.DataFrame:
        self.collapsed = (
            self.data.copy()
            .assign(
                level=get_lineage_level(self.collapsed[self.collapsed_col]))
        )
        threshold = self.__get_thresholds_based_on_pct(records_pct)
        logging.debug(
            'pct', records_pct, 'collapse_based_on_pct thr on pct', threshold
        )
        for level in range(self.collapsed['level'].max(), self.min_level, -1):
            logging.debug('level', level)
            self.collapsed = self.__collapse_lineages(level, threshold)
            logging.debug('coll on pct coll\n', self.collapsed)
        if convert_back:
            self.__convert_back_to_alias_and_drop_level()
        logging.debug('334', self.collapsed)
        return self.collapsed

    def collapse_recursively_to_at_least_n(
        self, n: int, thresholds: Union[list, tuple, set, pd.Series, None
        ] = None, percents: bool = True
    ) -> tuple:
        if not thresholds:
            if percents:
                thresholds = list(range(5, 35, 5))
            else:
                thresholds = list(range(1_000, 11_000, 1_000))
        else:
            thresholds = sorted(list(thresholds))

        local_collapsed = (
            self.data.copy()
            .assign(
                level=get_lineage_level(self.collapsed[self.collapsed_col]))
        )
        used_threshold = 0
        logging.debug('b4l coll', self.collapsed)
        logging.debug('b4l loc coll', local_collapsed)
        for threshold in thresholds:
            self.collapsed = (
                self.data.copy()
                .assign(
                    level=get_lineage_level(
                        self.collapsed[self.collapsed_col])
                )
            )
            logging.debug('\ncollapse_recursively_to_at_least_n: threshold=', threshold)
            if percents:
                self.collapse_based_on_pct(threshold, convert_back=False)
            else:
                self.collapse_based_on_threshold(threshold, convert_back=False)
            if self.collapsed[self.collapsed_col].nunique() > n:
                local_collapsed = self.collapsed.copy()
                used_threshold = threshold
            elif self.collapsed[self.collapsed_col].nunique() == n:
                used_threshold = threshold
                break
            else:
                self.collapsed = local_collapsed
                logging.debug('collapse n: inst coll=', self.collapsed)
                logging.debug('collapse n: lc=', local_collapsed)
                break
        logging.debug('collapse_recursively_to_at_least_n ut', used_threshold)
        logging.debug(f'return coll to n b4 tidy\n {self.collapsed}')
        return used_threshold, self.__convert_back_to_alias_and_drop_level()


def get_top_values(
    data: pd.DataFrame, values_to_pick: Union[str, list], counts_column: str,
    top_n: int
) -> list:
    """
    Returns list of top_n values_to_pick (a single column or a set of columns)
    in the dataframe based on the sum of counts_column.
    """
    return (
        data
        .groupby(values_to_pick)
        [counts_column]
        .sum()
        .nlargest(top_n)
        .index
        .to_list()
    )


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


def df_to_fasta(dataframe, ids_column, sequences_column, out_fas_path):
    """Take data frame with sequences and save fasta file."""
    with open(out_fas_path, 'w') as fas:
        fas.write(
            '\n'.join(
                ('>' + dataframe[ids_column] + '\n'
                 + dataframe[sequences_column])
                .to_list()
            )
        )


def setup_logging(log_filepath, task_name, level=logging.INFO):
    """
    Set up logging to the stdout and a log file at level given. Turns java
    logging level to error only to avoid verbose java logging.
    """
    logging.basicConfig(
        handlers=[logging.FileHandler(log_filepath),
                  logging.StreamHandler()],
        level=level,
        format=f'%(asctime)s : %(levelname)s : {task_name} : %(message)s')
    logging.getLogger("py4j.java_gateway").setLevel(logging.ERROR)


def sort_genomics_variant_call(variant_call):
    """
    Key function to sort variants list.
    e.g. sorted(list_of_variant_calls, key=sort_genomics_variant_call)
    """
    month_num = {
        'JAN': 1, 'FEB': 2, 'MAR': 3, 'APR': 4, 'MAY': 5, 'JUN': 6,
        'JUL': 7, 'AUG': 8, 'SEP': 9, 'OCT': 10, 'NOV': 11, 'DEC': 12
    }
    var_split = variant_call.split('-')
    if len(var_split) > 2:
        voc_vui = var_split[0]
        yy = var_split[1][:2]
        month = var_split[1][2:]
        num = var_split[2]
        var_split = [yy, month_num[month], num, voc_vui]
    elif variant_call == 'Omicron_Unassigned':
        var_split = ['21', 11, '00', 'REQ']
    elif variant_call.startswith('SIM'):
        var_split = ['_', '_', '00', 'SIM']
    else:
        var_split = ['_'] * 3 + var_split
    return var_split


def get_tree_file_dict_from_cog_id_list(list_of_cogids, path_to_tree_notebook="gen_tree_from_query_plus_bg_seqs"):
    """
    Takes a list of COG-IDs as input, runs the gen_tree_from_query_plus_bg_seqs notebook
    Returns a dictionary of the following format:
    
    {
        "png": "/17b2d21b-8a7f-413b-95e6-7d920dbb9574/out.png", 
        "pdf": "/17b2d21b-8a7f-413b-95e6-7d920dbb9574/out.pdf", 
        "ctree": "/17b2d21b-8a7f-413b-95e6-7d920dbb9574/out.ctree"
    }
    
    i.e. filetypes as keys and paths as values
    directories are a psuedorandom UUID
    """
    list_of_cogids = json.dumps(list_of_cogids)
    return_value = dbutils.notebook.run(path_to_tree_notebook, 0, {"id_list": list_of_cogids})
    return_value = json.loads(return_value)
    return return_value


def view_png(path_to_png_file, method="HTML"):
    """
    Displays a PNG file within a databricks notebook using the method specified.
    
    Currently only supports HTML, but in theory could support display via
    matplotlib, plotly, or other mechanisms
    """
    if method == "HTML":
        with open(path_to_png_file, "rb") as display_file:
            encoded_string_png = base64.b64encode(display_file.read()).decode('utf-8')
        
        return displayHTML(
            f"""<div style='overflow-y: scroll; height: 750px;'><img style='height: 750px; width: auto' src='data:image/png;base64, {encoded_string_png}' /></div>"""
        )
    else:
        raise Exception("Unknown method")


def get_list_of_required(yml, yaml_dict):
    requires = list()
    if 'requires' in yml.keys():
        parent = yml['requires']
        requires.append(yaml_dict[parent]['phe-label'])
        while 'requires' in yaml_dict[parent].keys():
            parent = yaml_dict[parent]['requires']
            requires.append(yaml_dict[parent]['phe-label'])
    return requires


def process_nested_variants(linelist_var):
    """
    This method takes in a linelist and produces a table with a single variant call for each sequence.
    Provide one with only central_sample_id and variant cols (i.e. no mutations or metadata such as Adm1)

    For seqs that are called as two variants that are not nested (e.g. BA.4 and BA.5) they will still be in the table
    twice and you need to deal with those as you see fit.
    """
    # replace low qc variant calls and calls that are redundant (e.g. omicron unassigned + BA.2)
    linelist_var = linelist_var.replace('Low_qc', np.nan)
    # remove alpha call for seqs called as V-21FEB-02
    if linelist_var[(linelist_var['VOC-21FEB-02'].notnull())].shape[0] > 0:
        linelist_var.loc[
            (linelist_var["VOC-21FEB-02"].isin(["Confirmed", "Probable", "Low_qc_high"])), "VOC-20DEC-01"] = np.nan
    # remove delta call for seqs called as AY.4.2
    if linelist_var[(linelist_var['VUI-21OCT-01'].notnull())].shape[0] > 0:
        linelist_var.loc[(linelist_var["VUI-21OCT-01"].isin(["Confirmed", "Probable"])), "VOC-21APR-02"] = np.nan

    # remove BA.2 call for BA.2.75 seqs because they are not nested variants
    if linelist_var[(linelist_var['V-22JUL-01'].notnull())].shape[0] > 0:
        linelist_var.loc[(linelist_var['V-22JUL-01'].isin(["Confirmed", "Probable", "Low_qc_high"])),
                         'VUI-22JAN-01'] = np.nan

    # remove Omicron_Unassigned call for XE seqs because they are not nested variants
    if linelist_var[(linelist_var['V-22APR-02'].notnull())].shape[0] > 0:
        linelist_var.loc[(linelist_var['V-22APR-02'].isin(["Confirmed", "Probable", "Low_qc_high"])),
                         'Omicron_Unassigned'] = np.nan

    # get all the yaml files from the git repo and remove all the calls for "required" definitions
    # this method removes the need for this bit of code to be updated every time a new variant is declared
    # if the relationships are described within the yaml files
    folder_id = uuid.uuid4()
    yamls = get_variant_yaml_git(f'/tmp/variant_definitions_{folder_id}')
    variant_yaml_dict = load_yaml_paths_glob_object_into_dict(yamls)
    for yml_id in variant_yaml_dict.keys():
        phe_label = variant_yaml_dict[yml_id]['phe-label']
        requires = get_list_of_required(variant_yaml_dict[yml_id], variant_yaml_dict)
        if linelist_var[(linelist_var[phe_label].notnull())].shape[0] > 0 and len(requires) > 0:
            linelist_var.loc[(linelist_var[phe_label].isin(["Confirmed", "Probable", "Low_qc_high"])),
                             requires] = np.nan

    # get list of unclassified sequences based on no result in edited table
    linelist_var["Unclassified"] = linelist_var.drop(["central_sample_id"], axis=1).isnull().values.all(axis=1)
    # switch to long format
    linelist_var_long = linelist_var.melt(id_vars=["central_sample_id"], var_name="variant", value_name="result",
                                          ignore_index=True)
    # replace True with "Confirmed" for Unclassified seqs
    linelist_var_long.loc[
        (linelist_var_long.variant == 'Unclassified') & (linelist_var_long.result == True), 'result'] = 'Confirmed'
    # filter to remove empty rows or unreviewed low_qc
    linelist_var_long = linelist_var_long[(linelist_var_long.result.isin(["Confirmed", "Probable", "Low_qc_high"]))]

    return linelist_var_long


def get_seq_metadata(date_cutoff):
    # get sequence metadata
    # samples that are in linkage will use specimen_date_sk, rest will use date submitted to climb
    sql = f"""SELECT * FROM 
              (SELECT genomics_genome_table.[COG-ID] AS central_sample_id, 
                      CASE WHEN Specimen_Date_SK IS NULL THEN Sample_date ELSE Specimen_Date_SK END AS specimen_date,
                      CASE WHEN Specimen_Date_SK IS NULL THEN 'CLIMB' ELSE 'Epilink2' END AS date_source, 
                      Published_date,
                      cdr_specimen_request_sk, finalid, Specimen_Number, 
                      Pillar AS collection_pillar, 
                      Adm1,
                      usher_lineage,
                      seq_org_code
              FROM genomics.genomics_genome_table 
              LEFT JOIN genomics.genomics_linkage ON genomics_linkage.cog_uk_id = genomics_genome_table.[COG-ID] COLLATE Latin1_General_100_CS_AS
              LEFT JOIN (SELECT [COG-ID], Cluster_value AS usher_lineage 
                         FROM genomics.genomics_cluster_table WHERE Cluster_key = 'usher_lineage') u 
              ON u.[COG-ID] = genomics_genome_table.[COG-ID]
              LEFT JOIN (SELECT [COG-ID], Cluster_value AS seq_org_code 
                         FROM genomics.genomics_cluster_table WHERE Cluster_key = 'sequencing_org_code') s 
              ON s.[COG-ID] = genomics_genome_table.[COG-ID]) t """

    if date_cutoff is not None:
        sql += "WHERE t.specimen_date >= '{date_cutoff}'"
    metadata = spark_read_sql(get_genomics_creds(), sql)
    metadata['unaliased_lineage'] = unalias_lineage(get_pango_aliases(), metadata.usher_lineage)

    return metadata


def get_sequence_metadata_table(mount_point, local_directory, date_cutoff='2022-01-03'):
    metadata = get_seq_metadata(date_cutoff)

    # get most recent linelist
    latest_line_list = get_latest_line_list_filename(0)

    read_from_storage(mount_point, join('line-list', latest_line_list), local_directory)

    linelist = pd.read_csv(join(local_directory, latest_line_list), low_memory=False)

    # filter to only include IDs that are also in the metadata column
    linelist = linelist[(linelist.central_sample_id.isin(metadata.central_sample_id.values.tolist()))]
    # drop columns that are not needed
    linelist_var = linelist.copy().drop(["published_date", "collection_pillar", "adm1", "E484K", "K417N", "Q493R",
                                         "Sanger_Provisional"], axis=1)
    linelist_long = process_nested_variants(linelist_var)

    merged = linelist_long.merge(metadata, on='central_sample_id', how='inner')
    merged['is_duplicate'] = merged.central_sample_id.duplicated(keep=False)

    return merged
