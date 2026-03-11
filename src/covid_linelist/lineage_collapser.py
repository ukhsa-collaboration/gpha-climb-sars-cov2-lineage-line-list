import logging
from typing import Union

import pandas as pd


class PangoAliasError(Exception):
    pass


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
        cols_to_aggregate, protect_lineages = self.validate_inputs(collapsed_col,
                                                                   cols_to_aggregate,
                                                                   dataframe,
                                                                   lineages_col,
                                                                   min_level,
                                                                   protect_lineages,
                                                                   totals_col)

        # create a copy of the dataframe and setup a copy of the column
        # containing the lineage info. Lineage aliases are replaced with the
        # corresponding pangolin lineages so collapser works on real levels
        self.data = dataframe.copy()
        self.lineage_col = lineages_col
        self.data[collapsed_col] = self.alias_to_lineage(
            self.data[self.lineage_col].copy(), reverse=False)
        self.totals_col = totals_col
        self.collapsed = self.data.copy()

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

    def validate_inputs(self, collapsed_col, cols_to_aggregate, dataframe, lineages_col, min_level, protect_lineages,
                        totals_col):
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
            if protect_lineages.empty:
                protect_lineages = None
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
        return cols_to_aggregate, protect_lineages

    def unalias_lineage(self, lineage_converter, lineage_series):
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
            type(self).pango_aliases = self.get_pango_aliases()

        lineage_converter = (
            type(self).pango_aliases
            if not reverse
            else self.reverse_alias_dict(type(self).pango_aliases)
        )

        replaced_series = self.unalias_lineage(lineage_converter, lineage_series)

        return replaced_series

    def __get_lineage_level(self) -> pd.Series:
        """
        # todo in the original code there are two methods with this name. this one is definitely recursive.
            this one calls the other. I am not sure if there is any recursion here

        Returns series of integers matching depth of lineage name for the
        collapsed_col if collapsed is True, otherwise for the lineages_col.

        Examples
        --------
        B.1                  2
        B.1.1.7     >>       4
        AY.1                 2
        """
        return self.get_lineage_level(self.collapsed[self.collapsed_col])

    def get_lineage_level(self, lineage_series: pd.Series) -> pd.Series:
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

    def reverse_alias_dict(self, alias_dict):
        """
        When changing lineages back to alias, dict needs to be sorted by lineage length to cope with nested alias
        e.g. B.1.1.529.5.2.1.5 == BA.5.2.1.5 == BF.5
        Therefore need to replace B.1.1.529.5.2.1 in list before replacing B.1.1.529
        """
        alias_swap = {v: k for k, v in alias_dict.items()}
        sorted_swap = {}
        for k in sorted(alias_swap, key=len, reverse=True):
            sorted_swap[k] = alias_swap[k]
        return sorted_swap

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
        self.collapsed['level'] = self.__get_lineage_level()
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
                level=self.get_lineage_level(self.collapsed[self.collapsed_col]))
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
                level=self.get_lineage_level(self.collapsed[self.collapsed_col]))
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
        thresholds = self.sort_or_set_thresholds(percents, thresholds)

        local_collapsed = (
            self.data.copy()
            .assign(
                level=self.get_lineage_level(self.collapsed[self.collapsed_col]))
        )
        used_threshold = 0
        logging.debug('b4l coll', self.collapsed)
        logging.debug('b4l loc coll', local_collapsed)
        for threshold in thresholds:
            self.collapsed = (
                self.data.copy()
                .assign(
                    level=self.get_lineage_level(
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

    def sort_or_set_thresholds(self, percents, thresholds):
        if not thresholds:
            if percents:
                thresholds = list(range(5, 35, 5))
            else:
                thresholds = list(range(1_000, 11_000, 1_000))
        else:
            thresholds = sorted(list(thresholds))
        return thresholds