import json
from datetime import date
from pathlib import Path

import pandas as pd
import pytest
import requests

from covid_linelist import linelist_utilities as ull

TEST_DATA_PATH = Path(Path(__file__).resolve().parent / "test_data")

END_DATE = date(2026, 3, 30)


@pytest.fixture
def test_lineage_df():
    lineage_df = pd.read_csv(
        Path(TEST_DATA_PATH / "test_data.csv"),
        parse_dates=["collection_date", "expected_period"],
        dayfirst=True,
    )
    lineage_df["collection_date"] = pd.to_datetime(
        lineage_df["collection_date"].dt.strftime("%Y%m%d")
    )
    return lineage_df

@pytest.fixture
def counts_df_one():
    counts_df = pd.read_csv("tests/test_data/totals_df.csv")
    return counts_df

@pytest.fixture
def expected_list_one():
    lineage_list = ["RE.1.1", "XBB.1.5", "XDV.1"]
    return lineage_list

@pytest.fixture
def counts_df_two():
    counts_df = pd.read_csv("tests/test_data/totals_df_v2.csv")
    return counts_df

@pytest.fixture
def expected_list_two():
    lineage_list = ["RE.1.1", "RF.1", "RT.1", "XBB.1", "XDV.1",  "XFN"]
    return lineage_list

@pytest.fixture
def pango_dict():
    pango_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json"
    do_filter = True
    with requests.get(pango_url) as url:
        pango_aliases_dict = json.loads(url.text)
        if do_filter:
            pango_aliases_dict = dict(
                [t for t in pango_aliases_dict.items() if not t[0].startswith('X') and t[0] not in ['A', 'B']])
            return pango_aliases_dict
        else:
            return pango_aliases_dict

def test_add_reporting_period_column(test_lineage_df):
    """
    This function adds a new column that rounds the collection date to the start date of whatever reporting period is
    defined. The reporting period size is given in days, and the reporting period start dates are given that many days
    apart. So with a period size of 7 days and an end range date of Monday 30th March, the reporting dates would be
    30th, 23rd, 16th, 9th in March, and so on, for as many weeks is defined by 'previous_weeks_to_include'.
    """
    ## test data has samples collected from 02/12/2025 - 31/03/2026, a range of <20 weeks
    ## Get Monday of current week as end of range for any time periods calculated.
    actual = ull.add_reporting_period_column(
        lineage_df=test_lineage_df,
        end_date=END_DATE,
        previous_weeks_to_include=20,
        period_size=7,  # days
    )
    # Print an example to explain the test before assertion
    eg_date = "2026-03-31"
    eg = (
        actual.loc[actual["collection_date"] == eg_date, "reporting_period"]
        .dt.strftime("%Y-%m-%d")
        .values[0]
    )
    print(
        f"Test add_reporting_period_column - expect sample with collection date {eg_date} to get reporting "
        f"period 2026-03-30, got {eg}"
    )

    assert actual["reporting_period"].equals(actual["expected_period"]), (
        f'Expected the column "reporting_period" created by the function to match "expected_period" in the test data,'
        f"but these rows do not match:\n{actual.loc[~actual['expected_period'] == actual['reporting_period']]}"
    )


def test_filter_lineage_df_to_date_cutoff(test_lineage_df):
    """
    This functions filters the lineage df based on date range which is between a specified filter_end_date and the
    number of previous_weeks_to_include. Date must not be str, must be datetime object
    Noting that the function rolls the date back to the Monday of the week, then goes back the provided number of weeks.
    """
    # The date 10 weeks ago is the 19th Jan 2026, there are 12 samples with a collection date later than that.
    filtered_lineage_df = ull.filter_lineage_df_to_date_cutoff(
        lineage_df=test_lineage_df, filter_end_date=END_DATE, previous_weeks_to_include=10
    )
    assert (l := len(filtered_lineage_df)) == 12, f"Expected 12 samples after filtering, got {l}"

    # The function rolls back to Monday before finding the date to filter to the provided number of weeks previous. When
    # providing a non-Monday date in the same week, should get the same return

    filtered_lineage_df_2 = ull.filter_lineage_df_to_date_cutoff(
        lineage_df=test_lineage_df, filter_end_date=date(2026, 4, 1), previous_weeks_to_include=10
    )
    assert filtered_lineage_df.equals(filtered_lineage_df_2), (
        "Expected these dataframes to match - check that the function rolls back to Monday before filtering."
    )


def test_lineage_counts_per_period(test_lineage_df):
    """Get the counts of lineages in a reporting period"""
    test_lineage_df = test_lineage_df.rename(columns={"expected_period": "reporting_period"})

    counts = ull.get_lineage_counts_per_period(test_lineage_df)

    assert counts.iloc[12]["seq_count"] == 3
    assert counts.iloc[10]["seq_count"] == 2

    print(
        "Expected 2 counts of XFG.1.1.1 in period 2026-03-16 (row 10) and 3 counts of XFG.1.1 in "
        "period 2026-03-23 (row 12). "
        f"\ngot:\n{counts}"
    )


def test_add_percentages_column():
    counts_by_period = pd.DataFrame(
        data=[
            ["2026-03-23", "KP.2.3", 1, 10.0],
            ["2026-03-23", "XFG.1.1", 3, 30.0],
            ["2026-03-23", "XFG.1.1.1", 5, 50.0],
            ["2026-03-23", "KP.2", 1, 10.0],
            ["2026-02-09", "KP.2.3", 1, 20.0],
            [
                "2026-02-09",
                "XFG.1.1",
                1,
                20.0,
            ],
            ["2026-02-09", "XFG.1.1.1", 3, 60.0],
        ],
        columns=["reporting_period", "lineage", "seq_count", "expected_pct"],
    )

    df = ull.add_percentages_column(counts_by_period)
    assert pd.Series.equals(df["pct_of_reporting_period"], df["expected_pct"])

@pytest.mark.parametrize(
    "counts_df, pango_aliases, expected_list",
    [
        ("counts_df_one", "pango_dict", "expected_list_one"),
        ("counts_df_two", "pango_dict", "expected_list_two")
    ],
)
def test_get_top_lineages_in_full_window(counts_df, pango_aliases, expected_list, request):
    counts_df = request.getfixturevalue(counts_df)
    pango_aliases = request.getfixturevalue(pango_aliases)
    expected_list = request.getfixturevalue(expected_list)

    already_protected = ["BA.3.2", "BA.3.2.2", "RE.1.1.2", "LF.1", "XFG.1"]
    collapse_limit = ["BA.2.86", "BA.2", "BA.3", "JN.1"]
    lineage_list = ull.get_top_lineages_in_full_window(counts_df, already_protected, collapse_limit, pango_aliases, 6)
    print(lineage_list)

    assert lineage_list == expected_list
