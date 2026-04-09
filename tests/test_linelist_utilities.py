from datetime import date
from pathlib import Path

import pandas as pd
import pytest
from covid_linelist import linelist_utilities as ull

test_data_path = Path(Path(__file__).resolve().parent / "test_data")


@pytest.fixture
def test_lineage_df():
    lineage_df = pd.read_csv(
        Path(test_data_path / "test_data.csv"),
        parse_dates=["collection_date", "expected_period"],
        dayfirst=True,
    )
    lineage_df["collection_date"] = pd.to_datetime(
        lineage_df["collection_date"].dt.strftime("%Y%m%d")
    )
    return lineage_df


def test_add_reporting_period_column(test_lineage_df):
    """
    This function adds a new column that rounds the collection date to the start date of whatever reporting period is
    defined. The reporting period size is given in days, and the reporting period start dates are given that many days
    apart. So with a period size of 7 days and an end range date of Monday 30th March, the reporting dates would be
    30th, 23rd, 16th, 9th in March, and so on, for as many weeks is defined by 'previous_weeks_to_include'.
    """
    ## test data has samples collected from 02/12/2025 - 31/03/2026, a range of <20 weeks
    ## Get Monday of current week as end of range for any time periods calculated.
    end_date = date(2026, 3, 30)
    actual = ull.add_reporting_period_column(
        lineage_df=test_lineage_df,
        end_date=end_date,
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
        f'but these rows do not match:\n{actual.loc[~actual["expected_period"]==actual["reporting_period"]]}'
    )
