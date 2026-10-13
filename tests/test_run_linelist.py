import argparse
import datetime as dt
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from covid_linelist import run_linelist

ROOT_PATH: Path = Path(__file__).parents[1]
TEST_DATA_PATH: Path = Path(ROOT_PATH / "tests" / "test_data")


@pytest.fixture
def expected_year_full_lineage_metadata():
    metadata = pd.read_csv(
        Path(TEST_DATA_PATH / "2026-10-13_expected_year_full_lineage_metadata.csv")
    )
    metadata = metadata.sort_values("central_sample_id")
    return metadata


@pytest.fixture
def expected_year_lineage_masking():
    metadata = pd.read_csv(Path(TEST_DATA_PATH / "2026-10-13_expected_year_lineage_masking.csv"))
    return metadata


def test_main(tmp_path, expected_year_full_lineage_metadata, expected_year_lineage_masking, caplog):
    arguments = argparse.Namespace(
        input=Path(TEST_DATA_PATH / "covid_linelist_validation_data_v3.csv"),
        output=tmp_path,
        reporting_date=run_linelist.date_to_monday_date(
            dt.datetime.strptime("2026-04-24", "%Y-%m-%d").date()
        ),
    )
    with patch("argparse.ArgumentParser.parse_args", return_value=arguments):
        run_linelist.main()

    # Read in the created outputs:
    today = dt.date.today()

    # First compare the actual_full_lineage_metadata:

    actual_full_lineage_metadata = pd.read_csv(
        Path(tmp_path / f"{today}_covid_linelist" / f"{today}_year_full_lineage_metadata.csv")
    )
    actual_full_lineage_metadata: pd.DataFrame = actual_full_lineage_metadata.sort_values(
        "central_sample_id"
    )

    compare_full_lineage: pd.DataFrame = expected_year_full_lineage_metadata.compare(
        actual_full_lineage_metadata
    )

    if not compare_full_lineage.empty:
        print("\n")
        print("=========FAIL=========")
        print("The output created is different to the output created in October 2026 using v3.1.0.")
        print("The differenences are:")
        print(compare_full_lineage)

    else:
        print("\n")
        print("=========PASS=========")
        print("Output 'year_full_lineage_metadata' is as expected.")

    # Second compare the year_lineage_masking:

    actual_year_lineage_masking: pd.DataFrame = pd.read_csv(
        Path(tmp_path / f"{today}_covid_linelist" / f"{today}_year_lineage_masking.csv")
    )
    compare_year_lineage_masking = expected_year_lineage_masking.compare(
        actual_year_lineage_masking
    )

    if not compare_year_lineage_masking.empty:
        print("=========FAIL=========")
        print("The output created is different to the output created in October 2026 using v3.1.0.")
        print("The differenences are:")
        print(compare_year_lineage_masking)

    else:
        print("=========PASS=========")
        print("Output 'year_lineage_masking' is as expected.")

    log_msgs = [
        "1 rows of data removed from the dataframe",
        "8 lineages identified to protect",
        "RE.1.1",
        "XEC.1",
        "XEK",
        "XDV.1",
        "XFG.1",
        "JN.1",
        "BA.3.2",
        "BA.3.2.2",
        "Protecting additional lineages: ['XBB.1.5']",
    ]
    not_found = []
    for msg in log_msgs:
        if msg not in caplog.text:
            not_found.append(msg)
    if not_found:
        print("=========FAIL=========")
        print(
            "The following messages or lineages were expected to be found in the log output "
            "but were not present:"
        )
        print(not_found)
    else:
        print("=========PASS=========")
        print("The log outputs contain the messages and lineages expected.")
