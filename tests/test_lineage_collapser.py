import json

import pandas as pd
import pytest
import requests

from covid_linelist.lineage_collapser import LineageCollapser


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

@pytest.fixture
def test_data():
    return pd.DataFrame(
        {
            "lineage": [
                "B.1.1",  # basic lineage
                "B.1.1.7",  # longer lineage
                "B.1.1.7.1",  # longer still
                "B.2",  # different base level
                "B.2.1",  # longer
                "XEC",  # recombinant base level
                "XEC.1",  # recombinant to stay
                "XEC.2",  # recombinant that should collapse
                "XEC.1.1.1",  # recombinant that should collapse
                "XEC.2.1.7.15",  # recombinant that should collapse
            ],
            # !!counts are rounded in the pct calculation!!
            "count": [
                10,  # B.1.1
                5,  # B.1.1.7
                3,  # B.1.1.7.1
                20,  # B.2
                3,  # B.2.1
                30,  # XEC - base
                15,  # XEC.1 - high count, should stay
                2,  # XEC.2 - low count, should collapse
                2,  # XEC.1.1.1 - low count, should collapse
                2,  # XEC.2.1.7.15 - low count, should collapse
            ],
        }
    )


@pytest.fixture
def collapser(test_data, pango_dict):
    return LineageCollapser(
        dataframe=test_data,
        lineages_col="lineage",
        totals_col="count",
        min_level=1,
        collapsed_col="collapsed",
        pango_aliases=pango_dict,
    )


def test_protect_recombinant(test_data, pango_dict):
    # if we protect XEC.2, it should remain when count is low
    protected = ["XEC.2"]
    lc = LineageCollapser(
        dataframe=test_data,
        lineages_col="lineage",
        totals_col="count",
        protect_lineages=protected,
        pango_aliases=pango_dict,
    )

    result = lc.collapse_based_on_threshold(threshold=5)
    assert "XEC.2" in result["collapsed"].values


def test_recombinant_collapse(collapser):
    # collapse anything under 5 which is:
    # B.1.1.7.1
    # B.2.1
    # XEC.2
    # XEC.1.1.1
    # XEC.2.1.7.15

    result = collapser.collapse_based_on_threshold(threshold=5)
    # XEC.2 should be collapsed to XEC
    assert "XEC.2" not in result["collapsed"].values
    # XEC.2.1.7.15 should be collapsed to XEC
    assert "XEC.2.1.7.15" not in result["collapsed"].values
    # XEC.1.1.1 should be collapsed to XEC
    assert "XEC.1.1.1" not in result["collapsed"].values
    # XEC.1 should stay
    assert "XEC.1" in result["collapsed"].values
    # Base XEC should also stay
    assert "XEC" in result["collapsed"].values


def test_percentage_with_recombinants(collapser):
    # total count is 90, so 5% is 4.5
    result = collapser.collapse_based_on_pct(records_pct=5)

    # XEC.2 should be collapsed to XEC
    assert "XEC.2" not in result["collapsed"].values
    # XEC.2.1.7.15 should be collapsed to XEC
    assert "XEC.2.1.7.15" not in result["collapsed"].values
    # XEC.1.1.1 should be collapsed to XEC.1
    assert "XEC.1.1.1" not in result["collapsed"].values
    # XEC.1 should stay
    assert "XEC.1" in result["collapsed"].values
    # XEC should stay
    assert "XEC" in result["collapsed"].values


def test_standard_lineages_collapse(collapser):
    result = collapser.collapse_based_on_threshold(threshold=5)

    # B.1.1 should stay
    assert "B.1.1" in result["collapsed"].values
    # B.1.1.7 should stay
    assert "B.1.1.7" in result["collapsed"].values
    # B.1.1.7.1 should collapse into B.1.1.7
    assert "B.1.1.7.1" not in result["collapsed"].values
    # B.2 should stay
    assert "B.2" in result["collapsed"].values
    # B.2.1 should collapse into B.2
    assert "B.2.1" not in result["collapsed"].values


def test_standard_percentage_collapse(collapser):
    result = collapser.collapse_based_on_pct(records_pct=5)

    # B.1.1 should stay
    assert "B.1.1" in result["collapsed"].values
    # B.1.1.7 should stay
    assert "B.1.1.7" in result["collapsed"].values
    # B.1.7.1 should collapse into B.1.1.7
    assert "B.1.1.7.1" not in result["collapsed"].values
    # B.2 should stay
    assert "B.2" in result["collapsed"].values
    # B.2.1 should collapse into B.2
    assert "B.2.1" not in result["collapsed"].values
