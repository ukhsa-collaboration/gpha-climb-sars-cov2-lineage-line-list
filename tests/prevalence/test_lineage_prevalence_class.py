import pytest
from src.prevalence.lineage_prevalence_class import LineageCollapser
from unittest.mock import patch
import pandas as pd


@pytest.fixture
def simple_aliases():
    """ A tiny deterministic alias map for testing. AY.4 → B.1.617.2.4 BA.1 → B.1.1.529.1 """
    return { "AY.4": "B.1.617.2.4", "BA.1": "B.1.1.529.1", }


@pytest.fixture
def df_basic():
    """ Simple DataFrame with a few lineages and counts. """
    return pd.DataFrame({
                            "lineage": ["B.1", "B.1.1", "B.1.1.4", "AY.4.1", "BA.1"],
                            "count": [10, 1, 4, 3, 2],
                            "region": ["X", "X", "X", "X", "X"],
                        })


# @pytest.fixture
# def setup_lineage_collapser():
#     min_level = 2
#     threshold = 5000
#     yield  LineageCollapser(week_counts, 'usher_lineage', 'seq_count', min_level=min_level, protect_lineages=to_protect,
#                           collapsed_col='lineage_clean')


def test_alias_to_lineage_expands_aliases(df_basic, simple_aliases):
    with patch.object(LineageCollapser, "get_pango_aliases", return_value=simple_aliases):
        lc = LineageCollapser(df_basic, "lineage", "count")
        expanded = lc.alias_to_lineage(df_basic["lineage"])
        assert "B.1.617.2.4.1" in expanded.values # AY.4.1 expands
        assert "B.1.1.529.1" in expanded.values # BA.1 expands

def test_reverse_alias_mapping(df_basic, simple_aliases):
    with patch.object(LineageCollapser, "get_pango_aliases", return_value=simple_aliases):
        lc = LineageCollapser(df_basic, "lineage", "count")
        expanded = lc.alias_to_lineage(df_basic["lineage"])
        reversed_back = lc.alias_to_lineage(expanded, reverse=True) # AY.4.1 should map back to AY.4.1 (prefix AY.4 restored)
        assert any(x.startswith("AY.4") for x in reversed_back)


def test_lineage_collapser():
    assert False


def test_validate_inputs(setup_lineage_collapser):
    assert False

def test_unalias_lineage(setup_lineage_collapser):
    assert False

def test_alias_to_lineage(setup_lineage_collapser):
    assert False

def test___get_lineage_level(setup_lineage_collapser):
    assert False

def test_get_pango_aliases(setup_lineage_collapser):
    assert False

def test_get_lineage_level(setup_lineage_collapser):
    assert False

def test_reverse_alias_dict(setup_lineage_collapser):
    assert False

def test___get_lineages_to_collapse(setup_lineage_collapser):
    assert False

def test___remove_highest_lineage_level(setup_lineage_collapser):
    assert False

def test___collapse_lineages(setup_lineage_collapser):
    assert False

def test___convert_back_to_alias_and_drop_level(setup_lineage_collapser):
    assert False

def test_collapse_based_on_threshold(setup_lineage_collapser):
    assert False

def test___get_thresholds_based_on_pct(setup_lineage_collapser):
    assert False

def test_collapse_based_on_pct(setup_lineage_collapser):
    assert False

def test_collapse_recursively_to_at_least_n(setup_lineage_collapser):
    assert False

def test_sort_or_set_thresholds(setup_lineage_collapser):
    assert False
