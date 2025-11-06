#!/usr/bin/env python3
import pandas as pd
import datetime
import sys
import json
import requests
import argparse
from glob import glob
from pango_aliasor.aliasor import Aliasor
from functools import cmp_to_key
from pathlib import Path


def compare_lineages(lineageA, lineageB):
    """
    Compares two pangolin lineages and returns
        -1 if a<b
        0 if a==b
        1 if a>b

    Can be used to sort a list as follows:
        compare_lineages_key = cmp_to_key(compare_lineages)
        sorted(["AY.4", "AY.5", "AY.4.1", "AY.4.10", "AY.4.2"], key=compare_lineages_key)
    """
    if not isinstance(lineageA, list):
        splitA = [int(x) if idx else x for idx, x in enumerate(lineageA.split("."))]
    else:
        splitA = lineageA
    if not isinstance(lineageB, list):
        splitB = [int(x) if idx else x for idx, x in enumerate(lineageB.split("."))]
    else:
        splitB = lineageB

    ziplist = zip(splitA, splitB)

    for a, b in ziplist:
        if a == b:
            continue
        if a > b:
            return 1
        if a < b:
            return -1
    if len(splitA) > len(splitB):
        return 1
    if len(splitA) < len(splitB):
        return -1
    return 0


def get_qualifying_lineages(start: None, end: None, threshold: float):
    ## read in latest metadata file
    print(f"{datetime.datetime.now()} Parsing cog_*_all_metadata.csv")
    meta_df = pd.read_csv(
        glob("/cephfs/covid/bham/results/msa/latest/alignments/cog_*_all_metadata.csv")[
            0
        ],
        usecols=["collection_date", "received_date", "usher_lineage"],
        low_memory=False,
    )

    print(f"{datetime.datetime.now()} Processing cog_*_all_metadata.csv")
    ## stage column
    meta_df["consensus_date"] = meta_df["collection_date"]
    ## fill in blanks with received date instead
    meta_df.loc[meta_df["collection_date"].isna(), "consensus_date"] = meta_df[
        "received_date"
    ]
    ## remove any empty dates
    meta_df = meta_df[meta_df["consensus_date"].notna()]
    ## drop unnecessary columns
    meta_df = meta_df.drop(["collection_date", "received_date"], axis=1)

    ## do filter if required
    if start is not None:
        rows_before = len(meta_df)
        formatted_date = datetime.datetime.strptime(start, "%Y%m%d").strftime("%Y-%m-%d")
        meta_df = meta_df[meta_df["consensus_date"] >= formatted_date]
        print(f"{datetime.datetime.now()} Filtered to dates >= {formatted_date} (rows {rows_before}->{len(meta_df)})")
    if end is not None:
        rows_before = len(meta_df)
        formatted_date = datetime.datetime.strptime(end, "%Y%m%d").strftime("%Y-%m-%d")
        meta_df = meta_df[meta_df["consensus_date"] <= formatted_date]
        print(f"{datetime.datetime.now()} Filtered to dates <= {formatted_date} (rows {rows_before}->{len(meta_df)})")

    ## yyyy-mm-dd to yyyy-mm
    meta_df["consensus_date"] = meta_df["consensus_date"].str.slice(0, 7)

    ## group by year-month
    print(f"{datetime.datetime.now()} Grouping lineages by month")
    group_df = (
        meta_df.groupby(["usher_lineage", "consensus_date"])
            .size()
            .reset_index()
            .rename({0: "lineage_count"}, axis=1)
    )

    ## sequences per month dataframe
    ## will be used as a denominator
    print(f"{datetime.datetime.now()} Generating denominators")
    denominator_df = (
        meta_df.groupby("consensus_date")
            .size()
            .reset_index()
            .rename({0: "sequence_count"}, axis=1)
    )

    ## merge the dataframes and calculate percentage prevalences
    print(f"{datetime.datetime.now()} Merging")
    merge_df = pd.merge(group_df, denominator_df, on="consensus_date", how="left")
    merge_df["prevalence_pct"] = (
            100 * merge_df["lineage_count"] / merge_df["sequence_count"]
    )

    threshold_pct = 100 * threshold  ## percent in month to qualify
    print(f"{datetime.datetime.now()} Using a threshold of {threshold_pct}%")

    qualifying_lineages = set(
        merge_df[merge_df["prevalence_pct"] >= threshold_pct]["usher_lineage"]
    )

    ## okay, now get the parent lineage of everything that isn't a qualifying lineage,
    ## sum the lineage counts and recalculate the prevalence_pct
    aliasor = Aliasor()

    def collapse_dataframe(inputDataframe):
        inputDataframe["usher_lineage"] = inputDataframe["usher_lineage"].apply(
            lambda x: aliasor.parent(x) if not x in qualifying_lineages else x
        )

        ## remove anything with a null parent
        inputDataframe = inputDataframe[inputDataframe["usher_lineage"] != ""]

        inputDataframe = inputDataframe.drop(
            ["sequence_count", "prevalence_pct"], axis=1
        )
        inputDataframe = (
            inputDataframe.groupby(["usher_lineage", "consensus_date"])
                .sum()
                .reset_index()
                .rename({0: "lineage_count"}, axis=1)
        )
        inputDataframe = pd.merge(
            inputDataframe, denominator_df, on="consensus_date", how="left"
        )
        inputDataframe["prevalence_pct"] = (
                100 * inputDataframe["lineage_count"] / inputDataframe["sequence_count"]
        )

        temp_qualifying_lineages = set(
            inputDataframe[inputDataframe["prevalence_pct"] >= threshold_pct][
                "usher_lineage"
            ]
        )

        qualifying_lineages.update(temp_qualifying_lineages)

        if len(set(inputDataframe["usher_lineage"]).difference(qualifying_lineages)):
            print(f"{datetime.datetime.now()} Recursing...")
            inputDataframe = collapse_dataframe(inputDataframe)
        else:
            return inputDataframe

    ## not a very interesting dataframe
    fully_collapsed_df = collapse_dataframe(merge_df)

    print(f"{datetime.datetime.now()} Prevalance lineage list generated ({len(qualifying_lineages)} lineage groups)")

    ## this is what we actually want
    return qualifying_lineages


def get_lineage_group(inputLineage, inputGroupDataframe):
    unaliased_input = aliasor.uncompress(inputLineage)

    candidate_list = filter(
        lambda x: unaliased_input.startswith(f"{x}.") or (unaliased_input == x),
        inputGroupDataframe["unaliased_lineage"],
    )

    filter_df = inputGroupDataframe[
        inputGroupDataframe["unaliased_lineage"].isin(candidate_list)
    ].sort_values("lineage_level", ascending=False)

    if len(filter_df):
        return list(filter_df["lineage"])[0]
    else:
        return None


def get_defined_variant_lineages():
    legacy_variants_dict = [
        ("VOC-20DEC-01", "B.1.1.7"),
        ("VOC-20DEC-02", "B.1.351"),
        ("VOC-21FEB-02", "B.1.1.7"),
        ("VOC-21JAN-02", "P.1"),
        ("V-21FEB-01", "A.23.1"),
        ("V-21FEB-03", "B.1.525"),
        ("V-21FEB-04", "B.1.1.318"),
        ("V-21JAN-01", "P.2"),
        ("V-21MAR-01", "B.1.324.1"),
        ("V-21MAR-02", "P.3"),
        ("V-21APR-01", "B.1.617.1"),
        ("VOC-21APR-02", "B.1.617.2"),
        ("V-21APR-03", "B.1.617.3"),
        ("V-21MAY-01", "AV.1"),
        ("V-21MAY-02", "C.36.3"),
        ("V-21JUN-01", "C.37"),
        ("V-21JUL-01", "B.1.621"),
        ("V-21OCT-01", "AY.4.2"),
        ("VOC-21NOV-01", "BA.1"),  ## not B.1.1.529 because I say so
        ("VOC-22JAN-01", "BA.2"),
        ("V-22APR-01", "XD"),
        ("V-22APR-02", "XE"),
        ("V-22APR-03", "BA.4"),
        ("V-22APR-04", "BA.5"),
        ("V-22JUL-01", "BA.2.75"),
        ("V-22SEP-01", "BA.4.6"),
        ("V-22OCT-01", "BQ.1"),
        ("V-22OCT-02", "XBB"),
        ("V-22DEC-01", "CH.1.1"),
        ("V-23JAN-01", "XBB.1.5"),
        ("V-23APR-01", "XBB.1.16"),
        ("V-23JUL-01", "EG.5.1"),
        ("V-23AUG-01", "BA.2.86"),
        ("V-23DEC-01", "JN.1"),
        ("SIM-BA3", "BA.3"),
    ]

    return set([x[1] for x in legacy_variants_dict])


def get_declared_lineages():
    unique_declared_lineages = set(
        pd.read_csv(
            "https://github.com/cov-lineages/pango-designation/raw/master/lineages.csv"
        )["lineage"]
    )

    return unique_declared_lineages


def get_root_lineages():
    alias_key_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json"

    alias_key = requests.get(alias_key_url).json()

    root_lineages = map(
        lambda x: x[0],
        filter(lambda x: (x[1] == "") or (isinstance(x[1], list)), alias_key.items()),
    )

    return set(root_lineages)