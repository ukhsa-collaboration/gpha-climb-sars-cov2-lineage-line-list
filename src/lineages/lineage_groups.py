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
    aliasor = Aliasor()
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


def generate_lineage_groups(start=None, end=None, filename=None, threshold=0.25):
    start_time = datetime.datetime.now()
    group_lineages = set()
    if filename is not None:
        try:
            print(f"{datetime.datetime.now()} Reading previously defined lineages")
            previously_generated_lineages = set(
                json.loads(
                    Path(filename).read_text()
                ).values()
            )
            print(
                f"{datetime.datetime.now()} Adding {len(previously_generated_lineages)} previously defined lineage groups")
            group_lineages.update(previously_generated_lineages)
        except Exception as e:
            print(
                f"{datetime.datetime.now()} Failed to parse lineages from {filename}\n{e}"
            )
            print(f"{datetime.datetime.now()} Run exited in {datetime.datetime.now() - start_time}")
            sys.exit()

        ## actual script
    prevalent_lineages = get_qualifying_lineages(end, start, threshold)  ## generated based on prevalence
    group_lineages.update(prevalent_lineages)
    variant_lineages = (
        get_defined_variant_lineages()
    )  ## defined as VOC/VUIs once upon a time
    group_lineages.update(variant_lineages)

    ## we additionally need to add all root lineages (A, B and recombinants)
    group_lineages.update(get_root_lineages())

    ## that's about it, the rest is just
    ## producing nice outputs
    print(f"{datetime.datetime.now()} Final lineage list generated ({len(group_lineages)} lineage groups)")

    ## begin "nice outputs"
    declared_lineages = get_declared_lineages()  ## every lineage we've heard of

    compare_lineages_key = cmp_to_key(compare_lineages)  ## key for sorting lineages

    aliasor = Aliasor()

    group_lineages_df = pd.DataFrame(
        [(x, aliasor.uncompress(x)) for x in group_lineages]
    ).rename({0: "lineage", 1: "unaliased_lineage"}, axis=1)

    group_lineages_df["lineage_level"] = group_lineages_df["unaliased_lineage"].apply(
        lambda x: len(x.split("."))
    )

    ## okay, now assign every lineage we've ever heard of to a group
    print(f"{datetime.datetime.now()} Assigning every lineage to a group")
    result_df = pd.DataFrame(
        [(x, aliasor.uncompress(x)) for x in declared_lineages]
    ).rename({0: "lineage", 1: "unaliased_lineage"}, axis=1)

    result_df["lineage_group"] = result_df["lineage"].apply(
        lambda x: get_lineage_group(x, group_lineages_df)
    )

    ## sort the dataframe
    ## it would be nice if we could use compare_lineages_key
    ## but this will do
    result_df = result_df.sort_values("unaliased_lineage")

    print(f"{datetime.datetime.now()} Writing out assignments CSV")
    result_df.to_csv(
        f"groups_{datetime.datetime.now().strftime('%Y%m%d-%H%M%S')}.csv", index=False
    )

    print(f"{datetime.datetime.now()} Writing out assignments JSON")
    with open(
            f"lineage_group_lookup_{datetime.datetime.now().strftime('%Y%m%d-%H%M%S')}.json",
            "w",
    ) as outfile:
        outfile.write(
            json.dumps(
                dict(result_df[["lineage", "lineage_group"]].to_records(index=False)),
                indent=4,
            )
        )

    print(f"{datetime.datetime.now()} Run completed in {datetime.datetime.now() - start_time}")