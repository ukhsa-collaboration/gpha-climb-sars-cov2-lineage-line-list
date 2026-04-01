
# <img src="https://upload.wikimedia.org/wikipedia/fi/6/69/UKHSA_Logo.svg.png" width="75" height="75"> GPHA CLIMB SARs-CoV-2 Lineage Line List 
---

## Table of Contents
- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#install)
- [Usage](#usage)
- [Config settings](<#config settings>)
- [Troubleshooting](#troubleshooting)

---

## Overview

This repo contains code to generate the covid linelist and lineage
groupings.

## Requirements
TODO: Add in section on accessing data from sftp

## Installation
Clone repo:
```bash
git clone https://github.com/ukhsa-collaboration/gpha-climb-sars-cov2-lineage-line-list.git
cd gpha-climb-sars-cov2-lineage-line-list
conda create -n linelist_env python
conda activate linelist_env
```
User installation:
```bash
pip install .
```
Developer installation:
```bash
pip install -e ".[dev]"
```

## Usage

```bash
python run_linelist.py --input /path/to/combined_pangolin_results.csv --output /path/to/results/folder/
```

## Config settings

### Global params
Global params are parameters that influence how the overall code runs or are used in both the recent and full reporting window.

| Parameter       | Description                              |
|-----------------|------------------------------------------|
| period_size     | Size of period to break the reporting window into **in days**. Currently set to 14 (2 weeks). |
| defined_variant | *Not yet populated*. Contains list of previously defined lineages to protect from collapsing/collapse no further than.  |

### Recent reporting window

Recent reporting window refers to the period of time in which percentage
based lineage prevalence is calculated. The recent reporting window is
split into  reporting periods (defined in global params as period_size)
and a lineage must reach the percentage prevalence threshold in at least
one of those periods to be protected in lineage grouping.

Parameters relating to the initial reporting window:

| Parameter          | Description                              |
|--------------------|------------------------------------------|
| max_lineages       | Maximum number of lineages to be protected. If more than max_lineages are identified to protect based on prevalence, these will be collapsed until max_lineages is reached. |
| min_lineages       | Minimum number of lineages to be protected. If less than min_lineages are identified to protect based on prevalence, lineages will be collapsed until at least min_lineages are identified to protect. |
| weeks_to_include   | Number of weeks to include in the recent reporting period. |
| percent_prevalence | Percent prevalence a lineage must reach to be protected. |

### Full reporting window

Full reporting window refers to the total period of time to include
samples from, excluding samples falling in the initial reporting window.
Returns the specified number of most prevalent lineages across the whole
reporting window. This is calculated as the full reporting window
max_lineages value minus the number of lineages identified to protect from
the recent reporting window. As a minimum, an additional two lineages
should always be protected in the full reporting window to account for
lineages that had high prevalence previously that have since decreased.

Parameters relating to the full reporting window:

| Parameter        | Description                              |
|------------------|------------------------------------------|
| weeks_to_include | Number of weeks to include in the full reporting period. |
| max_lineages     | Maximum number of lineages to be protected. |

## Troubleshooting
