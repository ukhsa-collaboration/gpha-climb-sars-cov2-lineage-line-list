GPHA CLIMB SARs-CoV2 Lineage Linelist

04.11.26: Debug

Codebase:

- DEBUG: update ONT folder identification method due to old method becoming obsolete due to Illumina changes to way they record their runs.
  - Old method: 8 digit date presence in run folder name
  - New method: 8 digit date and 4 digit time presence in run folder name
- updating of os.path instances -> pathlib
- removal of print statements
- removal of unused code

Features:

- VERSION: update to 3.1.1