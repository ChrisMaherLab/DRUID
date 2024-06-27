## DRUID: Deriving CircRNA Universal IDentifiers
##

Developed by the [Christopher Maher Lab](http://www.maherlab.com) at [Washington University in St. Louis](http://wustl.edu).

For additional details, see publication: PLACEHOLDER_LINK.

##

## Overview

DRUID (Deriving CircRNA Universal IDentifiers) is an open-source R package (with MIT license) that conveniently convert BED12 or BED16 coordinates of circular RNAs into Universal Identifiers (UID) according to nomenclature proposed by [Chen et al. (PMID: 36658223)](https://www.nature.com/articles/s41556-022-01066-9).

NOTE: For BED6 input coordinates, DRUID can only label the terminal exon blocks as intermediate full-length anatomy of circRNAs are unknown.

## Quick Start

Install CIRCUS in R with `devtools::install_github('ChrisMaherLab/DRUID')`

Note: The following pre-requisites need to be installed:

- `GenomicRanges (>= 1.50.2)`
- `genomation (>= 1.30.0)`
- `tidyr (>= 1.3.0)`
- `dplyr (>= 1.1.4)`
- `methods`
- `utils`
- `S4Vectors`

To load CIRCUS: `library(DRUID)`

To open help manual: `?DRUID`

## Input Format

4 input parameters are required:
| Variable | Description |
| --- | --- |
| ref_gpf_path | A txt file of two columns with no headers: first column with gene names, second column with transcript IDs. It should match all transcripts in ref_path. See `transcript2gene.txt` in test data. |
| ref_path | BED 12 file with no headers of reference transcripts for annotation. Must contain all 12 columns indicating exon composition of each transcript. 4th column should contain transcript IDs. See `annotation.bed` in test data. |
| bed6 | Boolean (T/F) indicator of whether circRNA BED file is BED12 or BED6. Note: When BED6 input is used, output circRNA UIDs will only contain information about the 5' and 3' terminal exon blocks (e.g., `circZFAND6(3,5)`) without intermediate exons (e.g., a BED12 input would result in `circZFAND6(3,4L,5)`). |
| bed_path | BED file with circRNAs to convert into circRNA UID. 4th column should contain an original identifier of choice (e.g., circRNA_1). Indicate whether this file is BED12 or BED6 with the bed6 parameter. See `bed6_coords.bed` and `bed12_coords.bed` in test data. |

## Test Cases

1. `test_bed12 <- DRUID(ref_gpf_path = system.file("extdata",'transcript2gene.txt',package="DRUID"),ref_path = system.file("extdata",'annotation.bed',package="DRUID"),bed_path = system.file("extdata",'bed12_coords.bed',package="DRUID"),bed6=FALSE)`

Expected output should look like `test_bed12.txt` in test data.

2. `test_bed6 <- DRUID(ref_gpf_path = system.file("extdata",'transcript2gene.txt',package="DRUID"),ref_path = system.file("extdata",'annotation.bed',package="DRUID"),bed_path = system.file("extdata",'bed6_coords.bed',package="DRUID"),bed6=TRUE)`

Expected output should look like `test_bed6.txt` in test data.
