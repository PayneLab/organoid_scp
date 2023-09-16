# Lung organoid single-cell proteomics

In collaboration with the Van Ry and Kelly labs at BYU.

## `oscutils` software package

A lot analysis functions for this project are contained within a Python package included in the repository called `oscutils` (Organoid Single-Cell UTILitieS). To install and use it, follow the instructions in [`oscutils/README.md`](https://github.com/PayneLab/organoid_scp/blob/main/oscutils/README.md).

## Note about old data

If you delve into the depths of this repo, you will notice that the data were ran twice each with both Proteome Discoverer and MetaMorpheus. However, the older run from each software has issues--specifically, the older PD run is missing one of the pseudo-bulk samples, and the older MM run included the blank, contaminated, no protein, and QC samples, which we found out later can mess up the normalization that MM does, so we ran it again without those samples.

So, even though those older data exist, you shouldn't use them for analysis. The `oscutils` package is set up to load the newer runs by default, and you should only use those.
