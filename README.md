# Replication Materials for "Explaining Disparities at Quantiles: An Augmented Kitagawa-Oaxaca-Blinder Decomposition"

## Software Requirements

Replication requires R. The required R packages are listed at the beginning of each script.

## Data

The `Data` folder contains two subfolders:

- `CPS Download`: raw CPS data files.
- `Data Cleaning`: the data-cleaning script, `aKOB_datacleaning.R`, and the cleaned analysis file, `CPS_UpdatedData.rds`.

The script `aKOB_datacleaning.R` processes the raw CPS files and produces `CPS_UpdatedData.rds`, which is used in the replication analyses.

## Replication

The main replication script is `aKOB_replication.R`. This script produces regression imputation, weighting, and doubly robust estimates for the union–nonunion wage-gap decomposition using both percentile rank and log-wage outcomes. The resulting estimates are used to construct the tables and figures reported in the main text and appendices.

Before running `aKOB_replication.R`, load the functions contained in the `Functions` folder. This folder includes subfolders with functions for producing results under different estimation strategies and outcome scales.
