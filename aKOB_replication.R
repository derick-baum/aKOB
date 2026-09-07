# Empirical illustration for the aKOB paper
# 
# This script produces regression imputation, weighting, and doubly robust estimates
# for the union–nonunion wage-gap decomposition using both
# percentile-rank and log-wage outcomes.
#
# Results are used to construct the tables and figures reported in the
# main text and appendices.

library(tidyverse)
library(spatstat)
library(DescTools)
library(Hmisc)
library(rlang)

# Load data and functions -------

CPS_UpdatedData <- readRDS("Data/Data Cleaning/CPS_UpdatedData.rds")
codes_path_RI <- list.files("Functions/Regression Imputation", full.names = TRUE)
codes_path_W <- list.files("Functions/Weighting", full.names = TRUE)
codes_path_DR <- list.files("Functions/Doubly Robust", full.names = TRUE)
codes_path_DR_nobootstrap <- list.files("Functions/Doubly Robust/No Bootstrapping", full.names = TRUE)

codes_source_RI <- walk(codes_path_RI, source)
codes_source_W <- walk(codes_path_W, source)
codes_source_DR <- walk(codes_path_DR[1:4], source)
codes_source_DR_nobootstrap <- walk(codes_path_DR_nobootstrap, source)

# Regression Imputation results -------

## Rank outcome -------

aKOB_riOLS_results_meanrank <- aKOB_riOLS_meanrank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  bootstrap_iters = 250
)

aKOB_riLogit_results_10percentrank <- aKOB_riLogit_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.10,
  bootstrap_iters = 250
)

aKOB_riLogit_results_25percentrank <- aKOB_riLogit_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.25,
  bootstrap_iters = 250
)

aKOB_riLogit_results_50percentrank <- aKOB_riLogit_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.50,
  bootstrap_iters = 250
)

aKOB_riLogit_results_75percentrank <- aKOB_riLogit_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.75,
  bootstrap_iters = 250
)

aKOB_riLogit_results_90percentrank <- aKOB_riLogit_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.90,
  bootstrap_iters = 250
)

aKOB_ri_results_rank <- list(
  aKOB_riOLS_results_meanrank,
  aKOB_riLogit_results_10percentrank,
  aKOB_riLogit_results_25percentrank,
  aKOB_riLogit_results_50percentrank,
  aKOB_riLogit_results_75percentrank,
  aKOB_riLogit_results_90percentrank
)
names(aKOB_ri_results_rank) <- c("meanrank", paste0(c(10, 25, 50, 75, 90), "percentrank"))

## Original outcome -------

aKOB_riOLS_results_meanoriginal <- aKOB_riOLS_meanoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  bootstrap_iters = 250
)

aKOB_riLogit_results_10percentoriginal <- aKOB_riLogit_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.10,
  bootstrap_iters = 250
)

aKOB_riLogit_results_25percentoriginal <- aKOB_riLogit_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.25,
  bootstrap_iters = 250
)

aKOB_riLogit_results_50percentoriginal <- aKOB_riLogit_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.50,
  bootstrap_iters = 250
)

aKOB_riLogit_results_75percentoriginal <- aKOB_riLogit_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.75,
  bootstrap_iters = 250
)

aKOB_riLogit_results_90percentoriginal <- aKOB_riLogit_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.90,
  bootstrap_iters = 250
)

aKOB_ri_results_original <- list(
  aKOB_riOLS_results_meanoriginal,
  aKOB_riLogit_results_10percentoriginal,
  aKOB_riLogit_results_25percentoriginal,
  aKOB_riLogit_results_50percentoriginal,
  aKOB_riLogit_results_75percentoriginal,
  aKOB_riLogit_results_90percentoriginal
)
names(aKOB_ri_results_original) <- c("meanoriginal", paste0(c(10, 25, 50, 75, 90), "percentoriginal"))

# Weighting results -------

## Rank outcome -------

aKOB_wLogit_results_meanrank <- aKOB_wLogit_meanrank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  bootstrap_iters = 250
)

aKOB_wLogit_results_10percentrank <- aKOB_wLogit_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.10,
  bootstrap_iters = 250
)

aKOB_wLogit_results_25percentrank <- aKOB_wLogit_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.25,
  bootstrap_iters = 250
)

aKOB_wLogit_results_50percentrank <- aKOB_wLogit_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.50,
  bootstrap_iters = 250
)

aKOB_wLogit_results_75percentrank <- aKOB_wLogit_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.75,
  bootstrap_iters = 250
)

aKOB_wLogit_results_90percentrank <- aKOB_wLogit_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.90,
  bootstrap_iters = 250
)

aKOB_wLogit_results_rank <- list(
  aKOB_wLogit_results_meanrank,
  aKOB_wLogit_results_10percentrank,
  aKOB_wLogit_results_25percentrank,
  aKOB_wLogit_results_50percentrank,
  aKOB_wLogit_results_75percentrank,
  aKOB_wLogit_results_90percentrank
)
names(aKOB_wLogit_results_rank) <-  c("meanrank", paste0(c(10, 25, 50, 75, 90), "percentrank"))

## Original outcome -------

aKOB_wLogit_results_meanoriginal <- aKOB_wLogit_meanoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  bootstrap_iters = 250
)

aKOB_wLogit_results_10percentoriginal <- aKOB_wLogit_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.10,
  bootstrap_iters = 250
)

aKOB_wLogit_results_25percentoriginal <- aKOB_wLogit_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.25,
  bootstrap_iters = 250
)

aKOB_wLogit_results_50percentoriginal <- aKOB_wLogit_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.50,
  bootstrap_iters = 250
)

aKOB_wLogit_results_75percentoriginal <- aKOB_wLogit_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.75,
  bootstrap_iters = 250
)

aKOB_wLogit_results_90percentoriginal <- aKOB_wLogit_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.90,
  bootstrap_iters = 250
)

aKOB_wLogit_results_original <- list(
  aKOB_wLogit_results_meanoriginal,
  aKOB_wLogit_results_10percentoriginal,
  aKOB_wLogit_results_25percentoriginal,
  aKOB_wLogit_results_50percentoriginal,
  aKOB_wLogit_results_75percentoriginal,
  aKOB_wLogit_results_90percentoriginal
)
names(aKOB_wLogit_results_original) <- c("meanoriginal", paste0(c(10, 25, 50, 75, 90), "percentoriginal"))

# DR results -------

## Rank outcome -------

aKOB_DR_results_meanrank <- aKOB_DR_meanrank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  bootstrap_iters = 250
)
aKOB_DR_results_10percentrank <- aKOB_DR_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.10,
  bootstrap_iters = 250
)
aKOB_DR_results_25percentrank <- aKOB_DR_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.25,
  bootstrap_iters = 250
)
aKOB_DR_results_50percentrank <- aKOB_DR_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.50,
  bootstrap_iters = 250
)
aKOB_DR_results_75percentrank <- aKOB_DR_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.75,
  bootstrap_iters = 250
)
aKOB_DR_results_90percentrank <- aKOB_DR_quantilerank_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage_rankalt,
  prob = 0.90,
  bootstrap_iters = 250
)

aKOB_DR_results_rank <- list(
  aKOB_DR_results_meanrank,
  aKOB_DR_results_10percentrank,
  aKOB_DR_results_25percentrank,
  aKOB_DR_results_50percentrank,
  aKOB_DR_results_75percentrank,
  aKOB_DR_results_90percentrank
)
names(aKOB_DR_results_rank) <- c("meanrank", paste0(c(10, 25, 50, 75, 90), "percentrank"))

## Original outcome -------

aKOB_DR_results_meanoriginal <- aKOB_DR_meanoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  bootstrap_iters = 250
)
aKOB_DR_results_10percentoriginal <- aKOB_DR_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.10,
  bootstrap_iters = 250
)
aKOB_DR_results_25percentoriginal <- aKOB_DR_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.25,
  bootstrap_iters = 250
)
aKOB_DR_results_50percentoriginal <- aKOB_DR_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.50,
  bootstrap_iters = 250
)
aKOB_DR_results_75percentoriginal <- aKOB_DR_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.75,
  bootstrap_iters = 250
)
aKOB_DR_results_90percentoriginal <- aKOB_DR_quantileoriginal_fun(
  data = CPS_UpdatedData,
  survey_weights = eweight,
  outcome = lwage,
  prob = 0.90,
  bootstrap_iters = 250
)

aKOB_DR_results_original <- list(
  aKOB_DR_results_meanoriginal,
  aKOB_DR_results_10percentoriginal,
  aKOB_DR_results_25percentoriginal,
  aKOB_DR_results_50percentoriginal,
  aKOB_DR_results_75percentoriginal,
  aKOB_DR_results_90percentoriginal
)
names(aKOB_DR_results_original) <- c("meanoriginal", paste0(c(10, 25, 50, 75, 90), "percentoriginal"))

## Rank outcome (No Bootstrapping) -------

probs <- seq(0.05, 0.95, 0.01)

aKOB_DR_results_rank_noboot <- map(
  probs,
  ~ aKOB_DR_quantilerank_nobootfun(
    data = CPS_UpdatedData,
    survey_weights = eweight,
    outcome = lwage_rankalt,
    prob = .x
  ),
  .progress = TRUE
)

names(aKOB_DR_results_rank_noboot) <- paste0(seq(5, 95, 1), "percentrank")

## Original outcome (No Bootstrapping) -------

aKOB_DR_results_original_noboot <- map(
  probs,
  ~ aKOB_DR_quantileoriginal_nobootfun(
    data = CPS_UpdatedData,
    survey_weights = eweight,
    outcome = lwage,
    prob = .x
  ),
  .progress = TRUE
)

names(aKOB_DR_results_original_noboot) <- paste0(seq(5, 95, 1), "percentoriginal")
