library(ipumsr)
library(haven)
library(tidyverse)
library(kdensity)
library(DescTools)
library(spatstat.geom)
library(Hmisc)

setwd("CPS Download")

UpdatedUnion_ddi <- read_ipums_ddi("cps_00004.xml")
UpdatedUnion_data <- read_ipums_micro(UpdatedUnion_ddi)

UpdatedUnion_dataremove <- UpdatedUnion_data %>%
  filter(MISH %in% c(4, 8), between(AGE, 16, 64), EMPSTAT %in% c(10, 12)) %>%
  mutate(
    valid_hourly =
      !is.na(HOURWAGE) &
      HOURWAGE > 0     &
      HOURWAGE < 999.99,
    valid_weekly =
      !is.na(EARNWEEK) &
      EARNWEEK > 0 &
      EARNWEEK < 9999.99 &
      !is.na(UHRSWORKORG) &
      UHRSWORKORG > 0 &
      UHRSWORKORG < 999
  ) %>%
  filter(valid_hourly | valid_weekly) %>%
  filter(SEX == 1)

UpdatedUnion_dataclean <- UpdatedUnion_dataremove %>%
  mutate(
    CPI99_recode = case_when(
      YEAR == 2015 ~ 0.704,
      YEAR == 2016 ~ 0.703,
      YEAR == 2017 ~ 0.694,
      YEAR == 2018 ~ 0.679,
      YEAR == 2019 ~ 0.663
    ),
    covered_union = case_when(UNION == 1 ~ 0, UNION == 2 |
                                UNION == 3 ~ 1, TRUE ~ NA_real_),
    union_recoded = fct_relevel(ifelse(covered_union == 1, "Union", "Nonunion"), "Nonunion"),
    age = AGE,
    educ_yrs = case_when(
      EDUC %in% c(0, 1, 2)       ~ 0,
      EDUC == 10                      ~ 2.5,
      EDUC == 11                      ~ 1,
      EDUC == 12                      ~ 2,
      EDUC == 13                      ~ 3,
      EDUC == 14                      ~ 4,
      EDUC == 20                      ~ 5.5,
      EDUC == 21                      ~ 5,
      EDUC == 22                      ~ 6,
      EDUC == 30                      ~ 7.5,
      EDUC == 31                      ~ 7,
      EDUC == 32                      ~ 8,
      EDUC == 40                      ~ 9,
      EDUC == 50                      ~ 10,
      EDUC == 60                      ~ 11,
      EDUC %in% c(70, 71, 72, 73)  ~ 12,
      EDUC %in% c(80, 81)            ~ 13,
      EDUC %in% c(90, 91, 92)       ~ 14,
      EDUC == 100                      ~ 15,
      EDUC %in% c(110, 111)            ~ 16,
      EDUC %in% c(120, 121)            ~ 17,
      EDUC == 122                      ~ 18,
      EDUC == 123                      ~ 18,
      EDUC == 124                      ~ 19,
      EDUC == 125                      ~ 20,
      TRUE                             ~ NA_real_
    ),
    educ_cat = factor(
      case_when(
        educ_yrs < 9 ~ "Elementary",
        educ_yrs >= 9 &
          educ_yrs < 12 ~ "HS_dropout",
        educ_yrs == 12 ~ "High_school",
        educ_yrs > 12 &
          educ_yrs <= 15 ~ "Some_college",
        educ_yrs == 16 ~ "College",
        educ_yrs > 16 ~ "Post_graduate"
      ),
      levels = c(
        "High_school",
        "Elementary",
        "HS_dropout",
        "Some_college",
        "College",
        "Post_graduate"
      )
    ),
    potexp = pmax(0, age - educ_yrs - 6),
    potexp_cat = fct_relevel(
      case_when(
        potexp < 5 ~ "Less_than_5",
        potexp >= 5 &
          potexp < 10 ~ "Between_5_9",
        potexp >= 10 &
          potexp < 15 ~ "Between_10_14",
        potexp >= 15 &
          potexp < 20 ~ "Between_15_19",
        potexp >= 20 &
          potexp < 25 ~ "Between_20_24",
        potexp >= 25 &
          potexp < 30 ~ "Between_25_29",
        potexp >= 30 &
          potexp < 35 ~ "Between_30_34",
        potexp >= 35 &
          potexp < 40 ~ "Between_35_39",
        potexp >= 40 ~ "More_than_40"
      ),
      "Between_20_24"
    ),
    hourlywage = case_when(
      valid_hourly ~ HOURWAGE,
      !valid_hourly &
        valid_weekly ~ EARNWEEK / UHRSWORKORG
    ),
    hourlywage_adj99 = hourlywage * CPI99_recode,
    loghourlywage_adj99 = log(hourlywage_adj99),
    marital_recoded = case_when(MARST == 1 |
                                  MARST == 2 ~ 1, between(MARST, 3, 7) ~ 0, TRUE ~ NA_real_),
    nonwhite_recoded = case_when(RACE == 100 ~ 0, between(RACE, 200, 830) ~ 1, TRUE ~ NA_real_),
    id = row_number()
  ) %>%
  relocate(id)

CPS_UpdatedData <- UpdatedUnion_dataclean %>%
  select(
    id,
    eweight = EARNWT,
    lwage = loghourlywage_adj99,
    union = union_recoded,
    covered = covered_union,
    education = educ_cat,
    experience = potexp_cat,
    marr = marital_recoded,
    nonwhite = nonwhite_recoded
  )

lwage_ecdf <- ewcdf(CPS_UpdatedData$lwage, weights = CPS_UpdatedData$eweight)

CPS_UpdatedData <- CPS_UpdatedData %>%
  mutate(
    lwage_rank = lwage_ecdf(lwage),
    lwage_rankalt = wtd.rank(signif(lwage, 10), weights = eweight) / sum(eweight)
  )

saveRDS(CPS_UpdatedData, "Data Cleaning/CPS_UpdatedData.rds")
