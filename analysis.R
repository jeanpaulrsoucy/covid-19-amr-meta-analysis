##### Antibiotic resistance associated with the COVID-19 pandemic: a rapid systematic review #####
##### Langford et al. #####
##### Script by Jean-Paul R. Soucy #####
##### https://github.com/jeanpaulrsoucy/covid-19-amr-meta-analysis #####

# load libraries
library(dplyr)
library(readr)
library(meta)
library(metafor)

# load functions
source("funs.R")

# create directory for outputs
dir.create("out", showWarnings = FALSE)

# load data
dat <- read_csv("data/amr.csv", show_col_types = FALSE)

# create variables
dat <- dat %>%
  mutate(
    changes_in_ipac_asp = factor(
      case_when(
        `Changes in IPAC procedures (Y/N)` == "Y" |
          `Presence of ASP (Y/N)?` == "Y" ~ "Yes",
        TRUE ~ "No or unknown"
      ), levels = c("No or unknown", "Yes")
    )
  )

# IRR forest plots

## MRSA
org <- "MRSA"
extract_raw(dat, org, "IRR") %>%
  run_meta() %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), ".png"),
              width = 9, height = 3)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "IRR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 5)

## VRE
org <- "VRE"
extract_raw(dat, org, "IRR") %>%
  run_meta() %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), ".png"),
              width = 9, height = 3)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "IRR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 5)

## Pseudomonas
org <- "Pseudomonas"
extract_raw(dat, org, "IRR") %>%
  run_meta() %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), ".png"),
              width = 9, height = 3)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "IRR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 5)

## Acinetobacter
org <- "Acinetobacter"
extract_raw(dat, org, "IRR") %>%
  run_meta() %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), ".png"),
              width = 9, height = 3)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "IRR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 5)

## ESBL
org <- "ESBL"
extract_raw(dat, org, "IRR") %>%
  run_meta() %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), ".png"),
              width = 9, height = 3)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "IRR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 5)

## CRE
org <- "CRE"
extract_raw(dat, org, "IRR") %>%
  run_meta() %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), ".png"),
              width = 9, height = 3)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "IRR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/irr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 5)

# RR forest plots

## MRSA
org <- "MRSA"
extract_raw(dat, org, "RR") %>%
  run_meta("denom") %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), ".png"),
              width = 11, height = 6)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "RR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 6)

## VRE
org <- "VRE"
extract_raw(dat, org, "RR") %>%
  run_meta("denom") %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), ".png"),
              width = 11, height = 5)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "RR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 5)

## Pseudomonas
org <- "Pseudomonas"
extract_raw(dat, org, "RR") %>%
  run_meta("denom") %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), ".png"),
              width = 11, height = 5)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "RR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 5)

## Acinetobacter
org <- "Acinetobacter"
extract_raw(dat, org, "RR") %>%
  run_meta("denom") %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), ".png"),
              width = 11, height = 3)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "RR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 3)

## ESBL
org <- "ESBL"
extract_raw(dat, org, "RR") %>%
  run_meta("denom") %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), ".png"),
              width = 11, height = 6)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "RR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 6)

## CRE
org <- "CRE"
extract_raw(dat, org, "RR") %>%
  run_meta("denom") %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), ".png"),
              width = 11, height = 5)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "RR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/rr_", clean_org_name(org), "_ipac_asp.png"),
              width = 11, height = 5)

# forest plots for organism groups

## Gram negative - IRR
org <- c("Pseudomonas", "Acinetobacter", "ESBL", "CRE")
extract_raw(dat, org, "IRR") %>%
  run_meta() %>%
  plot_forest(filename = paste0("out/irr_", "gram_negative", ".png"),
              width = 11, height = 5, include_org_in_study_name = TRUE)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "IRR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/irr_", "gram_negative", "_ipac_asp.png"),
              width = 11, height = 5, include_org_in_study_name = TRUE)

## Gram negative - RR
org <- c("Pseudomonas", "Acinetobacter", "ESBL", "CRE")
extract_raw(dat, org, "RR") %>%
  run_meta("denom") %>%
  plot_forest(filename = paste0("out/rr_", "gram_negative", ".png"),
              width = 12, height = 9, include_org_in_study_name = TRUE)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "RR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/rr_", "gram_negative", "_ipac_asp.png"),
              width = 12, height = 9, include_org_in_study_name = TRUE)

## Gram positive - IRR
org <- c("MRSA", "VRE")
extract_raw(dat, org, "IRR") %>%
  run_meta() %>%
  plot_forest(filename = paste0("out/irr_", "gram_positive", ".png"),
              width = 11, height = 5, include_org_in_study_name = TRUE)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "IRR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/irr_", "gram_positive", "_ipac_asp.png"),
              width = 11, height = 7, include_org_in_study_name = TRUE)

## Gram positive - RR
org <- c("MRSA", "VRE")
extract_raw(dat, org, "RR") %>%
  run_meta("denom") %>%
  plot_forest(filename = paste0("out/rr_", "gram_positive", ".png"),
              width = 11, height = 7, include_org_in_study_name = TRUE)
## plot by IPAC/ASP changes (will abort if all studies have the same value)
extract_raw(dat, org, "RR") %>%
  run_meta("ipac_asp", subgroup_abort = TRUE) %>%
  plot_forest(filename = paste0("out/rr_", "gram_positive", "_ipac_asp.png"),
              width = 11, height = 7, include_org_in_study_name = TRUE)
