# Antibiotic resistance associated with the COVID-19 pandemic: a systematic review and meta-analysis

## Purpose of this repository

This repository contains all data and code necessary to reproduce the analyses presented in the manuscript "[Antibiotic resistance associated with the COVID-19 pandemic: a systematic review and meta-analysis](https://doi.org/10.1016/j.cmi.2022.12.006)" by Langford et al (doi: [10.1016/j.cmi.2022.12.006](https://doi.org/10.1016/j.cmi.2022.12.006)). **This includes every figure and table except for figure 1 and table 1. Table 1 is a subset of the dataset found in the file `data/amr.csv`.** Note that the analysis script produces additional figures beyond those directly presented in the manuscript.

## Requirements

All code is written in the [R programming language](https://www.r-project.org/). The easiest way to run it is to use the [RStudio](https://rstudio.com/) IDE. An `.Rproj` file is included with this repository for ease of use with RStudio. The analysis script should run with any modern version of R.

The R packages required to reproduce the figures are listed at the top of the analysis script. They must be installed using `install.packages` or similar functionality within RStudio prior to running the script.

## Reproducing figures

Run `analysis.R` to create the output figures.

## Figures

- [Figure 2](out/irr_gram_positive.png)
- [Figure 3](out/rr_gram_positive_ipac_asp.png)
- [Figure 4](out/irr_gram_negative.png)
- [Figure 5](out/rr_gram_negative_ipac_asp.png)
