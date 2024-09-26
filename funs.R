##### Functions for: Antibiotic resistance associated with the COVID-19 pandemic: a rapid systematic review #####
##### Langford et al. #####
##### Script by Jean-Paul R. Soucy #####
##### https://github.com/jeanpaulrsoucy/covid-19-amr-meta-analysis #####

# define variable names for organisms
orgs <- list(
  "MRSA" = list(
    "num_pre" = "MRSA numerator PRE",
    "denom_pre" = "MRSA denominator PRE",
    "num_post" = "MRSA numerator post",
    "denom_post" = "MRSA denominator post",
    "eff_size" = "MRSA effect size"
  ),
  "VRE" = list(
    "num_pre" = "VRE numerator pre",
    "denom_pre" = "VRE denominator Pre",
    "num_post" = "VRE numerator Post",
    "denom_post" = "VRE denominator Post",
    "eff_size" = "VRE effect size"
  ),
  "Pseudomonas" = list(
    "num_pre" = "Pseudomonas-R numerator pre",
    "denom_pre" = "Pseudomonas-R denominator pre",
    "num_post" = "Pseudomonas-R numerator post",
    "denom_post" = "Pseudomonas-R denominator post",
    "eff_size" = "Pseudomonas effect size"
  ),
  "Acinetobacter" = list(
    "num_pre" = "A. baumanii numerator pre",
    "denom_pre" = "A. baumanii denominator pre",
    "num_post" = "A. baumanii numerator post",
    "denom_post" = "A. baumanii denominator post",
    "eff_size" = "A. baumanii effect size"
  ),
  "ESBL" = list(
    "num_pre" = "ESBL numerator pre",
    "denom_pre" = "ESBL denominator pre",
    "num_post" = "ESBL numerator post",
    "denom_post" = "ESBL denominator post",
    "eff_size" = "ESBL effect size"
  ),
  "CRE" = list(
    "num_pre" = "CRE numerator pre",
    "denom_pre" = "CRE denominator pre",
    "num_post" = "CRE numerator post",
    "denom_post" = "CRE denominator post",
    "eff_size" = "CRE effect size"
  )
)

# function to clean org names
clean_org_name <- function(org) {
  org <- gsub("\\s", "_", org)
  org <- gsub("\\.", "", org)
  tolower(org)
}

# extract raw data for specific organism(s) and summary measure
extract_raw <- function(dat, org, sm) {
  # check arguments
  match.arg(org, choices = names(orgs), several.ok = TRUE)
  match.arg(sm, c("IRR", "RR"))
  # filter based on SM
  if (sm == "IRR") {
    dat2 <- dat %>%
      filter(`AMR numerator` == "cases" &
               `AMR denominator` == "patient days")
  } else {
    dat2 <- dat %>%
      filter(`AMR numerator` == "cases" &
               `AMR denominator` %in% c(
                 "admissions", "discharges", "patients", "total HAIs", "total resistant and non-resistant"))
  }
  # extract data for each organism
  dat3 <- lapply(org, function(x) {
    # filter to studies with complete raw data for organism
    dat3 <- dat2 %>%
      filter(!is.na(!!sym(orgs[[x]][["num_pre"]])) &
               !is.na(!!sym(orgs[[x]][["num_post"]])) &
               !is.na(!!sym(orgs[[x]][["denom_pre"]])) &
               !is.na(!!sym(orgs[[x]][["denom_post"]]))
      )
    # format data
    dat3 <- dat3 %>%
      transmute(
        study = `Author, Year`,
        org = x,
        sm = sm,
        denom = droplevels(factor(case_when(
          sm == "IRR" ~ "Patient days",
          `AMR denominator` %in% c("admissions", "discharges", "patients") ~ "Patients",
          `AMR denominator` %in% c("total HAIs", "total resistant and non-resistant") ~ "Organisms",
          TRUE ~ "ERROR"
        ), levels = c("Patient days", "Patients", "Organisms", "ERROR"))),
        changes_in_ipac_asp,
        num_pre = as.integer(!!sym(orgs[[x]][["num_pre"]])),
        num_post = as.integer(!!sym(orgs[[x]][["num_post"]])),
        denom_pre = as.integer(!!sym(orgs[[x]][["denom_pre"]])),
        denom_post = as.integer(!!sym(orgs[[x]][["denom_post"]]))
      )
  })
  # return data
  bind_rows(dat3)
}

# run meta-analysis for specific organism and summary measure
run_meta <- function(dat, subgroup = "none", subgroup_abort = FALSE) {
  # read summary measure from input data
  sm = dat[["sm"]][1]
  # run meta-analysis
  if (sm == "IRR") {
    if (subgroup %in% c("ipac_asp") & length(unique(dat$changes_in_ipac_asp)) > 1) {
      metainc(
        data = dat,
        event.c = num_pre,
        time.c = denom_pre,
        label.c = "Pre-COVID",
        event.e = num_post,
        time.e = denom_post,
        label.e = "During COVID",
        sm = "IRR",
        method = "GLMM",
        subgroup = case_when(
          subgroup == "ipac_asp" ~ dat$changes_in_ipac_asp
        ),
        subgroup.name = case_when(
          subgroup == "ipac_asp" ~ "IPAC/ASP"#,
          # TRUE ~ "" # label takes up too much space
        )
      )
    } else {
      if (subgroup == "ipac_asp") {
        warning(paste("Ignoring subgroup 'ipac_asp' because all values are:", unique(dat$changes_in_ipac_asp)))
        if (subgroup_abort) {cat("Cannot produce subgroup analysis, aborting...", fill = TRUE); return(NULL)}
      }
      metainc(
        data = dat,
        event.c = num_pre,
        time.c = denom_pre,
        label.c = "Pre-COVID",
        event.e = num_post,
        time.e = denom_post,
        label.e = "During COVID",
        sm = "IRR",
        method = "GLMM")
    }
  } else {
    if ((subgroup == "ipac_asp" & length(unique(dat$changes_in_ipac_asp)) > 1) |
        (subgroup == "denom" & length(unique(dat$denom)) > 1)) {
      metabin(
        data = dat,
        event.c = num_pre,
        n.c = denom_pre,
        label.c = "Pre-COVID",
        event.e = num_post,
        n.e = denom_post,
        label.e = "During COVID",
        sm = "RR",
        method = "MH",
        method.tau = "PM",
        subgroup = case_when(
          subgroup == "ipac_asp" ~ as.character(dat$changes_in_ipac_asp),
          subgroup == "denom" ~ as.character(dat$denom)
        ),
        subgroup.name = case_when(
          subgroup == "ipac_asp" ~ "IPAC/ASP",
          subgroup == "denom" ~ "Denominator"#,
          # TRUE ~ "" # label takes up too much space
        )
      )
    } else {
      if (subgroup == "ipac_asp") {
        warning(paste("Ignoring subgroup 'ipac_asp' because all values are:", unique(dat$changes_in_ipac_asp)))
        if (subgroup_abort) {cat("Cannot produce subgroup analysis, aborting...", fill = TRUE); return(NULL)}
      } else if (subgroup == "denom") {
        warning(paste("Ignoring subgroup 'denom' because all values are:", unique(dat$denom)))
        if (subgroup_abort) {cat("Cannot produce subgroup analysis, aborting...", fill = TRUE); return(NULL)}
      }
      metabin(
        data = dat,
        event.c = num_pre,
        n.c = denom_pre,
        label.c = "Pre-COVID",
        event.e = num_post,
        n.e = denom_post,
        label.e = "During COVID",
        sm = "RR",
        method = "MH",
        method.tau = "PM"
      )
    }
  }
}

# forest plot of meta-analysis for specific organism and summary measure
plot_forest <- function(dat, filename, width, height, include_org_in_study_name = FALSE) {
  # abort if data is NULL
  if (is.null(dat)) {return(invisible(NULL))}
  # add study labels
  if (include_org_in_study_name) {
    dat$studlab <- paste0(dat$data$study, " (", dat$data$org, ")")
  } else {
    dat$studlab <- dat$data$study
  }
  # open graphics device
  png(file = filename, width = width, height = height, units = "in", res = 300)
  if (dat$sm == "IRR") {
    dat %>%
      forest(
        weight.study = "random",
        squaresize = 0.5,
        pooled.totals = TRUE,
        print.tau2 = TRUE,
        digits = 2,
        col.by = "black"
      )
  } else {
    dat %>%
      forest(
        weight.study = "random",
        squaresize = 0.5,
        pooled.totals = TRUE,
        print.tau2 = TRUE,
        digits = 2,
        col.by = "black"
      )
  }
  # close graphics device
  dev.off()
}
