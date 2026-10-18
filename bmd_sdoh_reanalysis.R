# =============================================================================
#  Social determinants of health (SDoH) and femoral neck BMD T-score
#  Re-analysis for submission to BMC Medical Informatics and Decision Making
#
#  Oct 2026 update: the last cycle is set in nhanes_cycles.R (environment
#  variable NHANES_LATEST). Default "P" = 2017-March 2020 pre-pandemic file in
#  place of 2017-2018 (results in results_2017_2020/); "J" reproduces the
#  submitted 2017-2018 analysis exactly (results in results/).
#
#  What this script does (in order)
#    1. Downloads NHANES 2005-2006, 2007-2008, 2009-2010, 2013-2014 and
#       2017-March 2020 (or 2017-2018) - all cycles with femur DXA - and caches
#       the .xpt files locally.
#    2. Harmonizes variables across cycles and builds the analytic sample
#       (adults >= 50 y with valid femoral neck BMD); records flow counts.
#    3. Survey-weighted descriptive table (Table 1) with missing counts.
#    4. Multiple imputation (mice), run separately by sex, EXCLUDING the
#       outcome so imputation cannot leak outcome information into CV.
#    5. Repeated 10-fold cross-validation (default 10 repeats), sex-stratified:
#         - Random forest (ranger), survey weights as case weights,
#           hyperparameters tuned INSIDE each training fold by OOB error
#         - Weighted linear regression with natural splines (benchmark)
#       Models: base (no SDoH), full (+ education, income),
#               extended (+ food security, insurance, usual care, partnered),
#               and race-free versions of base/full.
#       Paired differences (full - base, etc.) with Nadeau-Bengio corrected
#       resampled-t 95% CIs.
#    6. Held-out permutation importance for each predictor and for the SDoH
#       GROUP (education + income permuted jointly), with CIs and rank spread.
#    7. Secondary, using the WHO categories (normal T >= -1.0; low bone mass
#       -2.5 < T < -1.0; osteoporosis T <= -2.5):
#         - weighted AUC for osteoporosis vs not, and for low bone mass OR
#           osteoporosis vs normal, using predicted T-score as the risk score
#         - 3-category agreement (weighted accuracy, linear-weighted kappa)
#           between observed and predicted WHO category.
#    8. Accumulated local effects (ALE) plots for income and education.
#    9. Sensitivity analyses: complete-case; delta-adjusted MNAR for income;
#       2013-2018 only; excluding osteoporosis-medication users.
#   10. Survey-weighted linear model (svyglm, pooled over imputations) with
#       sex x education and sex x income interaction tests.
#
#  Outputs go to CONFIG$out_dir. Every table the manuscript needs is a CSV.
#
#  IMPORTANT - verify before running in production:
#    * Variable names/codes marked  # VERIFY  against the NHANES codebooks for
#      each cycle (harmonization across 2005-2018 is the main risk area).
#    * Run once with CONFIG$quick <- TRUE (minutes) before the full run.
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
pkgs <- c("haven", "dplyr", "tidyr", "purrr", "readr", "survey", "mitools",
          "mice", "ranger", "iml", "ggplot2")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) install.packages(missing_pkgs)
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(purrr); library(readr)
  library(survey); library(mitools); library(mice); library(ranger)
  library(splines); library(iml); library(ggplot2)
})
options(survey.lonely.psu = "adjust", dplyr.summarise.inform = FALSE)

# ---- 1. Configuration --------------------------------------------------------
# Project folder on OSC (note the spaces in "Random Forest ML" - keep the quotes)
PROJECT_DIR <- Sys.getenv("NHANES_BMD_DIR", "/users/PAS2982/summertea1109/Random Forest ML/NHANES_BMD")
if (dir.exists(PROJECT_DIR)) setwd(PROJECT_DIR)
message("Working directory: ", getwd())
source("nhanes_cycles.R")   # cycles, file names, pooled weights, results folder
EXPECTED_FILES <- if (USE_P) 70 else 71   # P: 14 files x 4 two-year cycles + 14 P_ files
if (!dir.exists("nhanes_cache") || length(list.files("nhanes_cache", "\\.xpt$")) < EXPECTED_FILES)
  warning("Expected at least ", EXPECTED_FILES, " .xpt files in ", file.path(getwd(), "nhanes_cache"),
          " - check the upload. Missing files will be downloaded if internet is available.")

# Use only the cores this session was given (on OSC, detectCores() reports the
# whole node, and over-subscribing it makes everything crawl)
detect_threads <- function() {
  for (v in c("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_JOB_CPUS_PER_NODE")) {
    n <- suppressWarnings(as.integer(sub("\\(.*", "", Sys.getenv(v))))
    if (!is.na(n) && n > 0) return(n)
  }
  n <- tryCatch(length(parallel::mcaffinity()), error = function(e) NA_integer_)
  if (!is.na(n) && n > 0) return(n)
  max(1L, min(4L, parallel::detectCores() - 1L))
}

CONFIG <- list(
  # Fast test run: sbatch --export=ALL,QUICK=TRUE run_osc.slurm  (or set TRUE here)
  quick        = toupper(Sys.getenv("QUICK", "FALSE")) == "TRUE",
  cache_dir    = "nhanes_cache",
  out_dir      = RESULTS_ROOT,     # "results" (2017-2018) or "results_2017_2020"  
  seed         = 20261001,
  m_imp        = 5,              # imputations
  k_folds      = 10,
  n_repeats    = 10,
  num_trees    = 500,
  tune_trees   = 300,
  grid         = expand.grid(mtry = c(2, 4, 6), min.node.size = c(10, 50)),
  n_perm       = 3,              # permutations per variable per fold
  # On OSC/SLURM use the cores you requested; elsewhere all cores but one
  threads      = detect_threads(),   # cores actually allocated (RStudio/OnDemand or sbatch)
  # NHANES III femoral neck reference, non-Hispanic white women aged 20-29
  # (Looker et al., Osteoporos Int 1998). ISCD recommends this reference for
  # T-scores in both sexes.                                            # VERIFY
  ref_mean     = 0.858,
  ref_sd       = 0.120,
  mnar_deltas  = c(-0.5, -1.0)   # shifts applied to imputed PIR (PIR units)
)
if (CONFIG$quick) {
  CONFIG$m_imp <- 2; CONFIG$k_folds <- 5; CONFIG$n_repeats <- 1
  CONFIG$num_trees <- 100; CONFIG$tune_trees <- 100; CONFIG$n_perm <- 1
  CONFIG$grid <- expand.grid(mtry = c(2, 4), min.node.size = 20)
}
dir.create(CONFIG$cache_dir, showWarnings = FALSE)
dir.create(CONFIG$out_dir, showWarnings = FALSE)
set.seed(CONFIG$seed)
message(sprintf("Mode: %s | cores used: %d | %d folds x %d repeats x %d imputations",
                if (CONFIG$quick) "QUICK TEST" else "FULL RUN", CONFIG$threads,
                CONFIG$k_folds, CONFIG$n_repeats, CONFIG$m_imp))
out <- function(f) file.path(CONFIG$out_dir, f)

# ---- 2. Download / read NHANES files ----------------------------------------
CYCLES <- CYC[, c("cycle", "suffix", "year")]   # from nhanes_cycles.R
COMPONENTS <- c("DEMO", "DXXFEM", "BMX", "SMQ", "ALQ", "PAQ", "MCQ", "OSQ",
                "RHQ", "DSQTOT", "RXQ_RX", "FSQ", "HIQ", "HUQ")

read_nhanes <- function(component, suffix, year, cycle) {
  fname <- nhanes_file(component, suffix)          # DEMO_J or P_DEMO
  path  <- file.path(CONFIG$cache_dir, paste0(fname, ".xpt"))
  if (!file.exists(path)) {
    urls <- c(
      sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s.xpt", year, fname),
      sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s.XPT", cycle, fname)
    )
    ok <- FALSE
    for (u in urls) {
      ok <- tryCatch({
        utils::download.file(u, path, mode = "wb", quiet = TRUE); TRUE
      }, error = function(e) FALSE, warning = function(w) FALSE)
      # CDC sometimes returns an HTML error page instead of a 404
      if (ok && file.size(path) < 2000) { unlink(path); ok <- FALSE }
      if (ok) break
    }
    if (!ok) { unlink(path); warning("Could not download ", fname); return(NULL) }
  }
  haven::zap_labels(haven::read_xpt(path))
}

# Variables kept from each component (any_of: absent ones are skipped)
KEEP <- list(
  DEMO   = c("SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "RIDEXPRG",
             "DMDEDUC2", "INDFMPIR", "DMDMARTL", "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
             "DMDMARTZ", "WTMECPRP"),                               # 2017-March 2020 names
  DXXFEM = c("SEQN", "DXXNKBMD"),
  BMX    = c("SEQN", "BMXWT", "BMXHT"),
  SMQ    = c("SEQN", "SMQ020", "SMQ040"),
  ALQ    = c("SEQN", "ALQ101", "ALQ110", "ALQ120Q", "ALQ120U",   # 2005-2014
             "ALQ111", "ALQ121",                                 # 2017-2018 and 2017-2020
             "ALQ130"),
  PAQ    = c("SEQN", "PAD200", "PAD320",                         # 2005-2006
             "PAQ605", "PAQ620", "PAQ650", "PAQ665"),            # 2007+
  MCQ    = c("SEQN", "MCQ160A", "MCQ190", "MCQ191", "MCQ195", "MCQ160M"),
  DSQ1   = c("SEQN", "DSDCOUNT"),                                 # 2005-2006 only
  OSQ    = c("SEQN", "OSQ010A", "OSQ010B", "OSQ010C", "OSQ130",
             "OSQ170", "OSQ200"),
  RHQ    = c("SEQN", "RHQ031", "RHQ131", "RHQ540"),
  DSQTOT = c("SEQN", "DSDCOUNT", "DSQTCALC", "DSQTVD"),
  FSQ    = c("SEQN", "FSDAD"),
  HIQ    = c("SEQN", "HIQ011"),
  HUQ    = c("SEQN", "HUQ030")
)

OSTEO_DRUGS <- paste(c("ALENDRONATE", "RISEDRONATE", "IBANDRONATE", "ZOLEDRONIC",
                       "DENOSUMAB", "RALOXIFENE", "TERIPARATIDE", "ABALOPARATIDE",
                       "ROMOSOZUMAB", "CALCITONIN"), collapse = "|")

read_cycle <- function(cycle, suffix, year) {
  message("Reading NHANES ", cycle)
  comps <- setdiff(COMPONENTS, "RXQ_RX")
  # 2005-2006 has no DSQTOT file; supplement count comes from DSQ1_D instead
  if (suffix == "D") comps <- c(setdiff(comps, "DSQTOT"), "DSQ1")
  parts <- lapply(comps, function(comp) {
    d <- read_nhanes(comp, suffix, year, cycle)
    if (is.null(d)) return(NULL)
    d %>% select(any_of(KEEP[[comp]]))
  })
  parts <- Filter(Negate(is.null), parts)
  merged <- Reduce(function(a, b) left_join(a, b, by = "SEQN"), parts)

  rx <- read_nhanes("RXQ_RX", suffix, year, cycle)
  if (!is.null(rx)) {
    rx_flag <- rx %>%
      mutate(drug = toupper(as.character(RXDDRUG))) %>%
      group_by(SEQN) %>%
      summarise(osteo_med_raw = as.integer(any(grepl(OSTEO_DRUGS, drug))))
    merged <- left_join(merged, rx_flag, by = "SEQN")
  }
  merged$cycle <- cycle
  merged
}

raw <- pmap_dfr(CYCLES, read_cycle)
saveRDS(raw, out("raw_pooled.rds"))

# ---- 3. Harmonization --------------------------------------------------------
col_or_na <- function(d, v) if (v %in% names(d)) d[[v]] else rep(NA_real_, nrow(d))
yn  <- function(x) dplyr::case_when(x == 1 ~ "Yes", x == 2 ~ "No", TRUE ~ NA_character_)
any_yes <- function(...) {
  m <- cbind(...)
  valid <- !is.na(m) & (m == 1 | m == 2)
  out <- ifelse(rowSums(m == 1, na.rm = TRUE) > 0, "Yes", "No")
  out[rowSums(valid) == 0] <- NA
  out
}
fct_yn <- function(x) factor(x, levels = c("No", "Yes"))

d <- raw
for (v in unique(unlist(KEEP))) if (!v %in% names(d)) d[[v]] <- NA_real_
if (!"osteo_med_raw" %in% names(d)) d$osteo_med_raw <- NA_integer_

# 2017-March 2020 names -> the names used for the other cycles
#   MEC weight: WTMECPRP (3.2-year weight; rescaled in pool_weight())
#   Marital status: DMDMARTZ 1 married/living with partner, 2 widowed/divorced/separated,
#   3 never married -> DMDMARTL codes 1, 2, 5 (only partnered vs unpartnered is used)
d <- d %>% mutate(
  WTMEC    = coalesce(WTMEC2YR, WTMECPRP),
  DMDMARTL = ifelse(cycle == "2017-2020",
                    case_when(DMDMARTZ == 1 ~ 1, DMDMARTZ == 2 ~ 2, DMDMARTZ == 3 ~ 5), DMDMARTL)
)

# Alcohol: harmonize to drinks per week, then NIAAA sex-specific categories
alq_days_2017 <- c(`1` = 7, `2` = 6, `3` = 3.5, `4` = 2, `5` = 1, `6` = 2.5 / 4.345,
                   `7` = 1 / 4.345, `8` = 9 / 52.14, `9` = 4.5 / 52.14, `10` = 1.5 / 52.14)
d <- d %>% mutate(
  # 2005-2014 items                                                    # VERIFY
  alq_q   = ifelse(ALQ120Q %in% c(777, 999), NA, ALQ120Q),
  alq_n   = ifelse(ALQ130 %in% c(777, 999), NA, ALQ130),
  perwk_o = case_when(alq_q == 0 ~ 0,
                      ALQ120U == 1 ~ alq_q,
                      ALQ120U == 2 ~ alq_q / 4.345,
                      ALQ120U == 3 ~ alq_q / 52.14),
  never_o = ALQ101 %in% 2 & ALQ110 %in% 2,
  dpw_old = case_when(never_o ~ 0, perwk_o == 0 ~ 0, TRUE ~ perwk_o * alq_n),
  # 2017-2018 and 2017-March 2020 items (same questions and codes)    # VERIFY
  never_n = ALQ111 %in% 2,
  dpw_new = case_when(never_n ~ 0,
                      ALQ121 %in% 0 ~ 0,
                      ALQ121 %in% 1:10 ~ unname(alq_days_2017[as.character(ALQ121)]) * alq_n),
  drinks_wk = ifelse(cycle %in% c("2017-2018", "2017-2020"), dpw_new, dpw_old)
)

d <- d %>% mutate(
  sex       = factor(ifelse(RIAGENDR == 1, "Male", "Female"), levels = c("Male", "Female")),
  age       = pmin(RIDAGEYR, 80),                 # 2005-06 top-coded at 85, later at 80
  race      = factor(case_when(RIDRETH1 %in% c(1, 2) ~ "Hispanic",
                               RIDRETH1 == 3 ~ "Non-Hispanic White",
                               RIDRETH1 == 4 ~ "Non-Hispanic Black",
                               RIDRETH1 == 5 ~ "Other/Multiracial"),
                     levels = c("Non-Hispanic White", "Non-Hispanic Black",
                                "Hispanic", "Other/Multiracial")),
  education = factor(case_when(DMDEDUC2 %in% c(1, 2) ~ "Less than high school",
                               DMDEDUC2 == 3 ~ "High school/GED",
                               DMDEDUC2 == 4 ~ "Some college/associate",
                               DMDEDUC2 == 5 ~ "College graduate or above"),
                     levels = c("Less than high school", "High school/GED",
                                "Some college/associate", "College graduate or above")),
  pir       = INDFMPIR,
  pir_cat   = factor(case_when(pir < 1.3 ~ "Low (<1.30)",
                               pir < 3.5 ~ "Middle (1.30-3.49)",
                               pir >= 3.5 ~ "High (>=3.50)"),
                     levels = c("Low (<1.30)", "Middle (1.30-3.49)", "High (>=3.50)")),
  weight    = BMXWT,
  height    = BMXHT,
  smoking   = factor(case_when(SMQ020 == 2 ~ "Never",
                               SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current",
                               SMQ020 == 1 & SMQ040 == 3 ~ "Former"),
                     levels = c("Never", "Former", "Current")),
  alcohol   = factor(case_when(drinks_wk == 0 ~ "None in past year",
                               sex == "Male"   & drinks_wk <= 14 ~ "Moderate",
                               sex == "Female" & drinks_wk <= 7  ~ "Moderate",
                               drinks_wk > 0 ~ "Heavy"),
                     levels = c("None in past year", "Moderate", "Heavy")),
  phys_active = fct_yn(ifelse(cycle == "2005-2006",
                              any_yes(PAD200, PAD320),
                              any_yes(PAQ605, PAQ620, PAQ650, PAQ665))),
  # Arthritis type: MCQ190 (2005-08), MCQ191 (2009-10), MCQ195 (2013+).
  # Rheumatoid arthritis is code 1 up to 2009-10 but code 2 from 2013 on
  # (code 1 = osteoarthritis there) - checked against CDC codebooks.
  arth_type = coalesce(MCQ190, MCQ191, MCQ195),
  ra_code   = ifelse(cycle %in% c("2013-2014", "2017-2018", "2017-2020"), 2, 1),
  ra        = fct_yn(case_when(MCQ160A == 2 ~ "No",
                               MCQ160A == 1 & arth_type == ra_code ~ "Yes",
                               MCQ160A == 1 & arth_type %in% 1:4 ~ "No")),
  thyroid   = fct_yn(yn(MCQ160M)),                  # checked: thyroid problem, all cycles
  prior_fx  = fct_yn(any_yes(OSQ010A, OSQ010B, OSQ010C)),                # hip, wrist, spine
  glucocorticoid = fct_yn(yn(OSQ130)),
  parent_hipfx   = fct_yn(any_yes(OSQ170, OSQ200)),
  osteo_med = fct_yn(case_when(osteo_med_raw == 1 ~ "Yes", osteo_med_raw == 0 ~ "No")),
  # Supplements: from 30-day dietary supplement totals (NOT dietary intake)
  # 2005-2006 has no nutrient totals: supplement users are left missing (imputed)
  calcium_supp = fct_yn(case_when(DSDCOUNT == 0 ~ "No",
                                  cycle == "2005-2006" ~ NA_character_,
                                  coalesce(DSQTCALC, 0) > 0 ~ "Yes",
                                  !is.na(DSDCOUNT) ~ "No")),
  vitd_supp    = fct_yn(case_when(DSDCOUNT == 0 ~ "No",
                                  cycle == "2005-2006" ~ NA_character_,
                                  coalesce(DSQTVD, 0) > 0 ~ "Yes",
                                  !is.na(DSDCOUNT) ~ "No")),
  # Female-specific (NA for men)
  postmeno  = fct_yn(ifelse(sex == "Female", yn(RHQ031) %>% recode(Yes = "No", No = "Yes"), NA)),
  ever_ht   = fct_yn(ifelse(sex == "Female", yn(RHQ540), NA)),
  ever_preg = fct_yn(ifelse(sex == "Female", yn(RHQ131), NA)),
  # Extended SDoH
  food_insecure = fct_yn(case_when(FSDAD %in% c(1, 2) ~ "No", FSDAD %in% c(3, 4) ~ "Yes")),
  insured       = fct_yn(yn(HIQ011)),
  routine_care  = fct_yn(case_when(HUQ030 %in% c(1, 3) ~ "Yes", HUQ030 == 2 ~ "No")),
  partnered     = fct_yn(case_when(DMDMARTL %in% c(1, 6) ~ "Yes",
                                   DMDMARTL %in% 2:5 ~ "No")),
  cycle     = factor(cycle),
  fn_bmd    = DXXNKBMD,
  tscore    = (DXXNKBMD - CONFIG$ref_mean) / CONFIG$ref_sd,
  # WHO categories (3 levels) and two binary splits used for AUC
  bmd_cat       = factor(case_when(tscore <= -2.5 ~ "Osteoporosis",
                                   tscore <  -1.0 ~ "Low bone mass (osteopenia)",
                                   !is.na(tscore) ~ "Normal"),
                         levels = c("Normal", "Low bone mass (osteopenia)", "Osteoporosis")),
  osteoporosis  = as.integer(tscore <= -2.5),        # osteoporosis vs not
  lowbmd_or_op  = as.integer(tscore < -1.0),         # low bone mass OR osteoporosis vs normal
  # pooled MEC weight: WTMEC2YR / 5 for five 2-year cycles; with 2017-March 2020,
  # each weight x (years in cycle / 11.2 total years) - see nhanes_cycles.R
  wt_pool   = pool_weight(WTMEC, cycle)
)

# Analytic sample + flow counts (for Figure 1)
flow <- tibble(step = "NHANES participants, 5 cycles", n = nrow(d))
d1 <- d %>% filter(RIDAGEYR >= 50);                     flow <- add_row(flow, step = "Aged >= 50 years", n = nrow(d1))
d2 <- d1 %>% filter(!(RIDEXPRG %in% 1));                flow <- add_row(flow, step = "Not pregnant", n = nrow(d2))
d3 <- d2 %>% filter(!is.na(fn_bmd));                    flow <- add_row(flow, step = "Valid femoral neck DXA", n = nrow(d3))
d4 <- d3 %>% filter(!is.na(WTMEC) & WTMEC > 0);   flow <- add_row(flow, step = "Positive MEC weight (analytic sample)", n = nrow(d4))
flow <- bind_rows(flow, d4 %>% count(sex) %>% transmute(step = paste0("  ", sex), n),
                  d4 %>% count(cycle) %>% transmute(step = paste0("  ", cycle), n))
write_csv(flow, out("fig1_flow_counts.csv"))
ana <- d4

# ---- 4. Predictor sets -------------------------------------------------------
BASE_COMMON <- c("age", "weight", "height", "race", "cycle", "smoking", "alcohol",
                 "phys_active", "thyroid", "ra", "calcium_supp", "vitd_supp",
                 "prior_fx", "glucocorticoid", "parent_hipfx", "osteo_med")
WOMEN_EXTRA <- c("postmeno", "ever_ht", "ever_preg")
SDOH        <- c("education", "pir")
SDOH_EXT    <- c("food_insecure", "insured", "routine_care", "partnered")
IMPUTE_VARS <- unique(c(BASE_COMMON, WOMEN_EXTRA, SDOH, SDOH_EXT))

model_sets <- function(sex) {
  base <- c(BASE_COMMON, if (sex == "Female") WOMEN_EXTRA)
  list(base        = base,
       full        = c(base, SDOH),
       extended    = c(base, SDOH, SDOH_EXT),
       base_norace = setdiff(base, "race"),
       full_norace = setdiff(c(base, SDOH), "race"))
}
COMPARISONS <- list(
  "SDoH (education + income)"           = c("full", "base"),
  "Extended SDoH (all six)"             = c("extended", "base"),
  "Extended items beyond edu + income"  = c("extended", "full"),
  "SDoH, race excluded from both"       = c("full_norace", "base_norace")
)

# ---- 5. Table 1: survey-weighted descriptives --------------------------------
des_all <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pool,
                     nest = TRUE, data = ana)

describe_var <- function(v) {
  x <- ana[[v]]
  miss <- tapply(is.na(x), ana$sex, sum)
  if (is.numeric(x)) {
    f  <- as.formula(paste0("~", v))
    wsd <- function(s) { k <- ana$sex == s & !is.na(x); w <- ana$wt_pool[k]
      m <- sum(w * x[k]) / sum(w); c(m, sqrt(sum(w * (x[k] - m)^2) / sum(w))) }
    sm <- wsd("Male"); sf <- wsd("Female")
    p  <- tryCatch(svyttest(as.formula(paste0(v, "~sex")), des_all)$p.value,
                   error = function(e) NA)
    tibble(variable = v, level = "mean (SD)",
           Male   = sprintf("%.2f (%.2f)", sm[1], sm[2]),
           Female = sprintf("%.2f (%.2f)", sf[1], sf[2]),
           p_value = p, missing_M = miss[["Male"]], missing_F = miss[["Female"]])
  } else {
    sub <- subset(des_all, !is.na(ana[[v]]))
    tab <- svytable(as.formula(paste0("~", v, "+sex")), sub)
    pct <- prop.table(tab, 2) * 100
    n_u <- table(ana[[v]], ana$sex)[rownames(pct), , drop = FALSE]
    p <- if (v %in% WOMEN_EXTRA) NA else
      tryCatch(svychisq(as.formula(paste0("~", v, "+sex")), sub)$p.value,
               error = function(e) NA)
    tibble(variable = v, level = rownames(pct),
           Male   = sprintf("%d (%.1f%%)", n_u[, "Male"],   pct[, "Male"]),
           Female = sprintf("%d (%.1f%%)", n_u[, "Female"], pct[, "Female"]),
           p_value = c(p, rep(NA, nrow(pct) - 1)),
           missing_M = c(miss[["Male"]], rep(NA, nrow(pct) - 1)),
           missing_F = c(miss[["Female"]], rep(NA, nrow(pct) - 1)))
  }
}
t1_vars <- c("tscore", "fn_bmd", "bmd_cat", "age", "weight",
             "height", "race", "education", "pir_cat", SDOH_EXT, "smoking", "alcohol",
             "phys_active", "calcium_supp", "vitd_supp", "thyroid", "ra", "prior_fx",
             "glucocorticoid", "parent_hipfx", "osteo_med", WOMEN_EXTRA)
ana_t1 <- ana
ana_bak <- ana; ana <- ana_t1
des_all <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pool,
                     nest = TRUE, data = ana)
table1 <- map_dfr(t1_vars, describe_var)
ana <- ana_bak
write_csv(table1, out("table1_weighted_descriptives.csv"))

# ---- 6. Multiple imputation (outcome excluded) --------------------------------
impute_sex <- function(dsex, sex) {
  vars <- unique(c(IMPUTE_VARS[!(IMPUTE_VARS %in% WOMEN_EXTRA) | sex == "Female"]))
  X <- as.data.frame(dsex[, vars])
  meth <- make.method(X)
  imp <- mice(X, m = CONFIG$m_imp, method = meth, seed = CONFIG$seed,
              printFlag = FALSE, maxit = 10)
  imp
}
build_imputed <- function(imp, dsex, delta = 0) {
  lapply(seq_len(imp$m), function(i) {
    X <- complete(imp, i)
    if (delta != 0) {
      miss_pir <- is.na(dsex$pir)
      X$pir[miss_pir] <- pmin(pmax(X$pir[miss_pir] + delta, 0), 5)
    }
    bind_cols(dsex %>% select(SEQN, tscore, osteoporosis, lowbmd_or_op, wt_pool,
                              SDMVPSU, SDMVSTRA, sex),
              X)
  })
}

# ---- 7. Modelling helpers ----------------------------------------------------
wmetrics <- function(y, yhat, w) {
  w <- w / sum(w); e <- y - yhat; ybar <- sum(w * y)
  mse <- sum(w * e^2)
  c(MSE = mse, RMSE = sqrt(mse), MAE = sum(w * abs(e)),
    R2 = 1 - mse / sum(w * (y - ybar)^2))
}
wmse <- function(y, yhat, w) sum(w * (y - yhat)^2) / sum(w)

# Weighted AUC: P(risk of a case > risk of a non-case), ties count 1/2
wauc <- function(risk, case, w) {
  ok <- !is.na(risk) & !is.na(case)
  risk <- risk[ok]; case <- case[ok]; w <- w[ok]
  if (length(unique(case)) < 2) return(NA_real_)
  u <- sort(unique(risk)); g <- factor(match(risk, u), levels = seq_along(u))
  wn <- tapply(w * (case == 0), g, sum); wn[is.na(wn)] <- 0
  wp <- tapply(w * (case == 1), g, sum); wp[is.na(wp)] <- 0
  below <- cumsum(wn) - wn
  sum(wp * (below + 0.5 * wn)) / (sum(wp) * sum(wn))
}

# WHO 3-category agreement between observed and predicted T-score
who3 <- function(t) ifelse(t <= -2.5, 3L, ifelse(t < -1.0, 2L, 1L))
wacc3 <- function(y, yhat, w) { o <- who3(y); q <- who3(yhat); sum(w * (o == q)) / sum(w) }
wkappa3 <- function(y, yhat, w) {          # linear-weighted kappa, survey-weighted
  o <- factor(who3(y), 1:3); q <- factor(who3(yhat), 1:3)
  P <- xtabs(w ~ o + q) / sum(w)
  A <- 1 - abs(outer(1:3, 1:3, "-")) / 2
  po <- sum(P * A); pe <- sum(outer(rowSums(P), colSums(P)) * A)
  if (pe >= 1) return(NA_real_)
  (po - pe) / (1 - pe)
}

# Nadeau & Bengio (2003) corrected resampled t for repeated CV
nb_ci <- function(x, n_train, n_test, level = 0.95) {
  x <- x[!is.na(x)]; J <- length(x)
  if (J < 2) return(c(est = mean(x), lo = NA, hi = NA, p = NA))
  m <- mean(x); se <- sqrt((1 / J + n_test / n_train) * var(x))
  tq <- qt(1 - (1 - level) / 2, df = J - 1)
  c(est = m, lo = m - tq * se, hi = m + tq * se, p = 2 * pt(-abs(m / se), df = J - 1))
}

# Choose mtry / min.node.size by out-of-bag error (training data only)
tune_rf <- function(train, preds, w, seed) {
  f <- reformulate(preds, "tscore")
  errs <- vapply(seq_len(nrow(CONFIG$grid)), function(i) {
    ranger(f, data = train, num.trees = CONFIG$tune_trees,
           mtry = min(CONFIG$grid$mtry[i], length(preds)),
           min.node.size = CONFIG$grid$min.node.size[i],
           case.weights = w, respect.unordered.factors = "order",
           num.threads = CONFIG$threads, seed = seed)$prediction.error
  }, numeric(1))
  as.list(CONFIG$grid[which.min(errs), ])
}
fit_rf <- function(train, preds, w, seed, params = NULL) {
  if (is.null(params)) params <- tune_rf(train, preds, w, seed)
  ranger(reformulate(preds, "tscore"), data = train, num.trees = CONFIG$num_trees,
         mtry = min(params$mtry, length(preds)),
         min.node.size = params$min.node.size,
         case.weights = w, respect.unordered.factors = "order",
         num.threads = CONFIG$threads, seed = seed)
}

lm_formula <- function(preds) {
  terms <- vapply(preds, function(p) {
    if (p %in% c("age", "weight", "height", "pir")) sprintf("ns(%s, df = 3)", p) else p
  }, character(1))
  as.formula(paste("tscore ~", paste(terms, collapse = " + ")))
}

perm_importance <- function(model, test, w, groups) {
  y <- test$tscore
  base <- wmse(y, predict(model, test, num.threads = CONFIG$threads)$predictions, w)
  vapply(groups, function(vars) {
    mean(replicate(CONFIG$n_perm, {
      tp <- test; idx <- sample(nrow(tp))
      tp[vars] <- tp[idx, vars, drop = FALSE]
      wmse(y, predict(model, tp, num.threads = CONFIG$threads)$predictions, w) - base
    }))
  }, numeric(1))
}

make_folds <- function(n) {
  folds <- list()
  for (r in seq_len(CONFIG$n_repeats)) {
    id <- sample(rep(seq_len(CONFIG$k_folds), length.out = n))
    for (k in seq_len(CONFIG$k_folds)) folds[[length(folds) + 1]] <-
      list(rep = r, fold = k, test = which(id == k))
  }
  folds
}

# One repeated-CV run for one sex. dat_list = list of imputed data frames
# Each finished fold is saved to results/checkpoints/, so if the session ends
# you can simply run the script again and it continues where it stopped.
run_cv <- function(dat_list, sex, sets = names(model_sets(sex)),
                   do_importance = TRUE, do_lm = TRUE, label = "primary") {
  ms <- model_sets(sex)[sets]
  set.seed(CONFIG$seed + sum(utf8ToInt(paste(label, sex))))   # same folds on every restart
  folds <- make_folds(nrow(dat_list[[1]]))
  imp_groups <- c(as.list(setNames(ms$full, ms$full)),
                  list(`SDoH group (education + income)` = SDOH))
  tune_set <- if ("full" %in% names(ms)) ms$full else ms[[length(ms)]]
  ck_dir <- out(file.path("checkpoints",
                          gsub("[^A-Za-z0-9]+", "_", paste(if (CONFIG$quick) "quick" else "full",
                                                           label, sex))))
  dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)
  res_list <- vector("list", length(folds)); t0 <- Sys.time(); n_new <- 0
  for (fi in seq_along(folds)) {
    ck <- file.path(ck_dir, sprintf("fold_%03d.rds", fi))
    if (file.exists(ck)) { res_list[[fi]] <- readRDS(ck); next }
    te <- folds[[fi]]$test; perf <- list(); imps <- list()
    # tune once per fold (imputation 1, full model) and reuse for every model
    tr1 <- dat_list[[1]][-te, ]
    params <- tune_rf(tr1, tune_set, tr1$wt_pool / mean(tr1$wt_pool), CONFIG$seed + fi)
    for (mi in seq_along(dat_list)) {
      dd <- dat_list[[mi]]; tr <- dd[-te, ]; ts <- dd[te, ]
      w_tr <- tr$wt_pool / mean(tr$wt_pool); w_ts <- ts$wt_pool
      for (mn in names(ms)) {
        rf <- fit_rf(tr, ms[[mn]], w_tr, seed = CONFIG$seed + fi, params = params)
        p  <- predict(rf, ts, num.threads = CONFIG$threads)$predictions
        perf[[length(perf) + 1]] <- tibble(
          sex = sex, analysis = label, rep = folds[[fi]]$rep, fold = folds[[fi]]$fold,
          imp = mi, learner = "Random forest", model = mn,
          !!!as.list(wmetrics(ts$tscore, p, w_ts)),
          AUC_osteoporosis = wauc(-p, ts$osteoporosis, w_ts),
          AUC_lowBMD_or_OP  = wauc(-p, ts$lowbmd_or_op, w_ts),
          Acc_WHO3 = wacc3(ts$tscore, p, w_ts),
          Kappa_WHO3 = wkappa3(ts$tscore, p, w_ts),
          n_train = nrow(tr), n_test = nrow(ts))
        # permutation importance on imputation 1 only (the slowest step)
        if (do_importance && mn == "full" && mi == 1) {
          vi <- perm_importance(rf, ts, w_ts, imp_groups)
          imps[[length(imps) + 1]] <- tibble(sex = sex, analysis = label,
            rep = folds[[fi]]$rep, fold = folds[[fi]]$fold, imp = mi,
            variable = names(vi), dMSE = unname(vi))
        }
        if (do_lm) {
          tr$.w <- w_tr   # weights must live in the data for lm()
          lmf <- tryCatch(lm(lm_formula(ms[[mn]]), data = tr, weights = .w),
                          error = function(e) NULL)
          if (!is.null(lmf)) {
            pl <- tryCatch(predict(lmf, ts), error = function(e) rep(NA, nrow(ts)))
            perf[[length(perf) + 1]] <- tibble(
              sex = sex, analysis = label, rep = folds[[fi]]$rep, fold = folds[[fi]]$fold,
              imp = mi, learner = "Linear (splines)", model = mn,
              !!!as.list(wmetrics(ts$tscore, pl, w_ts)),
              AUC_osteoporosis = wauc(-pl, ts$osteoporosis, w_ts),
              AUC_lowBMD_or_OP  = wauc(-pl, ts$lowbmd_or_op, w_ts),
              Acc_WHO3 = wacc3(ts$tscore, pl, w_ts),
              Kappa_WHO3 = wkappa3(ts$tscore, pl, w_ts),
              n_train = nrow(tr), n_test = nrow(ts))
          }
        }
      }
    }
    res_list[[fi]] <- list(perf = bind_rows(perf), imp = bind_rows(imps))
    saveRDS(res_list[[fi]], ck)
    n_new <- n_new + 1
    mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("  %s | %s: fold %d / %d done  (%.1f min so far, about %.0f min left)",
                    label, sex, fi, length(folds), mins,
                    (length(folds) - fi) * mins / n_new))
  }
  list(perf = bind_rows(lapply(res_list, `[[`, "perf")),
       imp  = bind_rows(lapply(res_list, `[[`, "imp")))
}

# Average over imputations within each fold, then summarize across folds
summarise_perf <- function(perf) {
  per_fold <- perf %>%
    group_by(sex, analysis, learner, model, rep, fold) %>%
    summarise(across(c(MSE, RMSE, MAE, R2, AUC_osteoporosis, AUC_lowBMD_or_OP, Acc_WHO3, Kappa_WHO3,
                       n_train, n_test), mean), .groups = "drop")
  overall <- per_fold %>%
    pivot_longer(c(MSE, RMSE, MAE, R2, AUC_osteoporosis, AUC_lowBMD_or_OP, Acc_WHO3, Kappa_WHO3),
                 names_to = "metric") %>%
    group_by(sex, analysis, learner, model, metric) %>%
    summarise(ci = list(nb_ci(value, mean(n_train), mean(n_test))), .groups = "drop") %>%
    mutate(est = map_dbl(ci, ~ .x[["est"]]), lo = map_dbl(ci, ~ .x[["lo"]]), hi = map_dbl(ci, ~ .x[["hi"]])) %>%
    select(-ci)
  deltas <- map_dfr(names(COMPARISONS), function(cn) {
    a <- COMPARISONS[[cn]][1]; b <- COMPARISONS[[cn]][2]
    if (!all(c(a, b) %in% per_fold$model)) return(NULL)   # e.g. sensitivity runs only base/full
    per_fold %>% filter(model %in% c(a, b)) %>%
      pivot_longer(c(MSE, RMSE, MAE, R2, AUC_osteoporosis, AUC_lowBMD_or_OP, Acc_WHO3, Kappa_WHO3),
                   names_to = "metric") %>%
      select(sex, analysis, learner, rep, fold, metric, model, value, n_train, n_test) %>%
      pivot_wider(names_from = model, values_from = value) %>%
      filter(!is.na(.data[[a]]), !is.na(.data[[b]])) %>%
      mutate(diff = .data[[a]] - .data[[b]]) %>%
      group_by(sex, analysis, learner, metric) %>%
      summarise(ci = list(nb_ci(diff, mean(n_train), mean(n_test))), .groups = "drop") %>%
      mutate(comparison = cn, est = map_dbl(ci, ~ .x[["est"]]), lo = map_dbl(ci, ~ .x[["lo"]]),
             hi = map_dbl(ci, ~ .x[["hi"]]), p = map_dbl(ci, ~ .x[["p"]])) %>% select(-ci)
  })
  list(per_fold = per_fold, overall = overall, deltas = deltas)
}

summarise_importance <- function(imp, n_train, n_test) {
  per_fold <- imp %>% group_by(sex, analysis, rep, fold, variable) %>%
    summarise(dMSE = mean(dMSE), .groups = "drop") %>%
    group_by(sex, analysis, rep, fold) %>%
    mutate(is_grp = variable == "SDoH group (education + income)",
           rank = rank(-ifelse(is_grp, -Inf, dMSE), ties.method = "average"),
           rank = ifelse(is_grp, NA, rank)) %>% ungroup()
  per_fold %>% group_by(sex, analysis, variable) %>%
    summarise(ci = list(nb_ci(dMSE, n_train, n_test)),
              median_rank = median(rank, na.rm = TRUE),
              rank_q1 = suppressWarnings(quantile(rank, .25, na.rm = TRUE)),
              rank_q3 = suppressWarnings(quantile(rank, .75, na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(dMSE = map_dbl(ci, ~ .x[["est"]]), lo = map_dbl(ci, ~ .x[["lo"]]), hi = map_dbl(ci, ~ .x[["hi"]])) %>%
    select(-ci) %>% arrange(sex, desc(dMSE))
}

# ---- 8. Primary analysis -----------------------------------------------------
results <- list(); imps_all <- list(); imputations <- list()
for (sx in c("Male", "Female")) {
  dsx <- ana %>% filter(sex == sx)
  imp_file <- out(sprintf("imputation_%s%s.rds", sx, if (CONFIG$quick) "_quick" else ""))
  if (file.exists(imp_file)) {
    message("Loading saved imputation: ", sx); imp <- readRDS(imp_file)
  } else {
    message("Imputing: ", sx, " (about 5-15 minutes, no messages while it runs)")
    imp <- impute_sex(dsx, sx); saveRDS(imp, imp_file)
  }
  imputations[[sx]] <- list(imp = imp, data = dsx)
  dl  <- build_imputed(imp, dsx)
  message("Primary CV: ", sx)
  res <- run_cv(dl, sx)
  results[[length(results) + 1]] <- res$perf
  imps_all[[length(imps_all) + 1]] <- res$imp
}
saveRDS(imputations, out("imputations.rds"))
perf_primary <- bind_rows(results)
imp_primary  <- bind_rows(imps_all)
write_csv(perf_primary, out("cv_raw_performance_primary.csv"))

S <- summarise_perf(perf_primary)
write_csv(S$overall, out("table_performance_by_model.csv"))
write_csv(S$deltas,  out("table_added_value_deltas.csv"))
ntr <- mean(perf_primary$n_train); nts <- mean(perf_primary$n_test)
VI <- summarise_importance(imp_primary, ntr, nts)
write_csv(VI, out("table2_permutation_importance.csv"))

# ---- 9. ALE plots (full RF model on imputation 1, per sex) -------------------
ale_list <- list()
for (sx in c("Male", "Female")) {
  dd <- build_imputed(imputations[[sx]]$imp, imputations[[sx]]$data)[[1]]
  preds <- model_sets(sx)$full
  rf <- fit_rf(dd, preds, dd$wt_pool / mean(dd$wt_pool), seed = CONFIG$seed)
  pr <- Predictor$new(rf, data = as.data.frame(dd[, preds]), y = dd$tscore,
                      predict.function = function(m, newdata) predict(m, newdata, num.threads = CONFIG$threads)$predictions)
  for (ft in c("pir", "education")) {
    fe <- FeatureEffect$new(pr, feature = ft, method = "ale", grid.size = 20)
    ale_list[[length(ale_list) + 1]] <- fe$results %>%
      rename(x = all_of(ft), ale = .value) %>%
      mutate(x = as.character(x), feature = ft, sex = sx)
  }
}
ale_df <- bind_rows(ale_list)
write_csv(ale_df, out("fig3_ale_values.csv"))
p_ale <- ggplot(ale_df %>% filter(feature == "pir"),
                aes(as.numeric(x), ale, colour = sex)) +
  geom_line(linewidth = 0.9) +
  labs(x = "Family poverty-income ratio", y = "ALE on predicted T-score", colour = NULL) +
  theme_minimal(base_size = 11)
ggsave(out("fig3a_ale_income.png"), p_ale, width = 6, height = 4, dpi = 300)
p_ale2 <- ggplot(ale_df %>% filter(feature == "education"),
                 aes(factor(x, levels = levels(ana$education)), ale, colour = sex, group = sex)) +
  geom_point(size = 2) + geom_line() +
  labs(x = "Education", y = "ALE on predicted T-score", colour = NULL) +
  theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(out("fig3b_ale_education.png"), p_ale2, width = 6, height = 4, dpi = 300)

# ---- 10. Sensitivity analyses ------------------------------------------------
# Only the SDoH contrast (full vs base) with RF + linear, no importance, to save time
sens <- list()
for (sx in c("Male", "Female")) {
  dsx <- imputations[[sx]]$data; imp <- imputations[[sx]]$imp
  ms  <- model_sets(sx)
  # S1 complete case
  cc <- dsx %>% select(SEQN, tscore, osteoporosis, lowbmd_or_op, wt_pool,
                       SDMVPSU, SDMVSTRA, sex, all_of(ms$full)) %>% tidyr::drop_na()
  sens[[length(sens) + 1]] <- run_cv(list(cc), sx, sets = c("base", "full"),
                                     do_importance = FALSE, label = "S1 complete case")$perf
  # S2 delta-adjusted MNAR for income
  for (dl in CONFIG$mnar_deltas) {
    sens[[length(sens) + 1]] <- run_cv(build_imputed(imp, dsx, delta = dl), sx,
                                       sets = c("base", "full"), do_importance = FALSE,
                                       label = sprintf("S2 MNAR delta %+.1f PIR", dl))$perf
  }
  # S3 recent cycles only (2013-2014 and 2017-2018 or 2017-March 2020)
  keep <- dsx$cycle %in% CYCLES_2013ON
  dl3 <- lapply(build_imputed(imp, dsx), function(x) {
    x <- x[keep, ]; x$cycle <- droplevels(x$cycle); x })
  sens[[length(sens) + 1]] <- run_cv(dl3, sx, sets = c("base", "full"),
                                     do_importance = FALSE,
                                     label = if (USE_P) "S3 2013-2020 only" else "S3 2013-2018 only")$perf
  # S4 exclude osteoporosis-medication users
  dl4 <- lapply(build_imputed(imp, dsx), function(x) x[x$osteo_med == "No", ])
  sens[[length(sens) + 1]] <- run_cv(dl4, sx, sets = c("base", "full"),
                                     do_importance = FALSE, label = "S4 no osteoporosis medication")$perf
}
perf_sens <- bind_rows(sens)
write_csv(perf_sens, out("cv_raw_performance_sensitivity.csv"))
S_sens <- summarise_perf(perf_sens)
write_csv(S_sens$deltas %>% filter(comparison == "SDoH (education + income)"),
          out("table_sensitivity_deltas.csv"))

# ---- 11. Associational model: svyglm with sex interactions ------------------
# Note: design built on the analytic subset (nest = TRUE); for strict domain
# estimation, define the design on all examined participants and subset().
assoc_vars <- c(BASE_COMMON, SDOH)
imp_combined <- lapply(seq_len(CONFIG$m_imp), function(i) {
  bind_rows(build_imputed(imputations$Male$imp,   imputations$Male$data)[[i]],
            build_imputed(imputations$Female$imp, imputations$Female$data)[[i]]) %>%
    mutate(pir_cat = cut(pir, c(-Inf, 1.3, 3.5, Inf), right = FALSE,
                         labels = levels(ana$pir_cat)))
})
des_mi <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pool,
                    nest = TRUE, data = imputationList(imp_combined))
f_assoc <- as.formula(paste("tscore ~ sex * (education + pir_cat) +",
                            paste(setdiff(BASE_COMMON, character(0)), collapse = " + ")))
fits <- with(des_mi, svyglm(f_assoc))
pooled <- MIcombine(fits)
assoc <- tibble(term = names(coef(pooled)), est = coef(pooled),
                se = sqrt(diag(vcov(pooled)))) %>%
  mutate(lo = est - 1.96 * se, hi = est + 1.96 * se, p = 2 * pnorm(-abs(est / se)))
write_csv(assoc, out("table_svyglm_sex_interactions.csv"))

# Joint Wald test of the sex x SDoH interaction terms (pooled over imputations)
int_terms <- grep("^sexFemale:", names(coef(pooled)), value = TRUE)
b <- coef(pooled)[int_terms]; V <- vcov(pooled)[int_terms, int_terms]
wald <- as.numeric(t(b) %*% solve(V) %*% b)
write_csv(tibble(test = "sex x (education + income) interactions",
                 chisq = wald, df = length(int_terms),
                 p = pchisq(wald, length(int_terms), lower.tail = FALSE)),
          out("table_interaction_wald.csv"))

# Missingness of PIR vs observed covariates (MAR-consistent check)
miss_tab <- ana %>% mutate(pir_missing = is.na(pir)) %>%
  group_by(sex, education, pir_missing) %>% summarise(n = n(), .groups = "drop")
write_csv(miss_tab, out("table_pir_missingness_by_education.csv"))

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("Done. Results in ", normalizePath(CONFIG$out_dir))
