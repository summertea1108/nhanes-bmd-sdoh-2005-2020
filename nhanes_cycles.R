# =============================================================================
#  nhanes_cycles.R - which NHANES releases the analysis uses, and how the
#  cycle-specific variables are read. Sourced by every analysis script.
#
#  NHANES_LATEST (environment variable) picks the last release in the series:
#    "P" (default)  2017-March 2020 pre-pandemic file (P_ prefix). It REPLACES
#                   2017-2018: the 2017-2018 participants are inside it again,
#                   under new SEQNs, so J and P must never be combined.
#    "J"            2017-2018, as in the submitted manuscript. Reproduces the
#                   submitted results exactly (same samples, seeds, folders).
#
#  Results go to results/ for "J" and to results_2017_2020/ for "P", so the
#  submitted results are never overwritten.
#
#  What differs in the 2017-March 2020 file (checked against the CDC codebooks
#  and the downloaded P_ files, Oct 2026):
#    * MEC weight is WTMECPRP (not WTMEC2YR). It covers 3.2 years, so when the
#      file is pooled with the earlier 2-year cycles each cycle's weight is
#      multiplied by (years in that cycle / total years) - CDC's guidance for
#      combining the pre-pandemic file with earlier cycles.
#    * Marital status is DMDMARTZ with 3 groups (1 married/living with partner,
#      2 widowed/divorced/separated, 3 never married). It is mapped onto the
#      DMDMARTL codes used elsewhere (1 -> 1, 2 -> 2, 3 -> 5); only the
#      partnered / unpartnered split is used, so nothing is lost.
#    * No household size (DMDHHSIZ): "unpartnered or living alone" can only use
#      marital status in 2017-2020. In 2017-2018 only 16 partnered adults
#      aged 50+ lived alone (34.1% vs 33.6% "alone" with and without it).
#    * No type of usual place of care (HUQ041): "poor access" is "no usual
#      place of care" only in 2017-2020 (the ER-as-usual-place part is lost;
#      9.3% vs 10.9% in 2017-2018).
#    * No housing questionnaire (no P_HOQ file): home tenure (renter) is not
#      available for 2017-2020. See SCORE_DROP below.
#    * Employment: OCQ380 is replaced by OCD383 (1 home/family, 2 school,
#      3 retired, 4 unable to work for health reasons/disabled, 7 other); the
#      "on layoff" category no longer exists.
#    * Femur DXA is limited to ages 50+ (same as our analytic sample).
#      2019-2020 scans used Hologic Horizon A scanners (cross-calibrated; all
#      2017-2020 scans analyzed with APEX 4.0).
# =============================================================================
NHANES_LATEST <- toupper(Sys.getenv("NHANES_LATEST", "P"))
if (!NHANES_LATEST %in% c("J", "P")) stop("NHANES_LATEST must be J (2017-2018) or P (2017-March 2020)")
USE_P <- NHANES_LATEST == "P"

CYC <- tibble::tribble(~suffix, ~year, ~cycle, ~years,
  "D", 2005, "2005-2006", 2,
  "E", 2007, "2007-2008", 2,
  "F", 2009, "2009-2010", 2,
  "H", 2013, "2013-2014", 2)
CYC <- if (USE_P) dplyr::bind_rows(CYC, tibble::tibble(suffix = "P", year = 2017, cycle = "2017-2020", years = 3.2)) else
                  dplyr::bind_rows(CYC, tibble::tibble(suffix = "J", year = 2017, cycle = "2017-2018", years = 2))
LATEST_CYCLE  <- CYC$cycle[nrow(CYC)]
CYCLES_2013ON <- c("2013-2014", LATEST_CYCLE)          # the "recent cycles only" sensitivity analysis
TOTAL_YEARS   <- sum(CYC$years)
RESULTS_ROOT  <- if (USE_P) "results_2017_2020" else "results"
EXPECTED_N    <- if (USE_P) NA_integer_ else 10981L      # analytic sample size check (J only)

# Pooled weight: each cycle's MEC weight x (years in cycle / total years).
# With five 2-year cycles this is WTMEC2YR / 5, written exactly as in the submitted code.
pool_weight <- function(wtmec, cycle) {
  if (all(CYC$years == 2)) return(wtmec / nrow(CYC))
  wtmec * CYC$years[match(as.character(cycle), CYC$cycle)] / TOTAL_YEARS
}

# Scored SDoH items that cannot be built in every cycle are dropped from the score and the block
SCORE_DROP <- if (USE_P) "renter" else character(0)
drop_unavailable <- function(x) setdiff(x, SCORE_DROP)

# File names: DEMO_J.xpt for 2-year cycles, P_DEMO.xpt for the pre-pandemic file
nhanes_file <- function(comp, suffix) if (suffix == "P") paste0("P_", comp) else paste0(comp, "_", suffix)
NOT_RELEASED <- list(P = c("HOQ"))                      # components with no 2017-March 2020 file
is_released <- function(comp, suffix) !(comp %in% NOT_RELEASED[[suffix]])

read_xpt_cached <- function(comp, s, y) {
  fname <- nhanes_file(comp, s)
  path <- file.path("nhanes_cache", paste0(fname, ".xpt"))
  if (!file.exists(path)) {
    u <- sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s.xpt", y, fname)
    ok <- tryCatch({ utils::download.file(u, path, mode = "wb", quiet = TRUE)
                      identical(rawToChar(readBin(path, "raw", 6)), "HEADER") },
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) { unlink(path)
      stop("Missing ", basename(path), ". Download it from ", u, " and upload it to nhanes_cache/.") }
  }
  haven::zap_labels(haven::read_xpt(path))
}
pick <- function(d, vars) { for (v in setdiff(vars, names(d))) d[[v]] <- NA_real_; d[, c("SEQN", vars)] }

# Extra SDoH variables used by sdoh_bundle.R, sdoh_v2_prediction.R and sdoh_v2_sensitivity.R
# (checked against CDC codebooks, Oct 2026). For J the result is identical to the submitted code.
build_extra <- function() {
  yesno <- function(cond_yes, cond_no) factor(dplyr::case_when(cond_yes ~ "Yes", cond_no ~ "No"), levels = c("No", "Yes"))
  extra <- purrr::pmap_dfr(CYC[, c("suffix", "year", "cycle")], function(suffix, year, cycle) {
    demo <- pick(read_xpt_cached("DEMO", suffix, year),
                 c("DMDHHSIZ", "DMDBORN", "DMDBORN2", "DMDBORN4", "SIALANG"))
    fsq  <- pick(read_xpt_cached("FSQ", suffix, year), c("FSQ165", "FSQ170", "FSQ171", "FSQ012"))
    huq  <- pick(read_xpt_cached("HUQ", suffix, year), c("HUQ040", "HUQ041"))
    hoq  <- if (is_released("HOQ", suffix)) pick(read_xpt_cached("HOQ", suffix, year), "HOQ065") else
              tibble::tibble(SEQN = demo$SEQN, HOQ065 = NA_real_)
    ocq  <- pick(read_xpt_cached("OCQ", suffix, year), c("OCD150", "OCQ380", "OCD383"))
    demo %>% dplyr::left_join(fsq, by = "SEQN") %>% dplyr::left_join(huq, by = "SEQN") %>%
      dplyr::left_join(hoq, by = "SEQN") %>% dplyr::left_join(ocq, by = "SEQN") %>% dplyr::mutate(cycle_x = cycle)
  })
  extra %>% dplyr::mutate(
    # SNAP in past 12 months (household): FSQ170 "authorized" (2005-06, asked of all), FSQ171 (2007-10)
    # and FSQ012 (2013 on) are asked ONLY if FSQ165 (ever received) = 1, so FSQ165 = 2 means "No"
    snap_raw = dplyr::case_when(cycle_x == "2005-2006" ~ FSQ170,
                                FSQ165 == 2 ~ 2,
                                cycle_x %in% c("2007-2008", "2009-2010") ~ FSQ171,
                                TRUE ~ FSQ012),
    snap = yesno(snap_raw == 1, snap_raw == 2),
    # Type of usual place: HUQ040 (2005-10, 7/9 missing) / HUQ041 (2013-18, 77/99 missing); ER = 3 in both.
    # Not asked in 2017-March 2020 (stays missing, so poor access = no usual place only).
    care_type = ifelse(cycle_x %in% c("2013-2014", "2017-2018", "2017-2020"), HUQ041, HUQ040),
    care_type = ifelse(care_type %in% c(7, 9, 77, 99), NA, care_type),
    # Home tenure: 1 owned/being bought; 2 rented; 3 other arrangement (no file for 2017-March 2020)
    renter = yesno(HOQ065 %in% c(2, 3), HOQ065 == 1),
    # Employment (sensitivity only): looking for work, laid off, or unable to work (health/disability);
    # working, retired, homemaker, student or "other" = No.
    # 2005-2018: OCQ380 (4 unable to work for health reasons, 5 on layoff, 6 disabled);
    # 2017-2020: OCD383 (4 unable to work for health reasons/disabled; no layoff category)
    unemployed = yesno(OCD150 == 3 | (OCD150 == 4 & (OCQ380 %in% c(4, 5, 6) | OCD383 %in% 4)),
                       OCD150 %in% c(1, 2) | (OCD150 == 4 & (OCQ380 %in% c(1, 2, 3, 7) | OCD383 %in% c(1, 2, 3, 7)))),
    # Context (not scored): US-born vs not (the only coding comparable across cycles)
    born_raw = dplyr::case_when(cycle_x == "2005-2006" ~ DMDBORN,
                                cycle_x %in% c("2007-2008", "2009-2010") ~ DMDBORN2,
                                TRUE ~ DMDBORN4),
    foreign_born = yesno(born_raw %in% c(2, 3, 4, 5), born_raw == 1),
    spanish_interview = yesno(SIALANG == 2, SIALANG == 1),
    lives_alone_raw = DMDHHSIZ == 1          # missing in 2017-2020 (not released)
  )
}
message(sprintf("NHANES cycles: %s | results in %s/%s", paste(CYC$cycle, collapse = ", "), RESULTS_ROOT,
                if (length(SCORE_DROP)) paste0(" | not available in every cycle, so dropped: ", paste(SCORE_DROP, collapse = ", ")) else ""))
