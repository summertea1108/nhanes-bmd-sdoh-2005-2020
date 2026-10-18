# =============================================================================
#  SDoH bundle study (follow-up to the education/income paper)
#
#  Aim 1  Association: cumulative SDoH disadvantage score (0-8; 0-7 with 2017-March 2020,
#         where home tenure is not available) vs femoral neck
#         T-score and osteoporosis, progressively adjusted; interactions with
#         sex and with race and ethnicity.
#  Aim 2  Prediction: does the whole SDoH block (or the score) improve
#         cross-validated prediction beyond clinical factors?
#  Aim 3  Race replacement: how much of the accuracy lost when race and
#         ethnicity are removed is recovered by adding the SDoH block?
#         Plus: how race coefficients change when SDoH and body size are added.
#
#  Needs (already on OSC from the main analysis):
#     <results folder>/imputations.rds  (analytic sample, harmonized variables)
#     <results folder> = results/ (NHANES_LATEST=J) or results_2017_2020/ (default, P)
#     nhanes_cache/DEMO_*, FSQ_*, HUQ_* (.xpt, 5 cycles)
#  NEW files needed in nhanes_cache/ (10 files) - downloaded automatically if
#  the session has internet; otherwise download them from the CDC website and
#  upload them to nhanes_cache/:
#     OCQ_D OCQ_E OCQ_F OCQ_H OCQ_J  HOQ_D HOQ_E HOQ_F HOQ_H HOQ_J  (.xpt)
#
#  Run in RStudio on OSC:
#     setwd("/users/PAS2982/summertea1109/Random Forest ML/NHANES_BMD")
#     Sys.setenv(QUICK = "TRUE"); source("sdoh_bundle.R")    # 10-20 min test first
#     Sys.setenv(QUICK = "FALSE"); source("sdoh_bundle.R")   # full run (several hours)
#  The cross-validation saves every fold, so if the session stops, run
#  source("sdoh_bundle.R") again and it resumes. Running the QUICK test twice
#  also checks that resuming works (the second run should load every fold).
#  The script STOPS after the harmonization check if any item has >20% missing
#  in any cycle; inspect 00_harmonization_check.csv (override only if you are
#  sure: Sys.setenv(SDOH_IGNORE_CHECK = "TRUE")).
#
#  Outputs go to <results folder>/sdoh_bundle/ (QUICK runs: <results folder>/sdoh_bundle_quick/).
#  FIRST look at 00_harmonization_check.csv: every item should have data in
#  every cycle with plausible "% yes".
# =============================================================================
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(purrr); library(readr)
  library(survey); library(mitools); library(mice); library(ranger); library(splines)
})
options(survey.lonely.psu = "adjust", dplyr.summarise.inform = FALSE)
PROJECT_DIR <- Sys.getenv("NHANES_BMD_DIR", "/users/PAS2982/summertea1109/Random Forest ML/NHANES_BMD")
if (dir.exists(PROJECT_DIR)) setwd(PROJECT_DIR)
source("nhanes_cycles.R")   # cycles, file names, results folder, extra SDoH variables

detect_threads <- function() {
  for (v in c("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_JOB_CPUS_PER_NODE")) {
    n <- suppressWarnings(as.integer(sub("\\(.*", "", Sys.getenv(v))))
    if (!is.na(n) && n > 0) return(n)
  }
  n <- tryCatch(length(parallel::mcaffinity()), error = function(e) NA_integer_)
  if (!is.na(n) && n > 0) return(n)
  max(1L, min(4L, parallel::detectCores() - 1L))
}
QUICK <- toupper(Sys.getenv("QUICK", "FALSE")) == "TRUE"
# SDOH_COMMON=TRUE: sensitivity analysis with definitions common to all survey periods
#   no usual place of care = HUQ030 only (emergency room as usual place not counted, 2005-2014)
#   not partnered          = marital status only (living alone not counted, 2005-2014)
COMMON <- toupper(Sys.getenv("SDOH_COMMON", "FALSE")) == "TRUE"
CFG <- list(seed = 20261009, k_folds = 10, n_repeats = 10, num_trees = 500, tune_trees = 300,
            grid = expand.grid(mtry = c(2, 4, 6), min.node.size = c(10, 50)),
            m_pred = 5, m_assoc = 20, threads = detect_threads())
if (QUICK) CFG <- modifyList(CFG, list(k_folds = 3, n_repeats = 2, num_trees = 100, tune_trees = 100,
                                       grid = expand.grid(mtry = c(2, 4), min.node.size = 20),
                                       m_pred = 2, m_assoc = 3))
OUT_DIR <- file.path(RESULTS_ROOT, paste0(if (QUICK) "sdoh_bundle_quick" else "sdoh_bundle", if (COMMON) "_common" else ""))
if (COMMON) message("SDOH_COMMON: common item definitions in all survey periods -> ", OUT_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)
message(sprintf("SDoH bundle | %s | cores: %d | %d x %d CV | m = %d / %d",
                if (QUICK) "QUICK TEST" else "FULL RUN", CFG$threads, CFG$k_folds,
                CFG$n_repeats, CFG$m_pred, CFG$m_assoc))

# =============================================================================
# 1. Read the extra SDoH variables (checked against CDC codebooks, Oct 2026)
# =============================================================================
# Cycles, file names and the extra variables come from nhanes_cycles.R
# (2017-March 2020 by default; NHANES_LATEST=J reproduces the submitted 2017-2018 analysis)
extra <- build_extra()
yesno <- function(cond_yes, cond_no) factor(case_when(cond_yes ~ "Yes", cond_no ~ "No"), levels = c("No", "Yes"))

# =============================================================================
# 2. Join to the analytic sample and define the scored items
# =============================================================================
imps <- readRDS(file.path(RESULTS_ROOT, "imputations.rds"))
BASE_COMMON <- c("age", "weight", "height", "race", "cycle", "smoking", "alcohol",
                 "phys_active", "thyroid", "ra", "calcium_supp", "vitd_supp",
                 "prior_fx", "glucocorticoid", "parent_hipfx", "osteo_med")
WOMEN_EXTRA <- c("postmeno", "ever_ht", "ever_preg")

prep <- function(sx) {
  imps[[sx]]$data %>%
    left_join(extra %>% select(SEQN, snap, care_type, renter, unemployed, foreign_born,
                               spanish_interview, lives_alone_raw), by = "SEQN") %>%
    mutate(
      uninsured   = yesno(HIQ011 == 2, HIQ011 == 1),
      poor_access = if (COMMON) yesno(HUQ030 == 2, HUQ030 %in% c(1, 3)) else
                    yesno(HUQ030 == 2 | care_type %in% 3, HUQ030 %in% c(1, 3)),
      # Unpartnered or living alone (DMDMARTL: 1 married, 6 living with partner)
      alone       = if (COMMON) yesno(DMDMARTL %in% 2:5, DMDMARTL %in% c(1, 6)) else
                    yesno(DMDMARTL %in% 2:5 | lives_alone_raw %in% TRUE,
                          DMDMARTL %in% c(1, 6) & !(lives_alone_raw %in% TRUE)),
      low_edu     = yesno(education == "Less than high school", !is.na(education))
    )
}
dat <- map(c(Male = "Male", Female = "Female"), prep)
n_total <- nrow(dat$Male) + nrow(dat$Female)
message("Analytic sample: ", n_total, " (", nrow(dat$Male), " men, ", nrow(dat$Female), " women)")
stopifnot(is.na(EXPECTED_N) || n_total == EXPECTED_N || QUICK)

# Scored items (adverse = "Yes"); low PIR is derived from imputed PIR after imputation
SCORE_ITEMS <- c("low_pir", "food_insecure", "snap", "low_edu", "uninsured", "poor_access", "renter", "alone")
CONTEXT     <- c("foreign_born", "spanish_interview")
NEW_IMPUTE  <- c("snap", "uninsured", "poor_access", "renter", "alone", "unemployed", CONTEXT, "food_insecure")
# The SDoH block used for prediction (individual items, education and PIR in full detail)
BLOCK <- c("education", "pir", "food_insecure", "snap", "uninsured", "poor_access", "renter", "alone", CONTEXT)
# Nativity and interview language partly stand in for ethnicity, so Aim 3 is also run without them
BLOCK_NOCX <- setdiff(BLOCK, CONTEXT)
# Items that cannot be built in every cycle (renter: no housing file for 2017-March 2020)
# are left out of the score, the imputation and the block, in all cycles
SCORE_ITEMS <- drop_unavailable(SCORE_ITEMS); NEW_IMPUTE <- drop_unavailable(NEW_IMPUTE)
BLOCK <- drop_unavailable(BLOCK); BLOCK_NOCX <- drop_unavailable(BLOCK_NOCX)
message("Scored items (", length(SCORE_ITEMS), "): ", paste(SCORE_ITEMS, collapse = ", "))

# ---- Harmonization check: % "Yes" and % missing by cycle (look at this first) ----
chk <- bind_rows(dat) %>%
  mutate(low_pir = yesno(pir < 1.3, pir >= 1.3)) %>%
  select(cycle, all_of(c(SCORE_ITEMS, "unemployed", CONTEXT))) %>%
  mutate(across(-cycle, as.character)) %>%
  pivot_longer(-cycle, names_to = "item") %>%
  group_by(item, cycle) %>%
  summarise(n = n(), pct_missing = round(100 * mean(is.na(value)), 1),
            pct_yes = round(100 * mean(value == "Yes", na.rm = TRUE), 1), .groups = "drop")
write_csv(chk, out("00_harmonization_check.csv"))
print(chk, n = Inf)
bad <- chk %>% filter(pct_missing > 20)
if (nrow(bad)) { print(bad, n = Inf)
  if (Sys.getenv("SDOH_IGNORE_CHECK") != "TRUE")
    stop("Some items have >20% missing in some cycles - inspect 00_harmonization_check.csv before continuing") }
message("Harmonization check passed (no item has >20% missing in any cycle).")

# =============================================================================
# 3. Multiple imputation (new SDoH items added to the original predictors)
# =============================================================================
impute <- function(sx, with_outcome, m, seed) {
  d <- dat[[sx]]
  vars <- c(BASE_COMMON, if (sx == "Female") WOMEN_EXTRA, "education", "pir", NEW_IMPUTE)
  X <- as.data.frame(d[, vars])
  if (with_outcome) X$tscore <- d$tscore
  mice(X, m = m, method = make.method(X), seed = seed, printFlag = FALSE, maxit = 10)
}
completed <- function(imp, sx) lapply(seq_len(imp$m), function(i) {
  X <- complete(imp, i); X$tscore <- NULL
  x <- bind_cols(dat[[sx]][, c("SEQN", "tscore", "osteoporosis", "wt_pool", "SDMVPSU", "SDMVSTRA", "sex")], X)
  x$low_pir <- factor(ifelse(x$pir < 1.3, "Yes", "No"), levels = c("No", "Yes"))
  x$low_edu <- factor(ifelse(x$education == "Less than high school", "Yes", "No"), levels = c("No", "Yes"))
  x$sdoh_score <- rowSums(sapply(SCORE_ITEMS, function(v) x[[v]] == "Yes"))
  x$score_cat  <- factor(pmin(x$sdoh_score, 3), levels = 0:3, labels = c("0", "1", "2", "3+"))
  x$score_trend <- pmin(x$sdoh_score, 3)
  x$sdoh_score_emp <- x$sdoh_score + (x$unemployed == "Yes")
  x
})
# Saved imputations and CV folds are tied to a fingerprint of the input data, so if the
# variable coding changes they are rebuilt instead of silently reused
fp <- function(x) { f <- tempfile(); saveRDS(x, f, compress = FALSE, version = 3); unname(tools::md5sum(f)) }
DAT_FP <- substr(fp(lapply(dat, function(d) as.data.frame(d[, c("SEQN", NEW_IMPUTE, "education", "pir")]))), 1, 8)
message("Data fingerprint: ", DAT_FP)
get_imp <- function(tag, with_outcome, m, seed) {
  f <- out(sprintf("imputation_%s_%s.rds", tag, DAT_FP))
  if (file.exists(f)) return(readRDS(f))
  message("Imputing (", tag, ") - several minutes per sex, no messages while it runs")
  r <- map(c(Male = "Male", Female = "Female"), ~ impute(.x, with_outcome, m, seed))
  saveRDS(r, f); r
}
imp_pred  <- get_imp("pred_no_outcome",  FALSE, CFG$m_pred,  CFG$seed)
imp_assoc <- get_imp("assoc_with_outcome", TRUE, CFG$m_assoc, CFG$seed + 7)
AL <- map(c(Male = "Male", Female = "Female"), ~ completed(imp_assoc[[.x]], .x))

# =============================================================================
# 4. Descriptives: prevalence of each item and score distribution (imputation 1
#    of the association imputation, survey-weighted), and item correlations
# =============================================================================
desc <- map_dfr(c("Male", "Female"), function(sx) {
  x <- AL[[sx]][[1]]
  des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pool, nest = TRUE, data = x)
  items <- map_dfr(c(SCORE_ITEMS, "unemployed", CONTEXT), function(v)
    tibble(sex = sx, measure = v, level = "Yes",
           pct = 100 * coef(svymean(as.formula(paste0("~I(", v, ' == "Yes")')), des))[[2]]))
  sc <- svymean(~score_cat, des)
  bind_rows(items, tibble(sex = sx, measure = "score_cat", level = levels(x$score_cat), pct = 100 * coef(sc)),
            tibble(sex = sx, measure = "score_mean", level = "mean", pct = coef(svymean(~sdoh_score, des))[[1]]))
})
write_csv(desc, out("01_item_prevalence_and_score.csv"))
cors <- map_dfr(c("Male", "Female"), function(sx) {
  x <- AL[[sx]][[1]]
  M <- sapply(SCORE_ITEMS, function(v) as.numeric(x[[v]] == "Yes"))
  C <- cov.wt(M, wt = x$wt_pool, cor = TRUE)$cor
  as_tibble(C, rownames = "item") %>% mutate(sex = sx, .before = 1)
})
write_csv(cors, out("02_item_correlations_weighted.csv"))

# =============================================================================
# 5. Aim 1 - associations (outcome-inclusive imputation, Rubin's rules)
# =============================================================================
wald_p <- function(b, V, tt) { w <- as.numeric(t(b[tt]) %*% solve(V[tt, tt]) %*% b[tt])
  c(chisq = w, df = length(tt), p = pchisq(w, length(tt), lower.tail = FALSE)) }
adj <- list("M0 age, cycle"                    = c("age", "cycle"),
            "M1 + race and ethnicity"           = c("age", "cycle", "race"),
            "M2 + weight, height"               = c("age", "cycle", "race", "weight", "height"),
            "M3 all base predictors"            = "ALL")
fit_pool <- function(des_imp, f, family = gaussian()) {
  p <- MIcombine(with(des_imp, svyglm(f, family = family))); list(b = coef(p), V = vcov(p))
}
tidy_terms <- function(fp, pattern, ...) {
  tt <- grep(pattern, names(fp$b), value = TRUE); se <- sqrt(diag(fp$V))[tt]
  tibble(term = tt, est = fp$b[tt], se = se, lo = fp$b[tt] - 1.96 * se, hi = fp$b[tt] + 1.96 * se,
         p = 2 * pnorm(-abs(fp$b[tt] / se)), ...)
}
aim1 <- list()
for (sx in c("Male", "Female")) {
  di <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pool, nest = TRUE,
                  data = imputationList(AL[[sx]]))
  for (mn in names(adj)) {
    cov <- adj[[mn]]; if (identical(cov, "ALL")) cov <- c(BASE_COMMON, if (sx == "Female") WOMEN_EXTRA)
    rhs <- paste(cov, collapse = " + ")
    # T-score: categories (joint Wald) and linear trend across 0/1/2/3+
    fc <- fit_pool(di, as.formula(paste("tscore ~ score_cat +", rhs)))
    ft <- fit_pool(di, as.formula(paste("tscore ~ score_trend +", rhs)))
    # Osteoporosis: odds ratios (quasibinomial for survey data)
    fo <- fit_pool(di, as.formula(paste("osteoporosis ~ score_cat +", rhs)), quasibinomial())
    fot <- fit_pool(di, as.formula(paste("osteoporosis ~ score_trend +", rhs)), quasibinomial())
    aim1[[length(aim1) + 1]] <- bind_rows(
      tidy_terms(fc, "^score_cat", outcome = "T-score (beta)",
                 p_joint = wald_p(fc$b, fc$V, grep("^score_cat", names(fc$b), value = TRUE))[["p"]]),
      tidy_terms(ft, "^score_trend$", outcome = "T-score (beta per category)"),
      tidy_terms(fo, "^score_cat", outcome = "Osteoporosis (log OR)",
                 p_joint = wald_p(fo$b, fo$V, grep("^score_cat", names(fo$b), value = TRUE))[["p"]]),
      tidy_terms(fot, "^score_trend$", outcome = "Osteoporosis (log OR per category)")) %>%
      mutate(sex = sx, model = mn, .before = 1)
  }
  message("Aim 1 done: ", sx)
}
aim1 <- bind_rows(aim1) %>%
  mutate(OR = ifelse(grepl("log OR", outcome), exp(est), NA), OR_lo = ifelse(grepl("log OR", outcome), exp(lo), NA),
         OR_hi = ifelse(grepl("log OR", outcome), exp(hi), NA))
write_csv(aim1, out("10_aim1_score_associations.csv"))

# Interactions (both sexes combined; female reproductive variables omitted)
ml <- AL$Male; fl <- AL$Female
comb <- lapply(seq_along(ml), function(i) { x <- bind_rows(ml[[i]], fl[[i]] %>% select(-any_of(WOMEN_EXTRA)))
  x$sex <- factor(x$sex, levels = c("Male", "Female")); x })
dc <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pool, nest = TRUE, data = imputationList(comb))
rhs_all <- paste(BASE_COMMON, collapse = " + ")
fsx <- fit_pool(dc, as.formula(paste("tscore ~ sex * score_trend +", rhs_all)))
frc <- fit_pool(dc, as.formula(paste("tscore ~ sex + race * score_trend +", setdiff(BASE_COMMON, "race") %>% paste(collapse = " + "))))
inter <- bind_rows(
  tidy_terms(fsx, "score_trend", test = "sex x score"),
  tidy_terms(frc, "score_trend", test = "race x score"))
wsx <- wald_p(fsx$b, fsx$V, grep(":score_trend|score_trend:", names(fsx$b), value = TRUE))
wrc <- wald_p(frc$b, frc$V, grep(":score_trend|score_trend:", names(frc$b), value = TRUE))
write_csv(inter, out("11_aim1_interaction_terms.csv"))
write_csv(tibble(test = c("sex x score (trend)", "race x score (trend)"),
                 chisq = c(wsx[["chisq"]], wrc[["chisq"]]), df = c(wsx[["df"]], wrc[["df"]]),
                 p = c(wsx[["p"]], wrc[["p"]])), out("12_aim1_interaction_wald.csv"))

# Aim 3b - race coefficients with and without the SDoH block and body size
race_tab <- list()
for (sx in c("Male", "Female")) {
  race_sets <- list(
    "R0 race + age, cycle"                       = c("age", "cycle"),
    "R1 + SDoH block"                            = c("age", "cycle", BLOCK),
    "R1b + SDoH block without nativity/language" = c("age", "cycle", BLOCK_NOCX),
    "R2 + weight, height"                        = c("age", "cycle", "weight", "height"),
    "R3 + SDoH block + weight, height"           = c("age", "cycle", BLOCK, "weight", "height"),
    "R4 + SDoH block + all base predictors"      = c(setdiff(BASE_COMMON, "race"),
                                                     if (sx == "Female") WOMEN_EXTRA, BLOCK))
  di <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pool, nest = TRUE,
                  data = imputationList(AL[[sx]]))
  for (rn in names(race_sets)) {
    fr <- fit_pool(di, as.formula(paste("tscore ~ race +", paste(unique(race_sets[[rn]]), collapse = " + "))))
    race_tab[[length(race_tab) + 1]] <- tidy_terms(fr, "^race", sex = sx, model = rn)
  }
}
write_csv(bind_rows(race_tab), out("13_aim3_race_coefficients.csv"))
message("Association analyses done.")

# =============================================================================
# 6. Aims 2 and 3 - cross-validated prediction (outcome-excluding imputation)
# =============================================================================
wmetrics <- function(y, yhat, w) { w <- w / sum(w); e <- y - yhat; ybar <- sum(w * y); mse <- sum(w * e^2)
  c(MSE = mse, RMSE = sqrt(mse), R2 = 1 - mse / sum(w * (y - ybar)^2)) }
wauc <- function(risk, case, w) {
  ok <- !is.na(risk) & !is.na(case); risk <- risk[ok]; case <- case[ok]; w <- w[ok]
  if (length(unique(case)) < 2) return(NA_real_)
  u <- sort(unique(risk)); g <- factor(match(risk, u), levels = seq_along(u))
  wn <- tapply(w * (case == 0), g, sum); wn[is.na(wn)] <- 0
  wp <- tapply(w * (case == 1), g, sum); wp[is.na(wp)] <- 0
  sum(wp * (cumsum(wn) - wn + 0.5 * wn)) / (sum(wp) * sum(wn)) }
nb_ci <- function(x, ntr, nts) { x <- x[!is.na(x)]; J <- length(x); m <- mean(x)
  if (J < 2) return(c(est = m, lo = NA, hi = NA, p = NA))
  se <- sqrt((1 / J + nts / ntr) * var(x)); tq <- qt(0.975, J - 1)
  c(est = m, lo = m - tq * se, hi = m + tq * se, p = 2 * pt(-abs(m / se), J - 1)) }
tune_rf <- function(tr, preds, w, seed) {
  f <- reformulate(preds, "tscore")
  errs <- vapply(seq_len(nrow(CFG$grid)), function(i) ranger(f, data = tr, num.trees = CFG$tune_trees,
    mtry = min(CFG$grid$mtry[i], length(preds)), min.node.size = CFG$grid$min.node.size[i],
    case.weights = w, respect.unordered.factors = "order", num.threads = CFG$threads,
    seed = seed)$prediction.error, numeric(1))
  as.list(CFG$grid[which.min(errs), ]) }
fit_rf <- function(tr, preds, w, seed, p) ranger(reformulate(preds, "tscore"), data = tr,
  num.trees = CFG$num_trees, mtry = min(p$mtry, length(preds)), min.node.size = p$min.node.size,
  case.weights = w, respect.unordered.factors = "order", num.threads = CFG$threads, seed = seed)
lm_formula <- function(preds) as.formula(paste("tscore ~", paste(vapply(preds, function(p)
  if (p %in% c("age", "weight", "height", "pir")) sprintf("ns(%s, df = 3)", p) else p, ""), collapse = " + ")))

model_sets <- function(sx) {
  base <- c(BASE_COMMON, if (sx == "Female") WOMEN_EXTRA)
  nor  <- setdiff(base, "race")
  list(base         = base,
       base_block   = c(base, BLOCK),
       base_score   = c(base, "sdoh_score"),
       base_blk_emp = c(base, BLOCK, "unemployed"),
       norace       = nor,
       norace_block = c(nor, BLOCK),
       norace_block_nocx = c(nor, BLOCK_NOCX))
}
COMPARISONS <- list(
  "Aim 2: SDoH block added"                      = c("base_block", "base"),
  "Aim 2: SDoH score added"                      = c("base_score", "base"),
  "Sensitivity: block + employment added"        = c("base_blk_emp", "base"),
  "Aim 3: loss from removing race"               = c("norace", "base"),
  "Aim 3: gain from SDoH block without race"     = c("norace_block", "norace"),
  "Aim 3: gain from SDoH without race, no nativity/language" = c("norace_block_nocx", "norace"),
  "Aim 3: no-race + SDoH vs model with race"     = c("norace_block", "base"))

perf <- list()
for (sx in c("Male", "Female")) {
  dl <- completed(imp_pred[[sx]], sx); ms <- model_sets(sx)
  set.seed(CFG$seed + sum(utf8ToInt(paste("sdoh_bundle", sx))))
  n <- nrow(dl[[1]]); folds <- list()
  for (r in seq_len(CFG$n_repeats)) { id <- sample(rep(seq_len(CFG$k_folds), length.out = n))
    for (k in seq_len(CFG$k_folds)) folds[[length(folds) + 1]] <- list(rep = r, fold = k, test = which(id == k)) }
  ck <- out(file.path("checkpoints", paste0(sx, "_", DAT_FP))); dir.create(ck, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time(); nnew <- 0
  for (fi in seq_along(folds)) {
    f_ck <- file.path(ck, sprintf("fold_%03d.rds", fi))
    if (file.exists(f_ck)) { perf[[length(perf) + 1]] <- readRDS(f_ck); next }
    te <- folds[[fi]]$test; tr1 <- dl[[1]][-te, ]
    params <- tune_rf(tr1, ms$base_block, tr1$wt_pool / mean(tr1$wt_pool), CFG$seed + fi)
    res <- list()
    for (mi in seq_along(dl)) {
      tr <- dl[[mi]][-te, ]; ts <- dl[[mi]][te, ]; w_tr <- tr$wt_pool / mean(tr$wt_pool); tr$.w <- w_tr
      for (mn in names(ms)) {
        p  <- predict(fit_rf(tr, ms[[mn]], w_tr, CFG$seed + fi, params), ts, num.threads = CFG$threads)$predictions
        pl <- predict(lm(lm_formula(ms[[mn]]), data = tr, weights = .w), ts)
        for (L in list(list("Random forest", p), list("Linear (splines)", pl)))
          res[[length(res) + 1]] <- tibble(sex = sx, rep = folds[[fi]]$rep, fold = folds[[fi]]$fold, imp = mi,
            learner = L[[1]], model = mn, !!!as.list(wmetrics(ts$tscore, L[[2]], ts$wt_pool)),
            AUC_osteoporosis = wauc(-L[[2]], ts$osteoporosis, ts$wt_pool), n_train = nrow(tr), n_test = nrow(ts))
      }
    }
    r <- bind_rows(res); saveRDS(r, f_ck); perf[[length(perf) + 1]] <- r; nnew <- nnew + 1
    mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("  SDoH bundle CV | %s: fold %d / %d (%.1f min, about %.0f min left)", sx, fi,
                    length(folds), mins, (length(folds) - fi) * mins / nnew))
  }
}
perf <- bind_rows(perf); write_csv(perf, out("20_cv_raw_performance.csv"))
pf <- perf %>% group_by(sex, learner, model, rep, fold) %>%
  summarise(across(c(R2, RMSE, AUC_osteoporosis, n_train, n_test), mean), .groups = "drop")

overall <- pf %>% pivot_longer(c(R2, RMSE, AUC_osteoporosis), names_to = "metric") %>%
  group_by(sex, learner, model, metric) %>%
  summarise(ci = list(nb_ci(value, mean(n_train), mean(n_test))), .groups = "drop") %>%
  mutate(est = map_dbl(ci, ~ .x[["est"]]), lo = map_dbl(ci, ~ .x[["lo"]]), hi = map_dbl(ci, ~ .x[["hi"]])) %>%
  select(-ci)
write_csv(overall, out("21_performance_by_model.csv"))

long <- pf %>% pivot_longer(c(R2, RMSE, AUC_osteoporosis), names_to = "metric") %>%
  select(sex, learner, rep, fold, metric, model, value, n_train, n_test) %>%
  pivot_wider(names_from = model, values_from = value)
deltas <- map_dfr(names(COMPARISONS), function(cn) {
  a <- COMPARISONS[[cn]][1]; b <- COMPARISONS[[cn]][2]
  long %>% mutate(diff = .data[[a]] - .data[[b]]) %>% group_by(sex, learner, metric) %>%
    summarise(ci = list(nb_ci(diff, mean(n_train), mean(n_test))), .groups = "drop") %>%
    mutate(comparison = cn, est = map_dbl(ci, ~ .x[["est"]]), lo = map_dbl(ci, ~ .x[["lo"]]),
           hi = map_dbl(ci, ~ .x[["hi"]]), p = map_dbl(ci, ~ .x[["p"]])) %>% select(-ci)
})
write_csv(deltas, out("22_added_value_deltas.csv"))

# Aim 3: share of the accuracy lost by removing race that the SDoH block recovers
# recovery = mean(norace_block - norace) / mean(base - norace); CI by bootstrapping
# repeats (whole CV repeats resampled, keeping folds together). Interpret with care
# when the loss is small (the ratio is then unstable).
set.seed(CFG$seed)
# NOTE: the interval only reflects how the CV splits fell (a Monte Carlo interval), not
# population uncertainty. Interpret the ratio only where the "loss from removing race"
# CI in 22_added_value_deltas.csv excludes zero.
recovery <- long %>% group_by(sex, learner, metric) %>% group_modify(function(g, key) {
  map_dfr(c(with_nativity_language = "norace_block", without_nativity_language = "norace_block_nocx"),
          function(m) {
    est <- mean(g[[m]] - g$norace) / mean(g$base - g$norace)
    reps <- unique(g$rep)
    boot <- if (length(reps) > 1) replicate(2000, {
      s <- sample(reps, replace = TRUE); gg <- bind_rows(lapply(s, function(r) g[g$rep == r, ]))
      mean(gg[[m]] - gg$norace) / mean(gg$base - gg$norace) }) else NA
    tibble(sdoh_block = m, recovery = est, mc_lo = quantile(boot, 0.025, na.rm = TRUE),
           mc_hi = quantile(boot, 0.975, na.rm = TRUE), loss_from_removing_race = mean(g$base - g$norace))
  })
}) %>% ungroup()
write_csv(recovery, out("23_aim3_recovery_of_race_information.csv"))

message("All SDoH bundle analyses done. Send back the CSV files in ", OUT_DIR,
        " (00_ to 23_; not the checkpoints or imputation .rds files).")
