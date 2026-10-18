# =============================================================================
#  SDoH and race: cross-validated prediction, version 2 (revision after audit)
#
#  What changed from sdoh_bundle.R (v1), and why
#   1. Four prediction strategies, each with its own information set:
#        A  clinical                       (clinical + background, no race, no social)
#        B  clinical + social block        (no race)
#        C  clinical + race and ethnicity  (no social)
#        D  clinical + race + social block
#      Key paired contrasts: D-C (social gain with race), B-A (social gain
#      without race), A-C (omitting race), B-C (social in place of race).
#      Secondary: B without nativity/language, C + disadvantage score, and the
#      Osteoporosis Self-assessment Tool (OST = 0.2 x [weight kg - age]) as a
#      simple benchmark.
#   2. Imputation is now done INSIDE each cross-validation fold, fitted on the
#      training participants only (mice "ignore" argument), and separately for
#      each information set. A strategy without race never uses race to fill
#      in its predictors; a clinical-only strategy never uses social variables.
#   3. Random forest tuning is done for each strategy separately (training data).
#   4. Performance is computed once per CV repetition from all held-out
#      predictions (each participant predicted once per repetition), not per fold.
#      This matters for AUC in men (about 134 cases in total).
#   5. Uncertainty: survey-design bootstrap (Rao-Wu subbootstrap of PSUs within
#      strata, 500 replicates) applied to the held-out predictions. The interval
#      describes sampling uncertainty for the population performance of the
#      fitted cross-validated predictions; it does not re-run imputation and
#      model fitting in each replicate (label it this way in the paper).
#      The spread across the 10 CV repetitions is reported as well.
#   6. No equivalence margin is used.
#
#  Association analyses (Aim 1, race coefficients) are NOT repeated: they use
#  outcome-inclusive imputation, which is appropriate and unchanged (v1 files 10-13).
#
#  Inputs (same as v1): <results folder>/imputations.rds, nhanes_cache/ and nhanes_cycles.R
#  <results folder> = results/ (NHANES_LATEST=J) or results_2017_2020/ (default, P)
#  Run (Slurm): sbatch run_sdoh_v2.slurm      QUICK test first:  QUICK=TRUE in the script
#  Outputs: <results folder>/sdoh_v2/ (QUICK: .../sdoh_v2_quick/). Folds are checkpointed,
#  so re-submitting resumes.
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
CFG <- list(seed = 20261010, k_folds = 10, n_repeats = 10, m = 3, maxit = 10, B = 500,
            num_trees = 500, tune_trees = 300,
            grid = expand.grid(mtry = c(2, 4, 6), min.node.size = c(10, 50)), threads = detect_threads())
if (QUICK) CFG <- modifyList(CFG, list(k_folds = 3, n_repeats = 2, m = 1, maxit = 5, B = 50,
                                       num_trees = 100, tune_trees = 100,
                                       grid = expand.grid(mtry = c(2, 4), min.node.size = 20)))
OUT_DIR <- file.path(RESULTS_ROOT, paste0(if (QUICK) "sdoh_v2_quick" else "sdoh_v2", if (COMMON) "_common" else ""))
if (COMMON) message("SDOH_COMMON: common item definitions in all survey periods -> ", OUT_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)
message(sprintf("SDoH v2 prediction | %s | cores: %d | %d x %d CV | m = %d per fold | %d bootstrap replicates",
                if (QUICK) "QUICK TEST" else "FULL RUN", CFG$threads, CFG$k_folds, CFG$n_repeats, CFG$m, CFG$B))

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

rm(imps); invisible(gc())
fp <- function(x) { f <- tempfile(); saveRDS(x, f, compress = FALSE, version = 3); unname(tools::md5sum(f)) }
DAT_FP <- substr(fp(lapply(dat, function(d) as.data.frame(d[, c("SEQN", NEW_IMPUTE, "education", "pir")]))), 1, 8)
message("Data fingerprint: ", DAT_FP)

# =============================================================================
# 3. Strategies and the information each one is allowed to use
# =============================================================================
CLIN <- function(sx) c(setdiff(BASE_COMMON, "race"), if (sx == "Female") WOMEN_EXTRA)
# Information sets = the variables an imputation model may see (outcome never included)
info_sets <- function(sx) { a <- CLIN(sx)
  list(A   = a,
       B   = c(a, BLOCK),
       Bnc = c(a, BLOCK_NOCX),
       C   = c(a, "race"),
       D   = c(a, "race", BLOCK),
       Dnc = c(a, "race", BLOCK_NOCX)) }
strategies <- function(sx) { a <- CLIN(sx)
  list(A           = list(set = "A",   preds = a),
       B           = list(set = "B",   preds = c(a, BLOCK)),
       C           = list(set = "C",   preds = c(a, "race")),
       D           = list(set = "D",   preds = c(a, "race", BLOCK)),
       B_nocontext = list(set = "Bnc", preds = c(a, BLOCK_NOCX)),
       C_score     = list(set = "Dnc", preds = c(a, "race", "sdoh_score")),
       OST         = list(set = "A",   preds = "ost", linear_only = TRUE)) }
CONTRASTS <- list(
  "D - C: social block added to model with race"       = c("D", "C"),
  "B - A: social block added to model without race"    = c("B", "A"),
  "A - C: race omitted"                                = c("A", "C"),
  "B - C: social block in place of race"               = c("B", "C"),
  "B_nocontext - A: social block without nativity/language, no race" = c("B_nocontext", "A"),
  "B_nocontext - C: social block without nativity/language in place of race" = c("B_nocontext", "C"),
  "C_score - C: disadvantage score added to model with race" = c("C_score", "C"),
  "A - OST: clinical model vs OST"                     = c("A", "OST"),
  "C - OST: clinical model with race vs OST"           = c("C", "OST"))
META <- c("SEQN", "tscore", "osteoporosis", "wt_pool")

# Derived variables are always rebuilt from the completed components
derive <- function(x) {
  yn <- function(z) factor(ifelse(z, "Yes", "No"), levels = c("No", "Yes"))
  if ("pir" %in% names(x))       x$low_pir <- yn(x$pir < 1.3)
  if ("education" %in% names(x)) x$low_edu <- yn(x$education == "Less than high school")
  if (all(SCORE_ITEMS %in% names(x))) x$sdoh_score <- rowSums(sapply(SCORE_ITEMS, function(v) x[[v]] == "Yes"))
  if (all(c("age", "weight") %in% names(x))) x$ost <- 0.2 * (x$weight - x$age)   # weight in kg
  x
}
# Imputation fitted on training rows only; test rows are filled in from the training models
run_mice <- function(d, vars, is_test, m, seed) {
  X <- as.data.frame(d[, vars])
  if (!anyNA(X)) return(rep(list(X), m))
  imp <- mice(X, m = m, method = make.method(X), ignore = is_test, seed = seed,
              printFlag = FALSE, maxit = CFG$maxit)
  L <- lapply(seq_len(m), function(i) complete(imp, i))
  left <- names(which(colSums(is.na(L[[1]])) > 0))
  if (length(left)) stop("Imputation left missing values in: ", paste(left, collapse = ", "))
  L
}

tune_rf <- function(tr, preds, w, seed) {
  f <- reformulate(preds, "tscore")
  errs <- vapply(seq_len(nrow(CFG$grid)), function(i) ranger(f, data = tr, num.trees = CFG$tune_trees,
    mtry = min(CFG$grid$mtry[i], length(preds)), min.node.size = CFG$grid$min.node.size[i],
    case.weights = w, respect.unordered.factors = "order", num.threads = 1,
    seed = seed)$prediction.error, numeric(1))
  as.list(CFG$grid[which.min(errs), ]) }
fit_rf <- function(tr, preds, w, seed, p) ranger(reformulate(preds, "tscore"), data = tr,
  num.trees = CFG$num_trees, mtry = min(p$mtry, length(preds)), min.node.size = p$min.node.size,
  case.weights = w, respect.unordered.factors = "order", num.threads = 1, seed = seed)
lm_formula <- function(preds) as.formula(paste("tscore ~", paste(vapply(preds, function(p)
  if (p %in% c("age", "weight", "height", "pir")) sprintf("ns(%s, df = 3)", p) else p, ""), collapse = " + ")))

# =============================================================================
# 4. Cross-validation: one checkpoint file per fold, folds run in parallel
# =============================================================================
make_folds <- function(sx) {
  set.seed(CFG$seed + sum(utf8ToInt(paste("sdoh_v2", sx))))
  n <- nrow(dat[[sx]]); fl <- list()
  for (r in seq_len(CFG$n_repeats)) { id <- sample(rep(seq_len(CFG$k_folds), length.out = n))
    for (k in seq_len(CFG$k_folds)) fl[[length(fl) + 1]] <- list(rep = r, fold = k, test = which(id == k)) }
  fl }
FOLDS <- map(c(Male = "Male", Female = "Female"), make_folds)
ck_dir <- function(sx) out(file.path("checkpoints", paste0(sx, "_", DAT_FP)))
walk(names(FOLDS), ~ dir.create(ck_dir(.x), recursive = TRUE, showWarnings = FALSE))
ck_file <- function(sx, fi) file.path(ck_dir(sx), sprintf("fold_%03d.rds", fi))

run_fold <- function(sx, fi) {
  t0 <- Sys.time(); fold <- FOLDS[[sx]][[fi]]
  d <- dat[[sx]]; is_te <- seq_len(nrow(d)) %in% fold$test
  sets <- info_sets(sx); st <- strategies(sx)
  comp <- list()
  for (s in unique(map_chr(st, "set")))
    comp[[s]] <- lapply(run_mice(d, sets[[s]], is_te, CFG$m, CFG$seed + 1000 * fi + match(s, names(sets))),
                        function(X) derive(bind_cols(d[, META], X)))
  res <- list()
  for (sn in names(st)) {
    s <- st[[sn]]; L <- comp[[s$set]]; rf <- !isTRUE(s$linear_only)
    tr1 <- L[[1]][!is_te, ]
    if (rf) params <- tune_rf(tr1, s$preds, tr1$wt_pool / mean(tr1$wt_pool), CFG$seed + fi)
    p_lm <- p_rf <- matrix(NA_real_, sum(is_te), length(L))
    for (i in seq_along(L)) {
      tr <- L[[i]][!is_te, ]; ts <- L[[i]][is_te, ]; tr$.w <- tr$wt_pool / mean(tr$wt_pool)
      p_lm[, i] <- predict(lm(lm_formula(s$preds), data = tr, weights = .w), ts)
      if (rf) p_rf[, i] <- predict(fit_rf(tr, s$preds, tr$.w, CFG$seed + fi, params), ts, num.threads = 1)$predictions
    }
    # predictions averaged over the m imputations of this fold
    res[[length(res) + 1]] <- tibble(sex = sx, rep = fold$rep, fold = fold$fold, SEQN = d$SEQN[is_te],
                                     strategy = sn, learner = "Linear (splines)", pred = rowMeans(p_lm))
    if (rf) res[[length(res) + 1]] <- tibble(sex = sx, rep = fold$rep, fold = fold$fold, SEQN = d$SEQN[is_te],
                                             strategy = sn, learner = "Random forest", pred = rowMeans(p_rf))
  }
  r <- bind_rows(res); saveRDS(r, ck_file(sx, fi))
  message(sprintf("  v2 CV | %s fold %d / %d done in %.1f min (%s)", sx, fi, length(FOLDS[[sx]]),
                  as.numeric(difftime(Sys.time(), t0, units = "mins")), format(Sys.time(), "%H:%M")))
  invisible(TRUE)
}

jobs <- bind_rows(map(names(FOLDS), ~ tibble(sex = .x, fi = seq_along(FOLDS[[.x]])))) %>%
  filter(!file.exists(ck_file(sex, fi)))
message(sprintf("%d folds to run (%d already saved) on %d cores", nrow(jobs),
                sum(lengths(FOLDS)) - nrow(jobs), CFG$threads))
if (nrow(jobs)) {
  rr <- parallel::mclapply(seq_len(nrow(jobs)), function(j)
          tryCatch(run_fold(jobs$sex[j], jobs$fi[j]), error = function(e) conditionMessage(e)),
          mc.cores = CFG$threads, mc.preschedule = FALSE)
  bad <- which(!vapply(rr, isTRUE, logical(1)))
  if (length(bad)) { for (b in bad) message("FAILED ", jobs$sex[b], " fold ", jobs$fi[b], ": ", rr[[b]])
    stop(length(bad), " folds failed - fix the error above, then resubmit (finished folds are kept)") }
}
P <- bind_rows(map(names(FOLDS), function(sx) map_dfr(seq_along(FOLDS[[sx]]), ~ readRDS(ck_file(sx, .x)))))
saveRDS(P, out("34_v2_heldout_predictions.rds"))

# =============================================================================
# 5. Performance per repetition, with survey-design bootstrap replicate weights
# =============================================================================
wm_mat <- function(y, yhat, W) {
  Wn <- sweep(W, 2, colSums(W), "/"); mse <- colSums(Wn * (y - yhat)^2)
  ybar <- colSums(Wn * y); sst <- colSums(Wn * sweep(matrix(y, length(y), ncol(W)), 2, ybar)^2)
  list(R2 = 1 - mse / sst, RMSE = sqrt(mse)) }
wauc_mat <- function(risk, case, W) {   # same estimator as v1, for every weight column at once
  g <- match(risk, sort(unique(risk)))
  wn <- rowsum(W * (case == 0), g, reorder = TRUE); wp <- rowsum(W * (case == 1), g, reorder = TRUE)
  cn <- apply(wn, 2, cumsum); if (!is.matrix(cn)) cn <- matrix(cn, nrow = nrow(wn))
  colSums(wp * (cn - 0.5 * wn)) / (colSums(wp) * colSums(wn)) }

counts <- list(); per_rep <- list(); M <- list()   # M[[sex]][[metric]][[strategy|learner]] = reps x (1 + B)
for (sx in names(FOLDS)) {
  M[[sx]] <- list(R2 = list(), RMSE = list(), AUC_osteoporosis = list())
  dd <- dat[[sx]] %>% mutate(stratum = interaction(cycle, SDMVSTRA, drop = TRUE))
  des <- svydesign(ids = ~SDMVPSU, strata = ~stratum, weights = ~wt_pool, nest = TRUE, data = dd)
  set.seed(CFG$seed + 99)
  rdes <- as.svrepdesign(des, type = "subbootstrap", replicates = CFG$B)
  W <- cbind(dd$wt_pool, as.matrix(weights(rdes, type = "analysis")))
  ratio <- colSums(W[, -1]) / sum(dd$wt_pool)
  if (any(ratio < 0.3 | ratio > 3)) stop("Replicate weights look wrong (total weight ratio out of range)")
  counts[[sx]] <- tibble(sex = sx, n = nrow(dd), n_osteoporosis = sum(dd$osteoporosis == 1),
                         weighted_prevalence_pct = 100 * weighted.mean(dd$osteoporosis == 1, dd$wt_pool),
                         n_strata = nlevels(dd$stratum), n_psu = nrow(unique(dd[, c("stratum", "SDMVPSU")])))
  Ps <- P %>% filter(sex == sx)
  for (key in unique(paste(Ps$strategy, Ps$learner, sep = "|"))) {
    sl <- strsplit(key, "|", fixed = TRUE)[[1]]
    for (r in seq_len(CFG$n_repeats)) {
      pr <- Ps %>% filter(strategy == sl[1], learner == sl[2], rep == r)
      yhat <- pr$pred[match(dd$SEQN, pr$SEQN)]
      if (anyNA(yhat) || nrow(pr) != nrow(dd)) stop("Repetition ", r, " of ", key, " does not predict every participant once")
      mm <- wm_mat(dd$tscore, yhat, W)
      mm$AUC_osteoporosis <- wauc_mat(-yhat, dd$osteoporosis, W)
      for (mt in names(mm)) M[[sx]][[mt]][[key]] <- rbind(M[[sx]][[mt]][[key]], mm[[mt]])
      per_rep[[length(per_rep) + 1]] <- tibble(sex = sx, strategy = sl[1], learner = sl[2], rep = r,
        R2 = mm$R2[1], RMSE = mm$RMSE[1], AUC_osteoporosis = mm$AUC_osteoporosis[1])
    }
  }
  message("Evaluation done: ", sx)
}
write_csv(bind_rows(counts), out("30_v2_event_counts.csv"))
write_csv(bind_rows(per_rep), out("35_v2_per_repetition.csv"))

# A matrix (reps x [full, B replicates]) -> estimate, SE and CI
summ <- function(X) {
  est <- mean(X[, 1]); boot <- colMeans(X[, -1, drop = FALSE])
  cv_sd <- if (nrow(X) > 1) sd(X[, 1]) else NA_real_
  se <- sqrt(var(boot) + ifelse(is.na(cv_sd), 0, cv_sd^2 / nrow(X)))
  tibble(est = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se,
         p = 2 * pnorm(-abs(est / se)), pct_lo = quantile(boot, 0.025), pct_hi = quantile(boot, 0.975),
         cv_sd_across_repetitions = cv_sd) }

perf <- map_dfr(names(M), function(sx) map_dfr(names(M[[sx]]), function(mt) map_dfr(names(M[[sx]][[mt]]), function(key) {
  sl <- strsplit(key, "|", fixed = TRUE)[[1]]
  summ(M[[sx]][[mt]][[key]]) %>% select(-p) %>% mutate(sex = sx, learner = sl[2], strategy = sl[1], metric = mt, .before = 1) })))
write_csv(perf, out("31_v2_performance_by_strategy.csv"))

contr <- map_dfr(names(M), function(sx) map_dfr(names(M[[sx]]), function(mt) map_dfr(c("Linear (splines)", "Random forest"), function(lr)
  map_dfr(names(CONTRASTS), function(cn) {
    a <- paste(CONTRASTS[[cn]][1], lr, sep = "|"); b <- paste(CONTRASTS[[cn]][2], lr, sep = "|")
    if (!all(c(a, b) %in% names(M[[sx]][[mt]]))) return(NULL)
    summ(M[[sx]][[mt]][[a]] - M[[sx]][[mt]][[b]]) %>% mutate(sex = sx, learner = lr, metric = mt, contrast = cn, .before = 1) }))))
write_csv(contr, out("32_v2_paired_contrasts.csv"))

# Share of the information lost by omitting race that the social block recovers: (B - A) / (C - A)
recov <- map_dfr(names(M), function(sx) map_dfr(names(M[[sx]]), function(mt) map_dfr(c("Linear (splines)", "Random forest"), function(lr)
  map_dfr(c("B", "B_nocontext"), function(bb) {
    g <- function(s) M[[sx]][[mt]][[paste(s, lr, sep = "|")]]
    num <- colMeans(g(bb) - g("A")); den <- colMeans(g("C") - g("A"))
    rat <- num[-1] / den[-1]
    tibble(sex = sx, learner = lr, metric = mt, social_block = bb, recovery = num[1] / den[1],
           pct_lo = quantile(rat, 0.025, na.rm = TRUE), pct_hi = quantile(rat, 0.975, na.rm = TRUE),
           race_information = den[1], race_information_lo = den[1] - 1.96 * sd(den[-1]),
           race_information_hi = den[1] + 1.96 * sd(den[-1])) }))))
write_csv(recov, out("33_v2_recovery_of_race_information.csv"))

writeLines(c(capture.output(sessionInfo()), "", paste("CFG:", capture.output(str(CFG)))), out("39_v2_session_info.txt"))
message("All v2 prediction analyses done. Send back files 30_ to 35_ and 39_ from ", OUT_DIR,
        " (not the checkpoints; 34_ is optional, it holds the held-out predictions).")
