# =============================================================================
#  SDoH and race: sensitivity and secondary prediction analyses, v2 procedure
#
#  Re-runs every prediction analysis that the manuscript reports outside the
#  main four-strategy comparison, using exactly the v2 procedure from
#  sdoh_v2_prediction.R (imputation inside each training fold with only the
#  variables each strategy may use; performance from all held-out predictions
#  per repetition; Rao-Wu survey bootstrap, 500 replicates, plus
#  between-repetition variance). Same fold seeds as the main v2 run.
#
#  Strategies in every scenario: A clinical, B + SDoH block, C + race,
#  D + race + SDoH block, AL1 clinical + education + income (PIR),
#  CL1 clinical + race + education + income.
#
#  SCENARIO (environment variable), one per Slurm array task:
#    main_extra  full sample; also D + unemployment
#    complete    complete cases only (no imputation needed)
#    mnar        imputed PIR shifted by -0.5 and -1.0 (missing not at random)
#    cyc1318     2013-2014 and the latest cycle only (2017-2018 or 2017-March 2020)
#    nomed       osteoporosis-medication users excluded, medication variable removed
#    hip         total hip T-score as the outcome
#
#  Inputs: <results folder>/imputations.rds, nhanes_cache/ and nhanes_cycles.R
#  <results folder> = results/ (NHANES_LATEST=J) or results_2017_2020/ (default, P)
#  Run:   sbatch run_sdoh_v2_sens.slurm          (all six scenarios as an array)
#  Quick: sbatch --export=ALL,QUICK=TRUE run_sdoh_v2_sens.slurm
#  Outputs: <results folder>/sdoh_v2_sens/<scenario>/40_ to 42_ and 49_ files.
#  Folds are checkpointed, so re-submitting resumes.
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
CFG <- list(seed = 20261010, k_folds = 10, n_repeats = 10, m = 3, maxit = 10, B = 500,
            num_trees = 500, tune_trees = 300,
            grid = expand.grid(mtry = c(2, 4, 6), min.node.size = c(10, 50)), threads = detect_threads())
if (QUICK) CFG <- modifyList(CFG, list(k_folds = 3, n_repeats = 2, m = 1, maxit = 5, B = 50,
                                       num_trees = 100, tune_trees = 100,
                                       grid = expand.grid(mtry = c(2, 4), min.node.size = 20)))
SCENARIO <- Sys.getenv("SCENARIO", "main_extra")
SCENARIOS <- c("main_extra", "complete", "mnar", "cyc1318", "nomed", "hip")
if (!SCENARIO %in% SCENARIOS) stop("SCENARIO must be one of: ", paste(SCENARIOS, collapse = ", "))
OUT_DIR <- file.path(RESULTS_ROOT, if (QUICK) "sdoh_v2_sens_quick" else "sdoh_v2_sens", SCENARIO)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)
message(sprintf("SDoH v2 sensitivity [%s] | %s | cores: %d | %d x %d CV | m = %d per fold | %d bootstrap replicates",
                SCENARIO, if (QUICK) "QUICK TEST" else "FULL RUN", CFG$threads, CFG$k_folds, CFG$n_repeats, CFG$m, CFG$B))

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
      poor_access = yesno(HUQ030 == 2 | care_type %in% 3, HUQ030 %in% c(1, 3)),
      # Unpartnered or living alone (DMDMARTL: 1 married, 6 living with partner)
      alone       = yesno(DMDMARTL %in% 2:5 | lives_alone_raw %in% TRUE,
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
# 3. Scenario: which participants, which outcome, which predictors
#    All scenarios use the v2 procedure: imputation inside each training fold
#    with only the variables each strategy may use, performance from all
#    held-out predictions per repetition, survey-design bootstrap CIs.
# =============================================================================
#   main_extra : full sample; education + income added (Layer 1), unemployment added
#   complete   : participants with no missing value in any predictor (no imputation)
#   mnar       : imputed PIR shifted by -0.5 and -1.0 (missing not at random)
#   cyc1318    : 2013-2014 and the latest cycle only
#   nomed      : users of osteoporosis medication excluded; medication variable removed
#   hip        : total hip T-score (NHANES III total femur reference) as the outcome
HIP_MEAN <- 0.942; HIP_SD <- 0.122
dat_full <- dat                           # full sample: used to build the survey replicate weights (domain estimation)
if (SCENARIO == "hip") {
  hip <- map_dfr(CYC$suffix, function(s) read_xpt_cached("DXXFEM", s, CYC$year[CYC$suffix == s]) %>%   # P_DXXFEM for 2017-2020
                   select(SEQN, DXXOFBMD)) %>% filter(!is.na(DXXOFBMD))
  dat <- map(dat, function(d) d %>% inner_join(hip, by = "SEQN") %>%
               mutate(tscore = (DXXOFBMD - HIP_MEAN) / HIP_SD, osteoporosis = as.integer(tscore <= -2.5)))
}
DROP_MED <- SCENARIO == "nomed"
CLIN <- function(sx) setdiff(c(setdiff(BASE_COMMON, "race"), if (sx == "Female") WOMEN_EXTRA),
                             if (DROP_MED) "osteo_med")
L1 <- c("education", "pir")
info_sets <- function(sx) { a <- CLIN(sx)
  s <- list(A = a, B = c(a, BLOCK), C = c(a, "race"), D = c(a, "race", BLOCK),
            AL1 = c(a, L1), CL1 = c(a, "race", L1))
  if (SCENARIO == "main_extra") s$Du <- c(a, "race", BLOCK, "unemployed")
  s }
strategies <- function(sx) { a <- CLIN(sx)
  st <- list(A   = list(set = "A",   preds = a),
             B   = list(set = "B",   preds = c(a, BLOCK)),
             C   = list(set = "C",   preds = c(a, "race")),
             D   = list(set = "D",   preds = c(a, "race", BLOCK)),
             AL1 = list(set = "AL1", preds = c(a, L1)),
             CL1 = list(set = "CL1", preds = c(a, "race", L1)))
  if (SCENARIO == "main_extra") st$D_unemp <- list(set = "Du", preds = c(a, "race", BLOCK, "unemployed"))
  if (SCENARIO == "mnar") for (dl in c("m05", "m10")) for (s in c("B", "D", "AL1", "CL1"))
    st[[paste0(s, "_", dl)]] <- list(set = paste0(s, "_", dl), preds = st[[s]]$preds)
  st }
CONTRASTS <- list(
  "D - C: social block added to model with race"    = c("D", "C"),
  "B - A: social block added to model without race" = c("B", "A"),
  "A - C: race omitted"                             = c("A", "C"),
  "B - C: social block in place of race"            = c("B", "C"),
  "CL1 - C: education and income added to model with race"    = c("CL1", "C"),
  "AL1 - A: education and income added to model without race" = c("AL1", "A"))
if (SCENARIO == "main_extra") CONTRASTS[["D_unemp - D: unemployment added to model with race and social block"]] <- c("D_unemp", "D")
if (SCENARIO == "mnar") for (dl in c("m05", "m10")) {
  lab <- if (dl == "m05") " [PIR -0.5]" else " [PIR -1.0]"
  CONTRASTS[[paste0("D - C", lab)]]   <- c(paste0("D_", dl), "C")
  CONTRASTS[[paste0("B - A", lab)]]   <- c(paste0("B_", dl), "A")
  CONTRASTS[[paste0("CL1 - C", lab)]] <- c(paste0("CL1_", dl), "C")
  CONTRASTS[[paste0("AL1 - A", lab)]] <- c(paste0("AL1_", dl), "A")
}
MNAR_SHIFT <- c(m05 = -0.5, m10 = -1.0)
# Imputation seed index per information set: A-D use the same index as sdoh_v2_prediction.R,
# so in the full-sample scenarios strategies A-D reproduce the main v2 results exactly (a built-in check)
SEED_IDX <- c(A = 1, B = 2, Bnc = 3, C = 4, D = 5, Dnc = 6, AL1 = 7, CL1 = 8, Du = 9)

# Participants kept in this scenario
keep_rows <- function(d, sx) {
  if (SCENARIO == "cyc1318") d <- d %>% filter(cycle %in% CYCLES_2013ON)   # 2013-2014 + latest cycle
  if (SCENARIO == "nomed")   d <- d %>% filter(osteo_med %in% "No")
  if (SCENARIO == "complete") { v <- unique(unlist(info_sets(sx))); d <- d[stats::complete.cases(d[, v]), ] }
  droplevels(d)
}
dat <- imap(dat, keep_rows)
message(sprintf("Scenario %s: n = %d men, %d women", SCENARIO, nrow(dat$Male), nrow(dat$Female)))
META <- c("SEQN", "tscore", "osteoporosis", "wt_pool")

# Derived variables are always rebuilt from the completed components
derive <- function(x) {
  yn <- function(z) factor(ifelse(z, "Yes", "No"), levels = c("No", "Yes"))
  if ("pir" %in% names(x))       x$low_pir <- yn(x$pir < 1.3)
  if ("education" %in% names(x)) x$low_edu <- yn(x$education == "Less than high school")
  if (all(SCORE_ITEMS %in% names(x))) x$sdoh_score <- rowSums(sapply(SCORE_ITEMS, function(v) x[[v]] == "Yes"))
  x
}
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
# 4. Cross-validation (same fold seeds as sdoh_v2_prediction.R; checkpointed)
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
  base_sets <- unique(sub("_m(05|10)$", "", map_chr(st, "set")))
  comp <- list()
  for (s in base_sets) {
    raw <- run_mice(d, sets[[s]], is_te, CFG$m, CFG$seed + 1000 * fi + SEED_IDX[[s]])
    comp[[s]] <- lapply(raw, function(X) derive(bind_cols(d[, META], X)))
    if (SCENARIO == "mnar" && "pir" %in% sets[[s]]) {
      miss <- is.na(d$pir)
      for (dl in names(MNAR_SHIFT)) comp[[paste0(s, "_", dl)]] <- lapply(raw, function(X) {
        X$pir[miss] <- pmin(pmax(X$pir[miss] + MNAR_SHIFT[[dl]], 0), 5)
        derive(bind_cols(d[, META], X)) })
    }
  }
  res <- list()
  for (sn in names(st)) {
    s <- st[[sn]]; L <- comp[[s$set]]
    tr1 <- L[[1]][!is_te, ]
    params <- tune_rf(tr1, s$preds, tr1$wt_pool / mean(tr1$wt_pool), CFG$seed + fi)
    p_lm <- p_rf <- matrix(NA_real_, sum(is_te), length(L))
    for (i in seq_along(L)) {
      tr <- L[[i]][!is_te, ]; ts <- L[[i]][is_te, ]; tr$.w <- tr$wt_pool / mean(tr$wt_pool)
      p_lm[, i] <- predict(lm(lm_formula(s$preds), data = tr, weights = .w), ts)
      p_rf[, i] <- predict(fit_rf(tr, s$preds, tr$.w, CFG$seed + fi, params), ts, num.threads = 1)$predictions
    }
    res[[length(res) + 1]] <- tibble(sex = sx, rep = fold$rep, fold = fold$fold, SEQN = d$SEQN[is_te],
                                     strategy = sn, learner = "Linear (splines)", pred = rowMeans(p_lm))
    res[[length(res) + 1]] <- tibble(sex = sx, rep = fold$rep, fold = fold$fold, SEQN = d$SEQN[is_te],
                                     strategy = sn, learner = "Random forest", pred = rowMeans(p_rf))
  }
  r <- bind_rows(res); saveRDS(r, ck_file(sx, fi))
  message(sprintf("  %s | %s fold %d / %d done in %.1f min (%s)", SCENARIO, sx, fi, length(FOLDS[[sx]]),
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

# =============================================================================
# 5. Performance per repetition; survey bootstrap replicate weights are built on
#    the full sample of each sex and restricted to the scenario's participants
#    (domain estimation), so subsets keep the full design.
# =============================================================================
wm_mat <- function(y, yhat, W) {
  Wn <- sweep(W, 2, colSums(W), "/"); mse <- colSums(Wn * (y - yhat)^2)
  ybar <- colSums(Wn * y); sst <- colSums(Wn * sweep(matrix(y, length(y), ncol(W)), 2, ybar)^2)
  list(R2 = 1 - mse / sst, RMSE = sqrt(mse), MAE = colSums(Wn * abs(y - yhat))) }
wauc_mat <- function(risk, case, W) {
  g <- match(risk, sort(unique(risk)))
  wn <- rowsum(W * (case == 0), g, reorder = TRUE); wp <- rowsum(W * (case == 1), g, reorder = TRUE)
  cn <- apply(wn, 2, cumsum); if (!is.matrix(cn)) cn <- matrix(cn, nrow = nrow(wn))
  colSums(wp * (cn - 0.5 * wn)) / (colSums(wp) * colSums(wn)) }

counts <- list(); M <- list()
for (sx in names(FOLDS)) {
  full <- dat_full[[sx]] %>% mutate(stratum = interaction(cycle, SDMVSTRA, drop = TRUE))
  des <- svydesign(ids = ~SDMVPSU, strata = ~stratum, weights = ~wt_pool, nest = TRUE, data = full)
  set.seed(CFG$seed + 99)
  rdes <- as.svrepdesign(des, type = "subbootstrap", replicates = CFG$B)
  Wf <- cbind(full$wt_pool, as.matrix(weights(rdes, type = "analysis")))
  dd <- dat[[sx]]; W <- Wf[match(dd$SEQN, full$SEQN), , drop = FALSE]
  if (anyNA(W)) stop("Scenario participants not found in the full sample")
  lowbm <- as.integer(dd$tscore < -1)
  counts[[sx]] <- tibble(scenario = SCENARIO, sex = sx, n = nrow(dd), n_osteoporosis = sum(dd$osteoporosis == 1),
                         n_low_bm_or_osteoporosis = sum(lowbm == 1),
                         weighted_prevalence_pct = 100 * weighted.mean(dd$osteoporosis == 1, dd$wt_pool))
  M[[sx]] <- list()
  Ps <- P %>% filter(sex == sx)
  for (key in unique(paste(Ps$strategy, Ps$learner, sep = "|"))) {
    sl <- strsplit(key, "|", fixed = TRUE)[[1]]
    for (r in seq_len(CFG$n_repeats)) {
      pr <- Ps %>% filter(strategy == sl[1], learner == sl[2], rep == r)
      yhat <- pr$pred[match(dd$SEQN, pr$SEQN)]
      if (anyNA(yhat) || nrow(pr) != nrow(dd)) stop("Repetition ", r, " of ", key, " does not predict every participant once")
      mm <- wm_mat(dd$tscore, yhat, W)
      mm$AUC_osteoporosis <- wauc_mat(-yhat, dd$osteoporosis, W)
      mm$AUC_low_bm_or_osteoporosis <- wauc_mat(-yhat, lowbm, W)
      for (mt in names(mm)) M[[sx]][[mt]][[key]] <- rbind(M[[sx]][[mt]][[key]], mm[[mt]])
    }
  }
  message("Evaluation done: ", sx)
}
summ <- function(X) {
  est <- mean(X[, 1]); boot <- colMeans(X[, -1, drop = FALSE])
  cv_sd <- if (nrow(X) > 1) sd(X[, 1]) else NA_real_
  se <- sqrt(var(boot, na.rm = TRUE) + ifelse(is.na(cv_sd), 0, cv_sd^2 / nrow(X)))
  tibble(est = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se, p = 2 * pnorm(-abs(est / se)),
         n_boot_ok = sum(is.finite(boot)), cv_sd_across_repetitions = cv_sd) }
perf <- map_dfr(names(M), function(sx) map_dfr(names(M[[sx]]), function(mt) map_dfr(names(M[[sx]][[mt]]), function(key) {
  sl <- strsplit(key, "|", fixed = TRUE)[[1]]
  summ(M[[sx]][[mt]][[key]]) %>% select(-p) %>%
    mutate(scenario = SCENARIO, sex = sx, learner = sl[2], strategy = sl[1], metric = mt, .before = 1) })))
contr <- map_dfr(names(M), function(sx) map_dfr(names(M[[sx]]), function(mt) map_dfr(c("Linear (splines)", "Random forest"), function(lr)
  map_dfr(names(CONTRASTS), function(cn) {
    a <- paste(CONTRASTS[[cn]][1], lr, sep = "|"); b <- paste(CONTRASTS[[cn]][2], lr, sep = "|")
    if (!all(c(a, b) %in% names(M[[sx]][[mt]]))) return(NULL)
    summ(M[[sx]][[mt]][[a]] - M[[sx]][[mt]][[b]]) %>%
      mutate(scenario = SCENARIO, sex = sx, learner = lr, metric = mt, contrast = cn, .before = 1) }))))
write_csv(bind_rows(counts), out(paste0("40_", SCENARIO, "_counts.csv")))
write_csv(perf,  out(paste0("41_", SCENARIO, "_performance.csv")))
write_csv(contr, out(paste0("42_", SCENARIO, "_contrasts.csv")))
writeLines(c(capture.output(sessionInfo()), "", paste("CFG:", capture.output(str(CFG)))),
           out(paste0("49_", SCENARIO, "_session_info.txt")))
message("Scenario ", SCENARIO, " done. Send back files 40_ to 42_ and 49_ from ", OUT_DIR, " (not the checkpoints).")
