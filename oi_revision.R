# =============================================================================
#  Revision analyses (Osteoporosis International pre-submission review)
#
#   Part 1 (~10-20 min): association analyses re-imputed WITH the outcome.
#          A separate mice imputation (m = 20) that includes femoral neck
#          T-score, used only for the association analyses:
#            - Table 4 (M0 education and PIR mutually adjusted ... M4)
#            - combined-sex survey regression with sex interactions (ESM 2)
#              and the joint Wald test
#          (The prediction analyses keep the outcome-excluding imputation.)
#   Part 2 (~1-2 h): sensitivity analysis S4 (osteoporosis-medication users
#          excluded) rerun for BOTH methods, with the medication variable
#          dropped from the base and full models. Checkpointed.
#   Part 3 (~1 min): age-eligible participants included vs excluded
#          (no valid femoral neck DXA), with reasons.
#
#  Oct 2026 update: cycles come from nhanes_cycles.R. Default NHANES_LATEST=P
#  (2017-March 2020 in place of 2017-2018; results_2017_2020/). NHANES_LATEST=J
#  reproduces the submitted 2017-2018 analysis (results/). Part 2 is not used in
#  the JBMR paper (the v2 "nomed" scenario replaces it), so it is skipped for P
#  unless RUN_PART2=TRUE.
#
#  Needs: <results folder>/imputations.rds, nhanes_cache/*.xpt and nhanes_cycles.R
#  Run in RStudio:
#     setwd("/users/PAS2982/summertea1109/Random Forest ML/NHANES_BMD")
#     source("oi_revision.R")
#  If the session stops, run source("oi_revision.R") again: Part 1 and Part 3
#  are quick, and Part 2 resumes from saved folds.
# =============================================================================
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(purrr); library(readr)
  library(survey); library(mitools); library(mice); library(ranger); library(splines)
})
options(survey.lonely.psu = "adjust", dplyr.summarise.inform = FALSE)
PROJECT_DIR <- Sys.getenv("NHANES_BMD_DIR", "/users/PAS2982/summertea1109/Random Forest ML/NHANES_BMD")
if (dir.exists(PROJECT_DIR)) setwd(PROJECT_DIR)
source("nhanes_cycles.R")   # cycles, file names, pooled weights, results folder
out <- function(f) file.path(RESULTS_ROOT, f)
RUN_PART2 <- toupper(Sys.getenv("RUN_PART2", if (USE_P) "FALSE" else "TRUE")) == "TRUE"

detect_threads <- function() {
  for (v in c("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_JOB_CPUS_PER_NODE")) {
    n <- suppressWarnings(as.integer(sub("\\(.*", "", Sys.getenv(v))))
    if (!is.na(n) && n > 0) return(n)
  }
  n <- tryCatch(length(parallel::mcaffinity()), error = function(e) NA_integer_)
  if (!is.na(n) && n > 0) return(n)
  max(1L, min(4L, parallel::detectCores() - 1L))
}
CFG <- list(seed = 20261001, k_folds = 10, n_repeats = 10, num_trees = 500, tune_trees = 300,
            grid = expand.grid(mtry = c(2, 4, 6), min.node.size = c(10, 50)),
            threads = detect_threads(), m_assoc = 20)
message("Cores used: ", CFG$threads)

imps <- readRDS(out("imputations.rds"))
BASE_COMMON <- c("age", "weight", "height", "race", "cycle", "smoking", "alcohol",
                 "phys_active", "thyroid", "ra", "calcium_supp", "vitd_supp",
                 "prior_fx", "glucocorticoid", "parent_hipfx", "osteo_med")
WOMEN_EXTRA <- c("postmeno", "ever_ht", "ever_preg")
add_pircat <- function(x) mutate(x, pir_cat = cut(pir, c(-Inf, 1.3, 3.5, Inf), right = FALSE,
                                                   labels = c("Low", "Middle", "High")))

# =============================================================================
# Part 1. Association analyses with an outcome-inclusive imputation
# =============================================================================
f_assoc <- out("imputations_assoc_with_outcome.rds")
if (file.exists(f_assoc)) {
  message("Part 1: loading saved outcome-inclusive imputation")
  assoc <- readRDS(f_assoc)
} else {
  assoc <- list()
  for (sx in c("Male", "Female")) {
    message("Part 1: imputing (with T-score) - ", sx)
    X <- imps[[sx]]$imp$data                       # same predictors, same missing values
    dsx <- imps[[sx]]$data
    X$tscore <- dsx$tscore                         # outcome now included (fully observed)
    meth <- make.method(X)
    imp2 <- mice(X, m = CFG$m_assoc, method = meth, seed = CFG$seed + 7,
                 printFlag = FALSE, maxit = 10)
    assoc[[sx]] <- list(imp = imp2, data = dsx)
  }
  saveRDS(assoc, f_assoc)
}
assoc_list <- function(sx) lapply(seq_len(assoc[[sx]]$imp$m), function(i) {
  X <- complete(assoc[[sx]]$imp, i); X$tscore <- NULL
  bind_cols(assoc[[sx]]$data[, c("SEQN", "tscore", "wt_pool", "SDMVPSU", "SDMVSTRA", "sex")], X) %>% add_pircat()
})

adj_sets <- list(
  "M0 education and PIR, mutually adjusted"  = character(0),
  "M1 + age, race and ethnicity, cycle"      = c("age", "race", "cycle"),
  "M2 + weight, height"                      = c("age", "race", "cycle", "weight", "height"),
  "M3 + smoking, alcohol, physical activity" = c("age", "race", "cycle", "weight", "height",
                                                 "smoking", "alcohol", "phys_active"),
  "M4 all base predictors"                   = "ALL")
wald_p <- function(b, V, tt) { w <- as.numeric(t(b[tt]) %*% solve(V[tt, tt]) %*% b[tt])
  c(chisq = w, df = length(tt), p = pchisq(w, length(tt), lower.tail = FALSE)) }
partA <- list()
for (sx in c("Male", "Female")) {
  des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pool, nest = TRUE,
                   data = imputationList(assoc_list(sx)))
  for (mn in names(adj_sets)) {
    cov <- adj_sets[[mn]]
    if (identical(cov, "ALL")) cov <- c(BASE_COMMON, if (sx == "Female") WOMEN_EXTRA)
    f <- as.formula(paste("tscore ~ education + pir_cat", if (length(cov)) paste("+", paste(cov, collapse = " + ")) else ""))
    pooled <- MIcombine(with(des, svyglm(f))); b <- coef(pooled); V <- vcov(pooled)
    terms <- grep("^education|^pir_cat", names(b), value = TRUE)
    partA[[length(partA) + 1]] <- tibble(sex = sx, model = mn, term = terms, est = b[terms],
        se = sqrt(diag(V))[terms]) %>%
      mutate(lo = est - 1.96 * se, hi = est + 1.96 * se, p = 2 * pnorm(-abs(est / se)),
             p_joint_education = wald_p(b, V, grep("^education", terms, value = TRUE))[["p"]],
             p_joint_income    = wald_p(b, V, grep("^pir_cat", terms, value = TRUE))[["p"]])
  }
  message("Part 1, Table 4 done: ", sx)
}
write_csv(bind_rows(partA), out("table_crude_vs_adjusted_v2.csv"))

# Combined-sex model with sex interactions (ESM 2). Adjustment set = all base
# predictors common to both sexes; the female reproductive variables are omitted
# because they are undefined for men.
ml <- assoc_list("Male"); fl <- assoc_list("Female")
comb <- lapply(seq_along(ml), function(i) {
  x <- bind_rows(ml[[i]], fl[[i]] %>% select(-any_of(WOMEN_EXTRA)))
  x$sex <- factor(x$sex, levels = c("Male", "Female")); x })
des_c <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pool, nest = TRUE,
                   data = imputationList(comb))
f_int <- as.formula(paste("tscore ~ sex * (education + pir_cat) +", paste(BASE_COMMON, collapse = " + ")))
pooled <- MIcombine(with(des_c, svyglm(f_int))); b <- coef(pooled); V <- vcov(pooled); se <- sqrt(diag(V))
write_csv(tibble(term = names(b), est = b, se = se, lo = b - 1.96 * se, hi = b + 1.96 * se,
                 p = 2 * pnorm(-abs(b / se))), out("table_svyglm_sex_interactions_v2.csv"))
w <- wald_p(b, V, grep("^sexFemale:", names(b), value = TRUE))
write_csv(tibble(test = "sex x (education + income)", chisq = w[["chisq"]], df = w[["df"]], p = w[["p"]]),
          out("table_interaction_wald_v2.csv"))
message("Part 1 done.")

# =============================================================================
# Part 3. Included vs excluded (run before the long Part 2)
# =============================================================================
rd <- function(f, s, vars) { p <- file.path("nhanes_cache", paste0(nhanes_file(f, s), ".xpt"))   # DEMO_J or P_DEMO
  x <- zap_labels(read_xpt(p)); for (v in setdiff(vars, names(x))) x[[v]] <- NA_real_; x[, vars] }
cyc <- setNames(CYC$cycle, CYC$suffix)   # from nhanes_cycles.R
sel <- map_dfr(names(cyc), function(s) {
  rd("DEMO", s, c("SEQN", "RIDSTATR", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "DMDEDUC2", "INDFMPIR",
                  "WTINT2YR", "WTMEC2YR", "WTINTPRP", "WTMECPRP", "SDMVPSU", "SDMVSTRA")) %>%
    left_join(rd("BMX", s, c("SEQN", "BMXWT", "BMXBMI")), by = "SEQN") %>%
    left_join(rd("DXXFEM", s, c("SEQN", "DXXNKBMD")), by = "SEQN") %>% mutate(cycle = cyc[[s]])
}) %>% filter(RIDAGEYR >= 50) %>%
  mutate(WTINT = coalesce(WTINT2YR, WTINTPRP), WTMEC = coalesce(WTMEC2YR, WTMECPRP),   # P file: *PRP weights
         included = !is.na(DXXNKBMD) & !is.na(WTMEC) & WTMEC > 0,
         group = factor(ifelse(included, "Included", "Excluded"), levels = c("Included", "Excluded")),
         reason = case_when(included ~ "Included",
                            RIDSTATR == 1 ~ "Interviewed only (no MEC examination)",
                            TRUE ~ "Examined, no valid femoral neck DXA"),
         female = as.integer(RIAGENDR == 2),
         race = factor(case_when(RIDRETH1 %in% c(1, 2) ~ "Hispanic", RIDRETH1 == 3 ~ "Non-Hispanic White",
                                 RIDRETH1 == 4 ~ "Non-Hispanic Black", TRUE ~ "Other or multiracial")),
         education = factor(case_when(DMDEDUC2 %in% c(1, 2) ~ "Less than high school", DMDEDUC2 == 3 ~ "High school/GED",
                                      DMDEDUC2 == 4 ~ "Some college/associate", DMDEDUC2 == 5 ~ "College graduate or above"),
                            levels = c("Less than high school", "High school/GED", "Some college/associate", "College graduate or above")),
         pir_cat = cut(INDFMPIR, c(-Inf, 1.3, 3.5, Inf), right = FALSE, labels = c("<1.30", "1.30-3.49", ">=3.50")),
         age80 = as.integer(RIDAGEYR >= 80),
         bmi_cat = cut(BMXBMI, c(-Inf, 25, 30, 40, Inf), right = FALSE, labels = c("<25", "25-29.9", "30-39.9", ">=40")),
         wt_int = pool_weight(WTINT, cycle))     # WTINT2YR / 5, or years-proportional with 2017-2020
write_csv(sel %>% count(reason), out("table_selection_reasons.csv"))
message(sprintf("Part 3 check: included = %d (expected %s), excluded = %d", sum(sel$included),
                if (USE_P) "12240" else "10981", sum(!sel$included)))
des_s <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_int, nest = TRUE, data = sel)
cont <- function(v, lab) {
  m <- svyby(as.formula(paste0("~", v)), ~group, des_s, svymean, na.rm = TRUE)
  s <- sapply(c("Included", "Excluded"), function(g) { k <- sel$group == g & !is.na(sel[[v]])
    x <- sel[[v]][k]; w <- sel$wt_int[k]; sqrt(sum(w * (x - sum(w * x) / sum(w))^2) / sum(w)) })
  p <- tryCatch(svyttest(as.formula(paste0(v, "~group")), des_s)$p.value, error = function(e) NA)
  tibble(variable = lab, level = "mean (SD)",
         Included = sprintf("%.1f (%.1f)", m[[v]][1], s[1]), Excluded = sprintf("%.1f (%.1f)", m[[v]][2], s[2]),
         missing_incl = sum(is.na(sel[[v]][sel$included])), missing_excl = sum(is.na(sel[[v]][!sel$included])), p = p) }
catv <- function(v, lab) {
  sub <- subset(des_s, !is.na(sel[[v]]))
  tb <- svytable(as.formula(paste0("~", v, "+group")), sub); pr <- prop.table(tb, 2)
  n <- table(sel[[v]], sel$group)
  p <- tryCatch(svychisq(as.formula(paste0("~", v, "+group")), sub)$p.value, error = function(e) NA)
  tibble(variable = lab, level = rownames(tb),
         Included = sprintf("%d (%.1f%%)", n[, "Included"], 100 * pr[, "Included"]),
         Excluded = sprintf("%d (%.1f%%)", n[, "Excluded"], 100 * pr[, "Excluded"]),
         missing_incl = sum(is.na(sel[[v]][sel$included])), missing_excl = sum(is.na(sel[[v]][!sel$included])),
         p = c(p, rep(NA, nrow(tb) - 1))) }
selt <- bind_rows(
  tibble(variable = "n", level = "", Included = format(sum(sel$included)), Excluded = format(sum(!sel$included)),
         missing_incl = NA, missing_excl = NA, p = NA),
  cont("RIDAGEYR", "Age, years"), catv("age80", "Aged 80 or older (1 = yes)"), catv("female", "Female (1 = yes)"),
  catv("race", "Race and ethnicity"), catv("education", "Education"), catv("pir_cat", "Family poverty-income ratio"),
  cont("BMXBMI", "BMI, kg/m2"), catv("bmi_cat", "BMI category"), catv("cycle", "Survey cycle"))
write_csv(selt, out("table_included_vs_excluded.csv"))
message("Part 3 done: ", out("table_included_vs_excluded.csv"))

if (!RUN_PART2) {
  message("Part 2 skipped (not used in the JBMR paper; set RUN_PART2=TRUE to run it). Send back: ",
          "table_crude_vs_adjusted_v2.csv, table_svyglm_sex_interactions_v2.csv, table_interaction_wald_v2.csv, ",
          "table_included_vs_excluded.csv, table_selection_reasons.csv from ", RESULTS_ROOT)
} else {
# =============================================================================
# Part 2. S4 without osteoporosis-medication users, medication variable dropped
# =============================================================================
wmetrics <- function(y, yhat, w) { w <- w / sum(w); e <- y - yhat; ybar <- sum(w * y); mse <- sum(w * e^2)
  c(MSE = mse, RMSE = sqrt(mse), MAE = sum(w * abs(e)), R2 = 1 - mse / sum(w * (y - ybar)^2)) }
wauc <- function(risk, case, w) {
  ok <- !is.na(risk) & !is.na(case); risk <- risk[ok]; case <- case[ok]; w <- w[ok]
  if (length(unique(case)) < 2) return(NA_real_)
  u <- sort(unique(risk)); g <- factor(match(risk, u), levels = seq_along(u))
  wn <- tapply(w * (case == 0), g, sum); wn[is.na(wn)] <- 0
  wp <- tapply(w * (case == 1), g, sum); wp[is.na(wp)] <- 0
  sum(wp * (cumsum(wn) - wn + 0.5 * wn)) / (sum(wp) * sum(wn)) }
nb_ci <- function(x, ntr, nts) { x <- x[!is.na(x)]; J <- length(x); m <- mean(x)
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
pred_list <- function(sx) {   # prediction analyses keep the ORIGINAL outcome-excluding imputation
  keep <- c("SEQN", "tscore", "osteoporosis", "wt_pool", "SDMVPSU", "SDMVSTRA", "sex")
  lapply(seq_len(imps[[sx]]$imp$m), function(i) bind_cols(imps[[sx]]$data[, keep], complete(imps[[sx]]$imp, i)))
}
perf <- list()
for (sx in c("Male", "Female")) {
  dl <- lapply(pred_list(sx), function(x) x[x$osteo_med == "No", ])
  base <- setdiff(c(BASE_COMMON, if (sx == "Female") WOMEN_EXTRA), "osteo_med"); full <- c(base, "education", "pir")
  set.seed(CFG$seed + sum(utf8ToInt(paste("S4v2", sx))))
  n <- nrow(dl[[1]]); folds <- list()
  for (r in seq_len(CFG$n_repeats)) { id <- sample(rep(seq_len(CFG$k_folds), length.out = n))
    for (k in seq_len(CFG$k_folds)) folds[[length(folds) + 1]] <- list(rep = r, fold = k, test = which(id == k)) }
  ck <- out(file.path("checkpoints", paste0("s4v2_", sx))); dir.create(ck, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time(); nnew <- 0
  for (fi in seq_along(folds)) {
    f_ck <- file.path(ck, sprintf("fold_%03d.rds", fi))
    if (file.exists(f_ck)) { perf[[length(perf) + 1]] <- readRDS(f_ck); next }
    te <- folds[[fi]]$test; tr1 <- dl[[1]][-te, ]
    params <- tune_rf(tr1, full, tr1$wt_pool / mean(tr1$wt_pool), CFG$seed + fi)
    res <- list()
    for (mi in seq_along(dl)) {
      tr <- dl[[mi]][-te, ]; ts <- dl[[mi]][te, ]; w_tr <- tr$wt_pool / mean(tr$wt_pool); tr$.w <- w_tr
      for (mn in c("base", "full")) {
        preds <- if (mn == "base") base else full
        p <- predict(fit_rf(tr, preds, w_tr, CFG$seed + fi, params), ts, num.threads = CFG$threads)$predictions
        pl <- predict(lm(lm_formula(preds), data = tr, weights = .w), ts)   # no tryCatch: errors must surface
        for (L in list(list("Random forest", p), list("Linear (splines)", pl)))
          res[[length(res) + 1]] <- tibble(sex = sx, rep = folds[[fi]]$rep, fold = folds[[fi]]$fold, imp = mi,
            learner = L[[1]], model = mn, !!!as.list(wmetrics(ts$tscore, L[[2]], ts$wt_pool)),
            AUC_osteoporosis = wauc(-L[[2]], ts$osteoporosis, ts$wt_pool), n_train = nrow(tr), n_test = nrow(ts))
      }
    }
    r <- bind_rows(res); saveRDS(r, f_ck); perf[[length(perf) + 1]] <- r; nnew <- nnew + 1
    mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("  S4 | %s: fold %d / %d done (%.1f min, about %.0f min left)", sx, fi,
                    length(folds), mins, (length(folds) - fi) * mins / nnew))
  }
}
perf <- bind_rows(perf); write_csv(perf, out("cv_raw_performance_s4v2.csv"))
pf <- perf %>% group_by(sex, learner, model, rep, fold) %>%
  summarise(across(c(R2, RMSE, AUC_osteoporosis, n_train, n_test), mean), .groups = "drop")
summ <- pf %>% pivot_longer(c(R2, RMSE, AUC_osteoporosis), names_to = "metric") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  group_by(sex, learner, metric) %>%
  summarise(ci = list(nb_ci(full - base, mean(n_train), mean(n_test))),   # paired differences first
            n = round(mean(n_train + n_test)), base = mean(base, na.rm = TRUE), full_mean = mean(full, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(delta = map_dbl(ci, ~ .x[["est"]]), lo = map_dbl(ci, ~ .x[["lo"]]),
         hi = map_dbl(ci, ~ .x[["hi"]]), p = map_dbl(ci, ~ .x[["p"]])) %>% select(-ci)
write_csv(summ, out("table_s4v2_deltas.csv"))
message("All revision analyses done. Send back: table_crude_vs_adjusted_v2.csv, table_svyglm_sex_interactions_v2.csv, ",
        "table_interaction_wald_v2.csv, table_included_vs_excluded.csv, table_selection_reasons.csv, ",
        "table_s4v2_deltas.csv, cv_raw_performance_s4v2.csv")
}
