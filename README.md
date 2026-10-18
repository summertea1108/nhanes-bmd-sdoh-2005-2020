# Race, social determinants and bone density prediction in US adults

Analysis code for:

Liu J, Jung J, Li Q, Wu L, Wu Q. *Race, social determinants and bone density prediction in US adults: implications for race-neutral risk tools.* Journal of Bone and Mineral Research (submitted).

Main analysis: NHANES 2005–2010, 2013–2014 and 2017–March 2020 (pre-pandemic file), adults aged ≥50 years with a valid femoral neck DXA scan (n = 12,240), seven-item social disadvantage score.
Sensitivity analysis: NHANES 2005–2010, 2013–2014 and 2017–2018 (n = 10,981), eight-item score including home tenure.

## Folder structure

```
analysis/              R scripts that run the analyses (survey regression, multiple imputation,
                       cross-validated prediction, survey-design bootstrap)
analysis/hpc/          Slurm job files used on the Ohio Supercomputer Center (Ascend)
reporting/             Python/JavaScript scripts that turn the result files into the
                       manuscript tables, figures and Supplementary Material
reporting/data/        Aggregate result files (CSV) produced by the R scripts
```

No individual-level data are included. All NHANES data are public and are downloaded
by the scripts from the National Center for Health Statistics
(https://wwwn.cdc.gov/nchs/nhanes/). Participant-level outputs (imputed datasets,
held-out predictions, checkpoints) are not included.

## Software

- R 4.5.2 with survey 4.5, mice 3.19.0, mitools 2.7, ranger 0.18.0, dplyr, tidyr, purrr,
  readr, haven, splines, iml, ggplot2 (full list in `reporting/data/results_2017_2020/sessionInfo.txt` and
  `reporting/data/results_2017_2020/sdoh_v2/39_v2_session_info.txt`)
- Python 3 with pandas, numpy, matplotlib, reportlab (tables, figures, supplement)
- Node.js with the `docx` package (editable Word tables)

Random seeds are set in each R script (configuration list `CFG`).

## Running the analysis

Set the project folder (the scripts default to the authors' OSC path):

```bash
export NHANES_BMD_DIR="/path/to/project"     # folder containing the R scripts
export NHANES_LATEST=P                       # P = 2005–March 2020 (main analysis, default)
                                             # J = 2005–2018 (eight-item sensitivity analysis)
```

Run in this order (each step reads the imputations created by step 1). Each script
has a quick test mode (`QUICK=TRUE`) that runs in minutes.

| Step | Script | Purpose | Main outputs (in `results_2017_2020/`, or `results/` with `NHANES_LATEST=J`) |
|---|---|---|---|
| 1 | `bmd_sdoh_reanalysis.R` | Downloads and harmonizes NHANES, builds the analytic sample, pooled weights and multiple imputations | `imputations.rds`, `table1_weighted_descriptives.csv`, `fig1_flow_counts.csv`, `table_pir_missingness_by_education.csv` |
| 2 | `oi_revision.R` | Education and income associations; sex interactions; included vs excluded participants | `table_crude_vs_adjusted_v2.csv`, `table_svyglm_sex_interactions_v2.csv`, `table_interaction_wald_v2.csv`, `table_included_vs_excluded.csv`, `table_selection_reasons.csv` |
| 3 | `sdoh_bundle.R` | Social disadvantage score: item checks, prevalence, associations, interactions, race and ethnicity coefficients | `sdoh_bundle/00_` to `23_` CSV files |
| 4 | `sdoh_v2_prediction.R` | Four prediction strategies, 10 × 10 cross-validation with within-fold imputation, Rao–Wu bootstrap | `sdoh_v2/30_` to `35_` CSV files |
| 5 | `sdoh_v2_sensitivity.R` | Sensitivity analyses (`SCENARIO` = `main_extra`, `complete`, `mnar`, `cyc1318`, `nomed`, `hip`) | `sdoh_v2_sens/<scenario>/40_` to `42_` CSV files |

`nhanes_cycles.R` is sourced by every script (survey cycles, file names, pooled weights,
cross-cycle variable harmonization). Pooled MEC weights: WTMEC2YR × 2/11.2 for each
2-year cycle and WTMECPRP × 3.2/11.2 for 2017–March 2020.

Sensitivity analysis with item definitions common to all survey periods (Table S13):
run steps 3 and 4 with `SDOH_COMMON=TRUE`; results are written to `sdoh_bundle_common/`
and `sdoh_v2_common/`.

On a Slurm cluster, the job files in `analysis/hpc/` run steps 1, 3, 4 and 5, for example:

```bash
sbatch --export=ALL,QUICK=TRUE run_sdoh_v2.slurm     # quick test
sbatch run_sdoh_v2.slurm                             # full run
sbatch --export=ALL,SDOH_COMMON=TRUE run_sdoh_v2.slurm
```

## Producing the tables, figures and supplement

The `reporting/` scripts read the CSV files in `reporting/data/`
(`results_2017_2020/` = main analysis, `results/` = eight-item 2005–2018 analysis,
`common/` = common-definition analysis).

```bash
cd reporting
python3 main_tables.py main.json           # Tables 1–6 (and figure legends) -> main.json
node build_tables.js main.json Tables.docx # editable Word tables
python3 figures.py out/                    # Figures 1–2 (PDF, EPS, 600-dpi TIFF)
python3 supp_si.py out/                    # Supplementary_Material.pdf (Fig. S1, Tables S1–S14)
```

## Map from result files to the manuscript

| Manuscript item | Result file(s) |
|---|---|
| Table 1 | `table1_weighted_descriptives.csv` |
| Table 2 | `sdoh_bundle/01_item_prevalence_and_score.csv` |
| Table 3 | `table_crude_vs_adjusted_v2.csv` |
| Table 4, Figure 1 | `sdoh_bundle/10_aim1_score_associations.csv`, `12_aim1_interaction_wald.csv` |
| Table 5 | `sdoh_v2/31_v2_performance_by_strategy.csv`, `30_v2_event_counts.csv` |
| Table 6 | `sdoh_bundle/13_aim3_race_coefficients.csv` |
| Figure 2 | `sdoh_v2/32_v2_paired_contrasts.csv` |
| Fig. S1 | `fig1_flow_counts.csv`, `table_selection_reasons.csv` |
| Tables S1–S2 | variable definitions (in `supp_si.py`; harmonization code in `nhanes_cycles.R` and `sdoh_bundle.R`) |
| Table S3 | `sdoh_bundle/00_harmonization_check.csv` |
| Table S4 | `table_svyglm_sex_interactions_v2.csv`, `table_interaction_wald_v2.csv`, `sdoh_bundle/11_aim1_interaction_terms.csv`, `12_aim1_interaction_wald.csv` |
| Table S5 | `sdoh_v2/31_v2_performance_by_strategy.csv` |
| Table S6 | `sdoh_v2/32_v2_paired_contrasts.csv` |
| Table S7 | `sdoh_v2_sens/main_extra/42_main_extra_contrasts.csv` |
| Table S8 | `sdoh_v2/33_v2_recovery_of_race_information.csv` |
| Tables S9–S10 | `sdoh_v2_sens/<scenario>/41_*` and `42_*` (complete, mnar, cyc1318, nomed, main_extra) |
| Table S11 | `sdoh_v2_sens/hip/42_hip_contrasts.csv` |
| Table S12 | `results/` (eight-item, 2005–2018) compared with `results_2017_2020/` |
| Table S13 | `common/sdoh_bundle_common/`, `common/sdoh_v2_common/` |
| Table S14 | `table_included_vs_excluded.csv` |

Text values not shown in a table (for example, PIR missingness by education) come from
`table_pir_missingness_by_education.csv` and `sdoh_bundle/02_item_correlations_weighted.csv`.
