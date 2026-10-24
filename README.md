# The double socioeconomic gradient of mental health care in Brazil

Analysis code for the manuscript:

> Dumay, HO. *The double socioeconomic gradient of mental health care in Brazil: access to care resources and unmet need among people aged 15 years or older with a diagnosed mental disorder in the 2019 National Health Survey.* Manuscript under review.

Archived version of this code: [Zenodo DOI, to be added after the first release].

The code estimates income gradients in mental health care among respondents of the 2019 Brazilian National Health Survey (Pesquisa Nacional de Saúde, PNS) who reported a professional diagnosis of depression or another mental illness (n = 10,956). It uses design-based survey regression throughout: the relative index of inequality (RII) comes from modified Poisson models, and the slope index of inequality (SII) comes from linear probability models, both with a weighted income ridit as the exposure.

## Data

The PNS 2019 microdata are public and are distributed by the Brazilian Institute of Geography and Statistics (IBGE). They are **not** included in this repository. The first script downloads them automatically with `PNSIBGE::get_pns()` into the `data_pnsibge/` folder (internet access and about 1 GB of disk space are required). The analyses used the IBGE release `PNS_2019_20220525.zip` (SHA-256 `a789f6b0cbf6057111939fb3d1fe62ad1d5e848ba9c1f7860d7b44daa4be6057`).

No individual-level data are written to this repository. The `results/` folder, which is created locally and excluded by `.gitignore`, contains only aggregate estimates.

## Requirements

- **R ≥ 4.1** with the packages `PNSIBGE`, `survey`, `dplyr`, `ggplot2` and `ranger`:

  ```r
  install.packages(c("PNSIBGE", "survey", "dplyr", "ggplot2", "ranger"))
  ```

- **Python ≥ 3.9** (figures only) with the packages listed in `requirements.txt`:

  ```bash
  pip install -r requirements.txt
  ```

The package versions used for the published results are recorded in `session_info/sessionInfo.txt`.

## How to reproduce

From the repository root:

```bash
Rscript run_all.R                  # all tables (CSV) in results/
python python/make_figures.py      # Figures 1, 2 and S1 in figures/
python python/graphical_abstract.py
```

`run_all.R` runs the scripts in the order below and writes `session_info/sessionInfo.txt`. Each script can also be run on its own, provided the scripts before it have been run.

All paths are relative to a project root, which defaults to the working directory. To keep the data and results elsewhere, set the environment variable `PNS_ROOT`, for example `Sys.setenv(PNS_ROOT = "D:/pns_project")` in R before sourcing the scripts.

## Pipeline

| Order | Script | What it does | Output folder |
|---|---|---|---|
| 1 | `R/01_prepare_cohort.R` | Builds the cohort (Q092 OR Q11006), harmonizes the care items of the depression and other-mental-illness questionnaire blocks, derives the diagnostic hierarchy and the within-cohort income ridit, and saves the survey design | `results/any_mental/`, `results/mental_care_design_anymh.rds` |
| 2 | `R/02_main_analysis.R` | Table 1; adjustment blocks A–E for the co-primary outcomes; shape tests; pre-specified and post hoc interactions; diagnostic-group models; sensitivity analyses; audit of fitted probabilities above 1 with log-binomial refits | `results/anymh/` |
| 3 | `R/03_income_bands.R` | Five per-capita income bands: prevalence, prevalence ratios, shape tests, and band-level results for the depression indicators | `results/income_bands_fine/` |
| 4 | `R/04_depression_cascade.R` | Exploratory depression care cascade (treatment, continuity, reaching a specialist, reasons for non-regular follow-up, perceived remission with current symptoms), each with a domain-specific income ridit | `results/depression_cascade_domain_rank/` |
| 5 | `R/05_detection_gap.R` | Exploratory detection gap among PHQ-9 screen-positive respondents | `results/detection_gap_domain_rank/` |
| 6 | `R/06_depression_model_checks.R` | Fitted-probability audit and log-binomial refits for the 38 exploratory depression models | `results/depression_prediction_check/` |
| 7 | `R/07_random_forest.R` | Exploratory random-forest check (supplementary; fixed seeds) | `results/forest_anymh/` |
| 8 | `R/08_referral_failure_blockD.R` | Referral failure re-estimated with the block-D covariates of the other Table 6 rows; the earlier model with public-sector care (script 04) is kept as a sensitivity model | `results/referral_failure_blockD/` |
| Figures | `python/make_figures.py`, `python/graphical_abstract.py` | Figures 1, 2, S1 and the graphical abstract, read from the CSVs above | `figures/` |

## Where each manuscript result comes from

| Manuscript item | Script(s) | Main CSV file(s) |
|---|---|---|
| Table 1 | 02, 03 | `anymh/01_tabela1_ponderada.csv`; income-band row: `income_bands_fine/A1_distribuicao_faixas_tabela1.csv` |
| Table 2 | 02 | `anymh/06_blocos_coorte.csv` |
| Table 3 | 03 (bands), 02 (quintiles) | `income_bands_fine/A2_prevalencia_por_faixa.csv`, `A3_PR_faixas_blocos_B_D.csv`, `A4_forma_faixas_novas.csv`; `anymh/03b_prevalencia_por_quintil.csv`, `07b_PR_quintis_blocoB.csv`, `08_forma_do_gradiente.csv` |
| Table 4 | 02, 03 | `anymh/04_prevalencia_por_grupo.csv`, `09_interacoes_pre_especificadas.csv`, `10_interacao_renda_x_grupo.csv`, `11_RII_por_grupo_modelo_interacao.csv`; `income_bands_fine/A5_grupo_x_faixa.csv` |
| Table 5 | 02 | `anymh/14_sensibilidades_coorte.csv`, `15_familia_de_modelo.csv`, `16_checagem_previsoes_e_logbinomial.csv` |
| Table 6 | 04, 05, 03, 08 | `depression_cascade_domain_rank/*_PR_income.csv` and `*_by_income.csv`; `detection_gap_domain_rank/03_lacuna_RII_renda.csv`; `income_bands_fine/B1_depressao_prevalencia_por_faixa.csv`; referral-failure row: `referral_failure_blockD/33_falha_encaminhamento_modelos.csv` (M0 minimal, M1 adjusted) |
| Figure 1 | 02 → `make_figures.py` | `anymh/03b_…`, `07b_…`, `08_…` |
| Figure 2 | 02, 04, 05, 08 → `make_figures.py` | see the docstring of `make_figures.py` |
| Graphical abstract | 03 → `graphical_abstract.py` | `income_bands_fine/A0_…`, `A1_…`, `B1_…` |
| Table S1 | 01 | `any_mental/audit_concordance_ambos.csv`, `audit_crosstab_*.csv` |
| Tables S2–S3 | 02 | `anymh/12_blocos_por_grupo.csv`, `13_sensibilidades_por_grupo.csv` |
| Table S4 | 05 | `detection_gap_domain_rank/01_…`, `02_…`, `03_lacuna_RII_renda.csv` |
| Table S5 | 02 | `anymh/15_familia_de_modelo.csv`, `16_checagem_previsoes_e_logbinomial.csv`, `17_perfil_previsoes_acima_de_1.csv` |
| Tables S6–S7, Figure S1 | 07 → `make_figures.py` | `forest_anymh/49b_*.csv`, `50b_*.csv`, `50c_*.csv` |
| Table S8 | 03 | `income_bands_fine/A0_contagens_classes_e_faixas.csv` |
| Table S9 | 06, 08 | `depression_prediction_check/32_checagem_previsoes_depressao.csv`; referral failure, block D: `referral_failure_blockD/33_falha_encaminhamento_modelos.csv` (log-binomial columns) |
| Table S10 | 04, 05 | `depression_cascade_domain_rank/cascade_stages.csv`, `critical_gap_PR_income.csv`; detection rows: `detection_gap_domain_rank/01_…`, `02_…` |

## Notes

- **Language.** Code comments, some variable labels and the names of the output CSV files are in Portuguese (for example, `desfecho` = outcome, `bloco` = adjustment block, `faixa` = income band, `renda` = income). The CSV column names follow the same convention.
- **Development names.** The scripts were renamed for this repository. Comments inside them refer to each other by their development names:

  | Repository name | Development name |
  |---|---|
  | `01_prepare_cohort.R` | `01c_prepare_any_mental.R` |
  | `02_main_analysis.R` | `30_anymh_analysis.R` |
  | `03_income_bands.R` | `31_faixas_renda_finas.R` |
  | `04_depression_cascade.R` | `16b_depression_cascade_domain_rank.R` |
  | `05_detection_gap.R` | `18b_detection_gap_domain_rank.R` |
  | `06_depression_model_checks.R` | `32_checagem_previsoes_depressao.R` |
  | `07_random_forest.R` | `15b_forest_anymh.R` |
  | `08_referral_failure_blockD.R` | `33_falha_encaminhamento_blocoD.R` |

  `30b_checagem_previsoes.R`, which re-ran only section 9 of `30_anymh_analysis.R`, is not included because `02_main_analysis.R` contains that section. Development notes in the comments were shortened, and one variable not used by any analysis was removed from `01_prepare_cohort.R`. Apart from these changes, the headers and the project-root line, the code is the version that produced the published results.
- **Survey design.** All estimates use the calibrated design returned by `PNSIBGE::pns_design()`, with `survey.lonely.psu = "adjust"` and `survey.adjust.domain.lonely = TRUE`. Analyses of subpopulations are domain analyses on the full design, not subsets of the data frame.

## License

Code: MIT License (see `LICENSE`). The PNS 2019 microdata are subject to IBGE's terms of use.

## Citation

If you use this code, please cite the article above and this repository (see `CITATION.cff`).
