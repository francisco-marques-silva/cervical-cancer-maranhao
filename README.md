# Cervical cancer care in Maranhão, Brazil: code and data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23306601.svg)](https://doi.org/10.5281/zenodo.23306601)

This repository contains the data and R code that reproduce all analyses of two studies of cervical cancer (ICD-10 C53) in women living in Maranhão, Brazil. Both studies use public data from the Brazilian Unified Health System (SUS):

1. **Distance to oncology hubs and access to cervical cancer care in Maranhão, Brazil: a population-based study combining individual and ecological analyses, 2013–2024.** It asks whether the distance from the municipality of residence to the nearest oncology hub is associated with late-stage disease, treatment delay (>60 days) and lack of registered treatment, at the individual level, and with case registration and mortality, at the municipal level.
2. **Desigualdades entre regiões de saúde no diagnóstico e no tratamento do câncer do colo do útero no Maranhão de 2020 a 2024: estudo ecológico** (Regional inequalities in the diagnosis and treatment of cervical cancer in Maranhão, 2020–2024: ecological study). It compares registration and mortality rates, timeliness indicators and the profile of deaths across the state's 19 health regions. The article is in Portuguese; the outputs in this repository are in English.

## Repository structure

```
.
├── run_all.R                     # reproduces every analysis from data/
├── R/
│   ├── 00_setup.R                # packages, paths, theme, constants
│   ├── functions.R               # age-standardised rates, Wilson CI, distances, EB smoothing
│   ├── 01_download_raw_data.R    # optional: rebuilds data/ from DATASUS and IBGE
│   ├── 02_prepare_data.R         # analytic datasets
│   ├── 03_article1_distance.R    # article 1: models, tables and figures
│   └── 04_article2_regions.R     # article 2: rates, tables and figures
├── data/                         # analysis inputs (see data/README.md)
│   └── maps/                     # GeoPackage boundaries (IBGE)
└── output/                       # created when the scripts run (not versioned)
```

## How to reproduce

The analyses were run with R 4.6.1 on Windows 11. From the repository root, run:

```r
# install.packages(c("tidyverse", "sf", "gtsummary", "flextable", "officer",
#                    "sandwich", "lmtest", "MASS", "patchwork", "scales"))
source("run_all.R")
```

`run_all.R` writes the following to `output/`:

- **Tables:** Word files (`output/tables/`).
- **Figures:** PNG, PDF and 600-dpi TIFF (`output/figures/`).
- **Key estimates:** CSV files with the main results, for example `output/article1_key_estimates.csv`.
- **Software environment:** `output/session_info.txt`.

The files in `data/` are enough to reproduce every result. To rebuild them from the original sources, run `R/01_download_raw_data.R`. It downloads about 250 MB and also requires the packages `microdatasus`, `brpop` and `geobr`. Set the environment variable `RAW_DIR` to reuse previously downloaded `.dbc` files.

### Package versions used

| Package | Version | Package | Version |
|---|---|---|---|
| tidyverse | 2.0.0 | sandwich | 3.1.3 |
| dplyr | 1.2.1 | lmtest | 0.9.40 |
| ggplot2 | 4.0.3 | MASS | 7.3.65 |
| sf | 1.1.3 | patchwork | 1.3.2 |
| gtsummary | 2.5.1 | scales | 1.4.0 |
| flextable | 0.10.0 | microdatasus* | 3.0.0 |
| officer | 0.7.6 | brpop* | 0.7.0 |
| geobr* | 2.1.0 | | |

\* only needed for `R/01_download_raw_data.R`.

## Data sources

| Source | Content | Access |
|---|---|---|
| Painel de Oncologia (PO-SUS), Ministry of Health | Cancer diagnoses and first treatments in the SUS, files by year of diagnosis (2013–2025) | DATASUS FTP: `ftp.datasus.gov.br/dissemin/publicos/PAINEL_ONCOLOGIA/DADOS/` |
| Sistema de Informações sobre Mortalidade (SIM) | Death certificates, Maranhão, 2013–2024 | DATASUS, via the `microdatasus` R package |
| Female population estimates | By municipality, year and 5-year age group | Ministry of Health/IBGE, via the `brpop` R package |
| Territorial data | Municipal boundaries, municipal seats and health regions | IBGE, via the `geobr` R package |

All sources are public and de-identified. The files in `data/` contain only the variables needed for the analyses. They do not include patient identifiers (the encrypted national health card number), dates of birth, exact dates of diagnosis, treatment or death, or facility codes. See `data/README.md` for the data dictionary.

## Main methods

- **Outcomes (PO-SUS):**
  - late stage (FIGO III–IV vs I–II), among women whose first treatment was chemotherapy and/or radiotherapy, because stage is only recorded for them in PO-SUS;
  - first treatment more than 60 days after diagnosis (Brazilian Law 12,732/2012);
  - no registered treatment.
- **Exposure:** great-circle distance from the seat of the municipality of residence to the seat of the nearest oncology hub. Hubs are municipalities that delivered ≥1% of first treatments: São Luís, Imperatriz and Teresina.
- **Individual-level models:** Poisson regression with cluster-robust variance (prevalence ratios), adjusted for age group, period and municipal population size, plus sensitivity analyses S1–S5.
- **Ecological models:** indirect standardisation with state age-specific rates; negative binomial models for standardised registration and mortality ratios; quasi-Poisson models for deaths per registered case.
- **Rates:** direct age standardisation to the WHO World Standard Population, with gamma confidence intervals.

## License

- Code: MIT License (see `LICENSE`).
- Derived data in `data/`: Creative Commons Attribution 4.0 International (CC BY 4.0). The original data are public information from the Brazilian Ministry of Health and IBGE.

## Citation

If you use this code or data, please cite the archived version:

> Vieira MOD, Marques Silva Neto F, Oliveira RS, Leite IRL. Cervical cancer care in Maranhão, Brazil: code and data (v1.0.0) [software]. Zenodo; 2026. https://doi.org/10.5281/zenodo.23306602

The DOI above refers to version 1.0.0, the version used in the articles. The concept DOI [10.5281/zenodo.23306601](https://doi.org/10.5281/zenodo.23306601) always resolves to the latest version. Citation metadata are in `CITATION.cff`.
