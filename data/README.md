# Data dictionary

All files are UTF-8 CSV or GeoPackage files, derived from public sources (see the main README) by `R/01_download_raw_data.R`. Municipality codes are 6-digit IBGE codes; Maranhão codes start with `21`.

## `oncology_panel_c53_ma.csv` — Painel de Oncologia (PO-SUS)

There is one row per patient and cancer site (first diagnosis): malignant cervical neoplasms (ICD-10 C53) in residents of Maranhão, diagnosed in 2013–2025. Original codes are kept as text.

| Variable | Description |
|---|---|
| `year_diagnosis` | Year of diagnosis |
| `mun_residence` | Municipality of residence |
| `mun_diagnosis` | Municipality of the facility where the diagnosis was recorded |
| `uf_treatment` | State (IBGE code) of the facility where the first treatment was recorded; empty if none |
| `mun_treatment` | Municipality of the facility where the first treatment was recorded; empty if none |
| `icd10` | Detailed diagnosis (ICD-10, 3 characters) |
| `diagnosis_category` | `01` = malignant neoplasm (Law 12,732/2012) |
| `sex` | `F` female, `M` male |
| `age` | Age at diagnosis in years (`999` = unknown) |
| `stage_code` | Recorded stage: `0`, `1` = I, `2` = II, `3` = III, `4` = IV, `5` = not applicable, `9` = unknown |
| `treatment_code` | First treatment: `1` surgery, `2` chemotherapy, `3` radiotherapy, `4` chemotherapy + radiotherapy, `5` no treatment information |
| `time_to_treatment_raw` | Interval between diagnosis and first treatment: sign (`+` after, `-` before the pathology report) followed by the number of days; `99999` = no treatment information |

In PO-SUS, stage is captured from chemotherapy and radiotherapy authorisations. Cases whose first treatment was surgery are coded `5` and cases without treatment `9`, so missing stage depends on the treatment modality.

## `mortality_c53_c55_ma.csv` — Sistema de Informações sobre Mortalidade (SIM)

There is one row per death of a woman resident in Maranhão, 2013–2024, with underlying cause C53 (cervix uteri), C54 (corpus uteri) or C55 (uterus, part unspecified).

| Variable | Description |
|---|---|
| `year` | Year of death |
| `age_years` | Age at death in years |
| `underlying_cause` | ICD-10 underlying cause (3 characters) |
| `mun_residence` | Municipality of residence |
| `race` | Race/skin colour: White, Black, Asian, Brown (parda), Indigenous, Unknown |
| `education` | Years of schooling: None, 1–3, 4–7, 8–11, 12 or more, Unknown |
| `marital_status` | Single, Married, Widowed, Legally separated/divorced, Consensual union, Unknown |
| `place_of_death` | Hospital, Other health facility, Home, Public road, Other, Indigenous village, Unknown |

## `female_population_ma.csv`

Female population by municipality (`mun_code`), year (2013–2024) and 5-year age group (`age5`, `00-04` … `80+`), from the Ministry of Health/IBGE estimates.

## `municipalities_ma.csv`

| Variable | Description |
|---|---|
| `mun_code` | Municipality code |
| `health_region_code`, `health_region` | Health region (19 in Maranhão) |
| `health_macroregion` | Health macroregion |
| `total_pop_mean_2013_2024` | Mean total population, 2013–2024 (used for municipal size categories) |

## `municipal_seats_brazil.csv`

Coordinates (SIRGAS 2000, longitude/latitude) of the seats of all Brazilian municipalities (IBGE, 2010), used to compute distances. The file covers every municipality because some women were treated outside Maranhão.

## `maps/`

| File | Content |
|---|---|
| `ma_municipalities.gpkg` | Municipal boundaries of Maranhão (IBGE, 2022) |
| `ma_health_regions.gpkg` | Health regions of Maranhão, dissolved from municipal boundaries |
| `brazil_states.gpkg` | State boundaries of Brazil (IBGE, 2020), used as map background |

## License

The derived data are released under CC BY 4.0. The original data are public information from the Brazilian Ministry of Health (DATASUS) and IBGE.
