# =============================================================================
# 02_prepare_data.R
# Builds the analytic datasets from the files in data/:
#   cases.rds          cervical cancer cases (PO-SUS), residents of Maranhão
#   deaths.rds         deaths from C53/C54/C55 (SIM), residents of Maranhão
#   population.rds     female population by municipality, year and age group
#   municipalities.rds municipalities with health region, size and distance to hubs
#   hubs.rds           oncology treatment hubs
# =============================================================================

source(file.path("R", "00_setup.R"))
source(file.path("R", "functions.R"))

chr <- cols(.default = col_character())
po_raw <- read_csv(file.path(dir_data, "oncology_panel_c53_ma.csv"), col_types = chr)
deaths_raw <- read_csv(file.path(dir_data, "mortality_c53_c55_ma.csv"),
                       col_types = cols(year = "i", age_years = "d", .default = "c"))
pop_raw <- read_csv(file.path(dir_data, "female_population_ma.csv"),
                    col_types = cols(year = "i", pop = "d", .default = "c"))
muns <- read_csv(file.path(dir_data, "municipalities_ma.csv"),
                 col_types = cols(total_pop_mean_2013_2024 = "d", .default = "c"))
seats <- read_csv(file.path(dir_data, "municipal_seats_brazil.csv"),
                  col_types = cols(seat_lon = "d", seat_lat = "d", .default = "c"))
num <- function(x) suppressWarnings(as.numeric(x))

# --- 1. Cases (PO-SUS) ---------------------------------------------------------------
cases <- po_raw |>
  filter(sex == "F", diagnosis_category %in% c("1", "01")) |>   # malignant neoplasms
  transmute(
    year = num(year_diagnosis),
    age = na_if(num(age), 999),
    mun_residence, mun_treatment = na_if(mun_treatment, ""),
    uf_treatment = na_if(uf_treatment, ""),
    stage_code, treatment_code,
    time_raw = time_to_treatment_raw
  ) |>
  mutate(
    # Days from diagnosis to first treatment; negative intervals (treatment
    # recorded before the pathology report, e.g. diagnostic surgery) set to 0
    days = if_else(time_raw %in% c("99999", "") | is.na(time_raw), NA_real_,
                   num(str_remove_all(time_raw, "[^0-9-]"))),
    negative = !is.na(days) & (str_detect(time_raw, "^-") | days < 0),
    days = if_else(negative, 0, abs(days)),
    treated = treatment_code %in% c("1", "2", "3", "4"),
    treatment = factor(treatment_code, levels = as.character(1:5),
                       labels = c("Surgery", "Chemotherapy", "Radiotherapy",
                                  "Chemoradiotherapy", "None registered")),
    # In PO-SUS, stage is captured from chemotherapy/radiotherapy authorisations:
    # surgery-first cases are coded "not applicable" (5) and untreated cases
    # "unknown" (9), so missing stage is structural.
    stage_group = case_when(
      stage_code %in% c("1", "2") ~ "I–II",
      stage_code %in% c("3", "4") ~ "III–IV",
      stage_code == "0" ~ "Inconsistent (0)",
      stage_code == "5" ~ "Not recorded (surgery)",
      TRUE ~ "Not recorded (no treatment)") |>
      factor(levels = c("I–II", "III–IV", "Not recorded (surgery)", "Inconsistent (0)",
                        "Not recorded (no treatment)")),
    stage_detail = factor(stage_code, levels = c("0", "1", "2", "3", "4", "5", "9"),
                          labels = c("0", "I", "II", "III", "IV", "Not applicable", "Unknown")),
    late_stage = case_when(stage_group == "III–IV" ~ 1L, stage_group == "I–II" ~ 0L,
                           TRUE ~ NA_integer_),
    # Sensitivity: primary surgery assumed stage I–II (indicated for FIGO IA–IIA)
    late_stage_s1 = case_when(!is.na(late_stage) ~ late_stage,
                              treatment_code == "1" ~ 0L, TRUE ~ NA_integer_),
    delay60 = if_else(treated, as.integer(days > 60), NA_integer_),
    no_treatment = as.integer(!treated),
    age_grp = age_group(age),
    age5 = age5_group(age),
    outside_mun = if_else(treated, mun_treatment != mun_residence, NA),
    outside_state = if_else(treated, uf_treatment != "21", NA)
  )

# --- 2. Oncology hubs and distances -----------------------------------------------------
# Hub = municipality (any state) that delivered >= 1% of first treatments of
# cervical cancer in residents of Maranhão, 2013–2024
hubs <- cases |>
  filter(treated, year <= 2024) |>
  count(mun_treatment) |>
  mutate(share = n / sum(n)) |>
  filter(share >= 0.01) |>
  left_join(seats, by = c("mun_treatment" = "mun_code"))
message("Hubs: ", paste(hubs$mun_name, collapse = ", "))

seats_ma <- seats |> filter(substr(mun_code, 1, 2) == "21")
dist_to <- function(lon, lat, h) sapply(seq_along(lon), \(i) min(haversine(lon[i], lat[i], h$seat_lon, h$seat_lat)))
municipalities <- muns |>
  left_join(seats_ma |> select(mun_code, mun_name, seat_lon, seat_lat), by = "mun_code") |>
  mutate(
    dist_hub_km = dist_to(seat_lon, seat_lat, hubs),
    # sensitivity: in-state hubs only (São Luís, Imperatriz)
    dist_hub_instate_km = dist_to(seat_lon, seat_lat, hubs |> filter(state == "MA")),
    nearest_hub = hubs$mun_name[sapply(seq_along(seat_lon), \(i)
      which.min(haversine(seat_lon[i], seat_lat[i], hubs$seat_lon, hubs$seat_lat)))],
    size = cut(total_pop_mean_2013_2024, c(0, 20e3, 50e3, 100e3, Inf), right = FALSE,
               labels = c("<20,000", "20,000–49,999", "50,000–99,999", "≥100,000"))
  )

cases <- cases |>
  left_join(municipalities |> select(mun_code, health_region_code, health_region,
                                     dist_hub_km, dist_hub_instate_km, size,
                                     lon_r = seat_lon, lat_r = seat_lat),
            by = c("mun_residence" = "mun_code")) |>
  left_join(seats |> select(mun_code, lon_t = seat_lon, lat_t = seat_lat),
            by = c("mun_treatment" = "mun_code")) |>
  mutate(dist_travelled_km = if_else(treated, haversine(lon_r, lat_r, lon_t, lat_t), NA_real_)) |>
  select(-lon_r, -lat_r, -lon_t, -lat_t)

# --- 3. Deaths (SIM) ------------------------------------------------------------------
deaths <- deaths_raw |>
  mutate(age5 = age5_group(age_years), age_grp = age_group(age_years)) |>
  left_join(municipalities |> select(mun_code, health_region_code, health_region),
            by = c("mun_residence" = "mun_code"))

# --- 4. Population ----------------------------------------------------------------------
population <- pop_raw |>
  mutate(age5 = factor(age5, levels = who_std$age5)) |>
  left_join(municipalities |> select(mun_code, health_region_code), by = "mun_code")

saveRDS(cases, file.path(dir_derived, "cases.rds"))
saveRDS(deaths, file.path(dir_derived, "deaths.rds"))
saveRDS(population, file.path(dir_derived, "population.rds"))
saveRDS(municipalities, file.path(dir_derived, "municipalities.rds"))
saveRDS(hubs, file.path(dir_derived, "hubs.rds"))
message("Analytic datasets saved to output/derived/.")
