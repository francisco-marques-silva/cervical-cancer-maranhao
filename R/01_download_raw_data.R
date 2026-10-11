# =============================================================================
# 01_download_raw_data.R  (optional)
# Rebuilds the input files in data/ from the original public sources:
#   - Oncology Panel (PO-SUS) microdata, DATASUS FTP, files by year of diagnosis
#   - Mortality Information System (SIM-DO), via the microdatasus package
#   - Female population by municipality, year and age group (brpop package)
#   - Municipal boundaries, municipal seats and health regions (IBGE, geobr)
# The versioned files in data/ are enough to reproduce all analyses; run this
# script only to rebuild them (requires internet; ~250 MB of downloads).
# Extra packages: microdatasus, brpop, geobr.
# =============================================================================

source(file.path("R", "00_setup.R"))
for (p in c("microdatasus", "brpop", "geobr")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Package required: ", p)
}
library(microdatasus)   # process functions need the package attached
dir.create(file.path(dir_raw, "po"), FALSE, TRUE)

# --- 1. Oncology Panel (PO-SUS) ----------------------------------------------------
# Cervical cancer (ICD-10 C53) in residents of Maranhão (state code 21).
# Patient identifiers (encrypted health card number), birth date, exact dates
# and facility codes are NOT exported.
ftp_po <- "ftp://ftp.datasus.gov.br/dissemin/publicos/PAINEL_ONCOLOGIA/DADOS/"
po <- map_dfr(2013:2025, function(y) {
  f <- file.path(dir_raw, "po", sprintf("POBR%d.dbc", y))
  if (!file.exists(f)) download.file(paste0(ftp_po, basename(f)), f, mode = "wb", quiet = TRUE)
  d <- read_dbc(f)
  names(d) <- toupper(names(d))
  d |>
    mutate(across(everything(), as.character)) |>
    filter(substr(DIAG_DETH, 1, 3) == "C53", UF_RESID == "21") |>
    transmute(
      year_diagnosis = ANO_DIAGN,
      mun_residence = MUN_RESID,
      mun_diagnosis = MUN_DIAG,
      uf_treatment = UF_TRATAM,
      mun_treatment = MUN_TRATAM,
      icd10 = DIAG_DETH,
      diagnosis_category = DIAGNOSTIC,
      sex = SEXO,
      age = IDADE,
      stage_code = ESTADIAM,
      treatment_code = TRATAMENTO,
      time_to_treatment_raw = TEMPO_TRAT
    )
})
write_csv(po, file.path(dir_data, "oncology_panel_c53_ma.csv"), na = "")
message("PO-SUS: ", nrow(po), " records")

# --- 2. Mortality Information System (SIM-DO) ------------------------------------
sim_raw <- fetch_datasus(year_start = 2013, year_end = 2024, uf = "MA",
                         information_system = "SIM-DO")
lab <- function(x, codes, labels) {
  factor(x, levels = codes, labels = labels) |> as.character() |> replace_na("Unknown")
}
deaths <- sim_raw |>
  mutate(cause3 = substr(CAUSABAS, 1, 3)) |>
  filter(SEXO %in% c("2", "F"), cause3 %in% c("C53", "C54", "C55")) |>
  transmute(
    year = as.integer(substr(DTOBITO, 5, 8)),
    # IDADE: first digit = unit (4 = years, 5 = 100 years or more)
    age_years = case_when(substr(IDADE, 1, 1) == "4" ~ as.numeric(substr(IDADE, 2, 3)),
                          substr(IDADE, 1, 1) == "5" ~ 100 + as.numeric(substr(IDADE, 2, 3)),
                          TRUE ~ NA_real_),
    underlying_cause = cause3,
    mun_residence = substr(CODMUNRES, 1, 6),
    race = lab(RACACOR, 1:5, c("White", "Black", "Asian", "Brown (parda)", "Indigenous")),
    education = lab(ESC, 1:5, c("None", "1–3 years", "4–7 years", "8–11 years", "12 years or more")),
    marital_status = lab(ESTCIV, 1:5, c("Single", "Married", "Widowed",
                                        "Legally separated/divorced", "Consensual union")),
    place_of_death = lab(LOCOCOR, 1:6, c("Hospital", "Other health facility", "Home",
                                         "Public road", "Other", "Indigenous village"))
  ) |>
  filter(substr(mun_residence, 1, 2) == "21")
write_csv(deaths, file.path(dir_data, "mortality_c53_c55_ma.csv"), na = "")
message("SIM: ", nrow(deaths), " deaths (C53–C55)")

# --- 3. Female population ------------------------------------------------------------
pop <- brpop::datasus2024_mun_female_pop() |>
  filter(substr(code_muni, 1, 2) == "21", year %in% 2013:2024, age_group != "Total") |>
  mutate(mun_code = substr(as.character(code_muni), 1, 6),
         age5 = age5_group(as.numeric(str_extract(age_group, "\\d+")))) |>
  group_by(mun_code, year, age5) |>
  summarise(pop = sum(pop), .groups = "drop")
write_csv(pop, file.path(dir_data, "female_population_ma.csv"))

# --- 4. Municipalities, health regions, seats and population size ---------------------
regions <- geobr::read_health_region(year = 2013, showProgress = FALSE) |>
  filter(abbrev_state == "MA")
seats <- geobr::read_municipal_seat(year = 2010, showProgress = FALSE)
xy <- st_coordinates(seats)
seats_df <- tibble(mun_code = substr(as.character(seats$code_muni), 1, 6),
                   mun_name = seats$name_muni, state = seats$abbrev_state,
                   seat_lon = xy[, 1], seat_lat = xy[, 2]) |>
  distinct(mun_code, .keep_all = TRUE)
write_csv(seats_df, file.path(dir_data, "municipal_seats_brazil.csv"))

size <- brpop::mun_pop_totals() |>
  filter(substr(code_muni, 1, 2) == "21", year %in% 2013:2024) |>
  mutate(mun_code = substr(as.character(code_muni), 1, 6)) |>
  group_by(mun_code) |>
  summarise(total_pop_mean_2013_2024 = round(mean(pop)), .groups = "drop")

muns <- regions |>
  st_drop_geometry() |>
  transmute(mun_code = as.character(code_muni6), health_region_code = code_health_region,
            health_region = name_health_region, health_macroregion = name_health_macroregion) |>
  distinct(mun_code, .keep_all = TRUE) |>
  left_join(size, by = "mun_code")
write_csv(muns, file.path(dir_data, "municipalities_ma.csv"))

# --- 5. Maps (GeoPackage) ------------------------------------------------------------
mun_sf <- geobr::read_municipality(code_muni = "MA", year = 2022, showProgress = FALSE) |>
  transmute(mun_code = substr(as.character(code_muni), 1, 6), mun_name = name_muni)
reg_sf <- regions |>
  st_transform(5880) |> st_make_valid() |>
  group_by(health_region_code = code_health_region, health_region = name_health_region) |>
  summarise(geometry = st_union(st_buffer(geometry, 50)), .groups = "drop") |>
  st_buffer(-50) |>
  st_transform(4674)
states_sf <- geobr::read_state(year = 2020, showProgress = FALSE) |>
  transmute(state = abbrev_state)
st_write(mun_sf, file.path(dir_data, "maps", "ma_municipalities.gpkg"), delete_dsn = TRUE, quiet = TRUE)
st_write(reg_sf, file.path(dir_data, "maps", "ma_health_regions.gpkg"), delete_dsn = TRUE, quiet = TRUE)
st_write(states_sf, file.path(dir_data, "maps", "brazil_states.gpkg"), delete_dsn = TRUE, quiet = TRUE)
message("Data files written to data/.")
