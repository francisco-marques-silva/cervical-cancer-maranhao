# =============================================================================
# 04_article2_regions.R
# Article 2 — Regional inequalities in the diagnosis and treatment of cervical
# cancer in Maranhão, 2020–2024: ecological study.
#   - Age-standardised rates of registered cases (PO-SUS) and mortality (SIM),
#     by year (2013–2024 as context) and by health region (2020–2024 pooled)
#   - Timeliness indicators (late stage, treatment >60 days) by health region
#     and by period (2020–2021 pandemic vs 2022–2024 recovery)
#   - Sociodemographic profile of deaths
# The article is written in Portuguese; tables and figures here are in English.
# =============================================================================

source(file.path("R", "00_setup.R"))
source(file.path("R", "functions.R"))

cases  <- readRDS(file.path(dir_derived, "cases.rds"))
deaths <- readRDS(file.path(dir_derived, "deaths.rds"))
pop    <- readRDS(file.path(dir_derived, "population.rds"))
hubs   <- readRDS(file.path(dir_derived, "hubs.rds"))
years <- 2020:2024
cs <- cases |> filter(year %in% years)
d53 <- deaths |> filter(underlying_cause == "C53")
res <- list()

# --- 1. Annual age-standardised rates, 2013–2024 -------------------------------------
cases_y <- cases |> filter(year <= 2024, !is.na(age5)) |> count(year, age5)
deaths_y <- d53 |> filter(!is.na(age5)) |> count(year, age5)
# Sensitivity: redistribute C55 (uterus, part unspecified) deaths proportionally
# between C53 and C54, by year and age group
deaths_corr <- deaths |> filter(!is.na(age5)) |>
  count(year, age5, underlying_cause) |>
  pivot_wider(names_from = underlying_cause, values_from = n, values_fill = 0) |>
  mutate(share53 = if_else(C53 + C54 > 0, C53 / (C53 + C54), 0.5), n = C53 + C55 * share53) |>
  select(year, age5, n)
pop_y <- pop |> filter(year <= 2024)
rates_year <- bind_rows(
  asr(cases_y, pop_y, year) |> mutate(measure = "Registered cases (PO-SUS)"),
  asr(deaths_y, pop_y, year) |> mutate(measure = "Deaths (SIM)"),
  asr(deaths_corr, pop_y, year) |> mutate(measure = "Deaths, C55 redistributed")
)
write_csv(rates_year, file.path("output", "article2_rates_by_year.csv"))
res$rates_year <- rates_year
res$rates_period <- bind_rows(
  asr(cs |> mutate(all = 1) |> filter(!is.na(age5)) |> count(all, age5),
      pop |> filter(year %in% years) |> mutate(all = 1), all) |> mutate(measure = "Registered cases"),
  asr(d53 |> filter(year %in% years, !is.na(age5)) |> mutate(all = 1) |> count(all, age5),
      pop |> filter(year %in% years) |> mutate(all = 1), all) |> mutate(measure = "Deaths"))

fig1 <- ggplot(rates_year |> filter(measure != "Deaths, C55 redistributed"),
               aes(year, rate, colour = measure)) +
  annotate("rect", xmin = 2019.5, xmax = 2024.5, ymin = -Inf, ymax = Inf, fill = "grey92") +
  annotate("text", x = 2022, y = Inf, label = "Study period", vjust = 1.5, size = 3, colour = col_text2) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = measure), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  scale_colour_manual(values = pal_cat[1:2], name = NULL) +
  scale_fill_manual(values = pal_cat[1:2], guide = "none") +
  scale_x_continuous(breaks = 2013:2024) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.08))) +
  labs(x = NULL, y = "Age-standardised rate per 100,000 women") +
  theme_paper(10)
save_fig(fig1, "article2_figure1_rates_by_year", w = 16, h = 9)

# --- 2. Health regions, 2020–2024 ------------------------------------------------------
pop_r <- pop |> filter(year %in% years)
reg <- asr(cs |> filter(!is.na(age5), !is.na(health_region_code)) |> count(health_region_code, age5),
           pop_r, health_region_code) |>
  rename_with(~ paste0("cases_", .x), -health_region_code) |>
  left_join(asr(d53 |> filter(year %in% years, !is.na(age5), !is.na(health_region_code)) |>
                  count(health_region_code, age5), pop_r, health_region_code) |>
              rename_with(~ paste0("deaths_", .x), -health_region_code), by = "health_region_code") |>
  left_join(cs |> group_by(health_region_code, health_region) |>
              summarise(n_staged = sum(!is.na(late_stage)), n_late = sum(late_stage, na.rm = TRUE),
                        n_treated = sum(treated), n_delay = sum(delay60, na.rm = TRUE), .groups = "drop"),
            by = "health_region_code") |>
  arrange(desc(cases_rate))
write_csv(reg, file.path("output", "article2_health_regions.csv"))
res$regions <- reg

ft_reg <- reg |>
  transmute(`Health region` = health_region,
            `Cases (n)` = cases_events,
            `Registration rate (95% CI)` = fmt_ci(cases_rate, cases_lcl, cases_ucl, 1),
            `Deaths (n)` = deaths_events,
            `Mortality rate (95% CI)` = fmt_ci(deaths_rate, deaths_lcl, deaths_ucl, 1),
            `Stage III–IV, % (n/N)` = sprintf("%.1f (%d/%d)", 100 * n_late / n_staged, n_late, n_staged),
            `Treatment >60 days, % (n/N)` = sprintf("%.1f (%d/%d)", 100 * n_delay / n_treated, n_delay, n_treated)) |>
  flextable() |>
  add_footer_lines("Rates per 100,000 women, age-standardised to the WHO World Standard Population; 95% CI by the gamma method. Stage III–IV computed among cases with recorded stage; treatment >60 days among cases with a registered treatment. Regions ordered by registration rate.") |>
  bold(part = "header") |> fontsize(size = 8, part = "all") |> autofit()
save_table(ft_reg, "article2_table3_health_regions", landscape = TRUE)

reg_sf <- st_read(file.path(dir_data, "maps", "ma_health_regions.gpkg"), quiet = TRUE) |>
  mutate(health_region_code = as.character(health_region_code)) |>
  select(-health_region) |>
  left_join(reg, by = "health_region_code")
neighbours <- st_read(file.path(dir_data, "maps", "brazil_states.gpkg"), quiet = TRUE) |>
  filter(state %in% c("PA", "TO", "PI"))
hubs_sf <- hubs |> st_as_sf(coords = c("seat_lon", "seat_lat"), crs = 4674)
bb <- st_bbox(reg_sf)
region_map <- function(var, title) {
  ggplot() +
    geom_sf(data = neighbours, fill = "grey96", colour = "grey80", linewidth = 0.2) +
    geom_sf(data = reg_sf, aes(fill = .data[[var]]), colour = "white", linewidth = 0.3) +
    geom_sf(data = hubs_sf, shape = 21, fill = "white", colour = "black", size = 2, stroke = 0.8) +
    scale_fill_gradientn(colours = pal_seq, name = "Age-standardised rate per 100,000",
                         guide = guide_colourbar(barwidth = 7, barheight = 0.4, title.position = "top")) +
    coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], expand = TRUE) +
    labs(title = title) +
    theme_void(base_size = 9) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 9),
          legend.title = element_text(size = 8))
}
save_fig(region_map("cases_rate", "A. Registered cases (PO-SUS)") +
           region_map("deaths_rate", "B. Mortality (SIM)"),
         "article2_figure2_maps", w = 16, h = 10)

# --- 3. Table 1: cases by period ---------------------------------------------------------
t1 <- cs |>
  mutate(
    period = factor(if_else(year <= 2021, "2020–2021", "2022–2024")),
    time_grp = case_when(!treated ~ "No registered treatment", delay60 == 0 ~ "≤60 days",
                         TRUE ~ ">60 days") |>
      factor(levels = c("≤60 days", ">60 days", "No registered treatment")),
    place = case_when(!treated ~ "No registered treatment",
                      !outside_mun ~ "Municipality of residence",
                      !outside_state ~ "Other municipality in Maranhão",
                      TRUE ~ "Other state") |>
      factor(levels = c("Municipality of residence", "Other municipality in Maranhão",
                        "Other state", "No registered treatment"))
  ) |>
  select(period, age, age_grp, stage_group, treatment, time_grp, place) |>
  tbl_summary(
    by = period,
    label = list(age ~ "Age (years), median (IQR)", age_grp ~ "Age group (years)",
                 stage_group ~ "Recorded stage", treatment ~ "First treatment",
                 time_grp ~ "Time to first treatment", place ~ "Place of first treatment"),
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})", all_categorical() ~ "{n} ({p})"),
    digits = list(all_categorical() ~ c(0, 1), all_continuous() ~ 0),
    missing = "no"
  ) |>
  add_overall(last = TRUE, col_label = "**Total**  \nN = {style_number(N)}") |>
  add_p(test = list(all_categorical() ~ "chisq.test", all_continuous() ~ "wilcox.test"),
        pvalue_fun = label_style_pvalue(digits = 3)) |>
  modify_header(label ~ "**Characteristic**") |>
  modify_footnote(everything() ~ NA)
ft1 <- as_flex_table(t1) |>
  add_footer_lines("IQR, interquartile range. n (%) unless otherwise stated. p value: Pearson chi-squared (categorical) or Mann-Whitney (age). In the Oncology Panel, stage is recorded only when the first treatment is chemotherapy and/or radiotherapy.") |>
  fontsize(size = 8, part = "all") |> autofit()
save_table(ft1, "article2_table1_cases_by_period")

# --- 4. Table 2: deaths ----------------------------------------------------------------------
t2 <- d53 |> filter(year %in% years) |>
  mutate(race = factor(race, levels = c("White", "Brown (parda)", "Black", "Asian", "Indigenous", "Unknown")),
         education = factor(education, levels = c("None", "1–3 years", "4–7 years", "8–11 years",
                                                  "12 years or more", "Unknown")),
         marital_status = fct_relevel(factor(marital_status), "Unknown", after = Inf),
         place_of_death = fct_relevel(fct_lump_n(factor(place_of_death), 3, other_level = "Other"),
                                      "Other", after = Inf)) |>
  select(age_years, age_grp, race, education, marital_status, place_of_death) |>
  tbl_summary(
    label = list(age_years ~ "Age (years), median (IQR)", age_grp ~ "Age group (years)",
                 race ~ "Race/skin colour", education ~ "Schooling (years of study)",
                 marital_status ~ "Marital status", place_of_death ~ "Place of death"),
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})", all_categorical() ~ "{n} ({p})"),
    digits = list(all_categorical() ~ c(0, 1), all_continuous() ~ 0),
    missing = "no"
  ) |>
  modify_header(label ~ "**Characteristic**", stat_0 ~ "**n (%)**  \nN = {style_number(N)}") |>
  modify_footnote(everything() ~ NA)
ft2 <- as_flex_table(t2) |>
  add_footer_lines("IQR, interquartile range. Underlying cause: ICD-10 C53.") |>
  fontsize(size = 8, part = "all") |> autofit()
save_table(ft2, "article2_table2_deaths")

# --- 5. Supplementary: timeliness indicators by year, 2013–2024 -------------------------------
ind <- cases |> filter(year <= 2024) |>
  group_by(year) |>
  summarise(n_staged = sum(!is.na(late_stage)), n_late = sum(late_stage, na.rm = TRUE),
            n_treated = sum(treated), n_delay = sum(delay60, na.rm = TRUE),
            n = n(), n_untreated = sum(!treated), .groups = "drop")
write_csv(ind, file.path("output", "article2_indicators_by_year.csv"))
ind_long <- bind_rows(
  ind |> transmute(year, indicator = "Stage III–IV (among staged)", x = n_late, n = n_staged),
  ind |> transmute(year, indicator = "Treatment started after 60 days", x = n_delay, n = n_treated)
) |> mutate(wilson(x, n))
figS1 <- ggplot(ind_long, aes(year, p, colour = indicator)) +
  annotate("rect", xmin = 2019.5, xmax = 2024.5, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_linerange(aes(ymin = lcl, ymax = ucl), linewidth = 0.5, alpha = 0.6) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  scale_colour_manual(values = pal_cat[1:2], name = NULL) +
  scale_x_continuous(breaks = 2013:2024) +
  scale_y_continuous(limits = c(0, 1), labels = label_percent(), expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "Proportion (95% CI)") +
  theme_paper(10)
save_fig(figS1, "article2_figureS1_indicators_by_year", w = 16, h = 9)

saveRDS(res, file.path(dir_derived, "results_article2.rds"))
message("Article 2 analysis complete.")
