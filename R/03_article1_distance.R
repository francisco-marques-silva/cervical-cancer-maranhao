# =============================================================================
# 03_article1_distance.R
# Article 1 — Distance to oncology hubs and access to cervical cancer care in
# Maranhão, Brazil, 2013–2024.
#   Individual level (PO-SUS): late stage (FIGO III–IV), treatment >60 days and
#   no registered treatment by distance from residence to the nearest hub
#   (Poisson regression with cluster-robust variance).
#   Ecological level (217 municipalities): standardised registration and
#   mortality ratios by distance (indirect standardisation, negative binomial).
# =============================================================================

source(file.path("R", "00_setup.R"))
source(file.path("R", "functions.R"))
library(MASS, exclude = "select")
library(splines)

cases  <- readRDS(file.path(dir_derived, "cases.rds"))
deaths <- readRDS(file.path(dir_derived, "deaths.rds"))
pop    <- readRDS(file.path(dir_derived, "population.rds"))
muns   <- readRDS(file.path(dir_derived, "municipalities.rds"))
hubs   <- readRDS(file.path(dir_derived, "hubs.rds"))
res <- list()

bands <- c("Hub municipality", "1–99 km", "100–199 km", "≥200 km")
dist_band <- function(x) cut(x, c(-Inf, 1, 100, 200, Inf), labels = bands, right = FALSE)

# --- 1. Study population ----------------------------------------------------------------
res$flow <- list(records_2013_2025 = nrow(cases), records_2013_2024 = sum(cases$year <= 2024))
d <- cases |>
  filter(year <= 2024, age >= 15) |>
  mutate(
    dist_cat = dist_band(dist_hub_km),
    dist_cat_instate = dist_band(dist_hub_instate_km),
    dist100 = dist_hub_km / 100,
    period = cut(year, c(2012, 2016, 2019, 2021, 2024),
                 labels = c("2013–2016", "2017–2019", "2020–2021", "2022–2024")),
    hub_used = case_when(!treated ~ NA_character_,
                         mun_treatment == hub_codes["sao_luis"] ~ "São Luís",
                         mun_treatment == hub_codes["imperatriz"] ~ "Imperatriz",
                         mun_treatment == hub_codes["teresina"] ~ "Teresina",
                         TRUE ~ "Other")
  )
res$flow <- c(res$flow, list(included = nrow(d), staged = sum(!is.na(d$late_stage)),
                             treated = sum(d$treated),
                             stage0 = sum(d$stage_group == "Inconsistent (0)")))

# --- 2. Table 1: characteristics by distance band -------------------------------------------
t1 <- d |>
  mutate(delay = factor(delay60, 0:1, c("≤60 days", ">60 days"))) |>
  select(dist_cat, age, age_grp, period, size, treatment, stage_group, days, delay) |>
  tbl_summary(
    by = dist_cat,
    label = list(age ~ "Age, years", age_grp ~ "Age group, years", period ~ "Period of diagnosis",
                 size ~ "Municipal population", treatment ~ "First treatment",
                 stage_group ~ "Recorded FIGO stage", days ~ "Days from diagnosis to treatment",
                 delay ~ "Time to treatment"),
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})", all_categorical() ~ "{n} ({p})"),
    digits = list(all_categorical() ~ c(0, 1), all_continuous() ~ 0),
    missing = "no"
  ) |>
  add_overall() |>
  modify_header(label ~ "**Characteristic**") |>
  modify_spanning_header(all_stat_cols(FALSE) ~ "**Distance from residence to nearest oncology hub**") |>
  modify_footnote(everything() ~ NA)
ft1 <- as_flex_table(t1) |>
  add_footer_lines("Median (interquartile range) or n (%). Time to treatment computed among women with a registered first treatment. Stage is recorded in PO-SUS only when the first treatment is chemotherapy and/or radiotherapy.") |>
  fontsize(size = 8, part = "all") |> autofit()
save_table(ft1, "article1_table1", landscape = TRUE)

# --- 3. Poisson models with cluster-robust variance ----------------------------------------
complete <- function(formula, data) data[complete.cases(data[, all.vars(formula)]), ]
pr_fit <- function(formula, data, expo = "dist_cat") {
  data <- complete(formula, data)
  m <- glm(formula, family = poisson(link = "log"), data = data)
  v <- vcovCL(m, cluster = data$mun_residence, type = "HC0")
  ct <- coeftest(m, vcov. = v); ci <- coefci(m, vcov. = v)
  tibble(term = rownames(ct), est = exp(ct[, 1]), lcl = exp(ci[, 1]), ucl = exp(ci[, 2]),
         p = ct[, 4], n = nobs(m)) |>
    filter(str_detect(term, expo))
}

d_late  <- d |> filter(!is.na(late_stage))
d_delay <- d |> filter(treated)
covs <- "age_grp + period + size"
f <- function(y, x = "dist_cat", extra = "") as.formula(paste(y, "~", x, "+", covs, extra))

models <- list(
  late_crude  = list(f = late_stage ~ dist_cat, data = d_late),
  late_adj    = list(f = f("late_stage"), data = d_late),
  delay_crude = list(f = delay60 ~ dist_cat, data = d_delay),
  delay_adj   = list(f = f("delay60"), data = d_delay),
  delay_adj2  = list(f = f("delay60", extra = "+ treatment + stage_group"),
                     data = d_delay |> mutate(treatment = fct_drop(treatment),
                                              stage_group = fct_drop(stage_group))),
  notr_crude  = list(f = no_treatment ~ dist_cat, data = d),
  notr_adj    = list(f = f("no_treatment"), data = d)
)
pr <- imap_dfr(models, ~ pr_fit(.x$f, .x$data) |> mutate(model = .y))
res$pr <- pr

res$trend <- bind_rows(
  late = pr_fit(f("late_stage", "dist100"), d_late, "dist100"),
  delay = pr_fit(f("delay60", "dist100"), d_delay, "dist100"),
  no_treatment = pr_fit(f("no_treatment", "dist100"), d, "dist100"),
  .id = "outcome")

prop_band <- bind_rows(
  d_late |> group_by(dist_cat) |> summarise(x = sum(late_stage), n = n()) |> mutate(outcome = "Late stage (III–IV)"),
  d_delay |> group_by(dist_cat) |> summarise(x = sum(delay60), n = n()) |> mutate(outcome = "Treatment >60 days"),
  d |> group_by(dist_cat) |> summarise(x = sum(no_treatment), n = n()) |> mutate(outcome = "No registered treatment")
) |> mutate(wilson(x, n))
res$prop_band <- prop_band

# --- 4. Sensitivity analyses ------------------------------------------------------------------
sens <- list(
  s1_late_surgery_early = pr_fit(f("late_stage_s1"), d |> filter(!is.na(late_stage_s1))),
  s2_untimely = pr_fit(f("untimely"), d |> mutate(untimely = as.integer(no_treatment == 1 | delay60 %in% 1))),
  s3_late_instate = pr_fit(f("late_stage", "dist_cat_instate"), d_late, "dist_cat_instate"),
  s3_delay_instate = pr_fit(f("delay60", "dist_cat_instate"), d_delay, "dist_cat_instate"),
  s4_delay_travel = pr_fit(f("delay60", "travel_cat"),
                           d_delay |> mutate(travel_cat = dist_band(dist_travelled_km)), "travel_cat"),
  s5_late_25_79 = pr_fit(f("late_stage"), d_late |> filter(age >= 25, age < 80)),
  s5_delay_25_79 = pr_fit(f("delay60"), d_delay |> filter(age >= 25, age < 80))
)
res$sensitivity <- bind_rows(sens, .id = "analysis")

# Effect modification by period (2013–2019 vs 2020–2024): multivariable Wald test
effect_mod <- function(y, data) {
  data <- data |> mutate(pandemic = year >= 2020)
  fm <- as.formula(paste(y, "~ dist_cat * pandemic + age_grp + size"))
  data <- complete(fm, data)
  m <- glm(fm, family = poisson, data = data)
  v <- vcovCL(m, cluster = data$mun_residence, type = "HC0")
  ii <- grep(":", names(coef(m)))
  wald <- as.numeric(t(coef(m)[ii]) %*% solve(v[ii, ii]) %*% coef(m)[ii])
  strata <- data |> group_split(pandemic) |>
    map_dfr(~ pr_fit(as.formula(paste(y, "~ dist_cat + age_grp + size")), .x) |>
              mutate(pandemic = unique(.x$pandemic)))
  list(p_interaction = pchisq(wald, length(ii), lower.tail = FALSE), strata = strata)
}
res$effect_modification <- list(late = effect_mod("late_stage", d_late),
                                delay = effect_mod("delay60", d_delay))

res$hub_delay <- d_delay |> group_by(hub_used) |>
  summarise(n = n(), prop_delay = mean(delay60), median_days = median(days),
            q1 = quantile(days, .25), q3 = quantile(days, .75))

# Natural splines; marginally standardised predicted proportions with delta-method CI
spline_pred <- function(y, data) {
  fm <- as.formula(paste(y, "~ ns(dist_hub_km, knots = c(50, 150)) + age_grp + period + size"))
  data <- complete(fm, data)
  m <- glm(fm, family = poisson, data = data)
  v <- vcovCL(m, cluster = data$mun_residence, type = "HC0")
  map_dfr(seq(0, 300, 5), function(k) {
    X <- model.matrix(delete.response(terms(m)), data |> mutate(dist_hub_km = k))
    p <- exp(X %*% coef(m))
    g <- colMeans(as.vector(p) * X)
    se <- sqrt(as.numeric(t(g) %*% v %*% g))
    tibble(dist = k, p = mean(p), lcl = mean(p) - 1.96 * se, ucl = mean(p) + 1.96 * se)
  })
}
outcome_levels <- c("Late stage (III–IV)", "Treatment >60 days", "No registered treatment")
res$spline <- bind_rows(
  spline_pred("late_stage", d_late) |> mutate(outcome = outcome_levels[1]),
  spline_pred("delay60", d_delay) |> mutate(outcome = outcome_levels[2]),
  spline_pred("no_treatment", d) |> mutate(outcome = outcome_levels[3])
)

# --- 5. Ecological analysis ---------------------------------------------------------------------
# Indirect standardisation with state age-specific rates (2013–2024)
pop_m <- pop |> filter(year <= 2024) |> group_by(mun_code, age5) |> summarise(pop = sum(pop), .groups = "drop")
cases_m <- d |> count(mun_code = mun_residence, age5)
deaths_m <- deaths |> filter(underlying_cause == "C53", !is.na(age5)) |> count(mun_code = mun_residence, age5)
state_rates <- pop_m |> group_by(age5) |> summarise(pop = sum(pop)) |>
  left_join(cases_m |> group_by(age5) |> summarise(c = sum(n)), by = "age5") |>
  left_join(deaths_m |> group_by(age5) |> summarise(o = sum(n)), by = "age5") |>
  mutate(across(c(c, o), ~ replace_na(.x, 0)), r_c = c / pop, r_o = o / pop)
eco <- pop_m |>
  left_join(state_rates |> select(age5, r_c, r_o), by = "age5") |>
  group_by(mun_code) |>
  summarise(E_c = sum(pop * r_c), E_o = sum(pop * r_o), pop = sum(pop), .groups = "drop") |>
  left_join(cases_m |> group_by(mun_code) |> summarise(O_c = sum(n)), by = "mun_code") |>
  left_join(deaths_m |> group_by(mun_code) |> summarise(O_o = sum(n)), by = "mun_code") |>
  mutate(across(c(O_c, O_o), ~ replace_na(.x, 0))) |>
  left_join(muns |> select(mun_code, dist_hub_km, nearest_hub, size), by = "mun_code") |>
  mutate(dist_cat = dist_band(dist_hub_km), dist100 = dist_hub_km / 100)
res$eco <- eco

eco_band <- eco |> group_by(dist_cat) |>
  summarise(n_mun = n(), O_c = sum(O_c), E_c = sum(E_c), O_o = sum(O_o), E_o = sum(E_o)) |>
  mutate(SRR = O_c / E_c, SMR = O_o / E_o, deaths_per_case = O_o / O_c)
res$eco_band <- eco_band

nb <- function(y, off, x, data) {
  m <- glm.nb(as.formula(paste(y, "~", x, "+ size + offset(log(", off, "))")), data = data)
  ci <- exp(confint.default(m))
  tibble(term = names(coef(m)), est = exp(coef(m)), lcl = ci[, 1], ucl = ci[, 2],
         p = summary(m)$coefficients[, 4]) |> filter(str_detect(term, x))
}
qp <- function(x, data) {
  m <- glm(as.formula(paste("O_o ~", x, "+ size + offset(log(O_c + 0.5))")), family = quasipoisson, data = data)
  ci <- exp(confint.default(m))
  tibble(term = names(coef(m)), est = exp(coef(m)), lcl = ci[, 1], ucl = ci[, 2],
         p = summary(m)$coefficients[, 4]) |> filter(str_detect(term, x))
}
eco_nonhub <- eco |> filter(dist_hub_km > 0)
eco_no_teresina <- eco_nonhub |> filter(nearest_hub != "Teresina")
res$eco_rr <- bind_rows(
  nb("O_c", "E_c", "dist_cat", eco) |> mutate(outcome = "Registered cases"),
  nb("O_o", "E_o", "dist_cat", eco) |> mutate(outcome = "Deaths"),
  qp("dist_cat", eco) |> mutate(outcome = "Deaths per registered case"))
res$eco_per100km <- bind_rows(
  nb("O_c", "E_c", "dist100", eco_nonhub) |> mutate(outcome = "Registered cases"),
  nb("O_o", "E_o", "dist100", eco_nonhub) |> mutate(outcome = "Deaths"),
  qp("dist100", eco_nonhub) |> mutate(outcome = "Deaths per registered case"))
res$eco_per100km_no_teresina <- bind_rows(
  nb("O_c", "E_c", "dist100", eco_no_teresina) |> mutate(outcome = "Registered cases"),
  nb("O_o", "E_o", "dist100", eco_no_teresina) |> mutate(outcome = "Deaths")) |>
  mutate(n_municipalities = nrow(eco_no_teresina))

# --- 6. Tables 2, 3 and S1 ----------------------------------------------------------------------
pr_col <- function(mod) {
  r <- pr |> filter(model == mod)
  c("1.00 (ref.)", fmt_ci(r$est, r$lcl, r$ucl))
}
block <- function(outcome, crude, adj, adj2 = NULL) {
  pb <- prop_band |> filter(outcome == !!outcome)
  tibble(Outcome = outcome, `Distance to nearest hub` = bands,
         `n/N (%)` = sprintf("%s/%s (%.1f)", format(pb$x, big.mark = ","), format(pb$n, big.mark = ","), 100 * pb$p),
         `Crude PR (95% CI)` = pr_col(crude), `Adjusted PR (95% CI)` = pr_col(adj),
         `Further adjusted PR (95% CI)` = if (is.null(adj2)) "—" else pr_col(adj2))
}
t2 <- bind_rows(block(outcome_levels[1], "late_crude", "late_adj"),
                block(outcome_levels[2], "delay_crude", "delay_adj", "delay_adj2"),
                block(outcome_levels[3], "notr_crude", "notr_adj"))
ft2 <- flextable(t2) |> merge_v(j = "Outcome") |> valign(j = 1, valign = "top") |>
  add_footer_lines("PR, prevalence ratio from Poisson regression with cluster-robust (municipality of residence) variance. Adjusted for age group, period of diagnosis and municipal population size. The further adjusted model also includes first treatment modality and recorded stage. Late stage analysed among women whose first treatment was chemotherapy and/or radiotherapy (stage recorded); time to treatment among women with a registered treatment; no registered treatment among all women.") |>
  bold(part = "header") |> fontsize(size = 8, part = "all") |> autofit() |>
  hline(i = c(4, 8), border = fp_border(color = "grey60"))
save_table(ft2, "article1_table2", landscape = TRUE)

rr_col <- function(out) {
  r <- res$eco_rr |> filter(outcome == out)
  c("1.00 (ref.)", fmt_ci(r$est, r$lcl, r$ucl))
}
t3 <- eco_band |>
  transmute(`Distance to nearest hub` = dist_cat, Municipalities = n_mun,
            `Registered cases, O/E (SRR)` = sprintf("%s/%.0f (%.2f)", format(O_c, big.mark = ","), E_c, SRR),
            `Adjusted RR, registered cases` = rr_col("Registered cases"),
            `Deaths, O/E (SMR)` = sprintf("%s/%.0f (%.2f)", format(O_o, big.mark = ","), E_o, SMR),
            `Adjusted RR, deaths` = rr_col("Deaths"),
            `Deaths per registered case` = sprintf("%.2f", deaths_per_case),
            `Adjusted ratio, deaths/registered` = rr_col("Deaths per registered case"))
ft3 <- flextable(t3) |>
  add_footer_lines("O, observed; E, expected from state age-specific rates (indirect standardisation, 5-year age groups, 2013–2024); SRR, standardised registration ratio; SMR, standardised mortality ratio. Rate ratios (RR) from negative binomial regression with log(E) offset, adjusted for municipal population size. Deaths per registered case modelled with quasi-Poisson regression and log(registered cases) offset. Deaths: ICD-10 C53 as underlying cause (SIM).") |>
  bold(part = "header") |> fontsize(size = 8, part = "all") |> autofit()
save_table(ft3, "article1_table3", landscape = TRUE)

labels_s <- c(s1_late_surgery_early = "S1. Late stage; primary surgery assumed stage I–II",
              s2_untimely = "S2. No treatment or treatment >60 days (composite)",
              s3_late_instate = "S3. Late stage; distance to in-state hubs only",
              s3_delay_instate = "S3. Treatment >60 days; distance to in-state hubs only",
              s4_delay_travel = "S4. Treatment >60 days; distance actually travelled",
              s5_late_25_79 = "S5. Late stage; women aged 25–79 years",
              s5_delay_25_79 = "S5. Treatment >60 days; women aged 25–79 years")
ts1 <- res$sensitivity |>
  mutate(Analysis = labels_s[analysis], band = str_remove(term, "dist_cat_instate|dist_cat|travel_cat"),
         est = fmt_ci(est, lcl, ucl)) |>
  select(Analysis, N = n, band, est) |>
  pivot_wider(names_from = band, values_from = est) |>
  mutate(`Hub municipality` = "1.00 (ref.)", .after = N)
ftS1 <- flextable(ts1) |>
  add_footer_lines("Adjusted prevalence ratios (95% CI); same covariates as the main models. In S4, the reference category is treatment in the municipality of residence.") |>
  bold(part = "header") |> fontsize(size = 8, part = "all") |> autofit()
save_table(ftS1, "article1_tableS1", landscape = TRUE)

# --- 7. Figures ---------------------------------------------------------------------------------
mun_sf <- st_read(file.path(dir_data, "maps", "ma_municipalities.gpkg"), quiet = TRUE) |>
  left_join(eco, by = "mun_code") |>
  mutate(srr_eb = eb_smooth(O_c, E_c),
         srr_cat = cut(srr_eb, c(0, 0.6, 0.8, 1.25, 1.6, Inf), right = FALSE,
                       labels = c("<0.60", "0.60–0.79", "0.80–1.24", "1.25–1.59", "≥1.60")))
neighbours <- st_read(file.path(dir_data, "maps", "brazil_states.gpkg"), quiet = TRUE) |>
  filter(state %in% c("PA", "TO", "PI"))
hubs_sf <- hubs |> st_as_sf(coords = c("seat_lon", "seat_lat"), crs = 4674) |>
  arrange(match(mun_treatment, hub_codes[c("imperatriz", "sao_luis", "teresina")]))
bb <- st_bbox(mun_sf)
base_map <- list(
  geom_sf(data = neighbours, fill = "grey96", colour = "grey80", linewidth = 0.2),
  theme_void(base_size = 9),
  theme(legend.position = "bottom", legend.title = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 9), legend.key.size = unit(0.35, "cm"))
)
m1 <- ggplot() + base_map +
  geom_sf(data = mun_sf, aes(fill = dist_cat), colour = "white", linewidth = 0.1) +
  geom_sf(data = hubs_sf, shape = 21, fill = "white", colour = "black", size = 2.2, stroke = 0.9) +
  geom_sf_label(data = hubs_sf, aes(label = mun_name), size = 2.5, nudge_x = c(0, 0.75, 0),
                nudge_y = c(0.42, 0.3, -0.42), fontface = "bold", label.size = 0, fill = alpha("white", 0.8)) +
  scale_fill_manual(values = pal_ord4, name = "Distance to nearest hub", drop = FALSE) +
  guides(fill = guide_legend(title.position = "top", nrow = 2)) +
  labs(title = "A. Distance to nearest oncology hub") +
  coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)])
m2 <- ggplot() + base_map +
  geom_sf(data = mun_sf, aes(fill = srr_cat), colour = "white", linewidth = 0.1) +
  geom_sf(data = hubs_sf, shape = 21, fill = "white", colour = "black", size = 2.2, stroke = 0.9) +
  scale_fill_manual(values = pal_div, name = "Standardised registration ratio (smoothed)", drop = FALSE) +
  guides(fill = guide_legend(title.position = "top", nrow = 2)) +
  labs(title = "B. Registered cases, observed/expected") +
  coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)])
save_fig(m1 + m2, "article1_figure1_maps", w = 17, h = 12)

d_f2 <- pr |> filter(model %in% c("late_adj", "delay_adj", "notr_adj")) |>
  mutate(band = factor(str_remove(term, "dist_cat"), levels = bands),
         outcome = recode(model, late_adj = outcome_levels[1], delay_adj = outcome_levels[2],
                          notr_adj = outcome_levels[3])) |>
  bind_rows(tibble(band = factor(bands[1], levels = bands), est = 1, lcl = 1, ucl = 1,
                   outcome = outcome_levels)) |>
  mutate(outcome = factor(outcome, levels = outcome_levels))
f2 <- ggplot(d_f2, aes(est, fct_rev(band))) +
  geom_vline(xintercept = 1, colour = col_ref, linewidth = 0.5) +
  geom_linerange(aes(xmin = lcl, xmax = ucl), colour = col_main, linewidth = 0.6) +
  geom_point(shape = 21, fill = "white", colour = col_main, size = 2.4, stroke = 1) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_x_log10(breaks = c(0.6, 0.8, 1, 1.25, 1.5), labels = label_number(accuracy = 0.01)) +
  labs(x = "Adjusted prevalence ratio (95% CI, log scale)", y = NULL) +
  theme_paper(9) + theme(panel.grid.major.y = element_blank(),
                         panel.grid.major.x = element_line(colour = "grey92"))
save_fig(f2, "article1_figure2_forest", w = 17, h = 7)

pois_ci <- function(O, E) tibble(lcl = qchisq(0.025, 2 * O) / 2 / E, ucl = qchisq(0.975, 2 * (O + 1)) / 2 / E)
d_f3 <- bind_rows(
  eco_band |> transmute(dist_cat, O = O_c, E = E_c, outcome = "Registered cases (SRR)"),
  eco_band |> transmute(dist_cat, O = O_o, E = E_o, outcome = "Deaths (SMR)")
) |> mutate(est = O / E, pois_ci(O, E),
            outcome = factor(outcome, levels = c("Registered cases (SRR)", "Deaths (SMR)")))
f3a <- ggplot(d_f3, aes(dist_cat, est, colour = outcome, group = outcome)) +
  geom_hline(yintercept = 1, colour = col_ref, linewidth = 0.5) +
  geom_line(linewidth = 0.6, position = position_dodge(0.3)) +
  geom_linerange(aes(ymin = lcl, ymax = ucl), linewidth = 0.6, position = position_dodge(0.3)) +
  geom_point(shape = 21, fill = "white", size = 2.4, stroke = 1, position = position_dodge(0.3)) +
  scale_colour_manual(values = pal_cat[1:2], name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.05))) +
  labs(x = "Distance from municipality to nearest hub", y = "Observed/expected (95% CI)",
       title = "A. Standardised ratios") +
  theme_paper(9)
d_f3b <- eco_band |> mutate(est = deaths_per_case,
                            lcl = est * exp(-1.96 * sqrt(1 / O_o + 1 / O_c)),
                            ucl = est * exp(1.96 * sqrt(1 / O_o + 1 / O_c)))
f3b <- ggplot(d_f3b, aes(dist_cat, est, group = 1)) +
  geom_line(colour = col_main, linewidth = 0.6) +
  geom_linerange(aes(ymin = lcl, ymax = ucl), colour = col_main, linewidth = 0.6) +
  geom_point(shape = 21, fill = "white", colour = col_main, size = 2.4, stroke = 1) +
  geom_text(aes(label = sprintf("%.2f", est)), nudge_x = 0.25, size = 2.8, colour = col_text2) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.05))) +
  labs(x = "Distance from municipality to nearest hub", y = "Deaths per registered case (95% CI)",
       title = "B. Deaths per registered case") +
  theme_paper(9)
save_fig(f3a + f3b, "article1_figure3_ecological", w = 17, h = 8)

fS1 <- ggplot(res$spline |> mutate(outcome = factor(outcome, levels = outcome_levels)), aes(dist, p)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = col_main, alpha = 0.15) +
  geom_line(colour = col_main, linewidth = 0.7) +
  geom_rug(data = d |> distinct(mun_residence, dist_hub_km) |> filter(dist_hub_km <= 300),
           aes(x = dist_hub_km), inherit.aes = FALSE, colour = "grey60", alpha = 0.5,
           length = unit(0.02, "npc")) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), labels = label_percent()) +
  scale_x_continuous(breaks = seq(0, 300, 100)) +
  labs(x = "Distance from residence to nearest oncology hub (km)", y = "Standardised predicted proportion") +
  theme_paper(9)
save_fig(fS1, "article1_figureS1_spline", w = 17, h = 7)

# Graphical abstract (531 x 1328 px or proportional; readable at 5 x 13 cm)
ga_data <- eco_band |> select(dist_cat, SRR, SMR) |>
  pivot_longer(-dist_cat, names_to = "measure", values_to = "ratio") |>
  mutate(measure = factor(recode(measure, SRR = "Registered", SMR = "Deaths"),
                          levels = c("Registered", "Deaths")), x = as.integer(dist_cat))
ga_plot <- ggplot(ga_data, aes(x, ratio, colour = measure)) +
  geom_hline(yintercept = 1, colour = col_ref, linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(shape = 21, fill = "white", size = 2, stroke = 1) +
  geom_text(data = ga_data |> filter(x == 4), aes(label = measure), hjust = 0, nudge_x = 0.15,
            size = 2.6, fontface = "bold") +
  scale_colour_manual(values = pal_cat[1:2], guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = c("Hub", "1–99", "100–199", "≥200 km"),
                     expand = expansion(add = c(0.2, 1.3))) +
  scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2)) +
  labs(x = "Distance to nearest oncology hub", y = "Observed/expected") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 7))
per100 <- res$eco_per100km |> filter(outcome == "Registered cases")
ga_text <- ggplot() +
  annotate("text", x = 0, y = 0.92, hjust = 0, vjust = 1, size = 2.7, fontface = "bold",
           label = "Cervical cancer, Maranhão, Brazil, 2013–2024") +
  annotate("text", x = 0, y = 0.70, hjust = 0, vjust = 1, size = 2.6, lineheight = 1.05,
           label = sprintf(paste0("Among %s registered women, distance was\n",
                                  "not associated with late stage or >60-day delay.\n\n",
                                  "Across 217 municipalities, registration fell\n",
                                  "%.0f%% per 100 km while mortality barely changed:\n",
                                  "deaths per registered case %.2f → %.2f."),
                           format(nrow(d), big.mark = ","), 100 * (1 - per100$est),
                           eco_band$deaths_per_case[1], eco_band$deaths_per_case[4])) +
  annotate("text", x = 0, y = 0.02, hjust = 0, vjust = 0, size = 2.6, fontface = "bold",
           colour = col_main, label = "Distance keeps women out of the care system.") +
  xlim(0, 1) + ylim(0, 1) + theme_void()
ga <- (ga_text | ga_plot) + plot_layout(widths = c(1.25, 1)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(dir_fig, "article1_graphical_abstract.tiff"), ga, width = 13.28, height = 5.31,
       units = "cm", dpi = 600, device = "tiff", compression = "lzw", bg = "white")
ggsave(file.path(dir_fig, "article1_graphical_abstract.png"), ga, width = 13.28, height = 5.31,
       units = "cm", dpi = 300, bg = "white")

# --- 8. Key estimates (machine-readable) -----------------------------------------------------
key <- bind_rows(
  pr |> mutate(analysis = paste0("PR_", model)),
  res$trend |> mutate(analysis = paste0("PR_per100km_", outcome)),
  res$sensitivity |> mutate(analysis = paste0("PR_", analysis)),
  res$eco_rr |> mutate(analysis = paste0("RR_band_", outcome)),
  res$eco_per100km |> mutate(analysis = paste0("RR_per100km_", outcome)),
  res$eco_per100km_no_teresina |> mutate(analysis = paste0("RR_per100km_noTeresina_", outcome))
) |> select(analysis, term, est, lcl, ucl, p, n)
write_csv(key, file.path("output", "article1_key_estimates.csv"))
saveRDS(res, file.path(dir_derived, "results_article1.rds"))
message("Article 1 analysis complete.")
