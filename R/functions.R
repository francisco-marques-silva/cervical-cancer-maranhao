# =============================================================================
# functions.R — helper functions
# =============================================================================

#' Directly age-standardised rates (WHO World Standard), with gamma 95% CI
#' (Fay & Feuer, 1997).
#' @param cases data frame with grouping columns, `age5` and `n`
#' @param pop   data frame with grouping columns, `age5` and `pop`
#' @param ...   grouping variables (unquoted)
asr <- function(cases, pop, ..., per = 1e5) {
  grp <- enquos(...)
  w <- who_std |> mutate(w = weight / sum(weight)) |> select(age5, w)
  pop |>
    group_by(!!!grp, age5) |>
    summarise(pop = sum(pop), .groups = "drop") |>
    left_join(cases |> group_by(!!!grp, age5) |> summarise(n = sum(n), .groups = "drop"),
              by = c(map_chr(grp, as_label), "age5")) |>
    mutate(n = replace_na(n, 0)) |>
    left_join(w, by = "age5") |>
    group_by(!!!grp) |>
    summarise(
      # age-specific terms are computed before `pop` is summarised
      rate = sum(w * n / pop),
      var = sum(w^2 * n / pop^2),
      wm = max(w / pop),
      events = sum(n),
      pop = sum(pop),
      crude = events / pop * per,
      .groups = "drop"
    ) |>
    mutate(
      lcl = if_else(events == 0, 0, qgamma(0.025, shape = rate^2 / var, scale = var / rate)),
      ucl = qgamma(0.975, shape = (rate + wm)^2 / (var + wm^2), scale = (var + wm^2) / (rate + wm)),
      across(c(rate, lcl, ucl), ~ .x * per)
    ) |>
    select(!!!grp, events, pop, crude, rate, lcl, ucl)
}

#' Wilson score 95% CI for a proportion
wilson <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  tibble(p = p, lcl = centre - half, ucl = centre + half)
}

#' Great-circle (haversine) distance in km
haversine <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  a <- sin((lat2 - lat1) * rad / 2)^2 +
    cos(lat1 * rad) * cos(lat2 * rad) * sin((lon2 - lon1) * rad / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}

#' Poisson–gamma empirical Bayes smoothing of standardised ratios
eb_smooth <- function(O, E) {
  r <- O / E; m <- sum(O) / sum(E)
  s2 <- max(sum(E * (r - m)^2) / sum(E) - m / mean(E), 1e-6)
  a <- m^2 / s2; b <- m / s2
  (O + a) / (E + b)
}

fmt <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = ",")
fmt_ci <- function(e, l, h, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), e, l, h)
