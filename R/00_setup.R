# =============================================================================
# 00_setup.R
# Packages, paths, graphical theme and constants shared by all scripts.
# All scripts must be run from the repository root (see run_all.R).
# =============================================================================

required <- c("tidyverse", "sf", "gtsummary", "flextable", "officer", "sandwich",
              "lmtest", "MASS", "splines", "patchwork", "scales")
missing <- setdiff(required, rownames(installed.packages()))
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(gtsummary)
  library(flextable)
  library(officer)
  library(sandwich)
  library(lmtest)
  library(patchwork)
  library(scales)
})

options(dplyr.summarise.inform = FALSE, timeout = 1200)
sf::sf_use_s2(FALSE)
theme_gtsummary_language("en", big.mark = ",")
set_flextable_defaults(font.family = "Times New Roman", font.size = 9, padding = 2)

# --- Paths ---------------------------------------------------------------------
dir_data    <- "data"                    # analysis inputs (versioned)
dir_raw     <- Sys.getenv("RAW_DIR", "data-raw")   # raw DATASUS files (not versioned)
dir_derived <- file.path("output", "derived")      # analytic datasets
dir_tab     <- file.path("output", "tables")
dir_fig     <- file.path("output", "figures")
for (d in c(dir_derived, dir_tab, dir_fig)) dir.create(d, FALSE, TRUE)

# --- Colour palette (colour-blind safe; fixed categorical order) ----------------
pal_cat     <- c("#2a78d6", "#eb6834", "#1baf7a")
col_main    <- "#2a78d6"
col_text2   <- "#52514e"
col_ref     <- "#c9c8c3"
pal_seq     <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
pal_ord4    <- c("#cde2fb", "#86b6ef", "#3987e5", "#184f95")
# Diverging: red (deficit) <-> neutral grey (about expected) <-> blue (excess)
pal_div     <- c("#b2182b", "#ef8a62", "#e6e5e1", "#86b6ef", "#1c5cab")

theme_paper <- function(base = 10) {
  theme_minimal(base_size = base, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = base + 1),
      plot.title.position = "plot",
      legend.position = "top",
      legend.justification = "left",
      strip.text = element_text(face = "bold", hjust = 0)
    )
}

save_fig <- function(p, name, w = 16, h = 10) {
  ggsave(file.path(dir_fig, paste0(name, ".png")), p, width = w, height = h,
         units = "cm", dpi = 300, bg = "white")
  ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = w, height = h,
         units = "cm", device = cairo_pdf)
  ggsave(file.path(dir_fig, paste0(name, ".tiff")), p, width = w, height = h,
         units = "cm", dpi = 600, bg = "white", device = "tiff", compression = "lzw")
}

save_table <- function(ft, name, landscape = FALSE) {
  orient <- if (landscape) "landscape" else "portrait"
  save_as_docx(ft, path = file.path(dir_tab, paste0(name, ".docx")),
               pr_section = prop_section(page_size(orient = orient)))
  saveRDS(ft, file.path(dir_tab, paste0(name, ".rds")))
}

# --- WHO World Standard Population (Ahmad et al., 2001), 17 five-year groups ---
who_std <- tibble(
  age5 = c("00-04", "05-09", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
           "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80+"),
  weight = c(8860, 8690, 8600, 8470, 8220, 7930, 7610, 7150, 6590, 6040, 5370,
             4550, 3720, 2960, 2210, 1520, 1545)
)
age5_group <- function(age) {
  cut(age, c(seq(0, 80, 5), Inf), right = FALSE, labels = who_std$age5)
}
age_group <- function(age) {
  cut(age, c(0, 30, 40, 50, 60, 70, Inf), right = FALSE,
      labels = c("<30", "30–39", "40–49", "50–59", "60–69", "≥70"))
}

# Hub codes (IBGE 6-digit municipality codes)
hub_codes <- c(sao_luis = "211130", imperatriz = "210530", teresina = "221100")
