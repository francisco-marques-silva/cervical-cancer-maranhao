# =============================================================================
# run_all.R — reproduces all analyses from the files in data/
# Run from the repository root:  Rscript run_all.R
# To rebuild data/ from the original DATASUS/IBGE sources first, run
#   Rscript R/01_download_raw_data.R
# =============================================================================

scripts <- c("R/02_prepare_data.R", "R/03_article1_distance.R", "R/04_article2_regions.R")
for (s in scripts) {
  message("\n>>> ", s)
  local(source(s, local = TRUE, encoding = "UTF-8"))
}

writeLines(capture.output(sessionInfo()), file.path("output", "session_info.txt"))
message("\nDone. Results are in output/.")
