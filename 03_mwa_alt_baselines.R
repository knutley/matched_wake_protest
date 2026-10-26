# =============================================================================
# 03_mwa_alt_baselines.R
# Left- and right-coded protests vs broader baselines of "typical" protest
#
# Project: Paper 2 -- Does Deployment Bias Become Institutionalised?
# Author:  Katelyn Nutley
#
# Why: centre-coded protest (02) is a narrow, low-attention category. Here each
# side is compared with two broader baselines:
#   Unknown    : events the Paper 1 classifier could not place (~45,000)
#   All others : everything that is neither left-only nor right-only
#                (centre + unknown + compound/counter-protest events, ~59,000)
# "All others" is the preferred baseline for typical protest: it is the largest,
# best matched, and does not depend on how the centre category was built.
# Unknown is a mixture that includes unclassified left- and right-coded events,
# which pulls estimates towards zero.
#
# Main radius only. 2 baselines x 2 sides x 3 windows x 2 outcomes = 24 calls.
# Centre results from 02 are added to the tables if output/centre_results.csv
# exists.
#
# Outputs: output/alt_baseline_results.csv, tab_2_baselines.tex,
#          tab_A6_baselines_volume.tex, fig_2_baselines.pdf
# =============================================================================

source("00_mwa_functions.R")

mwa_data <- load_mwa_data()
print(table(unknown_baseline = mwa_data$side_unknown))
print(table(others_baseline  = mwa_data$side_others))

baselines <- c(side_unknown = "Unknown", side_others = "All others")

alt <- bind_rows(lapply(names(baselines), function(g) {
  bind_rows(
    run_grid(mwa_data, g, "left",  "base", comparison = "Left vs baseline",  radii = MAIN_KM)$table,
    run_grid(mwa_data, g, "right", "base", comparison = "Right vs baseline", radii = MAIN_KM)$table
  ) %>% mutate(baseline = baselines[[g]])
}))

centre_file <- file.path(OUT_DIR, "centre_results.csv")
if (file.exists(centre_file)) {
  centre <- read.csv(centre_file) %>%
    filter(spat_window == MAIN_KM) %>%
    mutate(comparison = sub(" vs centre", " vs baseline", comparison))
  alt <- bind_rows(centre[, intersect(names(alt), names(centre))], alt)
} else {
  message("Run 02_mwa_centre_comparison.R first to include the centre baseline.")
}
alt$baseline <- factor(alt$baseline, levels = c("Centre", "Unknown", "All others"))
write.csv(alt, file.path(OUT_DIR, "alt_baseline_results.csv"), row.names = FALSE)

# --- Tables: one row per baseline x side, one column per window ------------------
wide_table <- function(tab) {
  tab <- tab[order(tab$comparison, tab$baseline, tab$t_window), ]
  cells <- fmt_pct_cell(tab$pct, tab$pct_lo, tab$pct_hi, tab$pvalue)
  key <- paste(tab$comparison, tab$baseline, sep = "|")
  rows <- unique(key)
  out <- t(sapply(rows, function(k) {
    tk <- tab[key == k, ]
    c(sub(" vs baseline", "", tk$comparison[1]), as.character(tk$baseline[1]),
      cells[key == k][match(WINDOWS, tk$t_window)],
      fmt_n(min(tk$control_post)))
  }))
  out[duplicated(out[, 1]), 1] <- ""      # print each comparison label once
  out
}

pol <- subset(alt, outcome == "Police-attended protests")
write_tex_table(
  wide_table(pol),
  file = file.path(OUT_DIR, "tab_2_baselines.tex"),
  caption = sprintf("Police-attended protests after left- and right-coded protests, by baseline (%d km)", MAIN_KM),
  label = "tab:baselines",
  header = paste("Side & Baseline &", paste(WINDOWS, "days", collapse = " & "), "& Baseline $N$"),
  align = "llcccr",
  notes = paste(TABLE_NOTE_PCT,
                "Baseline $N$: smallest number of matched baseline events across windows.",
                "Unknown: events the classifier could not place. All others: centre, unknown and",
                "compound (counter-protest) events."))

vol <- subset(alt, outcome == "All protests")
write_tex_table(
  wide_table(vol),
  file = file.path(OUT_DIR, "tab_A6_baselines_volume.tex"),
  caption = sprintf("Protest volume after left- and right-coded protests, by baseline (%d km)", MAIN_KM),
  label = "tab:baselines_volume",
  header = paste("Side & Baseline &", paste(WINDOWS, "days", collapse = " & "), "& Baseline $N$"),
  align = "llcccr",
  notes = TABLE_NOTE_PCT)

# --- Figure 2 ------------------------------------------------------------------
fig2 <- plot_pct(pol, "baseline", facet_var = "comparison",
                 ylab = "Police-attended protests:\ndifference from baseline (%)",
                 caption = sprintf("Matched wake analysis, %d km. Bars: 95%% intervals.", MAIN_KM),
                 colours = c("Centre" = "#b2182b", "Unknown" = "#2166ac", "All others" = "grey40"))
ggsave(file.path(OUT_DIR, "fig_2_baselines.pdf"), fig2, width = 9, height = 4.5)

sessionInfo()
