# =============================================================================
# 02_mwa_centre_comparison.R
# Left- and right-coded protests vs matched CENTRE-coded protests
#
# Project: Paper 2 -- Does Partisan Bias Become Institutionalised?
# Author:  Katelyn Nutley
#
# Why: every estimate in 01 compares right with left, so a gap could come from
# police standing down after the right OR stepping up after the left. Comparing
# each side with centre-coded protests (Paper 1's reference category) separates
# the two:
#   right vs centre < 0  -> lower attention after the right
#   left  vs centre > 0  -> higher attention after the left
#
# Note on the centre category: Paper 1 built it deliberately broad (including
# single-issue groups) so that the classification was conservative and
# minimised false-positive assignments to either pole. Centre-coded protest is
# therefore a low-attention category; 03 checks the results against broader
# baselines.
#
# Same specification as the main results. 2 sides x 3 windows x 2 outcomes
# = 12 matchedwake() calls, each covering 2, 4 and 6 km.
#
# Outputs: output/centre_results.csv, tab_A5_centre.tex, fig_A2_centre.pdf
# =============================================================================

source("00_mwa_functions.R")

mwa_data <- load_mwa_data()   # built by 01, or rebuilt here if missing
print(table(mwa_data$side_centre))   # base = centre, ~3,500 events

centre <- bind_rows(
  run_grid(mwa_data, "side_centre", "left",  "base", comparison = "Left vs centre")$table,
  run_grid(mwa_data, "side_centre", "right", "base", comparison = "Right vs centre")$table
)
centre$baseline <- "Centre"
write.csv(centre, file.path(OUT_DIR, "centre_results.csv"), row.names = FALSE)

# --- Table A5 ------------------------------------------------------------------
ct <- centre[order(centre$outcome != "Police-attended protests", centre$comparison,
                   centre$spat_window, centre$t_window), ]
n_pol  <- sum(ct$outcome == "Police-attended protests")
n_left <- sum(ct$outcome == "Police-attended protests" & ct$comparison == "Left vs centre")
n_vleft <- sum(ct$outcome == "All protests" & ct$comparison == "Left vs centre")

write_tex_table(
  data.frame(paste(ct$spat_window, "km"), paste(ct$t_window, "days"),
             fmt_pct_cell(ct$pct, ct$pct_lo, ct$pct_hi, ct$pvalue),
             fmt_n(ct$treatment_post), fmt_n(ct$control_post), fmt(ct$commonSupport_post)),
  file = file.path(OUT_DIR, "tab_A5_centre.tex"),
  caption = "Left- and right-coded protests vs matched centre-coded protests",
  label = "tab:centre",
  header = "Radius & Window & Difference (\\%) [95\\% CI] & Treated & Centre & Support (\\%)",
  align = "llcrrr",
  panel_rows = list(
    "Police-attended protests: left vs centre"  = 1,
    "Police-attended protests: right vs centre" = n_left + 1,
    "All protests: left vs centre"              = n_pol + 1,
    "All protests: right vs centre"             = n_pol + n_vleft + 1),
  notes = paste(TABLE_NOTE_PCT, "Treated: left- or right-coded events; Centre: matched centre-coded events."))

# --- Figure A2 -----------------------------------------------------------------
figA2 <- plot_pct(subset(centre, spat_window == MAIN_KM), "outcome", facet_var = "comparison",
                  ylab = "Difference vs matched centre-coded protests (%)",
                  caption = sprintf("Matched wake analysis, %d km. Bars: 95%% intervals.", MAIN_KM))
ggsave(file.path(OUT_DIR, "fig_A2_centre.pdf"), figA2, width = 9, height = 4.5)

sessionInfo()
