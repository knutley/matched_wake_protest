# =============================================================================
# 00_mwa_functions.R
# Shared settings and helper functions for the matched wake analysis (MWA)
#
# Project: Paper 2 -- MWA: Does Partisan Bias Become Institutionalised?
# Author:  Katelyn Nutley
#
# Sourced by:
#   01_mwa_setup_main.R         data, specification choice, main right-vs-left results
#   02_mwa_centre_comparison.R  left and right vs centre-coded protests
#   03_mwa_alt_baselines.R      left and right vs unknown / all other protests
#
# Data: the classified ACLED file from Paper 1 ("Who Is Spared?", under review).
# ACLED's terms of use do not allow redistribution; the build code is at
# https://github.com/knutley/selective_nonsecuritisation
# =============================================================================

# Java heap must be set BEFORE rJava/mwa load. If R has already started Java
# in this session, restart R first.
options(java.parameters = "-Xmx8g")

suppressPackageStartupMessages({
  library(rJava)
  library(mwa)       # Schutte & Donnay (2014); uses cem for matching
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(fixest)
})

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

DATA_PATH <- Sys.getenv(
  "MWA_DATA_PATH",
  "~/Documents/GitHub/selective_nonsecuritisation/data/combined/acled_merged_controls_rethresholded_v3.csv"
)
OUT_DIR <- "output"
dir.create(OUT_DIR, showWarnings = FALSE)

THRESHOLD <- 0.9          # police-presence classifier threshold (as in Paper 1)
MAIN_KM   <- 4            # main radius; see appendix for the choice
RADII     <- c(2, 4, 6)   # radii reported in the appendix grids (equally spaced)
WINDOWS   <- c(7, 30, 90) # cumulative windows: immediate, intermediate, longer

# Main matching specification: exact matching (CEM) on country, year and
# weekend, plus banded distance-to-government-building and same-day protest load
MATCH_COLS <- c("country", "year", "is_weekend", "dist_band", "load_band")

# Alternative matching sets, used only in the specification comparison
MATCH_COLS_FULL <- c("country", "year", "is_weekend", "incumbent_left", "incumbent_right",
                     "log_dist_govt_building", "log_dist_major_road", "log_protest_load")
MATCH_COLS_LEAN <- c("country", "year", "is_weekend", "log_dist_govt_building", "log_protest_load")

COMPOUND <- c("left_and_right", "centre_and_left", "centre_and_right")

OUTCOME_LABELS <- c(police_present = "Police-attended protests",
                    any_protest    = "All protests")

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

#' Build the MWA dataset from the Paper 1 file.
#' One row per protest event. Adds: timestamp/lat/lon (required by mwa),
#' grouping columns for each comparison, outcome indicators and matching
#' covariates (all numeric, as cem requires).
build_mwa_data <- function(path = DATA_PATH) {
  raw <- read_csv(path, show_col_types = FALSE)

  raw %>%
    # Same-day protest load across ALL events in the country (Paper 1's
    # exclusion restriction); computed before any filtering
    group_by(country, event_date) %>%
    mutate(protest_load = n() - 1L) %>%
    ungroup() %>%
    mutate(
      timestamp = paste(as.character(as.Date(event_date)), "00:00:00"),
      lat = latitude,
      lon = longitude,

      # Outcomes (dependent events)
      police_present = as.integer(police_presence_prob >= THRESHOLD),
      any_protest    = 1L,

      # Treatment / control groupings, one column per comparison
      side = case_when(                        # right vs left
        event_partisan_type_final == "right" ~ "right",
        event_partisan_type_final == "left"  ~ "left",
        TRUE ~ "other"),
      side_centre = case_when(                 # vs centre
        event_partisan_type_final %in% c("right", "left") ~ event_partisan_type_final,
        event_partisan_type_final == "centre" ~ "base",
        TRUE ~ "other"),
      side_unknown = case_when(                # vs unclassified events
        event_partisan_type_final %in% c("right", "left") ~ event_partisan_type_final,
        event_partisan_type_final == "unknown" ~ "base",
        TRUE ~ "other"),
      side_others = case_when(                 # vs everything else
        event_partisan_type_final %in% c("right", "left") ~ event_partisan_type_final,
        TRUE ~ "base"),                        # centre + unknown + compound

      # Matching covariates
      country_name           = country,
      country                = as.numeric(factor(country)),
      year                   = as.numeric(year),
      is_weekend             = as.numeric(is_weekend),
      incumbent_left         = as.numeric(incumbent_left),
      incumbent_right        = as.numeric(incumbent_right),
      log_dist_govt_building = log1p(dist_govt_building_m),
      log_dist_major_road    = log1p(dist_major_road_m),
      log_protest_load       = log1p(protest_load)
    ) %>%
    filter(!is.na(lat), !is.na(lon)) %>%
    mutate(
      dist_band = ntile(log_dist_govt_building, 4),   # quartiles
      load_band = ntile(log_protest_load, 3)          # quiet / normal / busy days
    ) %>%
    as.data.frame()
}

#' Load the cached dataset if 01 has already built it; otherwise build it.
load_mwa_data <- function(cache = file.path(OUT_DIR, "mwa_data.rds")) {
  if (file.exists(cache)) return(readRDS(cache))
  d <- build_mwa_data()
  saveRDS(d, cache)
  d
}

# -----------------------------------------------------------------------------
# Running and summarising MWA
# -----------------------------------------------------------------------------

#' One matchedwake() call for a single time window and one or more radii.
#' group_col holds the treatment value, the control value and anything else.
run_mwa <- function(data, group_col, treat, control, dependent, days,
                    radii = MAIN_KM, match_cols = MATCH_COLS,
                    tcm = FALSE, delete_sutva = FALSE) {
  stopifnot(all(sapply(data[match_cols], is.numeric)))
  spat <- if (length(radii) == 1) c(radii, radii, 1) else {
    step <- unique(diff(radii))
    stopifnot(length(step) == 1)            # radii must be equally spaced
    c(min(radii), max(radii), step)
  }
  message(sprintf("MWA: %s=%s vs %s | %s | %d days | %s km",
                  group_col, treat, control, dependent, days, paste(radii, collapse = ",")))
  matchedwake(
    data, t_window = c(days, days, 1), spat_window = spat,
    treatment = c(group_col, treat), control = c(group_col, control),
    dependent = c(dependent, 1), matchColumns = match_cols,
    t_unit = "days", TCM = tcm, deleteSUTVA = delete_sutva
  )
}

#' Tidy one matchedwake result: one row per (time window, radius).
#' pct = estimate as % of the mean post-period count around control events;
#' the SE is backed out of the reported estimate and p-value (normal approx.).
summarise_mwa <- function(res, ...) {
  keys <- c("t_window", "spat_window")
  base <- aggregate(dependent_post ~ t_window + spat_window,
                    data = subset(res$wakes, treatment == 0), FUN = mean)
  names(base)[3] <- "control_mean_post"
  out <- Reduce(function(a, b) merge(a, b, by = keys), list(
    res$estimates[, c(keys, "estimate", "pvalue")],
    res$matching[,  c(keys, "control_post", "treatment_post", "L1_pre", "L1_post", "commonSupport_post")],
    res$SUTVA[,     c(keys, "SO", "MO")],
    base
  ))
  out$se     <- abs(out$estimate) / qnorm(1 - out$pvalue / 2)
  out$pct    <- 100 * out$estimate / out$control_mean_post
  out$pct_lo <- 100 * (out$estimate - 1.96 * out$se) / out$control_mean_post
  out$pct_hi <- 100 * (out$estimate + 1.96 * out$se) / out$control_mean_post
  labels <- list(...)
  for (nm in names(labels)) out[[nm]] <- labels[[nm]]
  out
}

#' Run a comparison over windows x outcomes (each call covers all radii).
#' Returns list(table = tidy results, raw = matchedwake objects for reuse).
run_grid <- function(data, group_col, treat, control, comparison,
                     days = WINDOWS, radii = RADII,
                     outcomes = names(OUTCOME_LABELS), ...) {
  raw <- list(); tab <- list()
  for (d in days) for (o in outcomes) {
    key <- paste(o, d, sep = "_")
    raw[[key]] <- run_mwa(data, group_col, treat, control, o, d, radii, ...)
    tab[[key]] <- summarise_mwa(raw[[key]], comparison = comparison,
                                outcome = OUTCOME_LABELS[[o]])
  }
  list(table = bind_rows(tab), raw = raw)
}

#' Police presence per follow-up protest: Poisson model of post-period
#' police-attended events with the post-period count of all protests as
#' exposure (offset), controlling for pre-period counts. exp(coef) on
#' treatment is the rate ratio. With fe = TRUE, adds fixed effects for the
#' matching strata.
rate_model <- function(res_police, res_protest, km = MAIN_KM, fe = TRUE) {
  keys <- c("eventID", "t_window", "spat_window")
  wp <- subset(res_police$wakes,  spat_window == km)
  wa <- subset(res_protest$wakes, spat_window == km)
  fe_cols <- intersect(MATCH_COLS, names(wp))
  panel <- merge(
    wp[, c(keys, "treatment", "dependent_pre", "dependent_post", fe_cols)],
    wa[, c(keys, "dependent_pre", "dependent_post")],
    by = keys, suffixes = c("_police", "_protest")
  )
  panel <- subset(panel, dependent_post_protest > 0)
  rhs <- "treatment + log1p(dependent_pre_police) + log1p(dependent_pre_protest)"
  if (fe && length(fe_cols)) {
    fe_part <- if (all(c("country", "year") %in% fe_cols))
      paste(c("country^year", setdiff(fe_cols, c("country", "year"))), collapse = " + ")
    else paste(fe_cols, collapse = " + ")
    rhs <- paste(rhs, "|", fe_part)
  }
  m <- fepois(as.formula(paste("dependent_post_police ~", rhs)),
              offset = ~ log(dependent_post_protest), data = panel, vcov = "hetero")
  ci <- confint(m)["treatment", ]
  data.frame(t_window = unique(wp$t_window), spat_window = km, fixed_effects = fe,
             rate_ratio = unname(exp(coef(m)["treatment"])),
             rr_lo = exp(ci[[1]]), rr_hi = exp(ci[[2]]),
             pvalue = unname(pvalue(m)["treatment"]), n_cylinders = nobs(m))
}

# -----------------------------------------------------------------------------
# Output helpers: figures and LaTeX tables
# -----------------------------------------------------------------------------

#' Point-range plot of % differences by window.
plot_pct <- function(tab, colour_var, facet_var = NULL, ylab, caption,
                     colours = c("Police-attended protests" = "#b2182b", "All protests" = "grey45")) {
  p <- ggplot(tab, aes(x = factor(t_window), y = pct, colour = .data[[colour_var]])) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_pointrange(aes(ymin = pct_lo, ymax = pct_hi), position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = colours) +
    labs(x = "Days after protest", y = ylab, colour = NULL, caption = caption) +
    theme_minimal(base_size = 12) + theme(legend.position = "top")
  if (!is.null(facet_var)) p <- p + facet_wrap(as.formula(paste("~", facet_var)))
  p
}

stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "$^{***}$",
                     ifelse(p < 0.01, "$^{**}$", ifelse(p < 0.05, "$^{*}$", ""))))
fmt   <- function(x, d = 1) {               # LaTeX minus sign for negatives
  out <- formatC(abs(x), format = "f", digits = d, big.mark = ",")
  ifelse(is.na(x), "--", ifelse(x < 0 & round(abs(x), d) != 0, paste0("$-$", out), out))
}
fmt_n <- function(x) ifelse(is.na(x), "--", formatC(x, format = "d", big.mark = ","))
fmt_pct_cell <- function(pct, lo, hi, p)
  sprintf("%s%s [%s, %s]", fmt(pct), stars(p), fmt(lo), fmt(hi))

#' Write a booktabs/threeparttable LaTeX table.
#' body: character matrix or data frame (already formatted); header: one or
#' more LaTeX header lines (without trailing \\); notes: table note text.
write_tex_table <- function(body, file, caption, label, header, align, notes = NULL,
                            panel_rows = NULL) {
  body <- as.matrix(body)
  rows <- apply(body, 1, function(r) paste(paste(r, collapse = " & "), "\\\\"))
  if (!is.null(panel_rows)) {             # named list: row index -> panel title
    for (i in rev(seq_along(panel_rows))) {
      at <- panel_rows[[i]]
      title <- sprintf("\\addlinespace\\multicolumn{%d}{l}{\\textit{%s}} \\\\",
                       ncol(body), names(panel_rows)[i])
      rows <- append(rows, title, after = at - 1)
    }
  }
  lines <- c(
    "\\begin{table}[htbp]", "\\centering", "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}",
    "\\begin{threeparttable}",
    sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
    sprintf("\\begin{tabular}{%s}", align), "\\toprule",
    paste(header, "\\\\"), "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
    if (!is.null(notes)) c("\\begin{tablenotes}\\footnotesize", sprintf("\\item %s", notes), "\\end{tablenotes}"),
    "\\end{threeparttable}", "\\end{table}"
  )
  writeLines(lines, file)
  message("Wrote ", file)
}

TABLE_NOTE_PCT <- paste(
  "Cells give the difference in post-period counts between treatment and matched",
  "control events, as a percentage of the mean post-period count around control",
  "events, with 95\\% intervals in brackets. Windows are cumulative (days 1 to $t$).",
  "Matching: coarsened exact matching on country, year, weekend, distance to the",
  "nearest government building (quartiles) and same-day protest load (terciles).",
  "$^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.")
