# Author: Katelyn Nutley
# Date: 26-10-2026
# Title: MWA Setup 

# Note: The data for this article comes from a previous project, "Who Is Spared? Selective Non-Securitisation
# and the Partisan Policing of Protest in Western Europe" (under review). In this article, observational data 
# from ACLED was effectively classified and augmented. Given ACLED's UA, I cannot share their data, but the 
# process is replicable and the code is available here: https://github.com/knutley/selective_nonsecuritisation

# Required libraries

library(dplyr)
library(httr)
library(jsonlite)
library(readr)

# Need treatment, dependent, and control events. So, here's the research design walk-through: 
# https://matchedwake.com/aid_ex/ ; according to this, 

# Treatment Events = the events whose aftermath you want to study 
  # -> I'm broadly focused on partisanship, but only the right-coded protests achieved
  # significance, so I think they should probably be treatment
# Control Events = comparable events that lack the treatment feature, matched on 
# covariates so they acts as a quasi-counterfactual 
  # -> matched left-coded protests, then 
# Depdendent Events = (dependent after - dependent before) around treatment 
# vs when the same change around control (difference-in-differences)

# Events where police didn't show up aren't necessarily a separate category; they're just 
# protests that don't count as dependent events? (So, I'm going to remove them.)

# Here's the R package for MWA: https://cran.r-project.org/web/packages/mwa/mwa.pdf

install.packages("mwa")
install.packages("rJava")
options(java.parameters = "-Xmx8g")   # set BEFORE library(mwa); 100k events needs more than the 512 MB default
library(rJava)
library(mwa)

# Load data (again, see above for how to get the data -- alternatively, reach out and I'll share!)

THRESHOLD <- 0.9   # same presence threshold as Paper 1

og_data <- read_csv("~/Documents/GitHub/selective_nonsecuritisation/data/combined/acled_merged_controls_rethresholded_v3.csv",
                    show_col_types = FALSE)

mwa_data <- og_data %>%
  group_by(country, event_date) %>%
  mutate(protest_load = n() - 1L) %>%          # same-day protest load, all events
  ungroup() %>%
  mutate(
    timestamp = paste(as.character(as.Date(event_date)), "00:00:00"),
    lat = latitude,
    lon = longitude,
    side = case_when(
      event_partisan_type_final == "right" ~ "right",   # treatment
      event_partisan_type_final == "left"  ~ "left",    # control
      TRUE ~ "other"
    ),
    police_present = as.integer(police_presence_prob >= THRESHOLD),
    police_absent  = 1L - police_present,
    any_protest    = 1L,
    log_dist_govt_building = log1p(dist_govt_building_m),
    log_dist_major_road    = log1p(dist_major_road_m),
    log_protest_load       = log1p(protest_load)
  ) %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  as.data.frame()

# --- 3. Matching columns, converted to numeric BEFORE any matchedwake() call ---
match_cols      <- c("country", "year", "is_weekend", "incumbent_left", "incumbent_right",
                     "log_dist_govt_building", "log_dist_major_road", "log_protest_load")
match_cols_lean <- c("country", "year", "is_weekend", "log_dist_govt_building", "log_protest_load")

mwa_data <- mwa_data %>%
  mutate(across(all_of(match_cols), ~ as.numeric(if (is.character(.x) || is.factor(.x)) factor(.x) else .x)))
stopifnot(all(sapply(mwa_data[match_cols], is.numeric)))

# --- 4. Balance diagnostic: 30 days, 6 km -------------------------------------
res_bal <- matchedwake(
  mwa_data, t_window = c(30, 30, 1), spat_window = c(6, 6, 1),
  treatment = c("side", "right"), control = c("side", "left"),
  dependent = c("police_present", 1), matchColumns = match_cols_lean,
  t_unit = "days", TCM = TRUE, deleteSUTVA = FALSE
)

res_bal$matching[, c("L1_pre", "L1_post", "commonSupport_post")]

w <- res_bal$wakes
print(names(w))   # paste this if the next step errors

# The treatment flag column (1 = right-coded, 0 = left-coded)
treat_col <- intersect(c("treatment", "treat", "Treatment"), names(w))[1]
if (is.na(treat_col)) stop("Set treat_col by hand to the treatment-flag column printed above")

# Every numeric column except the flag, the event ID, the outcome counts and the overlap columns:
# that leaves the matching variables plus the prior-count columns TCM adds
skip <- grep("^(dependent|event|id$)|overlap|SO|MO", names(w), value = TRUE, ignore.case = TRUE)
bal_vars <- setdiff(names(w)[sapply(w, is.numeric)], c(treat_col, skip))
print(bal_vars)   # check this is the matching set (+ TCM counts)

imbalance(group = w[[treat_col]], data = w[, bal_vars])

# --- 5. Balance fixes: run these three and compare L1 / common support ---------
# (A) Same as above but WITHOUT the prior-count matching (TCM)
res_A <- matchedwake(
  mwa_data, t_window = c(30, 30, 1), spat_window = c(6, 6, 1),
  treatment = c("side", "right"), control = c("side", "left"),
  dependent = c("police_present", 1), matchColumns = match_cols_lean,
  t_unit = "days", TCM = FALSE, deleteSUTVA = FALSE
)

# (B) Coarsen the continuous matching variables into bands first
mwa_data <- mwa_data %>%
  mutate(
    dist_band = ntile(log_dist_govt_building, 4),   # quartiles of distance to govt building
    load_band = ntile(log_protest_load, 3)          # quiet / normal / busy days
  )
match_cols_binned <- c("country", "year", "is_weekend", "dist_band", "load_band")

res_B <- matchedwake(
  mwa_data, t_window = c(30, 30, 1), spat_window = c(6, 6, 1),
  treatment = c("side", "right"), control = c("side", "left"),
  dependent = c("police_present", 1), matchColumns = match_cols_binned,
  t_unit = "days", TCM = FALSE, deleteSUTVA = FALSE
)

# (C) Bands, with TCM back on
res_C <- matchedwake(
  mwa_data, t_window = c(30, 30, 1), spat_window = c(6, 6, 1),
  treatment = c("side", "right"), control = c("side", "left"),
  dependent = c("police_present", 1), matchColumns = match_cols_binned,
  t_unit = "days", TCM = TRUE, deleteSUTVA = FALSE
)

compare <- rbind(
  original      = cbind(res_bal$matching, res_bal$estimates[, c("estimate", "pvalue")]),
  A_noTCM       = cbind(res_A$matching,   res_A$estimates[,   c("estimate", "pvalue")]),
  B_bands_noTCM = cbind(res_B$matching,   res_B$estimates[,   c("estimate", "pvalue")]),
  C_bands_TCM   = cbind(res_C$matching,   res_C$estimates[,   c("estimate", "pvalue")])
)
print(compare[, c("L1_pre", "L1_post", "commonSupport_post",
                  "treatment_post", "control_post", "estimate", "pvalue")], digits = 3)


res_bal$matching[, c("L1_pre", "L1_post", "commonSupport_post")]

w <- res_bal$wakes
print(names(w))   # paste this if the next step errors

# The treatment flag column (1 = right-coded, 0 = left-coded)
treat_col <- intersect(c("treatment", "treat", "Treatment"), names(w))[1]
if (is.na(treat_col)) stop("Set treat_col by hand to the treatment-flag column printed above")

# Every numeric column except the flag, the event ID, the outcome counts and the overlap columns:
# that leaves the matching variables plus the prior-count columns TCM adds
skip <- grep("^(dependent|event|id$)|overlap|SO|MO", names(w), value = TRUE, ignore.case = TRUE)
bal_vars <- setdiff(names(w)[sapply(w, is.numeric)], c(treat_col, skip))
print(bal_vars)   # check this is the matching set (+ TCM counts)

imbalance(group = w[[treat_col]], data = w[, bal_vars])

res_B_radius <- matchedwake(
  mwa_data, t_window = c(30, 30, 1), spat_window = c(2, 10, 2),   # 2, 4, 6, 8, 10 km
  treatment = c("side", "right"), control = c("side", "left"),
  dependent = c("police_present", 1), matchColumns = match_cols_binned,
  t_unit = "days", TCM = FALSE, deleteSUTVA = FALSE
)

# approximate baseline: average police-present events after left-coded protests, per radius
base <- aggregate(dependent_post ~ spat_window, data = subset(res_B_radius$wakes, treatment == 0), FUN = mean)
names(base)[2] <- "control_mean_post"

radius_tab <- Reduce(function(a, b) merge(a, b, by = c("t_window", "spat_window")),
                     list(res_B_radius$estimates[, c("t_window", "spat_window", "estimate", "pvalue")],
                          res_B_radius$matching[, c("t_window", "spat_window", "L1_post", "commonSupport_post", "treatment_post", "control_post")],
                          res_B_radius$SUTVA[, c("t_window", "spat_window", "SO", "MO")]))
radius_tab <- merge(radius_tab, base, by = "spat_window")
radius_tab$pct_of_baseline <- round(100 * radius_tab$estimate / radius_tab$control_mean_post, 1)
print(radius_tab, digits = 3)

# Author: Katelyn Nutley
# Title: MWA next steps -- main specification (B settings, 4 km)
#
# Run AFTER mwa_balance_check.R sections 1-3 in the same session
# (needs mwa_data, the numeric conversion, and the bands below).
#
#   Step 1: time windows 7 / 30 / 90 days at 4 km  -> widening or narrowing?
#   Step 2: protest-volume check (per-protest rate model) at 30 days, 4 km
#   Step 3: paper figure (gap as % of baseline, with 95% intervals)

library(dplyr)
library(fixest)
library(ggplot2)

# --- Main specification ---------------------------------------------------------
if (!"dist_band" %in% names(mwa_data)) {
  mwa_data <- mwa_data %>%
    mutate(dist_band = ntile(log_dist_govt_building, 4),
           load_band = ntile(log_protest_load, 3))
}
match_cols_binned <- c("country", "year", "is_weekend", "dist_band", "load_band")
MAIN_KM <- 4

run_B <- function(days, km, dependent_col) {
  matchedwake(
    mwa_data, t_window = c(days, days, 1), spat_window = c(km, km, 1),
    treatment = c("side", "right"), control = c("side", "left"),
    dependent = c(dependent_col, 1), matchColumns = match_cols_binned,
    t_unit = "days", TCM = FALSE, deleteSUTVA = FALSE
  )
}

summarise_B <- function(res, label) {
  base <- mean(res$wakes$dependent_post[res$wakes$treatment == 0])
  est  <- res$estimates$estimate
  p    <- res$estimates$pvalue
  se   <- abs(est) / qnorm(1 - p / 2)        # back out the SE from estimate and p-value
  data.frame(
    days = res$estimates$t_window, km = res$estimates$spat_window, outcome = label,
    estimate = est, pvalue = p, se = se,
    control_mean_post = base,
    pct = 100 * est / base,
    pct_lo = 100 * (est - 1.96 * se) / base,
    pct_hi = 100 * (est + 1.96 * se) / base,
    treatment_n = res$matching$treatment_post, control_n = res$matching$control_post,
    common_support = res$matching$commonSupport_post,
    MO = res$SUTVA$MO
  )
}

# --- Step 1: 7 / 30 / 90 days ----------------------------------------------------
runs <- list()
for (d in c(7, 30, 90)) {
  runs[[paste0("police_", d)]]  <- run_B(d, MAIN_KM, "police_present")
  runs[[paste0("protest_", d)]] <- run_B(d, MAIN_KM, "any_protest")
}

time_tab <- bind_rows(
  lapply(c(7, 30, 90), function(d) summarise_B(runs[[paste0("police_", d)]],  "Police-attended protests")),
  lapply(c(7, 30, 90), function(d) summarise_B(runs[[paste0("protest_", d)]], "All protests"))
)
print(time_tab, digits = 3)
write.csv(time_tab, "mwa_main_time_windows.csv", row.names = FALSE)

# --- Step 2: per-protest rate (is it police, or just fewer protests?) -----------
rate_model <- function(d) {
  wp <- runs[[paste0("police_", d)]]$wakes
  wa <- runs[[paste0("protest_", d)]]$wakes
  panel <- merge(
    wp[, c("eventID", "treatment", "dependent_pre", "dependent_post")],
    wa[, c("eventID", "dependent_pre", "dependent_post")],
    by = "eventID", suffixes = c("_police", "_protest")
  )
  m <- fepois(
    dependent_post_police ~ treatment + log1p(dependent_pre_police) + log1p(dependent_pre_protest),
    offset = ~ log(dependent_post_protest),
    data   = subset(panel, dependent_post_protest > 0),
    vcov   = "hetero"
  )
  ci <- confint(m)["treatment", ]
  data.frame(days = d,
             rate_ratio = exp(coef(m)["treatment"]),
             lo = exp(ci[[1]]), hi = exp(ci[[2]]),
             pvalue = pvalue(m)["treatment"],
             n_cylinders = nobs(m))
}
rate_tab <- bind_rows(lapply(c(7, 30, 90), rate_model))
print(rate_tab, digits = 3)
# rate_ratio = police presence PER FOLLOW-UP PROTEST after right- vs left-coded protests.
# 0.85 means police attend 15% fewer of the protests that follow. ~1 means the gap was protest volume.
write.csv(rate_tab, "mwa_rate_model.csv", row.names = FALSE)

# --- Step 3: paper figure ---------------------------------------------------------
fig <- ggplot(time_tab, aes(x = factor(days), y = pct, colour = outcome)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_pointrange(aes(ymin = pct_lo, ymax = pct_hi),
                  position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("Police-attended protests" = "#b2182b", "All protests" = "grey45")) +
  labs(x = "Days after protest (cumulative window)",
       y = "Right vs left-coded difference\n(% of left-coded baseline)",
       colour = NULL,
       caption = paste0("Matched wake analysis, ", MAIN_KM,
                        " km radius. Points are estimates; bars are 95% intervals.")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
print(fig)
ggsave("mwa_main_figure.pdf", fig, width = 7, height = 4.5)

####### 

# Author: Katelyn Nutley
# Title: MWA next steps -- main specification (B settings, 4 km)
#
# Run AFTER mwa_balance_check.R sections 1-3 in the same session
# (needs mwa_data, the numeric conversion, and the bands below).
#
#   Step 1: time windows 7 / 30 / 90 days at 4 km  -> widening or narrowing?
#   Step 2: protest-volume check (per-protest rate model) at 30 days, 4 km
#   Step 3: paper figure (gap as % of baseline, with 95% intervals)

library(dplyr)
library(fixest)
library(ggplot2)

# --- Main specification ---------------------------------------------------------
if (!"dist_band" %in% names(mwa_data)) {
  mwa_data <- mwa_data %>%
    mutate(dist_band = ntile(log_dist_govt_building, 4),
           load_band = ntile(log_protest_load, 3))
}
match_cols_binned <- c("country", "year", "is_weekend", "dist_band", "load_band")
MAIN_KM <- 4

run_B <- function(days, km, dependent_col) {
  matchedwake(
    mwa_data, t_window = c(days, days, 1), spat_window = c(km, km, 1),
    treatment = c("side", "right"), control = c("side", "left"),
    dependent = c(dependent_col, 1), matchColumns = match_cols_binned,
    t_unit = "days", TCM = FALSE, deleteSUTVA = FALSE
  )
}

summarise_B <- function(res, label) {
  base <- mean(res$wakes$dependent_post[res$wakes$treatment == 0])
  est  <- res$estimates$estimate
  p    <- res$estimates$pvalue
  se   <- abs(est) / qnorm(1 - p / 2)        # back out the SE from estimate and p-value
  data.frame(
    days = res$estimates$t_window, km = res$estimates$spat_window, outcome = label,
    estimate = est, pvalue = p, se = se,
    control_mean_post = base,
    pct = 100 * est / base,
    pct_lo = 100 * (est - 1.96 * se) / base,
    pct_hi = 100 * (est + 1.96 * se) / base,
    treatment_n = res$matching$treatment_post, control_n = res$matching$control_post,
    common_support = res$matching$commonSupport_post,
    MO = res$SUTVA$MO
  )
}

# --- Step 1: 7 / 30 / 90 days ----------------------------------------------------
runs <- list()
for (d in c(7, 30, 90)) {
  runs[[paste0("police_", d)]]  <- run_B(d, MAIN_KM, "police_present")
  runs[[paste0("protest_", d)]] <- run_B(d, MAIN_KM, "any_protest")
}

time_tab <- bind_rows(
  lapply(c(7, 30, 90), function(d) summarise_B(runs[[paste0("police_", d)]],  "Police-attended protests")),
  lapply(c(7, 30, 90), function(d) summarise_B(runs[[paste0("protest_", d)]], "All protests"))
)
print(time_tab, digits = 3)
write.csv(time_tab, "mwa_main_time_windows.csv", row.names = FALSE)

# --- Step 1b: full main table, 2 / 4 / 6 km x 7 / 30 / 90 days x both outcomes ---
# One matchedwake() call per (days, outcome) covers all three radii: 6 calls in total.
summarise_multi <- function(res, label) {
  keys <- c("t_window", "spat_window")
  base <- aggregate(dependent_post ~ t_window + spat_window,
                    data = subset(res$wakes, treatment == 0), FUN = mean)
  names(base)[3] <- "control_mean_post"
  out <- Reduce(function(a, b) merge(a, b, by = keys), list(
    res$estimates[, c(keys, "estimate", "pvalue")],
    res$matching[,  c(keys, "treatment_post", "control_post", "commonSupport_post")],
    res$SUTVA[,     c(keys, "MO")],
    base
  ))
  out$se     <- abs(out$estimate) / qnorm(1 - out$pvalue / 2)
  out$pct    <- 100 * out$estimate / out$control_mean_post
  out$pct_lo <- 100 * (out$estimate - 1.96 * out$se) / out$control_mean_post
  out$pct_hi <- 100 * (out$estimate + 1.96 * out$se) / out$control_mean_post
  out$outcome <- label
  out
}

run_B_radii <- function(days, dependent_col) {
  matchedwake(
    mwa_data, t_window = c(days, days, 1), spat_window = c(2, 6, 2),   # 2, 4, 6 km
    treatment = c("side", "right"), control = c("side", "left"),
    dependent = c(dependent_col, 1), matchColumns = match_cols_binned,
    t_unit = "days", TCM = FALSE, deleteSUTVA = FALSE
  )
}

main_tab <- bind_rows(lapply(c(7, 30, 90), function(d) bind_rows(
  summarise_multi(run_B_radii(d, "police_present"), "Police-attended protests"),
  summarise_multi(run_B_radii(d, "any_protest"),    "All protests")
)))
main_tab <- main_tab[order(main_tab$outcome, main_tab$spat_window, main_tab$t_window), ]
print(main_tab[, c("outcome", "spat_window", "t_window", "estimate", "pvalue",
                   "pct", "pct_lo", "pct_hi", "treatment_post", "commonSupport_post", "MO")], digits = 3)
write.csv(main_tab, "mwa_main_table.csv", row.names = FALSE)

# Figure: one panel per radius
fig_grid <- ggplot(main_tab, aes(x = factor(t_window), y = pct, colour = outcome)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_pointrange(aes(ymin = pct_lo, ymax = pct_hi), position = position_dodge(width = 0.4)) +
  facet_wrap(~ paste(spat_window, "km")) +
  scale_colour_manual(values = c("Police-attended protests" = "#b2182b", "All protests" = "grey45")) +
  labs(x = "Days after protest (cumulative window)",
       y = "Right vs left-coded difference\n(% of left-coded baseline)", colour = NULL,
       caption = "Matched wake analysis. Points are estimates; bars are 95% intervals.") +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
print(fig_grid)
ggsave("mwa_main_figure_by_radius.pdf", fig_grid, width = 9, height = 4.5)

# --- Step 2: per-protest rate (is it police, or just fewer protests?) -----------
rate_model <- function(d) {
  wp <- runs[[paste0("police_", d)]]$wakes
  wa <- runs[[paste0("protest_", d)]]$wakes
  panel <- merge(
    wp[, c("eventID", "treatment", "dependent_pre", "dependent_post")],
    wa[, c("eventID", "dependent_pre", "dependent_post")],
    by = "eventID", suffixes = c("_police", "_protest")
  )
  m <- fepois(
    dependent_post_police ~ treatment + log1p(dependent_pre_police) + log1p(dependent_pre_protest),
    offset = ~ log(dependent_post_protest),
    data   = subset(panel, dependent_post_protest > 0),
    vcov   = "hetero"
  )
  ci <- confint(m)["treatment", ]
  data.frame(days = d,
             rate_ratio = exp(coef(m)["treatment"]),
             lo = exp(ci[[1]]), hi = exp(ci[[2]]),
             pvalue = pvalue(m)["treatment"],
             n_cylinders = nobs(m))
}
rate_tab <- bind_rows(lapply(c(7, 30, 90), rate_model))
print(rate_tab, digits = 3)
# rate_ratio = police presence PER FOLLOW-UP PROTEST after right- vs left-coded protests.
# 0.85 means police attend 15% fewer of the protests that follow. ~1 means the gap was protest volume.
write.csv(rate_tab, "mwa_rate_model.csv", row.names = FALSE)

# --- Step 3: paper figure ---------------------------------------------------------
fig <- ggplot(time_tab, aes(x = factor(days), y = pct, colour = outcome)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_pointrange(aes(ymin = pct_lo, ymax = pct_hi),
                  position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("Police-attended protests" = "#b2182b", "All protests" = "grey45")) +
  labs(x = "Days after protest (cumulative window)",
       y = "Right vs left-coded difference\n(% of left-coded baseline)",
       colour = NULL,
       caption = paste0("Matched wake analysis, ", MAIN_KM,
                        " km radius. Points are estimates; bars are 95% intervals.")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
print(fig)
ggsave("mwa_main_figure.pdf", fig, width = 7, height = 4.5)


