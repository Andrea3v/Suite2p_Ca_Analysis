## this is based on peaks detected on deltaFoverF (dFF, lowpass filtered) calcium traces -> report amplitude/prominence
# FWHM is extracted from '[]_burst_events.csv' where the FWHM is calc on the burst peaks but reports the amplitude
# from dFF, ignore FWHM reported here...

library(readr)
library(dplyr)
library(purrr)
library(fs)
# install.packages("see")
library(see)
# install.packages("tidyverse")  # if not already installed")
library(tidyverse)
# install.packages("patchwork")  # only if you want them stacked
library(patchwork)
library(rstatix)   # wilcox_test, adjust_pvalue, add_xy_position
library(ggpubr)    # stat_pvalue_manual

collate_calcium_peaks <- function(path, pattern = "detected_peaks_calcium\\.csv$") {
  files <- list.files(
    path        = path,
    pattern     = pattern,
    full.names  = TRUE,
    recursive   = FALSE
  )
  
  if (length(files) == 0) {
    stop("No files matching '", pattern, "' found in: ", path)
  }
  
  message("Found ", length(files), " file(s). Reading...")
  
  files |>
    set_names(\(f) sub("_?detected_peaks_calcium\\.csv$", "", basename(f))) |>
    map_dfr(read_csv, show_col_types = FALSE, .id = "source")
}

# --- load all data ---
data_dir <- r"(C:\Users\andre\Documents\Boaz_Ca_Analysis\Data_Jun2026)"
df <- collate_calcium_peaks(data_dir)

df <- df |>
  tidyr::separate(
    col   = source,
    into  = c("age", "location", NA),   # NA silently drops the date portion
    sep   = "_",
    extra = "drop",    # any tokens beyond the 3rd are dropped
    remove = FALSE     # keep the original source column
  )
df <- df |>
  mutate(
    age_num   = as.numeric(str_extract(age, "[0-9]+")),   # extracts 14, 15, 18, 19
    age_group = case_when(
      age_num <= 15 ~ "w14-15",
      age_num >= 18 ~ "w18-19",
      .default = NA_character_
    )
  )

glimpse(df)

#### plot

# metrics <- c("PeakHeight", "Prominence", "FWHM_s")
metrics <- c("PeakHeight", "Prominence")
df_long <- df |>
  select(source, Cell, age, location, age_group, all_of(metrics)) |>
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") |>
  mutate(
    metric    = factor(metric,    levels = metrics),
    location  = factor(location,  levels = c("apex", "base")),
    age_group = factor(age_group, levels = c("w14-15", "w18-19"))
  )

# ---- Plot 1: average within each (source, Cell) first ----
df_cell_means <- df_long |>
  group_by(source, Cell, metric) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

# remove outlies for better visualization 

# --- outlier helper (Tukey 1.5x IQR) ---
filter_iqr <- function(x, k = 1.5) {
  q   <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- diff(q)
  x >= q[1] - k * iqr & x <= q[2] + k * iqr
}

# Clean raw peaks, per metric
df_long_clean <- df_long |>
  group_by(metric) |>
  filter(filter_iqr(value)) |>
  ungroup()

# How many got dropped?
df_long |>
  count(metric, name = "before") |>
  left_join(count(df_long_clean, metric, name = "after"), by = "metric") |>
  mutate(removed = before - after,
         pct_removed = round(100 * removed / before, 1)) |>
  print()

# set metrics labels for plot
metric_labels <- c(
  PeakHeight = "Peak height (ΔF/F0)",
  Prominence = "Prominence (ΔF/F0)"
  # FWHM_s     = "FWHM (s)"
)

# Recompute per-cell means from cleaned data
df_cell_means_clean <- df_long_clean |>
  group_by(source, Cell, metric, age_group, location) |>
  summarise(value = median(value, na.rm = TRUE), .groups = "drop")

# ── Statistics (Task 2) ───────────────────────────────────────────────────────
# Design: 2×2 location × age_group, non-normal data, n = 5 experiments
# Test:   Wilcoxon rank-sum (Mann-Whitney U), BH FDR correction
#
# Two bracket sets per facet:
#   A) within-x  — age_group effect AT each location   (short dodge-level bracket)
#   B) across-x  — location effect WITHIN each age_group (bracket spanning x-axis)
# Computed per metric so y-positions respect each facet's free y-scale.

dodge_w <- 0.8
age_pal <- c("w14-15" = "#66C2A5", "w18-19" = "#FC8D62")  # Set2 col 1 & 2

run_wilcox_stats <- function(data, dodge = 0.8) {
  purrr::map_dfr(levels(data$metric), function(m) {
    d <- dplyr::filter(data, metric == m)

    # A) within each location: w14-15 vs w18-19
    st_A <- d |>
      dplyr::group_by(location) |>
      rstatix::wilcox_test(value ~ age_group) |>
      rstatix::adjust_pvalue(method = "BH") |>
      rstatix::add_significance() |>
      rstatix::add_xy_position(x = "location", dodge = dodge,
                               data = d, fun = "max") |>
      dplyr::mutate(metric = factor(m, levels = levels(data$metric)),
                    comparison_type = "within_x")

    # B) within each age_group: apex vs base
    st_B_raw <- d |>
      dplyr::group_by(age_group) |>
      rstatix::wilcox_test(value ~ location) |>
      rstatix::adjust_pvalue(method = "BH") |>
      rstatix::add_significance() |>
      rstatix::add_xy_position(x = "location", dodge = dodge,
                               group = "age_group",
                               data = d, fun = "max") |>
      dplyr::mutate(metric = factor(m, levels = levels(data$metric)),
                    comparison_type = "across_x")

    # stack B above A so brackets don't overlap
    y_ceil <- max(st_A$y.position, na.rm = TRUE)
    y_step <- 0.10 * diff(range(d$value, na.rm = TRUE))
    st_B   <- st_B_raw |>
      dplyr::arrange(age_group) |>
      dplyr::mutate(y.position = y_ceil + dplyr::row_number() * y_step)

    dplyr::bind_rows(st_A, st_B)
  })
}

stat_p1 <- run_wilcox_stats(df_cell_means_clean, dodge = dodge_w)
stat_p2 <- run_wilcox_stats(df_long_clean,        dodge = dodge_w)

# ── shared theme ─────────────────────────────────────────────────────────────
plot_theme <- theme_minimal(base_size = 12) +
  theme(strip.placement   = "outside",
        strip.background  = element_blank(),
        strip.text        = element_text(face = "bold", angle = 90),
        legend.position   = "bottom",
        plot.subtitle     = element_text(size = 8, face = "plain"),
        plot.title        = element_text(size = 8, face = "bold"),
        panel.grid        = element_blank(),
        axis.line.y       = element_line(colour = "grey30", linewidth = 0.4),
        axis.ticks.y      = element_line(colour = "grey30", linewidth = 0.4),
        axis.ticks.length.y = unit(3, "pt"))

# ---- Plot 1: per-cell medians  (x = location, fill = age_group) ----

p1 <- ggplot(df_cell_means_clean,
             aes(x = location, y = value, fill = age_group)) +
  geom_violin(aes(colour = age_group),
              position  = position_dodge(dodge_w),
              alpha = 0.4, trim = TRUE, linewidth = 0.5) +
  geom_boxplot(width = 0.12, position = position_dodge(dodge_w),
               outlier.shape = NA, linewidth = 0.4, colour = "grey25") +
  stat_pvalue_manual(stat_p1, label = "p.adj.signif",
                     tip.length = 0.01, hide.ns = TRUE,
                     size = 3.5, bracket.size = 0.4) +
  scale_fill_manual(values  = age_pal, name = "Age group") +
  scale_colour_manual(values = age_pal, name = "Age group") +
  facet_wrap(~ metric, scales = "free", nrow = 1,
             strip.position = "left",
             labeller = labeller(metric = metric_labels)) +
  plot_theme +
  labs(title    = "Per-cell medians",
       subtitle = "Wilcoxon rank-sum, BH-corrected | *p<0.05  **p<0.01  ***p<0.001",
       x = NULL, y = NULL)

# ---- Plot 2: every peak pooled  (x = location, fill = age_group) ----

p2 <- ggplot(df_long_clean,
             aes(x = location, y = value, fill = age_group)) +
  geom_violin(position  = position_dodge(dodge_w),
              alpha = 0.4, trim = TRUE, linewidth = 0.5) +
  geom_boxplot(width = 0.12, position = position_dodge(dodge_w),
               outlier.shape = NA, linewidth = 0.4, colour = "grey25") +
  stat_pvalue_manual(stat_p2, label = "p.adj.signif",
                     tip.length = 0.01, hide.ns = TRUE,
                     size = 3.5, bracket.size = 0.4) +
  scale_fill_manual(values = age_pal, name = "Age group") +
  facet_wrap(~ metric, scales = "free", nrow = 1,
             strip.position = "left",
             labeller = labeller(metric = metric_labels)) +
  plot_theme +
  labs(title    = "All peaks pooled",
       subtitle = "\u26a0 n = individual peaks \u2014 brackets exploratory (pseudoreplication)",
       x = NULL, y = NULL)
# Combined view
p3 <- p1 / p2
p3
file_name_comb <- file.path(data_dir, "_PeakAmpl.png")
ggsave(file_name_comb, p3, width = 9, height = 8, dpi = 300)
# Save if desired
# ggsave("violins_cell_means.png", p1, width = 9, height = 4, dpi = 300)
# ggsave("violins_all_peaks.png",  p2, width = 9, height = 4, dpi = 300)
