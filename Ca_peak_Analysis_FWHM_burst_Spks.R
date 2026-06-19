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
library(RColorBrewer)
library(rstatix)
library(ggpubr)


collate_calcium_peaks <- function(path, pattern = "_burst_events\\.csv$") {
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
    into  = c("age", "location", NA),
    sep   = "_",
    extra = "drop",
    remove = FALSE
  ) |>
  mutate(
    age_num   = as.numeric(str_extract(age, "[0-9]+")),
    age_group = case_when(
      age_num <= 15 ~ "w14-15",
      age_num >= 18 ~ "w18-19",
      .default = NA_character_
    ),
    location  = factor(location,  levels = c("apex", "base")),
    age_group = factor(age_group, levels = c("w14-15", "w18-19"))
  )

# --- outlier helper (Tukey 1.5x IQR) — defined here so it's available everywhere ---
filter_iqr <- function(x, k = 1.5) {
  q   <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- diff(q)
  x >= q[1] - k * iqr & x <= q[2] + k * iqr
}

glimpse(df)

burst_stats <- df |>
  group_by(source) |>
  mutate(recording_duration_s = max(BurstTime_s, na.rm = TRUE)) |>
  group_by(source, Cell, location, age_group, recording_duration_s) |>
  summarise(n_bursts = n(), .groups = "drop") |>
  mutate(burst_frequency_Hz = n_bursts / recording_duration_s)

glimpse(burst_stats)

#### plot

my_cols <- brewer.pal(8, "Dark2")[c(3, 4, 5)]

# plot 0: 

metrics_burst <- c(
  # "n_bursts", 
  "burst_frequency_Hz")

burst_long <- burst_stats |>
  select(source, Cell, location, age_group, all_of(metrics_burst)) |>
  pivot_longer(all_of(metrics_burst), names_to = "metric", values_to = "value") |>
  mutate(
    metric    = factor(metric,    levels = metrics_burst),
    location  = factor(location,  levels = c("apex", "base")),
    age_group = factor(age_group, levels = c("w14-15", "w18-19"))
  )

# Same IQR cleaner as before
burst_long_clean <- burst_long |>
  group_by(metric) |>
  filter(filter_iqr(value)) |>
  ungroup()

burst_metric_labels <- c(
  # n_bursts           = "Bursts per cell (count)",
  burst_frequency_Hz = "Burst frequency (Hz)"
)

p0 <- ggplot(burst_long_clean, aes(x = metric, y = value, fill = metric)) +
  geom_violinhalf(alpha = 0.5, trim = TRUE) +
  geom_jitter(aes(color = metric), width = 0.08, alpha = 0.7, size = 2) +
  facet_wrap(~ metric, scales = "free", nrow = 1,
             strip.position = "left",
             labeller = labeller(metric = burst_metric_labels)) +
  scale_fill_manual(values = my_cols[1], guide = "none") +
  scale_color_manual(values = my_cols[1], guide = "none") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        strip.placement  = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", angle = 90),
        plot.subtitle = element_text(size = 8, face = "plain"),
        plot.title     = element_text(size = 8, face = "bold")) +
  labs(title = "Burst activity per cell",
       subtitle = paste0("n = ", nrow(burst_stats),
                         " cells across ", dplyr::n_distinct(burst_stats$source),
                         " recordings"),
       x = NULL, y = NULL)

# metrics <- c("SpksAmplitude", "FWHM_s")
metrics <- c("FWHM_s") # SpksAmplitude already calc from raw calcium traces

df_long <- df |>
  select(source, Cell, location, age_group, all_of(metrics)) |>
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") |>
  mutate(
    metric    = factor(metric,    levels = metrics),
    location  = factor(location,  levels = c("apex", "base")),
    age_group = factor(age_group, levels = c("w14-15", "w18-19"))
  )

fwhm_metric_labels <- c(FWHM_s = "Burst FWHM (s)")

# Clean raw bursts, per metric
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

# Per-cell medians from cleaned data
df_cell_means_clean <- df_long_clean |>
  group_by(source, Cell, metric, age_group, location) |>
  summarise(value = median(value, na.rm = TRUE), .groups = "drop")

# ── Statistics ────────────────────────────────────────────────────────────────
dodge_w <- 0.8
age_pal <- c("w14-15" = "#66C2A5", "w18-19" = "#FC8D62")

run_wilcox_stats <- function(data, dodge = 0.8) {
  purrr::map_dfr(levels(data$metric), function(m) {
    d <- dplyr::filter(data, metric == m)
    st_A <- d |>
      dplyr::group_by(location) |>
      rstatix::wilcox_test(value ~ age_group) |>
      rstatix::adjust_pvalue(method = "BH") |>
      rstatix::add_significance() |>
      rstatix::add_xy_position(x = "location", dodge = dodge,
                               data = d, fun = "max") |>
      dplyr::mutate(metric = factor(m, levels = levels(data$metric)),
                    comparison_type = "within_x")
    st_B_raw <- d |>
      dplyr::group_by(age_group) |>
      rstatix::wilcox_test(value ~ location) |>
      rstatix::adjust_pvalue(method = "BH") |>
      rstatix::add_significance() |>
      rstatix::add_xy_position(x = "location", dodge = dodge,
                               group = "age_group", data = d, fun = "max") |>
      dplyr::mutate(metric = factor(m, levels = levels(data$metric)),
                    comparison_type = "across_x")
    y_ceil <- max(st_A$y.position, na.rm = TRUE)
    y_step <- 0.10 * diff(range(d$value, na.rm = TRUE))
    st_B   <- st_B_raw |>
      dplyr::arrange(age_group) |>
      dplyr::mutate(y.position = y_ceil + dplyr::row_number() * y_step)
    dplyr::bind_rows(st_A, st_B)
  })
}

stat_p0 <- run_wilcox_stats(burst_long_clean,   dodge = dodge_w)
stat_p1 <- run_wilcox_stats(df_cell_means_clean, dodge = dodge_w)
stat_p2 <- run_wilcox_stats(df_long_clean,        dodge = dodge_w)

# ── shared theme ─────────────────────────────────────────────────────────────
plot_theme <- theme_minimal(base_size = 12) +
  theme(strip.placement     = "outside",
        strip.background    = element_blank(),
        strip.text          = element_text(face = "bold", angle = 90),
        legend.position     = "bottom",
        plot.subtitle       = element_text(size = 8, face = "plain"),
        plot.title          = element_text(size = 8, face = "bold"),
        panel.grid          = element_blank(),
        axis.line.y         = element_line(colour = "grey30", linewidth = 0.4),
        axis.ticks.y        = element_line(colour = "grey30", linewidth = 0.4),
        axis.ticks.length.y = unit(3, "pt"))

# ---- p0: burst frequency per cell ----
p0 <- ggplot(burst_long_clean,
             aes(x = location, y = value, fill = age_group)) +
  geom_violin(aes(colour = age_group),
              position = position_dodge(dodge_w),
              alpha = 0.4, trim = TRUE, linewidth = 0.5) +
  geom_boxplot(width = 0.12, position = position_dodge(dodge_w),
               outlier.shape = NA, linewidth = 0.4, colour = "grey25") +
  stat_pvalue_manual(stat_p0, label = "p.adj.signif",
                     tip.length = 0.01, hide.ns = TRUE,
                     size = 3.5, bracket.size = 0.4) +
  scale_fill_manual(values  = age_pal, name = "Age group") +
  scale_colour_manual(values = age_pal, name = "Age group") +
  facet_wrap(~ metric, scales = "free", nrow = 1,
             strip.position = "left",
             labeller = labeller(metric = burst_metric_labels)) +
  plot_theme +
  labs(title    = "Burst frequency per cell",
       subtitle = "Wilcoxon rank-sum, BH-corrected | *p<0.05  **p<0.01  ***p<0.001",
       x = NULL, y = NULL)

# ---- p1: per-cell medians (FWHM) ----
p1 <- ggplot(df_cell_means_clean,
             aes(x = location, y = value, fill = age_group)) +
  geom_violin(aes(colour = age_group),
              position = position_dodge(dodge_w),
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
             labeller = labeller(metric = fwhm_metric_labels)) +
  plot_theme +
  labs(title    = "Per-cell medians (FWHM)",
       subtitle = "Wilcoxon rank-sum, BH-corrected | *p<0.05  **p<0.01  ***p<0.001",
       x = NULL, y = NULL)

# ---- p2: all bursts pooled (FWHM) ----
p2 <- ggplot(df_long_clean,
             aes(x = location, y = value, fill = age_group)) +
  geom_violin(position = position_dodge(dodge_w),
              alpha = 0.4, trim = TRUE, linewidth = 0.5) +
  geom_boxplot(width = 0.12, position = position_dodge(dodge_w),
               outlier.shape = NA, linewidth = 0.4, colour = "grey25") +
  stat_pvalue_manual(stat_p2, label = "p.adj.signif",
                     tip.length = 0.01, hide.ns = TRUE,
                     size = 3.5, bracket.size = 0.4) +
  scale_fill_manual(values = age_pal, name = "Age group") +
  facet_wrap(~ metric, scales = "free", nrow = 1,
             strip.position = "left",
             labeller = labeller(metric = fwhm_metric_labels)) +
  plot_theme +
  labs(title    = "All bursts pooled (FWHM)",
       subtitle = "\u26a0 n = individual bursts \u2014 brackets exploratory (pseudoreplication)",
       x = NULL, y = NULL)

# Combined view
p3 <- (p0 | p1 | p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p3
file_name_comb <- file.path(data_dir, "violins_combined_burst.png")
ggsave(file_name_comb, p3, width = 14, height = 6, dpi = 300)

