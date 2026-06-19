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


collate_calcium_peaks <- function(path, pattern = "_synchrony_metrics\\.csv$") {
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
    set_names(\(f) sub("_?synchrony_metrics\\.csv$", "", basename(f))) |>
    map_dfr(read_csv, show_col_types = FALSE, .id = "source")
}

# --- load all data ---
data_dir <- r"(C:\Users\andre\Documents\Boaz_Ca_Analysis\Data_Jun2026)"
df_sync <- collate_calcium_peaks(data_dir)

df_sync <- df_sync |>
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

glimpse(df_sync)

# --- average synchrony metrics per recording ---

sync_means <- df_sync |>
  group_by(source, Metric, location, age_group) |>
  summarise(value = mean(Value, na.rm = TRUE), .groups = "drop") |>
  rename(metric = Metric) |>
  mutate(metric = factor(metric,
                         levels = c("MedianPearsonR", "PopSyncIndex_GolombRinzel")))

sync_labels <- c(
  MedianPearsonR            = "Median Pearson R",
  PopSyncIndex_GolombRinzel = "Pop. sync index (Golomb–Rinzel)"
)

# --- plot experiment averages ---

my_cols <- brewer.pal(8, "Dark2")[c(1, 4, 5)]

# ── Statistics ────────────────────────────────────────────────────────────────
# n = 5 recordings only — Wilcoxon tests shown where group n ≥ 2; tryCatch
# guards against groups with only 1 observation.
dodge_w <- 0.8
age_pal <- c("w14-15" = "#66C2A5", "w18-19" = "#FC8D62")

run_wilcox_stats <- function(data, dodge = 0.8) {
  purrr::map_dfr(levels(data$metric), function(m) {
    d <- dplyr::filter(data, metric == m)
    st_A <- tryCatch({
      d |>
        dplyr::group_by(location) |>
        rstatix::wilcox_test(value ~ age_group) |>
        rstatix::adjust_pvalue(method = "BH") |>
        rstatix::add_significance() |>
        rstatix::add_xy_position(x = "location", dodge = dodge,
                                 data = d, fun = "max") |>
        dplyr::mutate(metric = factor(m, levels = levels(data$metric)),
                      comparison_type = "within_x")
    }, error = function(e) NULL)
    st_B_raw <- tryCatch({
      d |>
        dplyr::group_by(age_group) |>
        rstatix::wilcox_test(value ~ location) |>
        rstatix::adjust_pvalue(method = "BH") |>
        rstatix::add_significance() |>
        rstatix::add_xy_position(x = "location", dodge = dodge,
                                 group = "age_group", data = d, fun = "max") |>
        dplyr::mutate(metric = factor(m, levels = levels(data$metric)),
                      comparison_type = "across_x")
    }, error = function(e) NULL)
    if (is.null(st_A) && is.null(st_B_raw)) return(tibble::tibble())
    if (is.null(st_B_raw)) return(st_A)
    if (is.null(st_A))     return(st_B_raw)
    y_ceil <- max(st_A$y.position, na.rm = TRUE)
    y_step <- 0.10 * diff(range(d$value, na.rm = TRUE))
    st_B   <- st_B_raw |>
      dplyr::arrange(age_group) |>
      dplyr::mutate(y.position = y_ceil + dplyr::row_number() * y_step)
    dplyr::bind_rows(st_A, st_B)
  })
}

stat_p4 <- run_wilcox_stats(sync_means, dodge = dodge_w)

plot_theme <- theme_minimal(base_size = 12) +
  theme(strip.placement     = "outside",
        strip.background    = element_blank(),
        strip.text          = element_text(face = "bold", angle = 90),
        legend.position     = "bottom",
        plot.subtitle       = element_text(size = 8, face = "plain"),
        plot.title          = element_text(size = 11),
        panel.grid          = element_blank(),
        axis.line.y         = element_line(colour = "grey30", linewidth = 0.4),
        axis.ticks.y        = element_line(colour = "grey30", linewidth = 0.4),
        axis.ticks.length.y = unit(3, "pt"))

p4 <- ggplot(sync_means, aes(x = location, y = value, fill = age_group)) +
  geom_boxplot(width = 0.35, position = position_dodge(dodge_w),
               alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(colour = age_group),
              position = position_jitterdodge(dodge.width  = dodge_w,
                                              jitter.width = 0.08),
              alpha = 0.9, size = 2.5) +
  stat_pvalue_manual(stat_p4, label = "p.adj.signif",
                     tip.length = 0.01, hide.ns = TRUE,
                     size = 3.5, bracket.size = 0.4) +
  scale_fill_manual(values  = age_pal, name = "Age group") +
  scale_colour_manual(values = age_pal, name = "Age group") +
  facet_wrap(~ metric, scales = "free", nrow = 1,
             strip.position = "left",
             labeller = labeller(metric = sync_labels)) +
  plot_theme +
  labs(title    = "Synchrony metrics (per-experiment means)",
       subtitle = paste0("\u26a0 n = ", dplyr::n_distinct(sync_means$source),
                         " recordings \u2014 interpret stat brackets with caution"),
       x = NULL, y = NULL)

# save as png
file_name_comb <- file.path(data_dir, "Sync_Metrics.png")
ggsave(p4, filename = file_name_comb, width = 8, height = 4, dpi = 300)
ggsave(file.path(data_dir, "Sync_Metrics.svg"), p4, width = 8, height = 4)
