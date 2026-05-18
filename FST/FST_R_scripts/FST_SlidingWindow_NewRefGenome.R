# =============================================================================
# Sex-biased F_ST Manhattan Plots — Crab and Wave Ecotypes
# Littorina saxatilis, Koster populations, new reference genome
#
# Weighted and mean F_ST calculated per ecotype (male vs. female) using
# VCFtools (Weir & Cockerham estimator) in 50 kb windows, 10 kb steps.
# Windows below the genome-wide median are excluded before plotting.
#
# Written by Amin Ghane
# Last update : 2026-05
# =============================================================================

library(tidyverse)
library(stringr)

setwd("/Users/aminghane/Master MER+/MyThesis/00_Koster_SexChr_WGS/Littorina_saxatilis_sex_determination/FST/FST_in_out_files")


# =============================================================================
# 0. Shared inputs
# =============================================================================

scaffold_lengths <- read_tsv(
  "Lsax_genome_CLR_HiC_curated_freeze_1_2023_02_17.fasta.fai",
  col_names = c("scaf", "length")
)

# Shared plot dimensions for ggsave (A5 landscape)
PLOT_W  <- 210    # mm
PLOT_H  <- 148.5  # mm
PLOT_DPI <- 300


# =============================================================================
# 1. Crab ecotype
# =============================================================================

# --- 1.1 Load and filter -----------------------------------------------------
# "w" = windowed / weighted
# Pre-filter row counts (before FST >= 0 filter): see lab notebook
Fst_w_saxa_crab <- read_tsv(
  "sex_saxatilis_crab_fst_vcftools_window50kb_step10kb.windowed.weir.fst"
) %>%
  replace_na(list(WEIR_AND_COCKERHAM_FST = 0)) %>%
  rename(scaf = CHROM, base = BIN_START) %>%
  mutate(base = as.numeric(base)) %>%
  filter(WEIGHTED_FST >= 0) %>%
  filter(MEAN_FST    >= 0)

# --- 1.2 Quantile cut-offs ---------------------------------------------------
quantile(Fst_w_saxa_crab$WEIGHTED_FST, c(0.50, 0.95, 0.99), na.rm = TRUE)
# 50%: 0.0105822 | 95%: 0.0554242 | 99%: 0.0893876
mean(Fst_w_saxa_crab$WEIGHTED_FST)
# mean: 0.01671537

CRAB_MEDIAN_FST <- 0.0105822   # genome-wide median (50th percentile)

# --- 1.3 Join scaffold lengths and retain windows above median ---------------
Fst_W_Saxa_crab <- left_join(scaffold_lengths, Fst_w_saxa_crab) %>%
  filter(WEIGHTED_FST >= CRAB_MEDIAN_FST) %>%
  na.omit()

# --- 1.4 Extract linkage group number from scaffold name ---------------------
Fst_W_Saxa_crab$CHR <- str_extract(Fst_W_Saxa_crab$scaf, "(?<=_)[0-9]+") %>%
  as.numeric()

# --- 1.5 Restrict to super-scaffolds (LG 1–17) and compute cumulative pos. --
Fst_W_Saxa_crab_SUPER <- subset(Fst_W_Saxa_crab, CHR <= 17)

fst_w_plot_crab <- Fst_W_Saxa_crab_SUPER %>%
  group_by(CHR) %>%
  summarise(chr_len = max(base)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(Fst_W_Saxa_crab_SUPER, ., by = "CHR") %>%
  arrange(CHR, base) %>%
  mutate(basecum = base + tot)

axisdf_w_crab <- fst_w_plot_crab %>%
  group_by(CHR) %>%
  summarize(center = (max(basecum) + min(basecum)) / 2)

# --- 1.6 Manhattan plot: weighted F_ST, all LGs ------------------------------
fst_w_crab_all <- ggplot(fst_w_plot_crab, aes(x = basecum, y = WEIGHTED_FST)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.3, size = 0.5) +
  scale_color_manual(values = rep(c("black", "#E67E22"), 22)) +
  scale_x_continuous(labels = axisdf_w_crab$CHR, breaks = axisdf_w_crab$center) +
  scale_y_continuous(expand = c(0, 0)) +
  ylim(CRAB_MEDIAN_FST, 0.3) +
  labs(title = "Weighted FST in Crab ecotype",
       x = "Super-scaffold (LG)",
       y = expression(Weighted ~ italic(F)[ST] ~ (Male:Female))) +
  theme_bw() +
  theme(
    legend.position   = "none",
    axis.title.y      = element_text(size = 10, angle = 90, hjust = .5, vjust = .5),
    plot.title        = element_text(hjust = 0.5),
    panel.border      = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

fst_w_crab_all_smooth <- fst_w_crab_all +
  geom_smooth(method = "gam", se = TRUE, level = 0.95)

# --- 1.7 Per-LG zoom: LG3 ---------------------------------------------------
fst_w_crab_CHR3 <- ggplot(
  subset(fst_w_plot_crab, CHR == 3),
  aes(x = basecum, y = WEIGHTED_FST)
) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.4, size = 0.5) +
  scale_color_manual(values = rep("#E67E22", 22)) +
  scale_x_continuous(labels = axisdf_w_crab$CHR, breaks = axisdf_w_crab$center) +
  scale_y_continuous(expand = c(0, 0)) +
  ylim(CRAB_MEDIAN_FST, 0.22) +
  labs(title = "Weighted FST in Crab ecotype — LG3",
       x = "Super-scaffold (LG)",
       y = expression(Weighted ~ italic(F)[ST] ~ (Male:Female))) +
  theme_bw() +
  theme(
    legend.position    = "none",
    axis.title.y       = element_text(size = 10, angle = 90, hjust = .5, vjust = .5),
    plot.title         = element_text(hjust = 0.5),
    panel.border       = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

fst_w_crab_CHR3_smooth <- fst_w_crab_CHR3 +
  geom_smooth(method = "gam", se = TRUE, level = 0.95)

# --- 1.8 Manhattan plot: mean F_ST, all LGs ----------------------------------
fst_mean_crab_all <- ggplot(fst_w_plot_crab, aes(x = basecum, y = MEAN_FST)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.4, size = 0.5) +
  scale_color_manual(values = rep(c("black", "#E67E22"), 22)) +
  scale_x_continuous(labels = axisdf_w_crab$CHR, breaks = axisdf_w_crab$center) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Mean FST in Crab ecotype",
       x = "Super-scaffold (LG)",
       y = expression(italic(Mean ~ F)[ST] ~ between ~ males ~ and ~ females ~ (Crab))) +
  theme_bw() +
  theme(
    legend.position    = "none",
    axis.title.y       = element_text(size = 10, angle = 90, hjust = .5, vjust = .5),
    plot.title         = element_text(hjust = 0.5),
    panel.border       = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

fst_mean_crab_all_smooth <- fst_mean_crab_all +
  geom_smooth(method = "gam", se = TRUE, level = 0.95)


# =============================================================================
# 2. Wave ecotype
# =============================================================================

# --- 2.1 Load and filter -----------------------------------------------------
# Pre-filter row counts: 112,025 total; 40,787 after FST >= 0 filter
Fst_w_saxa_wave <- read_tsv(
  "sex_saxatilis_wave_fst_vcftools_window50kb_step10kb.windowed.weir.fst"
) %>%
  replace_na(list(WEIR_AND_COCKERHAM_FST = 0)) %>%
  rename(scaf = CHROM, base = BIN_START) %>%
  mutate(base = as.numeric(base)) %>%
  filter(WEIGHTED_FST >= 0) %>%
  filter(MEAN_FST    >= 0)

# --- 2.2 Quantile cut-offs ---------------------------------------------------
quantile(Fst_w_saxa_wave$WEIGHTED_FST, c(0.50, 0.95, 0.99), na.rm = TRUE)
# 50%: 0.01221590 | 95%: 0.04832190 | 99%: 0.08109203
mean(Fst_w_saxa_wave$WEIGHTED_FST)
# mean: 0.01701703

WAVE_MEDIAN_FST <- 0.01221590   # genome-wide median (50th percentile)

# --- 2.3 Join scaffold lengths and retain windows above median ---------------
Fst_W_Saxa_wave <- left_join(scaffold_lengths, Fst_w_saxa_wave) %>%
  filter(WEIGHTED_FST >= WAVE_MEDIAN_FST) %>%
  na.omit()

# --- 2.4 Extract linkage group number ----------------------------------------
Fst_W_Saxa_wave$CHR <- str_extract(Fst_W_Saxa_wave$scaf, "(?<=_)[0-9]+") %>%
  as.numeric()

# --- 2.5 Restrict to super-scaffolds and compute cumulative position ---------
Fst_W_Saxa_wave_SUPER <- subset(Fst_W_Saxa_wave, CHR <= 17)

fst_w_plot_wave <- Fst_W_Saxa_wave_SUPER %>%
  group_by(CHR) %>%
  summarise(chr_len = max(base)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(Fst_W_Saxa_wave_SUPER, ., by = "CHR") %>%
  arrange(CHR, base) %>%
  mutate(basecum = base + tot)

axisdf_w_wave <- fst_w_plot_wave %>%
  group_by(CHR) %>%
  summarize(center = (max(basecum) + min(basecum)) / 2)

# --- 2.6 Manhattan plot: weighted F_ST, all LGs ------------------------------
fst_w_wave_all <- ggplot(fst_w_plot_wave, aes(x = basecum, y = WEIGHTED_FST)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.3, size = 0.5) +
  scale_color_manual(values = rep(c("black", "#007fff"), 22)) +
  scale_x_continuous(labels = axisdf_w_wave$CHR, breaks = axisdf_w_wave$center) +
  scale_y_continuous(expand = c(0, 0)) +
  ylim(WAVE_MEDIAN_FST, 0.32) +
  labs(title = "Weighted FST in Wave ecotype",
       x = "Super-scaffold (LG)",
       y = expression(Weighted ~ italic(F)[ST] ~ (Male:Female))) +
  theme_bw() +
  theme(
    legend.position    = "none",
    axis.title.y       = element_text(size = 10, angle = 90, hjust = .5, vjust = .5),
    plot.title         = element_text(hjust = 0.5),
    panel.border       = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

fst_w_wave_all_smooth <- fst_w_wave_all +
  geom_smooth(color = "#98FB98", method = "gam", se = TRUE, level = 0.95)

# --- 2.7 Per-LG zooms: LG9 and LG8 ------------------------------------------
make_wave_lg_plot <- function(lg_num) {
  ggplot(subset(fst_w_plot_wave, CHR == lg_num),
         aes(x = basecum, y = WEIGHTED_FST)) +
    geom_point(aes(color = as.factor(CHR)), alpha = 0.4, size = 0.5) +
    scale_color_manual(values = rep("#007fff", 22)) +
    scale_x_continuous(labels = axisdf_w_wave$CHR, breaks = axisdf_w_wave$center) +
    scale_y_continuous(expand = c(0, 0)) +
    ylim(WAVE_MEDIAN_FST, 0.32) +
    labs(title = paste0("Weighted FST in Wave ecotype — LG", lg_num),
         x = "Super-scaffold (LG)",
         y = expression(Weighted ~ italic(F)[ST] ~ (Male:Female))) +
    theme_bw() +
    theme(
      legend.position    = "none",
      axis.title.y       = element_text(size = 10, angle = 90, hjust = .5, vjust = .5),
      plot.title         = element_text(hjust = 0.5),
      panel.border       = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

fst_w_wave_CHR9 <- make_wave_lg_plot(9)
fst_w_wave_CHR8 <- make_wave_lg_plot(8)

fst_w_wave_CHR9_smooth <- fst_w_wave_CHR9 +
  geom_smooth(color = "#98FB98", method = "gam", se = TRUE, level = 0.95)
fst_w_wave_CHR8_smooth <- fst_w_wave_CHR8 +
  geom_smooth(color = "#98FB98", method = "gam", se = TRUE, level = 0.95)


# =============================================================================
# 3. Save all plots
# =============================================================================

plots_to_save <- list(
  fst_weighted_plot_crab_all         = fst_w_crab_all,
  fst_weighted_plot_crab_all_smooth  = fst_w_crab_all_smooth,
  fst_weighted_plot_crab_CHR3_smooth = fst_w_crab_CHR3_smooth,
  fst_mean_plot_crab_all             = fst_mean_crab_all,
  fst_weighted_plot_wave_all         = fst_w_wave_all,
  fst_weighted_plot_wave_all_smooth  = fst_w_wave_all_smooth,
  fst_weighted_plot_wave_CHR9_smooth = fst_w_wave_CHR9_smooth,
  fst_weighted_plot_wave_CHR8_smooth = fst_w_wave_CHR8_smooth
)

out_dir <- "/Users/aminghane/Master MER+/MyThesis/00_Koster_SexChr_WGS/Figures/FST"

for (fname in names(plots_to_save)) {
  ggsave(
    filename = paste0(fname, ".png"),
    plot     = plots_to_save[[fname]],
    path     = out_dir,
    device   = png,
    width    = PLOT_W,
    height   = PLOT_H,
    units    = "mm",
    dpi      = PLOT_DPI
  )
}