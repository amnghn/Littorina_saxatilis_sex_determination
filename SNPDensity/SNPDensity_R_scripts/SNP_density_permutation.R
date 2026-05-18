### Amin Ghane; ISTA
### Littorina saxatilis sex determination project
### SNP density analysis — Crab and Wave ecotypes, new reference genome
### 2023-07-05

### ============================================================
### Required packages
### ============================================================
library(tidyverse)
library(ggpubr)
library(ggthemes)

### ============================================================
### Paths — adjust as needed
### ============================================================
WORK_DIR   <- "/Users/aminghane/Master MER+/MyThesis/00_Koster_SexChr_WGS/Littorina_saxatilis_sex_determination/SNPDensity/SNPDensity_in_out_files"
FIG_DIR    <- "/Users/aminghane/Master MER+/MyThesis/00_Koster_SexChr_WGS/Figures/SNPDensity"
GENOME_FAI <- "/Users/aminghane/Master MER+/MyThesis/00_Koster_SexChr_WGS/Littorina_saxatilis_sex_determination/SNPDensity/SNPDensity_in_out_files/Lsax_genome_CLR_HiC_curated_freeze_1_2023_02_17.fasta.fai"

### ============================================================
### Helper functions
### ============================================================

# Read all per-individual .snpden files from a directory into one wide tibble.
# Files must be named {Sex}_{ecotype}_{SampleID}_*.vcf.snpden
# Returns a tibble with columns: scaf, base, <SampleID_1>, <SampleID_2>, ...
read_snpden_files <- function(dir, pattern) {
  files <- list.files(path = dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files found matching pattern: ", pattern, " in ", dir)

  # Build backbone from first file
  first      <- read_delim(files[1], delim = "\t", col_names = TRUE, show_col_types = FALSE)
  parts      <- strsplit(basename(files[1]), "_")[[1]]
  col_id     <- paste(parts[1], parts[3], sep = "_")
  backbone   <- first %>%
    unite(LOCATION, c("CHROM", "BIN_START"), sep = ":") %>%
    dplyr::rename(!!col_id := "VARIANTS/KB") %>%
    select(-SNP_COUNT)

  # Loop remaining files
  for (i in seq(2, length(files))) {
    dat    <- read_delim(files[i], delim = "\t", col_names = TRUE, show_col_types = FALSE)
    parts  <- strsplit(basename(files[i]), "_")[[1]]
    col_id <- paste(parts[1], parts[3], sep = "_")
    dat    <- dat %>%
      unite(LOCATION, c("CHROM", "BIN_START"), sep = ":") %>%
      dplyr::rename(!!col_id := "VARIANTS/KB") %>%
      select(-SNP_COUNT)
    backbone <- full_join(backbone, dat, by = "LOCATION")
  }

  backbone %>% separate(LOCATION, c("scaf", "base"), sep = ":")
}

# Calculate per-window mean SNP density for males and females, and their difference.
# Returns a tibble with: scaf, base, mean_Males, mean_Females, mean_MvF_dif
calc_true_snp_density <- function(snpden_males, snpden_females) {
  snpden <- full_join(snpden_males, snpden_females, by = c("scaf", "base"))
  snpden$base <- as.numeric(snpden$base)
  snpden <- snpden %>%
    replace_na(list(Males = 0, Females = 0)) %>%
    replace(is.na(.), 0)

  male_cols   <- snpden[, grepl("Male",   names(snpden))]
  female_cols <- snpden[, grepl("Female", names(snpden))]

  males_mean   <- bind_cols(snpden %>% select(scaf, base), male_cols   %>% mutate(mean_Males   = rowMeans(.))) %>% select(scaf, base, mean_Males)
  females_mean <- bind_cols(snpden %>% select(scaf, base), female_cols %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf, base, mean_Females)

  list(
    snpden      = snpden,
    male_n      = ncol(male_cols),
    female_n    = ncol(female_cols),
    true_density = full_join(males_mean, females_mean, by = c("scaf", "base")) %>%
      mutate(mean_MvF_dif = mean_Males - mean_Females)
  )
}

# Run 1000-iteration permutation test on SNP density data.
# Shuffles individual columns, recalculates male/female means, stores per-permutation MvF diff.
# Returns perm_backbone (wide tibble with p1..p1000 columns).
run_permutation <- function(snpden, male_n, female_n, n_perm = 1000) {
  run_one_perm <- function(snpden, male_n, female_n) {
    perm_idx      <- sample(3:ncol(snpden))
    snpden_perm   <- snpden %>% select(1:2, all_of(perm_idx))
    males_perm    <- bind_cols(
      snpden_perm %>% select(1:2),
      snpden_perm %>% select(3:(male_n + 2)) %>% mutate(mean_Males = rowMeans(.))
    ) %>% select(scaf, base, mean_Males)
    females_perm  <- bind_cols(
      snpden_perm %>% select(1:2),
      snpden_perm %>% select((male_n + 3):(male_n + 2 + female_n)) %>% mutate(mean_Females = rowMeans(.))
    ) %>% select(scaf, base, mean_Females)
    full_join(males_perm, females_perm, by = c("scaf", "base")) %>%
      mutate(mean_MvF_dif = mean_Males - mean_Females) %>%
      select(scaf, base, mean_MvF_dif)
  }

  perm_backbone <- run_one_perm(snpden, male_n, female_n) %>% rename(p1 = mean_MvF_dif)

  for (i in 2:n_perm) {
    col_id        <- paste0("p", i)
    perm_upgrade  <- run_one_perm(snpden, male_n, female_n) %>%
      dplyr::rename(!!col_id := "mean_MvF_dif")
    perm_backbone <- full_join(perm_backbone, perm_upgrade, by = c("scaf", "base"))
  }
  perm_backbone
}

# Calculate empirical p-values: proportion of permutation MvF diffs more extreme than observed.
# Uses a one-sided test conditional on the direction of the observed difference.
calc_pvalues <- function(true_density, perm_backbone) {
  perm_with_true <- full_join(
    true_density %>% select(scaf, base, mean_MvF_dif),
    perm_backbone,
    by = c("scaf", "base")
  )
  perm_cols <- perm_with_true[, grep("^p", names(perm_with_true))]
  n_perm    <- ncol(perm_cols)

  perm_with_true %>%
    mutate(Pvalue = case_when(
      mean_MvF_dif > 0  ~ (rowSums(perm_cols > mean_MvF_dif) + 1) / (n_perm + 1),
      mean_MvF_dif < 0  ~ (rowSums(perm_cols < mean_MvF_dif) + 1) / (n_perm + 1),
      mean_MvF_dif == 0 ~ 1
    ))
}

# Extract the chromosome number from scaffold names (e.g. "SUPER_1" → 1).
add_chr_number <- function(df) {
  df %>% mutate(CHR = as.numeric(str_extract(scaf, "(?<=_)[0-9]+")))
}

# Manhattan-style plot of mean MvF SNP density difference or p-value across super-scaffolds.
plot_snpden <- function(data_super, y_var, y_label, title, point_color, y_limits = NULL) {
  p <- ggplot(data_super, aes(x = base, y = .data[[y_var]])) +
    geom_point(aes(color = as.factor(CHR)), alpha = 1, size = 0.7) +
    scale_x_discrete(name = "SUPER scaffolds", breaks = data_super$CHR) +
    scale_color_manual(values = rep(c("black", point_color), 22)) +
    facet_grid(cols = vars(CHR), scales = "free_x", space = "free_x", switch = "x") +
    labs(title = title) +
    theme_igray() +
    theme(
      axis.text.x        = element_text(size = 8),
      axis.text.y        = element_text(size = 8),
      axis.title.y       = element_text(size = 10, angle = 90, hjust = .5, vjust = .5),
      plot.title         = element_text(hjust = 0.5),
      legend.position    = "none",
      panel.border       = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    ylab(y_label)
  if (!is.null(y_limits)) p <- p + ylim(y_limits)
  p
}

### ============================================================
### Genome index (scaffold lengths)
### ============================================================
saxatilis_index <- read_tsv(GENOME_FAI, col_names = FALSE, show_col_types = FALSE) %>%
  rename(scaf = X1, length = X2)

### ============================================================
### CRAB ecotype
### ============================================================
cat("--- Processing Crab ecotype ---\n")
crab_dir <- file.path(WORK_DIR, "individual_SNP_density_crab")

snpden_males_crab   <- read_snpden_files(crab_dir, "Male_crab")
snpden_females_crab <- read_snpden_files(crab_dir, "Female_crab")

crab <- calc_true_snp_density(snpden_males_crab, snpden_females_crab)
cat("  Crab male_n =", crab$male_n, "| female_n =", crab$female_n, "\n")
cat("  Crab MvF dif > 0:", sum(crab$true_density$mean_MvF_dif > 0),
    "| < 0:", sum(crab$true_density$mean_MvF_dif < 0),
    "| == 0:", sum(crab$true_density$mean_MvF_dif == 0), "\n")

## -- Run permutation (or load pre-computed result) --
# Uncomment the block below to re-run the permutation test (~1000 iterations, slow):
# perm_backbone_crab     <- run_permutation(crab$snpden, crab$male_n, crab$female_n)
# perm_with_true_p_crab  <- calc_pvalues(crab$true_density, perm_backbone_crab)
# write_tsv(perm_with_true_p_crab, file.path(crab_dir, "SNPdensity_perm_with_true_p_saxatilis_crab.txt"))
# write_tsv(perm_with_true_p_crab %>% select(scaf, base, mean_MvF_dif, Pvalue),
#           file.path(crab_dir, "SNPdensity_SexFindR_saxatilis_crab.txt"))

## Load pre-computed permutation results:
perm_with_true_p_crab <- read_delim(
  file.path(crab_dir, "SNPdensity_SexFindR_saxatilis_crab.txt"),
  show_col_types = FALSE
)

## Subset to p <= 0.01 and annotate
perm_p01_crab <- perm_with_true_p_crab %>%
  filter(Pvalue <= 0.01) %>%
  add_chr_number()

proportion_p01_crab <- perm_p01_crab %>%
  count(scaf) %>%
  left_join(saxatilis_index, by = "scaf") %>%
  mutate(proportion = (n * 10000) / length) %>%
  add_chr_number()

## Separate super-scaffolds (CHR 1–17) from small scaffolds
perm_p01_crab_SUPER <- filter(perm_p01_crab, CHR <= 17)
perm_p01_crab_other <- filter(perm_p01_crab, CHR > 17) %>%
  mutate(CHR = "small_scaffolds")

### ============================================================
### WAVE ecotype
### ============================================================
cat("--- Processing Wave ecotype ---\n")
wave_dir <- file.path(WORK_DIR, "individual_SNP_density_wave")

snpden_males_wave   <- read_snpden_files(wave_dir, "Male_wave")
snpden_females_wave <- read_snpden_files(wave_dir, "Female_wave")

wave <- calc_true_snp_density(snpden_males_wave, snpden_females_wave)
cat("  Wave male_n =", wave$male_n, "| female_n =", wave$female_n, "\n")
cat("  Wave MvF dif > 0:", sum(wave$true_density$mean_MvF_dif > 0),
    "| < 0:", sum(wave$true_density$mean_MvF_dif < 0),
    "| == 0:", sum(wave$true_density$mean_MvF_dif == 0), "\n")

## -- Run permutation (or load pre-computed result) --
# Uncomment the block below to re-run the permutation test (~1000 iterations, slow):
# perm_backbone_wave     <- run_permutation(wave$snpden, wave$male_n, wave$female_n)
# perm_with_true_p_wave  <- calc_pvalues(wave$true_density, perm_backbone_wave)
# write_tsv(perm_with_true_p_wave, file.path(wave_dir, "SNPdensity_perm_with_true_p_saxatilis_wave.txt"))
# write_tsv(perm_with_true_p_wave %>% select(scaf, base, mean_MvF_dif, Pvalue),
#           file.path(wave_dir, "SNPdensity_SexFindR_saxatilis_wave.txt"))

## Load pre-computed permutation results:
perm_with_true_p_wave <- read_delim(
  file.path(wave_dir, "SNPdensity_SexFindR_saxatilis_wave.txt"),
  show_col_types = FALSE
)

## Subset to p <= 0.01 and annotate
perm_p01_wave <- perm_with_true_p_wave %>%
  filter(Pvalue <= 0.01) %>%
  add_chr_number()

proportion_p01_wave <- perm_p01_wave %>%
  count(scaf) %>%
  left_join(saxatilis_index, by = "scaf") %>%
  mutate(proportion = (n * 10000) / length) %>%
  add_chr_number()

## Separate super-scaffolds from small scaffolds
perm_p01_wave_SUPER <- filter(perm_p01_wave, CHR <= 17)
perm_p01_wave_other <- filter(perm_p01_wave, CHR > 17) %>%
  mutate(CHR = "small_scaffolds")

### ============================================================
### Plots
### ============================================================

## Crab — mean MvF difference
snpden_crab_plot_Mdif <- plot_snpden(
  perm_p01_crab_SUPER,
  y_var     = "mean_MvF_dif",
  y_label   = "10kb SNP density (Male mean − Female mean)",
  title     = "SNP density in Crab ecotype",
  point_color = "#E67E22",
  y_limits  = c(-4.09, 4.5)
)

## Crab — p-value
snpden_crab_plot_Pvalue <- plot_snpden(
  perm_p01_crab_SUPER,
  y_var     = "Pvalue",
  y_label   = "P-value",
  title     = "SNP density in Crab ecotype",
  point_color = "#E67E22"
)

## Wave — mean MvF difference
snpden_wave_plot_Mdif <- plot_snpden(
  perm_p01_wave_SUPER,
  y_var     = "mean_MvF_dif",
  y_label   = "10kb SNP density (Male mean − Female mean)",
  title     = "SNP density in Wave ecotype",
  point_color = "#007fff",
  y_limits  = c(-4.09, 4.5)
)

## Wave — p-value
snpden_wave_plot_Pvalue <- plot_snpden(
  perm_p01_wave_SUPER,
  y_var     = "Pvalue",
  y_label   = "P-value",
  title     = "SNP density in Wave ecotype",
  point_color = "#007fff"
)

### ============================================================
### Save figures
### ============================================================
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

ggsave(snpden_crab_plot_Mdif,
       filename = "snpden_plot_crab_MvF_dif.png",
       path     = FIG_DIR,
       device   = png,
       width    = 210,
       height   = 148.5,
       units    = "mm",
       dpi      = 700)

ggsave(snpden_wave_plot_Mdif,
       filename = "snpden_plot_wave_MvF_dif.png",
       path     = FIG_DIR,
       device   = png,
       width    = 210,
       height   = 148.5,
       units    = "mm",
       dpi      = 700)

ggsave(snpden_crab_plot_Pvalue,
       filename = "snpden_plot_crab_Pvalue.png",
       path     = FIG_DIR,
       device   = png,
       width    = 210,
       height   = 148.5,
       units    = "mm",
       dpi      = 700)

ggsave(snpden_wave_plot_Pvalue,
       filename = "snpden_plot_wave_Pvalue.png",
       path     = FIG_DIR,
       device   = png,
       width    = 210,
       height   = 148.5,
       units    = "mm",
       dpi      = 700)

cat("Figures saved to:", FIG_DIR, "\n")