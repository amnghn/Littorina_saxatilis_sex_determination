# Sex-associated genomic variation in *Littorina saxatilis* — Koster analyses

Bioinformatics code for identifying sex-linked genomic regions in the Crab and Wave
ecotypes of *Littorina saxatilis* collected from the Koster archipelago, southwest Sweden.
Analyses are part of a collaborative manuscript on sex determination in *L. saxatilis*
(Ghane et al., *in prep*), building on Hearn et al. (2022).

---

## Table of contents

- [Sex-associated genomic variation in *Littorina saxatilis* — Koster analyses](#sex-associated-genomic-variation-in-littorina-saxatilis--koster-analyses)
  - [Table of contents](#table-of-contents)
  - [Background](#background)
  - [Repository structure](#repository-structure)
  - [Data](#data)
  - [Methods overview](#methods-overview)
    - [1. F\_ST (fixation index)](#1-f_st-fixation-index)
    - [2. SNP density](#2-snp-density)
    - [3. GWAS](#3-gwas)
    - [4. Genomic coverage](#4-genomic-coverage)
  - [Software requirements](#software-requirements)
  - [Usage](#usage)
    - [F\_ST (SLURM cluster)](#f_st-slurm-cluster)
    - [F\_ST plots (R)](#f_st-plots-r)
  - [Output files](#output-files)
  - [Citation](#citation)
  - [Contact](#contact)

---

## Background

*Littorina saxatilis* has a female heterogametic (ZW/ZZ) sex determination system in the
Crab ecotype, associated with chromosomal inversions on linkage group 12 (Hearn et al.,
2022; Koch et al., 2021). Whether a comparable system exists in the Wave ecotype remains
unclear. To investigate sex-linked genomic regions in both ecotypes, we applied the
SexFindR workflow (Grayson et al., 2022), which combines four complementary approaches:
genomic coverage, SNP density, F_ST, and GWAS.

Analyses were performed on 98 whole-genome-sequenced individuals (Crab: 15M / 33F;
Wave: 15M / 19F; mixed: 30M / 52F) mapped to a new *L. saxatilis* reference genome v2.0
(De Jode, unpublished; N50 = 298,711 bp; ~1.35 Gbp).

---

## Repository structure

```
.
├── 01_FST/
│   ├── FST_window.sh            # SLURM: sliding-window F_ST with VCFtools
│   └── FST_plot_crab_wave.R     # Manhattan plots (weighted & mean F_ST)
│
├── 02_SNP_density/
│   └── ...                      # SNP density scripts (VCFtools + SexFindR R scripts)
│
├── 03_GWAS/
│   └── ...                      # PLINK conversion + GEMMA association analysis
│
├── 04_coverage/
│   └── ...                      # DifCover male–female coverage comparison
│
├── input/
│   ├── snpDen_males_crab.txt    # Sample list: Crab males
│   ├── snpDen_females_crab.txt  # Sample list: Crab females
│   ├── snpDen_males_wave.txt    # Sample list: Wave males
│   └── snpDen_females_wave.txt  # Sample list: Wave females
│
└── README.md
```

> **Note:** The filtered VCF (`filtered_with_bcftools_7_3.vcf.gz`) and the genome index
> (`.fai`) are not hosted here due to file size. See [Data](#data) for access details.

---

## Data

| File | Description | Access |
|------|-------------|--------|
| `filtered_with_bcftools_7_3.vcf.gz` | Filtered biallelic SNPs (1,910,753 SNPs; 98 individuals) | Available on request |
| `Lsax_genome_CLR_HiC_curated_freeze_1_2023_02_17.fasta` | *L. saxatilis* reference genome v2.0 | A. De Jode, pers. comm. |
| `Lsax_genome_...fasta.fai` | FASTA index (samtools faidx) | Available on request |

**Variant filtering criteria applied (VCFtools):** biallelic sites only; no indels;
MAF > 0.05; missing data < 0.1; quality score ≥ 30; mean depth 5–25×; per-genotype
depth 4–25×. Final dataset: 1,910,753 SNPs.

---

## Methods overview

### 1. F_ST (fixation index)

Male–female F_ST was calculated per ecotype using VCFtools (v0.1.16) with the
Weir & Cockerham (1984) estimator. A sliding-window approach was used
(window = 50 kb, step = 10 kb). Windows with F_ST below the genome-wide median
were excluded before plotting. Results are visualised as Manhattan plots across the
17 super-scaffolds (linkage groups).

- **Script:** `01_FST/FST_window.sh` (SLURM), `01_FST/FST_plot_crab_wave.R`
- **Key parameters:** `--fst-window-size 50000`, `--fst-window-step 10000`
- **Cut-offs used:** Crab median = 0.0106; Wave median = 0.0122

### 2. SNP density

Per-sample SNP density was calculated in non-overlapping 10 kb windows (VCFtools
`--SNPdensity`). Mean densities for males and females were compared, and a 1,000-
permutation test (sex labels shuffled) was used to calculate window-level p-values.
Windows with p ≤ 0.01 were retained as candidates.

- **Scripts:** `02_SNP_density/`
- **Reference:** Grayson et al. (2022)

### 3. GWAS

Sex was treated as a binary phenotype (case = female, control = male). VCFtools was
used to convert the VCF to PLINK format, and GEMMA (v0.98.1; default parameters) was
used for association testing. Only sites with likelihood ratio test p-value < 0.05 were
retained for Manhattan plotting (qqman R package).

- **Scripts:** `03_GWAS/`

### 4. Genomic coverage

DifCover was used to compare male–female read coverage across the genome in 20
male–female pairs (10 Crab, 10 Wave). Coverage windows were parameterised with
target valid bases = 100 kb, min/max coverage = 1/20×, and enrichment score threshold
P = 0.74. Coverage was not analysed separately per ecotype.

- **Scripts:** `04_coverage/`

---

## Software requirements

| Software | Version | Reference |
|----------|---------|-----------|
| VCFtools | 0.1.16 | Danecek et al. (2011) |
| BCFtools | 1.14 | Danecek et al. (2021) |
| GEMMA | 0.98.1 | Zhou & Stephens (2012) |
| PLINK | 1.09b | Purcell et al. (2007) |
| DifCover | — | Symonová et al. |
| R | ≥ 4.0 | R Core Team |
| tidyverse | — | Wickham et al. |
| qqman | — | Turner (2018) |
| mgcv | — | Wood (2017) |

R package versions used in this analysis are recorded in `session_info.txt`
(generated with `sessionInfo()`).

---

## Usage

### F_ST (SLURM cluster)

```bash
sbatch 01_FST/FST_window.sh
```

The script expects the VCF and population list files to be present in the working
directory (`/proj/snic2022-6-238/Amin/FST_newRef/`). Adjust the `#SBATCH` directives
and paths at the top of the script for your own cluster environment.

### F_ST plots (R)

```r
# Set working directory to the folder containing the VCFtools output files
setwd("path/to/FST/results")
source("01_FST/FST_plot_crab_wave.R")
```

Plots are saved as 300 dpi PNG files (A5 landscape) in the same directory.

---

## Output files

| File | Description |
|------|-------------|
| `*.windowed.weir.fst` | VCFtools windowed F_ST output |
| `fst_weighted_plot_crab_all.png` | Genome-wide weighted F_ST, Crab |
| `fst_weighted_plot_crab_all_smooth.png` | + GAM smooth |
| `fst_weighted_plot_crab_CHR3_smooth.png` | LG3 zoom, Crab |
| `fst_mean_plot_crab_all.png` | Genome-wide mean F_ST, Crab |
| `fst_weighted_plot_wave_all.png` | Genome-wide weighted F_ST, Wave |
| `fst_weighted_plot_wave_all_smooth.png` | + GAM smooth |
| `fst_weighted_plot_wave_CHR9_smooth.png` | LG9 zoom, Wave |
| `fst_weighted_plot_wave_CHR8_smooth.png` | LG8 zoom, Wave |

---

## Citation

If you use this code, please cite:

> Ghane A, *et al.* (*in prep*). Sex determination in *Littorina saxatilis*: evidence
> from Koster populations.

And the underlying workflow:

> Grayson P, *et al.* (2022). SexFindR: a computational workflow to identify and
> characterise sex-linked regions. *Molecular Ecology Resources.*

---

## Contact

**Amin Ghane** — Institute of Science and Technology Austria (ISTA)
aminghane22@gmail.com