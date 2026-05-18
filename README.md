# Sex-associated genomic variation in *Littorina saxatilis* — Koster analyses

Bioinformatics code accompanying the manuscript on sex determination in *Littorina saxatilis*
(Ghane et al., *in prep*). Analyses were performed on whole-genome sequencing data from
Crab and Wave ecotype individuals collected from the Koster archipelago, southwest Sweden.

---

## Background

*Littorina saxatilis* has a female heterogametic (ZW/ZZ) sex determination system in the
Crab ecotype, associated with chromosomal inversions on linkage group 12 (Hearn et al., 2022).
Whether a comparable system exists in the Wave ecotype remains unclear. To identify
sex-linked genomic regions in both ecotypes, we applied the SexFindR workflow
(Grayson et al., 2022), which combines four complementary approaches: genomic coverage,
SNP density, F_ST, and GWAS. Analyses were performed on 98 individuals (Crab: 15M / 33F;
Wave: 15M / 19F) mapped to the *L. saxatilis* reference genome v2.0 (De Jode, unpublished).

---

## Analyses

### 1. F_ST

Male–female F_ST was calculated per ecotype using the Weir & Cockerham estimator in
VCFtools, using a sliding-window approach (50 kb windows, 10 kb steps). Windows below
the genome-wide median F_ST were excluded, and results are visualised as Manhattan plots
across the 17 linkage groups.

- `FST/FST_BASH_scripts/FST_SlidingWindow_NewRefGenome.sh` — runs VCFtools on the
  cluster (SLURM)
- `FST/FST_R_scripts/FST_SlidingWindow_NewRefGenome.R` — produces Manhattan plots
  for weighted and mean F_ST (Crab and Wave)

### 2. SNP density

Per-sample SNP density was calculated in non-overlapping 10 kb windows. Mean densities
for males and females were compared, and a 1,000-permutation test (sex labels shuffled)
was used to assign window-level p-values. Windows with p ≤ 0.01 were retained as
candidates.

- `SNP_density/` *(scripts to be added)*

### 3. GWAS

Sex was treated as a binary phenotype and a genome-wide association study was performed
using GEMMA. Only sites with a likelihood ratio test p-value < 0.05 were used for
Manhattan plotting.

- `GWAS/` *(scripts to be added)*

### 4. Genomic coverage

Male–female read coverage was compared across the genome in 20 male–female pairs using
DifCover. The absence of consistent coverage differences confirmed that large heteromorphic
sex chromosome regions are not present in *L. saxatilis*, consistent with its karyotype.

- `coverage/` *(scripts to be added)*

---

## Data

The filtered VCF and reference genome are not hosted here due to file size. Please get
in touch if you need access to the raw data.

---

## Related publications

Hearn, K. E., Koch, E. L., Stankowski, S., Butlin, R. K., Faria, R., Johannesson, K., &
Westram, A. M. (2022). Differing associations between sex determination and sex-linked
inversions in two ecotypes of *Littorina saxatilis*. *Evolution Letters*, 6(5), 358–374.
https://doi.org/10.1002/evl3.295

Reeve, J., Butlin, R. K., Koch, E. L., Stankowski, S., & Faria, R. (2023). Chromosomal
inversion polymorphisms are widespread across the species ranges of rough periwinkles
(*Littorina saxatilis* and *L. arcana*). *Molecular Ecology*, 33(24), e17160.
https://doi.org/10.1111/mec.17160

Grayson, P., et al. (2022). SexFindR: a computational workflow to identify and characterise
sex-linked regions. *Molecular Ecology Resources*.

---

## Contact

**Amin Ghane** — University of Vienna
amin.ghane@univie.ac.at