#!/bin/bash

# Written by Amin Ghane
# Last Update in 2026-05
### SNP density calculation per individual, using the new reference genome

### ============================================================
### SLURM job parameters
### ============================================================
#SBATCH -A naiss2023-5-264
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -J SNP_density
#SBATCH --output=%x_%A.log
#SBATCH --error=%x_%A.err
#SBATCH --mail-type=FAIL,TIME_LIMIT_80
#SBATCH --mail-user=aminghane22@gmail.com

echo "Start time: $(date)"

### ============================================================
### Required software
### ============================================================
module load bioinfo-tools
module load bcftools   # tested with bcftools/1.17
module load vcftools   # tested with vcftools/0.1.16
module list

### ============================================================
### Arguments
### ============================================================
# $1  sample name (or comma-separated sample list) to subset
# $2  output VCF path (uncompressed)
# $3  input VCF/BCF path (can be bgzipped)

SAMPLE=$1
OUT_VCF=$2
IN_VCF=$3

### ============================================================
### Usage
### ============================================================
# Submit one job per individual from a sample list.
# Example — Koster Crab ecotype:
#
#   for file in $(cat snpDen_females_crab.txt); do
#       sbatch SNP_density.sh \
#           "$file" \
#           individual_SNP_density_crab/Female_crab_${file%%.*}.vcf \
#           filtered_with_bcftools_7_3.vcf.gz
#       sleep 0.1
#   done
#
#   for file in $(cat snpDen_males_crab.txt); do
#       sbatch SNP_density.sh \
#           "$file" \
#           individual_SNP_density_crab/Male_crab_${file%%.*}.vcf \
#           filtered_with_bcftools_7_3.vcf.gz
#       sleep 0.1
#   done
#
# Example — Koster Wave ecotype:
#
#   for file in $(cat snpDen_females_wave.txt); do
#       sbatch SNP_density.sh \
#           "$file" \
#           individual_SNP_density_wave/Female_wave_${file%%.*}.vcf \
#           filtered_with_bcftools_7_3.vcf.gz
#       sleep 0.1
#   done
#
#   for file in $(cat snpDen_males_wave.txt); do
#       sbatch SNP_density.sh \
#           "$file" \
#           individual_SNP_density_wave/Male_wave_${file%%.*}.vcf \
#           filtered_with_bcftools_7_3.vcf.gz
#       sleep 0.1
#   done

### ============================================================
### Step 1 — Subset and trim alleles for target sample(s)
### ============================================================
# bcftools view options used:
#   -a / --trim-alt-alleles  remove ALT alleles absent in the subset genotypes
#                            (ALT is set to "." if none remain; record is kept)
#   -s / --samples           comma-separated sample ID(s) to retain
#   -o / --output            output file

bcftools view -a -s "$SAMPLE" -o "$OUT_VCF" "$IN_VCF"

### ============================================================
### Step 2 — Compress and index the subset VCF
### ============================================================
bgzip -c "$OUT_VCF" > "${OUT_VCF}.gz"
bcftools index "${OUT_VCF}.gz"

### ============================================================
### Step 3 — Calculate SNP density in 10 kb windows
### ============================================================
vcftools --vcf "$OUT_VCF" --SNPdensity 10000 --out "$OUT_VCF"

echo "End time: $(date)"