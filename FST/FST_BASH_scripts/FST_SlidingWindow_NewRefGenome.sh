#!/bin/bash
# =============================================================================
# Sex-biased FST Analysis in Littorina saxatilis (Koster populations)
# Sliding-window FST between males and females, calculated per ecotype
# using the Weir & Cockerham estimator (VCFtools).
#
# Written by Amin Ghane
# Last Update in 2026-05
# =============================================================================

# --- SLURM job configuration -------------------------------------------------
#SBATCH -A naiss2023-5-264
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -J FST_SlidingWindow_NewRefGenome
#SBATCH --output=/proj/snic2022-6-238/Amin/FST_newRef/%x_%A.log
#SBATCH --error=/proj/snic2022-6-238/Amin/FST_newRef/%x_%A.err
#SBATCH --mail-type=FAIL,TIME_LIMIT_80
#SBATCH --mail-user=aminghane22@gmail.com

# --- Parameters --------------------------------------------------------------
VCF=/proj/snic2022-6-238/Amin/FST_newRef/filtered_with_bcftools_7_3.vcf.gz
WINDOW_SIZE=10000
WINDOW_STEP=100
OUTDIR=/proj/snic2022-6-238/Amin/FST_newRef

# --- Environment setup -------------------------------------------------------
echo "Job started  : $(date)"
echo "Running on   : $(hostname)"
echo "Working dir  : ${OUTDIR}"

module load bioinfo-tools
module load vcftools
module list

# --- Analysis ----------------------------------------------------------------
cd "${OUTDIR}"

# Crab ecotype: male vs. female FST in sliding windows
echo "[1/2] Computing FST — Crab ecotype..."
vcftools --gzvcf "${VCF}" \
    --weir-fst-pop snpDen_males_crab.txt \
    --weir-fst-pop snpDen_females_crab.txt \
    --fst-window-size "${WINDOW_SIZE}" \
    --fst-window-step "${WINDOW_STEP}" \
    --out sex_saxatilis_crab_fst_vcftools_window

# Wave ecotype: male vs. female FST in sliding windows
echo "[2/2] Computing FST — Wave ecotype..."
vcftools --gzvcf "${VCF}" \
    --weir-fst-pop snpDen_males_wave.txt \
    --weir-fst-pop snpDen_females_wave.txt \
    --fst-window-size "${WINDOW_SIZE}" \
    --fst-window-step "${WINDOW_STEP}" \
    --out sex_saxatilis_wave_fst_vcftools_window

echo "Job finished : $(date)"