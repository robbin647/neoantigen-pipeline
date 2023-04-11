#!/bin/sh
#SBATCH --job-name=RNA-SSNV-OV001-T_L_OV001-M_P1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/data3/yangtianxia/SOFTWARE/vep/auto_scripts/snakemake.20230325.log
#SBATCH --error=/data3/yangtianxia/SOFTWARE/vep/auto_scripts/snakemake.20230325.log

echo `hostname`;
echo ======="Start at $(date)"========= && \
cd /data3/yangtianxia/SOFTWARE/vep/auto_scripts && \
. /data1/apps/software/Anaconda3/2022.05/etc/profile.d/conda.sh && \
conda activate "/data3/yangtianxia/SOFTWARE/vep" && \

snakemake -s workflow/Snakefile --configfile configs/config.yaml --cores 2 --rerun-incomplete && \
echo ======="Finish at $(date)"=========

