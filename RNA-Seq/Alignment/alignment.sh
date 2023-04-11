#!/bin/sh
#SBATCH --job-name=RNA-SSNV-Prepare-Bam
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/yangtianxia/my_code/RNA-SSNV/rna_preprocess.20230328.log
#SBATCH --error=/home/yangtianxia/my_code/RNA-SSNV/rna_preprocess.20230328.log

echo `hostname` 
echo ======="Start at $(date)"========= && \
cd /data3/yangtianxia/SOFTWARE/RNA-SSNV && \
. /home/yangtianxia/my_code/RNA-SSNV/lmod_init_bash.sh && \
. env/bin/activate && \
snakemake --cores 6 \
          -s rules/RNA-Preparation-tsv-Snakefile.smk \
	  --configfile configs/rna_preprocess.yaml --rerun-incomplete && \
echo ======="Finish at $(date)"=========
