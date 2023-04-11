#!/bin/sh
#SBATCH --job-name=RNA-SSNV-Mutect
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/yangtianxia/my_code/RNA-SSNV/run_Mutect.slurm.log
#SBATCH --error=/home/yangtianxia/my_code/RNA-SSNV/run_Mutect.slurm.log

echo `hostname` 
echo ======="Start at $(date)"========= && \
cd /data3/yangtianxia/SOFTWARE/RNA-SSNV && \
. /home/yangtianxia/my_code/RNA-SSNV/lmod_init_bash.sh && \
. env/bin/activate && \
snakemake --cores 4 \
  -s rules/RNA-Somatic-tsv-Snakefile.smk \
  --configfile configs/project_config.yaml \
  --rerun-incomplete && \
echo ======="Finish at $(date)"=========
