#!/bin/bash

#SBATCH --job-name=SeCNV_OV004-T_R2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=9-00:00:00
#SBATCH --output=/home/yangtianxia/my_code/secnv.OV004-T_R2.20221027.log
#SBATCH --error=/home/yangtianxia/my_code/secnv.OV004-T_R2.20221027.log

SAMPLE_ID="OV004-T_R2"

bash ~/my_code/se_cnv.sh --sample-id $SAMPLE_ID --input-file /data6/yangtianxia/neoantigen/DATA/SC_MATRIX/$SAMPLE_ID/$SAMPLE_ID.average.tsv --out-dir /data6/yangtianxia/neoantigen/DATA/SC_MATRIX/$SAMPLE_ID/SeCNV 


