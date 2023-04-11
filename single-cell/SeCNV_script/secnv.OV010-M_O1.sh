#!/bin/bash

#SBATCH --job-name=SeCNV_OV010-M_O1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=9-00:00:00
#SBATCH --output=/home/yangtianxia/my_code/secnv.OV010-M_O1.20221027.log
#SBATCH --error=/home/yangtianxia/my_code/secnv.OV010-M_O1.20221027.log

SAMPLE_ID="OV010-M_O1"

bash ~/my_code/se_cnv.sh --sample-id $SAMPLE_ID --input-file /data6/yangtianxia/neoantigen/DATA/SC_MATRIX/$SAMPLE_ID/$SAMPLE_ID.average.tsv --out-dir /data6/yangtianxia/neoantigen/DATA/SC_MATRIX/$SAMPLE_ID/SeCNV 


