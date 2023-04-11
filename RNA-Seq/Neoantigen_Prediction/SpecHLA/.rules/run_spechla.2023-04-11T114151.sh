
echo `hostname`; \ 
echo ========Start at $(date)=========== && \
cd /data3/yangtianxia/SOFTWARE/spechla && \
. /data1/apps/software/Anaconda3/2022\.05/etc/profile.d/conda.sh && \
conda activate "/data3/yangtianxia/SOFTWARE/spechla/spechla_env" && \
conda env list && \
    
echo [`date`] Working on OV005-T_R2 && \
mkdir -p /data5/yangtianxia/spechla//OV005-T_R2 && \
mkdir -p /data5/yangtianxia/spechla//OV005-T_R2/extract_hla && \
bash /data3/yangtianxia/SOFTWARE/spechla/script/ExtractHLAread.sh -s OV005-T_R2 \
     -b /data4/yangtianxia/RNA-SSNV/bam_20230228/OV005-T_R2/OV005-T_R2.bam \
     -r hg38 \
     -o /data5/yangtianxia/spechla//OV005-T_R2/extract_hla && \
mkdir -p /data5/yangtianxia/spechla//OV005-T_R2/spec_hla && \
bash /data3/yangtianxia/SOFTWARE/spechla/script/whole/SpecHLA.sh -n OV005-T_R2 \
     -1 /data5/yangtianxia/spechla//OV005-T_R2/extract_hla/OV005-T_R2_extract_1.fq.gz -2 /data5/yangtianxia/spechla//OV005-T_R2/extract_hla/OV005-T_R2_extract_2.fq.gz \
     -o /data5/yangtianxia/spechla//OV005-T_R2/spec_hla \
     -u 1 \
     -j 1 && \
echo [`date`] Done with OV005-T_R2 && \
        
echo ========Finish at $(date)=============