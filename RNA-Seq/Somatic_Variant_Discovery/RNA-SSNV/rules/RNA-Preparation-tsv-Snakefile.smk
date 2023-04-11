# RNA相关
# 从原始RNA fastq.gz文件起始，到可用以突变检测的RNA bam文件为终止
# 应用STAR完成alignment和mapping操作，应用GATK相关工具完成对RNA reads的前处理事宜

"""
MarkDuplicate过程中出现问题"Mapped mate should have mate reference name"
https://github.com/cbrueffer/tophat-recondition/issues/1
"""
import pandas as pd
import gzip
import os

from snakemake.utils import min_version

min_version("3.2")

workdir: config["workDir"]

# Prepare corresponding sample info dataframe
samples = pd.read_table(config["samples"]).set_index("file_name", drop=False)
# remove duplicate case fastq files
# samples.drop_duplicates(subset=["aliquots_id"], keep='first', inplace=True)

# Prepare lists to iterate over
aliquots_ids = set(samples["aliquots_id"])

##### 2023-03-10 Find PU information from source FASTQ files #####
# Platform Unit (PU) information for each sample
platform_units = {}
def extract_flowcell_lanes(fastq_file: str)-> tuple:
    """Extract flowcell id and lane number from the first line of a FASTQ/FASTQ.GZ file"""
    # add code to decompress fastq_file as a gzip file 
    # if it ends with .gz
    if fastq_file.endswith(".gz"):
        with gzip.open(fastq_file, 'rt') as f:
            first_line = f.readline().strip()
            flowcell_id = first_line.split(":")[2]
            lane_number = first_line.split(":")[3]
            return (flowcell_id, lane_number)
    elif fastq_file.endswith(".fastq"):
        with open(fastq_file, 'rt') as f:
            first_line = f.readline().strip()
            flowcell_id = first_line.split(":")[2]
            lane_number = first_line.split(":")[3]
            return (flowcell_id, lane_number)
    else:
        raise ValueError("File type not supported. File must be .fastq or .fastq.gz")
# iterate through each aliquots_id
for aliquot in aliquots_ids:
    fastq_file_basename = samples.loc[samples['aliquots_id'] == aliquot, ].head(1).index[0]
    fastq_file_path = os.path.join(config["sampleDir"], fastq_file_basename)
    flowcell_id, lane_number = extract_flowcell_lanes(fastq_file_path)
    platform_units[aliquot] = "{}.{}".format(flowcell_id, lane_number)
##################################################

rule all:
    input:
        expand("STAR-2-Pass/{aliquots_id}/{aliquots_id}.bam", aliquots_id = aliquots_ids),
        expand("STAR-2-Pass/{aliquots_id}/{aliquots_id}.bam.bai", aliquots_id = aliquots_ids)
        # expand("STAR-2-Pass/{aliquots_id}/{aliquots_id}-PU-info.txt", aliquots_id = aliquots_ids)

##### Test Find PU information from source FASTQ files ######
#rule test_pu:
#    params:
#        RG_PU = lambda wildcards: platform_units[wildcards.aliquots_id]
#    output:
#        result = "STAR-2-Pass/{aliquots_id}/{aliquots_id}-PU-info.txt"
#    shell:
#        """
#        if [ ! -f {output.result} ]; then
#            touch {output.result} 
#        fi
#        printf $"{wildcards.aliquots_id}\t{params.RG_PU}" > {output.result}
#        """

# specific for STAR-2.6.0c
rule STAR_2_pass:
    input:
        fastqs = lambda wildcards : [os.path.join(config["sampleDir"], samples["file_name"][file_name])
                                        for file_name in samples.loc[samples["aliquots_id"]==wildcards.aliquots_id, ].index]
        # fastq1 = lambda wildcards: os.path.join(config["sampleDir"], wildcards.aliquots_id+"m", wildcards.aliquots_id+"_1.fq.gz"),
        # fastq2 = lambda wildcards: os.path.join(config["sampleDir"], wildcards.aliquots_id+"m", wildcards.aliquots_id+"_2.fq.gz")
    params:
        STAR_genome_dir = config["ref"]["genome_dir"],
        RG_ID = "{aliquots_id}",
        RG_SM = "{aliquots_id}",
        RG_PL = "ILLUMINA",
        RG_LB = "LIB1",
        ##### Find PU information from source FASTQ files #####
        RG_PU = lambda wildcards: platform_units[wildcards.aliquots_id],
        #########################################
        output_prefix = "STAR-2-Pass/{aliquots_id}/{aliquots_id}-"
    output:
        bam = "STAR-2-Pass/{aliquots_id}/{aliquots_id}.bam"
    threads: 8
    shell:
        """
        STAR \
        --readFilesIn {input.fastqs} \
        --outSAMattrRGline ID:{params.RG_ID} SM:{params.RG_SM} PL:{params.RG_PL} LB:{params.RG_LB} PU:{params.RG_PU} \
        --alignIntronMax 1000000 \
        --alignIntronMin 20 \
        --alignMatesGapMax 1000000 \
        --alignSJDBoverhangMin 1 \
        --alignSJoverhangMin 8 \
        --alignSoftClipAtReferenceEnds Yes \
        --chimJunctionOverhangMin 15 \
        --chimMainSegmentMultNmax 1 \
        --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
        --chimSegmentMin 15 \
        --genomeDir {params.STAR_genome_dir} \
        --genomeLoad NoSharedMemory \
        --limitSjdbInsertNsj 1200000 \
        --outFileNamePrefix {params.output_prefix} \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS nM NM ch \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --quantMode TranscriptomeSAM GeneCounts \
        --readFilesCommand zcat \
        --runThreadN {threads} \
        --twopassMode Basic && mv {params.output_prefix}Aligned.sortedByCoord.out.bam {output.bam}
        """

rule bam_index:
    input:
        bam = "STAR-2-Pass/{aliquots_id}/{aliquots_id}.bam"
    output:
        index = "STAR-2-Pass/{aliquots_id}/{aliquots_id}.bam.bai"
    threads: 2
    shell:
        """
        samtools index -@ {threads} {input.bam} {output.index} && \
        printf "\n[$(date)]\n=====Index for BAM {input.bam} was created!====\n"
        """
