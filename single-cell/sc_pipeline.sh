#!/bin/sh

MY_DIR=$( dirname -- "${BASH_SOURCE[0]}" )
sample_id=
WORK_DIR=

while [[ $# -gt 0 ]]; do
  case $1 in
    --id)
      sample_id="$2"
      shift 
      shift 
      ;;
    --work-dir)
      WORK_DIR="$2"
      shift
      shift
      ;;
    -*|--*)
      echo "Unknown option $1"
      echo "Usage $0 --id <SAMPLE_ID> --work-dir <WORK_DIR>"
      exit 1
      ;;
  esac
done

if [ -z $WORK_DIR -o -z $sample_id ]; then
    printf "Required parameters not all set!"
    exit 1
fi

. /home/yangtianxia/my_code/RNA-SSNV/lmod_init_bash.sh && \
module purge && module reset && \
module load 'GCC/11.2.0' 'OpenMPI/4.1.1' 'Python/3.9.6' 'SciPy-bundle/2021.10' 'matplotlib' 'scikit-learn/1.0.1' 'scVelo/0.2.4' 'R/4.2.0' && \
export PYTHONPATH="$PYTHONPATH:/home/yangtianxia/SOFTWARE/pyensembl"; 
export PYTHONPATH="$PYTHONPATH:/home/yangtianxia/SOFTWARE/rpy2";
python $MY_DIR/mtx_restore.py "$sample_id" $WORK_DIR && \
python $MY_DIR/mtx_average.py \
        --id "$sample_id" \
        --in $WORK_DIR/$sample_id.tsv \
        --out $WORK_DIR/$sample_id.average.tsv \
        --tmp /data4/yangtianxia/TMP && \
python $MY_DIR/SeCNV_script/SeCNV_RNA.py \
              --input_file $WORK_DIR/$sample_id.average.tsv \
              --output_dir $WORK_DIR/SeCNV \
              --KR_norm --save_figures --ensembl 78 && \
python $MY_DIR/rm_column.py --list $WORK_DIR/SeCNV/normal_list.csv \
            --in $WORK_DIR/$sample_id.average.tsv \
            --out $WORK_DIR/$sample_id.average.filtered.csv;
echo "Finished with $sample_id. See output in $WORK_DIR"