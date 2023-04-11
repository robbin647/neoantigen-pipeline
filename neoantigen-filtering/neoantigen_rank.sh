#!/bin/bash

MATRIX_CORRELATE=/root/software/scikit-bio/matrix_correlate.py
SC_TREE=/root/software/scikit-bio/rank_sc_tree.py
IN_SC=
IN_PEP=
SAMPLE_ID=
WORK_DIR=

# get command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --sc)
      IN_SC="$2"
      shift
      shift
      ;;
    --pep)
      IN_PEP="$2"
      shift
      shift
      ;;
    --id)
      SAMPLE_ID="$2"
      shift
      shift
      ;;
    --out)
      WORK_DIR="$2"
      shift
      shift
      ;;
    -*|--*)
      echo "Usage: $0 --sc <SINGLE_CELL> --pep <PEPTIDE> --id <SAMPLE_ID> --out <WORK_DIR>"
      exit 1
      ;;
  esac
done

echo "IN_SC=$IN_SC"
echo "IN_PEP=$IN_PEP"
echo "SAMPLE_ID=$SAMPLE_ID"
echo "WORK_DIR=$WORK_DIR"
# create working directory if not exist
if [ ! -d $WORK_DIR ]; then
  mkdir -p $WORK_DIR
fi;

source /root/software/scikit-bio/env/bin/activate && \
python $MATRIX_CORRELATE --sc $IN_SC --peptide $IN_PEP \
                         --out "$WORK_DIR/$SAMPLE_ID.best_pep_cov_mtx.tsv" \
                         --filtered-aggr "$WORK_DIR/$SAMPLE_ID.neoantigen.filtered.tsv" && \
python $SC_TREE --in "$WORK_DIR/$SAMPLE_ID.best_pep_cov_mtx.tsv" \
                --sample-id $SAMPLE_ID \
                --top-n-out "$WORK_DIR/$SAMPLE_ID.top_n_neoantigens.txt" \
                --newick-out "$WORK_DIR/$SAMPLE_ID.newick" && \
printf "Done. Newick tree and filtered neoantigen list were written to %s\n" $WORK_DIR
