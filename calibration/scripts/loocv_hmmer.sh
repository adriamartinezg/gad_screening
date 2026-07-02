CONDA_ENV="base"
REPS_DIR="clusters" 
GOLD_POS="gold/positives.faa"
GOLD_NEG="gold/negatives.faa"
OUTPUT_LOOCV="loocv"
PARSER="scripts/parser.sh"

mkdir -p "$OUTPUT_LOOCV"

conda activate $CONDA_ENV

for gene in "$REPS_DIR"/*_rep_seq.fasta; do
    gen_name=$(basename "$gene" _clustered_rep_seq.fasta)
    echo ">>> Iniciando LOOCV para el gen: $gen_name"
    
    gen_out="$OUTPUT_LOOCV/$gen_name"
    mkdir -p "$gen_out"
