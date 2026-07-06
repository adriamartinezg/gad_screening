#!/bin/bash

usage() {
    echo "Usage: $0 [-m models] [-r results] [-t tc.csv] [-s file] [-x taxon] [-f]"
    exit 1
}

MODELS_DIR=""
RESULTS_DIR=""
TC_FILE=""
SINGLE_FILE=""
FINAL_ONLY=false
TAXON_NAME="dataset"

while getopts "m:r:t:s:x:fh" opt; do
    case $opt in
        m) MODELS_DIR=$OPTARG ;;
        r) RESULTS_DIR=$OPTARG ;;
        t) TC_FILE=$OPTARG ;;
        s) SINGLE_FILE=$OPTARG ;;
        x) TAXON_NAME=$OPTARG ;;
        f) FINAL_ONLY=true ;;
        *) usage ;;
    esac
done

RAW_DIR="$RESULTS_DIR/hmm_raw_$TAXON_NAME"
SUMMARY_DIR="$RESULTS_DIR/summary"
LOCAL_MODELS="$RESULTS_DIR/active_models"

mkdir -p "$RAW_DIR" "$SUMMARY_DIR" "$LOCAL_MODELS"

# ---------- cargar TC ----------
if [[ "$FINAL_ONLY" == false ]]; then
    sed 's/,/\t/g' "$TC_FILE" > summary.csv
    TC_GADA=$(awk -F'\t' '$1=="gadA"{print $2}' summary.csv)
    TC_GADB=$(awk -F'\t' '$1=="gadB"{print $2}' summary.csv)
    TC_GADC=$(awk -F'\t' '$1=="gadC"{print $2}' summary.csv)

    for gene in gadA gadB gadC; do
        cp "$MODELS_DIR/$gene.hmm" "$LOCAL_MODELS/"
    done
fi

# ---------- modo SINGLE ----------
if [[ -n "$SINGLE_FILE" ]]; then
    sample_id=$(basename "$SINGLE_FILE" | sed 's/\.[^.]*$//')

    for gene in gadA gadB gadC; do
        current_tc=$(eval echo \$TC_${gene^^})
        hmmsearch --noali \
          --domtblout "$RAW_DIR/${sample_id}_${gene}.domtblout" \
          -T "$current_tc" \
          "$LOCAL_MODELS/${gene}.hmm" "$SINGLE_FILE" >/dev/null
    done
    exit 0
fi

# ---------- reporte final ----------
FINAL_REPORT="$SUMMARY_DIR/GAD_${TAXON_NAME}_classification.tsv"
echo -e "Genome_ID\tgadA\tgadB\tgad_unassigned\tgadC\tSystem_Status" > "$FINAL_REPORT"

ls "$RAW_DIR"/*_gadA.domtblout | xargs -n1 basename | sed 's/_gadA.domtblout//' | sort -u | while read sample; do
    awk -v tc="$TC_GADA" '$1!~/^#/ && $14>=tc{print $1}' "$RAW_DIR/${sample}_gadA.domtblout" | sort -u > A
    awk -v tc="$TC_GADB" '$1!~/^#/ && $14>=tc{print $1}' "$RAW_DIR/${sample}_gadB.domtblout" | sort -u > B

    count_A=$(comm -23 A B | wc -l)
    count_B=$(comm -13 A B | wc -l)
    count_un=$(comm -12 A B | wc -l)

    count_C=$(awk -v tc="$TC_GADC" '$1!~/^#/ && $14>=tc{print $1}' "$RAW_DIR/${sample}_gadC.domtblout" | sort -u | wc -l)

    total=$((count_A+count_B+count_un))

    if [[ $total -gt 0 && $count_C -gt 0 ]]; then status="COMPLETE"
    elif [[ $total -gt 0 ]]; then status="PARTIAL_ENZYME_ONLY"
    elif [[ $count_C -gt 0 ]]; then status="PARTIAL_TRANSPORTER_ONLY"
    else status="ABSENT"; fi

    echo -e "$sample\t$count_A\t$count_B\t$count_un\t$count_C\t$status" >> "$FINAL_REPORT"
done

rm -f A B
