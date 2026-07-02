#!/bin/bash

while getopts "i:o:" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
  esac
done

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Uso: $0 -i input.domtbl -o output.tsv"
    exit 1
fi

# Escribir cabecera
echo -e "SeqID\tQuery\tE-value_full_sequence\tScore_full_sequence\tcE-value_best_domain\tiE-value_best_domain\tScore_best_domain\tCoverage\tDescription" > "$OUTPUT"

awk '!/^#/ {
    seqid=$1; tlen=$3; query=$4; qlen=$6;
    seqEvalue=$7; seqscore=$8;
    dom_cEvalue=$12; dom_iEvalue=$13; domscore=$14;
    hmmfrom=$16; hmmto=$17;

    # Cálculo de coverage
    cov = (((hmmto - hmmfrom) + 1) / qlen) * 100;

    desc = ""; for (i=23; i<=NF; i++) desc = desc $i " ";

    # Solo guardar si es el mejor score (columna 14) para este ID
    if (!(seqid in max_score) || domscore > max_score[seqid]) {
        max_score[seqid] = domscore;
        line[seqid] = seqid "\t" query "\t" seqEvalue "\t" seqscore "\t" dom_cEvalue "\t" dom_iEvalue "\t" domscore "\t" cov "\t" desc;
    }
}
END {
    for (id in line) print line[id];
}' "$INPUT" >> "$OUTPUT"
