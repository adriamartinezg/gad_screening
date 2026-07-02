#!/bin/bash
set -e

FASTA_FILE=""
OUTPUT_TSV=""

while getopts "f:o:h" opt; do
    case $opt in
        f) FASTA_FILE="$OPTARG" ;;
        o) OUTPUT_TSV="$OPTARG" ;;
        h) echo "Uso: $0 -f secuencias.fasta -o database.tsv"; exit 0 ;;
    esac
done

if [ ! -f "$FASTA_FILE" ]; then
    echo "Error: No se encuentra el archivo FASTA: $FASTA_FILE"
    exit 1
fi

# 1. Limpieza de temporales
rm -rf biolib_results

# 2. Ejecución única en BioLib
echo "Enviando $FASTA_FILE a BioLib (DeepTMHMM)..."
biolib run DTU/DeepTMHMM --fasta "$FASTA_FILE"

# 3. Verificación de resultados
RESULT_FILE="biolib_results/TMRs.gff3"
if [ ! -f "$RESULT_FILE" ]; then
    echo "Error: BioLib no generó resultados."
    exit 1
fi

# 4. Crear la Base de Datos (TSV con ID, Longitud y Número de TMs)
# Parseamos el GFF3 para extraer: ID, Length y el conteo de TMRs
echo "Generando base de datos en $OUTPUT_TSV..."
awk '/^#/ {
    seqid = $2; 
    for (i=3; i<=NF; i++) {
        if ($i == "Length:") l[seqid] = $(i+1); 
        if ($i == "TMRs:") t[seqid] = $(i+1)
    }
} 
END {
    # Imprimimos cabecera
    print "ID\tseq_length\tn_tm";
    for (s in l) {
        print s "\t" l[s] "\t" (t[s]?t[s]:0)
    }
}' "$RESULT_FILE" > "$OUTPUT_TSV"

echo "¡Listo! Base de datos creada con $(wc -l < "$OUTPUT_TSV") entradas."
