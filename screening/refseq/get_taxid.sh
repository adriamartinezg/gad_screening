#!/bin/bash

PATTERN="results_bacteria_classes/summary/GAD_*_classification.tsv"
FILES=$(ls $PATTERN 2>/dev/null)

if [ -z "$FILES" ]; then
    echo "Error: No se encontraron archivos con el patrón $PATTERN"
    exit 1
fi

echo "Iniciando procesamiento de múltiples clases..."

Rscript genomes_status.R $FILES

echo "Proceso terminado."