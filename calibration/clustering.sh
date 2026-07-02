#!/bin/bash

CONDA_ENV="base"
INPUT_DIR="."
OUTPUT_DIR="./results"

while getopts "e:i:o:h" opt; do
    case ${opt} in
        e ) CONDA_ENV="$OPTARG" ;;
        i ) INPUT_DIR="$OPTARG" ;;
        o ) OUTPUT_DIR="$OPTARG" ;;
        h )
            echo "Usage: $0 [-e conda environment] [-i input directory] [-o output directory]"
            exit 0
            ;;
        \? )
            echo "Invalid option. Use -h for help."
            exit 1
            ;;
    esac
done

mkdir -p "$OUTPUT_DIR"
TMP_DIR="tmp_mmseqs_global"
mkdir -p "$TMP_DIR"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"

shopt -s nullglob
faa_files=("$INPUT_DIR"/*.faa)

if [ ${#faa_files[@]} -eq 0 ]; then
    echo "Error: No .faa files found in $INPUT_DIR"
    exit 1
fi

for file in "${faa_files[@]}"; do
    # 1. Obtener solo el nombre del archivo
    base_name=$(basename "$file" .faa)
    
    prefix_out="$OUTPUT_DIR/${base_name}_clustered"

    echo ">>> Procesando: $base_name ..."

    # Ejecutar MMseqs2
    mmseqs easy-cluster "$file" "$prefix_out" "$TMP_DIR" \
        --min-seq-id 0.7 \
        -c 0.8 \
        --cov-mode 0 \
        -v 3

    echo "Finalizado: $base_name"
    echo "Archivos generados en $OUTPUT_DIR con prefijo ${base_name}_clustered"
    echo "------------------------------------------------"
done

# Limpieza final
rm -rf "$TMP_DIR"
echo "Proceso completo."