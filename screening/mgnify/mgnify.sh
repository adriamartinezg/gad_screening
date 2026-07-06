#!/bin/bash

# Configuración
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes"
BASE_DIR=$(pwd)
MODELS="$BASE_DIR/model"           # Ruta completa a carpeta de modelos
TC_FILE="$BASE_DIR/thresholds.csv" # Ruta completa al CSV
CLASSIFIER="$BASE_DIR/scripts/classifier.sh"
RESULTS_BASE="$BASE_DIR/results_mgnify_meta"
BIOMES=(
    "barley-rhizosphere" "chicken-gut" "cow-rumen" "honeybee-gut" 
    "human-gut" "human-oral" "human-skin" "human-vaginal" 
    "maize-rhizosphere" "marine-eukaryotes" "marine" "marine_sediment" 
    "mouse-gut" "non-model-fish-gut" "pig-gut" "sheep-rumen" 
    "soil" "tomato-rhizosphere" "zebrafish-fecal"
)

for BIOME in "${BIOMES[@]}"; do
    echo "----------------------------------------------------------"
    echo "PROCESSING BIOME: $BIOME"
    echo "----------------------------------------------------------"
    
    VERSION=$(curl -s "$BASE_URL/$BIOME/" | grep -oE 'v[0-9]+\.[0-9]+(\.[0-9]+)?' | sort -V | tail -n 1)
    
    if [ -z "$VERSION" ]; then
        echo "[!] Not any version was detected for biome $BIOME. Next..."
        continue
    fi
    echo "[+] Detected version: $VERSION"

    BIOME_DIR="$RESULTS_BASE/$BIOME"
    mkdir -p "$BIOME_DIR/faa_files"
    
    METADATA_URL="$BASE_URL/$BIOME/$VERSION/genomes-all_metadata.tsv"
    echo "[+] Downloading metadata..."
    wget -N -P "$BIOME_DIR" "$METADATA_URL"

    CATALOGUE_URL="$BASE_URL/$BIOME/$VERSION/species_catalogue/"
    BLOCKS=$(curl -s "$CATALOGUE_URL" | grep -oE 'MGYG[0-9]+/' | uniq)

    for BLOCK in $BLOCKS; do
            BLOCK_URL="${CATALOGUE_URL}${BLOCK}"
            SPECIES_IDS=$(curl -s "$BLOCK_URL" | grep -oE 'MGYG[0-9]+/' | uniq)
            
            for S_ID_RAW in $SPECIES_IDS; do
                S_ID=$(echo "$S_ID_RAW" | sed 's/\///')
                FAA_FILE="$BIOME_DIR/faa_files/${S_ID}.faa"

                if [ ! -f "$FAA_FILE" ]; then
                    echo "    [Download] $S_ID..."
                    echo "${BLOCK_URL}${S_ID_RAW}genome/${S_ID}.faa"
                    wget -O "$FAA_FILE" "${BLOCK_URL}${S_ID_RAW}genome/${S_ID}.faa"
                fi

            bash $CLASSIFIER -m "$MODELS" -r "$BIOME_DIR/hits" -t "$TC_FILE" -s "$FAA_FILE"
            rm "$FAA_FILE"
            done
        done
        echo "--- Finished bioma: $BIOME ---"
    done

echo "[+] Obtaining final report for $BIOME..."
    $CLASSIFIER -m "$MODELS" -r "$BIOME_DIR" -t "$TC_FILE" -x "$BIOME" -f
    
    echo "--> BIOME FINISHED. Report in: $BIOME_DIR/summary/"
done
