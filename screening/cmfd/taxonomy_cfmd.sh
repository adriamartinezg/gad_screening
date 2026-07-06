#!/bin/bash

mkdir -p results_cfmd/taxonomy

CFMD_LIST="results_cfmd/taxonomy/cFMD_mags_list.tsv"
if [ ! -f "$CFMD_LIST" ]; then
    echo "Descargando lista de taxonomía de cFMD..."
    wget "https://raw.githubusercontent.com/SegataLab/cFMD/refs/heads/main/cFMD_mags_list.tsv" \
     -O "$CFMD_LIST"
fi

awk -F'\t' '$9 == "COMPLETE" || $9 == "ABSENT" {print $3 "\t" $5 "\t" $9}' GAD_SYSTEM_MAG_LEVEL_RESULTS.tsv > mag_ids_all.tmp

OUTPUT="results_cfmd/taxonomy/final_taxonomy_results_all.tsv"

echo -e "MAG_id\tCategory\tSystem_Status\tSpecies\tKuperkingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tCompleteness\tContamination" > "$OUTPUT"

echo "Asignando taxonomía a los MAGs..."
while read -r mag_id category status; do
    line=$(grep -w "$mag_id" "$CFMD_LIST")
    if [ -n "$line" ]; then
        specie=$(echo "$line" | cut -f14)
        skingdom=$(echo "$line" | cut -f8)
        phylum=$(echo "$line" | cut -f9)
        class=$(echo "$line" | cut -f10)
        order=$(echo "$line" | cut -f11)
        family=$(echo "$line" | cut -f12)
        genus=$(echo "$line" | cut -f13)
        completeness=$(echo "$line" | cut -f18)
        contamination=$( echo "$line" | cut -f19)

        echo -e "$mag_id\t$category\t$status\t$specie\t$skingdom\t$phylum\t$class\t$order\t$family\t$genus\t$completeness\t$contamination" >> "$OUTPUT"
    fi
done < mag_ids_all.tmp

[ -f mag_ids_all.tmp ] && rm mag_ids_all.tmp

echo "¡Proceso finalizado con éxito!"
echo "Archivos generados en results_cfmd/taxonomy/"

