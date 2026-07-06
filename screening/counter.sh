#!/bin/bash

OUTPUT="RESUMEN_PHYLUM.tsv"

echo -e "Phylum\tTotal_Genomes\tCOMPLETE\tPARTIAL_ENZYME\tPARTIAL_TRANSPORTER\tABSENT\t%_COMPLETE\t%_PARTIAL\t%_ABSENT" > "$OUTPUT"

echo "Procesando archivos..."

for file in GAD_system_*.tsv; do
    phylum=$(echo "$file" | sed 's/GAD_system_//;s/\.tsv//')
    
    counts=$(tail -n +2 "$file" | awk '
        {
            total++
            if ($NF == "COMPLETE") comp++
            else if ($NF == "PARTIAL_ENZYME_ONLY") p_enz++
            else if ($NF == "PARTIAL_TRANSPORTER_ONLY") p_tra++
            else if ($NF == "ABSENT") abs++
        }
        END {
            # Cálculo porcentajes
            if (total > 0) {
                p_comp = (comp / total) * 100
                p_part = ((p_enz + p_tra) / total) * 100
                p_abs = (abs / total) * 100
            } else {
                p_comp = p_part = p_abs = 0
            }
            
            # Resultados separados por tabulador
            printf "%d\t%d\t%d\t%d\t%d\t%.2f%%\t%.2f%%\t%.2f%%", 
            total, comp, p_enz, p_tra, abs, p_comp, p_part, p_abs
        }')

    echo -e "$phylum\t$counts" >> "$OUTPUT"
done

echo "--------------------------------------------------"
column -t -s $'\t' "$OUTPUT"
echo "--------------------------------------------------"
echo "Resumen guardado en: $OUTPUT"