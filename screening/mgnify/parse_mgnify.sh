#!/bin/bash

# Nombre del archivo final (se guardará en la raíz tfm_app)
OUTPUT="MGNIFY_FINAL_GAD.tsv"

# 1. Escribir la cabecera
echo -e "Genome_ID\tBiome\tStatus\tSpecies\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tCompleteness\tContamination" > "$OUTPUT"

# 2. Definir la raíz de los datos
DATA_ROOT="results_mgnify_meta"

# 3. Iterar por cada carpeta de bioma dentro de results_mgnify_meta
for BIOME_PATH in "$DATA_ROOT"/*/; do
    # Extraer solo el nombre de la carpeta (ej: barley-rhizosphere)
    BIOME_NAME=$(basename "$BIOME_PATH")
    
    # Definir rutas exactas según tu estructura
    META_FILE="${BIOME_PATH}genomes-all_metadata.tsv"
    HITS_FILE="${BIOME_PATH}hits/summary/GAD_dataset_classification.tsv"

    echo "Procesando bioma: $BIOME_NAME..."

    # Verificar que existen ambos archivos
    if [[ -f "$META_FILE" && -f "$HITS_FILE" ]]; then
        
        # 4. Unir con AWK
        awk -F'\t' -v biome="$BIOME_NAME" '
        # Leer Metadata primero
        NR==FNR {
            if (FNR == 1) {
                for(i=1;i<=NF;i++){
                    if($i=="Genome") c_id=i;
                    if($i=="Completeness") c_comp=i;
                    if($i=="Contamination") c_cont=i;
                    if($i=="Lineage") c_lin=i;
                }
                next
            }
            # Guardar: Comp | Cont | Lin
            meta[$c_id] = $c_comp "\t" $c_cont "\t" $c_lin
            next
        }
        # Leer Clasificación
        {
            if (FNR == 1) next;
            id = $1; status = $6
            if (id in meta) {
                split(meta[id], a, "\t")
                comp = a[1]; cont = a[2]; lin = a[3]
                
                # Parsear Linaje GTDB
                split(lin, t, ";")
                d="Unassigned"; p="Unassigned"; c="Unassigned"; o="Unassigned"; f="Unassigned"; g="Unassigned"; s="Unassigned"
                for (i in t) {
                    if (t[i] ~ /^d__/) d = substr(t[i], 4)
                    if (t[i] ~ /^p__/) p = substr(t[i], 4)
                    if (t[i] ~ /^c__/) c = substr(t[i], 4)
                    if (t[i] ~ /^o__/) o = substr(t[i], 4)
                    if (t[i] ~ /^f__/) f = substr(t[i], 4)
                    if (t[i] ~ /^g__/) g = substr(t[i], 4)
                    if (t[i] ~ /^s__/) s = substr(t[i], 4)
                }
                # Rellenar vacíos
                if(d=="")d="Unassigned"; if(p=="")p="Unassigned"; if(c=="")c="Unassigned"
                if(o=="")o="Unassigned"; if(f=="")f="Unassigned"; if(g=="")g="Unassigned"; if(s=="")s="Unassigned"
                
                print id "\t" biome "\t" status "\t" s "\t" d "\t" p "\t" c "\t" o "\t" f "\t" g "\t" comp "\t" cont
            }
        }
        ' "$META_FILE" "$HITS_FILE" >> "$OUTPUT"
    else
        echo "   [!] Faltan archivos en $BIOME_NAME."
        [[ ! -f "$META_FILE" ]] && echo "       Falta: $META_FILE"
        [[ ! -f "$HITS_FILE" ]] && echo "       Falta: $HITS_FILE"
    fi
done

echo "----------------------------------------------------------"
echo "PROCESO FINALIZADO. Archivo generado: $OUTPUT"