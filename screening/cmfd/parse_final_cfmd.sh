#!/bin/bash

# --- CONFIGURACIÓN ---
RESULTS_DIR="results_cfmd"
METADATA_FILE="$RESULTS_DIR/cFMD_metadata.tsv"
FINAL_TABLE="GAD_SYSTEM_MAG_LEVEL_RESULTS.tsv"
HITS_BASE="$RESULTS_DIR/hits"

echo -e "dataset_name\tsample_id\tMAG_id\tmacrocategory\tcategory\ttype\tgadEnzyme_count\tgadTransporter_count\tstatus" > "$FINAL_TABLE"

# Lista de datasets (puedes dejar tu array tal cual)
DATASETS=('AlvarezOrdonezA_xxxx' 'ArikanM_2020' 'BertuzziAS_2018' 'ChaconVargasK_2020' 'CM_INJERA' 'CrovadoreJ_2017' 'DeFilippisF_xxxx' 'DeRoosJ_2020' 'DuR_2020' 'DuruIC_2018' 'EinsonEJ_2018' 'EscobarZepedaA_2016' 'FerrocinoI_2018' 'HeilCS_2018' 'HellmannSL_2020' 'KastmanEK_2016' 'KawaiT_2012' 'KumarJ_2019' 'LandisEA_2021' 'LeechJ_2020' 'LeechJ_xxxx' 'LeonardSR_2016' 'LiZ_2018' 'LiZ_2019' 'LordanR_2019' 'MASTER_WP4_CSIC_1' 'MASTER_WP4_CSIC_2' 'MASTER_WP4_FFoQSI_1' 'MASTER_WP4_FFoQSI_2' 'MASTER_WP4_FFoQSI_3' 'MASTER_WP4_MATIS_1' 'MASTER_WP4_TEAGASC_1' 'MASTER_WP4_UNINA_1' 'MASTER_WP4_UNINA_2' 'McHughAJ_2020' 'MilaniC_2019' 'MortensenS_xxxx' 'PasolliE_2020' 'PatroJN_2016' 'PfeferT_xxxx' 'PorcellatoD_2016' 'PothakosV_2020' 'QuigleyL_2016' 'RippF_2014' 'SalvettiE_2016' 'SomervilleV_2019' 'SternesPR_2017' 'SulaimanJ_2014' 'VerceM_2019' 'WalshAM_2016' 'WalshAM_2017' 'WalshAM_2020' 'WalshL_xxxx' 'WolfeBE_2014' 'XieM_2019' 'YaoG_2017' 'YasirM_2020' 'YulandiA_2020' 'ZhaoCC_2020' 'ShangpliangH_2023_a' 'ShangpliangH_2023_b' 'KharnaiorP_2023' 'YouL_2022' 'LimaC_2020' 'DecadtH_2024' 'YasirM_2022' 'FontanaF_2023' 'FranciosaI_2021' 'MotaGutierrezJ_2021' 'SaakC_2023' 'YangC_2021' 'LopezSanchezR_2023' 'GonzalezOrozcoB_2023' 'SalgadoTS_2021' 'FalardeauJ_2023' 'SequinoG_2024_a' 'SequinoG_2024_b' 'MagliuloR_2024' 'YapM_2020' 'OlgaP_2019' 'QuijadaN_2022' 'YuY_2022' 'TomarS_2023' 'AlmeidaO_2020' 'CM_UNINA_FFOOD')

for DS in "${DATASETS[@]}"; do
    RAW_DIR="$HITS_BASE/$DS/hmm_raw_dataset"
    [ ! -d "$RAW_DIR" ] && continue

    # Generamos la lista de MAGs únicos detectados en la carpeta
    # sed elimina el sufijo para dejar solo el nombre del MAG
    MAG_LIST=$(ls "$RAW_DIR"/*.domtblout 2>/dev/null | sed -E 's/_(gadA|gadB|gadC)\.domtblout//' | sort -u)

    for MAG_ID in $MAG_LIST; do
        
        # 1. Obtenemos nombre limpio y SampleID
        CLEAN_MAG_ID=$(basename "$MAG_ID")
        SAMP_ID=$(echo "$CLEAN_MAG_ID" | sed -E 's/.*__(.*)__bin.*/\1/')
        
        # 2. Búsqueda exacta: Buscamos la fila donde Col1==DS y Col2==SAMP_ID
        # Usamos awk para filtrar por las dos columnas simultáneamente
        META_LINE=$(awk -F'\t' -v d="$DS" -v s="$SAMP_ID" '$1 == d && $2 == s {print; exit}' "$METADATA_FILE")

        # 3. Extraemos los campos (columna 3=macro, 4=category, 5=type)
        MACRO=$(echo "$META_LINE" | awk -F'\t' '{print $3}')
        CAT=$(echo "$META_LINE" | awk -F'\t' '{print $4}')
        TYPE=$(echo "$META_LINE" | awk -F'\t' '{print $5}')
        # Calcular conteos por MAG
        TOTAL_ENZ=0
        for suffix in "_gadA" "_gadB"; do
            f="$MAG_ID$suffix.domtblout"
            [ -f "$f" ] && COUNT=$(grep -vc "^#" "$f") && TOTAL_ENZ=$((TOTAL_ENZ + COUNT))
        done

        TOTAL_TRA=0
        f="${MAG_ID}_gadC.domtblout"
        [ -f "$f" ] && TOTAL_TRA=$(grep -vc "^#" "$f")

        # Definir status
        if [ $TOTAL_ENZ -gt 0 ] && [ $TOTAL_TRA -gt 0 ]; then ST="COMPLETE"
        elif [ $TOTAL_ENZ -gt 0 ]; then ST="PARTIAL_ENZYME"
        elif [ $TOTAL_TRA -gt 0 ]; then ST="PARTIAL_TRANSPORTER"
        else ST="ABSENT"; fi
	
        echo -e "$DS\t$SAMP_ID\t$CLEAN_MAG_ID\t$MACRO\t$CAT\t$TYPE\t$TOTAL_ENZ\t$TOTAL_TRA\t$ST" >> "$FINAL_TABLE"
    done
done
