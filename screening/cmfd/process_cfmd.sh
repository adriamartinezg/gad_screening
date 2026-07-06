#!/bin/bash

# --- CONFIGURACIÓN ---
MODELS_DIR="model"
THRESHOLDS="thresholds.csv"
RESULTS_DIR="results_cfmd"
METADATA_URL="https://github.com/SegataLab/cFMD/blob/86872566a49ebff286291b7102b6c1666bfed1ff/cFMD_metadata.tsv"

DATASETS=('LiZ_2018' 'LiZ_2019' 'LordanR_2019' 'MASTER_WP4_CSIC_1' 'MASTER_WP4_CSIC_2' 'MASTER_WP4_FFoQSI_1' 'MASTER_WP4_FFoQSI_2' 'MASTER_WP4_FFoQSI_3' 'MASTER_WP4_MATIS_1'
'MASTER_WP4_TEAGASC_1' 'MASTER_WP4_UNINA_1' 'MASTER_WP4_UNINA_2' 'McHughAJ_2020' 'MilaniC_2019' 'MortensenS_xxxx' 'PasolliE_2020' 'PatroJN_2016' 'PfeferT_xxxx' 'PorcellatoD_2016' 'PothakosV_2020' 'QuigleyL_2016' 'RippF_2014'
'SalvettiE_2016' 'SomervilleV_2019' 'SternesPR_2017' 'SulaimanJ_2014' 'VerceM_2019' 'WalshAM_2016' 'WalshAM_2017' 'WalshAM_2020' 'WalshL_xxxx' 'WolfeBE_2014' 'XieM_2019' 'YaoG_2017' 'YasirM_2020' 'YulandiA_2020' 'ZhaoCC_2020' 'ShangpliangH_2023_a' 'ShangpliangH_2023_b' 'KharnaiorP_2023' 'YouL_2022' 'LimaC_2020' 'DecadtH_2024' 'YasirM_2022' 'FontanaF_2023' 'FranciosaI_2021'
'MotaGutierrezJ_2021' 'SaakC_2023' 'YangC_2021' 'LopezSanchezR_2023' 'GonzalezOrozcoB_2023' 'SalgadoTS_2021' 'FalardeauJ_2023' 'SequinoG_2024_a' 'SequinoG_2024_b' 'MagliuloR_2024'
'YapM_2020' 'OlgaP_2019' 'QuijadaN_2022' 'YuY_2022' 'TomarS_2023' 'AlmeidaO_2020' 'CM_UNINA_FFOOD')

mkdir -p "$RESULTS_DIR/faa_all"
mkdir -p "$RESULTS_DIR/hits"
mkdir -p "$RESULTS_DIR/summary"

echo "### Downloading cFMD Metadata ###"
wget "$METADATA_URL" -O "$RESULTS_DIR/cFMD_metadata.tsv"

for DS in "${DATASETS[@]}"; do
    echo "==========================================="
    echo "### PROCESSING DATASET: $DS ###"
    echo "==========================================="

    echo "### Downloading dataset's MAGs ###"
    # 1. Descargar (Usando el script proporcionado)
    bash scripts/download_mags.sh "$DS"

    # Prokka
    echo "### Step 1: Running Prokka on MAGs ###"
    # 1. Descomprimimos (igual que antes)
    bzip2 -d "${DS}_mags"/*.bz2 2>/dev/null

    # 2. Bucle de predicción
    shopt -s nullglob
    for MAG in "${DS}_mags"/*.{fna,fa,fasta}; do
        
        # Verificamos que el archivo tenga contenido
        [ -s "$MAG" ] || continue

        # Sacamos el ID limpio para el nombre del archivo
        FULL_NAME=$(basename "$MAG")
        MAG_ID="${FULL_NAME%.*}"

        echo "### Running Prodigal on: $MAG_ID ###"

        # 3. Lanzamos Prodigal (Sustituye a Prokka)
        # -i: entrada (nucleótidos)
        # -a: salida de proteínas (aminoácidos)
        # -p meta: optimizado para MAGs/fragmentos
        # -q: modo silencioso
        prodigal -i "$MAG" \
                 -a "$RESULTS_DIR/faa_all/${MAG_ID}.faa" \
                 -p meta \
                 -q

        # 4. Verificación simple
        if [ -f "$RESULTS_DIR/faa_all/${MAG_ID}.faa" ]; then
             echo "   [OK] Proteínas generadas."
        else
             echo "   [ERROR] Algo falló con $MAG_ID"
        fi
    done

    # 3. Clasificador HMM
    echo "### Running GAD Classifier ###"
    # Concatenamos todas las proteínas del dataset para pasar el clasificador una sola vez
    cat "$RESULTS_DIR/faa_all"/*.faa > "${DS}_temp_all.faa"
    
    bash scripts/classifier.sh -m "$MODELS_DIR" -t "$THRESHOLDS" \
         -r "$RESULTS_DIR/hits/$DS" -s "${DS}_temp_all.faa"

    # 4. Limpieza
    echo "### Cleaning up raw genomic data ###"
    rm -rf "${DS}_mags"
    rm "${DS}_temp_all.faa"
    # Opcional: rm "$RESULTS_DIR/faa_all"/*.faa 
    # (Si ya tienes los resultados de HMMer, no necesitas los .faa de 2GB)
done

echo "==========================================="
echo "### ALL DATASETS PROCESSED ###"
echo "==========================================="
