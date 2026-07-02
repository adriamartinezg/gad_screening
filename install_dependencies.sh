#!/bin/bash

set -e

ENV_NAME="gad_system_env"

echo "=========================================================="
echo "Installing bioinformatic dependencies for the GAD system screening"
echo "=========================================================="

echo "--> Creating Conda environment and installin packages..."
conda create -n $ENV_NAME -y \
    -c conda-forge -c bioconda -c defaults \
    python=3.10 \
    hmmer \
    scikit-learn \
    matplotlib \
    mafft \
    trimal \
    skani \
    gtdbtk \
    prodigal \
    mmseqs2 \
    drep \
    ncbi-datasets

echo "--> Activating the environment to install Python libraries..."
# Usamos este comando para activar conda de forma segura dentro de un script script
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

echo "--> Installing pybiolib via pip3..."
pip3 install --upgrade pybiolib

echo "=========================================================="
echo "Succeed in installation! "
echo "=========================================================="
echo "Please, activate your environment:"
echo "  conda activate $ENV_NAME"
echo ""
echo "To check DeepTMHMM execute:"
echo "  biolib run DTU/DeepTMHMM --help"
echo "=========================================================="
