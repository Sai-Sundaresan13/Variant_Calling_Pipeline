#!/bin/bash

ENV_NAME="variant_pipeline"  # Change environment name if desired

# Update system
sudo apt-get update && sudo apt-get upgrade -y

# Install basic dependencies
sudo apt-get install -y \
    wget \
    curl \
    unzip \
    zip \
    default-jdk \
    python3 \
    python3-pip \
    build-essential


# Install Conda (if not already installed)
if ! command -v conda &> /dev/null; then
    echo "Installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$HOME/miniconda"
    export PATH="$HOME/miniconda/bin:$PATH"
    echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc
    conda init bash
    source ~/.bashrc
else
    echo "Conda already installed, skipping..."
fi


# Add required channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge


# Create dedicated environment
echo "Creating conda environment: ${ENV_NAME}..."
conda create -n "${ENV_NAME}" -y


# Install all tools into the environment
echo "Installing tools into ${ENV_NAME}..."
conda install -n "${ENV_NAME}" -y -c bioconda -c conda-forge \
    fastqc \
    multiqc \
    trimmomatic \
    hisat2 \
    samtools \
    picard \
    gatk4 \
    snpeff \
    bcftools


# Verify installations
echo ""
echo " Verifying Installations "
conda run -n "${ENV_NAME}" bash << 'EOF'
echo -n "FastQC:       "; fastqc --version 2>&1 | head -1
echo -n "MultiQC:      "; multiqc --version 2>&1 | head -1
echo -n "Trimmomatic:  "; trimmomatic -version 2>&1 | head -1
echo -n "HISAT2:       "; hisat2 --version 2>&1 | head -1
echo -n "SAMtools:     "; samtools --version 2>&1 | head -1
echo -n "Picard:       "; picard MarkDuplicates --version 2>&1 | head -1
echo -n "GATK:         "; gatk --version 2>&1 | head -1
echo -n "BCFtools:     "; bcftools --version 2>&1 | head -1
EOF



echo ""
echo "Installation complete!"
echo ""
echo "To use the pipeline, activate the environment first by running:"
echo ""
echo "    conda activate ${ENV_NAME}"
echo ""
echo "Then run your analysis pipeline as normal."
