
# Nextflow Variant Calling Pipeline

## 1. Overview
This repository contains a modular Nextflow DSL2 workflow for basic NGS variant calling. The pipeline follows best practices by keeping all modules inside the workflow file, placing the workflow in the workflows/ directory, keeping configuration in nextflow.config, and providing a Conda environment for reproducibility.

## 2. What the Pipeline Does
The pipeline executes the following stages in order: input FASTQ reads and reference genome, quality control using FastQC, read alignment to the reference using BWA, alignment processing using Samtools, and variant calling using BCFtools to generate a final VCF output.

## 3. Repository Structure
The repository is organized in a clean DSL2 format. The main.nf file is kept minimal, and the workflow orchestration is stored in workflows/.


nextflow-pipeline/

├── main.nf

├── nextflow.config

├── README.md

├── workflows/

│   └── workflow.nf

├── modules/

│   ├── download.nf

│   ├── fastqc.nf

│   ├── align.nf

│   └── variant.nf

├── data/

├── reference/

└── results/

## 4. Requirements 
This pipeline was updated to satisfy required workflow structure and reproducibility rules. Modules are not imported in main.nf; instead, module imports are done only inside the workflow file. The workflow is placed inside the workflows/ folder. The pipeline does not mention nextflow.enable.dsl = 2 because Nextflow 25+ runs DSL2 by default. Tool binaries are defined in nextflow.config using absolute paths. Modules reference tool paths from nextflow.config, avoiding hard-coded commands. The Conda environment YAML file is included in the repository. The README provides required commands such as clone, setup, run, and output locations.

## 5. Installation and Setup (Commands)
Follow these commands step-by-step to use this repository on a new system.

### 5.1 Clone the repository: 

git clone https://github.com/shewaledipika0-netizen/nextflow-pipeline.git

cd nextflow-pipeline

### 5.2 Create the Conda environment

conda env create -f environment.yml

conda activate nextflow_variant_env

Note: If you are using mamba, you can replace conda with mamba for faster installation.

### 5.3 Verify tools are installed correctly

which nextflow

which fastqc

which bwa

which samtools

which bcftools

## 6. Running the Pipeline (Commands)
Run the workflow using the main entry script. This will automatically call the workflow stored inside workflows/.

nextflow run main.nf

nextflow run main.nf -resume

## 7. Output Files
All pipeline outputs are stored inside the results/ folder. The pipeline generates FastQC HTML reports, aligned BAM output, and a VCF file containing detected variants.

## 8. Git Commands Used During Development
The following commands were used to update and push repository changes to GitHub:

git status

git add .

git commit -m "Update workflow structure, config, and README"

git push origin main

## 9. Common Issues and Fixes

Issue 1: bcftools: command not found

Cause: BCFtools is not installed in the active environment or its path is not correctly set in nextflow.config.

Fix: Install via Conda and confirm the binary path:

which bcftools

Issue 2: Git warning about embedded repository (submodule warning)

Cause: A folder inside the project was accidentally initialized as a separate Git repository.

Fix: Remove the nested repository and push again:

git rm --cached <folder_name>

rm -rf <folder_name>

git add .

git commit -m "Remove nested git repository"

git push origin main

## 10. Conclusion
This project demonstrates a clean modular Nextflow DSL2 implementation of a basic variant calling pipeline. It follows recommended workflow structure, uses configuration-driven tool paths, and provides a reproducible Conda environment. The repository can be cloned and executed successfully on a fresh system using the documented commands.
