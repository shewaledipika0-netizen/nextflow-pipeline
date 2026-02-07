# Nextflow Variant Calling Pipeline

## 1. Overview
This repository contains a modular Nextflow DSL2 workflow for basic NGS variant calling. The pipeline follows best practices by keeping all modules inside the workflow file, placing the workflow in the `workflows/` directory, keeping configuration in `nextflow.config`, and providing a Conda environment for reproducibility.

## 2. What the Pipeline Does
The pipeline executes the following stages in order:

- Input FASTQ reads + reference genome  
- Quality Control (FastQC)  
- Read alignment to reference (BWA)  
- Alignment processing (Samtools)  
- Variant calling (BCFtools) → VCF output  

## 3. Repository Structure
The repository is organized in a clean DSL2 format. The `main.nf` file is kept minimal, and the workflow orchestration is stored in `workflows/`.

```text
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

## 4. Requirements / Best Practices Implemented
This pipeline was updated to satisfy the following requirements:

Modules are NOT imported in main.nf (imports are done only inside the workflow).

Workflow is placed inside the workflows/ folder.

No usage of nextflow.enable.dsl = 2 (Nextflow 25+ runs DSL2 by default).

Tool binaries are defined in nextflow.config using absolute paths.

Modules reference tool paths from nextflow.config (no hard-coded commands).

Conda environment .yml file is included in the repository.

README includes all required commands: clone, setup, run, outputs.

## 5. Installation and Setup (Commands)
Follow these commands step-by-step to use this repository on a new system.

### 5.1 Clone the repository
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
Optional: Run again using resume (saves time if already executed):

nextflow run main.nf -resume

## 7. Output Files
All pipeline outputs are stored inside the results/ folder. The pipeline produces:

FastQC HTML reports (quality control)

Aligned BAM file

VCF file containing variants

## 8. Git Commands Used During Development
The following commands were used to update and push the repository changes to GitHub:

git status
git add .
git commit -m "Update workflow structure, config, and README"
git push origin main

## 9. Common Issues and Fixes
Issue: bcftools: command not found
Cause: BCFtools is not installed in the active environment OR its path is not correctly set in nextflow.config.

Fix (recommended):
Install via Conda and confirm using:

which bcftools
Issue: Git warning about embedded repository (submodule warning)
Cause: A folder inside your project was accidentally initialized as a separate git repository.

Fix:

git rm --cached <folder_name>
rm -rf <folder_name>
git add .
git commit -m "Remove nested git repository"
git push origin main

## 10. Conclusion
This project demonstrates a clean modular Nextflow DSL2 implementation of a basic variant calling pipeline. It follows recommended workflow structure, uses configuration-driven tool paths, and provides a reproducible Conda environment. The repository can be cloned and executed successfully on a fresh system using the documented commands.
