# SBayesAPP

SBayesAPP (Summary-data-based Bayesian method leveraging biological Annotations to quantify Pleiotropy and Polygenicity) is a Bayesian framework for dissecting the shared genetic architecture between complex traits using GWAS summary statistics and functional annotations.

This method estimates:
- Annotation-stratified genetic covariance matrix
- Annotation-stratified SNP effect covariance matrix
- Annotation-stratified Polygenic proportions 

SBayesAPP enables researchers to identify whether coheritability enrichment arises from many shared variants with weak effects or few variants with strong pleiotropic effects.

---

## 📖 Overview

- 🔍 Distinguishes between pleiotropy, shared polygenicity, and LD-driven genetic correlation
- 📊 Supports genome-wide and annotation-stratified inference
- 🧬 Integrates functional annotations (e.g., single-cell expression data)
- 🧪 Benchmarking via simulation and real data (e.g., T2D and fasting glucose)

---

## 🛠 Installation

SBayesAPP is written in Julia and uses some Python and R scripts for preprocessing and visualization.

### Dependencies
- Julia 1.8.x
- Python ≥3.8
- R ≥4.1
- Required packages are listed in [`install_instructions.md`](./install_instructions.md)

---

## 📂 Repository Structure

```text
SBayesAPP/
├── STSBayesC_code/             # Main Julia code for SBayesAPP
├── Simulation/                 # Scripts and data for simulation studies
├── example_data/              # Example input files
├── results/                   # Example output files
├── scripts/                   # Pipeline and preprocessing scripts
├── doc/                       # Documentation and figures
├── README.md
├── install_instructions.md
└── CITATION.cff
