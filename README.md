# SBayesAPP

SBayesAPP (Summary-data-based Bayesian method leveraging biological Annotations to quantify Pleiotropy and Polygenicity) is a Bayesian framework for dissecting the shared genetic architecture between complex traits using GWAS summary statistics and functional annotations.

This method estimates:
- Annotation-stratified genetic covariance matrix
- Annotation-stratified SNP effect covariance matrix
- Annotation-stratified Polygenic proportions 

SBayesAPP enables researchers to identify whether coheritability enrichment arises from many shared variants with weak effects or few variants with strong pleiotropic effects.

---

## 🛠 Installation

SBayesAPP is written in Julia and uses some R scripts for preprocessing and visualization.

### Dependencies
- Julia 1.11.x
- R ≥4.1
- Required packages are listed in [`install_instructions.md`](./install_instructions.md)

---

## Data 
The total data for T2D and FG (cell type analysis) is provided in https://drive.google.com/drive/folders/1nR_wAj9Hwk1LCrRcFVEr9J-qXiNeaaRp. 

## 📂 Repository Structure

```text
SBayesAPP/
├── SBayesAPP_code/             # Main Julia code for SBayesAPP
├── Simulation/                 # Scripts for simulation studies
├── example_data/              # Example input files
├── results/                   # Example output files
├── scripts/                   # Pipeline and preprocessing scripts
├── doc/                       # Documentation and figures
├── README.md
├── install_instructions.md
└── CITATION.cff
