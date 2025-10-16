# SBayesAPP

SBayesAPP (Summary-data-based Bayesian method leveraging biological Annotations to quantify Pleiotropy and Polygenicity) is a Bayesian framework for dissecting the shared genetic architecture between complex traits using GWAS summary statistics and functional annotations.

This method estimates:
- Annotation-stratified genetic covariance matrix
- Annotation-stratified SNP effect covariance matrix
- Annotation-stratified polygenic proportions 

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
The dataset for T2D and FG (cell type analysis) is provided in https://drive.google.com/drive/folders/1nR_wAj9Hwk1LCrRcFVEr9J-qXiNeaaRp. 

## 📂 Repository Structure

```text
SBayesAPP/
├── Example/                    # demo script to run SBayesAPP on example data (no MPI required)
├── RealData/                   # code for real data analysis used in the manuscript
├── Simulation/                 # code for simulated data analysis used in the manuscript
├── README.md
├── install_instructions.md
└── CITATION.cff
