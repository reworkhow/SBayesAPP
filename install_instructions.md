# 🛠️ Installation Instructions for SBayesAPP

SBayesAPP is currently implemented in **Julia** and designed to run efficiently on high-performance computing (HPC) environments using **MPI** (Message Passing Interface). A future version with **OpenMP support** will be released to simplify setup and improve accessibility for single-node or multi-threaded users.

---

## 📌 Requirements

### 🧬 Software Dependencies

| Tool          | Version      | Notes                                     |
|---------------|--------------|-------------------------------------------|
| **Julia**     | 1.11.x        | Tested with 1.11.4                       |
| **MPI**       | OpenMPI or MPICH | Used for distributed parallelism      |
| **R**         | ≥4.1         | (Optional) for post-analysis visualization|

SBayesAPP is currently implemented using **MPI-compatible Julia** for distributed parallel computing on HPC systems.
> 🔧 **Note**: You must have an **MPI-compatible Julia setup** to run the full MPI version of SBayesAPP.
> However, for ease of testing and smaller-scale use, we provide a non-MPI equivalent of the core code in the `example_data/` folder. This code is functionally identical but does **not** require MPI and can be run directly in a single Julia session. This allows users to explore and validate the method locally before deploying to HPC environments.

## 📦 Julia Package Setup

1. **Launch Julia**:
   ```bash
   julia

2. **Enter package manager** (press `]` key), then run:

   ```julia
   add CSV
   add DataFrames
   add DelimitedFiles
   add Distributions
   add LinearAlgebra
   add ProgressMeter
   add Random
   add Statistics
   add Dates
   add JLD2
   ```

3. If your are using MPI implementation, also add:

   ```julia
   add MPI
   ```
---

## 📦 Configure MPI on server 
Since you’re on a server, it’s likely that an MPI implementation like OpenMPI or MPICH is already installed. You can check this by running:

   You may need to set the environment variable before building:

   ```bash
   export JULIA_MPI_BINARY=system
   julia -e 'using Pkg; Pkg.build("MPI")'
   ```

---

## 🚀 Running the MPI Version

Make sure you are on an HPC system with `mpirun` or `srun` available.

Example run command:

```bash
mpirun -np 20 julia --project sbayesapp_mpi.jl \
    --sumstats1 trait1.txt \
    --sumstats2 trait2.txt \
    --annot annotation_matrix.txt \
    --ld_dir path_to_ld_blocks/ \
    --out results/
```

You may also submit this using `sbatch` on SLURM-based systems.

---

## 📦 Optional: Python & R Environment

**Python** (for munge\_sumstats):

```bash
pip install pandas numpy scipy
```

**R** (for visualization):

```r
install.packages(c("ggplot2", "reshape2", "dplyr", "pheatmap"))
```

---

## 🔄 Future Update: OpenMP Version

We are actively developing an **OpenMP-compatible** version of SBayesAPP. This upcoming release will:

* Eliminate the need for MPI setup
* Support multi-threading within a single machine
* Be more accessible to users without HPC environments

Stay tuned in the [GitHub Issues](https://github.com/reworkhow/S-MT-Bayes/issues) or watch this repository for updates.

---

## 💬 Help

If you encounter issues during installation or setup:

* Check open and closed [Issues](https://github.com/reworkhow/S-MT-Bayes/issues)
* Or contact the authors listed in the [README](./README.md)

```

---

Let me know if you're using a specific cluster (e.g., SLURM, PBS) and I can tailor the MPI `sbatch` submission example. I can also help prepare a Julia `Project.toml` for package reproducibility.
```

