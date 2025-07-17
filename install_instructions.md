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

---

Let me know if you'd like a minimal command-line example for the non-MPI version too.



> 🔧 You must have an **MPI-compatible Julia setup** to run the MPI-version code. However, code provided in the example_data folder is exactly as the MPI-version code but without need to use MPI. 

---

## 📦 Julia Package Setup

1. **Launch Julia**:
   ```bash
   julia

