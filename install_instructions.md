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

> 🔧 You must have an **MPI-compatible Julia setup**, e.g., with `MPICH` or `OpenMPI` loaded on your HPC.

---

## 📦 Julia Package Setup

1. **Launch Julia**:
   ```bash
   julia

