## Code Explanation

* **`generate_annot_gz.sh`**
  Generates `SimData.1.annot.gz` by running `Rscript annot.gz4sldsc.R`.
  **Output:** `SimData.1.annot.gz`

* **`compute_ldsc_4_sldsc.sh`**
  Computes LD scores using the annotation file `SimData.1.annot.gz`.
  **Output:** LD score files (e.g., `SimData.1.l2.ldscore.gz` + index files)

* **`sldsc.sh`**
  Runs partitioned heritability (binary annotation) with LDSC using the LD scores above.
  **Output:** Heritability partition results (e.g., `.log`)
