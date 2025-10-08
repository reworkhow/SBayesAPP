***Example data and outputs***: https://drive.google.com/drive/folders/1kgnrb8hnYfxZ5rCkh6Ay-nIfKdD6EsY8?usp=drive_link
- `SBayesAPP_nonMPI.jl`: Single-node (no MPI) SBayesAPP script.
- `submit_example.sh`: SLURM script to run SBayesAPP_nonMPI.jl via sbatch on an HPC cluster. 
- `data_preprocess`: Optional preprocessing scripts and a SLURM submission example. Not required for running the main example here because preprocessed inputs are already provided (see links above).

***Quick start (no preprocessing)***
- Download the preprocessed inputs (ST_res and SBayesAPP_input_first10blks) from the link above.
- Edit `submit_example.sh` to point to your Julia path and input/output locations.
- Submit the job: sbatch submit_example.sh

***(Optional) Run the preprocessing example***
If you want to see how preprocessed inputs are generated, use the scripts in data_preprocess/ with the example inputs: https://drive.google.com/drive/folders/1fy-27MNSjTo6YgFnlOWWHbGHD9XhrESM?usp=drive_link.
