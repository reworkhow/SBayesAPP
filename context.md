I have a repo here. The current repo structure centers on the Julia package under `src/` and the package-backed non-MPI launcher `scripts/run_nonmpi.jl`. The old MPI path and the legacy non-MPI source entrypoint have been removed.

This repo also contains `script/` where `run.sh` provides a convenience wrapper for the non-MPI package runner, and additional shell helpers can target longer cluster runs. The example inputs live under `example/`, including `SBayesAPP_input_first10blks/`, `ST_res/`, and `annotation_df.txt`. The shell wrapper reads from the example data and writes to an example output directory. Please note that I just construct this repo, so there might be a chance that the relative path is not accurate. Please don't directly refer to the path included in any files.

