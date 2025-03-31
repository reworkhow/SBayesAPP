#!/bin/bash

# Parent folder
parent_folder="h2_0.1_pi0.9"

# Output file
output_file="$parent_folder/MT_input_result.csv"
echo "subfolder_name,T1_varg,T1_pi,T2_varg,T2_pi" > "$output_file"

# Loop through subfolders
for subfolder in $parent_folder/h2_trait1.0.5.h2_trait2.0.2_pleioPercent*_sampleSize300000_seed*; do
  if [ -d "$subfolder" ]; then
    t1_varg=$(cat "$subfolder/Trait1/mean_varg_total.txt")
    t1_pi=$(cat "$subfolder/Trait1/mean_pi.txt")
    t2_varg=$(cat "$subfolder/Trait2/mean_varg_total.txt")
    t2_pi=$(cat "$subfolder/Trait2/mean_pi.txt")
    echo "$(basename "$subfolder"),$t1_varg,$t1_pi,$t2_varg,$t2_pi" >> "$output_file"
  fi
done

