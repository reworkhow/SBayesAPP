#!/bin/bash

# Define the parent folder
parent_folder="MTSBayesCC_SNPorder_ShuffleTraits_estGscale_STinput2/"

# Output file
output_file="${parent_folder}/gcov_ratios_concatenated.csv"
echo "Subfolder,gcov_ratio_c1,gcov_ratio_c2" > $output_file

# Loop through each subfolder and process the estGenr_whq.txt file
find "$parent_folder" -type f -name "estGenr_whq.txt" | while read file; do
  # Extract the subfolder name
  subfolder=$(dirname "$file" | xargs basename)
  
  # Read the gcov ratios
  gcov_ratio_c1=$(sed -n '1p' "$file")
  gcov_ratio_c2=$(sed -n '2p' "$file")
  
  # Append to the output file
  echo "$subfolder,$gcov_ratio_c1,$gcov_ratio_c2" >> $output_file
done

