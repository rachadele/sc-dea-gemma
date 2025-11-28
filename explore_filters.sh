#!/bin/bash
# set to stop on error
set -e
# This script is used to explore different variance filter and min_cells thresholds for DEA analysis

variance_filters=(0.01 0.02 0.05)
min_cells_values=(10 20 50 100)
levels=("class" "subclass")

for var_filter in "${variance_filters[@]}"; do
  for min_cells in "${min_cells_values[@]}"; do
	for level in "${levels[@]}"; do

		echo "Running DEA with variance filter: $var_filter, min_cells: $min_cells, level: $level"
		nextflow run main.nf \
	  --variance_filter $var_filter \
	  --min_cells $min_cells \
	  --level $level \
	  -resume
	  
	  done
	done
done