import pandas as pd
import numpy as np
import sys
import os
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Map contrast metadata to DEA results")
    parser.add_argument("--contrast", type=str, help="contrast ID", default="635382")
    parser.add_argument("--dea_results", type=str, help="Path to DEA results file", default="/space/grp/rschwartz/rschwartz/dea-granularity/class_results/dea_results/GSE133283_635382_de_res.tsv")
    parser.add_argument("--dea_meta", type=str, help="Path to DEA metadata file", default="/space/grp/rschwartz/rschwartz/dea-granularity/class_results/meta/GSE133283_meta.tsv")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

        

def main():
	args = parse_arguments()
	# Load DEA results
	dea_results_df = pd.read_csv(args.dea_results, sep="\t")

	# Load DEA metadata
	dea_meta_df = pd.read_csv(args.dea_meta, sep="\t")

	# Filter metadata for the specified contrast
	contrast_meta = dea_meta_df[dea_meta_df['result.ID'] == args.contrast]

	if contrast_meta.empty:
		print(f"No metadata found for contrast ID: {args.contrast}", file=sys.stderr)
		sys.exit(1)

	# Merge DEA results with contrast metadata
	merged_df = dea_results_df.merge(contrast_meta, left_on='contrast', right_on='result.ID', how='left')

	# Save the merged results
	output_file = f"mapped_{os.path.basename(args.dea_results)}"
	merged_df.to_csv(output_file, sep="\t", index=False)
	print(f"Merged DEA results saved to: {output_file}")
 
if __name__ == "__main__":
	main()