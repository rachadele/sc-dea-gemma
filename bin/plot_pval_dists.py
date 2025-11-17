#!/bin/python
import pandas as pd
import numpy as np
import sys
import os
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Map contrast metadata to DEA results")
    parser.add_argument("--dea_results_mapped", type=str, help="Path to DEA results file", default="/space/grp/rschwartz/rschwartz/dea-granularity/work/bf/f12e751e72af5cc120a2508013480c/GSE280569_636866.tsv")
    parser.add_argument("--experiment", type=str, help="Experiment accession ID", default="GSE280569")
    parser.add_argument("--result_id", type=str, help="Result ID to filter metadata", default="636866")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
    
def split_de_contrasts(dea_results_df):
    # split by experimental_factor, cell_type

 
def main():
	args = parse_arguments()

	# Load DEA results
	dea_results_mapped = pd.read_csv(args.dea_results_mapped, sep="\t")

	
if __name__ == "__main__":
    main()