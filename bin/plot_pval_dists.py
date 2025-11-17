#!/bin/python
import pandas as pd
import numpy as np
import sys
import os
import argparse
import matplotlib.pyplot as plt
import re
import textwrap

def parse_arguments():
    parser = argparse.ArgumentParser(description="Map contrast metadata to DEA results")
    parser.add_argument("--dea_results_mapped", type=str, help="Path to DEA results file", default="/space/grp/rschwartz/rschwartz/dea-granularity/class_results/dea_meta_mapped/GSE280569/GSE280569_636886_no_ID_mapped.tsv")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args


def check_unique(vals, name):
    unique_vals = vals[~pd.isna(vals)]
    if len(unique_vals) > 1:
        raise ValueError(f"Warning: Multiple unique values found for {name}: {unique_vals}")
    return unique_vals[0] if len(unique_vals) > 0 else 'unknown'

def main():

  args = parse_arguments()

  # Extract GSE, result_id, and contrast_id from the filename
  # Expecting filename like GSE_resultID_contrastID(_mapped).tsv
  base = os.path.basename(args.dea_results_mapped)
  # Remove _mapped and extension
  base = re.sub(r'_mapped\\.tsv$', '', base)
  base = re.sub(r'\\.tsv$', '', base)
  # Split by _
  parts = base.split('_mapped')[0].split('_')
  experiment = parts[0]
  result_id = parts[1]
  contrast_id = "_".join(parts[2:])
  

  # Load DEA results
  dea_results_mapped = pd.read_csv(args.dea_results_mapped, sep="\t")

  # Identify all contrast pvalue columns
  pval_cols = [col for col in dea_results_mapped.columns if re.match(r"contrast(_[\w\d]+)?_pvalue$", col)]

  title_dict = {}

  for column in ["factor_category", "experimental_factor", "cell_type"]:
      unique_value = check_unique(dea_results_mapped[column].unique(), column)
      title_dict[column] = unique_value
  
# make title
  title_parts = [f"{experiment}"] + [f"{value}" for key, value in title_dict.items() if value != "unknown"]
  title_str = " | ".join(title_parts)

  # For every contrast, plot p-value distribution as a single histogram (no grouping)
  for pval_col in pval_cols:
      m = re.match(r"contrast_([\w\d]+)_pvalue$", pval_col)
      this_contrast_id = m.group(1) if m else 'no_ID'

      plt.figure(figsize=(10,4))
      pvals = dea_results_mapped[pval_col].dropna()
      plt.hist(pvals, bins=50, range=(0,1), color='blue', alpha=0.7)
      # Main title: smaller font, concise, wrap if too long
      #title = f"{experiment} | {result_id} | {this_contrast_id}"
      wrapped_title = '\n'.join(textwrap.wrap(title_str, width=80))
      plt.title(wrapped_title, fontsize=10)
      plt.xlabel('p-value')
      plt.ylabel('Frequency')
      fname = f"pval_dist_{experiment}_{result_id}_{this_contrast_id}.png"
      plt.tight_layout()
      plt.savefig(fname)
      plt.close()
      print(f"Saved {fname}")
          
if __name__ == "__main__":
    main()