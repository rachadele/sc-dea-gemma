#!/bin/python
import pandas as pd
import numpy as np
import sys
import os
import argparse
import matplotlib.pyplot as plt
import re
import textwrap
import math

def parse_arguments():
      parser = argparse.ArgumentParser(description="Map contrast metadata to DEA results")
      parser.add_argument("--dea_results_mapped", type=str, help="Path to DEA results file", default="/space/grp/rschwartz/rschwartz/dea-granularity/work/8c/2fc36b9f7cb8e6eb5e3fb3cbc02284/GSE280569_oligodendrocyte_precursor_cell_637201_mapped.tsv")
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

  # Improved extraction: GSE, cell type (may have underscores), result_id (numeric, last part before _mapped)
  base = os.path.basename(args.dea_results_mapped)
  base = re.sub(r'_mapped\.tsv$', '', base)
  parts = base.split('_')
  experiment = parts[0]
  # Find the last numeric part as result_id
  result_id_idx = None
  for i in range(len(parts)-1, 0, -1):
    if re.match(r'^\d+$', parts[i]):
      result_id_idx = i
      break
  if result_id_idx is not None and result_id_idx > 1:
    cell_type = '_'.join(parts[1:result_id_idx])
    result_id = parts[result_id_idx]
  else:
    cell_type = 'unknown'
    result_id = parts[-1] if len(parts) > 1 else 'unknown'

  # Load DEA results
  dea_results_mapped = pd.read_csv(args.dea_results_mapped, sep="\t")
  # check if file is empty
  if dea_results_mapped.empty:
      print(f"DEA results file {args.dea_results_mapped} is empty. Exiting.")
      # write empty plot file
      outpath = f"{experiment}_{cell_type}_hist.png"
      plt.figure()
      plt.savefig(outpath)
      print(f"Saved empty plot to {outpath}")
      sys.exit(0)
  # Identify all contrast pvalue columns
  pval_cols = [col for col in dea_results_mapped.columns if re.match(r"contrast(_[\w\d]+)?_pvalue$", col)]

  contrast_ids = []

  for pval_col in pval_cols:
      m = re.match(r"contrast_([\w\d]+)_pvalue$", pval_col)
      if m:
          contrast_ids.append(m.group(1))
  
  if not contrast_ids:
    contrast_ids = ['no_ID']
    
          
  # Build all combinations
  n = len(contrast_ids)
  ncols = min(3, n)
  nrows = math.ceil(n / ncols)
  fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 4*nrows), squeeze=False)
  for idx, contrast_id in enumerate(contrast_ids):
      # extract experimental factors
      factor_value_col = f"contrast_{contrast_id}_experimental_factor"
      exp_factor = dea_results_mapped[factor_value_col].dropna().unique()
      factor_category_col = f"contrast_{contrast_id}_factor_category"
      factor_category = dea_results_mapped[factor_category_col].dropna().unique()
      
      row, col = divmod(idx, ncols)
      ax = axes[row][col]
      # Find the pval column for this contrast
      pval_col = f"contrast_{contrast_id}_pvalue"
      # Subset data for this cell type
      if pval_col in dea_results_mapped.columns:
        pvals = dea_results_mapped[pval_col].dropna()
      else:
        pvals = dea_results_mapped["contrast_pvalue"].dropna()
      ax.hist(pvals, bins=50, range=(0,1), color='blue', alpha=0.7)
      ax.set_title(f"{contrast_id} | {exp_factor} | {factor_category}", fontsize=9)
      ax.set_xlabel('p-value')
      ax.set_ylabel('Frequency')
  # Hide any unused subplots
  for idx in range(len(contrast_ids), nrows*ncols):
      row, col = divmod(idx, ncols)
      fig.delaxes(axes[row][col])
  # Main title
  title_str = f"{experiment} | {cell_type}"
  wrapped_title = '\n'.join(textwrap.wrap(title_str, width=80))
  fig.suptitle(wrapped_title, fontsize=12)
  plt.tight_layout(rect=[0, 0.03, 1, 0.95])
  fname = f"{experiment}_{cell_type}_hist.png"
  plt.savefig(fname)
  plt.close()
  print(f"Saved {fname}")
          
if __name__ == "__main__":
    main()