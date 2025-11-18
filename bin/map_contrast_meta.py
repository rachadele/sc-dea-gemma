#!/bin/python
import pandas as pd
import numpy as np
import sys
import os
import argparse
import re
def parse_arguments():
    parser = argparse.ArgumentParser(description="Map contrast metadata to DEA results")
    parser.add_argument("--dea_results", type=str, help="Path to DEA results file", default="/space/grp/rschwartz/rschwartz/dea-granularity/class_results/dea_results_files/GSE213364/GSE213364_637157.tsv")
    parser.add_argument("--dea_meta", type=str, help="Path to DEA metadata file", default="/space/grp/rschwartz/rschwartz/dea-granularity/class_results/meta_files/GSE213364_meta.tsv")
    parser.add_argument("--experiment", type=str, help="Experiment accession ID", default="GSE213364")
    parser.add_argument("--result_id", type=str, help="Result ID to filter metadata", default="637157")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

def split_de_contrasts(dea_results_df):      
  # extract unique contrast IDs from contrast_{id}_* columns and return a single data frame for each contrast
  contrast_cols = [col for col in dea_results_df.columns if col.startswith("contrast_")]
  contrast_ids = set()
  no_id_cols = []
  for col in contrast_cols:
    parts = col.split("_")
    if len(parts) > 2:
      contrast_id = str(parts[1])  # e.g. contrast_636875_log2fc
      contrast_ids.add(contrast_id)
    else:
      # e.g. contrast_log2fc (no id)
      no_id_cols.append(col)
  contrast_dfs = {}
  # Handle contrast IDs with explicit id
  for contrast_id in contrast_ids:
    contrast_specific_cols = [col for col in contrast_cols if col.startswith(f"contrast_{contrast_id}_")]
    selected_cols = ["Probe", "NCBIid", "gene_ensembl_id", "GeneSymbol", "GeneName", "pvalue", "corrected_pvalue", "rank"] + contrast_specific_cols
    contrast_df = dea_results_df[selected_cols].copy()
    contrast_dfs[str(contrast_id)] = contrast_df
  # Handle columns with no id (map to 'None')
  if no_id_cols:
    selected_cols = ["Probe", "NCBIid", "gene_ensembl_id", "GeneSymbol", "GeneName", "pvalue", "corrected_pvalue", "rank"] + no_id_cols
    contrast_df = dea_results_df[selected_cols].copy()
    contrast_dfs['no_ID'] = contrast_df
  return contrast_dfs

def split_meta(dea_meta_df, result_id):
  # Filter metadata for the specified contrast and result_id
  # change NaNs to None
  # make a string type
  
  dea_meta_df["result_ID"] = dea_meta_df["result_ID"].astype(str)
  dea_meta_filtered = dea_meta_df[dea_meta_df["result_ID"] == result_id]
  # create dict of contrast.id: contrast metadata
  contrast_meta_dict = {}
  for _, row in dea_meta_filtered.iterrows():
    if np.isnan(row["contrast_ID"]):
        contrast_id = "no_ID"
    else:
        contrast_id = str(int(row["contrast_ID"]))
    contrast_meta_dict[contrast_id] = row.to_dict()
  # pop result.ID contrast.id and experiment.ID from each contrast metadata
  for contrast_id in contrast_meta_dict:
    contrast_meta_dict[contrast_id].pop("result_ID", None)
    contrast_meta_dict[contrast_id].pop("experiment_ID", None)
    contrast_meta_dict[contrast_id].pop("contrast_ID", None)
  return contrast_meta_dict

def clean_celltype(ct):
  if pd.isna(ct):
    return "unknown"
  return str(ct).replace(' ', '_').replace('/', '_')

def main():
  args = parse_arguments()
  # Load DEA results
  dea_results_df = pd.read_csv(args.dea_results, sep="\t")
  # Load DEA metadata
  dea_meta_df = pd.read_csv(args.dea_meta, sep="\t")
  # check if columns have no content
  if dea_results_df.isna().all().all():
      Warning("DEA metadata file is empty. Exiting without mapping.")
      # write fake tsv to trick nextflow
      outpath = f"{args.experiment}_{args.result_id}_no_res.tsv"
      pd.DataFrame().to_csv(outpath, sep="\t", index=False)
      sys.exit(0)
  result_id = str(args.result_id)
  
  de_results_dict = split_de_contrasts(dea_results_df)
  contrast_meta_dict = split_meta(dea_meta_df, result_id)
  de_results_combined = pd.DataFrame()
  # map contrast metadata to each de_results_dict df based on dict keys

  for contrast_id, de_results_contrast_df in de_results_dict.items():
    if contrast_id in contrast_meta_dict.keys():
      print(f"Mapping metadata for contrast ID: {contrast_id}")
      contrast_meta = contrast_meta_dict[contrast_id]
      de_results_contrast_df[f"contrast_{contrast_id}_factor_category"] = contrast_meta.get("factor_category", None)
      de_results_contrast_df[f"contrast_{contrast_id}_experimental_factor"] = contrast_meta.get("experimental_factor", None)
      # Clean cell type name for both column and output
      cell_type_val = contrast_meta.get("cell_type", None)
      cell_type_clean = clean_celltype(cell_type_val)
      de_results_contrast_df[f"contrast_{contrast_id}_cell_type"] = cell_type_clean
      
      if de_results_combined.empty:
          de_results_combined = de_results_contrast_df
      else:
          # merge on gene identifiers
          de_results_combined = pd.merge(de_results_combined, de_results_contrast_df, on=["Probe", "NCBIid", "gene_ensembl_id", "GeneSymbol", "GeneName", "pvalue", "corrected_pvalue", "rank"], how="outer")




  # check that cell type is the same across all contrasts
  cell_types = set()
  
  for col in [c for c in de_results_combined.columns if c.endswith("_cell_type")]:
    cts = de_results_combined[col].unique()
    cell_types = set(list(cell_types) + list(cts))
    
  if len(cell_types) > 1:
      raise ValueError(f"Warning: Multiple unique cell types found across contrasts: {cell_types}")
  cell_type = list(cell_types)[0] if cell_types else "unknown"
  de_results_combined.to_csv(f"{args.experiment}_{cell_type}_{result_id}_mapped.tsv", sep="\t", index=False)
  
  
  
if __name__ == "__main__":
  main()