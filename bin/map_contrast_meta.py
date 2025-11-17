#!/bin/python
import pandas as pd
import numpy as np
import sys
import os
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Map contrast metadata to DEA results")
    parser.add_argument("--dea_results", type=str, help="Path to DEA results file", default="/space/grp/rschwartz/rschwartz/dea-granularity/class_results/dea_results_files/GSE213364/GSE213364_636937.tsv")
    parser.add_argument("--dea_meta", type=str, help="Path to DEA metadata file", default="/space/grp/rschwartz/rschwartz/dea-granularity/class_results/meta_files/GSE213364_meta.tsv")
    parser.add_argument("--experiment", type=str, help="Experiment accession ID", default="GSE280569")
    parser.add_argument("--result_id", type=str, help="Result ID to filter metadata", default="636927")
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
 
  # map contrast metadata to each de_results_dict df based on dict keys
  for contrast_id, de_results_contrast_df in de_results_dict.items():
      if contrast_id in contrast_meta_dict.keys():
        print(f"Mapping metadata for contrast ID: {contrast_id}")
        contrast_meta = contrast_meta_dict[contrast_id]
        de_results_contrast_df["factor_category"] = contrast_meta.get("factor_category", None)
        de_results_contrast_df["experimental_factor"] = contrast_meta.get("experimental_factor", None)
        de_results_contrast_df["cell_type"] = contrast_meta.get("cell_type", None)
       #breakpoint
        # Save the updated DEA results with contrast metadata in current directory
        outpath = f"{args.experiment}_{result_id}_{contrast_id}_mapped.tsv"
        de_results_contrast_df.to_csv(outpath, sep="\t", index=False)
 
if __name__ == "__main__":
  main()