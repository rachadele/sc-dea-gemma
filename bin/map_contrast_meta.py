#!/bin/python
import pandas as pd
import numpy as np
import sys
import os
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Map contrast metadata to DEA results")
    parser.add_argument("--dea_results", type=str, help="Path to DEA results file", default="/space/grp/rschwartz/rschwartz/dea-granularity/work/90/2578ee14a2ba7da16a8e94f02ea5fb/GSE157827_636863_dea_results.tsv")
    parser.add_argument("--dea_meta", type=str, help="Path to DEA metadata file", default="/space/grp/rschwartz/rschwartz/dea-granularity/work/90/2578ee14a2ba7da16a8e94f02ea5fb/GSE157827_meta.tsv")
    parser.add_argument("--experiment", type=str, help="Experiment accession ID", default="GSE157827")
    parser.add_argument("--result_id", type=str, help="Result ID to filter metadata", default="636863")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

def split_de_contrasts(dea_results_df):      
  # extract unique contrast IDs from contrast_{id}_* columns and return a single data frame for each contrast
  contrast_cols = [col for col in dea_results_df.columns if col.startswith("contrast_")]
  contrast_ids = set()
  for col in contrast_cols:
    contrast_id = str(col.split("_")[1])  # ensure string type
    contrast_ids.add(contrast_id)
  contrast_dfs = {}
  for contrast_id in contrast_ids:
    # contrast columns follow the pattern contrast_[id]_*
    contrast_specific_cols = [col for col in contrast_cols if col.startswith(f"contrast_{contrast_id}_")]
    # select Probe, NCBIid, gene_ensembl_id, GeneSymbol, GeneName, pvalue, correctedpvalue, rank,contrast specific columns
    selected_cols = ["Probe", "NCBIid", "gene_ensembl_id", "GeneSymbol", "GeneName", "pvalue", "corrected_pvalue", "rank"] + contrast_specific_cols
    contrast_df = dea_results_df[selected_cols].copy()
    contrast_dfs[str(contrast_id)] = contrast_df  # ensure string key
  return contrast_dfs

def split_meta(dea_meta_df, result_id):
  # Filter metadata for the specified contrast and result_id
  # make a string type
  dea_meta_df["result.ID"] = dea_meta_df["result.ID"].astype(str)
  dea_meta_filtered = dea_meta_df[dea_meta_df["result.ID"] == result_id]
  # create dict of contrast.id: contrast metadata
  contrast_meta_dict = {}
  for _, row in dea_meta_filtered.iterrows():
    contrast_id = str(row["contrast.id"])  # ensure string key
    contrast_meta_dict[contrast_id] = row.to_dict()
  # pop result.ID contrast.id and experiment.ID from each contrast metadata
  for contrast_id in contrast_meta_dict:
    contrast_meta_dict[contrast_id].pop("result.ID", None)
    contrast_meta_dict[contrast_id].pop("experiment.ID", None)
    contrast_meta_dict[contrast_id].pop("contrast.id", None)
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
  # if result not in contrast metadata result.ID, exit
  if result_id not in dea_meta_df["result.ID"].astype(str).values:
      Warning(f"Result ID {result_id} not found in DEA metadata. Exiting without mapping.")
      # write fake tsv to trick nextflow
      outpath = f"{args.experiment}_{args.result_id}_no_meta.tsv"
      pd.DataFrame().to_csv(outpath, sep="\t", index=False)
      sys.exit(0)
      
  contrast_meta_dict = split_meta(dea_meta_df, result_id)
 
  # map contrast metadata to each de_results_dict df based on dict keys
  for contrast_id, de_results_contrast_df in de_results_dict.items():
      if contrast_id in contrast_meta_dict.keys():
        print(f"Mapping metadata for contrast ID: {contrast_id}")
        contrast_meta = contrast_meta_dict[contrast_id]
        for key, value in contrast_meta.items():
          de_results_contrast_df[key] = value
        # Save the updated DEA results with contrast metadata in current directory
        outpath = f"{args.experiment}_{result_id}_{contrast_id}_mapped.tsv"
        de_results_contrast_df.to_csv(outpath, sep="\t", index=False)
 
if __name__ == "__main__":
  main()