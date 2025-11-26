import gemmapy
import argparse
import pandas as pd
import sys

def parse_arguments():
  parser = argparse.ArgumentParser(description='Download DEA results from Gemma and save as TSV files.')
  parser.add_argument('--username', type=str, help='Gemma username', default="raschwar")
  parser.add_argument('--password', type=str, help='Gemma password', default="7nddtt")
  parser.add_argument('--experiment', type=str, help='Gemma experiment accession', default="GSE213364")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

def extract_factors(df):
  # factors column stores data frames where the "value" column has the factor level
  # extrac t the "value" column from each data frame in the factors column
  # use lambda function to check if item is a data frame
  
  levels = df["experimental_factors"].apply(lambda x: x["summary"][0] if isinstance(x, pd.DataFrame) else [])
  subset_factor = df["subset_factor"].apply(lambda x: x["summary"][0] if isinstance(x, pd.DataFrame) else [])

  #baseline_factors = df["baseline_factors"].apply(lambda x: x["summary"][0] if isinstance(x, pd.DataFrame) else [])
  df["experimental_factor"] = levels
  df["cell_type"] = subset_factor
  
  return df

def main():
  args = parse_arguments()
  experiment = args.experiment
  
  # Authenticate if credentials are provided
  if args.username and args.password:
      api = gemmapy.GemmaPy(auth = (args.username, args.password))
  else:
    raise ValueError("Username and password must be provided for authentication.")

  # Fetch DEA meta using gemmapy
  analyses = api.get_dataset_differential_expression_analyses(experiment)
  df = pd.DataFrame(analyses)

  # if analyses is empty, raise error
  if df.isna().all().all():    
    # write empty tsv file
    outpath = f"{experiment}_no_meta.tsv"
    pd.DataFrame().to_csv(outpath, sep="\t", index=False)
    Warning(f"Empty DEA meta for {experiment} written to {outpath}")
    sys.exit(0)
    
  # Convert to DataFrame if possible
  df = extract_factors(df)

  # drop data frame columns
  cols_to_drop = [col for col in df.columns if isinstance(df[col].iloc[0], pd.DataFrame)]
  df = df.drop(columns=cols_to_drop)
  
  output_file = f"{experiment}_meta.tsv"
  df.to_csv(output_file, sep='\t', index=False)
  print(f"DEA meta for {experiment} written to {output_file}")


if __name__ == "__main__":
    main()
