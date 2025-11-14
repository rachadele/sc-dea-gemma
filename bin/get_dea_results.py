import gemmapy
import argparse
import pandas as pd
import os


def parse_arguments():
    parser = argparse.ArgumentParser(description='Download DEA results from Gemma and save as TSV files.')
    parser.add_argument('--username', type=str, help='Gemma username', default="raschwar")
    parser.add_argument('--password', type=str, help='Gemma password', default="7nddtt")
    parser.add_argument('--experiment', type=str, help='Gemma experiment accession', default="GSE280569")
    known_args, _ = parser.parse_known_args()
    return known_args
    
    
def main():
    args = parse_arguments()
    # Authenticate if credentials are provided
    if args.username and args.password:
        api = gemmapy.GemmaPy(auth=(args.username, args.password))
    else:
        raise ValueError("Username and password must be provided for authentication.")

    # Fetch DEA results using gemmapy
    results = api.get_differential_expression_values(args.experiment, readableContrasts=True)

    # Save each result set as a separate TSV file
    for result_id, df in results.items():
        if not isinstance(df, pd.DataFrame):
            df = pd.DataFrame(df)
        filename = f"{args.experiment}_{result_id}.tsv"
        df.to_csv(filename, sep='\t', index=False)
        print(f"Saved: {filename}")

if __name__ == "__main__":
    main()
