import requests
import argparse
import os 
import pandas as pd

def get_fasta_file(url, output_file): 
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_file, 'w') as f:
            f.write(response.text)
        print(f"FASTA file saved to {output_file}")



if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Download a FASTA file from a given URL.")
    arg_parser.add_argument("--kinase_file", type=str, help="The csv containing kinase data")
    arg_parser.add_argument("--init_protein_fasta", type=str, help="The base protein fasta file being compared against")

    args = arg_parser.parse_args()

    df = pd.read_csv(args.kinase_file)

    for id in df['UniprotID']:
        url = f"https://rest.uniprot.org/uniprotkb/{id}.fasta"
        output_file = f"fasta_files/{id}.fasta"
        get_fasta_file(url, output_file)

    for fasta_file in df['UniprotID'].apply(lambda x: f"fasta_files/{x}.fasta"):
        with open(fasta_file, 'a') as f: 
            with open(args.init_protein_fasta, 'r') as init_f: 
                for line in init_f: 
                    f.write(line)


