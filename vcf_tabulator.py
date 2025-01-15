import argparse
import os
import pandas as pd
from funcs import tabulate_vcf, summarise_variants


parser = argparse.ArgumentParser(description='Find UPD events in NGS Trio Data')
parser.add_argument('--vcf', type=str, required=True,
				help='The path to the VCF file. Must be bgzipped and tabixed.')

args = parser.parse_args()

vcf = args.vcf

file_path = "./df_all_variants_binaries.csv"

if os.path.exists(file_path):
    print("All variants csv file exists, loading ...")
    df = pd.read_csv("./df_all_variants_binaries.csv", index_col=0)
else:
    print("All variants csv file does not exist, tabulating VCF")
    df = tabulate_vcf(vcf)
    df.to_csv(file_path) 
    print(f"All variants csv file saved to {file_path}")

summary_path = "./df_path_variants_conseq.csv"
summary_df = summarise_variants(df)
summary_df.to_csv(summary_path)
print(f"Summary table csv file saved to {summary_path}")