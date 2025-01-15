import argparse
import pandas as pd
from funcs import tabulate_vcf, summarise_variants


pd.set_option('display.max_rows', 500)

parser = argparse.ArgumentParser(description='Find UPD events in NGS Trio Data')
parser.add_argument('--vcf', type=str, required=True,
				help='The path to the VCF file. Must be bgzipped and tabixed.')

args = parser.parse_args()

vcf = args.vcf

pd.set_option('display.max_rows', 500)

### First time you run the script, unhash the following the 3 lines. Hash again for subsequent runs ###
# df = tabulate_vcf(vcf)
# df.to_csv("./df_all_variants_binaries.csv") 
# print("dataframe saved to file.")

df = pd.read_csv("./df_all_variants_binaries.csv", index_col=0)

summary_df = summarise_variants(df)
summary_df.to_csv("./df_path_variants_conseq.csv")