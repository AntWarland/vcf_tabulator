import argparse
# import csv
# import sys
import pandas as pd
# from pyvariantfilter.family import Family
# from pyvariantfilter.family_member import FamilyMember
# from pyvariantfilter.variant_set import VariantSet
from funcs import tabulate_vcf, summarise_variants
# from pysam import VariantFile

pd.set_option('display.max_rows', 500)

parser = argparse.ArgumentParser(description='Find UPD events in NGS Trio Data')
parser.add_argument('--vcf', type=str, required=True,
				help='The path to the VCF file. Must be bgzipped and tabixed.')

args = parser.parse_args()

vcf = args.vcf

pd.set_option('display.max_rows', 500)
df = tabulate_vcf(vcf)
df.to_csv("./df_all_variants_binaries.csv") 
print("dataframe saved to file.")


# summary_df = summarise_variants(df)