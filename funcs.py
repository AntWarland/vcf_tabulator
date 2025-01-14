import io
import pandas as pd
from pysam import VariantFile
import re
# from pyvariantfilter.variant import Variant
# import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
# import scipy


def tabulate_vcf(vcf):

    vcf_in = VariantFile(vcf)
    counter = 0
    gene_symbols = []
    gene_ids = []
    clnsigs = []
    mc1s = []
    mc2s = []
    poses = []
    paths = []

    for rec in vcf_in.fetch():
        # if counter < 1000:
            var_id = rec.id
            
            try:
                gene_symbol = rec.info['GENEINFO'].split(":",2)[0]
                gene_id = rec.info['GENEINFO'].split(":",2)[1]
            except KeyError:
                gene_symbol = "NA"
                gene_id = "NA"
                # print(f"{var_id} geneinfo invalid")
            
            try:
                clnsig = rec.info['CLNSIG'][0]                
            except KeyError:
                clnsig = "NA"
                # print(f"{var_id} clnsig invalid")

            try:
                mc1 = rec.info['MC'][0]
                mc1 = re.findall(r"\|(\S*)$", mc1)[0]
                if len(rec.info['MC']) > 1 :
                    mc2 = rec.info['MC'][1]
                    mc2 = re.findall(r"\|(\S*)$", mc2)[0]
                else:
                    mc2 = "NA"
            except KeyError:
                mc1 = "NA"
                mc2 = "NA"
                # print(f"{var_id} variant type invalid")

            pos = rec.pos         
            # if clnsig == "L"
            # path = ""


            gene_symbols.append(gene_symbol)
            gene_ids.append(gene_id)
            clnsigs.append(clnsig)
            mc1s.append(mc1)
            mc2s.append(mc2)
            poses.append(pos)
            # print(f"rec.info.CLNSIG = {rec.info['CLNSIG']}")
        # counter += 1
    big_dict = {'gene_symbol': gene_symbols,
                'gene_id': gene_ids,
                'clnsig': clnsigs,
                'mc1': mc1s,
                'mc2': mc2s,
                'pos': poses}
    
    df = pd.DataFrame(big_dict)
    # print(df['clnsig'].unique())
    df_clnsig = df["clnsig"].str.split('|', expand=True).rename(columns=lambda x: f"clnsig_{x+1}")
    print(df.head())
    print(df_clnsig.head()) 
    df = pd.merge(df, df_clnsig, left_index=True, right_index=True)
    return df


def summarise_variants(variants_df):
      variants_df["path"] = ""


      grouped_df = variants_df.groupby()