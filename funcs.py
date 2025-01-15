import io
import pandas as pd
from pysam import VariantFile
import re
import numpy as np
import math
# from pyvariantfilter.variant import Variant
# import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
# import scipy


def get_t(dictionary, col_to_search, df, empty_col):
    """
    'Translates' values in a dataframe to a new column according to a given dictionary
    """
    for key in dictionary:
        df.loc[df[col_to_search].str.contains(key, regex=False, na=False), empty_col] = dictionary[key]

def tabulate_vcf(vcf):
    """
    - Read in vcf file, iterate through rows and pull out the relevant info (gene info, clinical significance, consequence of variant, etc.)
    - Store info in lists and construct dataframe using dictionary, format if necessary.
    - Format significance column to create one binary column that reads Pathogenic/Likely_pathogenic or Other
    - Create binary consequence columns i.e. missense: True or False
    - Slice df to remove unneeded columns created during formatting
    - return sliced df
    """
    vcf_in = VariantFile(vcf)
    counter = 0
    gene_symbols = []
    gene_ids = []
    clnsigs = []
    mc1s = []
    mc2s = []
    poses = []

    for rec in vcf_in.fetch():
            
            try:
                gene_symbol = rec.info['GENEINFO'].split(":",2)[0]
                gene_id = rec.info['GENEINFO'].split(":",2)[1]
            except KeyError:
                gene_symbol = "NA"
                gene_id = "NA"
            
            try:
                clnsig = rec.info['CLNSIG'][0]                
            except KeyError:
                clnsig = "NA"
                # print(f"{var_id} clnsig invalid")
            
            # MC : note maximum number of elements observed in the rec.info['MC'] tuple is 2 #
            try:
                # slice the first element in the consequence tuple #
                mc1 = rec.info['MC'][0]
                # filter out noise around the consequence string #
                mc1 = re.findall(r"\|(\S*)$", mc1)[0]

                if len(rec.info['MC']) > 1 :
                    mc2 = rec.info['MC'][1]
                    mc2 = re.findall(r"\|(\S*)$", mc2)[0]
                else:
                    mc2 = "NA"
            except KeyError:
                mc1 = np.NaN
                mc2 = np.NaN
                # print(f"{var_id} variant type invalid")

            pos = rec.pos         

            gene_symbols.append(gene_symbol)
            gene_ids.append(gene_id)
            clnsigs.append(clnsig)
            mc1s.append(mc1)
            mc2s.append(mc2)
            poses.append(pos)
            # print(f"rec.info.CLNSIG = {rec.info['CLNSIG']}")
            counter += 1
    print(counter)

    df_dict = {'gene_symbol': gene_symbols,
                'gene_id': gene_ids,
                'clnsig': clnsigs,
                'mc1': mc1s,
                'mc2': mc2s,
                'pos': poses}
    
    df = pd.DataFrame(df_dict)

    print(f"All variants = {df.shape[0]}")

    df_clnsig = df["clnsig"].str.split('|', expand=True).rename(columns=lambda x: f"clnsig_{x+1}")
        # clnsig split creates four colums # 
    
    # merge df with clnsig columns #
    df = pd.merge(df, df_clnsig, left_index=True, right_index=True)
    
    # mc_unique = pd.unique(np.append(df['mc1'].unique(), df['mc2'].unique()))
    """ unique entries for variant consequences
    ['intron_variant' 'missense_variant' 'synonymous_variant'
    'splice_donor_variant' 'inframe_deletion' 'nonsense' 'frameshift_variant'
    'inframe_insertion' 'splice_acceptor_variant' 'no_sequence_alteration'
    '5_prime_UTR_variant' 'NA' '3_prime_UTR_variant'
    'initiator_codon_variant' 'inframe_indel' 'non-coding_transcript_variant'
    'stop_lost' 'genic_downstream_transcript_variant'
    'genic_upstream_transcript_variant']
    """

    print(f"rows in df after merge = {df.shape[0]}")
    # clnsig_unique = pd.unique(np.concatenate((df['clnsig_1'].unique(), df['clnsig_2'].unique(), df['clnsig_3'].unique(), df['clnsig_4'].unique()), axis=0))
    """ unique entries for clnsig  
    clnsig_unique unique ['Uncertain_significance' 'Likely_benign' 'Benign'
    'Conflicting_classifications_of_pathogenicity' 'Benign/Likely_benign'
    'Pathogenic' 'Likely_pathogenic' 'Pathogenic/Likely_pathogenic'
    'not_provided' 'risk_factor' 'NA' 'Affects'
    'no_classification_for_the_single_variant' 'association' 'drug_response'
    'Uncertain_risk_allele' 'other'
    'no_classifications_from_unflagged_records' 'Likely_risk_allele'
    'Pathogenic/Likely_pathogenic/Pathogenic' 'protective'
    'Uncertain_significance/Uncertain_risk_allele'
    'Likely_pathogenic/Likely_risk_allele' 'association_not_found'
    'Pathogenic/Likely_pathogenic/Likely_risk_allele' 'Pathogenic/Pathogenic'
    'Pathogenic/Likely_risk_allele' 'confers_sensitivity' None]
    """

    t_dict = {'Pathogenic': 'Pathogenic/Likely_pathogenic',
              'Likely_pathogenic': 'Pathogenic/Likely_pathogenic',
              'Pathogenic/Likely_pathogenic': 'Pathogenic/Likely_pathogenic',
              'Pathogenic/Likely_pathogenic/Pathogenic': 'Pathogenic/Likely_pathogenic',
              'Likely_pathogenic/Likely_risk_allele':'Pathogenic/Likely_pathogenic',
              'Pathogenic/Pathogenic':'Pathogenic/Likely_pathogenic',
              'Pathogenic/Likely_pathogenic/Likely_risk_allele': 'Pathogenic/Likely_pathogenic',
              'Pathogenic/Likely_risk_allele': 'Pathogenic/Likely_pathogenic'}

    clnsig_cols = ['clnsig_1', 'clnsig_2','clnsig_3','clnsig_4']
    for col in clnsig_cols:
        get_t(t_dict, f"{col}", df, f"{col}_path")
    
    clnsig_path_cols = ['clnsig_1_path', 'clnsig_2_path','clnsig_3_path','clnsig_4_path']
    
    df["path"] = "Other"

    counter2= 0
    for i, row in df.iterrows():
        for col in clnsig_path_cols:
            if pd.notna(row[col]):
                df.loc[i, "path"] = row[col]
        counter2+=1
        print(f"df['path'] row {counter2} = {df.loc[i, 'path']}")
    
    # Assuming df and clnsig_path_cols are already defined
    # df['path'] = df[clnsig_path_cols].bfill(axis=1).iloc[:, 0]

    # If you need to print the updated rows
    # for i, path in enumerate(df['path']):
    #     print(f"df['path'] row {i+1} = {path}")

    conseq_cols = ["mc1", "mc2"]
    conseqs = ['missense', 'lof','in_frame','non_coding', 'other_consequence']
    
    counter2=0
    for conseq in conseqs:
        print(f"Now starting conseq {conseq}")
        df[conseq] = "False"
        for i, row in df.iterrows():
            for col in conseq_cols:
                if pd.notna(row[col]):
                    df.loc[i, conseq] = "True"
            counter2+=1
            print(f"df['path'] row {counter2} = {df.loc[i, 'path']}")

    
    # for conseq in conseqs:
    #     print(f"Now starting conseq {conseq}")
    #     df[conseq] = df[conseq_cols].notna().any(axis=1).astype(str)
        print(pd.crosstab(df["path"], df[conseq], margins=True))
        print(pd.crosstab(df["path"], df[conseq], margins=True, normalize=True))
    # If you need to print the updated rows
    # for i, row in df.iterrows():
    #     print(f"df['{conseq}'] row {i+1} = {row[conseq]}")


    df_consice = df[['gene_symbol', 'gene_id', 'path', 'clnsig', 'missense', 'lof','in_frame','non_coding', 'other_consequence']]
    print(df_consice.head())
    return df_consice


def summarise_variants(variants_df):
      variants_df["path"] = ""
      grouped_df = variants_df.groupby()