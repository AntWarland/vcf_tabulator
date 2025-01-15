import pandas as pd
from pysam import VariantFile
import numpy as np


def translate_column(dictionary, col_to_search, df, empty_col):
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
    mcs = []
    poses = []

    for rec in vcf_in.fetch():
            
            try:
                gene_symbol = rec.info['GENEINFO'].split(":",2)[0]
                gene_id_pre = rec.info['GENEINFO'].split(":",2)[1]
                gene_id = gene_id_pre.split("|",2)[0]
            except KeyError:
                gene_symbol = "NA"
                gene_id = "NA"
            
            try:
                clnsig = rec.info['CLNSIG'][0]          
            except KeyError:
                clnsig = "NA"
            
            # MC : note maximum number of elements observed in the rec.info['MC'] tuple is 2 #
            try:
                if len(rec.info['MC']) > 1 :
                    mc = str(rec.info['MC'][0]) + " " + str(rec.info['MC'][0])
                    # mc2 = re.findall(r"\|(\S*)$", mc2)[0]
                else:
                    mc = str(rec.info['MC'][0])
            except KeyError:
                mc = np.NaN
                # print(f"{var_id} variant type invalid")

            pos = rec.pos         

            gene_symbols.append(gene_symbol)
            gene_ids.append(gene_id)
            clnsigs.append(clnsig)
            mcs.append(mc)
            poses.append(pos)
            counter += 1
    print(counter)

    df_dict = {'gene_symbol': gene_symbols,
                'gene_id': gene_ids,
                'clnsig': clnsigs,
                'mc': mcs,
                # 'mc2': mc2s,
                'pos': poses}
    
    df = pd.DataFrame(df_dict)

    print(f"All variants = {df.shape[0]}")    
    
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
      

    clnsig_dict = {'Pathogenic': True,
              'Likely_pathogenic': True,
              'Pathogenic/Likely_pathogenic': True,
              'Pathogenic/Likely_pathogenic/Pathogenic': True,
              'Likely_pathogenic/Likely_risk_allele':True,
              'Pathogenic/Pathogenic':True,
              'Pathogenic/Likely_pathogenic/Likely_risk_allele': True,
              'Pathogenic/Likely_risk_allele': True}

    df["clnsig_path"] = False
    translate_column(clnsig_dict, "clnsig", df, "clnsig_path")

    missense_dict = {'missense_variant': True}
    lof_dict = {'nonsense': True,
                'frameshift_variant': True,
                'splice_donor_variant': True,
                'splice_acceptor_variant': True}
    if_indels_dict = {'ininframe_insertion': True,
                      'inframe_deletion': True,
                      'inframe_indel': True}
    non_coding_dict = {'intron_variant': True,
                       '5_prime_UTR_variant': True,
                       '3_prime_UTR_variant': True,
                       'genic_upstream_transcript_variant': True}
    other_dict = {'synonymous_variant': True,
                  'no_sequence_alteration': True,
                  'initiator_codon_variant': True,
                  'non-coding_transcript_variant': True,
                  'stop_lost': True,
                  'genic_downstream_transcript_variant': True,
                  }
    conseq_dict_list = [missense_dict, lof_dict, if_indels_dict, non_coding_dict, other_dict]
    conseqs = ['missense', 'lof','in_frame','non_coding', 'other']

    for conseq, conseq_dict in zip(conseqs, conseq_dict_list):
        print(f"Now starting {conseq} variants")
        df[conseq] = False
        translate_column(conseq_dict, "mc", df, f"{conseq}")
        # print(pd.crosstab(df["clnsig_path"], df[conseq], margins=True))
        # print(pd.crosstab(df["clnsig_path"], df[conseq], margins=True, normalize=True))


    df_consice = df[['gene_symbol', 'gene_id', 'clnsig_path', 'missense', 'lof','in_frame','non_coding', 'other']]
    print(df_consice.head())
    return df_consice


def summarise_variants(variants_df):
      df = variants_df.loc[variants_df['clnsig_path']==True, :]
      grouped_df = df.groupby(['gene_symbol', 'gene_id']).sum()
      print(grouped_df.head())
      return grouped_df