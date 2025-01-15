Input : VCF file (gzipped and tabixed)

- Reads in VCF
- pulls out the pathogenic/likely_pathogenic status and variant consequence status
- produces a summary table with the headings:
- gene_symbol, gene_id, clnsig_path, missense, lof, in_frame, non_coding, other

Where:
    gene_symbol = gene symbol
    gene_id = gene_id e.g. HGNC or entrez
    clinsig_path = number of variants with a pathogenic or likely_pathogenic clinical significance
    missense = number of missense variants which are pathogenic or likely pathogenic
    lof = number of loss of function (stop gained, frameshift, canonical splicing etc) which are pathogenic or likely pathogenic
    in_frame = number of in frame indels which are pathogenic or likely pathogenic
    non_coding = number of non coding variants (intronic, UTR, upstream etc) that are pathogenic or likely pathogenic
    other = number of variants with other conseqeunces which are  pathogenic or likely pathogenic

Purpose:
    - so we can display to users what type of variants cause disease in each gene e.g. some genes don't have loss of function as a mechanism of disease

example run command:

python clinvar_tabulator.py --vcf ~/Documents/projects/clinvar_project/clinvar.vcf.gz