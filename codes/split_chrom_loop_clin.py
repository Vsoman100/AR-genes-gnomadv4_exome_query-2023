import pandas as pd
import argparse as ag

##############################
parser = ag.ArgumentParser(description='Enter Date')
# Add arguments
parser.add_argument('--date', help='Date of Pipeline Run')

# Parse the command-line arguments
args = parser.parse_args()

# Access argument values
date = args.date
##############################

# Opens TSV file of all Clinvar autosomal variants that are pathogenic/
# likely pathogenic associated with gene from list
clinvarVariants = open(f"output/parse_clin_out/{date}_selected_autorec_plp_nonconf_tsv2.tsv", "r")

# Adds empty list to be filled with variants from each chromosome
clinChrom = [[] for _ in range(22)]

# Reads through Clinvar variants
for line in clinvarVariants:
    # Splits Clinvar variant info
    clinData = line.split("\t")

    # Excludes header
    if clinData[1] != "CHROM":
        # Adds variant data to corresponding chromosome list
        clinChrom[int(clinData[1])-1].append(line.split('\t'))

clinvarVariants.close()

# Reads through chromosome
for chrom in clinChrom:
    # Reads through variant in chromosome list
    for var in chrom:
        # Removes new line from GNOMAD_ID
        if var[-1].__contains__("\n"):
            var[-1] = var[-1].replace("\n", "")
        # Removes number for each variant
        var.pop(0)
    # Makes dataframe of Clinvar data sorted by chromosome
    tempExtract = pd.DataFrame(chrom, columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'ALLELEID', 'CLNDISDB',
                                                      'CLNHGVS', 'CLNSIG', 'CLNVC', 'MC', 'AF_ESP', 'AF_EXAC', 'AF_TGP',
                                                      'RS', 'CLNSIGCONF', 'CLNDISDBINCL', 'GENEINFO', 'GENE',
                                                      'GNOMAD_ID'])
    # Puts variant data in chromosome file
    tempExtract.to_csv(f"{date}_clin_" + str(chrom[0][0]) + ".tsv", sep="\t")
