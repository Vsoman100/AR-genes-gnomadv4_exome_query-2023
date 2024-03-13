import pandas as pd
import numpy as np
import argparse as ag

##############################
parser = ag.ArgumentParser(description='Enter Date, clinvar vcf filename, and gene list name')
# Add arguments
parser.add_argument('--date', help='Date of Pipeline Run')
parser.add_argument('--clinvar', help='Clinvar vcf file name')
parser.add_argument('--genelist', help='Genelist')

# Parse the command-line arguments
args = parser.parse_args()

# Access argument values
date = args.date
clinvar = args.clinvar
genelist = args.genelist

##############################

clinvar = open(clinvar, 'r') # Downloaded May 31, 2023

list_non_parsed_plp = []

list_of_abbreviated_plp = []

desired_info_fields = ['ALLELEID', 'CLNDISDB', 'CLNHGVS', 'CLNSIG', 'CLNVC', 'MC', 'AF_ESP', 'AF_EXAC', 'AF_TGP', 'RS', 'CLNSIGCONF', 'CLNDISDBINCL', 'GENEINFO']

########################################################

for line in clinvar:

    plp_non_conf = False # pathogenic or likely pathogenic with no conflicts

    # Remove header lines
    if line.startswith("#"):
        continue

    # Get all fields into separate list item
    split_line = line.split("\t")

    # Remove \n from end of info item
    split_line[-1] = split_line[-1].strip("\n")

    # Split info field into component parts
    split_info_field = split_line[-1].split(";")

    # Remove un-split info field
    split_line.pop()

    for info_item in split_info_field:
        split_line.append(info_item)
        # Find if the line is an autosomal recessive plp with no conflicts
    for field in split_line:
        if field[0:7] == "CLNSIG=":
            if ("athogenic" in field) and ("onflict" not in field):
                plp_non_conf = True

        # If plp with no conflicts, add it to the list
    if plp_non_conf:
            list_non_parsed_plp.append(split_line)


print("Done extracting AR non-conflicting PLP to list_non_parsed_plp")

########################################################

# Move through plps in list of plp non conflicting entries
for plp in list_non_parsed_plp:

    # Take CHROM, POS, ID, REF, ALT
    abbreviated_plp = [plp[0], plp[1], plp[2], plp[3], plp[4]]

    # Move through list of desired info fields
    for desired_info_field in desired_info_fields:

        # Initialize variable to track whether the variant contains the field
        field_exists = False

        # Move through the info fields in the current variant
        for existing_info_field in plp:

            # If the current info field is the desired info field
            if existing_info_field.split("=")[0] == desired_info_field:

                # Specify that the field exists
                field_exists = True

                # Add the value of the current info field to the variant
                abbreviated_plp.append(existing_info_field.split("=")[1])

                continue
            # If the desired info field isn't found in the current plp
        if field_exists == False:
            # Add an empty string to the plp
            abbreviated_plp.append("")
    # Exclude X, Y, or MT genes to leave only variants of autosomal genes
    if abbreviated_plp[0] not in ["X", "Y", "MT"]:

        # Append the abbreviated plp variant to the list
        list_of_abbreviated_plp.append(abbreviated_plp)
# Create a dataframe for the abbreviated plp variants
df = pd.DataFrame(list_of_abbreviated_plp, columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'ALLELEID', 'CLNDISDB', 'CLNHGVS', 'CLNSIG', 'CLNVC', 'MC', 'AF_ESP', 'AF_EXAC', 'AF_TGP', 'RS', 'CLNSIGCONF', 'CLNDISDBINCL', 'GENEINFO'])

# Convert the dataframe of plp variants to a csv
df.to_csv(f"{date}_output_autorec_plp_nonconf_2.csv")

print("Done generating df of all abbreviated plps")

####################################################

# Get a list of the genes from Guo Table S2 with 4 additional desired genes
genes_list = open(genelist, 'r')
list_of_desired_genes = genes_list.read().split("\n")
print("Done opening list of genes of interest")

####################################################

list_of_extracted_rows = []
list_of_unfound = []

# Iterate through genes in the list from the file of desired genes
for gene_of_interest in list_of_desired_genes:

    # Set current gene to unfound
    gene_been_found = False
 # Update all of the genes according to names in Clinvar file.
    # Note that all of the genes are the same, just labeled differently
    # as verified using genecards.org and gnomAD
    if gene_of_interest == "G6PC":
        gene_of_interest = "G6PC1"
    elif gene_of_interest == "MUT":
        gene_of_interest = "MMUT"
    elif gene_of_interest == "GARS":
        gene_of_interest = "GARS1"
    elif gene_of_interest == "C5orf42":
        gene_of_interest = "CPLANE1"
    elif gene_of_interest == "WISP3":
        gene_of_interest = "CCN6"

    # Iterate through the abbreviated plps
    for row in list_of_abbreviated_plp:

        # Split the different genes in the GENEINFO field
        row_genes_with_number = row[-1].split("|")

        # Create list to hold the names of the genes without their numerical ID
        row_genes_without_number = []

        # Iterate through each of the genes (including numerical ID)
        for row_gene in row_genes_with_number:
            # Add the gene name without numerical ID to the list
            row_genes_without_number.append(row_gene.split(":")[0])

        # If the current gene is present in the list of genes without num ID
        if row_genes_without_number != [""] and gene_of_interest in row_genes_without_number:
            # Set the current gene to found

            gene_been_found = True

            # Add the name of the gene that this allele was selected based on
            # to the end of the row for that allele
            row.append(gene_of_interest)

            # Correct for the issue where some alleles had the gene name added
            # to the end of the row twice
            if len(row) == 20:
                row.pop()

            # Add the ID to cross with gnomAD to the selected row
            row.append(str(row[0]) + "-" + str(row[1]) + "-" + row[3] + "-" + row[4])

            # Add this new row to the list of selected variants
            list_of_extracted_rows.append(row)

    # If the gene has not been found, add it to the list of unfound genes
    if gene_been_found == False:
        list_of_unfound.append(gene_of_interest)

# Create a dataframe to hold the alleles for the selected genes
extracted_rows_df = pd.DataFrame(list_of_extracted_rows, columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'ALLELEID', 'CLNDISDB', 'CLNHGVS', 'CLNSIG', 'CLNVC', 'MC', 'AF_ESP', 'AF_EXAC', 'AF_TGP', 'RS', 'CLNSIGCONF', 'CLNDISDBINCL', 'GENEINFO', 'GENE', 'GNOMAD_ID'])

extracted_rows_df.to_csv(f"{date}_selected_autorec_plp_nonconf_2.csv")
extracted_rows_df.to_csv(f"{date}_selected_autorec_plp_nonconf_tsv2.tsv", sep="\t")

print("Done generating df abbreviated plps for genes of interest")

####################################################

# Create a list for all of the genes represented in the list of alleles
list_of_extracted_genes = sorted(extracted_rows_df["GENE"].unique())

# Create a comma-separated list of these genes in string format
str_list_of_extracted_genes = ", ".join(list_of_extracted_genes)

# Write this gene list to a text file
final_list_genes = open(f'{date}_final_gene_list2.txt', 'w')
final_list_genes.write(str_list_of_extracted_genes)

####################################################

clinvar.close()
genes_list.close()
final_list_genes.close()
