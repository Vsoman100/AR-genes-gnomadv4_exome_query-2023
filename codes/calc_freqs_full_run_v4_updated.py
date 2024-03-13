import pandas as pd
import sys

pop_file = sys.argv[3]
pops = open(pop_file, "r")

desired_info_fields = ["AC", "AN", "AF", "grpmax"]
num_pops = 0
pop_list = []

# Allele count, allele number, allele frequency, and number of homozygotes for each population
for pop in pops:
    pop = pop.strip()
    pop_list.append(pop)
    #print(pop)
    num_pops += 1
    desired_info_fields.append("AC_"+pop)
    desired_info_fields.append("AN_"+pop)
    desired_info_fields.append("AF_"+pop)
    desired_info_fields.append("nhomalt_"+pop)

desired_info_fields.append("AC_grpmax")
desired_info_fields.append("AN_grpmax")
desired_info_fields.append("AF_grpmax")
desired_info_fields.append("nhomalt_grpmax")

pops.close()


# Finds input file of extracted gnomAD data by chromosome
input_file = sys.argv[1]
# Chromosome number
chrom = sys.argv[2]
# Opens TSV gnomAD data for reading
gnomadVariants = open(input_file, "r")

# Initializes lists for later use in storing data
combinedOut = []
intermediateOut =[]

# Positions for each population corresponding for alelle count, allele number, and number of homozygotes for VCF calculations
initial_triplet = [9,10,12]
VCFvariables = [initial_triplet]
for i in range(1, num_pops):
    x = VCFvariables[-1][-1]
    new_triplet = [x+1, x+2, x+4]
    VCFvariables.append(new_triplet)
    

# Function that calculates VCF for inputted variant and adds to intermediate list
def calcVCF(var):
    for vars in VCFvariables:
        if (int(var[vars[1]])) != 0:
            intermediateOut.append((int(var[vars[0]])-int(var[vars[2]]))/(0.5*int(var[vars[1]])))
        else: 
            intermediateOut.append(0)

# Loops through variant in gnomAD file
for var in gnomadVariants:
    # Splits data by tab
    varDelim = var.split("\t")
    # Replace 'chrNum' with just 'Num' of chromosome
    varDelim[0] = varDelim[0].replace("chr", "")
    # Skips variants that don't have extracted gnomAD data
    if varDelim[3].strip() == "NA":
        continue
    intermediateOut = []
    # Adds chromosome, position, ref seq, alt seq, and gene to intermediate list
    for i in range(0, 5):
        intermediateOut.append(varDelim[i])
    # Splits gnomAD data by semicolon
    gnomadExtract = varDelim[5].split(";")
    # Loops through desired data fields
    for field in desired_info_fields:
        infoFound = False
        # Nested loop that reads through all provided gnomAD fields
        for item in gnomadExtract:
            # Separates field label from data
            fieldExtract = item.split("=")
            # If desired field is found, add data to intermediate list
            if fieldExtract[0] == field:
                infoFound = True
                intermediateOut.append(fieldExtract[1])
        # If desired field is not found, NA is added to intermediate list
        if infoFound is False:
            intermediateOut.append("NA")
    # VCF for all populations are calculated and added to list
    #print(intermediateOut)
    calcVCF(intermediateOut)
    # Variant's data is added to overall list
    combinedOut.append(intermediateOut)

gnomadVariants.close()
# Dataframe is created with variant data

columns_vcf = ['CHROM', 'POS', 'GENE', "POS", "ALT", "AC", "AN", "AF", "grpmax"]
for pop in pop_list:
    columns_vcf.append("AC_"+pop)
    columns_vcf.append("AN_"+pop)
    columns_vcf.append("AF_"+pop)
    columns_vcf.append("nhomalt_"+pop)
    
columns_vcf.append("AC_grpmax")
columns_vcf.append("AN_grpmax")
columns_vcf.append("AF_grpmax")
columns_vcf.append("nhomalt_grpmax")

for pop in pop_list:
    columns_vcf.append("VCF_"+pop)
    

finalExtractdf = pd.DataFrame(combinedOut, columns=columns_vcf)
# TSV file is made with variant data 
finalExtractdf.to_csv(chrom + "_VCF.tsv", sep="\t")

##############################################

# Dictionary is created
gene_gcf_dict = {}

vcf_cols = len(columns_vcf) - num_pops
# Loops through all populations
for i in range (0,num_pops):
    # Nested loop reads through variant
    for var in combinedOut:
        # Only proceeds if VCF is nonzero
        # change from 53 to 49 bc "oth" is gone
        if float(var[vcf_cols+i]) != 0.0:
            # If gene is in dictionary, population (organized by index) has 1-VCF multiplied by current value and adds 1 to variant count
            if var[2] in gene_gcf_dict:
                gene_gcf_dict[var[2]][2*i] *= (1-float(var[vcf_cols+i]))
                gene_gcf_dict[var[2]][1+(2*i)] += 1
            # If gene is not in dictionary, gene values are initialized with 1 and 0 for each population (1 so multiplication by 1-VCF is accurate
            # and 0 because there are 0 contributing varaints to begin with
            # remove 1 set to get rid of oth
            else:
                gene_gcf_dict[var[2]] = [1, 0] * num_pops
                gene_gcf_dict[var[2]][2 * i] *= (1-float(var[vcf_cols+i]))
                gene_gcf_dict[var[2]][1 + (2 * i)] += 1

intermediateOut = []
combinedOut = []
# List of populations as they appear in this data
# remove oth
gcf_index = pop_list
# Loops through variants in dictionary
for key, value in gene_gcf_dict.items():
    # Initializes list
    intermediateOut = []
    # Begins with maximum GCF as 0
    maxCount = [0, 0, '']
    # Adds variant data to list
    intermediateOut.append(key)
    # Loops through VCF data for each variant
    for index, gcf in enumerate(value):
        # Only looks at 1-VCF values (not contributing varaints)
        if index % 2 == 0:
            # Adds final GCF calculation to list 
            intermediateOut.append(1-float(gcf))
            # Finds maximum GCF information
            if (1-float(gcf) > maxCount[0]):
                maxCount[0] = 1-float(gcf)
                maxCount[1] = value[index+1]
                maxCount[2] = gcf_index[int((index/2))]
        # Contributing variant counts are added without being changed
        else:
            intermediateOut.append(gcf)
    # Adds maxium GCF data to list
    for element in maxCount:
        intermediateOut.append(element)
    # Makes final list containing GCF for each population 
    combinedOut.append(intermediateOut)


columns_gcf = ["Gene"]
for pop in pop_list:
    columns_gcf.append("GCF_"+pop)
    columns_gcf.append("contributing_vars")
columns_gcf.append("max_GCF")
columns_gcf.append("contributing_vars")
columns_gcf.append("max_pop")


finalGCF = pd.DataFrame(combinedOut, columns=columns_gcf)
# Makes TSV file
finalGCF.to_csv(chrom + "_GCF.tsv", sep="\t")


