
# Loads tabix tool
module load htslib/1.13
module load bcftools

basepath_exome="$2"
mkdir -p OUT/info_output

# Loops through all chromosomes
for chrom in {1..22}
do

echo "Iterating through chr ${chrom}"

# Inputs Clinvar variants by chromosome TSV file
input_file="output/split_chrom/${1}_clin_${chrom}.tsv"
# Outputs TSV file for each chromosome
output_file="output/database_query/${chrom}_database_search.tsv"
info_output_file="OUT/info_output/${chrom}_info_output_exome.txt"

>"$info_output_file"
>"$output_file"

vcf_file_exome_curr=${basepath_exome}/*chr${chrom}.vcf.bgz

#Check if vcf files are indexed
if [ ! -e ${vcf_file_exome_curr}.tbi ]; then
        echo "Index file not found for chromosome ${chrom}. Running bcftools index."
        tabix -p ${vcf_file_exome_curr}
else
        echo "Exome Index file already exists for chromosome ${chrom}. Skipping."
fi


echo "Querying $(ls $vcf_file_exome_curr)"

# Reads through Clinvar variants
while IFS= read -r line 
do
  # Storing Clinvar columns (chromosome, position, reference sequence, alternate sequence, and gene) in bash variables
  chr=chr$(echo "$line" | cut -f 1)
  pos=$(echo "$line" | cut -f 2)
  ref=$(echo "$line" | cut -f 3)
  alt=$(echo "$line" | cut -f 4)
  gene=$(echo "$line" | cut -f 5)

# Extracts gnomAD annotations if position, reference, and alternate are the same

##
#tabix -p ${basepath}/*chr${chrom}.vcf.bgz
##

info=$(tabix ${vcf_file_exome_curr} "${chr}:${pos}-${pos}" | awk -F '\t' -v pos=$pos -v ref=$ref -v alt=$alt '{OFS = "\t"}{if($2 == pos && $4 == ref && $5 == alt)print $4,$5,$8}' | head -n 1)
echo "$info" >> "$info_output_file"

# If data was found, it stores the data, otherwise NA is inputted 
if [ ! -z "$info" ]
  then
      	echo -e "${chr}\t${pos}\t${gene}\t${info}" >> "$output_file"
  else
      	echo -e "${chr}\t${pos}\t${gene}\tNA" >> "$output_file"
  fi

# Extracts chromosome, position, reference, alternate, and gene from Clinvar data
done < <(tail -n+2 "$input_file" | cut -f 2,3,5,6,20)
done
