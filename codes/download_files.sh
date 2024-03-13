wget --timestamping 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20230527.vcf.gz' ./
gunzip clinvar_20230527.vcf.gz

for chrom in {1..22}
do
   wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr${chrom}.vcf.bgz"
   wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr${chrom}.vcf.bgz.tbi"
done
