# AR-genes-gnomadv4_database_query-2023
Pipeline to query ClinVar vcf file based list of autosomal recessive genes to identify pathogenic/likely pathogenic variants. Then query those variants against gnomad v4 database to extract population frequencies for given variants. Then calculate VCF and GCF values based on the input varints.  

Gnomad v4 incorporates both genome and exome data, and since the genome data is very similar to that of v3, we are focusing on only the exome data for this version of the pipeline. There is also a slight change in included populations compared to the v3 version of the pipeline.
