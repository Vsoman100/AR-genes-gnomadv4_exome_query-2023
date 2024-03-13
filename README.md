# AR-genes-gnomadv4_exome_query-2023
Pipeline to query ClinVar vcf file based list of autosomal recessive genes to identify pathogenic/likely pathogenic variants. Then query those variants against gnomad v4 database to extract population frequencies for given variants. Then calculate VCF and GCF values based on the input varints.  

Gnomad v4 incorporates both genome and exome data, and since the genome data is very similar to that of v3, we are focusing on only the exome data for this version of the pipeline. There is also a slight change in included populations compared to the v3 version of the pipeline.

The entire pipeline can be run by running the `GCF_pipeline_master.job` file. Variables to modify before run:
  -date: date of pipeline run
  -clinvar_vcf: Full path to your clinvar vcf file. Currently set to our vcf version which was downlaoded on 10/7/23
  -pops_list: list of populations to extract frequency info from gnomAD. Current list contains all pops except for "Remaining Individuals" (RMI) population and "Amish" (AMI)
  -vcf_base_path_exome: Full path to directory containing downloaded gnoMAD exome vcf files
