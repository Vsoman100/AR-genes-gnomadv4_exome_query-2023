# AR-genes-gnomadv4_exome_query-2023

## **Overview**

This repository can be downloaded with `git clone https://github.com/Vsoman100/AR-genes-gnomadv4_exome_query-2023.git`. 

This pipeline takes a filtered list of ClinVar variants as input (list created by Dr. Mahmoud Aarabi and Chandran Lab). The list is filtered from the clinVar vcf file from Oct 7 2023, and contains variants that are pathogenic/likely pathogenic without conflicting interpretation and that are associated with genes linked to autosomal recessive condition. Next it queries vcf files from a variant frequency database of interest (gnomAD in our case), and extracts population variant frequency data for matching variants that are found. Tabix is used to rapidly retrieve sorted data from vcf files. It then uses that variant frequency data to calculate variant carrier frequency (VCF) and gene carrier frequency (GCF). Equations for those can be found below:

<img width="136" alt="Screenshot 2023-11-09 at 4 49 48 PM" src="https://github.com/Vsoman100/AR-genes-database_query-2023/assets/42780677/ed0a8747-f6de-4a3c-a7b6-acba61714845"> \
<img width="193" alt="Screenshot 2023-11-09 at 4 49 52 PM" src="https://github.com/Vsoman100/AR-genes-database_query-2023/assets/42780677/3ce127f7-0667-4337-b572-bac206a3e6cf">

Python version 3.11.9 (conda env provided) and bash version 4.4.19 were used in this pipeline. Computing resources from the University of Pittsburgh's Center for Research Computing were used. All code was run with through the Pittsburgh Supercomputing Center with 4 cores per node. A complete run filled 2.3 TB (including vcf file sizes) of spaces and took about 30 hours to complete (not including vcf file download time).


Pipeline to query ClinVar vcf file based list of autosomal recessive genes to identify pathogenic/likely pathogenic variants. Then query those variants against gnomad v4 database to extract population frequencies for given variants. Then calculate VCF and GCF values based on the input varints.  

Gnomad v4 incorporates both genome and exome data, and since the genome data is very similar to that of v3, we are focusing on only the exome data for this version of the pipeline. There is also a slight change in included populations compared to the v3 version of the pipeline.

The entire pipeline can be run by running the `GCF_pipeline_master.job` file. Variables to modify before run:
  -date: date of pipeline run \n
  -clinvar_vcf: Full path to your clinvar vcf file. Currently set to our vcf version which was downlaoded on 10/7/23
  -pops_list: list of populations to extract frequency info from gnomAD. Current list contains all pops except for "Remaining Individuals" (RMI) population and "Amish" (AMI)
  -vcf_base_path_exome: Full path to directory containing downloaded gnoMAD exome vcf files
