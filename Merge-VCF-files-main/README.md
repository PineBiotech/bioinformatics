# Merge-VCF-files
Find intersection of two VCF files

Usage:
python merge_vcf_dict_RLN.py -i1 diploid_samples.vcf -i2 clinVar_pathogenic_tp53_edt.vcf -o common_dip_clvar_tp53.vcf

The python script merges two VCF files based on the chromosomal number, position and the reference nucleotide. 
While merging, the output file is populated with the variants, qscore from the first vcf file.
