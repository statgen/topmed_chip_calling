
bcftools annotate merged/merged_ad_vaf.bcf -c 'FILTER,INFO' -a filter/merged.updated_info.filtered.sites.vcf.gz -Ob -o filter/merged.updated_info.filtered.ad_vaf.bcf

#bcftools annotate merged_with_u2af1.sorted.germage_u2af1_thresh.bcf -c 'FILTER,INFO' -a merged_with_u2af1.sorted.germage_u2af1_thresh.updated_info.filtered.sites.vcf.gz -Ob -o merged_with_u2af1.sorted.germage_u2af1_thresh.updated_info.filtered.ad_vaf.bcf
