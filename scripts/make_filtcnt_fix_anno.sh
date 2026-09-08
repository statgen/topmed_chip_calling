#!/bin/bash

f=$1

(q='%CHROM\t%POS\t%REF\t%ALT\t%FILTCNT_artifact_in_normal\t%FILTCNT_base_quality\t%FILTCNT_clustered_events\t%FILTCNT_contamination\t%FILTCNT_duplicate_evidence\t%FILTCNT_fragment_length\t%FILTCNT_germline_risk\t%FILTCNT_mapping_quality\t%FILTCNT_multiallelic\t%FILTCNT_panel_of_normals\t%FILTCNT_read_position\t%FILTCNT_str_contraction\t%FILTCNT_strand_artifact\t%FILTCNT_t_lod\n'; echo -ne "#$q" | tr -d '%'; bcftools query -f "$q" $f | awk -F'\t' -v 'OFS=\t' '{ for (i=5; i<=NF; i++) { if ($i == ".") { $i="0" }} print}')

# bcftools annotate -a filtcnt_fix_anno.tsv.gz -c 'CHROM,POS,REF,ALT,NS,FILTCNT_artifact_in_normal,FILTCNT_base_quality,FILTCNT_clustered_events,FILTCNT_contamination,FILTCNT_duplicate_evidence,FILTCNT_fragment_length,FILTCNT_germline_risk,FILTCNT_mapping_quality,FILTCNT_multiallelic,FILTCNT_panel_of_normals,FILTCNT_read_position,FILTCNT_str_contraction,FILTCNT_strand_artifact,FILTCNT_t_lod' merged_ad_vaf.germ_age_info.min_alt_3.min_vaf_2e-2.sites.vcf.gz | awk -F'\t' -v 'OFS=\t' '{ if ($0 !~ /^#/) { $7="." } print}' | bgzip > merged_ad_vaf.germ_age_info.min_alt_3.min_vaf_2e-2.sites.filt_cnt_fix.filt_reset.vcf.gz
