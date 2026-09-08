#!/bin/bash

vcf=$1
filt_count_anno=$2
chip_flags_anno=$3

bcftools annotate $vcf -a $filt_count_anno -c 'CHROM,POS,REF,ALT,FILTCNT_artifact_in_normal,FILTCNT_base_quality,FILTCNT_clustered_events,FILTCNT_contamination,FILTCNT_duplicate_evidence,FILTCNT_fragment_length,FILTCNT_germline_risk,FILTCNT_mapping_quality,FILTCNT_multiallelic,FILTCNT_panel_of_normals,FILTCNT_read_position,FILTCNT_str_contraction,FILTCNT_strand_artifact,FILTCNT_t_lod' | bcftools annotate -H '##INFO=<ID=KNOWN_CHIP,Number=0,Type=Flag,Description="Variant is a known missense or LOF CHIP mutation">' -H '##INFO=<ID=LOF,Number=0,Type=Flag,Description="Variant is nonsense, nonstop, or frameshif and in loss of function gene list">' -H '##INFO=<ID=SPLICE,Number=0,Type=Flag,Description="Is a predicted splicing variant and in splicing gene list">' -H '##INFO=<ID=EXCEPTION,Number=0,Type=Flag,Description="Variant meets special gene-specific criteria">' -a $chip_flags_anno -c 'CHROM,POS,REF,ALT,KNOWN_CHIP,LOF,SPLICE,EXCEPTION' -Ov

