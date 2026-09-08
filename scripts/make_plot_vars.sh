
f=$1

q='%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%DP\t%ECNT\t%IN_PON\t%N_ART_LOD\t%POP_AF\t%P_CONTAM\t%P_GERMLINE\t%RU\t%STR\t%TLOD\t%ANNOVAR_DATE\t%Func.refGene\t%Gene.refGene\t%GeneDetail.refGene\t%ExonicFunc.refGene\t%AAChange.refGene\t%cosmic70\t%ALLELE_END\t%FILTCNT_artifact_in_normal\t%FILTCNT_base_quality\t%FILTCNT_clustered_events\t%FILTCNT_contamination\t%FILTCNT_duplicate_evidence\t%FILTCNT_fragment_length\t%FILTCNT_germline_risk\t%FILTCNT_mapping_quality\t%FILTCNT_multiallelic\t%FILTCNT_panel_of_normals\t%FILTCNT_read_position\t%FILTCNT_str_contraction\t%FILTCNT_strand_artifact\t%FILTCNT_t_lod\t%NS\t%NC\t%NCA\t%AGE_P\t%GERM_P\t%AGE_M_DIFF\t%VAF_M\t%KNOWN_CHIP\t%LOF\t%SPLICE\t%EXCEPTION\n'

echo -ne "$q" | tr -d '%'
bcftools query -f "$q" $f
