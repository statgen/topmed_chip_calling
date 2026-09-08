#!/bin/bash

f=$1

(q='%CHROM\t%POS\t%REF\t%ALT\t%GERM_P\t%AGE_P\t%NS\t%Func.refGene\t%Gene.refGene\t%ExonicFunc.refGene\t%GeneDetail.refGene\t%AAChange.refGene\n'; echo -ne "$q" | tr -d '%'; bcftools query -f "$q" $f | perl -ne '@a = split; my $prfx = join "\t", @a[0..9]; if (@a[10] eq "." && @a[11] eq "." ) { print "$prfx\t.\t.\n"; } else { if (@a[10] ne ".") { foreach my $gd (split("\\\\x3b", @a[10])) { print "$prfx\t$gd\t.\n"; }} if (@a[11] ne "." ) { foreach my $aa (split(",", @a[11])) { print "$prfx\t.\t$aa\n"; }}}')

#' if ($gdetail ne ".") { foreach my $gd (split("\\\\x3b", $gdetail)) { print "$prfx\t$gd\t.\n"; }} if ($aachange ne "." ) { foreach my $aa (split(",", $aachange)) { print "$prfx\t.\t$aa\n"; }}')

#(q='%CHROM\t%POS\t%REF\t%ALT\t%GERM_P\t%AGE_P\t%NS\t%FILTCNT_artifact_in_normal\t%FILTCNT_base_quality\t%FILTCNT_clustered_events\t%FILTCNT_contamination\t%FILTCNT_duplicate_evidence\t%FILTCNT_fragment_length\t%FILTCNT_germline_risk\t%FILTCNT_mapping_quality\t%FILTCNT_multiallelic\t%FILTCNT_panel_of_normals\t%FILTCNT_read_position\t%FILTCNT_str_contraction\t%FILTCNT_strand_artifact\t%FILTCNT_t_lod\t%Func.refGene\t%Gene.refGene\t%ExonicFunc.refGene\t%GeneDetail.refGene\t%AAChange.refGene\n'; echo -ne "$q" | tr -d '%'; bcftools query -f "$q" $f | perl -ne 'my ($chrom, $pos, $ref, $alt, $germp, $agep, $gdetail, $aachange) = split; my $prfx = "$chrom\t$pos\t$ref\t$alt\t$germp\t$agep"; if ($gdetail ne ".") { foreach my $gd (split("\\\\x3b", $gdetail)) { print "$prfx\t$gd\t.\n"; }} if ($aachange ne "." ) { foreach my $aa (split(",", $aachange)) { print "$prfx\t.\t$aa\n"; }}') 


# bcftools annotate -H '##INFO=<ID=KNOWN_CHIP,Number=0,Type=Flag,Description="Variant is a known missense or LOF CHIP mutation">' -H '##INFO=<ID=LOF,Number=0,Type=Flag,Description="Variant is nonsense, nonstop, or frameshif and in loss of function gene list">' -H '##INFO=<ID=SPLICE,Number=0,Type=Flag,Description="Is a predicted splicing variant and in splicing gene list">' -a chip_flags.tsv.gz -c 'CHROM,POS,REF,ALT,KNOWN_CHIP,LOF,SPLICE' merged_ad_vaf.germ_age_info.min_alt_3.min_vaf_2e-2.sites.filt_cnt_fix.filt_reset.vcf.gz
