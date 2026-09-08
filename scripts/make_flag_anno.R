library(readr)
library(stringr)

args = commandArgs(trailingOnly=T)

df = as.data.frame(read_tsv(args[1])) #"merged_ad_vaf.germ_age_info.min_alt_3.min_vaf_2e-2.filter_vars.tsv.gz"))

df$var_id = paste(df$CHROM, df$POS, df$REF, df$ALT, sep=":")
df_out = df[!duplicated(df$var_id), c("CHROM","POS","REF","ALT")]
rownames(df_out) = paste(df_out$CHROM, df_out$POS, df_out$REF, df_out$ALT, sep=":")

known_chip_variants = as.data.frame(read_tsv("whitelist_filter_files/NEJM_2017_genes_01262020.MLL_fix.flatten.cleaned.additional_missense_v3.tsv"))
known_chip_keys = paste(known_chip_variants$Accession, known_chip_variants$AAChange)
lof_genes = as.data.frame(read_tsv("whitelist_filter_files/NEJM_2017_genes_01262020_nocr_mll_fix.lof.tsv"))
splice_genes = as.data.frame(read_tsv("whitelist_filter_files/NEJM_2017_genes_01262020_nocr_mll_fix.splice.tsv"))

is_known = paste(str_split_i(df$AAChange.refGene, ":", 2), str_split_i(df$AAChange.refGene, ":", 5)) %in% known_chip_keys
is_lof = str_split_i(df$AAChange.refGene, ":", 2) %in% lof_genes$Accession & grepl("fs|X|\\*", str_split_i(df$AAChange.refGene, ":", 5)) & df$ExonicFunc.refGene %in% c("frameshift_deletion","frameshift_insertion","stopgain","stoploss","startloss")
is_splice = str_split_i(df$GeneDetail.refGene, ":", 1) %in% splice_genes$Accession & grepl("splicing", df$Func.refGene)

df_out$KNOWN_CHIP = as.integer(rownames(df_out) %in% df[is_known,]$var_id)
df_out$LOF = as.integer(rownames(df_out) %in% df[is_lof,]$var_id)
df_out$SPLICE = as.integer(rownames(df_out) %in% df[is_splice,]$var_id)

# exceptions

# ASXL1 NM_015338 Frameshift/nonsense/splice-site in exon 11-12
is_asxl1_exception_lof = str_split_i(df$AAChange.refGene, ":", 2) %in% c("NM_015338") & str_split_i(df$AAChange.refGene, ":", 3) %in% c("exon11","exon12") & grepl("fs|X|\\*", str_split_i(df$AAChange.refGene, ":", 5)) & df$ExonicFunc.refGene %in% c("frameshift_deletion","frameshift_insertion","stopgain","stoploss","startloss")
is_asxl1_exception_splice = str_split_i(df$GeneDetail.refGene, ":", 1) %in% c("NM_015338") & str_split_i(df$GeneDetail.refGene, ":", 2) %in% c("exon11","exon12") & grepl("splicing", df$Func.refGene)

# ASXL2 NM_018263 Frameshift/nonsense/splice-site in exon 11-12
is_asxl2_exception_lof = str_split_i(df$AAChange.refGene, ":", 2) %in% c("NM_018263") & str_split_i(df$AAChange.refGene, ":", 3) %in% c("exon11","exon12") & grepl("fs|X|\\*", str_split_i(df$AAChange.refGene, ":", 5)) & df$ExonicFunc.refGene %in% c("frameshift_deletion","frameshift_insertion","stopgain","stoploss","startloss")
is_asxl2_exception_splice = str_split_i(df$GeneDetail.refGene, ":", 1) %in% c("NM_018263") & str_split_i(df$GeneDetail.refGene, ":", 2) %in% c("exon11","exon12") & grepl("splicing", df$Func.refGene)

# PPM1D NM_003620 Frameshift/nonsense in exon 5 or 6
is_ppm1d_exception = str_split_i(df$AAChange.refGene, ":", 2) %in% c("NM_003620") & str_split_i(df$AAChange.refGene, ":", 3) %in% c("exon5","exon6") & grepl("fs|X|\\*", str_split_i(df$AAChange.refGene, ":", 5)) & df$ExonicFunc.refGene %in% c("frameshift_deletion","frameshift_insertion","stopgain","stoploss","startloss")

# TET2 NM_001127208 missense mutations in catalytic domains (p.1104-1481 and 1843-2002)
aapos = as.integer(substr(str_split_i(df$AAChange.refGene, ":", 5), 4, 7))
is_tet2_exception = str_split_i(df$AAChange.refGene, ":", 2) %in% c("NM_001127208") & df$ExonicFunc.refGene %in% c("nonsynonymous_SNV") & nchar(str_split_i(df$AAChange.refGene, ":", 5)) == 8 & ((aapos >= 1104 & aapos <= 1481) | (aapos >= 1843 & aapos <= 2002))

# CBL  CBL RING finger missense p.381-421 NM_005188
aapos = as.integer(substr(str_split_i(df$AAChange.refGene, ":", 5), 4, 6))
is_cbl_exception = str_split_i(df$AAChange.refGene, ":", 2) %in% c("NM_005188") & df$ExonicFunc.refGene %in% c("nonsynonymous_SNV") & nchar(str_split_i(df$AAChange.refGene, ":", 5)) == 7 & (aapos >= 381 & aapos <= 421)

# CBLB RING finger missense p.372-412 NM_170662
# aapos = as.integer(substr(str_split_i(df$AAChange.refGene, ":", 5), 4, 6))
is_cblb_exception = str_split_i(df$AAChange.refGene, ":", 2) %in% c("NM_170662") & df$ExonicFunc.refGene %in% c("nonsynonymous_SNV") & nchar(str_split_i(df$AAChange.refGene, ":", 5)) == 7 & (aapos >= 372 & aapos <= 412)

df_out$EXCEPTION = as.integer(rownames(df_out) %in% df[is_asxl1_exception_lof | is_asxl1_exception_splice | is_asxl2_exception_lof | is_asxl2_exception_splice | is_ppm1d_exception | is_tet2_exception | is_cbl_exception | is_cblb_exception,]$var_id)

nms = colnames(df_out)
nms[1] = "#CHROM"
colnames(df_out) = nms
write_tsv(df_out, file=stdout())
