library(readr)

args <- commandArgs(trailingOnly = TRUE)

df = as.data.frame(read_tsv(args[1], col_names=c("#CHROM","POS","REF","ALT","MUTATION","DP","ALT_DP")))

stopifnot(df$MUTATION == c("D14G","S34Y","S34F","R35Q","R35L","R156H","Q157P","Q157R","D14G","S34Y","S34F","R35Q","R35L","R156H","Q157P","Q157R"))

df_out = df[1:8,1:2]
df_out$ID = "."
df_out$REF = df$REF[1:8]
df_out$ALT = df$ALT[1:8]
df_out$QUAL = "."
df_out$FILTER = "PASS"
df_out$INFO = "."
df_out$FORMAT = "GT:AD"

alt = df$ALT_DP[1:8] + df$ALT_DP[9:16]
dp = df$DP[1:8] + df$DP[9:16]

df_out[,args[2]] = paste0(ifelse(alt > 0, "0/1", "0/0"), ":", dp - alt, ",", alt)

write_tsv(df_out[alt > 0,], file=stdout())
