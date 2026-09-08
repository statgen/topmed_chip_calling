
{ 
  if ($0 ~ /INFO=<ID=DP,/) { 
    print "##INFO=<ID=DP,Number=1,Type=Float,Description=\"Mean approximate read depth; some reads may have been filtered\">" 
  } else if ($0 ~ /INFO=<ID=ECNT,/) { 
    print "##INFO=<ID=ECNT,Number=1,Type=Float,Description=\"Mean number of events in this haplotype\">" 
  } else if ($0 ~ /INFO=<ID=NLOD,/) {
    print "##INFO=<ID=NLOD,Number=A,Type=Float,Description=\"Mean normal LOD score\">"  
  } else if ($0 ~ /INFO=<ID=POP_AF,/) {
    print "##INFO=<ID=POP_AF,Number=A,Type=Float,Description=\"Mean population allele frequencies of alt alleles\">"
  } else if ($0 ~ /INFO=<ID=N_ART_NLOD,/) {
    print "##INFO=<ID=N_ART_LOD,Number=A,Type=Float,Description=\"Mean log odds of artifact in normal with same allele fraction as tumor\">"
  } else if ($0 ~ /INFO=<ID=P_CONTAM,/) {
    print "##INFO=<ID=P_CONTAM,Number=1,Type=Float,Description=\"Mean posterior probability for alt allele to be due to contamination\">"
  } else if ($0 ~ /INFO=<ID=P_GERMLINE,/) {
    print "##INFO=<ID=P_GERMLINE,Number=A,Type=Float,Description=\"Mean posterior probability for alt allele to be germline variants\">"
  } else if ($0 ~ /INFO=<ID=TLOD,/) {
    print "##INFO=<ID=TLOD,Number=A,Type=Float,Description=\"Mean tumor LOD score\">"
  } else {
    print
  } 
}

