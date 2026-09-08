BEGIN { OFS="\t"; FS="\t" } 
{
  if ($0 ~ /^##/) { 
    print
  } else if ($0 ~ /^#CHROM/) {
    print "##INFO=<ID=FILTCNT_artifact_in_normal,Number=1,Type=Integer,Description=\"artifact_in_normal filter count\">"
    print "##INFO=<ID=FILTCNT_base_quality,Number=1,Type=Integer,Description=\"base_quality filter count\">"
    print "##INFO=<ID=FILTCNT_clustered_events,Number=1,Type=Integer,Description=\"clustered_events filter count\">"
    print "##INFO=<ID=FILTCNT_contamination,Number=1,Type=Integer,Description=\"contamination filter count\">"
    print "##INFO=<ID=FILTCNT_duplicate_evidence,Number=1,Type=Integer,Description=\"duplicate_evidence filter count\">"
    print "##INFO=<ID=FILTCNT_fragment_length,Number=1,Type=Integer,Description=\"fragment_length filter count\">"
    print "##INFO=<ID=FILTCNT_germline_risk,Number=1,Type=Integer,Description=\"germline_risk filter count\">"
    print "##INFO=<ID=FILTCNT_mapping_quality,Number=1,Type=Integer,Description=\"mapping_quality filter count\">"
    print "##INFO=<ID=FILTCNT_multiallelic,Number=1,Type=Integer,Description=\"multiallelic filter count\">"
    print "##INFO=<ID=FILTCNT_panel_of_normals,Number=1,Type=Integer,Description=\"panel_of_normals filter count\">"
    print "##INFO=<ID=FILTCNT_read_position,Number=1,Type=Integer,Description=\"read_position filter count\">"
    print "##INFO=<ID=FILTCNT_str_contraction,Number=1,Type=Integer,Description=\"str_contraction filter count\">"
    print "##INFO=<ID=FILTCNT_strand_artifact,Number=1,Type=Integer,Description=\"strand_artifact filter count\">"
    print "##INFO=<ID=FILTCNT_t_lod,Number=1,Type=Integer,Description=\"t_lod filter count\">"
    print $0
  } else if ($7 == "PASS") {
    print
  } else {
    n=split($7, a, ";")
    for (i=1; i<=n; i++) { $8=$8";FILTCNT_"a[i]"=1" }
    print 
  } 
}
