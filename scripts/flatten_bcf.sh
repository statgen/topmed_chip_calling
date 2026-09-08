#!/bin/bash
set -euo pipefail

f=$1

bcftools view -h $f | bcftools annotate -h <(printf '##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample ID">\n##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele fraction">\n') | awk '{if ($0 ~ /^##/) { print } else { print  $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8 }}'
bcftools query $f -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO[\tSAMPLE=%SAMPLE;VAF=%VAF]\n' -i 'VAF>0' | awk '{ for (i=9; i <= NF; i++) { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8";"$i; }}'
