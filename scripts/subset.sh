#!/bin/bash
set -euo pipefail

#merged_with_u2af1.sorted.germage_u2af1_thresh.updated_info.filtered.pass_only.carriers_only.vaf.bcf
#merged_with_u2af1.sorted.germage_u2af1_thresh.updated_info.filtered.sites.bcf
#merged_with_u2af1.sorted.germage_u2af1_thresh.updated_info.filtered.pass_only.vaf.bcf
#merged_with_u2af1.sorted.germage_u2af1_thresh.updated_info.filtered.ad_vaf.bcf


#bcftools view ../merged_with_u2af1.sorted.germage_u2af1_thresh.updated_info.filtered.ad_vaf.bcf -S release_ids.txt -Ou | bcftools annotate -x 'INFO/NCA,INFO/VAF_M,INFO/AGE_M_DIFF' -Ou | bcftools +fill-tags -Ou -- -t 'NC:1=int(COUNT(VAF>0))' | bcftools view -i 'NC >= 1' -Ob -o topmed_chip.freeze3.pass_and_fail.ad_vaf.bcf

#bcftools annotate topmed_chip.freeze3.pass_and_fail.ad_vaf.bcf -x 'FMT/AD' -Ou | bcftools view -i 'FILTER="PASS"' -Ob -o topmed_chip.freeze3.pass_only.vaf.bcf

#bcftools view -G topmed_chip.freeze3.pass_and_fail.ad_vaf.bcf -Ob -o topmed_chip.freeze3.pass_and_fail.sites.bcf

bash ../flatten_bcf.sh topmed_chip.freeze3.pass_only.vaf.bcf | bcftools view -Ob -o topmed_chip.freeze3.pass_only.carriers_only.vaf.bcf
