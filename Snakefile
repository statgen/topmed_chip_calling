import os
import math


configfile: "config.yml"
singularity: "chip-gcp.sif"


mem_step_size = 16000

sample_ids = None

def get_ids():
    if sample_ids is None:
        with open(config["ids_file"]) as f:
            sample_ids = f.read().splitlines()
    return sample_ids

rule unfiltered_vcf:
    input:
        config["cram_expr"]
    output:
        "single_sample/m2/{sample_id}.vcf.gz"
    resources:
        mem_mb = 16000
    shell:
        """
        tmp_dir=`mktemp -d`
        tmp_out_vcf=$tmp_dir/$(basename {output})
        tmp_in={input}
        set -eu
        # build the reference sequence cache
        export REF_CACHE={config[ref_cache]} #./md5/%2s/%2s/%s

        #export GATK_LOCAL_JAR="{config[gatk4_jar]}"
        #export PATH="/MitoHPC/bin:$PATH"


        #touch bamout.bam
        #echo "" > normal_name.txt

        gatk --java-options "-Xmx{resources.mem_mb}m" GetSampleName -R {config[ref_fasta]} -I $tmp_in -O $tmp_dir/tumor_name.txt -encode
        tumor_command_line="-I $tmp_in -tumor `cat $tmp_dir/tumor_name.txt`"

        gatk --java-options "-Xmx{resources.mem_mb}m" Mutect2 \
          -R {config[ref_fasta]} \
          $tumor_command_line \
          --germline-resource {config[gnomad]} \
          -pon {config[pon]} \
          -L {config[intervals]} \
          -O $tmp_out_vcf
        rc=$?

        if [[ $rc == 0 ]]; then
          mv ${{tmp_out_vcf}}* $(dirname {output})/
          rc=$?
        fi

        rm -r $tmp_dir
        exit $rc
        """

rule contamination_table:
    input:
        config["cram_expr"]
    output:
        pileups="single_sample/contam/{sample_id}.pileup.tsv",
        contamination_table="single_sample/contam/{sample_id}.contam.tsv",
        maf_segments="single_sample/contam/{sample_id}.segments.tsv"
    resources:
        mem_mb = 16000
    shell:
      """
      set -eu
      tmp_dir=`mktemp -d`
      tmp_pileups_table=$tmp_dir/$(basename {output.pileups})
      tmp_contam_table=$tmp_dir/$(basename {output.contamination_table})
      tmp_segments_table=$tmp_dir/$(basename {output.maf_segments})

      export REF_CACHE={config[ref_cache]}
      #export GATK_LOCAL_JAR="{config[gatk4_jar]}"
      #export PATH="/MitoHPC/bin:$PATH"


      gatk --java-options "-Xmx{resources.mem_mb}m" GetPileupSummaries \
        -R {config[ref_fasta]} \
        -I {input} \
        -L {config[intervals]} \
        -V {config[variants_for_contamination]} \
        -O $tmp_pileups_table

      gatk --java-options "-Xmx{resources.mem_mb}m" CalculateContamination \
        -I $tmp_pileups_table \
        -O $tmp_contam_table \
        --tumor-segmentation $tmp_segments_table
      
      mv $tmp_pileups_table {output.pileups}
      mv $tmp_contam_table {output.contamination_table}
      mv $tmp_segments_table {output.maf_segments}
      rm -r $tmp_dir
      """


rule filtered_vcf:
    input:
        unfiltered_vcf = rules.unfiltered_vcf.output,
        contamination_table = rules.contamination_table.output.contamination_table,
        maf_segments = rules.contamination_table.output.maf_segments
    output:
        "single_sample/filter/{sample_id}.filtered.vcf.gz"
    resources:
        mem_mb = 16000
    shell:
        """
        set -eu
        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})

        export REF_CACHE={config[ref_cache]}
        #export GATK_LOCAL_JAR="{config[gatk4_jar]}"
        #export PATH="/MitoHPC/bin:$PATH"

        gatk --java-options "-Xmx{resources.mem_mb}m" FilterMutectCalls \
          -V {input.unfiltered_vcf} \
          -O $tmp_out \
          --contamination-table {input.contamination_table} \
          --tumor-segmentation {input.maf_segments}

        mv $tmp_out {output}
        mv ${{tmp_out}}.tbi {output}.tbi
        rm -r $tmp_dir
        """


rule annovar_vcf:
    input:
        rules.filtered_vcf.output
    output:
        txt = "single_sample/annovar/{sample_id}.annovar_out.hg38_multianno.txt",
        vcf = "single_sample/annovar/{sample_id}.annovar_out.hg38_multianno.vcf"
    resources:
        mem_mb = mem_step_size
    shell:
        """
        set -eu
        tmp_dir=`mktemp -d`
        out_prefix=$tmp_dir/{wildcards.sample_id}.annovar_out

        perl /annovar/table_annovar.pl {input} /annovar/humandb/ -buildver hg38 -out $out_prefix -remove -protocol refGene,cosmic70 -operation g,f -nastring . -vcfinput

        mv ${{out_prefix}}.hg38_multianno.txt {output.txt}
        mv ${{out_prefix}}.hg38_multianno.vcf {output.vcf}
        rm -r $tmp_dir
        """

rule all_annovar:
    input: [rules.annovar_vcf.output[0].format(sample_id=id) for id in get_ids()]


rule region_counts:
    singularity: "docker:jweinstk/pileup_region"
    input:
        config["cram_expr"]
    output:
        "single_sample/region_counts/{sample_id}.counts.tsv"
    resources:
        mem_mb = mem_step_size
    shell:
        """
        set -eu
        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})

        pileup_region u2af1_vars.txt {input} {config[ref_fasta]} > $tmp_out

        mv $tmp_out {output}
        rm -r $tmp_dir
        """


rule merged_batch:
    input:
        lambda wc: ["single_sample/annovar/{sid}.annovar_out.hg38_multianno.vcf".format(sid=s) for s in get_ids()[(int(wc.batch_beg)-1):int(wc.batch_end)]]
    output:
        "merged_batch/merged_batch.{batch_beg}_{batch_end}.bcf"
    params:
        ref = "resources/ref/hs38DH.fa"
    shell:
        """
        set +e
        set -uo pipefail

        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/out/$(basename {output})
        mkdir $tmp_dir/in/
        mkdir $tmp_dir/out/

        for f in {input}; do
          awk -f update_merge_info_headers.awk $f | awk -f set_filter_count_fields.awk | bcftools annotate -x '^FMT/AD,^FMT/GT,INFO/RPA,INFO/NLOD' -Ou | bcftools norm --check-ref e --fasta-ref {params.ref}  --multiallelics - -Ob -o $tmp_dir/in/$(basename $f .vcf).bcf \
          && bcftools index $tmp_dir/in/$(basename $f .vcf).bcf
          rc=$?

          if [[ $rc != 0 ]]; then
            rm -r $tmp_dir
            exit $rc
          fi
        done

        bcftools merge $tmp_dir/in/*.bcf --filter-logic + -m none --info-rules DP:sum,ECNT:sum,N_ART_LOD:sum,POP_AF:sum,P_CONTAM:sum,P_GERMLINE:sum,TLOD:sum,FILTCNT_artifact_in_normal:sum,FILTCNT_base_quality:sum,FILTCNT_clustered_events:sum,FILTCNT_contamination:sum,FILTCNT_duplicate_evidence:sum,FILTCNT_fragment_length:sum,FILTCNT_germline_risk:sum,FILTCNT_mapping_quality:sum,FILTCNT_multiallelic:sum,FILTCNT_panel_of_normals:sum,FILTCNT_read_position:sum,FILTCNT_str_contraction:sum,FILTCNT_strand_artifact:sum,FILTCNT_t_lod:sum -Ou | bcftools +fill-tags -Ou -- -t NS | bcftools annotate -x 'FMT/GT' -Ob -o $tmp_out && bcftools index $tmp_out
        rc=$?

        if [[ $rc == 0 ]]; then
          mv $tmp_out {output} && mv $tmp_out.csi {output}.csi
          rc=$?
        fi

        rm -r $tmp_dir
        exit $rc
        """

rule merged_ad_bcf:
    input: [rules.merged_batch.output[0].format(batch_beg=i*1000+1, batch_end=min((i+1)*1000, len(get_ids()))) for i in range(0, math.ceil(len(get_ids())/1000))]
    output:
        "merged/merged_ad.bcf"
    shell:
        """
        set -euo pipefail
        
        q='%CHROM\t%POS\t%REF\t%ALT\t%NS\t%DP\t%ECNT\t%N_ART_LOD\t%POP_AF\t%P_CONTAM\t%P_GERMLINE\t%TLOD\n';

        bcftools merge {input} --filter-logic + -m none --info-rules DP:sum,ECNT:sum,N_ART_LOD:sum,POP_AF:sum,P_CONTAM:sum,P_GERMLINE:sum,TLOD:sum,FILTCNT_artifact_in_normal:sum,FILTCNT_base_quality:sum,FILTCNT_clustered_events:sum,FILTCNT_contamination:sum,FILTCNT_duplicate_evidence:sum,FILTCNT_fragment_length:sum,FILTCNT_germline_risk:sum,FILTCNT_mapping_quality:sum,FILTCNT_multiallelic:sum,FILTCNT_panel_of_normals:sum,FILTCNT_read_position:sum,FILTCNT_str_contraction:sum,FILTCNT_strand_artifact:sum,FILTCNT_t_lod:sum,NS:sum -Ob -o {output}.tmp &&
        (echo -ne "#$q" | tr -d '%'; bcftools query -f "$q" {output}.tmp | awk -F'\t' -v 'OFS=\t' '{{ for (i=6; i<=NF; i++) {{ $i = $i / $5 }}; print }}') | bgzip > {output}.mean_anno.tsv.gz &&
        tabix -s1 -b2 -e2 {output}.mean_anno.tsv.gz &&
        bcftools annotate -a {output}.mean_anno.tsv.gz -c 'CHROM,POS,REF,ALT,-,DP,ECNT,N_ART_LOD,POP_AF,P_CONTAM,P_GERMLINE,TLOD' {output}.tmp -Ob -o {output} &&
        bcftools index {output}
        """


rule merged_u2af1_batch:
    input:
        lambda wc: ["single_sample/region_counts/{sid}.counts.tsv".format(sid=s) for s in get_ids()[(int(wc.batch_beg)-1):int(wc.batch_end)]]
    output:
        "merged_u2af1_batch/merged_u2af1_batch.{batch_beg}_{batch_end}.bcf"
    params:
        ref = "resources/ref/hs38DH.fa"
    shell:
        """
        set +e
        set -uo pipefail

        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/out/$(basename {output})
        mkdir $tmp_dir/in/
        mkdir $tmp_dir/out/


        for f in {input}; do
          id=$(basename $f .counts.tsv)
          cols="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$id\n"

          (cat mutect2_headers.txt; echo -ne "$cols"; cat $f | perl make_u2af1_vcf.pl) | bcftools sort | bcftools norm --check-ref e --fasta-ref {params.ref}  --multiallelics - -Ob -o $tmp_dir/in/$id.bcf \
          && bcftools index $tmp_dir/in/$id.bcf
          rc=$?

          if [[ $rc != 0 ]]; then
            rm -r $tmp_dir
            exit $rc
          fi
        done

        bcftools merge $tmp_dir/in/*.bcf --filter-logic + -m none -Ou | bcftools +fill-tags -Ob -o $tmp_out -- -t NS && bcftools index $tmp_out
        # bcftools merge $tmp_dir/in/*.bcf --filter-logic + -m none -Ou | bcftools +fill-tags -Ou -- -t NS | bcftools annotate -x 'FMT/GT' -Ob -o $tmp_out && bcftools index $tmp_out
        rc=$?

        if [[ $rc == 0 ]]; then
          mv $tmp_out {output} && mv $tmp_out.csi {output}.csi
          rc=$?
        fi

        rm -r $tmp_dir
        exit $rc
        """

rule merged_u2af1_vcf:
    input: [rules.merged_u2af1_batch.output[0].format(batch_beg=i*1000+1, batch_end=min((i+1)*1000, len(get_ids()))) for i in range(0, math.ceil(len(get_ids())/1000))]
    output:
        "merged_u2af1/merged_u2af1_ad.vcf.gz"
    shell:
        """
        set -euo pipefail


        bcftools merge {input} --filter-logic + -m none --info-rules NS:sum -Oz -o {output} &&
        bcftools index {output}
        """


"""
singularity exec chip_with_annovar.sif perl /annovar/table_annovar.pl merged_u2af1/merged_u2af1_ad.vcf.gz /annovar/humandb/ -buildver hg38 -out merged_u2af1/merged_u2af1_ad.annovar -remove -protocol refGene,cosmic70 -operation g,f -nastring . -vcfinput; echo sleeping $?; sleep 30d

bcftools reheader -s nhlbi.6945.freeze10.exclude_and_keep.convert.tab merged_u2af1/merged_u2af1_ad.annovar.hg38_multianno.vcf | bcftools view -S merged/merged_ad.nwd_ids.sample_ids.txt -Oz -o merged_u2af1/merged_u2af1_ad.annovar.hg38_multianno.nwd_ids.vcf.gz

bcftools annotate -x 'INFO/AC,INFO/AN,FMT/GT' merged_u2af1/merged_u2af1_ad.annovar.hg38_multianno.nwd_ids.vcf.gz -Oz -o merged_u2af1/merged_u2af1_ad.annovar.hg38_multianno.nwd_ids.no_gt.vcf.gz

bcftools view -t '^chr21:43092956-43107570' merged/merged_ad.mean_anno.nwd_ids.bcf -Ob -o merged/merged_ad.mean_anno.nwd_ids.exclude_u2af1.bc

bcftools concat merged/merged_ad.mean_anno.nwd_ids.exclude_u2af1.bcf merged_u2af1/merged_u2af1_ad.annovar.hg38_multianno.nwd_ids.no_gt.vcf.gz -Ob -o merged_with_u2af1.bcf -a

bcftools sort merged_with_u2af1.bcf -Ob -o merged_with_u2af1.sorted.bcf

"""

rule vaf_calls:
    input:
        "merged/merged_ad.bcf"
    output:
        "merged/merged_ad_vaf.bcf"
    shell:
        """
        set -euo pipefail
        
        chfilter {input} {output}
        bcftools index --force {output}
        """

rule sites_only:
    input:
        rules.vaf_calls.output
    output:
        "merged/merged_sites.bcf"
    shell:
        """
        set -euo pipefail 

        bcftools view -G {input} -Ob -o {output}
        bcftools index --force {output}
        """

rule filt_count_anno:
    input:
        rules.sites_only.output
    output:
        "filter/filt_counts.tsv.gz"
    shell:
        """
        set -euo pipefail
        
        scripts/make_filtcnt_fix_anno.sh {input} | bgzip > {output}
        tabix -s1 -b2 -e2 --force {output} 
        """

rule chip_flags_anno:
    input:
        rules.sites_only.output
    output:
        "filter/chip_flags.tsv.gz"
    shell:
        """
        set -euo pipefail

        scripts/extract_mutations_by_accession.sh {input} | bgzip > filter/mut_by_acc.tsv.gz
        Rscript scripts/make_flag_anno.R filter/mut_by_acc.tsv.gz | bgzip > {output}
        tabix -s1 -b2 -e2 --force {output}
        """

rule updated_info_vcf:
    input:
        vcf = rules.sites_only.output,
        filt_counts = rules.filt_count_anno.output,
        chip_flags = rules.chip_flags_anno.output
    output:
        "filter/merged.updated_info.sites.vcf.gz"
    shell:
        """
        set -euo pipefail

        scripts/annotate_info.sh {input} | bgzip > {output}
        bcftools index --force {output}
        """

rule filtered_sites_vcf:
    input:
        rules.updated_info_vcf.output
    output:
        "merged.updated_info.filtered.sites.vcf.gz"
    shell:
        """
        set -euo pipefail

        bcftools filter {input} -i '1==1' -m x | bcftools filter -m + -s OFF_TARGET -i 'KNOWN_CHIP=1 || LOF=1 || SPLICE=1 || EXCEPTION=1' | bcftools filter -m + -s artifact_in_normal -e '(FILTCNT_artifact_in_normal / NS) >= 0.50' | bcftools filter -m + -s base_quality -e '(FILTCNT_base_quality / NS) >= 0.20' | bcftools filter -m + -s clustered_events -e '(FILTCNT_clustered_events / NS) >= 0.40' | bcftools filter -m + -s contamination -e '(FILTCNT_contamination / NS) >= 0.50' | bcftools filter -m + -s duplicate_evidence -e '(FILTCNT_duplicate_evidence / NS) >= 0.50' | bcftools filter -m + -s fragment_length -e '(FILTCNT_fragment_length / NS) >= 0.50' | bcftools filter -m + -s germline_risk -e '(FILTCNT_germline_risk / NS) >= 0.90' | bcftools filter -m + -s mapping_quality -e '(FILTCNT_mapping_quality / NS) >= 0.50' | bcftools filter -m + -s multiallelic -e '(FILTCNT_multiallelic / NS) >= 0.50' | bcftools filter -m + -s panel_of_normals -e '(FILTCNT_panel_of_normals / NS) >= 0.50' | bcftools filter -m + -s read_position -e '(FILTCNT_read_position / NS) >= 0.50' | bcftools filter -m + -s str_contraction -e '(FILTCNT_str_contraction / NS) >= 0.75' | bcftools filter -m + -s strand_artifact -e '(FILTCNT_strand_artifact / NS) >= 0.50' | bcftools filter -m + -s t_lod -e '(FILTCNT_t_lod / NS) >= 0.85' | bcftools filter -m + -s GERMLINE_PROB -e 'P_GERMLINE != "." && P_GERMLINE > -1.30103 && (GERM_P == "." || GERM_P > -1.30103)' -Oz -o {output}
        bcftools index --force {output}
        """
