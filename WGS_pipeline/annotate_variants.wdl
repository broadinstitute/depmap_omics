version 1.0
# Filters out low quality variants, annotate if vairants overlap with segmental duplication or repeatmasker,
# add SnpEff, VEP and open cravat annotation

import "remove_filtered.wdl" as removeFiltered
import "../sandbox/hgvs/hgvs.wdl" as hgvs
import "opencravat_dm.wdl" as openCravat
import "mask_variants.wdl" as mask_variants

workflow annotateVariants {

    input {
        File input_vcf
        String sample_id
        String bcftools_exclude_string = 'FILTER~"weak_evidence"||FILTER~"map_qual"||FILTER~"strand_bias"||FILTER~"slippage"||FILTER~"clustered_events"||FILTER~"base_qual"'
        String vep_pick_order = "mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length,ensembl,refseq"
        Int hgvs_boot_disk_size=100
        Int hgvs_disk_space=200
        String hgvs_vep_data="gs://cds-vep-data/homo_sapiens_vep_110_GRCh38.tar.gz"
        String hgvs_docker_image="us.gcr.io/cds-docker-containers/hgvs"
        Int oc_boot_disk_size=600
        Int oc_disk_space=600
        Int oc_mem=64
    }

    call removeFiltered.RemoveFiltered as RemoveFiltered {
        input:
            sample_id=sample_id,
            input_vcf=input_vcf,
            bcftools_exclude_string=bcftools_exclude_string
    }

    call mask_variants.run_mask_variants as mask_variants {
        input:
            vcf=RemoveFiltered.output_vcf,
            sample_id=sample_id,
    }

    call hgvs.HgvsWorkflow as HgvsWorkflow{
        input:
            input_vcf=mask_variants.mask_annotated,
            sample_id=sample_id,
            boot_disk_size=hgvs_boot_disk_size,
            disk_space=hgvs_disk_space,
            vep_pick_order=vep_pick_order
    }

    call openCravat.opencravat as open_cravat {
        input:
            vcf=HgvsWorkflow.vep,
            boot_disk_size=oc_boot_disk_size,
            disk_space=oc_disk_space,
            memory=oc_mem
    }

    output {
        File hgvs_maf = HgvsWorkflow.maf
        File hgvs_oc_vcf = open_cravat.oc_main_file
        File oc_error_file=open_cravat.oc_error_file
        File oc_log_file=open_cravat.oc_log_file
    }
}
