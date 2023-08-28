version 1.0

import "remove_filtered.wdl" as removeFiltered
import "../sandbox/hgvs/hgvs.wdl" as hgvs
import "opencravat_dm.wdl" as openCravat

workflow annotateVariants {

    input {
        File sample_id
        File input_vcf
        String sample_id
        String bcftools_exclude_string = 'FILTER~"weak_evidence"||FILTER~"map_qual"||FILTER~"strand_bias"||FILTER~"slippage"||FILTER~"clustered_events"||FILTER~"base_qual"'
    }

    call removeFiltered.RemoveFiltered as RemoveFiltered {
        input:
            sample_id=sample_name,
            input_vcf=input_vcf,
            bcftools_exclude_string=bcftools_exclude_string
    }

    call hgvs.HgvsWorkflow as HgvsWorkflow{
        input:
            input_vcf=RemoveFiltered.output_vcf,
            sample_id=sample_id,
    }

    call openCravat.opencravat as open_cravat {
        input:
            vcf=HgvsWorkflow.vep,
    }

    output {
        File hgvs_maf = HgvsWorkflow.output_maf
        File hgvs_oc_vcf = open_cravat.oc_main_file
    }
}
