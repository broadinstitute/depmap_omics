version 1.0

import "../sandbox/hgvs/hgvs.wdl" as hgvs
import "opencravat_dm.wdl" as openCravat

workflow annotateVariants {

    input {
        File input_vcf
        String sample_id
    }

    call hgvs.HgvsWorkflow as HgvsWorkflow{
        input:
            input_vcf=input_vcf,
            sample_id=sample_id,
    }

    call openCravat.opencravat as open_cravat {
        input:
            vcf=HgvsWorkflow.vep,
    }

    output {
        hgvs_oc_vcf = open_cravat.oc_main_file
    }
