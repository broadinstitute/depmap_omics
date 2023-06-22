version 1.0


workflow run_echtvar_cosmic {
    input {
        String sample_id
        File input_vcf
    }

    call echtvar_cosmic {
        input:
            input_vcf=input_vcf,
            sample_id=sample_id
    }

    output {
        File annotated_vcf_gnomad=echtvar_cosmic.annotated_vcf_gnomad
        File annotated_vcf_gnomad_cosmic=echtvar_cosmic.annotated_vcf_gnomad_cosmic
    }
}

task echtvar_cosmic {
    input {
        String sample_id
        File input_vcf
        File ref_fasta="gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
        File gnomad="gs://echtvar-dabases/gnomad.v3.1.2.echtvar.v2.zip"
        File cosmic="gs://echtvar-dabases/cosmic_coding_ev.zip"

        String docker_image="depmapomics:test"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=60
        Int cpu = 4
        Int mem = 32
    }

    command {
        bcftools norm -m- ${input_vcf} -w 10000 -f ${ref_fasta} -O b -o ${sample_id}_bcftools_normalized.bcf

        ~/bin/echtvar anno -e ${gnomad} ${sample_id}_bcftools_normalized.bcf ${sample_id}_echtvar_gnomad.bcf

        ~/bin/echtvar anno -e ${cosmic} ${sample_id}_echtvar_gnomad.bcf ${sample_id}_echtvar_gnomad_cosmic.bcf
    }

    runtime {
        disks: "local-disk ~{disk_space} HDD"
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {
        File annotated_vcf_gnomad = "${sample_id}_echtvar_gnomad.bcf"
        File annotated_vcf_gnomad_cosmic = "${sample_id}_echtvar_gnomad_cosmic.bcf"
    }
}
