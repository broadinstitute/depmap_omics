version 1.0


workflow run_vcf2maf {
    input {
        File input_vcf
        String sample_id
        File vep_data = "gs://cds-vep-data/homo_sapiens_vep_109_GRCh38.tar.gz"
        File fasta = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz"
        File fai = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz.fai"
        File gzi = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz.gzi"
    }

    call vcf2maf {
        input:
            input_vcf=input_vcf,
            sample_id=sample_id,
            fasta=fasta,
            vep_data=vep_data
    }

    output {
        #File vcf=vcf2maf.output_vcf
        File maf=vcf2maf.output_maf
    }
}

# Standard interface to run vcf to maf
task vcf2maf {
    input {
        File input_vcf
        File fasta
        File vep_data
        String sample_id

        String docker_image="us.gcr.io/cds-docker-containers/vcf2maf:test"
        String assembly="GRCh38"
        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=20
        Int cpu = 6
        Int mem = 32
    }

    command {
        tar -C /tmp -xvzf ~{vep_data} 
        ls /tmp
        chmod 777 /tmp/homo_sapiens
        ls /tmp/homo_sapiens

        perl /tmp/vcf2maf/vcf2maf.pl \
         --input-vcf ~{input_vcf} \
         --output-maf ~{sample_id}.maf \
         --ref ~{fasta} \
         --vep-path /opt/conda/envs/vep/bin/ \
         --vep-data /tmp \
         --ncbi-build ~{assembly}
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
        #File output_vcf = "~{sample_id}.vep.vcf"        
        File output_maf = "~{sample_id}.maf"        
    }
}