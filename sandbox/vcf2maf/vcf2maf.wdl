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
            fai=fai,
            gzi=gzi,
            vep_data=vep_data
    }

    output {
        File maf=vcf2maf.output_maf
    }
}

# Standard interface to run vcf to maf
task vcf2maf {
    input {
        File input_vcf
        File fasta
        File fai
        File gzi
        File vep_data
        String sample_id

        String docker_image="us.gcr.io/cds-docker-containers/vcf2maf:test"
        String assembly="GRCh38"
        Int preemptible=2
        Int boot_disk_size=60
        Int disk_space=60
        Int cpu = 4
        Int mem = 80
    }

    command {
        tar -C /tmp -xvzf ~{vep_data} 
        ls /tmp
        chmod 777 /tmp/homo_sapiens
        ls /tmp/homo_sapiens
        cp ~{fai} /tmp
        cp ~{fasta} /tmp
        cp ~{gzi} /tmp
        du -sh /tmp/Homo_sapiens_assembly38.fasta.gz*

        IS_GZ=`echo ~{input_vcf} | grep -Pic ".gz"`
        if [ "$IS_GZ" -eq "1" ]; 
        then
           zcat ~{input_vcf} > ~{sample_id}.vcf
           perl /tmp/vcf2maf/vcf2maf.pl \
            --input-vcf ~{sample_id}.vcf \
            --output-maf ~{sample_id}.maf \
            --ref /tmp/Homo_sapiens_assembly38.fasta.gz \
            --vep-path /opt/conda/envs/vep/bin/ \
            --vep-data /tmp \
            --ncbi-build ~{assembly}
        else
           perl /tmp/vcf2maf/vcf2maf.pl \
            --input-vcf ~{input_vcf} \
            --output-maf ~{sample_id}.maf \
            --ref /tmp/Homo_sapiens_assembly38.fasta.gz \
            --vep-path /opt/conda/envs/vep/bin/ \
            --vep-data /tmp \
            --ncbi-build ~{assembly}
        fi
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
        File output_maf = "~{sample_id}.maf"        
    }
}
