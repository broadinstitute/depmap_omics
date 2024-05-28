version 1.0


workflow VEP_SV_Workflow {
    input {
        File input_vcf
        String sample_id
        File pLi = "gs://cds-vep-data/pLI_values.txt"
        File LoF = "gs://cds-vep-data/LoFtool_scores.txt"

        File vep_data = "gs://cds-vep-data/homo_sapiens_vep_110_GRCh38.tar.gz"
        File fasta = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz"
        File fai = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz.fai"
        File gzi = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz.gzi"
        Int boot_disk_size=60
        Int disk_space=60
        Int cpu = 10
        Int mem = 80
    }

    call annotate_sv_vep {
        input:
            input_vcf=input_vcf,
            sample_id=sample_id,
            fasta=fasta,
            fai=fai,
            gzi=gzi,
            vep_data=vep_data,
            pLi=pLi,
            LoF=LoF,
            alphamis=alphamis,
            alphamis_idx=alphamis_idx,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu = cpu,
            mem = mem,
    }

    output { 
        File vep=annotate_sv_vep.output_vep_vcf
    }
}

# Standard interface to run vcf to maf
task annotate_sv_vep {
    input {
        File input_vcf
        File fasta
        File fai
        File gzi
        File pLi
        File LoF
        File vep_data
        String sample_id
        String docker_image="us.gcr.io/cds-docker-containers/hgvs"
        String assembly="GRCh38"
        Int preemptible=2
        Int boot_disk_size=60
        Int disk_space=60
        Int cpu = 10
        Int mem = 80
    }

    command {

        bcftools norm -m- -w 10000 -f ~{fasta} -O z -o ~{sample_id}.norm.vcf.gz ~{input_vcf}

        mkdir -p /tmp/Plugins
        wget -c https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/111/AlphaMissense.pm -O /tmp/Plugins/AlphaMissense.pm

        tar -C /tmp -xvzf ~{vep_data} 
        ls /tmp
        chmod 777 /tmp/homo_sapiens
        ls /tmp/homo_sapiens
        cp ~{fai} /tmp
        cp ~{fasta} /tmp
        cp ~{gzi} /tmp
        du -sh /tmp/Homo_sapiens_assembly38.fasta.gz*
       
        cp ~{pLi} ~{LoF} /tmp

        vep --species homo_sapiens --cache --assembly ~{assembly} --no_progress --no_stats --everything --dir /tmp \
            --input_file ~{sample_id}.norm.vcf.gz \
            --output_file ~{sample_id}.norm.snpeff.clinvar.vep.vcf \
            --plugin pLI,/tmp/pLI_values.txt --plugin LoFtool,/tmp/LoFtool_scores.txt \
            --force_overwrite --offline --fasta /tmp/Homo_sapiens_assembly38.fasta.gz --fork ~{cpu} --vcf \
            --pick --max_sv_size 50000000

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
        File output_vep_vcf = "~{sample_id}.norm.snpeff.clinvar.vep.vcf"
    }
}

