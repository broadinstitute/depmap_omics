version 1.0


workflow HgvsWorkflow {
    input {
        File input_vcf
        String sample_id
        File pLi = "gs://cds-vep-data/pLI_values.txt"
        File LoF = "gs://cds-vep-data/LoFtool_scores.txt"

        File clinvar_data = "gs://snpsift_data/clinvar-latest.vcf.gz"
        File clinvar_data_tbi = "gs://snpsift_data/clinvar-latest.vcf.gz.tbi"

        File snpeff = "gs://snpsift_data/snpEff.jar"
        File snpeff_config = "gs://snpsift_data/snpEff.config"
        File snpsift = "gs://snpsift_data/SnpSift.jar"

        File alphamis = "gs://cds-vep-data/AlphaMissense_hg38.tsv.gz"
        File alphamis_idx = "gs://cds-vep-data/AlphaMissense_hg38.tsv.gz.tbi"

        File vep_data = "gs://cds-vep-data/homo_sapiens_vep_110_GRCh38.tar.gz"
        File fasta = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz"
        File fai = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz.fai"
        File gzi = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz.gzi"
        Int boot_disk_size=60
        Int disk_space=60
        Int cpu = 10
        Int mem = 80
    }

    call annotate_hgvs_task {
        input:
            input_vcf=input_vcf,
            sample_id=sample_id,
            fasta=fasta,
            fai=fai,
            gzi=gzi,
            vep_data=vep_data,
            snpeff=snpeff,
            snpsift=snpsift,
            snpeff_config=snpeff_config,
            clinvar_data=clinvar_data,
            clinvar_data_tbi=clinvar_data_tbi,
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
        File maf=annotate_hgvs_task.output_maf   
        File vep=annotate_hgvs_task.output_vep_vcf
    }
}

# Standard interface to run vcf to maf
task annotate_hgvs_task {
    input {
        File input_vcf
        File fasta
        File fai
        File alphamis
        File alphamis_idx
        File gzi
        File pLi
        File LoF
        File vep_data
        File snpeff
        File snpeff_config
        File snpsift
        File clinvar_data
        File clinvar_data_tbi
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

        cp ~{snpeff_config} ~{snpeff} .

        java -Xmx8g -jar ~{snpeff} GRCh38.mane.1.0.ensembl ~{sample_id}.norm.vcf.gz > ~{sample_id}.norm.snpeff.vcf

        mkdir -p ./db/GRCh38/clinvar/
        cp ~{clinvar_data} ~{clinvar_data_tbi} ./db/GRCh38/clinvar/
        java -jar ~{snpsift} Annotate -clinvar -db ~{clinvar_data} ~{sample_id}.norm.snpeff.vcf > ~{sample_id}.norm.snpeff.clinvar.vcf

        mkdir -p /tmp/Plugins
        wget -c https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/111/AlphaMissense.pm -O /tmp/Plugins/AlphaMissense.pm

        tar -C /tmp -xzf ~{vep_data}
        chmod 777 /tmp/homo_sapiens
        cp ~{fai} /tmp
        cp ~{fasta} /tmp
        cp ~{gzi} /tmp

        cp ~{pLi} ~{LoF} ~{alphamis} ~{alphamis_idx} /tmp

        vep --species homo_sapiens --cache --assembly ~{assembly} --no_progress --no_stats --everything --dir /tmp --input_file ~{sample_id}.norm.snpeff.clinvar.vcf \
            --output_file ~{sample_id}.norm.snpeff.clinvar.vep.vcf \
            --plugin pLI,/tmp/pLI_values.txt --plugin LoFtool,/tmp/LoFtool_scores.txt \
            --plugin AlphaMissense,file=/tmp/AlphaMissense_hg38.tsv.gz \
            --force_overwrite --offline --fasta /tmp/Homo_sapiens_assembly38.fasta.gz --fork ~{cpu} --vcf \
            --pick 

        perl /vcf2maf/vcf2maf.pl \
            --input-vcf ~{sample_id}.norm.snpeff.clinvar.vep.vcf \
            --output-maf ~{sample_id}.maf \
            --ref /tmp/Homo_sapiens_assembly38.fasta.gz \
            --vep-path /opt/conda/envs/vep/bin/ \
            --vep-data /tmp/ \
            --ncbi-build ~{assembly} --inhibit-vep
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
        File output_vep_vcf = "~{sample_id}.norm.snpeff.clinvar.vep.vcf"
    }
}

