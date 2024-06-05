version 1.0


workflow VEP_SV_Workflow {
    input {
        File input_vcf
        String sample_id

        File vep_data = "gs://cds-vep-data/homo_sapiens_vep_110_GRCh38.tar.gz"
        File fasta = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz"
        File fai = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz.fai"
        File gzi = "gs://cds-vep-data/Homo_sapiens_assembly38.fasta.gz.gzi"
        File gnomad = "gs://gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz"
        File gnomad_idx = "gs://gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz.tbi"

        Int max_sv_size = 50000000
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
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu = cpu,
            mem = mem,
            gnomad=gnomad,
            gnomad_idx=gnomad_idx,
            max_sv_size=max_sv_size
    }

    output { 
        File vep_annotated_sv = annotate_sv_vep.output_vep_vcf
        File vep_sv_stats = annotate_sv_vep.output_vep_stats
    }
}

# Standard interface to run vcf to maf
task annotate_sv_vep {
    input {
        File input_vcf
        File fasta
        File fai
        File vep_data
        String sample_id
        File gzi
        File gnomad
        File gnomad_idx
        File sv_plugin="gs://cds-vep-data/StructuralVariantOverlap.pm"
        String docker_image="ensemblorg/ensembl-vep:release_112.0"
        String assembly="GRCh38"
        Int preemptible=2
        Int boot_disk_size=60
        Int disk_space=60
        Int cpu = 10
        Int mem = 80
        Int max_sv_size = 50000000
        Int sv_match_percentage=80
    }

    String gnomad_basename = basename(gnomad)
    String plugin_basename = basename(sv_plugin)

    command {
        set -euo pipefail

        ln -s ~{fasta} genome_reference.fasta

        mkdir -p /tmp/vep_cache

        cp ~{gnomad} /tmp/vep_cache
        cp ~{gnomad_idx} /tmp/vep_cache 

        tar -C /tmp/vep_cache -xvzf ~{vep_data}
        chmod 777 /tmp/vep_cache/homo_sapiens
        cp ~{fai} /tmp/vep_cache
        cp ~{fasta} /tmp/vep_cache
        cp ~{gzi} /tmp/vep_cache
        du -sh /tmp/vep_cache/Homo_sapiens_assembly38.fasta.gz*


        perl /opt/vep/src/ensembl-vep/vep --force_overwrite \
            --input_file ~{input_vcf} \
            --vcf \
            --output_file ~{sample_id}_sv_vep_annotated.vcf \
            --stats_file ~{sample_id}_sv_vep_stats.txt \
            --stats_text \
            --cache \
            --dir_cache /tmp/vep_cache \
            --fasta genome_reference.fasta \
            --fork ~{cpu} \
            --numbers --offline --hgvs --shift_hgvs 0 --terms SO --symbol \
            --sift b --polyphen b --total_length --ccds --canonical --biotype \
            --protein --xref_refseq --mane --pubmed --af --max_af --af_1kg --af_gnomadg \
            --max_sv_size ~{max_sv_size}

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
        File output_vep_vcf = "~{sample_id}_sv_vep_annotated.vcf"
        File output_vep_stats = "~{sample_id}_sv_vep_stats.txt"
    }
}

