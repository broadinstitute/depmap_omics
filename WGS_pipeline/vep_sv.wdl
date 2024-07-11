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

    call pr_sr_filter {
        input:
            input_vcf=input_vcf,
            sample_id=sample_id,
    }

    call annotate_sv_vep {
        input:
            input_vcf=pr_sr_filter.output_filtered_vcf,
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

    call vcf2bedpe {
        input:
            input_vcf=annotate_sv_vep.output_vep_vcf,
            output_vcf_basename=sample_id,

    }

    call gnomad_filter {
        input:
            input_vcf=annotate_sv_vep.output_vep_vcf,
            sample_id=sample_id,
    }

    call bedpe_to_depmap {
        input:
            input_bedpe=vcf2bedpe.output_bedpe,
            sample_id=sample_id
    }

    output { 
        File vep_annotated_sv = annotate_sv_vep.output_vep_vcf
        File bedpe = vcf2bedpe.output_bedpe
        File gnomad_filtered_no_rescue = gnomad_filter.output_filtered_vcf
        File expanded_sv_bedpe=bedpe_to_depmap.expanded_bedpe
        File expanded_filtered_sv_bedpe=bedpe_to_depmap.expanded_filtered_bedpe
    }
}


task pr_sr_filter {
    input {
        File input_vcf
        String sample_id
        String docker_image="dceoy/bcftools"
        Int preemptible=2
        Int boot_disk_size=10
        Int disk_space=10
        Int cpu = 2
        Int mem = 10
    }

    command <<<
        set -euo pipefail

        bcftools filter --include 'SUM(FORMAT/PR[0:1]+FORMAT/SR[0:1]) > 3' ~{input_vcf} > ~{sample_id}.pr_sr_filtered.vcf
        bcftools view -i 'FILTER="PASS"|FILTER="MaxDepth"' ~{sample_id}.pr_sr_filtered.vcf > ~{sample_id}.pr_sr_filtered.passed.vcf
    >>>

    runtime {
        disks: "local-disk ~{disk_space} HDD"
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {     
        File output_filtered_vcf = "~{sample_id}.pr_sr_filtered.passed.vcf"
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
        Int gnomad_distance=0
        String sv_match_type="overlap" # overlap (default), within, surrounding or exact
        Int sv_match_reciprocal=0

    }

    String gnomad_basename = basename(gnomad)
    String plugin_basename = basename(sv_plugin)

    command {
        set -euo pipefail

        ln -s ~{fasta} genome_reference.fasta

        mkdir -p /tmp/vep_cache

        mkdir -p /tmp/Plugins
        curl https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/StructuralVariantOverlap.pm -o /tmp/Plugins/StructuralVariantOverlap.pm
        cp ~{gnomad} /tmp/Plugins
        cp ~{gnomad_idx} /tmp/Plugins 

        tar -C /tmp/vep_cache -xzf ~{vep_data}
        chmod 777 /tmp/vep_cache/homo_sapiens
        cp ~{fai} /tmp/vep_cache
        cp ~{fasta} /tmp/vep_cache
        cp ~{gzi} /tmp/vep_cache
        du -sh /tmp/vep_cache/Homo_sapiens_assembly38.fasta.gz*

        perl /opt/vep/src/ensembl-vep/vep --force_overwrite \
            --input_file ~{input_vcf} \
            --output_file ~{sample_id}_sv_vep_annotated.vcf \
            --cache \
            --fasta genome_reference.fasta \
            --dir_cache /tmp/vep_cache --dir_plugins /tmp/Plugins \
            --fork 10 \
            --vcf \
            --dont_skip \
            # --pick \
            --numbers --offline --hgvs --shift_hgvs 0 --terms SO --symbol --mane \
            --total_length --ccds --canonical --biotype --protein \
            --max_sv_size ~{max_sv_size} \
            --plugin StructuralVariantOverlap,file=/tmp/Plugins/~{gnomad_basename},same_type=1,overlap_cutoff=80,reciprocal=0,same_type=1,fields=AC%AF

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
    }
}

# stolen from https://github.com/hall-lab/sv-pipeline/blob/master/scripts/SV_Tasks.wdl
task vcf2bedpe {
    input {
        File input_vcf
        String output_vcf_basename
        Int preemptible_tries=2
    }

    command {
        set -eo pipefail

        bgzip ~{input_vcf} 
        zcat ${input_vcf}.gz \
            | svtools vcftobedpe \
            > ${output_vcf_basename}.bedpe
    }

    runtime {
        docker: "halllab/svtools:v0.5.1"
        cpu: "1"
        memory: "3 GB"
        disks: "local-disk " +  3*ceil( size(input_vcf, "GB")) + " HDD"
        preemptible: preemptible_tries
    }

    output {
        File output_bedpe = "${output_vcf_basename}.bedpe"
    }
}

task gnomad_filter {
    input {
        File input_vcf
        String sample_id
        Float gnomad_cutoff=0.001
        String docker_image="dceoy/bcftools"
        Int preemptible=2
        Int boot_disk_size=10
        Int disk_space=10
        Int cpu = 2
        Int mem = 10
    }

    command <<<
        set -euo pipefail

        bcftools view -h ~{input_vcf} > ~{sample_id}.pr_sr_filtered.vep_annotated.gnomad_filtered.vcf

        bcftools view ~{input_vcf} | awk -F"\t" '{
            split($8, info, ";");
            for (i=1; i<=length(info); i++) {
                if (info[i] ~ /^CSQ=/) {
                    split(substr(info[i], 5), csq, "|");
                    if (csq[length(csq)-2] == "") {
                            print $0
                    } else {
                        split(csq[length(csq)-2], af, "&");
                        max=0
                        for (j=1; j<=length(af); j++) {
                            if (af[j] >= max) {
                                max = af[j]
                            }
                        }
                        if (max <= ~{gnomad_cutoff}) {
                            print $0
                        }
                    }
                }
            }
        }' >> ~{sample_id}.pr_sr_filtered.vep_annotated.gnomad_filtered.vcf
    >>>

    runtime {
        disks: "local-disk ~{disk_space} HDD"
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {     
        File output_filtered_vcf = "~{sample_id}.pr_sr_filtered.vep_annotated.gnomad_filtered.vcf"
    }
}


task bedpe_to_depmap {
    input {
        File input_bedpe
        String sample_id

        String docker_image="us-docker.pkg.dev/depmap-omics/public/bedpe_to_depmap:test"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=40
        Int cpu = 4
        Int mem = 32
    }

    command {
        python -u /home/bedpe_to_depmap.py \
              ~{input_bedpe} \
              ~{sample_id} 
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
        File expanded_bedpe = "~{sample_id}.svs.expanded.bedpe"
        File expanded_filtered_bedpe = "~{sample_id}.svs.expanded.filtered.bedpe"
    }
}
