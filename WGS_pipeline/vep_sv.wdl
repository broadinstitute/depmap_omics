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

    call gnomad_filter {
        input:
            input_vcf=annotate_sv_vep.output_vep_vcf,
            sample_id=sample_id,
    }

    call vcf2bedpe {
        input:
            input_vcf=gnomad_filter.output_filtered_vcf,
            output_vcf_basename=sample_id,

    }

    output { 
        File vep_annotated_sv = annotate_sv_vep.output_vep_vcf
        File vep_sv_stats = annotate_sv_vep.output_vep_stats
        File gnomad_filtered_sv = gnomad_filter.output_filtered_vcf
        File gnomad_filtered_bedpe = vcf2bedpe.output_bedpe
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
        Int gnomad_distance=500
        String sv_match_type="overlap" # overlap (default), within, surrounding or exact

    }

    String gnomad_basename = basename(gnomad)
    String plugin_basename = basename(sv_plugin)

    command {
        set -euo pipefail

        ln -s ~{fasta} genome_reference.fasta

        mkdir -p /tmp/vep_cache

        cp ~{gnomad} /tmp/vep_cache
        cp ~{gnomad_idx} /tmp/vep_cache 

        tar -C /tmp/vep_cache -xzf ~{vep_data}
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
            --dont_skip \
            --dir_cache /tmp/vep_cache \
            --fasta genome_reference.fasta \
            --fork ~{cpu} \
            --pick \
            --numbers --offline --hgvs --shift_hgvs 0 --terms SO --symbol \
            --sift b --polyphen b --total_length --ccds --canonical --biotype \
            --protein --xref_refseq --mane --pubmed --af --max_af --af_1kg --af_gnomadg \
            --custom file=/tmp/vep_cache/~{gnomad_basename},short_name=gnomAD_SV,format=vcf,type=~{sv_match_type},overlap_cutoff=~{sv_match_percentage},distance=~{gnomad_distance},same_type=1,fields=AC%AF%FILTER,summary_stats=min \
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

        bcftools filter --include 'SUM(FORMAT/PR[0:1]+FORMAT/SR[0:1]) > 3' ~{input_vcf} > ~{sample_id}.vep_annotated.pr_sr_filtered.vcf

        bcftools view -h ~{input_vcf} > ~{sample_id}.vep_annotated.pr_sr_filtered.gnomad_filtered.vcf

        bcftools view --with-header ~{sample_id}.vep_annotated.pr_sr_filtered.vcf | awk -F"\t" '{
            split($8, info, ";");
            for (i=1; i<=length(info); i++) {
                if (info[i] ~ /^CSQ=/) {
                    split(substr(info[i], 5), csq, "|");
                    if (csq[length(csq)-2] == "") {
                            print $0
                    } else {
                        split(csq[length(csq)-2], af, "&");
                        min=1
                        for (j=1; j<=length(af); j++) {
                            if (af[j] <= min) {
                                min = af[j]
                            }
                        }
                        if (min <= ~{gnomad_cutoff}) {
                            print $0
                        }
                    }
                }
            }
        }' >> ~{sample_id}.vep_annotated.gnomad_filtered.vcf
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
        File output_filtered_vcf = "~{sample_id}.vep_annotated.gnomad_filtered.vcf"
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