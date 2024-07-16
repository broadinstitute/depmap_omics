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

    call reannotate_genes {
        input:
            input_bedpe=vcf2bedpe.output_bedpe,
            sample_id=sample_id
    }

    call bedpe_to_depmap {
        input:
            input_bedpe=vcf2bedpe.output_bedpe,
            gene_annotation=reannotate_genes.output_reannotated_bedpe,
            sample_id=sample_id
    }

    call intersect_dels_and_dups {
        input:
            dels_bed=bedpe_to_depmap.dels,
            dups_bed=bedpe_to_depmap.dups,
            sample_id=sample_id
    }

    output { 
        File vep_annotated_sv = annotate_sv_vep.output_vep_vcf
        File bedpe = vcf2bedpe.output_bedpe
        File gnomad_filtered_no_rescue = gnomad_filter.output_filtered_vcf
        File expanded_sv_bedpe=bedpe_to_depmap.expanded_bedpe
        File expanded_filtered_sv_bedpe=bedpe_to_depmap.expanded_filtered_bedpe
        File reannotate_genes_bedpe=reannotate_genes.output_reannotated_bedpe
        File sv_del_genes = intersect_dels_and_dups.del_genes
        File sv_dup_genes = intersect_dels_and_dups.dup_genes
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
        File gtf="gs://depmap-omics-24q4-rna-ref/gencode.v38.primary_assembly.CORRECTED_MISSING_IDs.annotation.sorted.gtf.gz"
        File gtf_index="gs://depmap-omics-24q4-rna-ref/gencode.v38.primary_assembly.CORRECTED_MISSING_IDs.annotation.sorted.gtf.gz.tbi"
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
            --gtf ~{gtf} \
            --dir_cache /tmp/vep_cache --dir_plugins /tmp/Plugins \
            --fork 10 \
            --vcf \
            --dont_skip \
            --pick \
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

task reannotate_genes {
    input {
        File input_bedpe
        String sample_id
        File gtf_bed="gs://ccleparams/gencode.v38.primary_assembly.CORRECTED_MISSING_IDs.annotation.GENES_ONLY.bed"
        String docker_image="biocontainers/bedtools:v2.28.0_cv2"
        Int preemptible=2
        Int boot_disk_size=10
        Int disk_space=10
        Int cpu = 2
        Int mem = 10
    }

    command <<<
        set -euo pipefail

        sed '/^#/d' ~{input_bedpe} |\
        cut -f1-3,13,16 |\
        bedtools intersect -a stdin -b ~{gtf_bed} -wao | \
        sed 's/;//g' | sed 's/"//g' | \
        awk -F"\t" '{ \
                    if ($5 != ".") { \
                    split($4,arr,":");  \
                    print $1"\t"$2"\t"$3"\t"$4"\t"arr[1]":"arr[2]":"arr[3]":"arr[4]":"arr[5]":"arr[6]":"arr[7]"\t.\t"$(NF-1) \
                    } \
                    else { \
                    print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t.\t"$9 \
                    } \
        }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | \
        bedtools groupby -g 1,2,3,4,5 -c 7 -o distinct | sort -k5,5 > gene_overlaps.A.bed


        sed '/^#/d' ~{input_bedpe} | \
        cut -f4,5,6,13,16 | \
        bedtools intersect -a stdin -b ~{gtf_bed} -wao | \
        sed 's/;//g' | sed 's/"//g' | \
        awk -F"\t" '{ \
                    if ($5 != ".") { \
                    split($4,arr,":");  \
                print $1"\t"$2"\t"$3"\t"$4"\t"arr[1]":"arr[2]":"arr[3]":"arr[4]":"arr[5]":"arr[6]":"arr[7]"\t.\t"$(NF-1) \
                } \
                else { \
                print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t.\t"$9 \
                } \
        }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | \
        bedtools groupby -g 1,2,3,4,5 -c 7 -o distinct | sort -k5,5 > gene_overlaps.B.bed

        join -1 5 -2 5 gene_overlaps.A.bed gene_overlaps.B.bed  | sed 's/ /\t/g' | cut -f2- > ~{sample_id}.SV.gene_overlaps.txt

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
        File output_reannotated_bedpe = "~{sample_id}.SV.gene_overlaps.txt"
    }
}


task bedpe_to_depmap {
    input {
        File input_bedpe
        File gene_annotation
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
              ~{gene_annotation} \
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
        File dels = "~{sample_id}_dels.bed"
        File dups = "~{sample_id}_dups.bed"
    }
}

task intersect_dels_and_dups {
    input {
        File dels_bed
        File dups_bed
        String sample_id
        File gtf_bed="gs://ccleparams/gencode.v38.primary_assembly.CORRECTED_MISSING_IDs.annotation.GENES_ONLY.HUGO_ONLY.bed"
        String docker_image="biocontainers/bedtools:v2.28.0_cv2"
        Int preemptible=2
        Int boot_disk_size=10
        Int disk_space=10
        Int cpu = 2
        Int mem = 10
    }

    command <<<
        set -euo pipefail

        bedtools intersect -a ~{dels_bed} -b ~{gtf_bed} -wao | bedtools groupby -g 1,2,3 -c 7 -o distinct > ~{sample_id}.del_genes.bed
        bedtools intersect -a ~{dups_bed} -b ~{gtf_bed} -wao | bedtools groupby -g 1,2,3 -c 7 -o distinct > ~{sample_id}.dup_genes.bed
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
        File del_genes = "~{sample_id}.del_genes.bed"
        File dup_genes = "~{sample_id}.dup_genes.bed"
    }
}