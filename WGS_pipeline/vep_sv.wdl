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
        File gtf_bed = "gs://ccleparams/gencode.v38.primary_assembly.CORRECTED_MISSING_IDs.annotation.GENES_ONLY.bed"
    }

    call pr_sr_filter {
        input:
            input_vcf = input_vcf,
            sample_id = sample_id
    }

    call annotate_sv_vep {
        input:
            input_vcf = pr_sr_filter.output_filtered_vcf,
            sample_id = sample_id,
            fasta = fasta,
            fai = fai,
            gzi = gzi,
            vep_data = vep_data,
            gnomad = gnomad,
            gnomad_idx = gnomad_idx
    }

    call vcf2bedpe {
        input:
            input_vcf = annotate_sv_vep.output_vep_vcf,
            output_vcf_basename = sample_id

    }

    call gnomad_filter_qc {
        input:
            input_vcf = annotate_sv_vep.output_vep_vcf,
            sample_id = sample_id
    }

    call reannotate_genes {
        input:
            input_bedpe = vcf2bedpe.output_bedpe,
            sample_id = sample_id,
            gtf_bed = gtf_bed
    }

    call bedpe_to_depmap {
        input:
            input_bedpe = vcf2bedpe.output_bedpe,
            gene_annotation = reannotate_genes.output_reannotated_bedpe,
            del_annotation = reannotate_genes.annotated_overlap_del,
            dup_annotation = reannotate_genes.annotated_overlap_dup,
            sample_id = sample_id
    }

    output { 
        File vep_annotated_sv = annotate_sv_vep.output_vep_vcf
        File bedpe = vcf2bedpe.output_bedpe
        File expanded_sv_bedpe = bedpe_to_depmap.expanded_bedpe
        File expanded_filtered_sv_bedpe = bedpe_to_depmap.expanded_filtered_bedpe
        File reannotate_genes_bedpe = reannotate_genes.output_reannotated_bedpe
    }
}


task pr_sr_filter {
    # drop SVs where PR_alt + SR_alt < 3 OR FILTER is not PASS or MaxDepth
    input {
        File input_vcf
        String filter_statement = "SUM(FORMAT/PR[0:1]+FORMAT/SR[0:1]) > 3"
        String view_statement = 'FILTER="PASS"|FILTER="MaxDepth"'
        String sample_id

        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(10 * size(input_vcf, "GiB")) + 10

    command <<<
        set -euo pipefail

        bcftools filter --include '~{filter_statement}' ~{input_vcf} > ~{sample_id}.pr_sr_filtered.vcf
        bcftools view -i '~{view_statement}' ~{sample_id}.pr_sr_filtered.vcf > ~{sample_id}.pr_sr_filtered.passed.vcf
    >>>

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/bcftools:production"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {     
        File output_filtered_vcf = "~{sample_id}.pr_sr_filtered.passed.vcf"
    }

    meta {
        allowNestedInputs: true
    }
}


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
        File sv_plugin = "gs://cds-vep-data/StructuralVariantOverlap.pm"
        String assembly = "GRCh38"
        Int max_sv_size = 50000000

        Int cpu = 2
        Int mem_gb = 8
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
        Int boot_disk_size = 60
    }

    Int disk_space = ceil(
        10 * size(input_vcf, "GiB") + size([fasta, vep_data, gnomad, gzi], "GiB")
    ) + 10

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

        perl /opt/vep/src/ensembl-vep/vep --force_overwrite \
            --input_file ~{input_vcf} \
            --output_file ~{sample_id}_sv_vep_annotated.vcf \
            --cache \
            --fasta genome_reference.fasta \
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
        docker: "ensemblorg/ensembl-vep:release_112.0"
        bootDiskSizeGb: boot_disk_size
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {     
        File output_vep_vcf = "~{sample_id}_sv_vep_annotated.vcf"
    }

    meta {
        allowNestedInputs: true
    }
}

# stolen from https://github.com/hall-lab/sv-pipeline/blob/master/scripts/SV_Tasks.wdl
# converts SV vcf to bedpe
task vcf2bedpe {
    input {
        File input_vcf
        String output_vcf_basename

        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(3 * size(input_vcf, "GiB")) + 10

    command <<<
        set -eo pipefail

        bgzip "~{input_vcf}"
        zcat "~{input_vcf}.gz" \
            | svtools vcftobedpe \
            > "~{output_vcf_basename}.bedpe"
    >>>

    runtime {
        docker: "halllab/svtools:v0.5.1"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {
        File output_bedpe = "~{output_vcf_basename}.bedpe"
    }

    meta {
        allowNestedInputs: true
    }
}


task gnomad_filter_qc {
    # !!!!! this task is only for QC evaluations, NOT where the actual gnomad filtering happens !!!!!
    # !!!!! the actual gnomad filtering with rescue enabled happens in bedpe_to_depmap !!!!!
    input {
        File input_vcf
        String sample_id
        Float gnomad_cutoff = 0.001

        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(3 * size(input_vcf, "GiB")) + 10

    command <<<
        set -euo pipefail

        bcftools view -h ~{input_vcf} > ~{sample_id}.pr_sr_filtered.vep_annotated.gnomad_filtered.vcf

        # parse the SV_overlap_AF component (overlapping structural variant allele
        # frequency) in the CSQ annotation and remove variants where the highest
        # observed frequency is less than gnomad_cutoff
        awk -F"\t" '{
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
        }' ~{input_vcf} >> ~{sample_id}.pr_sr_filtered.vep_annotated.gnomad_filtered.vcf
    >>>

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/bcftools:production"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {     
        File output_filtered_vcf = "~{sample_id}.pr_sr_filtered.vep_annotated.gnomad_filtered.vcf"
    }

    meta {
        allowNestedInputs: true
    }
}

task reannotate_genes {
    # since VEP doesn't correctly annotate genes at breakpoints, we have to intersect them ourselves
    input {
        File input_bedpe
        String sample_id
        File gtf_bed

        Int mem_gb = 8
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(3 * size(input_bedpe, "GiB")) + ceil(size(gtf_bed, "GiB")) + 10

    command <<<
        set -euo pipefail

        # annotate genes at breakpoint A
        sed '/^#/d' ~{input_bedpe} | \
            cut -f1-3,13,16 | \
            bedtools intersect -a stdin -b ~{gtf_bed} -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            awk -F"\t" '{
                if ($5 != ".") {
                    split($4,arr,":");
                    print $1"\t"$2"\t"$3"\t"$4"\t"arr[1]":"arr[2]":"arr[3]":"arr[4]":"arr[5]":"arr[6]":"arr[7]"\t.\t"$(NF-1)
                }
                else {
                    print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t.\t"$9
                }
            }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | \
            bedtools groupby -g 1,2,3,4,5 -c 7 -o distinct | sort -k5,5 \
            > gene_overlaps.A.bed # collapse all genes to one row

        # annotate genes at breakpoint B
        sed '/^#/d' ~{input_bedpe} | \
            cut -f4,5,6,13,16 | \
            bedtools intersect -a stdin -b ~{gtf_bed} -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            awk -F"\t" '{
                if ($5 != ".") {
                    split($4,arr,":");
                    print $1"\t"$2"\t"$3"\t"$4"\t"arr[1]":"arr[2]":"arr[3]":"arr[4]":"arr[5]":"arr[6]":"arr[7]"\t.\t"$(NF-1)
                }
                else {
                    print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t.\t"$9
                }
            }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | \
            bedtools groupby -g 1,2,3,4,5 -c 7 -o distinct | sort -k5,5 \
            > gene_overlaps.B.bed # collapse all genes to one row

        # join breakpoint A annotations and breakpoint B annotations
        join -1 5 -2 5 gene_overlaps.A.bed gene_overlaps.B.bed | \
            sed 's/ /\t/g' | \
            cut -f2- > ~{sample_id}.SV.gene_overlaps.txt

        # subset DEL and DUP variants
        awk -F"\t" '{ if ($11 == "DEL") print }' ~{input_bedpe} > ~{sample_id}.DEL.bedpe
        awk -F"\t" '{ if ($11 == "DUP") print }' ~{input_bedpe} > ~{sample_id}.DUP.bedpe

        # for DEL, not just look at genes at the breakpoints, but also genes that lie between breakpoints
        # so here intersect (start_a, end_b) with GTF
        sed '/^#/d'  ~{sample_id}.DEL.bedpe | \
            cut -f1,2,6,13,16 | \
            bedtools intersect -a stdin -b ~{gtf_bed} -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            awk -F"\t" '{ \
                if ($5 != ".") {
                    split($4,arr,":");
                    print $1"\t"$2"\t"$3"\t"$4"\t"arr[1]":"arr[2]":"arr[3]":"arr[4]":"arr[5]":"arr[6]":"arr[7]"\t.\t"$(NF-1)
                }
                else {
                    print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t.\t"$9
                }
            }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | \
            bedtools groupby -g 1,2,3,4,5 -c 7 -o distinct | sort -k5,5 \
            > ~{sample_id}.DEL.overlap.bed

        # for DUP, not just look at genes at the breakpoints, but also genes that lie between breakpoints
        # so here intersect (start_a, end_b) with GTF
        sed '/^#/d'  ~{sample_id}.DUP.bedpe | \
            cut -f1,2,6,13,16 | \
            bedtools intersect -a stdin -b ~{gtf_bed} -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            awk -F"\t" '{
                if ($5 != ".") {
                    split($4,arr,":");
                    print $1"\t"$2"\t"$3"\t"$4"\t"arr[1]":"arr[2]":"arr[3]":"arr[4]":"arr[5]":"arr[6]":"arr[7]"\t.\t"$(NF-1)
                }
                else {
                    print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t.\t"$9
                }
            }' | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | \
            bedtools groupby -g 1,2,3,4,5 -c 7 -o distinct | sort -k5,5 \
            > ~{sample_id}.DUP.overlap.bed
    >>>

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/bedtools:production"
        disks: "local-disk ~{disk_space} SSD"
        memory: "~{mem_gb} GB"
        cpu: cpu
        preemptible: preemptible
    }

    output {     
        File output_reannotated_bedpe = "~{sample_id}.SV.gene_overlaps.txt"
        File annotated_overlap_del = "~{sample_id}.DEL.overlap.bed"
        File annotated_overlap_dup = "~{sample_id}.DUP.overlap.bed"
    }

    meta {
        allowNestedInputs: true
    }
}


task bedpe_to_depmap {
    # see python module for details, but in short, it includes:
    # expanding INFO annotation column, incorporating gene re-annotation generated above,
    # gnomad filering, rescuing, keeping only columns of interest
    input {
        File input_bedpe
        File gene_annotation
        File del_annotation
        File dup_annotation
        String sample_id

        Int mem_gb = 8
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(
        3 * size(input_bedpe, "GiB") + size([gene_annotation, del_annotation, dup_annotation], "GiB")
    ) + 10

    command <<<
        python -u /home/bedpe_to_depmap.py \
            ~{input_bedpe} \
            ~{gene_annotation} \
            ~{del_annotation} \
            ~{dup_annotation} \
            ~{sample_id}
    >>>

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/bedpe_to_depmap:production"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {
        File expanded_bedpe = "~{sample_id}.svs.expanded.reannotated.bedpe"
        File expanded_filtered_bedpe = "~{sample_id}.svs.expanded.reannotated.filtered.bedpe"
    }

    meta {
        allowNestedInputs: true
    }
}
