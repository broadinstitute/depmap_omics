version 1.0

workflow post_GaTSV_SV_Workflow {
    input {
        File input_bedpe
        String sample_id

        File gtf_bed = "gs://ccleparams/gencode.v38.primary_assembly.CORRECTED_MISSING_IDs.annotation.GENES_ONLY.bed"
    }

    call reannotate_genes {
        input:
            input_bedpe = input_bedpe,
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

        # tweaked the task in vep_sv to work with GaTSV output format
        # annotate genes at breakpoint A
        sed '/^#/d' ~{input_bedpe} | \
            cut -f1-3,7 | \
            bedtools intersect -a stdin -b ~{gtf_bed} -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            bedtools groupby -g 1,2,3,4 -c 8 -o distinct | sort -k4,4 \
            > gene_overlaps.A.bed # collapse all genes to one row

        # annotate genes at breakpoint B
        sed '/^#/d' ~{input_bedpe} | \
            cut -f4-7 | \
            bedtools intersect -a stdin -b ~{gtf_bed} -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            bedtools groupby -g 1,2,3,4 -c 8 -o distinct | sort -k4,4 \
            > gene_overlaps.B.bed # collapse all genes to one row

        # join breakpoint A annotations and breakpoint B annotations
        join -1 4 -2 4 gene_overlaps.A.bed gene_overlaps.B.bed | \
            sed 's/ /\t/g' > gatsv_test.SV.gene_overlaps.txt

        # subset DEL and DUP variants
        awk -F"\t" '{ if ($8 == "DEL") print }' ~{input_bedpe} > ~{sample_id}.DEL.bedpe
        awk -F"\t" '{ if ($8 == "DUP") print }' ~{input_bedpe} > ~{sample_id}.DUP.bedpe

        # for DEL, not just look at genes at the breakpoints, but also genes that lie between breakpoints
        # so here intersect (start_a, end_b) with GTF
        sed '/^#/d'  ~{sample_id}.DEL.bedpe | \
            cut -f1,2,6,7 | \
            bedtools intersect -a stdin -b ~{gtf_bed} -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            bedtools groupby -g 1,2,3,4 -c 8 -o distinct | sort -k4,4 \
            > ~{sample_id}.DEL.overlap.bed # collapse all genes to one row

        # for DUP, not just look at genes at the breakpoints, but also genes that lie between breakpoints
        # so here intersect (start_a, end_b) with GTF
        sed '/^#/d'  ~{sample_id}.DUP.bedpe | \
            cut -f1,2,6,7 | \
            bedtools intersect -a stdin -b ~{gtf_bed} -wao | \
            sed 's/;//g' | sed 's/"//g' | \
            bedtools groupby -g 1,2,3,4 -c 8 -o distinct | sort -k4,4 \
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
        File cosmic_fusion_pairs = "gs://cds-cosmic/cosmic_fusion_gene_pairs_v100.csv"
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
            ~{cosmic_fusion_pairs} \
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
