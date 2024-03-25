version 1.0

workflow annotate_variants_wf {
    
    input {
        File input_vcf
        File input_vcf_index
        String base_name

        File ref_fasta
        File ref_fasta_index

        File? bam
        File? bam_index

        # annotation options
        Boolean incl_snpEff = true
        Boolean incl_dbsnp = true
        Boolean incl_gnomad = true
        Boolean incl_rna_editing = true
        Boolean include_read_var_pos_annotations = true
        Boolean incl_repeats = true
        Boolean incl_homopolymers = true
        Boolean incl_splice_dist = true
        Boolean incl_blat_ED = true
        Boolean incl_cosmic = true
        Boolean incl_cravat = true

        File? dbsnp_vcf
        File? dbsnp_vcf_index

        File? gnomad_vcf
        File? gnomad_vcf_index

        File? rna_editing_vcf
        File? rna_editing_vcf_index

        File? cosmic_vcf
        File? cosmic_vcf_index

        File? ref_splice_adj_regions_bed
        
        File? cravat_lib_tar_gz
        String? cravat_lib_dir

        File? repeat_mask_bed

        String? genome_version

        String docker = "trinityctat/ctat_mutations:latest"
        String plugins_path = "/usr/local/src/ctat-mutations/plugins"
        String scripts_path = "/usr/local/src/ctat-mutations/src"

        Int preemptible
        Int cpu
    }


    Boolean vcf_input = defined(vcf)

    parameter_meta {

        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}
        bam:{help:"Previously aligned bam file. When VCF is provided, the output from ApplyBQSR should be provided as the bam input."}
        bai:{help:"Previously aligned bam index file"}
        vcf:{help:"Previously generated vcf file to annotate and filter. When provided, the output from ApplyBQSR should be provided as the bam input."}
        sample_id:{help:"Sample id"}

        # resources
        ref_fasta:{help:"Path to the reference genome to use in the analysis pipeline."}
        ref_fasta_index:{help:"Index for ref_fasta"}
        gtf:{help:"Annotations GTF."}


        dbsnp_vcf:{help:"dbSNP VCF file for the reference genome."}
        dbsnp_vcf_index:{help:"dbSNP vcf index"}

        gnomad_vcf:{help:"gnomad vcf file w/ allele frequencies"}
        gnomad_vcf_index:{help:"gnomad VCF index"}

        rna_editing_vcf:{help:"RNA editing VCF file"}
        rna_editing_vcf_index:{help:"RNA editing VCF index"}

        cosmic_vcf:{help:"Coding Cosmic Mutation VCF annotated with Phenotype Information"}
        cosmic_vcf_index:{help:"COSMIC VCF index"}

        repeat_mask_bed:{help:"bed file containing repetive element (repeatmasker) annotations (from UCSC genome browser download)"}

        ref_splice_adj_regions_bed:{help:"For annotating exon splice proximity"}

        cravat_lib_tar_gz:{help:"CRAVAT resource archive"}
        cravat_lib_dir:{help:"CRAVAT resource directory (for non-Terra / local use where ctat genome lib is installed already)"}

        genome_version:{help:"Genome version for annotating variants using Cravat and SnpEff", choices:["hg19", "hg38"]}

        include_read_var_pos_annotations :{help: "Add vcf annotation that requires variant to be at least 6 bases from ends of reads."}

        plugins_path:{help:"Path to plugins"}
        scripts_path:{help:"Path to scripts"}

        docker:{help:"Docker or singularity image"}
    }


    # leftnorm and split multiallelics
    call left_norm_vcf {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            base_name = base_name,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            scripts_path = scripts_path,
            docker = docker,
            preemptible = preemptible
    }


    if (incl_snpEff) {

        call snpEff {
            input:
                input_vcf = left_norm_vcf.vcf,
                input_vcf_index = left_norm_vcf.vcf_index,
                base_name = base_name,
                scripts_path = scripts_path,
                plugins_path = plugins_path,
                genome_version = select_first([genome_version]),
                docker = docker,
                preemptible = preemptible,
                cpu = cpu
        }
    }


    if (incl_dbsnp) {

        call annotate_dbsnp {
            input:
                input_vcf =  select_first([snpEff.vcf, left_norm_vcf.vcf]),
                input_vcf_index = select_first([snpEff.vcf_index, left_norm_vcf.vcf_index]),
                dbsnp_vcf = select_first([dbsnp_vcf]),
                dbsnp_vcf_index = select_first([dbsnp_vcf_index]),
                base_name = base_name,
                docker = docker,
                preemptible = preemptible
        }
    }

    if (incl_gnomad) {

        call annotate_gnomad {
            input:
                input_vcf = select_first([annotate_dbsnp.vcf, snpEff.vcf, left_norm_vcf.vcf]),
                input_vcf_index = select_first([annotate_dbsnp.vcf_index, snpEff.vcf_index, left_norm_vcf.vcf_index]),
                gnomad_vcf = select_first([gnomad_vcf]),
                gnomad_vcf_index = select_first([gnomad_vcf_index]),
                base_name = base_name,
                docker = docker,
                preemptible = preemptible
        }
    }
                

    if (incl_rna_editing) {

        call annotate_RNA_editing {
            input:
                input_vcf = select_first([annotate_gnomad.vcf, annotate_dbsnp.vcf, snpEff.vcf, left_norm_vcf.vcf]),
                input_vcf_index = select_first([annotate_gnomad.vcf_index, annotate_dbsnp.vcf_index, snpEff.vcf_index, left_norm_vcf.vcf_index]),
                rna_editing_vcf = select_first([rna_editing_vcf]),
                rna_editing_vcf_index = select_first([rna_editing_vcf_index]),
                base_name = base_name,
                scripts_path = scripts_path,
                docker = docker,
                preemptible = preemptible
        }
    }

    if (include_read_var_pos_annotations) {
        
        call annotate_PASS_reads {
            input:
                input_vcf = select_first([annotate_RNA_editing.vcf, annotate_gnomad.vcf, annotate_dbsnp.vcf, snpEff.vcf, left_norm_vcf.vcf]),
                input_vcf_index = select_first([annotate_RNA_editing.vcf_index, annotate_gnomad.vcf_index, annotate_dbsnp.vcf_index, snpEff.vcf_index, left_norm_vcf.vcf_index]),
                bam = select_first([bam]),
                bam_index = select_first([bam_index]),
                base_name = base_name,
                scripts_path = scripts_path,
                docker = docker,
                preemptible = preemptible,
                cpu = cpu
        }
    }


    if (incl_repeats) {
    
        call annotate_repeats {
            input:
                input_vcf = select_first([annotate_PASS_reads.vcf, annotate_RNA_editing.vcf, annotate_gnomad.vcf, annotate_dbsnp.vcf, snpEff.vcf, left_norm_vcf.vcf]),
                input_vcf_index = select_first([annotate_PASS_reads.vcf_index, annotate_RNA_editing.vcf_index, annotate_gnomad.vcf_index, annotate_dbsnp.vcf_index, snpEff.vcf_index, left_norm_vcf.vcf_index]),
                repeat_mask_bed = select_first([repeat_mask_bed]),
                base_name = base_name,
                scripts_path = scripts_path,
                docker = docker,
                preemptible = preemptible,
                cpu = cpu
       }
    }


    if (incl_homopolymers) {

        call annotate_homopolymers_n_entropy {
            input:
                input_vcf = select_first([annotate_repeats.vcf, annotate_PASS_reads.vcf, annotate_RNA_editing.vcf, annotate_gnomad.vcf, annotate_dbsnp.vcf, snpEff.vcf, left_norm_vcf.vcf]), 
                input_vcf_index = select_first([annotate_repeats.vcf_index, annotate_PASS_reads.vcf_index, annotate_RNA_editing.vcf_index, annotate_gnomad.vcf_index, annotate_dbsnp.vcf_index, snpEff.vcf_index, left_norm_vcf.vcf_index]),
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                base_name = base_name,
                scripts_path = scripts_path,
                docker = docker,
                preemptible = preemptible
          }
    }
                
    if (incl_splice_dist) {
        call annotate_splice_distance {
            input:
                input_vcf = select_first([annotate_homopolymers_n_entropy.vcf, annotate_repeats.vcf, annotate_PASS_reads.vcf, annotate_RNA_editing.vcf, annotate_gnomad.vcf, annotate_dbsnp.vcf, snpEff.vcf, left_norm_vcf.vcf]),
                input_vcf_index = select_first([annotate_homopolymers_n_entropy.vcf_index, annotate_repeats.vcf_index, annotate_PASS_reads.vcf_index, annotate_RNA_editing.vcf_index, annotate_gnomad.vcf_index, annotate_dbsnp.vcf_index, snpEff.vcf_index, left_norm_vcf.vcf_index]),
                ref_splice_adj_regions_bed = select_first([ref_splice_adj_regions_bed]),
                base_name = base_name,
                scripts_path = scripts_path,
                docker = docker,
                preemptible = preemptible
         }
    }


    if (incl_blat_ED) {
        call annotate_blat_ED {
            input:
                input_vcf = select_first([annotate_splice_distance.vcf, annotate_homopolymers_n_entropy.vcf, annotate_repeats.vcf, annotate_PASS_reads.vcf, annotate_RNA_editing.vcf, annotate_gnomad.vcf, annotate_dbsnp.vcf, snpEff.vcf, left_norm_vcf.vcf]),
                input_vcf_index = select_first([annotate_splice_distance.vcf_index,  annotate_homopolymers_n_entropy.vcf_index, annotate_repeats.vcf_index, annotate_PASS_reads.vcf_index, annotate_RNA_editing.vcf_index, annotate_gnomad.vcf_index, annotate_dbsnp.vcf_index, snpEff.vcf_index, left_norm_vcf.vcf_index]),
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                base_name = base_name,
                scripts_path = scripts_path,
                docker = docker,
                preemptible = preemptible,
                cpu = cpu
        }

    }


    if (incl_cosmic) {
        call annotate_cosmic_variants {
            input:
                input_vcf = select_first([annotate_blat_ED.vcf, annotate_splice_distance.vcf, annotate_homopolymers_n_entropy.vcf, annotate_repeats.vcf, annotate_PASS_reads.vcf, annotate_RNA_editing.vcf, annotate_gnomad.vcf, annotate_dbsnp.vcf, snpEff.vcf, left_norm_vcf.vcf]),
                input_vcf_index = select_first([annotate_blat_ED.vcf_index, annotate_splice_distance.vcf_index,  annotate_homopolymers_n_entropy.vcf_index, annotate_repeats.vcf_index, annotate_PASS_reads.vcf_index, annotate_RNA_editing.vcf_index, annotate_gnomad.vcf_index, annotate_dbsnp.vcf_index, snpEff.vcf_index, left_norm_vcf.vcf_index]), 
                cosmic_vcf = select_first([cosmic_vcf]),
                cosmic_vcf_index = select_first([cosmic_vcf_index]),
                base_name = base_name,
                scripts_path = scripts_path,
                docker = docker,
                preemptible = preemptible
        }

    }

    if (incl_cravat) {
         
        call open_cravat {
            input:
                input_vcf = select_first([annotate_cosmic_variants.vcf, annotate_blat_ED.vcf, annotate_splice_distance.vcf, annotate_homopolymers_n_entropy.vcf, annotate_repeats.vcf, annotate_PASS_reads.vcf, annotate_RNA_editing.vcf, annotate_gnomad.vcf, annotate_dbsnp.vcf, snpEff.vcf, left_norm_vcf.vcf]),
                input_vcf_index = select_first([annotate_cosmic_variants.vcf_index, annotate_blat_ED.vcf_index, annotate_splice_distance.vcf_index,  annotate_homopolymers_n_entropy.vcf_index, annotate_repeats.vcf_index, annotate_PASS_reads.vcf_index, annotate_RNA_editing.vcf_index, annotate_gnomad.vcf_index, annotate_dbsnp.vcf_index, snpEff.vcf_index, left_norm_vcf.vcf_index]),
                cravat_lib_tar_gz = cravat_lib_tar_gz,
                cravat_lib_dir = cravat_lib_dir,
                genome_version = select_first([genome_version]),
                base_name = base_name,
                scripts_path = scripts_path,
                docker = docker,
                preemptible = preemptible
        }
    }


    call rename_vcf {
        input:
            input_vcf = select_first([open_cravat.vcf, annotate_cosmic_variants.vcf, annotate_blat_ED.vcf, annotate_splice_distance.vcf, annotate_homopolymers_n_entropy.vcf, annotate_repeats.vcf, annotate_PASS_reads.vcf, annotate_RNA_editing.vcf, annotate_gnomad.vcf, annotate_dbsnp.vcf, snpEff.vcf, left_norm_vcf.vcf]),
            input_vcf_index = select_first([open_cravat.vcf_index, annotate_cosmic_variants.vcf_index, annotate_blat_ED.vcf_index, annotate_splice_distance.vcf_index,  annotate_homopolymers_n_entropy.vcf_index, annotate_repeats.vcf_index, annotate_PASS_reads.vcf_index, annotate_RNA_editing.vcf_index, annotate_gnomad.vcf_index, annotate_dbsnp.vcf_index, snpEff.vcf_index, left_norm_vcf.vcf_index]),
            base_name = base_name,
            docker=docker,
            preemptible=preemptible

    }


   output {
    
        File vcf = rename_vcf.vcf
        File vcf_index = rename_vcf.vcf_index
    
    }
}


task left_norm_vcf {
    input {
        File input_vcf
        File input_vcf_index
        String base_name
        File ref_fasta
        File ref_fasta_index
        String scripts_path

        String docker
        Int preemptible
        Int cpu = 1
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
        
    }

    command <<<
        set -ex

        echo "####### leftnorm and split multiallelics ##########"
      
        # leftnorm and split multiallelics
        bcftools norm \
        -f ~{ref_fasta} \
        -m -any \
        -o ~{base_name}.norm.vcf \
        ~{input_vcf}

        ~{scripts_path}/groom_vcf.py ~{base_name}.norm.vcf ~{base_name}.norm.groom.vcf

        bcftools sort -T . ~{base_name}.norm.groom.vcf > ~{base_name}.norm.groom.sorted.vcf
        bgzip -c ~{base_name}.norm.groom.sorted.vcf > ~{base_name}.norm.groom.sorted.vcf.gz
        tabix ~{base_name}.norm.groom.sorted.vcf.gz

   >>>

    output {
        File vcf = "~{base_name}.norm.groom.sorted.vcf.gz"
        File vcf_index = "~{base_name}.norm.groom.sorted.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }


}


task snpEff {

    input {
        File input_vcf
        File input_vcf_index
        String base_name
        String scripts_path
        String base_name
        String plugins_path
        String genome_version

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
   
    }


    command <<<

        set -ex


        echo "######### SnpEFF #########"
      
        bgzip -cd ~{input_vcf} | \
            java -Xmx3500m -jar ~{plugins_path}/snpEff.jar \
            -nostats -noLof -no-downstream -no-upstream -noLog \
            ~{genome_version} > ~{base_name}.snpeff.tmp.vcf

        ~{scripts_path}/update_snpeff_annotations.py \
            ~{base_name}.snpeff.tmp.vcf \
            ~{base_name}.snpeff.vcf

        bgzip -c ~{base_name}.snpeff.vcf > ~{base_name}.snpeff.vcf.gz
        tabix ~{base_name}.snpeff.vcf.gz
           

    >>>

    output {
        File vcf = "~{base_name}.snpeff.vcf.gz"
        File vcf_index = "~{base_name}.snpeff.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_dbsnp {

    input {
        File input_vcf
        File input_vcf_index
        File dbsnp_vcf
        File dbsnp_vcf_index
        String base_name

        String docker
        Int preemptible
        Int cpu = 1
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
    }

    command <<<
        set -ex

        echo "####### Annotate dbSNP ########"
      
        bcftools annotate \
            --output-type z \
            --annotations ~{dbsnp_vcf} \
            --columns "INFO/OM,INFO/PM,INFO/SAO,INFO/RS" \
            --output ~{base_name}.dbsnp.vcf.gz \
            ~{input_vcf}
            tabix ~{base_name}.dbsnp.vcf.gz

    >>>

    output {
        File vcf = "~{base_name}.dbsnp.vcf.gz"
        File vcf_index = "~{base_name}.dbsnp.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_gnomad {

    input {
        File input_vcf
        File input_vcf_index
        File gnomad_vcf
        File gnomad_vcf_index
        String base_name

        String docker
        Int preemptible
        Int cpu = 1
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
    }

    command <<<
        set -ex

        echo "####### Annotate gnomAD ########"

      
        bcftools annotate \
            --output-type z \
            --annotations ~{gnomad_vcf} \
            --columns "INFO/gnomad_RS,INFO/gnomad_AF" \
            --output ~{base_name}.gnomad.vcf.gz \
            ~{input_vcf}

            tabix ~{base_name}.gnomad.vcf.gz

    >>>

    output {
        File vcf = "~{base_name}.gnomad.vcf.gz"
        File vcf_index = "~{base_name}.gnomad.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_RNA_editing {

    input {
        File input_vcf
        File input_vcf_index
        File rna_editing_vcf
        File rna_editing_vcf_index
        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu = 1
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
    }

    command <<<
        set -ex


         echo "######### Annotate RNA Editing #############"
      
         bcftools annotate \
            --output-type z \
            --annotations ~{rna_editing_vcf} \
            --columns "INFO/RNAEDIT" \
            --output ~{base_name}.rna_editing.gz \
            ~{input_vcf}

        #must groom for gatk compat
        ~{scripts_path}/groom_vcf.py ~{base_name}.rna_editing.gz ~{base_name}.rna_editing.groom.vcf


        bgzip -c ~{base_name}.rna_editing.groom.vcf > ~{base_name}.rna_editing.groom.vcf.gz
        tabix ~{base_name}.rna_editing.groom.vcf.gz


    >>>

    output {
        File vcf = "~{base_name}.rna_editing.groom.vcf.gz"
        File vcf_index = "~{base_name}.rna_editing.groom.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_PASS_reads {

    input {
        File input_vcf
        File input_vcf_index
        File bam
        File bam_index
        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(bam, "GB") * 4) + (size(input_vcf, "GB") * 10) + 20)
        
    }


    command <<<

        set -ex

        echo "######## Annotate PASS Reads #########"
      
        samtools index ~{bam}

        ~{scripts_path}/annotate_PASS_reads.py \
            --vcf ~{input_vcf}  \
            --bam ~{bam} \
            --output_vcf ~{base_name}.annot_pass_reads.vcf \
            --threads ~{cpu}

            bgzip -c ~{base_name}.annot_pass_reads.vcf > ~{base_name}.annot_pass_reads.vcf.gz
            tabix ~{base_name}.annot_pass_reads.vcf.gz

    >>>


    output {
        File vcf = "~{base_name}.annot_pass_reads.vcf.gz"
        File vcf_index = "~{base_name}.annot_pass_reads.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_repeats {
    input {
        File input_vcf
        File input_vcf_index

        File repeat_mask_bed

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 20)
       
    }

    command <<<

        set -ex

        echo "####### Annotate Repeats #########"

      
        ~{scripts_path}/annotate_repeats.py \
            --input_vcf ~{input_vcf} \
            --repeats_bed ~{repeat_mask_bed} \
            --output_vcf ~{base_name}.annot_repeats.vcf

            bgzip -c ~{base_name}.annot_repeats.vcf > ~{base_name}.annot_repeats.vcf.gz

            tabix ~{base_name}.annot_repeats.vcf.gz
    >>>


    output {
        File vcf = "~{base_name}.annot_repeats.vcf.gz"
        File vcf_index = "~{base_name}.annot_repeats.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_homopolymers_n_entropy {

    input {
        File input_vcf
        File input_vcf_index

        File ref_fasta
        File ref_fasta_index

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu = 1
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
    }

    command <<<

        set -ex

         echo "########## Annotate Entropy and Homopolymers ###############"
      
         ~{scripts_path}/annotate_entropy_n_homopolymers.py \
            --input_vcf ~{input_vcf} \
            --ref_genome_fa ~{ref_fasta} \
            --tmpdir $TMPDIR \
            --output_vcf ~{base_name}.homopolymer.vcf

            bgzip -c ~{base_name}.homopolymer.vcf > ~{base_name}.homopolymer.vcf.gz
            tabix ~{base_name}.homopolymer.vcf.gz

    >>>


    output {
        File vcf = "~{base_name}.homopolymer.vcf.gz"
        File vcf_index = "~{base_name}.homopolymer.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_splice_distance {

    input {
        File input_vcf
        File input_vcf_index

        File ref_splice_adj_regions_bed

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu = 1
        Int disk = ceil((size(input_vcf, "GB") * 3) + 20)

    }


    command <<<

        set -ex


        echo "########## Annotate Splice Distance ##############"
      
        ~{scripts_path}/annotate_DJ.py \
            --input_vcf ~{input_vcf} \
            --splice_bed ~{ref_splice_adj_regions_bed} \
            --temp_dir $TMPDIR \
            --output_vcf ~{base_name}.splice_distance.vcf

        bgzip -c ~{base_name}.splice_distance.vcf > ~{base_name}.splice_distance.vcf.gz
        tabix ~{base_name}.splice_distance.vcf.gz

    >>>


    output {
        File vcf = "~{base_name}.splice_distance.vcf.gz"
        File vcf_index = "~{base_name}.splice_distance.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}    


task annotate_blat_ED {

    input {
        File input_vcf
        File input_vcf_index

        File ref_fasta
        File ref_fasta_index

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 200)

    }


    command <<<

        set -ex
      
        echo "########### Annotate BLAT ED #############"
      
        ~{scripts_path}/annotate_ED.py \
            --input_vcf ~{input_vcf} \
            --output_vcf ~{base_name}.blat_ED.vcf \
            --reference ~{ref_fasta} \
            --temp_dir $TMPDIR \
            --threads ~{cpu}

            bgzip -c ~{base_name}.blat_ED.vcf > ~{base_name}.blat_ED.vcf.gz
            tabix ~{base_name}.blat_ED.vcf.gz

    >>>


    output {
        File vcf = "~{base_name}.blat_ED.vcf.gz"
        File vcf_index = "~{base_name}.blat_ED.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "40G"
        preemptible: preemptible
        cpu : cpu
    }

}    




task annotate_cosmic_variants {

    input {
        File input_vcf
        File input_vcf_index

        File cosmic_vcf
        File cosmic_vcf_index

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu = 1
        Int disk = ceil((size(input_vcf, "GB") * 3) + 20)

    }


    command <<<

        set -ex

        echo "############# Annotate COSMIC Variants ################"
      
        bcftools annotate \
            --annotations ~{cosmic_vcf} \
            --columns "INFO/COSMIC_ID,INFO/TISSUE,INFO/TUMOR,INFO/FATHMM,INFO/SOMATIC" \
            --output ~{base_name}.annot_cosmic.tmp.vcf \
            ~{input_vcf}


            #must groom for gatk compat
            ~{scripts_path}/groom_vcf.py \
                ~{base_name}.annot_cosmic.tmp.vcf ~{base_name}.annot_cosmic.vcf
      
            bgzip ~{base_name}.annot_cosmic.vcf
      
            tabix ~{base_name}.annot_cosmic.vcf.gz

    >>>

    output {
        File vcf = "~{base_name}.annot_cosmic.vcf.gz"
        File vcf_index = "~{base_name}.annot_cosmic.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}    




task open_cravat {

    input {
        File input_vcf
        File input_vcf_index

        # must specify cravat_lib_dir or cravat_lib_tar_gz
        File? cravat_lib_tar_gz #providing the tar.gz file with the cravat resources
        String? cravat_lib_dir  #path to existing cravat lib dir in the ctat genome lib
        String genome_version

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu = 1
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50 + if(defined(cravat_lib_tar_gz))then 100 else 0)
        

    }


    command <<<

        set -ex

        echo "########### Annotate CRAVAT #############"

        cravat_lib_dir="~{cravat_lib_dir}"
        # cravat
        if [ "$cravat_lib_dir" == "" ]; then
            
            if [ "~{cravat_lib_tar_gz}" == "" ]; then
                 echo "Error, must specify cravat_lib_tar_gz or cravat_lib path"
                 exit 1
            fi
            
            #use the provided tar.gz cravat lib      
            cravat_lib_dir="~{cravat_lib_tar_gz}"

            mkdir cravat_lib_dir
            compress="pigz"

            if [[ $cravat_lib_dir == *.bz2 ]] ; then
                compress="pbzip2"
            fi

            tar -I $compress -xf $cravat_lib_dir -C cravat_lib_dir --strip-components 1
            cravat_lib_dir="cravat_lib_dir"

        fi
        
        export TMPDIR=/tmp # https://github.com/broadinstitute/cromwell/issues/3647

        ~{scripts_path}/annotate_with_cravat.py \
            --input_vcf ~{input_vcf} \
            --genome ~{genome_version} \
            --cravat_lib_dir $cravat_lib_dir \
            --output_vcf ~{base_name}.cravat.tmp.vcf

            #must groom for gatk compat
            ~{scripts_path}/groom_vcf.py \
                ~{base_name}.cravat.tmp.vcf ~{base_name}.cravat.groom.vcf

            bcftools sort -T . ~{base_name}.cravat.groom.vcf > ~{base_name}.cravat.vcf
            bgzip -c ~{base_name}.cravat.vcf > ~{base_name}.cravat.vcf.gz
            tabix ~{base_name}.cravat.vcf.gz

        

    >>>

    output {
        File vcf = "~{base_name}.cravat.vcf.gz"
        File vcf_index = "~{base_name}.cravat.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}    


task rename_vcf {
    input {
        File input_vcf
        File input_vcf_index
        String base_name

        String docker        
        Int preemptible
        Int cpu = 1
        Int disk = ceil((size(input_vcf, "GB") * 2)) 
    }


    command <<<

        set -ex

        echo "####### Final step: Renaming Vcf ########"
      
        mv ~{input_vcf} ~{base_name}.vcf.gz
        mv ~{input_vcf}.tbi ~{base_name}.vcf.gz.tbi

     >>>



    output {
        File vcf = "~{base_name}.vcf.gz"
        File vcf_index = "~{base_name}.vcf.gz.tbi"
    }
	
    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}

