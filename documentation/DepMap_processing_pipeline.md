![schema](documentation/architecture_diagram_white.png)

We are using a set of key tools to process the sequencing output:
- __star (from docker image `gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10`)__:
  - [https://www.ncbi.nlm.nih.gov/pubmed/23104886](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
  - [https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
  - aligns RNAseq bam files for downstream processing
- __rsem (from docker image `gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10`)__: 
  - [https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)
  - quantifes gene and isoform abundances from RNAseq data
- __star fusion (from docker image `trinityctat/starfusion:1.7.0`)__: 
  - [https://github.com/STAR-Fusion/STAR-Fusion/wiki](https://github.com/STAR-Fusion/STAR-Fusion/wiki)
  - [http://biorxiv.org/content/early/2017/03/24/120295](http://biorxiv.org/content/early/2017/03/24/120295)
  - outputs fusion predictions from RNAseq data
- __mutect__: 
  - [https://software.broadinstitute.org/cancer/cga/mutect](https://software.broadinstitute.org/cancer/cga/mutect)
  - [https://youtu.be/rN-cLrb5aGs](https://youtu.be/rN-cLrb5aGs)
  - [https://www.nature.com/articles/nbt.2514](https://www.nature.com/articles/nbt.2514)
- __mutect2 (from `broadinstitute/gatk:4.2.6.0`)__:
  - [https://gatk.broadinstitute.org/hc/en-us/articles/360036490432-Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360036490432-Mutect2)
  - outputs mutation calls from RNAseq/WES/WGS data
- __gatk cnv__:
  - [https://software.broadinstitute.org/gatk/documentation/article?id=11682](https://software.broadinstitute.org/gatk/documentation/article?id=11682)
  - outputs relative segment and copy number from WES/WGS data
- __strelka__:
  - [https://www.nature.com/articles/s41592-018-0051-x](https://www.nature.com/articles/s41592-018-0051-x)
  - [https://github.com/Illumina/strelka](https://github.com/Illumina/strelka)
- __PureCN__:
  - [https://github.com/lima1/PureCN](https://github.com/lima1/PureCN)
  - computes absolute copy number, as well as features including loss of heterozygosity (LOH), LOH fraction, ploidy estimate, Whole Genome Doubling (WGD), Chromasomal Instability (CIN) from WES/WGS data
- __MSIsensor2__:
  - [https://github.com/niu-lab/msisensor2](https://github.com/niu-lab/msisensor2)
  - computes Microsatellite Instability (MSI) score from WES/WGS data
- __Manta__:
  - [https://github.com/Illumina/manta](https://github.com/Illumina/manta)
  - calls structural variants from WES/WGS data

The following flowchart provides another good overview of what the pipeline is doing.

![](documentation/updated-flowchart.png)

Note that the input references and parameters used in generating DepMap data for the following pipelines in any given quarter can be found in `data/*quarter*/`.

#### Copy Numbers and Somatic Mutations

We are running the following workflows in this order to generate copy number and mutation datasets:

[WGS_pipeline](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/WGS_pipeline:master?tab=info) runs several sub-processes to generate relative and absolute copy number segments, mutation MAF data, structural variant calls, and various genomic features including loss of heterozygosity (LOH), LOH fraction, ploidy estimate, Whole Genome Doubling (WGD), Chromasomal Instability (CIN), and MSI score.

[WGS_aggregate](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/WGS_aggregate:master?tab=info) aggregates CN segments and MAFs into their respective files.

The outputs to be downloaded will be saved under the sample set that you ran. The outputs we use for the release are:

*   combined_seg_file
*   filtered_CGA_MAF_aggregated
*   unfiltered_CGA_MAF_aggregated

* Note that additional files are used for QC

#### Expression and Fusion

We are generating both expression and fusion datasets with RNAseq data. Specifically, we use the [GTEx pipeline](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md) to generate the expression dataset, and [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) to generate gene fusion calls. This task also contains a flag that lets you specify if you want to delete the intermediates (fastqs) that can be large and might cost a lot to store. Run the following tasks on all samples that you need, in this order:

[RNA_pipeline](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/RNA_pipeline:master?tab=info) imports and runs several sub-processes to generate RNA expression and fusion data matrices.

[RNA_aggregate](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/RNA_aggregate:master?tab=info) aggregates expression and fusion data files into their respective aggregated file.

The outputs to be downloaded will be saved under the sample set that you ran. The outputs we use for the release are:

*   rsem_genes_expected_count
*   rsem_genes_tpm
*   rsem_transcripts_tpm
*   fusions_star

__Finally, we save the workflow configurations used in the pipeline runs__

**Remarks:**
- for the copy number pipeline we have parametrized both an XX version and an XY version, we recommend using the XY version as it covers the entire genome
- for the mutation pipeline we are working on Tumor-Normal pairs which explain some of the back and forth between the two workspace data table. (workflows works as well with or without matched normals.)
- for the expression pipeline, we have an additional set of workflows to call mutations from RNAseq, this might not be relevant to your need.

### 3. Downloading and Postprocessing (often called **on local** in the notebooks) <a name="downloading-postprocessing"></a>

This step will do a set of tasks:
- clean the workspaces by deleting large useless files, including unmapped bams.
- retrieve from the workspace interesting QC results.
- copy realigned bam files to our own data storage bucket.
- download the results.
- remove all duplicate samples from our downloaded file (keeping only the latest version of each sample).

_You would only be interested here at minima in the result downloading_
 
...and postprocessing tasks. The main postprocessing steps for each pipeline are as followed:

#### CN

`copynumbers.py` contains the main postprocessing function (`postProcess()` and their wrappers for internal use) responsible for postprocessing segments and creating gene-level CN files.

#### Mutation

`mutations.py` contains `postProcess()` (and its wrappers for internal use, including one for filtered and one for unfiltered mutation data), which is responsible for postprocessing aggregated MAF files and generating various types of mutation datasets.

#### RNA

`expressions.py` contains the main postprocessing function (`postProcess()` and their wrappers for internal use) responsible for postprocessing aggregated expression data from RSEM, which removes duplicates, renames genes, filters and log transforms entries, and generates protein-level expression data files.

##### Fusion

Functions responsible for postprocessing aggregated fusion data can be found in `fusions.py`. We want to apply filters to the fusion table to reduce the number of artifacts in the dataset. Specifically, we filter the following:

* Remove fusions involving mitochondrial chromosomes, or HLA genes, or immunoglobulin genes
* Remove red herring fusions (from STAR-Fusion annotations column)
* Remove recurrent in CCLE (>= 25 samples)
* Remove fusion with (SpliceType=" INCL_NON_REF_SPLICE" and LargeAnchorSupport="No" and FFPM < 0.1)
* Remove fusions with FFPM < 0.05 (STAR-Fusion suggests using 0.1, but looking at the translocation data, this looks like it might be too aggressive)


**Remarks:**
- in the RNAseq pipeline we have an additional sub-pipeline at the end of the notebook to process the fusion calls from starFusion
- to get the exact same results as in CCLE, be sure to run `genecn = genecn.apply(lambda x: np.log2(1+x))` to the genecn dataframe in the CNV pipeline (present at the end of the validation steps).
- we do not yet have integrated our germline calling in the mutation pipeline but you can still find the haplotypeCaller\|DeepVariant workflows and their parameters


### 4. QC, Grouping and Uploading to the Portal (internal use only) <a name="qc-grouping-uploading"></a>

We then perform the following QC tasks for each dataset. These tasks should not be very interesting for external user as they revolve around manual checks of the data, comparison to previous releases, etc.

#### CN

Once the CN files are saved, we load them back in python and do some validations, in brief:

- mean, max, var...
- to previous version: same mean, max, var...
- checkAmountOfSegments: flag any samples with a very high number of segments
- checkGeneChangeAccrossAll: flag any genes which stay at a similar value across all samples

#### Mutation

__Compare to previous release (broad only)__

We compare the results to the previous releases MAF. Namely:

- Count the total number of mutations per cell line, split by type (SNP, INS, DEL)
- Count the total number of mutations observed by position (group by chromosome, start position, end position and count the number of mutations)
- Look at specific differences between the two MAFs (join on DepMap_ID, Chromosome, Start position, End position, Variant_Type). This is done for WES only


#### RNA

Once the expression files are saved, we do the following validations:
- mean, max, var...
- to previous version: same mean, max, var...
- we QC on the amount of genes with 0 counts for each samples


> After QC, we are also preparing the data to be released to different groups, removing the samples per access category: Blacklist\|Internal\|DepMapConsortium\|Public.

We are then uploading the data to a server called taiga where it will be used in the depmap portal