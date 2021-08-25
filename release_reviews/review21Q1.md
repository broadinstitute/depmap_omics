Based on this message I want to remove :

https://broadinstitute.slack.com/archives/C01B2CLLQEQ/p1604509300225800?thread_ts=1604502738.219300&cid=C01B2CLLQEQ

'ACH-000561',
 'ACH-001131',
 'ACH-001189',
 'ACH-002217',
 'ACH-002315',
 'ACH-002341',
 'ACH-002390',
 'ACH-002391',
 'ACH-002393',
 'ACH-002394',
 'ACH-002395',
 'ACH-002396'

## 25/12
- Mutations
  - sample failed QC, 
  - need to add wgs and legacy
  - waiting for the WGS data. 
  - Need to get results from CN pipeline (QC, duplicates, misslabelling)
  - and change them in our workspaces). 
  - Need to check for additional issues in legacy datasets. 
  - reupload them afterwarrd
- RNA
  - finished. things to discuss.
    - some new genes, some genes got rermoved. 
    - I don't really know why. likely because of filtering... 
    - should we do filtering?
- CN
  - Sample failed on Jeremie's standard QC. 
    - Running chekc for misslabelling (similarly to RNAsq, based on correlation with previous releases). 
    - Then need to modify them everywhere and in mutations as well
  - remapping, need to check plots with new PON, 
    - then launch with a smoothing of 5
  - Need to check samples on legacy (misslabellings)
  - upddate workspace and tracker with cn data

## 4/01

- Finished a bunch of stuff during the holidays
- Found some issues we might want to discuss
- I will never be able to write down every changes in the pipeline but we have the log files
- we really can be pretty happy with quality of WGS CN
- 

- mutations
  - wgs: complex line kept failing (not enough compute) <-- waiting for it
  - wes: looking for issues (misslabelling, not found, failed qc lines in legacies)
  - expected to finish by today
- cn:
  - need to make sure of the quality of the merge file (final qc)
  - expected finish by today
- fusion:
  - done 
- expression:
  - genes added, genes removed... (genes renaming got updated? genes set to 0 when all var/mean==0)
  - otherwise: done
- upload to taiga:
  - need to make the readme
  - need to test the new upload function (based on feedback from last quarter)
  - expected finish by tomorrow



### not found 
ACH-002217 Sanger
ACH-002335 Chordoma
ACH-002378 Sanger

### failed qc 
ACH-001956, ACH-001955', 'ACH-001957' (chordoma)

### wrong new
ACH-001011
ACH-001108
ACH-001187

### wrong (from 20Q4)
ACH-001189
ACH-002303
ACH-002315
ACH-001675

### found version for:
'CDS-b5ElTm':"ACH-000157", 
"CDS-up4Vo5":"ACH-000662",
"CDS-CWA37D": "ACH-000825", 
"CDS-CCAK2f":"ACH-001328",
"CDS-2jBQ8n":"ACH-000757",
"CDS-T8W6P4": "ACH-000398",
"CDS-9TDVpH": "ACH-000685",
"CDS-dQKiht": "ACH-000375",
"CDS-Ckptje": "ACH-002291",
"CDS-ljFuDX": "ACH-001339",
"CDS-5x4qLj": "ACH-000608",
"CDS-UxJcOY": "ACH-000561",
"CDS-TUYedU": "ACH-000261",
"CDS-RLVrVE": "ACH-001523",
"CDS-6liik0": "ACH-000561",
"CDS-b5ElTm": "ACH-000157",
"CDS-u9hZ60": "ACH-000077",
"CDS-NUlX3d": "ACH-000458",
"CDS-2HO10g": "ACH-000278"

### duplicates:
RI-1 RI1
{'ACH-000261', 'ACH-000398'}
REH Reh
{'ACH-000960', 'ACH-000473'}
RERF-LC-AI RERFLCAI
{'ACH-000261', 'ACH-000960'}
Hs 852.T HS852T
{'ACH-000274', 'ACH-001523'}
RT112 RT-112
{'ACH-000398', 'ACH-000473'}
CH157MN CH-157MN
{'ACH-000511', 'ACH-000025'}
CORL47 COR-L47
{'ACH-000695', 'ACH-000662'}
CALU1 Calu-1
{'ACH-000608', 'ACH-000511'}
COR-L23 CORL23
{'ACH-001339', 'ACH-000662'}

## R

mutations
we donot care about this 53rd line


expression
check duplicates
add in readmes http://useast.ensembl.org/Homo_sapiens/Info/Annotation

taiga:
write down readmes
list of things:



## expression:
- updated list of genes. reminder, updates now come from updates in ensembl's biomart: http://useast.ensembl.org/Homo_sapiens/Info/Annotation
- removed duplicate genes
- new better ssGSEA geneset file now look for: now just using the "All gene set" file from msigDB http://www.gsea-msigdb.org/gsea/downloads.jsp.
- renaming files: 
  - expression_transcripts_tpm: CCLE_RNAseq_transcripts
  - expression_genes_tpm: CCLE_expression_full
  - expression_genes_expected_count: CCLE_RNAseq_reads
  - expression_proteincoding_genes_tpm: CCLE_expression

## fusion:
- updated the fusion pipeline to starfusion:1.7.0

## mutations:
- renaming files:
  - all_somatic_mutations_withlegacy: CCLE_mutations

## copy number:
- stopped reporting X chromosome amplification status as it is biased (we are using a male PON)
- improved and solved reporting of WGS CN data (solved issues related to gene copy number set to 0, solved logtransform applied twice to some lines, solved X chromosome being reported with a Copy ratio of 2 when from female samples (now male have copy ratio of .5 and female of 1. But for Y chromosome, males have a copy ratio of 1))
- removed 'chr' chromosome prefix in seg CN file
- re-extended some segments
- resolved wrong source annotation in segment file
- removed duplicate lines
- renaming files:
  - all_merged_genes_cn: CCLE_gene_cn
  - all_merged_segments: CCLE_segment_cn

## CN and mutations
- removed lines that should have been removed in 20Q4 ACH-001189, ACH-002303, ACH-002315, ACH-002341
- removed a bunch of normal lines: 'ACH-000043', 'ACH-000063', 'ACH-000071', 'ACH-000079', 'ACH-000083', 'ACH-000088', 'ACH-000119', 'ACH-000125', 'ACH-000131', 'ACH-000133', 'ACH-000134', 'ACH-000135', 'ACH-000154', 'ACH-000165', 'ACH-000175', 'ACH-000180', 'ACH-000184', 'ACH-000185', 'ACH-000194', 'ACH-000199', 'ACH-000214', 'ACH-000224', 'ACH-000229', 'ACH-000230', 'ACH-000240', 'ACH-000275', 'ACH-000284', 'ACH-000306', 'ACH-000340', 'ACH-000413', 'ACH-000526', 'ACH-000529', 'ACH-000531', 'ACH-000539', 'ACH-000540', 'ACH-000742', 'ACH-000797', 'ACH-000850', 'ACH-001093', 'ACH-001207', 'ACH-001767', 'ACH-002342', 'ACH-002343', 'ACH-002344', 'ACH-002346', 'ACH-002347', 'ACH-002350', 'ACH-002358', 'ACH-002359', 'ACH-002360', 'ACH-002361', 'ACH-002362', 'ACH-002363', 'ACH-002364', 'ACH-002365', 'ACH-002366', 'ACH-002367', 'ACH-002368', 'ACH-002369', 'ACH-002370', 'ACH-002371', 'ACH-002372', 'ACH-002373', 'ACH-002374', 'ACH-002375', 'ACH-002376', 'ACH-002377', 'ACH-002379', 'ACH-002380', 'ACH-002381', 'ACH-002382', 'ACH-002383', 'ACH-002384', 'ACH-003210'
- removing 3 duplicates linked to crispr engineered lines in Achilles: ACH-003000, ACH-002875, ACH-002874

## TODO:

- rerun the CN with the new set
- rerun the mutations with the new set
- check correlation with rnaseq of failed qc lines
