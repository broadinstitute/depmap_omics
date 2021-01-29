# 20Q4 post mortem

- [20Q4 post mortem](#20q4-post-mortem)
  - [Recap:](#recap)
  - [Issues:](#issues)
  - [Solutions](#solutions)
    - [Create a review process for updates to the pipelines.](#create-a-review-process-for-updates-to-the-pipelines)
    - [Changes in the way we code](#changes-in-the-way-we-code)
    - [Have direct communications.](#have-direct-communications)
  - [Current CCLE Issues:](#current-ccle-issues)
  - [TODO:](#todo)
  - [postponed changes](#postponed-changes)
  - [current changes:](#current-changes)
    - [Made a document explainnig the pipelines.](#made-a-document-explainnig-the-pipelines)
    - [Made a document explainnig time spent on things.](#made-a-document-explainnig-time-spent-on-things)
  - [Time management](#time-management)
    - [data loading:](#data-loading)
      - [+](#ulliliul)
      - [-](#ulliliul-1)
    - [Terra processing](#terra-processing)
      - [+](#ulliliul-2)
      - [-](#ulliliul-3)
    - [Post processing](#post-processing)
      - [+](#ulliliul-4)
      - [-](#ulliliul-5)
    - [Taiga Upload](#taiga-upload)
      - [+](#ulliliul-6)
      - [-](#ulliliul-7)
    - [Conclusion](#conclusion)

## Recap:

In 19Q3 I was given a set of very approximative code and guidelines and a bunch of issues to solve.
I had some supervision by Neekesh related mostly to the biology and some ideas of why things were the way they were. There was no one to explain how things were done before and why we are doing them this way, except Allie for some of the pipeline.
So I made mistakes for the first 2 releases. Mostly regarding how I understood the data should be delivered and in what format.
We had also many people coming to us with issues regarding the data.


I built a pipeline from this and tried to make it so that it could run while solving blocking issues that it had.
Out of all of this we laid out a plan with Neekesh, in what to do to solve them and what to do next for ccle.
This plan was set as series of tasks to be done in Asana.
I explained the plan to Javad (maybe too early) and worked on it since 19Q4.
A lot of the changes were not seen until 20Q2 as we postponed them until we couldn't.


I am working on the pipeline only during the release: 3-4 week process. I have only half my time on ccle so 6-8 weeks which is 2 months. I have a third of my time in doing things related to science and research: 2-3 weeks / Quarters.
When I was hired I was promised by Neekesh that the pipeline would not be more than a week every 3 months.
Now I end up doing 2/3 of my time at the Broad doing soft engineering and being very stressed by deadlines etc.

-> Not what I came to the US for.

## Issues:

- We have had 2 issues during this release:
  - Some samples were added from 2 Quarters ago.
  - Discrepancy between the way we checked samples
  - Issues with managing the changelog changes (too many got added / removed...)
  - Some columns/files got added to the releases and should not have been.
  - Modifications in log-transform/filenames

-> This has led to some distrust.

- At the end this release I have seen people taking notes for me.
- Going over the details of everything I have to do.
- Wanting to redo some of my analysis.

I don't think that is useful to do and it overlooks the real issues.
Hear me here: I am extremely convinced that anyone else would face at least as much issues given the state of things.

- The first one is a misunderstanding from what happened before and highlights the many good changes that undertakes but might lead to confusion between what happened before and what happens now.
- The second ones highlight the issues in communication between our teams. Not helped with covid.
- The other ones also show the state of the pipeline. it is yet hard to make different version of the data

## Solutions

### Create a review process for updates to the pipelines.

- everything is always on the Asana task. and changes are well in advance in the sub-project.
- We should have a meeting a couple weeks before the start of the release to discuss some non postponable changes.
- There is a latest dataset on taiga. the idea was for all non required changes to go there and be used as the latest version of the pipeline (tried on by achilles and the portal team) 

### Changes in the way we code
  
- Create two version of the pipeline: dev/prod (soft eng). it is so that we can always manage the updates and keep the previous working version. (yes it should have been done sooner. but I started with a file on google drive)
- Portal team and ccle team need to work on the same sample subletting. we can't have two versions. For code you don't create two function that do the same things. you create one that works. Most of the issues we found were because oof pipelines doing things differently. We need to work on the same file/set of function to do this. Not have each our own code and not talking to each other nor understanding how it works.

### Have direct communications.

- Becky can't always be in between the portal team, the CCLE team, Achilles team. some information will get lost and it would be more efficient to have meetings between us.
- I need to be present with you Javad when you present changes to the pipeline. you told me that at many occasions you could not explain the detail of things. that is normal. I need to be present to explain why things cannot be done differently and understand why the other group need things the way they are. The more intermediates the more chances for confusion.

## Current CCLE Issues:

- Binary matrices are bad
- Star-fusion out of date (even more now than before). not the same data for all of them
- Many ,way-too-different copy-number data.
- Around a dozen duplicates to be removed
- Around a dozen samples are likely mislabelled
- Portal hotspot needs to be removed and showcased as stars.
- Never know what to add in the log file. who is it for. is there an actual changelog?
- LogFoldChange just some data. users get confused
- No easy way to remove files from the terra pipelines. storage costs getting pretty high.
- Duplicates: 

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

- normals ACH-002375 ACH-002361 ACH-002383 ACH-002373 ACH-002382 ACH-002374 ACH-002380 ACH-002360 ACH-002370 ACH-002381 ACH-002377 ACH-002372 ACH-002378 ACH-002376 ACH-002363 ACH-002369 ACH-002365 ACH-002379 ACH-002362 ACH-002364 ACH-002384 ACH-002266 ACH-002371 ACH-002367 ACH-002366 ACH-002368 ACH-002343 ACH-002154 ACH-002150 ACH-002265 ACH-002153 ACH-002151 ACH-002152 ACH-002346 ACH-002347 ACH-002351 ACH-002344 ACH-002350 ACH-002359 ACH-002342

## TODO:

- Starfusion/star needs to be updated
- Create a wdl update script
- Create 2 erase script for files generated by our pipelines (when a job fails & unused files)
- make a dev and prod version of the pipeline.

## postponed changes

- Renaming most files to better acknowledge the new files that we will have either this release or the next one:
  - all_merged_genes_cn: CCLE_gene_cn
  - all_merged_segments: CCLE_segment_cn
  - all_somatic_mutations_withlegacy: CCLE_mutations
  - expression_transcripts_tpm: CCLE_RNAseq_transcripts
  - expression_genes_tpm: CCLE_expression_full
  - expression_genes_expected_count: CCLE_RNAseq_reads
  - expression_proteincoding_genes_tpm: CCLE_expression
- Removed log2+1 transform for all datasets that had it, so that it is consitent everywhere and does not lead to misinterpretations
- Unmerged datasets (CN WGS, Mutations non, legacy only)
- Full CN reprocessing
- New binary matrix release
- Binary matrix for 
  - all_somatic_mutations_boolmatrix_fordepmap_damaging
  - all_somatic_mutations_boolmatrix_fordepmap_othercons
  - all_somatic_mutations_boolmatrix_fordepmap_othernoncons
  - all_somatic_mutations_boolmatrix_fordepmap_hotspot
- Additional columns in mutations
- Unfiltered mutations
- GSEA file
- Rationalize READMEs (will now be in the datasets on taiga as well, will stop using taiga ReadMes)
- Make sure that ACH-003000, ACH-002875, ACH-002874 are never released

## current changes:

- Latest dataset: If everything goes well, what is in the latest dataset will be in the next dataset on the next release. It should be used to prepare your pipelines for the new version. Every changes made are explained in the readmes on taiga.
- Unfiltered mutations (the raw somatic mutations calls from our pipeline. no quality filtering.
- Rerun achilles on new CN dataset
- Rerun TDA on new mutation matrices and new CN
- Try to release some more files in the portals (unfiltered mutations for ex?) 
- Issue with mutation annotations in the portal

### Made a document explainnig the pipelines.

### Made a document explainnig time spent on things.

## Time management

### data loading:

_1-4 days_

#### +
+ if any changes in the release process
+ if more samples

#### -
- if more sm-id and check on cellosaurus in advance _(GP's decision, Ops Fix)_
- if check for duplicates in advances _(GP's decision, Ops Fix)_

Best: 20mns

### Terra processing

_3-7 days_

#### +
+ change in the bam files
+ adding/updating pipelines
+ change in who get the samples from
+ unknown changes in Terra
+ unknown changes in groups we use pipelines from

#### -
- productionalization (_everything set up in asana. 5-6 months of work minimum. requires a computationalist_)

Best: 
1 day

### Post processing

_2-7 days_

#### +
+ new data generated
+ issues to fix
+ quality of the samples

#### -
- finish auto QCs (CN/Mutations) _(are we all ok with automated QCs, 1-2 months of work, some thinking to do)_
- do pre-releases (_not enough bandwidth_)

Best: 
1 day

### Taiga Upload

_1-2 days_

#### +
+ changes in the release process
+ changes in the subsetting process
+ new datasets

#### -
- do the subset on the portal's side (_to decide_)

Best:
2 hours

### Conclusion

Most important factor: changing things.
Most relevant todo: get productionalization, work on the input and output processes
If everything is decided: possible to run the pipeline without much supervisions, with more control from OPs' side, all in ~2days
