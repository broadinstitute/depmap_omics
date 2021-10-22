Based on this message I want to remove :

https://broadinstitute.slack.com/archives/C01B2CLLQEQ/p1604509300225800?thread_ts=1604502738.219300&cid=C01B2CLLQEQ

new
ACH-001011
ACH-001108
ACH-001187
already decided
ACH-001189
ACH-002303
ACH-002315
ACH-002341

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
