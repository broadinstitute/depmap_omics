# Report 21Q3

## RNA

switched ACH-001096 - ACH-001196
dropped {'ACH-001741', 'ACH-001163', 'ACH-001316', 'ACH-000010', 'ACH-001173'}
ACH-001709 WES
## WGS

1. CGA pipeline does not call cache
2. can't merge with current CN segment merging tool (need to better parse file path (going from the end))


## TODO:

- check why did not propagate from renaming json in mutation to fusions
- drop 001709 in mutations
- check if still add 001709 and 000195 even if no data (Sanger mutation only)
-

# decisions

## took already

removed blacklist status: (?)
CDS-VnMBYD (ACH-000164) wgs 
CDS-Rd4nMx (ACH-000454) wgs
CDS-dgxjAa (ACH-000550) wgs
CDS-TGTiB8 (ACH-000739) wgs
CDS-YSRYLi (ACH-000928) wgs
CDS-ZJh6UN (ACH-001239) wgs
all (ACH-000632)

changed version: (V)
ACH-000015 CDS-DjghjI (WES Sanger) v1 to v2 and CDS-DAiFZk v2 to v1
ACH-001172 CDS-o86jP7 (WES) v2 to v1 and CDS-rQjmI6 (WES) v1 to v2
ACH-001172 CDS-FFtZ3f (RNA) v2 to v1 and CDS-Y7UJg0 (RNA) v1 to v2
ACH-000831  CDS-pEGMSA (Sanger) v1 to v2 and CDS-aLoqXc (Broad) v2 to v1
ACH-000833 CDS-utyylb (rna): v1 -> v2 and CDS-pvNwpR (rna) : v2 -> v3 and CDS-tTwJLo (rna) : v3 -> v1
ACH-001164 CDS-Fauu92 (rna) v2 to v1 and CDS-V5yMiU (rna) v1 to v2

failed QC (WGS) (?)
'CDS-AqZLna' → ACH-002512
'CDS-dIijHP' → ACH-000309
'CDS-ovDo8V' → ACH-000890
--> decision to keep

(RNA) (?)
'CDS-VvIG7p'
--> decision to keep

keeping failed QC chordoma lines (?)
'CDS-Rl87Z1': 'ACH-001956',
'CDS-mys9Dm': 'ACH-001955',
'CDS-TzQAjG': 'ACH-001957'
and broad WES
'ACH-001709': 'CDS-TuKZau'


Removed ​​the following from the tracker: (N)
all (ACH-001078)
all (ACH-002184)
all (ACH-001146)
all (ACH-002022)
all (ACH-002707)
all (ACH-001790)
CDS-P4ZoH4 (ACH-002184)

changed arxspan ids (?)
(ISSUE: done without changing the other cell line specific informations) 
changed CDS-cXQOkP to ACH-001543
Changed CDS-uNAcAR to ACH-000201, (and changed CDS-eaMPk2 (ACH-000201) to version to 2 )

removed from legacy: (V)
{'ACH-001096', 'ACH-001196', 'ACH-000961', 'ACH-000511', 'ACH-000375', 'ACH-000278', 'ACH-001709', 'ACH-002475', 'ACH-001063', 'ACH-000641', 'ACH-000090', 'ACH-000658','ACH-001146', 'ACH-002260', 'ACH-000010', 'ACH-001078'}


## NOT TOOK yet
stop blacklisting:
ACH-000901


## on the code (already shown above)

additional blacklisted (should already be blacklisted in the sample tracker but they got removed from the tracker so caused an issue in the pipeline):

CN (shoulld be blacklisted)  todropwes+=['CDS-GNOJc5', 'CDS-P4ZoH4', 'CDS-ywnbJT']

todropwgs+=['CDS-ymRIxH']

to remove from wes legacy

WES_legacy_blacklist = ['ACH-001096', 'ACH-001196', 'ACH-000961', 'ACH-000511',
                          'ACH-000375', 'ACH-000278', 'ACH-001709', 'ACH-002475',
                          'ACH-001063', 'ACH-000641', 'ACH-000090', 'ACH-000658']
