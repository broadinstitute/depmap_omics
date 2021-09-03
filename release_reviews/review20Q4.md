# review of the 20Q4 release process:

1. managing all the different tasks and seting up WGS, getting access to Terra
2. loading the data for WES-WGS, setting up the task with the group, finding issues. 
3. creating functions to retrieve data that had no arxspan ids, rewritinig the loading function for WGS lines
4. finishing relabelling, solving some issues with ways bam date was computed for some bams, solving workspaces which had bad naming. creating pairs, sets, rewriring the pipelines for samples from ccle. Running the first workflows
5. Issues with BAM contigs. Now continuing running the pipelines, going to RNAseq in the meantime 
6. debugging small issues during the week end and continuing processes.
7. - look at other data from CNV pipeline (why are we using this output; there is so many others that seem good).
	 - Change CNN2D into CNN1D and then rerun CNN with CNN2D parameter
	 -  
8. issues with most aggregation tasks on WGS. rewriting part of the WGS pipeline for getting the CN data and removing failed data, uploading remmaped version, getting QC. 
	- getting older version of the data, finding issues with some lines
	- making a legacy dataset
	- re-owing the post processing and QCing and comparison to previous versions
9. debugging issues with gene level CN compute.
	- comparison of old to new, XY to XX, versions to current version
	- making the different datasets
	- separating the achilles tasks from ccle's
	- investigating issues with correlations
10: 
	- solving issues related to wrong versions



## NOW:

### Objectives
- WGS: release everything. 
- WES: on my version of averaging only if correlation if above 0.8.
-> if we have >5% that are below 0.75 replace them
- set as latest version the most highly correlated of the two.
- release new version only.

### tasks CN
- create new taiga folders:
- put all on taiga
- do the same for WGS
- add WGS to Achilles version 

### tasks mut

- separate legacy mutations from the rest
- make mutation postp in python
- rerun MANTA in WGS mode
- get merge of candidates SV and indels from manta for unfiltered SV data


- debug fusion docker image (see wdl)


##For later

- add funcotator to CN data
- add VEP & funcotator to SVs and indels
- make a google spreadsheet merging tool
- try to use CNN_2D for germline mutations

## For ASAP
- remove the toremove set in the WES workspaces
oldtoremove = ['CDS-jqOvtj',
 'CDS-FRxdcH',
 'CDS-6Yy3Yj',
 'CDS-CuJ0f8',
 'CDS-rLRUbG',
 'CDS-PdUZxY',
 'CDS-eUqT7L',
 'CDS-KbbgMb',
 'CDS-6da3hu',
 'CDS-fXMRF9',
 'CDS-CMenCH',
 'CDS-MLJbT2',
 'CDS-QVhVDT',
 'CDS-XevQNc',
 'CDS-0pZb0j',
 'CDS-6l3V79',
 'CDS-MnF3x8',
 'CDS-ihI7Dp',
 'CDS-34hKv3',
 'CDS-TyWjJs',
 'CDS-4sr6RL',
 'CDS-M8xDMS',
 'CDS-TpDBjm',
 'CDS-W80jkV',
 'CDS-agZcmk',
 'CDS-cYWYp7',
 'CDS-IJnjkY',
 'CDS-0aJ4Yh',
 'CDS-txTRwz',
 'CDS-gIMBax',
 'CDS-1p2nnc',
 'CDS-KQDgIV',
 'CDS-Eq9UNX',
 'CDS-3M6Pq9',
 'CDS-qZsCuJ',
 'CDS-0lfqVz',
 'CDS-o4dXGr',
 'CDS-uQ8nnX',
 'CDS-iqPqOr',
 'CDS-Dkl8OF',
 'CDS-Hj3xAa',
 'CDS-3WygAj',
 'CDS-oHu1Ik',
 'CDS-X3c4UY',
 'CDS-PYw8ID',
 'CDS-Sp18uD',
 'CDS-leGxSD',
 'CDS-SJq3p4',
 'CDS-no7ysz',
 'CDS-UnDaBI',
 'CDS-eowEZF',
 'CDS-HNytLD',
 'CDS-KYkMDa',
 'CDS-OgPf0h',
 'CDS-OCkOqy',
 'CDS-QU7ftt',
 'CDS-iEULQm',
 'CDS-ODmXrP',
 'CDS-YMIv9D',
 'CDS-5rD8XC',
 'CDS-QXBhht',
 'CDS-9XPgHB',
 'CDS-Ig6N9S',
 'CDS-UtrDTK',
 'CDS-nby0QM',
 'CDS-49azaP',
 'CDS-9qDPiX',
 'CDS-KgRznV',
 'CDS-picEuX',
 'CDS-L0pDPl',
 'CDS-kxNZ5S',
 'CDS-1djAlo',
 'CDS-YYLKZ0',
 'CDS-pXMN9C',
 'CDS-gRA4SM',
 'CDS-QHp4h4',
 'CDS-B0qAaq',
 'CDS-1b1Hxk',
 'CDS-5wYxZS',
 'CDS-cyuMYb',
 'CDS-XQkXf4',
 'CDS-7PFldq',
'CDS-3EBt51',
 'CDS-UV1pVE',
 'CDS-WedVJA',
 'CDS-WfjTcJ',
 'CDS-bntBUl',
 'CDS-cAEii6',
 'CDS-d18Xie',
 'CDS-dpub1O',
 'CDS-yPSmxb']

- add Baits, origin, parent cell line info in our release pipelines
- clean workspaces better

- ask back/find back deleted WES 
'CDS-519oCd','gs://cclebams/wes/C836.JHOS-4.1.bam',
 'CDS-5LolJ4', 'gs://cclebams/wes/C836.MJ.2.bam',
 'CDS-8sQWae', 'gs://cclebams/wes/C836.OVMANA.2.bam',
 'CDS-AJMYsd', 'gs://cclebams/wes/C836.Panc_02.13.2.bam',
 'CDS-EzZEgz', 'gs://cclebams/wes/C836.SNU-738.2.bam',
 'CDS-NE61ys', 'gs://cclebams/wes/C836.TO_175.T.1.bam',
 'CDS-Sl300T', 'gs://cclebams/wes/C836.YD-8.1.bam',
 'CDS-T10Uph', 'gs://cclebams/wes/CDS-8sQWae.bam',
 'CDS-VS9XDY', 'gs://cclebams/wes/CDS-AJMYsd.bam',
 'CDS-cmV75B', 'gs://cclebams/wes/CDS-EzZEgz.bam',
 'CDS-d2xe4x', 'gs://cclebams/wes/CDS-T10Uph.bam',
 'CDS-dGaHMd', 'gs://cclebams/wes/CDS-VS9XDY.bam',
 'CDS-ezLdbO', 'gs://cclebams/wes/CDS-cmV75B.bam',
 'CDS-fb9VZf', 'gs://cclebams/wes/SANGER_C32_SKIN.bam',
 'CDS-klYQtA', 'gs://cclebams/wes/SANGER_HN_UPPER_AERODIGESTIVE_TRACT.bam',
 'CDS-u9hZ60' 'gs://cclebams/wes/SANGER_NCIH1755_LUNG.bam',

- remap hg38 removed WES, and readd size, hash and prevCRC32 for these lines
{'CDS-0M43vk',
 'CDS-0c0cBO',
 'CDS-0lrlra',
 'CDS-0nxrUZ',
 'CDS-0pS7rH',
 'CDS-13Wsxa',
 'CDS-1CmPe1',
 'CDS-2C2P2k',
 'CDS-2HO10g',
 'CDS-2aFKN1',
 'CDS-49azaP',
 'CDS-4IhPBG',
 'CDS-4Z3flU',
 'CDS-4c306A',
 'CDS-519oCd',
 'CDS-53dHII',
 'CDS-54TzK8',
 'CDS-5J9sBu',
 'CDS-5LolJ4',
 'CDS-5x4qLj',
 'CDS-5zTVEn',
 'CDS-6PZKz8',
 'CDS-6l3V79',
 'CDS-6vimSZ',
 'CDS-7Kh68B',
 'CDS-7PFldq',
 'CDS-7f41z4',
 'CDS-7v1juQ',
 'CDS-865Noq',
 'CDS-8GqFo5',
 'CDS-8KOsf2',
 'CDS-8jIuUa',
 'CDS-8m3WyO',
 'CDS-8sQWae',
 'CDS-9VWmlK',
 'CDS-9XPgHB',
 'CDS-9u5DMn',
 'CDS-9xyTMv',
 'CDS-9zidMf',
 'CDS-A58ERr',
 'CDS-AJMYsd',
 'CDS-AaVRuK',
 'CDS-AcA4C8',
 'CDS-B1yOQi',
 'CDS-B2rI2R',
 'CDS-BXilIo',
 'CDS-BksHHX',
 'CDS-BsqCMO',
 'CDS-C2RlCj',
 'CDS-CRPZeK',
 'CDS-CZstO2',
 'CDS-CgJtOU',
 'CDS-DledZB',
 'CDS-EEm0MQ',
 'CDS-EefiWg',
 'CDS-Eqr06D',
 'CDS-EzZEgz',
 'CDS-F7DKFi',
 'CDS-FInR9b',
 'CDS-Fr9AXu',
 'CDS-Fyjj8I',
 'CDS-GNOJc5',
 'CDS-GUJUKu',
 'CDS-GZ2d8P',
 'CDS-HPFNQJ',
 'CDS-HbJiG3',
 'CDS-Hj3xAa',
 'CDS-HjGCvC',
 'CDS-HpdIND',
 'CDS-HyNpTt',
 'CDS-I6G2pp',
 'CDS-IQuj9W',
 'CDS-Igadct',
 'CDS-JNpgXf',
 'CDS-JPXwPO',
 'CDS-JuKIq7',
 'CDS-KhywM7',
 'CDS-L4Qt3X',
 'CDS-L7Z07J',
 'CDS-LVeuLY',
 'CDS-Lck6JE',
 'CDS-LtJKmx',
 'CDS-Mzr7zM',
 'CDS-N6nU6v',
 'CDS-NE61ys',
 'CDS-NvYPYt',
 'CDS-OApL5r',
 'CDS-OEAx8K',
 'CDS-OHN1rp',
 'CDS-PETJ2c',
 'CDS-PUjlMd',
 'CDS-Pg9VCG',
 'CDS-PlugSJ',
 'CDS-Q3NFuZ',
 'CDS-Q7a4HW',
 'CDS-QO3CaX',
 'CDS-QkW8kw',
 'CDS-RHwFnF',
 'CDS-RMGYeV',
 'CDS-RQfFDW',
 'CDS-Rb0qOn',
 'CDS-RhqRTK',
 'CDS-RwI5Zp',
 'CDS-SI20Wi',
 'CDS-Sl300T',
 'CDS-Sqq5b3',
 'CDS-T10Uph',
 'CDS-TDblpN',
 'CDS-TNKqSC',
 'CDS-TQfoys',
 'CDS-TSDUCK',
 'CDS-UC4ACX',
 'CDS-Uh1A6J',
 'CDS-UtrDTK',
 'CDS-VS9XDY',
 'CDS-VXUtqm',
 'CDS-VdAsYO',
 'CDS-VseK7Y',
 'CDS-VuyYNI',
 'CDS-WAPQGk',
 'CDS-WHZolj',
 'CDS-WHdVXR',
 'CDS-WmVwp8',
 'CDS-XWBhIr',
 'CDS-XiAraT',
 'CDS-YP761S',
 'CDS-YnodyM',
 'CDS-ZAVZwo',
 'CDS-Zvb39e',
 'CDS-aXqwpM',
 'CDS-adi8ww',
 'CDS-b7i84A',
 'CDS-bntBUl',
 'CDS-bnzaOP',
 'CDS-cFqDEh',
 'CDS-cKMeDY',
 'CDS-cMQr77',
 'CDS-cmV75B',
 'CDS-cqP17c',
 'CDS-d2xe4x',
 'CDS-dGaHMd',
 'CDS-dOucLc',
 'CDS-dQKiht',
 'CDS-dbb9ty',
 'CDS-doGTxL',
 'CDS-e8OmCc',
 'CDS-eHF4ky',
 'CDS-eZQaqH',
 'CDS-eb9th6',
 'CDS-eyELts',
 'CDS-ezLdbO',
 'CDS-fE2Hcg',
 'CDS-fUKMM7',
 'CDS-fb9VZf',
 'CDS-gB4IzD',
 'CDS-gDwmIQ',
 'CDS-gE7gGi',
 'CDS-gHBkmw',
 'CDS-gPkQ5J',
 'CDS-hMjWFh',
 'CDS-hvqBN5',
 'CDS-iP4meB',
 'CDS-iXhz1w',
 'CDS-iYXN4r',
 'CDS-iymIlk',
 'CDS-jbtmKO',
 'CDS-k6MXSN',
 'CDS-kNHP4Y',
 'CDS-kUHdCP',
 'CDS-khPv3M',
 'CDS-klYQtA',
 'CDS-ldZUrx',
 'CDS-liTQH6',
 'CDS-mH7O2b',
 'CDS-mazUYU',
 'CDS-meB55t',
 'CDS-mhHwwD',
 'CDS-miDyKo',
 'CDS-ms8hJc',
 'CDS-n3RkdQ',
 'CDS-nYIBWR',
 'CDS-nkHgKT',
 'CDS-oQeSKR',
 'CDS-oRM8DN',
 'CDS-oi7eGQ',
 'CDS-pcGs7k',
 'CDS-picEuX',
 'CDS-q0NBsY',
 'CDS-qDNbzG',
 'CDS-qitYY1',
 'CDS-rDkmLA',
 'CDS-rsq3pa',
 'CDS-rvqgUM',
 'CDS-rxFkEN',
 'CDS-sEOmEs',
 'CDS-sVHzxI',
 'CDS-sieIuO',
 'CDS-tBTWSF',
 'CDS-tDo7ZF',
 'CDS-tGFzkR',
 'CDS-tJlDZV',
 'CDS-tVykaB',
 'CDS-tbW4cG',
 'CDS-u9hZ60',
 'CDS-uOnYpA',
 'CDS-uqhYTt',
 'CDS-vDTrlN',
 'CDS-vGT14g',
 'CDS-vUJJJp',
 'CDS-vYJEbO',
 'CDS-vcC6c0',
 'CDS-vdjaWG',
 'CDS-vupX1R',
 'CDS-wNFQSF',
 'CDS-wQyAqu',
 'CDS-wRs9eP',
 'CDS-wSV3OM',
 'CDS-waif9Z',
 'CDS-wkBjwZ',
 'CDS-wr5fSI',
 'CDS-x0Y44f',
 'CDS-xB0mQ4',
 'CDS-xCdJnR',
 'CDS-xiWEJn',
 'CDS-yLlkxa',
 'CDS-z6YDCo',

- reindex lost index WES (do gsutil ls on all bai files)
- Why do we say WES when we don't have WES.

- CDS-XPrQp9 is ACH-001183

- Merge lines that are exactly the same ones
[('ACH-000004',2-3
 ('ACH-000007',2-3
 ('ACH-000019',2-3
 ('ACH-000047',2-3
 ('ACH-000200',2-3
 ('ACH-000217',2-3
 ('ACH-000247',2-3
 ('ACH-000304',2-3
 ('ACH-000312',2-3
 ('ACH-000337',2-3
 ('ACH-000434',2-3
 ('ACH-000452',2-3
 ('ACH-000466',
  {1: 0.8482086730250371, 2: 0.8482086730250371, 3: 0.17625973821393717}),
 ('ACH-000514',
  {2: 0.896832622636499, 1: 0.6656716042615397, 3: 0.896832622636499}),
 ('ACH-000544',
  {2: 0.9345038866482839, 1: 0.7471091890352364, 3: 0.9345038866482839}),
 ('ACH-000557',
  {1: 0.817295785095991, 2: 0.817295785095991, 3: 0.14246048034733536}),
 ('ACH-000577',
  {2: 0.9321428128338566, 1: 0.8971449135860029, 3: 0.9321428128338566}),
 ('ACH-000596',
  {2: 0.9743359998983518, 1: 0.8686719715398205, 3: 0.9743359998983518}),
 ('ACH-000655',
  {2: 0.789934172750995, 1: 0.6330873174573377, 3: 0.789934172750995}),
 ('ACH-000672',
  {2: 0.9552560167698771, 1: 0.8250994080714341, 3: 0.9552560167698771}),
 ('ACH-000698',
  {2: 0.8965929815985778, 1: 0.7640744125107671, 3: 0.8965929815985778}),
 ('ACH-000740',
  {2: 0.9185670766468564, 1: 0.8239993288063965, 3: 0.9185670766468564}),
 ('ACH-000767',
  {2: 0.8403397440303697, 1: 0.7746706486003392, 3: 0.8403397440303697}),
 ('ACH-000800',
  {2: 0.9083513826565347, 1: 0.7716243079907239, 3: 0.9083513826565347}),
 ('ACH-000823',
  {2: 0.91066140842814, 1: 0.8448973049162971, 3: 0.91066140842814}),
 ('ACH-000837',
  {2: 0.7040218847948029, 1: 0.72967955505462, 3: 0.7040218847948029}),
 ('ACH-000866',
  {2: 0.6571284229701372, 1: 0.25945318745529355, 3: 0.6571284229701372}),
 ('ACH-001075',
  {2: 0.7301852475512532, 1: 0.4658622928321053, 3: 0.7301852475512532}),
 ('ACH-001129',
  {2: 0.8863396266123929, 1: 0.852985394639509, 3: 0.852985394639509}),
 ('ACH-001151',
  {2: 0.03968715803945472, 1: 0.8717807211442776, 3: 0.03968715803945472}),
 ('ACH-001190',
  {2: 0.8539697921208084, 1: 0.725892742242607, 4: 0.8539697921208084}),
 ('ACH-001333',
  {2: 0.7067981212273347, 1: 0.06967689262994049, 3: 0.7067981212273347}),
 ('ACH-001336',
  {2: 0.8960620737573596, 1: 0.79928494127311, 3: 0.8960620737573596}),
 ('ACH-001341',
  {2: 0.9077597944745668, 1: 0.8498669306518626, 3: 0.9077597944745668}),
 ('ACH-001345',
  {2: 0.6367572702292515, 1: -0.10679906674707781, 3: 0.6367572702292515}),
 ('ACH-001360',
  {2: 0.925375697750159, 1: 0.7933324784501316, 3: 0.925375697750159}),
 ('ACH-001368',
  {2: 0.7641925801439161, 1: 0.4003030833367517, 3: 0.7641925801439161}),
 ('ACH-001373',
  {3: 0.9755659962415464,
   2: 0.9720263251667639,
   1: 0.8910183348945192,
   4: 0.9755659962415464}),
 ('ACH-001374',
  {2: 0.6072688939068169, 1: 0.011929775387380263, 3: 0.6072688939068169}),
 ('ACH-001401',
  {2: 0.9257412538154829, 1: 0.7844764808093024, 3: 0.9257412538154829}),
 ('ACH-001402',
  {2: 0.9662469505742706, 1: 0.9066334074108684, 3: 0.9662469505742706}),
 ('ACH-001418',
  {2: 0.9401133945270204, 1: 0.8689310043704049, 3: 0.9401133945270204}),
 ('ACH-001443',
  {2: 0.9202540211940454, 1: 0.8431159318525633, 3: 0.9202540211940454}),
 ('ACH-001496',
  {2: 0.9196081155073421, 1: 0.8348312747913919, 3: 0.9196081155073421}),
 ('ACH-001497',
  {2: 0.9355132686608456, 1: 0.8009847601911252, 3: 0.9355132686608456}),
 ('ACH-001500',
  {2: 0.9368175705155969, 1: 0.786239933999951, 3: 0.9368175705155969}),
 ('ACH-001517',
  {2: 0.901567503853154, 1: 0.753937358039363, 3: 0.901567503853154}),
 ('ACH-001525',
  {2: 0.90084507431811, 1: 0.7996988222355362, 3: 0.90084507431811}),
 ('ACH-001526',
  {2: 0.9282332750686934, 1: 0.8459846806687525, 3: 0.9282332750686934}),
 ('ACH-001530',
  {2: 0.8040588835425967, 1: 0.6352248929938014, 3: 0.8040588835425967}),
 ('ACH-001542',
  {2: 0.920379261709659, 1: 0.7935286074774694, 3: 0.920379261709659}),
 ('ACH-001549',
  {2: 0.9518898206744432, 1: 0.8369508905042846, 3: 0.9518898206744432}),
 ('ACH-001630',
  {2: 0.9465769352906812, 1: 0.7944170899035772, 3: 0.9465769352906812}),
 ('ACH-001638',
  {2: 0.8327515519793249, 1: 0.34862334911983683, 3: 0.8327515519793249}),
 ('ACH-001642',
  {2: 0.7221686848532586, 1: 0.4953043680014088, 3: 0.7221686848532586}),
 ('ACH-001650',
  {2: 0.8317862067475097, 1: 0.05648450497066899, 3: 0.8317862067475097}),
 ('ACH-001654',
  {2: 0.9524797570871436, 1: 0.8600279093012886, 3: 0.9524797570871436}),
 ('ACH-001655',
  {2: 0.8638769601541457, 1: 0.7724983130789185, 3: 0.8638769601541457}),
 ('ACH-001674',
  {2: 0.8414666641804616, 1: 0.672633700707711, 3: 0.8414666641804616}),
 ('ACH-001702',
  {2: 0.800196442501435, 1: 0.6811619549706996, 3: 0.800196442501435}),
 ('ACH-001715',
  {2: 0.9521698698216043, 1: 0.822900711003482, 3: 0.9521698698216043}),
 ('ACH-001794',
  {1: 0.8378041545446324, 2: 0.35683526310261676, 3: 0.8378041545446324}),
 ('ACH-001861',
  {2: 0.8272623487299418, 1: 0.6741365308264543, 3: 0.8272623487299418}),
 ('ACH-002446',
  {1: 0.9088843662141951, 2: 0.9088843662141951, 3: 0.8484310922357836}),
 ('ACH-001398',
  {2: 0.7860358137582468, 1: 0.6338261302739719, 3: 0.7860358137582468}),
 ('ACH-001453',
  {2: 0.857926178982342, 1: 0.6837877405158376, 3: 0.857926178982342}),
 ('ACH-001498',
  {2: 0.6772635053741624, 1: -0.1852267397593894, 3: 0.6772635053741624}),
 ('ACH-001617',
  {2: 0.5984826612145459, 1: 0.17532580322480878, 3: 0.5984826612145459}),
 ('ACH-001653',
  {2: 0.8407073971972692, 3: 0.8403345124927907, 1: 0.8407073971972692}),
 ('ACH-001670',
  {2: 0.8901028330106359, 1: 0.7990857937751111, 3: 0.8901028330106359}),
 ('ACH-001819',
  {2: 0.9489571558297776, 1: 0.802459440416392, 3: 0.9489571558297776}),
 ('ACH-002029',
  {2: 0.8698124730675842, 1: 0.7765321768541431, 3: 0.8698124730675842}),
 ('ACH-002059',
  {2: 0.7852957441863467, 1: 0.19317501411009305, 3: 0.7852957441863467}),
 ('ACH-002062',
  {2: 0.7713690661944266, 1: 0.5697895015503149, 3: 0.7713690661944266}),
 ('ACH-001421',
  {1: 0.7740459556625691,
   2: 0.7740459556625691,
   3: 0.7740459556625691,
   4: 0.7740459556625691}),
 ('ACH-001711',
  {2: 0.8584477973935452,
   3: 0.8584477973935452,
   1: 0.6450068591484311,
   4: 0.8584477973935452,
   5: 0.8584477973935452}),
 ('ACH-001679',
  {2: 0.6410030046666684, 3: 0.6442528562808794, 4: 0.6410030046666684})]

- take care of dup lines
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
- rename lines that had duplicate namings
MDAMB435S MDA-MB-435S
{'ACH-000884'}
U-178 U178
{'ACH-000208'}
SW1573 SW 1573
{'ACH-000677'}
KYSE140 KYSE-140
{'ACH-000823'}
RI-1 RI1
{'ACH-000419'}
NMC-G1 NMCG1
{'ACH-000200'}
DMS53 DMS 53
{'ACH-000698'}
NCIH2882 NCI-H2882
{'ACH-000700'}
HEL 92.1.7 HEL9217
{'ACH-000005'}
LN235 LN-235
{'ACH-000591'}
HUT78 HuT 78
{'ACH-000509'}
HEC-265 HEC265
{'ACH-000946'}
WSU-NHL WSUNHL
{'ACH-001709'}
PA-1 PA1
{'ACH-001374'}
NCIH647 NCI-H647
{'ACH-000378'}
NCI-H1915 NCIH1915
{'ACH-000434'}
SNU840 SNU-840
{'ACH-000280'}
SNU-761 SNU761
{'ACH-000537'}
JHH2 JHH-2
{'ACH-000577'}
SNU-1077 SNU1077
{'ACH-000302'}
NCI-H1339 NCIH1339
{'ACH-000921'}
NCIH2286 NCI-H2286
{'ACH-000912'}
HCC827GR5 HCC827 GR5
{'ACH-000029'}
LN-464 LN464
{'ACH-000676'}
SCCOHT-1 SCCOHT1
{'ACH-001279'}
NCI-H1092 NCIH1092
{'ACH-000514'}
SK-MEL-2 SKMEL2
{'ACH-001190'}
SNU886 SNU-886
{'ACH-000316'}
LN-340 LN340
{'ACH-000634'}
RBE RBE_
{'ACH-001856'}
OVCAR-5 OVCAR5
{'ACH-001151'}
NCI-H446 NCIH446
{'ACH-000800'}
NCI-H3122 NCIH3122
{'ACH-000337'}
C33A C-33 A
{'ACH-001333'}
SKRC31 SK-RC-31
{'ACH-001194'}
RERFLCAD2 RERF-LC-Ad2
{'ACH-000774'}
OCIMY7 OCI-My7
{'ACH-000436'}
MONOMAC1 MONO-MAC-1
{'ACH-001129'}
SNU-601 SNU601
{'ACH-000736'}
SNU-1196 SNU1196
{'ACH-000461'}
KCIMOH1 KCI-MOH1
{'ACH-001098'}
SUPT11 SUP-T11
{'ACH-000122'}
F36P F-36P
{'ACH-000487'}
JHOM1 JHOM-1
{'ACH-000237'}
RERFLCAD1 RERF-LC-Ad1
{'ACH-000791'}
OCUM1 OCUM-1
{'ACH-000247'}
NCIH322 NCI-H322
{'ACH-000837'}
IA-LM IALM
{'ACH-000672'}
TE-8 TE8
{'ACH-000452'}
SNU-216 SNU216
{'ACH-000466'}
SNU-719 SNU719
{'ACH-000898'}
NCIH2887 NCI-H2887
{'ACH-000251'}
YD15 YD-15
{'ACH-000836'}
IMR32 IMR-32
{'ACH-000310'}
KPL1 KPL-1
{'ACH-000028'}
A253 A-253
{'ACH-000740'}
NCIH526 NCI-H526
{'ACH-000767'}
LN-443 LN443
{'ACH-000673'}
KE97 KE-97
{'ACH-000167'}
DLD1 DLD-1
{'ACH-001061'}
CASKI Ca Ski
{'ACH-001336'}
BC3C BC-3C
{'ACH-000593'}
PA-TU-8988S PATU8988S
{'ACH-000022'}
NCIH292 NCI-H292
{'ACH-001075'}
NCI-H1048 NCIH1048
{'ACH-000866'}
HEC1B HEC-1-B
{'ACH-000941'}
NCI-H2077 NCIH2077
{'ACH-000010'}
WM115 WM-115
{'ACH-000304'}
AML193 AML-193
{'ACH-000557'}
JHUEM7 JHUEM-7
{'ACH-000993'}
OUMS-23 OUMS23
{'ACH-000296'}
HCC366 HCC-366
{'ACH-000840'}
SNU-8 SNU8
{'ACH-000460'}

## look and get all homonyms from arxspan

# Other:

ChangeLog:
- [CN] we changed the way we compute gene level copy number. From an average of all segments within gene to a weighted average based on segment length within genes.
- [CN] removing ACH-001189, ACH-002303, ACH-002315, ACH-002341 as they are duplicates from other cell lines.
- [CN] removed logic for CN values being clamped at a certain max value
- [CN] adding amplification status in segmented CN (+,-,0) as a new columns (directly from the GATK pipeline)
- [CN] WGS copy numbers. WGS has replaced the lines that had either SNP array or WES. WGS provides data outside of the exonic regions (whereas it is infered from exonic region in WES)
- [CN] Y chromosome in CN. We are now reporting Y chromosome copy number data in some of our samples. Doing so we had to change our Panel Of Normals (changing the CN ratios slightly) for these samples. This will introduce a batch effect in this release.
- [CN] renaming some segment level sequencing type from Broad WES to Broad SNP	as it seems they came from SNP array (based on the baam files we have)

- [Mutations] changing mutation columns from wgs to legacy_wgs_exon_only to reflect that these are WGS mutations on exonic region called with out legacy pipeline. CGA_WES_AC now contains exonic mutation using the CGA WES pipeline. We are using WGS in this pipeline if no WES is available and we are treating WGS as if they were WES.
- [Mutations] solved issues with some RNA mutations being reported twice
- [Mutations] removed some 0:0 allelic ratio for some mutations

- [Expression] adding new genes that were removed before (ERCC genes for examples) (some of them don't have an hgnc name so they will only have their ensembl_id)
- [Expression] removing some wrong genes (not listed in hgnc and not expressed anywhere)
- [Expression] removing 10 lines that did not pass QC (based on https://github.com/getzlab/rnaseqc/blob/master/Metrics.md, using QCs of minmapping: 0.7, minendmapping: 0.66, minefficiency: 0.6, maxendmismatch: 0.02, maxmismatch: 0.02, minhighqual: 0.7, minexon: 0.66, maxambiguous: 0.1, maxsplits: 0.1, maxalt: 0.5, maxchim: 0.2, minreads: 20000000, minlength: 80, maxgenes: 35000, mingenes: 10000, see RNAseq Readmes for moree information.


- [Expression] fully adding the reprocessed samples. we had halted the release of differrent samples in 20Q2 and 20Q3.
- [Expression] two new datasets: expression_proteincoding_genes_expected_count: the subset of only protein coding genes for expected counts. expression_transcripts_expected_count: the expected counts at the transcript level.
- [Expression] missing a line for expression_proteincoding_genes_expected_count & expression_transcripts_expected_count. This is a missing line that was added from our legcy datasets (datasets that have unreproducible data), the two datasets are new for this release and thus don't contain this line

- [Fusions] No new cell lines available due to an issue in the processing pipeline.
- [Fusions] better post processing. some lines that were initially dropped got added back in the fusion dataset


## TODO:

- remove duplicates in legacies
- give missing lines to Becky (not WES, lost, not in WES/RNA, failed QC)

