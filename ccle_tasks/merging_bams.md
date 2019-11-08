```python
 ccnames = {}
    ...: for k, val in cclsizes.items():
    ...: if val in ccnames:
    ...: ccnames[val].append(k)
    ...: else:
    ...: ccnames[val] = [k]
    ...:   # get all bai in tsv
```

```python
for k, val in names.items():
    ...: if len(val) > 1:
    ...: for v in val:
    ...: if 'rna' in v:
    ...: if v not in rnabams:
    ...: toremove.append(v)
    ...: elif 'wes' in v:
    ...: if v not in wesbam:
    ...: toremove.append(v)
```
```python
# get ls of all files folder
samples = os.popen('gsutil -m ls -al ' + folder + '**.bai').read().split('\n')
# compute size filepath

sizes = {'gs://' + val.split('gs://')[1].split('#')[0]: int(val.split("2019-")[0]) for val in samples[:-2]}
names = {}
for k, val in sizes.items():
  if val in names:
    names[val].append(k)
  else:
    names[val] = [k]
# get all bai in tsv
samp = dm.WorkspaceManager(workspace).get_samples()
```

```python
for k, val in samp.iloc[i + 1:].iterrows():
if val[bamcol] != 'NA' and val[baicol] != 'NA' and val[bamcol] is not None:
    # if bai has duplicate size
	new = val[bamcol]
	print(new)
	prev = new
	if newgs not in new:
		for prevgs in prevgslist:
			new = new.replace(prevgs, newgs)
if flag_non_matching:
  if 'gs://' == prev[:5]:
    if new == prev:
      flaglist.append(prev)
      continue
      code = os.system('gsutil ls ' + new)
      if code==0:
        print('just changing name')
        samp.loc[k,bamcol]=new
        samp.loc[k,baicol]=new.replace('.bam','.bai') if os.system(
'gsutil ls ' +
          new.replace('.bam','.bai')) == 0 else new.replace('.bam',
'.bam.bai')
      elif code == 256:
        print('not found')
        bai = new.replace('.bam','.bai') if os.system('gsutil ls '
+
          new.replace('.bam','.bai')) == 0 else new.replace('.bam',
'.bam.bai')
        for va in names[sizes[bai]]:
          # for all duplicate size
          # if ls bam of bai duplicate size work
          # mv bam to bampath in folder
          if '.bam' in va:
            if os.system('gsutil ls ' + va.split('.bam.bai')[0] + "
.bam") == 0:
              if va.split('.bam.bai')[0] + ".bam" in previous:
                print(str(k)+' was also found in wes... skipping '+
str(va))
                continue
              print("refound")
              samp.loc[k,bamcol]=va.split('.bam.bai')[0] + ".bam"
              samp.loc[k,baicol]=va
              break
          elif os.system('gsutil ls ' + va.split('.bai')[0] + ".bam
") == 0:
            if va.split('.bai')[0] + ".bam" in previous:
              print(str(k)+' was also found in wes... skipping '+st
r(va))
              continue
            print("refound")
            samp.loc[k,bamcol]=va.split('.bai')[0] + ".bam"
            samp.loc[k,baicol]=va
            break
      elif code == signal.SIGINT:
        print('Awakened')
        break
    else:
      print("no data for " + str(k))

``` 
# when fails
```python
for i in range(2000):
	if samp.iloc[i].name == k:
	  break
```

```python
for val in flaglist:
   if "ccle_bams" not in val:
       v = "gs://ccle_bams/wes/"+val.split('/')[-1]
       if v not in previous:
           v = v.replace('.bam','.*')
           code = os.system("gsutil -m cp "+ v +" gs://ccle_bams/rna/")
   if code == signal.SIGINT:
       print('awakened')
       break
```


# missing bams

"gs://ccle_bams/wes/DepMap_CellLine_WES_Batch3_June2019/RP-1561/Exome/R256/v1/R256.bam",
"gs://ccle_bams/wes/chordomafoundationsequencedata/BAM_Files/ExomeSeq/ExomeSeq/U-CH2_ExomeSeq_GS001438.cleaned.rg.sorted.markedDupes.realigned.recal.bam",
"gs://ccle_bams/wes/SANGER_J82_URINARY_TRACT.bam",
"gs://ccle_bams/wes/SANGER_KYSE520_OESOPHAGUS.bam",


# missing versions

