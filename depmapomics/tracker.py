#tracker.py
from genepy.utils import helper as h
import ipdb
import pandas as pd
from depmapomics import loading
from gsheets import Sheets
from taigapy import TaigaClient
tc = TaigaClient()

def merge(tracker, new, old, arxspid, cols):
  """
  given a tracker a a new and old arxspan id, will merge the two cells lines in the tracker
  """
  #loc = tracker[tracker[arxspid]==old].index


def updateFromTracker(samples, ccle_refsamples, arxspan_id='arxspan_id', participant_id='participant_id', toupdate={"sex":[],
  "primary_disease":[],
  "cellosaurus_id":[],
  "age":[],
  "primary_site":[],
  "subtype":[],
  "subsubtype":[],
  "origin":[],
  "parent_cell_line":[],
  "matched_normal":[],
  "comments":[],
  "mediatype":[],
  "condition":[],
  'stripped_cell_line_name':[],
  "participant_id":[]}):
    # If I have a previous samples I can update unknown data directly
    index = []
    notfound = []
    for k, val in samples.iterrows():
        dat = ccle_refsamples[ccle_refsamples[arxspan_id] == val[arxspan_id]]
        if len(dat) > 0:
            index.append(k)
            for k, v in toupdate.items():
                toupdate[k].append(dat[k].tolist()[0])
        else:
            notfound.append(k)
    # doing so..
    import pdb;pdb.set_trace()
    for k, v in toupdate.items():
        samples.loc[index, k] = v
    len(samples.loc[notfound][participant_id]
        ), samples.loc[notfound][participant_id].tolist()
    return samples, notfound    


def removeOlderVersions(names, refsamples, arxspan_id="arxspan_id", 
                        version="version", priority=None):
  """
  Given a dataframe containing ids, versions, sample_ids and you dataset df indexed by the same ids, will set it to your sample_ids using the latest version available for each sample

  Args:
  -----
    refsamples: df[id, version, arxspan_id,...] the reference metadata
    names: list[id] only do it on this set of samples
    arxspan_id: the name of the id field
    version: the name of the version field

  Returns:
  --------
    the subsetted dataframe

  """
  lennames = len(names)
  res = {}
  refsamples = refsamples.loc[names].copy()
  if lennames > len(refsamples):
    print(set(names) - set(refsamples.index))
    ipdb.set_trace()
    raise ValueError('we had some ids in our dataset not registered in this refsample dataframe')
  for arxspan in set(refsamples[arxspan_id]):
    allv = refsamples[refsamples[arxspan_id] == arxspan]
    for k, val in allv.iterrows():
      if priority is None:
        if val[version] == max(allv.version.values):
          res[k] = arxspan
          break
      else:
        if val[version] == max(allv.version.values):
          res[k] = arxspan
        if val[priority] == 1:
          res[k] = arxspan
          break

  print("removed " + str(lennames - len(res)) + " duplicate samples")
  # remove all the reference metadata columns except the arxspan ID
  return res



def changeCellLineNameInNew(ref, new, datatype, dupdict, toupdate=['stripped_cell_line_name',
                                                                      'arxspan_id', "patient_id",
                                                                      "sex", "primary_disease",
                                                                      "cellosaurus_id", "age",
                                                                      "primary_site", "subtype",
                                                                      "subsubtype"]):
  """
  Rename a/some line in a DF and takes care of corresponding metadata and versions from the sample tracker

  Args:
  -----
    new: change the cell line name in this dataframe
    dupdict: dict(tochange,newname)
    datatype: str for a ref with many datatype (to get the right version number)

  Returns:
  --------
    the updated dataframe
  """
  for k, v in dupdict.items():
    new.loc[new[new.arxspan_id == k].index, toupdate] = ref[ref.arxspan_id == v][toupdate].values[0]
    new.loc[new[new.arxspan_id == v].index, 'version'] = len(ref[(ref.arxspan_id == v) & (ref.datatype == datatype)]) + 1
  return new


def changeCellLineName(ref, datatype, dupdict, toupdate=["stripped_cell_line_name",
                                                         "participant_id",
                                                         "cellosaurus_id",
                                                         "sex",
                                                         "arxspan_id",
                                                         "matched_normal",
                                                         "age",
                                                         "primary_site",
                                                         "primary_disease",
                                                         "subtype",
                                                         "subsubtype",
                                                         "origin"]):
  """
  Rename a/some line in our sample tracker and takes care of corresponding metadata and versions from the sample tracker

  Args:
  -----
    dupdict: dict(tochange,newname)
    datatype: str for a ref with many datatype (to get the right version number)

  Returns:
  --------
    the updated dataframe
  """
  for k, v in dupdict.items():
    try:
      ref.loc[k, toupdate] = ref[ref.arxspan_id == v][toupdate].values[0]
      ref.loc[k, 'version'] = len(ref[(ref.arxspan_id == v) & (ref.datatype == datatype)]) + 1
    except IndexError:
      raise IndexError(str(v)+" not found in tracker")
  return ref



def updatePairs(workspaceID, tracker, removeDataFiles=True, ):
  """
  looks at the current sample tracker and updates the pairs in Terra

  It will add and remove them based on what information of match normal is available in the sample tracker. if an update happens it will remove the data files for the row.
  """

def cleanVersions(tracker, samplecol='arxspan_id', dryrun=False, datatypecol='datatype', versioncol="version"):
  """
  updates and sorts the versions of samples in the sample tracker:

  checks that we get 1,2,3 instead of 2,4,5 when samples are renamed or removed

  Args:
  -----
    tracker: dataframe of the sample tracker
    samplecol: str colname of samples in the trackerr
    dryrun: bool whether or not to apply it or just print if there is issues
    datatypecol: str colname of the datatype values
    versioncol: str colname of the version values

  Returns:
  -------
    an updated sample tracker
  """
  tracker = tracker.copy()
  tracker['samval']=tracker[samplecol]+tracker[datatypecol]
  for v in set(tracker['samval']):
      sams = tracker[tracker['samval']==v]
      vs = sams[versioncol].tolist()
      if max(vs)==len(vs):
          continue
      print("found issue")
      if dryrun:
        continue
      vs.sort()
      rn = {}
      for i,v in enumerate(vs):
          rn.update({v:i+1})
      for k, val in sams.iterrows():
          tracker.loc[k, versioncol] = rn[val[versioncol]]
  tracker = tracker.drop(columns='samval')
  return tracker


def setRightName(tracker, name='stripped_cell_line_name', signs=['-','_','.',' ']):
  """
  cell line name need to be all uppercase and not contain special signs will take carre of ones that are not

  BE CARREFUL, it does not solve issues of duplicate lines (diff sample id, same new name)
  see findLikelyDup for that

  Args:
  -----
    tracker: dataframe of the sample tracker

  Returns:
  -----
    an updated sample tracker
  """
  new = []
  for val in tracker[name]:
      for s in signs:
          val=val.replace(s,'')
      new.append(val.upper())
  tracker[name]=new
  return tracker


def findLikelyDup(tracker, name='stripped_cell_line_name', signs=['-', '_', '.', ' '], arxspid='arxspan_id', looksub=True):
  """
  find cell lines that are likely to be duplicates

  will return ,  as well,
  and

  Args:
  -----
    tracker: dataframe of the sample tracker
    looksub: bool, look if a name if within another name (can flag many derivatives)

  Returns:
  --------
    a list[tuples(str,str)] of likly duplicate names as tuples (rh13, RH-13)
    a list[tuples(str,str)] of associated arxspan ids
    a dict[str:set(str)] of arxspan ids that have multiple cell line names associated
  """
  names = set(tracker[name])
  simi = []
  arxsp = []
  issues = {}
  for i, name1 in enumerate(names):
    h.showcount(i, len(names))
    n1 = name1
    for s in signs:
      name1 = name1.replace(s, '')
    name1 = name1.upper()
    for name2 in names-set([n1]):
      n2 = name2
      for s in signs:
        name2 = name2.replace(s, '')
      name2 = name2.upper()
      if name1 == name2:
        if (looksub and (name1 in name2 or name2 in name1) and abs(len(name1)-len(name2)) < 2) or not looksub:
          if (n1, n2) not in simi and (n2, n1) not in simi:
            simi.append((n1, n2))
            arxsp.append(
                (tracker[tracker[name] == n1][arxspid][0], tracker[tracker[name] == n2][arxspid][0]))
  for val in set(tracker[name]):
    v = set(tracker[tracker[name] == val][arxspid])
    if len(v) > 1:
      issues.update({val: v})
  return simi, arxsp, issues


def resolveIssues(tracker, issus, arxspid, cols):
  """
  given a dict of names: [arxp ids] will try to find back the right name for the right
  arxspan id by looking at their rfequncy of occurance along the tracker

  if we have rh12: [ACH-00001,ACH-0002]
  and rh12 is associated 1 time with ach-00002 and 3 with ach-00001
  and rh13 is assocated 2 time with ach-00002, then it associates:
  ach-00001 : rh12
  ach-00002 : rh13
  """
  #for val in issus:


def retrieveFromCellLineName(noarxspan, ccle_refsamples, datatype, extract={}, my_id='~/.client_secret.json',
                            stripped_cell_line_name="stripped_cell_line_name", arxspan_id="arxspan_id", 
                            mystorage_id="~/.storage.json",
                            depmappvlink="https://docs.google.com/spreadsheets/d/1uqCOos-T9EMQU7y2ZUw4Nm84opU5fIT1y7jet1vnScE"):
  """
  Given a List of samples with no arxspan ids, will try to retrieve an arxspan id and associated data from trackers

  For now the sample tracker and paquita's depmap pv are used.

  Args:
  -----
    noarxspan: dataframe of samples with missing arxspan ids
    ccle_refsamples: dataframe of the sample tracker
    datatype: str the datatype we are interested in (wes/rna/..) see the ones in the sample tracker
    extract: dict(str:str) see the extract in the "resolveFromWorkspace" function
    stripped_cell_line_name: str colname where the cell line name is stored
    arxspan_id: str colname wherre the sample id is stored in both noarxspan and ccle_refsamples
    depmappvlink: str the url to the depmap_pv google sheet

  Returns:
  --------
    a new dataframe with filled annotations when possible
  """
  sheets = Sheets.from_files(my_id, mystorage_id)

  # find back from cell line name in ccle ref samples
  noarxspan.arxspan_id = [ccle_refsamples[ccle_refsamples[stripped_cell_line_name] == i].arxspan_id[0] if i in ccle_refsamples[stripped_cell_line_name].tolist() else 0 for i in noarxspan[arxspan_id]]
  a = [ccle_refsamples[ccle_refsamples[stripped_cell_line_name] == i][arxspan_id][0]
        if i in ccle_refsamples[stripped_cell_line_name].tolist() else 0 for i in noarxspan[stripped_cell_line_name]]
  noarxspan[arxspan_id] = [i if i != 0 else a[e]
                            for e, i in enumerate(noarxspan.arxspan_id)]

  # get depmap pv
  depmap_pv = sheets.get(depmappvlink).sheets[0].to_frame(header=2)
  depmap_pv = depmap_pv.drop(depmap_pv.iloc[:1].index)

  # find back from depmapPV
  signs = ['-', '_', '.', ' ']
  for k, val in noarxspan[noarxspan[arxspan_id] == 0].iterrows():
    val = val[stripped_cell_line_name].upper()
    for s in signs:
      val = val.replace(s, '')
    a = depmap_pv[depmap_pv['CCLE_name'].str.contains(
        val) | depmap_pv['Stripped Cell Line Name'].str.contains(val) | depmap_pv['Aliases'].str.contains(val)]
    if len(a) > 0:
      noarxspan.loc[k, arxspan_id] = a['DepMap_ID'].values[0]
  noarxspan[arxspan_id] = noarxspan[arxspan_id].astype(str)
  new_noarxspan= loading.resolveFromWorkspace(noarxspan[noarxspan[arxspan_id].str.contains('ACH-')], refsamples=ccle_refsamples[ccle_refsamples['datatype'] == datatype], match=[
                                    'ACH', 'CDS'], participantslicepos=10, accept_unknowntypes=True, extract=extract)
  return pd.concat([new_noarxspan, noarxspan[~noarxspan[arxspan_id].str.contains('ACH-')]])


def updateSamplesSelectedForRelease(refsamples, releaseName, samples):
  """
  given a list of samples, a release name and our sample tracker, will set these samples as 1 for this releasename and the rest at 0
  """
  refsamples[releaseName] = '0'
  refsamples.loc[samples, releaseName] = '1'
  return refsamples
