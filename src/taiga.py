#taiga.py


def addToVirtual(virtualname, folderfrom = None, files=[]):
  """
  will add some files from a taiga folder to a taiga virtual dataset folder and preserve the previous files

  Args:
  ----
    virtualname: the taiga virtual dataset name
    folderfrom: the taiga folder wher the files are
    files: a list(tuples(newfilename,prevfilename)) can be from the folder
  """
  file_dict = {}
  if len(files)>0:
    assert type(files[0]) is tuple
  if folderfrom is not None:
    versiona = max([int(i['name']) for i in tc.get_dataset_metadata(folderfrom)['versions']])
  versionb = max([int(i['name']) for i in tc.get_dataset_metadata(virtualname)['versions']])
  keep = [(i['name'], i['underlying_file_id']) for i in tc.get_dataset_metadata(virtualname, version=versionb)
          ['datasetVersion']['datafiles'] if 'underlying_file_id' in i]

  for i, val in enumerate(files):
    if "/" in val[1]:
      print("assuming "+val[1]+" to be a local file")
      file_dict.update({val[1]: val[0]})
      files.pop(i)    
    else:
      files[i] = (val[0], folderfrom + '.' + str(versiona) + '/' + val[1])
  print(files)
  tc.update_dataset(dataset_permaname=virtualname, add_taiga_ids=files,
                    upload_file_path_dict=file_dict, add_all_existing_files=True)
