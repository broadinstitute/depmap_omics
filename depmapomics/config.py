import os.getenv

def fileToList(filename):
  """
  loads an input file with a\\n b\\n.. into a list [a,b,..]
  """
  with open(filename) as f:
    return [val[:-1] for val in f.readlines()]


if fileToList("../.git/HEAD")[0].split('/')[-1] == "master" or os.getenv('DEPMAP_ENV') == "PROD":
  from depmapomics.config_prod import *
else:
  from depmapomics.config_dev import *
