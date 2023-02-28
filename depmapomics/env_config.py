import os

def _file_to_list(filename):
    """
    loads an input file with a\\n b\\n.. into a list [a,b,..]
    """
    with open(filename) as f:
        return [val[:-1] for val in f.readlines()]


# if we are on master, will use the master config file except if stipulated otherwise in the DEPMAP_ENV environment variable
if (
    _file_to_list("/".join(__file__.split("/")[:-2]) + "/.git/HEAD")[0].split("/")[-1]
    == "master"
    or os.getenv("DEPMAP_ENV") == "PROD"
):
    from depmapomics.config_prod import *
else:
    from depmapomics.config_dev import *
