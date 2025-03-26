import json
import subprocess
import os
import pandas as pd
import random
import string
import io

prevshowcount = 100


def showcount(i, size):
    """
    pretty print of i/size%, to put in a for loop
    """
    global prevshowcount
    a = 1 + int(100 * (i / size))
    if a != prevshowcount:
        print(str(a) + "%", end="\r")
        prevshowcount = a


def fileToList(filename):
    """
    loads an input file with a\\n b\\n.. into a list [a,b,..]
    """
    with open(filename) as f:
        return [val[:-1] for val in f.readlines()]


def listToFile(l, filename):
    """
    loads a list with [a,b,..] into an input file a\\n b\\n..
    """
    with open(filename, "w") as f:
        for item in l:
            f.write("%s\n" % item)


def generateGeneNames(
    # v104 corresponds to GENCODE v38 https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=genes&hgta_track=wgEncodeGencodeV38&hgta_table=wgEncodeGencodeCompV38&hgta_doSchema=describe+table+schema
    ensemble_server="http://may2021.archive.ensembl.org/biomart",
    useCache=False,
    cache_folder="/".join(__file__.split("/")[:-3]) + "/",
    attributes=[],
    default_attr=[
        "ensembl_gene_id",
        "external_gene_name",
        "hgnc_symbol",
        "gene_biotype",
        "entrezgene_id",
    ],
):
    """generate a genelist dataframe from ensembl's biomart

    Args:
        ensemble_server ([type], optional): [description]. Defaults to ENSEMBL_SERVER_V.
        useCache (bool, optional): [description]. Defaults to False.
        cache_folder ([type], optional): [description]. Defaults to CACHE_PATH.

    Raises:
        ValueError: [description]

    Returns:
        [type]: [description]
    """
    assert cache_folder[-1] == "/"

    cache_folder = os.path.expanduser(cache_folder)
    createFoldersFor(cache_folder)
    cachefile = os.path.join(cache_folder, ".biomart.csv")
    if useCache & os.path.isfile(cachefile):
        print("fetching gene names from biomart cache")
        res = pd.read_csv(cachefile)
    else:
        print("downloading gene names from biomart")
        res = _fetchFromServer(ensemble_server, default_attr + attributes)
        res.to_csv(cachefile, index=False)

    res.columns = default_attr + attributes
    if type(res) is not type(pd.DataFrame()):
        raise ValueError("should be a dataframe")

    if "external_gene_name" in default_attr and "hgnc_symbol" in default_attr:
        res = res[~(res["external_gene_name"].isna() & res["hgnc_symbol"].isna())]
        res.loc[res[res.hgnc_symbol.isna()].index, "hgnc_symbol"] = res[
            res.hgnc_symbol.isna()
        ]["external_gene_name"]

    return res


def _fetchFromServer(ensemble_server, attributes):
    # deferring the import until last possible moment because it's unclear
    # to me whether we are actually fetching data from biomart. (I don't think we should)
    from biomart import BiomartServer

    server = BiomartServer(ensemble_server)
    ensmbl = server.datasets["hsapiens_gene_ensembl"]
    res = pd.read_csv(
        io.StringIO(
            ensmbl.search({"attributes": attributes}, header=1).content.decode()
        ),
        sep="\t",
    )
    return res


def createFoldersFor(filepath):
    """
    will recursively create folders if needed until having all the folders required to save the file in this filepath
    """
    prevval = ""
    for val in os.path.expanduser(filepath).split("/")[:-1]:
        prevval += val + "/"
        if not os.path.exists(prevval):
            os.mkdir(prevval)


def dups(lst):
    """
    shows the duplicates in a list
    """
    seen = set()
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in lst if x in seen or seen.add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)


def dictToFile(d, filename):
    """
    turn a dict into a json file
    """
    with open(filename, "w") as json_file:
        json.dump(d, json_file)


def askif(quest):
    """
    asks a y/n question to the user about something and returns true or false given his answer
    """
    print(quest)
    inp = input()
    if inp in ["yes", "y", "Y", "YES", "oui", "si"]:
        return 1
    elif inp in ["n", "no", "nope", "non", "N"]:
        return 0
    else:
        return askif("you need to answer by yes or no")


def parrun(cmds, cores, add=[]):
    """
    runs a set of commands in parallel using the "&" command

    Args:
      cmds: the list of commands
      cores: number of parallel execution
      add: an additional list(len(cmds)) of command to run in parallel at the end of each parallel run
    """
    count = 0
    exe = ""
    if len(add) != 0 and len(add) != len(cmds):
        raise ValueError("we would want them to be the same size")
    else:
        addexe = ""
    fullres = []
    for i, cmd in enumerate(cmds):
        count += 1
        exe += cmd
        if len(add) != 0:
            addexe += add[i]
        if count < cores and i < len(cmds) - 1:
            exe += " & "
            if len(add) != 0:
                addexe += " & "
        else:
            count = 0
            res = subprocess.run(exe, capture_output=True, shell=True)
            if res.returncode != 0:
                raise ValueError("issue with the command: " + str(res.stderr))
            exe = ""
            if len(add) != 0:
                res = subprocess.run(addexe, capture_output=True, shell=True)
                if res.returncode != 0:
                    raise ValueError("issue with the command: " + str(res.stderr))
                addexe = ""
            fullres.append(res.stdout.decode("utf-8"))
    return fullres


def randomString(stringLength=6, stype="all", withdigits=True):
    """
    Generate a random string of letters and digits

    Args:
    -----
      stringLength: the amount of char
      stype: one of lowercase, uppercase, all
      withdigits: digits allowed in the string?

    Returns:
    -------
      the string
    """
    if stype == "lowercase":
        lettersAndDigits = string.ascii_lowercase
    elif stype == "uppercase":
        lettersAndDigits = string.ascii_uppercase
    else:
        lettersAndDigits = string.ascii_letters
    if withdigits:
        lettersAndDigits += string.digits
    return "".join(random.choice(lettersAndDigits) for i in range(stringLength))
