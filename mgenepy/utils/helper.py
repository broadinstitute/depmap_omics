import json
import subprocess
import os

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
