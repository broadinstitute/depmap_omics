import subprocess
import os

def mkdir_if_not_exist(inputdir):
    if not os.path.exists(inputdir):
        os.makedirs(inputdir)
    return inputdir

def run_unix_cmd(cmd, verbose=True):
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if verbose:
        print(output.decode())
    if p.returncode != 0:
        raise Exception(error.decode())
    return output.decode()


def gsutil_cp(path1, path2, payer_project_id=None, make_dir=False, verbose=True, path2isfile=False):
    if path2isfile:
        path_arg = ''
    else:
        path_arg = ' -r'
    if make_dir:
        mkdir_if_not_exist(path2)

    if payer_project_id==None:
        cmd = 'gsutil -m cp{:s} {:s} {:s}'.format(path_arg, path1, path2)
    else:
        cmd = 'gsutil -u {:s} -m cp{:s} {:s} {:s}'.format(payer_project_id, path_arg, path1, path2)
    if verbose:
        print('copying {:s} -> {:s}'.format(path1, path2))
    output = run_unix_cmd(cmd, verbose=verbose)
    return output
