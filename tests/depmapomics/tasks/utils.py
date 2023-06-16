import os
import subprocess
import json
from typing import Dict
import shutil
import pandas as pd
import pandas.testing as testing
import hashlib
import gzip
from google.cloud.storage import Client
import re


def run_wdl(tmpdir, wdl_path, inputs, destination):
    """
    Run a WDL script given the provided inputs, and write the result in the destination directory.
    """
    dest_dir = str(tmpdir.join("work"))

    tmp_inputs = tmpdir.join("inputs.json")
    tmp_inputs.write(json.dumps(inputs))
    tmp_outputs = tmpdir.join("output.json")

    # run the WDL script, writing out to the tmp dir
    subprocess.check_call(
        [
            "miniwdl",
            "run",
            wdl_path,
            "-i",
            str(tmp_inputs),
            "-o",
            str(tmp_outputs),
            "-d",
            dest_dir,
        ]
    )

    with open(str(tmp_outputs), "rt") as fd:
        output = json.load(fd)

    # rewrite the outputs file so that it the files it references all are relative to the output directory
    # so that we can upload these to the cloud, and in general, move them around without worrying about
    # the working directory that the paths are relative to, changing.
    write_outputs(output["dir"], output["outputs"], destination)


def _file_sha256_hash(fd, bufsize=100_000):
    "Return the sha256 hash of the given file"
    hasher = hashlib.sha256()
    while True:
        data = fd.read(bufsize)
        if not data:
            break
        hasher.update(data)
    return hasher.hexdigest()


def assert_gziped_files_equal(expected_file, new_file):
    "throw an assertion if the two ungziped files are different. (Compare the uncompressed versions because there is more than one compressed encoding of a file)"

    def get_ungzip_hash(path):
        with gzip.open(_localize(path), "rb") as fd:
            return _file_sha256_hash(fd)

    assert get_ungzip_hash(expected_file) == get_ungzip_hash(
        new_file
    ), f"file {expected_file} is different then {new_file}"


def assert_parquet_files_equal(expected_file, new_file):
    testing.assert_frame_equal(
        pd.read_parquet(_localize(expected_file)), pd.read_parquet(new_file)
    )


def assert_output_files_match(expected_file, new_file):
    "Compare two files and throw an assertion if they don't match. Based on the file extension, we may do different comparisions."
    # Open question: How much information do we want to provide if there's a mismatch
    if expected_file.endswith(".csv.gz"):
        assert_gziped_files_equal(expected_file, new_file)
    elif expected_file.endswith(".parquet"):
        assert_parquet_files_equal(expected_file, new_file)
    else:
        raise Exception("Did not know how to compare {expected_file} and {new_file}")


cache_root = "gcs_cache"

_gcs_client = None


def _get_client():
    "Get a gcs client, creating it if necessary"
    global _gcs_client
    if _gcs_client is None:
        _gcs_client = Client()
    return _gcs_client


def _localize(path):
    "Given a path (local or GCS), return a path to a local file which we can read with those contents"
    m = re.match("gs://([^/]+)/(.+)", path)
    if m is None:
        return path
    else:
        assert m is not None
        client = _get_client()
        bucket_name, key = m.groups()
        bucket = client.bucket(bucket_name)
        blob = bucket.blob(key)
        local_path = f"{cache_root}/{bucket}/{key}"
        ensure_parent_dir_exists(local_path)
        blob.download_to_filename(local_path)
        return local_path


def _exists_local_or_gs(path):
    "Returns true if the path (local or GCS) exists"
    m = re.match("gs://([^/]+)/(.+)", path)
    if m is None:
        return os.path.exists(path)
    else:
        client = _get_client()
        assert m is not None
        bucket_name, key = m.groups()
        bucket = client.bucket(bucket_name)
        blob = bucket.blob(key)
        return blob.exists()


def assert_output_dirs_match(expected_output_dir, new_output_dir):
    """
    Compare a directory (either local, or a in GCS) with expected results, with the results stored
    in a local directory.
    """

    def read_json(path):
        with open(_localize(os.path.join(path, "output.json")), "rt") as fd:
            return json.load(fd)

    expected_output = read_json(expected_output_dir)
    new_output = read_json(new_output_dir)

    missing_outputs = set(expected_output.keys()).difference(new_output.keys())
    extra_outputs = set(new_output.keys()).difference(expected_output.keys())

    assert (
        len(missing_outputs) == 0
    ), f"{expected_output_dir} contained outputs named {missing_outputs} but {new_output_dir} did not"
    assert (
        len(missing_outputs) == 0
    ), f"{expected_output_dir} did not contain outputs named {extra_outputs} but {new_output_dir} did"

    # if we've reached here, they have the same keys
    for name, expected_value in expected_output.items():

        def assert_value_matches(expected_value, new_value):
            if _exists_local_or_gs(os.path.join(expected_output_dir, expected_value)):
                # check if these are files which should be compared
                assert _exists_local_or_gs(os.path.join(new_output_dir, new_value))

                assert_output_files_match(
                    os.path.join(expected_output_dir, expected_value),
                    os.path.join(new_output_dir, new_value),
                )
            else:
                # if not, compare these as strings
                assert (
                    expected_value == new_value
                ), f"Mismatch in {name}: expected {expected_value} but was {new_value}"

        new_value = new_output[name]

        assert type(new_value) == type(expected_value)
        if isinstance(new_value, list):
            assert len(expected_value) == len(new_value)
            for expected_v, new_v in zip(expected_value, new_value):
                assert_value_matches(expected_v, new_v)
        else:
            assert_value_matches(expected_value, new_value)


def ensure_parent_dir_exists(path):
    "Given a path of a file we want to write, creates the parent directory if it doesn't exist"
    parent_dir = os.path.dirname(path)
    if os.path.exists(parent_dir):
        assert os.path.isdir(parent_dir)
    else:
        os.makedirs(parent_dir)


def write_outputs(output_dir: str, outputs: dict, dest_dir: str):
    """
    Write outputs as a json file named output.json in the dest_dir directory. Also
    copying any files referenced to a subdirectoryof dest_dir and writing the updated
    path into the json file. (In the json file the
    paths will be relative to dest_dir)
    """

    def is_file_path(value: str):
        return value.startswith(output_dir) and os.path.exists(value)

    def rewrite_path(dest_dir, x: str):
        if is_file_path(x):
            dest_path = os.path.join(dest_dir, x[len(output_dir) + 1 :])
            ensure_parent_dir_exists(dest_path)
            shutil.copy(x, dest_path)
            return os.path.relpath(dest_path, dest_dir)
        else:
            return x

    # go through the outputs and identify those which are file paths vs those that are strings
    new_outputs = {}
    for name, value in outputs.items():
        if isinstance(value, list):
            value = [rewrite_path(dest_dir, x) for x in value]
        else:
            assert isinstance(value, str)
            value = rewrite_path(dest_dir, value)
        new_outputs[name] = value

    dest_outputs = os.path.join(dest_dir, "output.json")
    ensure_parent_dir_exists(dest_outputs)
    with open(dest_outputs, "wt") as fd:
        fd.write(json.dumps(new_outputs))