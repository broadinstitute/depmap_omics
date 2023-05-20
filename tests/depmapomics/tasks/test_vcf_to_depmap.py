import os
import subprocess
import json
from dataclasses import dataclass
from typing import Dict


@dataclass
class WDLResult:
    original_output: Dict
    files: Dict[str, str]
    properties: Dict[str, str]


def run_wdl(tmpdir, wdl_path, inputs):
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
        original_outputs = json.load(fd)

    def is_file_path(value):
        return value.startswith(dest_dir) and os.path.exists(value)

    # go through the outputs and identify those which are file paths vs those that are strings
    files = {}
    properties = {}
    for name, value in original_outputs["outputs"].items():
        if isinstance(value, list):
            if all([is_file_path(x) for x in value]):
                files[name] = [x[len(dest_dir) + 1 :] for x in value]
            else:
                properties[name] = value
        elif isinstance(value, str):
            if is_file_path(value):
                files[name] = value[len(dest_dir) + 1 :]
            else:
                properties[name] = value

    return WDLResult(
        original_output=original_outputs, files=files, properties=properties
    )


def _file_sha256_hash(fd, bufsize=100_000):
    import hashlib

    hasher = hashlib.sha256()
    while True:
        data = fd.read(bufsize)
        if not data:
            break
        hasher.update(data)
    return hasher.hexdigest()


import gzip


def assert_gziped_files_equal(expected_file, new_file):
    def get_ungzip_hash(path):
        with gzip.open(path, "rb") as fd:
            return _file_sha256_hash(fd)

    assert get_ungzip_hash(expected_file) == get_ungzip_hash(
        new_file
    ), f"file {expected_file} is different then {new_file}"


import pandas as pd
import pandas.testing as testing


def assert_parquet_files_equal(expected_file, new_file):
    testing.assert_frame_equal(
        pd.read_parquet(expected_file), pd.read_parquet(new_file)
    )


def assert_output_files_match(expected_file, new_file):
    # Open question: How much information do we want to provide if there's a mismatch
    if expected_file.endswith(".csv.gz"):
        assert_gziped_files_equal(expected_file, new_file)
    elif expected_file.endswith(".parquet"):
        assert_parquet_files_equal(expected_file, new_file)
    else:
        raise Exception("Did not know how to compare {expected_file} and {new_file}")


def assert_output_dirs_match(expected_output_dir, new_output_dir):
    def read_json(path):
        with open(os.path.join(path, "output.json"), "rt") as fd:
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
            if os.path.exists(os.path.join(expected_output_dir, expected_value)):
                # check if these are files which should be compared
                assert os.path.exists(os.path.join(new_output_dir, new_value))

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


def test_vcs_to_depmap(tmpdir):
    result = run_wdl(
        tmpdir,
        "WGS_pipeline/vcf_to_depmap.wdl",
        {"input_vcf": "tests/CDS-R99tiF.subset.vcf.gz", "sample_id": "CDS-R99tiF"},
    )

    output = result.original_output

    write_outputs(output["dir"], output["outputs"], "test_vcs_to_depmap/latest")

    assert_output_dirs_match("test_vcs_to_depmap/expected", "test_vcs_to_depmap/latest")
    # print(result)

    # def name_with_extension(name, orig_name):
    #     if "." not in orig_name:
    #         return name
    #     extension = orig_name[orig_name.index(".") + 1 :]
    #     return f"{name}.{extension}"

    # import shutil

    # os.path.join()
    # for name, path in result.files.items():
    #     shutil.copy(
    #         path,
    #         os.path.join(
    #             "vcs_to_depmap", "outputs", "last", name_with_extension(name, path)
    #         ),
    #     )


import shutil


def write_outputs(output_dir: str, outputs: dict, dest_dir: str):
    """
    Write outputs as a json file named output.json in the dest_dir directory. Also
    copying any files referenced to a subdirectoryof dest_dir and writing the updated
    path into the json file. (In the json file the
    paths will be relative to dest_dir)
    """

    def ensure_parent_dir_exists(path):
        parent_dir = os.path.dirname(path)
        if os.path.exists(parent_dir):
            assert os.path.isdir(parent_dir)
        else:
            os.makedirs(parent_dir)

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
    print("writing ", dest_outputs)
    ensure_parent_dir_exists(dest_outputs)
    with open(dest_outputs, "wt") as fd:
        fd.write(json.dumps(new_outputs))


def test_write_outputs(tmpdir):
    # create a mock output directory with a single file in it
    output_dir = tmpdir.join("sample")
    output_dir.mkdir()
    file_path = output_dir.join("file")
    file_path.write("file")

    dest_dir = tmpdir.join("dest")
    write_outputs(
        str(output_dir),
        {"file_path": str(file_path), "str_value": "string"},
        str(dest_dir),
    )

    # now verify dest_dir contains the expected files
    dest_outputs = json.loads(dest_dir.join("output.json").read())
    assert set(dest_outputs.keys()) == {"file_path", "str_value"}
    assert dest_outputs["str_value"] == "string"
    assert dest_dir.join(dest_outputs["file_path"]).read() == "file"
