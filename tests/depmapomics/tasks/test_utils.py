from .utils import write_outputs
import json


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
