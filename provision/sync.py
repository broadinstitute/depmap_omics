import dalmatian
import firecloud.api
from dataclasses import dataclass
from dump_workspace import ParsedWDL, dumb_wdl_parse, rewrite_with_relative_imports
from typing import Dict
import json
import os
import hashlib
import re



@dataclass
class MethodConfig:
    name: str
    inputs : dict
    outputs: dict
    prerequisites : dict
    entity_root_type : str
    wdl: ParsedWDL
    wdl_sha256: str
    imports: Dict[str, ParsedWDL]


def _get_sha256(bytes):
    return hashlib.sha256(bytes).hexdigest()


def _find_configs(src_dir):
    config_files = []
    for dirpath, dirnames, filenames in os.walk(src_dir):
        for filename in filenames:
            if filename == "config.json":
                config_files.append(os.path.join(dirpath, filename))


def _read_config(src_dir, get_hosted_path):
    config_name = os.path.basename(src_dir)
    with open(os.path.join(src_dir, "method.wdl"), "rt") as fd:
        wdl = fd.read()
    with open(os.path.join(src_dir, "config.json"), "rt") as fd:
        config = json.load(fd)

    parsed = dumb_wdl_parse(wdl)

    imports = {}

    def get_input_name(filename):
        with open(filename, "rb") as fd:
            wdl = fd.read()

        # Currently doesn't support recursive imports at this time
        parsed = dumb_wdl_parse(wdl)
        assert len(parsed.get_imports()) == 0

        new_filename = os.path.join(_get_sha256(wdl), os.path.basename(filename))
        url = get_hosted_path(new_filename)
        imports[url] = parsed
        return url

    new_wdl = rewrite_with_relative_imports(parsed, get_input_name)
    return MethodConfig(
        name=config_name,
        config=config,
        wdl=new_wdl,
        imports=imports,
        wdl_sha256=_get_sha256(new_wdl.str()),
    )


def read_export(src_dir, get_hosted_path):
    config_files = _find_configs(os.path.join(src_dir, "configs"))

    configs = []
    for config_file in config_files:
        config = _read_config(os.path.dirname(config_file), get_hosted_path)
        configs.append(config)

    return configs



def _get_wdl_by_sha256(namespace, name):
    r = firecloud.api.list_repository_methods(namespace=namespace, name=name)
    assert r.status_code == 200
    for rec in r.json():
        m = re.match("SHA256:(\\S+)", rec.get("snapshotComment", ""))
        if m:
            sha256 = m.group(1)
            yield sha256, rec["snapshotId"]


class WDLCache:
    def __init__(self) -> None:
        self.scanned = set()
        self.wdl_lookup = {}  # keyed by (name, wdl_sha256)

    def get(self, method_namespace, method_name, wdl_sha256):
        method_key = (method_namespace, method_name)
        if method_key not in self.scanned:
            for sha256, version in _get_wdl_by_sha256(method_namespace, method_name):
                self.wdl_lookup[(method_namespace, method_name, sha256)] = version
            self.scanned.add(method_key)

import tempfile

def _sync_config(
    ws_namespace: str,
    ws_name: str,
    method_namespace: str,
    config: MethodConfig,
    cache: WDLCache,
):
    method_name = config.name
    method_version = cache.get(method_namespace, method_name, config.wdl_sha256)
    if method_version is None:
        synopsis = ""  # figure out where this should come from
        with tempfile.NamedTemporaryFile(mode="wt") as t:
            t.write(config.wdl.str())
            t.flush()
            method_version = _update_method(
                method_namespace, method_name, synopsis, t.name
            )

    _update_config(
        ws_namespace,
        ws_name,
        method_namespace,
        method_name,
        method_version,
        config.inputs,
        config.outputs,
        config.prerequisites,
        config.entity_root_type,
    )


def sync_workspace(ws_namespace, ws_name, method_namespace, export):
    cache = WDLCache()

    for config in export:
        _sync_config(method_namespace, config, cache)


def _update_method(namespace, name, synopsis, wdl):
    r = firecloud.api.update_repository_method(namespace, name, synopsis, wdl)
    assert r.status_code == 201
    return r.json()["snapshotId"]


def _update_config(
    ws_namespace,
    ws_name,
    method_namespace,
    method_name,
    method_version,
    inputs,
    outputs,
    prerequisites,
    entity_root_type,
):
    config = {
        "namespace": method_namespace,
        "name": method_name,
        "rootEntityType": entity_root_type,
        "methodRepoMethod": {
            "sourceRepo": "agora",
            "methodName": method_name,
            "methodNamespace": method_namespace,
            "methodVersion": method_version,
        },
        "methodConfigVersion": 2,
        "inputs": inputs,
        "outputs": outputs,
        "prerequisites": prerequisites,
        "deleted": False,
    }

    r = firecloud.api.create_workspace_config(
        ws_namespace, ws_name, config
    )
    assert r.status_code in [ 200, 201], f"status: {r.status_code}\nbody: {r.content}"
    # _update_method()


# update_config(config) method of dalmatian.wmanager.WorkspaceManager instance
#     Create or update a method configuration (separate API calls)

#     json_body = {
#        'namespace': config_namespace,
#        'name': config_name,
#        'rootEntityType' : entity,
#        'methodRepoMethod': {'methodName':method_name, 'methodNamespace':method_namespace, 'methodVersion':version},
#        'inputs':  {},
#        'outputs': {},
#        'prerequisites': {},
#        'deleted': False
#     }
