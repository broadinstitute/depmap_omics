import dalmatian
import os
import json
import requests
import urllib.parse

import re
import io
import firecloud.api

from dataclasses import dataclass
from typing import Union, Dict, List, Callable

IMPORT_PATTERN = re.compile('import\\s+"([^"]+)"\\s+as\\s+(\\S+)')


@dataclass
class TextSpan:
    text: str

    def str(self):
        return self.text


@dataclass
class Import:
    path: str
    name: str

    def str(self):
        return f'import "{self.path}" as {self.name}'


@dataclass
class ParsedWDL:
    statements: List[Union[TextSpan, Import]]

    def str(self):
        buffer = io.StringIO()
        for stmt in self.statements:
            buffer.write(stmt.str())
        return buffer.getvalue()


@dataclass
class MethodConfig:
    config: dict
    # synopsis : str
    wdl: ParsedWDL
    imports: Dict[str, ParsedWDL]


# def _write_config(dest_dir, config):

#     filename = os.path.join(dest_dir, "configs", config['namespace'], config['name'], "config.json")
#     parent = os.path.dirname(filename)
#     os.makedirs(parent, exist_ok=True)
#     with open(filename, "wt") as fd:
#         fd.write(json.dumps(min_config, indent=1, sort_keys=True))


def _get_wdl(method_config):
    method_repo_method = method_config["methodRepoMethod"]
    repo = method_repo_method.get("sourceRepo", "agora")
    if repo == "agora":
        wdl_name = f"agora:{method_repo_method['methodNamespace']}/{method_repo_method['methodName']}/{method_repo_method['methodVersion']}"
        print(f"Fetching {wdl_name}")
        wdl = _get_agora_method_wdl(
            method_repo_method["methodNamespace"],
            method_repo_method["methodName"],
            method_repo_method["methodVersion"],
        )
    else:
        assert repo == "dockstore"
        wdl_name = f"dockstore:{method_repo_method['methodPath']}/{method_repo_method['methodVersion']}"
        print(f"Fetching {wdl_name}")
        wdl = _fetch_dockstore_wdl(
            method_repo_method["methodPath"], method_repo_method["methodVersion"]
        )
    return wdl, wdl_name


def _get_all_configs(wm):
    configs = wm.list_configs()
    full_configs = []
    for config in configs:
        method_config = wm.get_config(f"{config['namespace']}/{config['name']}")
        imports = {}
        wdl, wdl_name = _get_wdl(method_config)
        parsed_wdl = dumb_wdl_parse(wdl)

        stack = []

        def fetch_import(path):
            print(f"Fetching {path} (imported from {stack[-1]})")
            resp = requests.get(path)
            name = os.path.basename(path)
            assert (
                name not in imports
            ), f"Found multiple files named {name} while traversing imports"
            parsed = dumb_wdl_parse(resp.content.decode("utf8"))
            stack.append(path)
            imports[name] = rewrite_with_relative_imports(parsed, fetch_import)
            return name

        stack.append(wdl_name)
        parsed_wdl = rewrite_with_relative_imports(parsed_wdl, fetch_import)
        stack.pop()
        full_configs.append(
            MethodConfig(config=method_config, wdl=parsed_wdl, imports=imports)
        )
    return full_configs


# def export_configs(wm, dest_dir):
#     full_configs = _get_all_configs(wm)

#     for config in full_configs:
#         _write_config(dest_dir, config)


def export_attributes(wm, dest_dir):
    attributes = wm.get_attributes()
    filename = os.path.join(dest_dir, "attributes.json")
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "wt") as fd:
        fd.write(json.dumps(attributes, indent=1))


def _get_agora_method(namespace, name, snapshot):
    r = firecloud.api.get_repository_method(namespace, name, snapshot)
    return r.json()


def _get_agora_method_wdl(namespace, name, snapshot):
    m = _get_agora_method(namespace, name, snapshot)
    return m["payload"]


def _write_dockstore_method(dest_dir, method):
    assert method["sourceRepo"] == "dockstore"
    wdl = _fetch_dockstore_wdl(method["methodPath"], method["methodVersion"])

    namespace, name = method["methodPath"].split("/")[-2:]

    filename = os.path.join(dest_dir, "methods", namespace, name, "method.json")
    parent = os.path.dirname(filename)
    os.makedirs(parent, exist_ok=True)
    with open(filename, "wt") as fd:
        fd.write(json.dumps(method, indent=1))
    filename = os.path.join(dest_dir, "methods", namespace, name, "method.wdl")
    with open(filename, "wt") as fd:
        fd.write(wdl)


def _write_agora_method(dest_dir, method):
    if "namespace" in method:
        namespace = method["namespace"]
    else:
        namespace = method["methodNamespace"]
    if "name" in method:
        name = method["name"]
    else:
        name = method["methodName"]


def _write_config(config: MethodConfig, dest_dir: str):
    method = config.config

    keys_to_keep = {"prerequisites", "rootEntityType", "inputs", "outputs"}
    keys_to_drop = {
        "methodConfigVersion",
        "methodRepoMethod",
        "name",
        "namespace",
        "deleted",
    }

    extra_keys = set(method.keys()).difference(keys_to_keep.union(keys_to_drop))
    assert len(extra_keys) == 0, f"Found extra unknown fields {extra_keys} in {config}"

    min_method = {k: method[k] for k in keys_to_keep}

    if "namespace" in method:
        namespace = method["namespace"]
    else:
        namespace = method["methodNamespace"]
    if "name" in method:
        name = method["name"]
    else:
        name = method["methodName"]

    filename = os.path.join(dest_dir, "configs", namespace, name, "config.json")
    parent = os.path.dirname(filename)
    os.makedirs(parent, exist_ok=True)
    with open(filename, "wt") as fd:
        fd.write(json.dumps(min_method, indent=1, sort_keys=True))

    filename = os.path.join(dest_dir, "configs", namespace, name, "method.wdl")
    with open(filename, "wt") as fd:
        fd.write(config.wdl.str())

    for imported_name, wdl in config.imports.items():
        filename = os.path.join(dest_dir, "configs", namespace, name, imported_name)
        parent = os.path.dirname(filename)
        os.makedirs(parent, exist_ok=True)
        with open(filename, "wt") as fd:
            fd.write(wdl.str())


def export_configs(wm, dest_dir):
    full_configs = _get_all_configs(wm)

    for config in full_configs:
        _write_config(config, dest_dir)


def _fetch_dockstore_wdl(method_path: str, method_version: int) -> ParsedWDL:
    response = requests.get(
        "https://dockstore.org/api/workflows/path/entry/{}/published".format(
            urllib.parse.quote_plus(method_path)
        )
    )
    x = response.json()

    repo_path = x["path"]
    assert repo_path.startswith("github.com/")
    repo_path = repo_path[len("github.com/") :]
    wdl_path = x["workflow_path"]
    assert wdl_path.startswith("/")
    wdl_path = wdl_path[1:]
    url = f"https://raw.githubusercontent.com/{repo_path}/{method_version}/{wdl_path}"

    response = requests.get(url)
    return response.content.decode("utf8")


def dumb_wdl_parse(wdl):
    matches = IMPORT_PATTERN.finditer(wdl)
    buffer = []
    prev = 0
    for match in matches:
        buffer.append(TextSpan(wdl[prev : match.start()]))
        path, name = match.groups()
        buffer.append(Import(path, name))
        prev = match.end()
    buffer.append(TextSpan(wdl[prev:]))
    return ParsedWDL(buffer)


def rewrite_with_relative_imports(
    parsed_wdl: ParsedWDL, get_import_path: Callable[[str], str]
):
    new_buffer = []

    for stmt in parsed_wdl.statements:
        if isinstance(stmt, Import):
            resp = requests.get(stmt.path)

            new_path = get_import_path(stmt.path)
            new_buffer.append(Import(new_path, stmt.name))
        else:
            new_buffer.append(stmt)

    return ParsedWDL(new_buffer)
