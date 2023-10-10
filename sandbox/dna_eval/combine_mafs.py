""" MAF aggregation scripts
"""

import numpy as np
import json
import uuid

import dalmatian
from dataclasses import dataclass
import pandas as pd
import pickle
import os
import glob
import multiprocessing as mp
import argparse

from google.cloud import bigquery
from tqdm import tqdm

from google.api_core.exceptions import BadRequest
import queue

namespace = "broad-firecloud-ccle"
workspaces = ["DepMap_WES_CN_hg38", "DepMap_WGS_CN"]


@dataclass
class Transfer:
    """A data class for transfer data state."""

    srcs: str
    cds_id: str


def get_transfers(workspace, dest_dataset='', ids=None):
    """Get all transfers from Terra workspaces."""
    wm = dalmatian.WorkspaceManager(f"{namespace}/{workspace}")

    sample = wm.get_entities("sample")
    sample = sample.reset_index()

    transfers = []
    i = 0
    for rec in sample.to_dict("records"):
        if 'depmap_maf_23q4' in rec:
            if not isinstance(rec['depmap_maf_23q4'], float):
                if ids is None:
                    transfers.append(Transfer(rec['depmap_maf_23q4'], rec["sample_id"]))
                elif rec["sample_id"] in ids:
                    transfers.append(Transfer(rec['depmap_maf_23q4'], rec["sample_id"]))
    return transfers


def create_ext_table(client, srcs, dest_table, schema, column_types: dict):
    """Create bigquery table from external parquet files.
    """
    table = bigquery.Table(dest_table)

    external_config = bigquery.ExternalConfig(
        bigquery.external_config.ExternalSourceFormat.PARQUET)
    external_config.source_uris = srcs
    external_config.schema = schema
    table.external_data_configuration = external_config

    client.create_table(table, exists_ok=True)
    if schema is None and column_types is not None:
        # bigqueries autodetection of schema is not reliable. In particular if the first parquet file in srcs has
        # rows with all null values, it will guess the type is integer which will cause a problem if the a later
        # file contains rows with strings as the values. To avoid this problem, fetch the schema from the table
        # to find out the column names, but then look up those column names in our map where we've defined the
        # column types instead of guessing the column types
        new_schema = []
        fetched_table = client.get_table(dest_table)
        for field in fetched_table.schema:
            field_type = column_types[field.name]
            # if field_type != field.field_type:
            #     print(
            #         f"overriding type on {field.name}: {field.field_type} -> {field_type}"
            #     )
            new_schema.append(bigquery.SchemaField(field.name, field_type, "NULLABLE"))

        # now recreate the ext table with this new schema
        client.delete_table(table)
        table = create_ext_table(client, srcs, dest_table, new_schema, column_types)

    return table


def cleanup_uris(uris):
    result = []
    for uri in uris:
        result.extend([x.strip() for x in uri.split(",")])
    return list(set(result))

@dataclass
class Failure:
    transfer: Transfer
    exception: str
    query: str
    
@dataclass
class Success:
    transfer: Transfer

def do_transfer(transfer):
    try:
        df = pd.read_csv(transfer.srcs, index_col=0)
        df['cds_id'] = transfer.cds_id
        assert df.index.name == '0915'
    except BadRequest as ex:
        return Failure(transfer, str(ex), transfer.cds_id)
    return Success(df)

def start_worker(q, results, dest_table):
    while True:
        try:
            transfer = q.get(block=False)
        except queue.Empty:
            return
        result = do_transfer(
            transfer
        )
        results.put(result)


def read_in_parallel(transfers, parallism, dest_table):
    failures = []
    q = mp.Queue()
    results = mp.Queue()
    processes = [mp.Process(target=start_worker, args=(q, results, dest_table)) for i in range(parallism)]

    # fill up the queue with transfers
    for transfer in transfers:
        q.put(transfer)
    
    for p in processes:
        p.start()

    success_dfs = []
    for _ in tqdm(transfers):
        result = results.get()
        if isinstance(result, Failure):
            failures.append(result)
        else:
            success_dfs.append(result)
        
    for p in processes:
        p.join()
    return failures, success_dfs


def unify_schema(a, b):
    if a == b:
        print("No change")
        return

    # compare common fields
    fields = set(a).intersection(b)
    for field in fields:
        field_type = a[field]
        if a[field] != b[field]:
            # this is often the case of it guesses INTEGER but the real type is something else
            if a[field] == "INTEGER":
                field_type = b[field]
            print(f"{field}: {a[field]} != {b[field]}. (Picking {field_type})")
            a[field] = field_type
    # add unique fields in b into a
    for field in set(b).difference(a):
        assert field not in a
        a[field] = b[field]
    # print(f"oc_brca1_func_assay__class -> {a.get('oc_brca1_func_assay__class')}")


def unify_schemas(schemas: list, overrides: dict):
    # return a map of field name -> field type. First definition of field type wins, but disagreements are
    # printed and one can add it to the provided overrides dict to decide what it should be
    unified = dict(overrides)
    for schema in schemas:
        unify_schema(unified, {field.name: field.field_type for field in schema})
    return unified


def get_schema_from_parquet(client, temp_table_name, filename):
    table = bigquery.Table(temp_table_name)

    external_config = bigquery.ExternalConfig(
        bigquery.external_config.ExternalSourceFormat.PARQUET
    )
    external_config.source_uris = filename
    table.external_data_configuration = external_config

    print(f"creating table {temp_table_name} for {filename}")
    client.create_table(table, exists_ok=True)
    # read it back, and now schema will be populated
    fetched_table = client.get_table(temp_table_name)
    client.delete_table(table)

    return fetched_table.schema


def dump_schema(client, transfers):
    def get_schemas():
        for transfer in transfers:
            print(f"Reading files for {transfer.dest_table}...")
            for src in transfer.srcs:
                schema = get_schema_from_parquet(
                    client, transfer.dest_table + "_tmp_pscan", src
                )
                yield schema

    column_types = unify_schemas(
        get_schemas(),
        {},
    )

    filename = "column_types.json"
    print(f"dumping column types to {filename}")
    with open(filename, "wt") as fd:
        fd.write(json.dumps(column_types, indent=2, sort_keys=True))


# Concatenate table adding CDS_ID to the dest table
def concatenate_tables(dest_table, client, transfers, parallelism=10):
    # create table for first transfer to copy schema from that one
    failures, success_df = read_in_parallel(
        transfers=transfers, parallism=parallelism, dest_table=dest_table
    )
    print(failures)
    print(len(success_df))
    print(f"{len(failures)} failures")
    return success_df


def main():
    pickled_transfers = "transfers.pickled"
    ids = None
    #ids = ['CDS-00Nrci', 'CDS-099jzP', 'CDS-00rz9N', 'CDS-01bI6z']

    if os.path.exists(pickled_transfers):
        print(f"Reading cached transfers from {pickled_transfers}. Delete this file if you want to get a fresh snapshot of all files from Terra")
        transfers = pickle.load(open(pickled_transfers,"rb"))
    else:
        transfers = []
        for workspace in workspaces:
            transfers.extend(get_transfers(workspace, ids=ids))
        with open(pickled_transfers, "wb") as fd:
            pickle.dump(transfers, fd)
    success_df = concatenate_tables(f"mafs_latest", None, transfers, parallelism=15)
    df = pd.concat([s.transfer for s in success_df], axis=0)
    df.to_csv("23Q4_mutation_maf_latest.tsv", sep='\t')


if __name__ == "__main__":
    main()
