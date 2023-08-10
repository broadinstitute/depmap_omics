""" BigQuery batch search script for DepMap mutation data
"""

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
# workspaces = ["DepMap_WES_CN_hg38", "DepMap_WGS_CN"]
workspaces = ["DEV_DepMap_WES_CN_hg38", "DEV_DepMap_WGS_CN"]

import json

with open("column_types.json", "rt") as fd:
    COLUMN_TYPES = json.load(fd)


@dataclass
class Transfer:
    """A data class for transfer data state."""

    srcs: str
    dest_table: str
    cds_id: str


def get_transfers(workspace, dest_dataset="", ids=None):
    """Get all transfers from Terra workspaces."""
    wm = dalmatian.WorkspaceManager(f"{namespace}/{workspace}")

    sample = wm.get_entities("sample")
    sample = sample.reset_index()
    print(sample.head())
    print(sample.loc[:, "sample_id"].isin(ids).sum())

    transfers = []
    for rec in sample.to_dict("records"):
        if isinstance(rec["parquet_with_hgvs"], list):
            dest_table = f"{dest_dataset}.{rec['sample_id'].replace('-', '_')}".lower()
            if ids is None:
                transfers.append(
                    Transfer(rec["full_file"], dest_table, rec["sample_id"])
                )
            elif isinstance(ids, list) and rec["sample_id"] in ids:
                transfers.append(
                    Transfer(rec["parquet_with_hgvs"], dest_table, rec["sample_id"])
                )
    return transfers


def create_ext_table(client, srcs, dest_table, schema, column_types: dict):
    """Create bigquery table from external parquet files."""
    table = bigquery.Table(dest_table)

    external_config = bigquery.ExternalConfig(
        bigquery.external_config.ExternalSourceFormat.PARQUET
    )
    external_config.source_uris = srcs
    external_config.schema = schema
    table.external_data_configuration = external_config

    # print(f"creating table {dest_table}")
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


def do_transfer(client, transfer, dest_table, dest_column_names, expected_schema):
    uris = cleanup_uris(
        transfer.srcs
    )  # at least one row has a entry which looks like a string containing a comma seperated list instead of a real list
    parq_table = create_ext_table(
        client, uris, transfer.dest_table, expected_schema, COLUMN_TYPES
    )
    try:
        column_names = [x.name for x in client.get_table(transfer.dest_table).schema]

        assert dest_column_names[0] == "CDS_ID"
        dest_column_names = dest_column_names[1:]

        # assert "civic_description" in dest_column_names
        # # hack: civic_id is of type INTEGER in some parquet files. Add a cast to prevent an exception
        # dest_column_names = [x if x != 'civic_description' else 'CAST(civic_description AS STRING)' for x in dest_column_names]

        selection = ", ".join(dest_column_names)

        if "CDS_ID" not in column_names:
            selection = f"'{transfer.cds_id}' CDS_ID, {selection}"

        append_stmt = f"insert into {dest_table} select {selection} from {transfer.dest_table}"  # where hugo_symbol != '' and hugo_symbol is not NULL"
        job = client.query(append_stmt)
        try:
            job.result()  # wait for completion
        except BadRequest as ex:
            return Failure(transfer, str(ex), append_stmt)
    finally:
        # clean up temp external table
        client.delete_table(parq_table)
    return Success(transfer)


def start_worker(q, results, dest_table):
    client = bigquery.Client()
    dest_schema = client.get_table(dest_table).schema
    column_names = [x.name for x in dest_schema]
    expected_schema = dest_schema[1:]  # everything except the first column

    while True:
        try:
            transfer = q.get(block=False)
        except queue.Empty:
            return
        result = do_transfer(
            client, transfer, dest_table, column_names, expected_schema
        )
        results.put(result)


def copy_in_parallel(transfers, parallism, dest_table):
    failures = []
    q = mp.Queue()
    results = mp.Queue()
    processes = [
        mp.Process(target=start_worker, args=(q, results, dest_table))
        for i in range(parallism)
    ]

    # fill up the queue with transfers
    for transfer in transfers:
        q.put(transfer)

    for p in processes:
        p.start()

    for _ in tqdm(transfers):
        result = results.get()
        if isinstance(result, Failure):
            breakpoint()
            print(result)
            print(result.query)
            failures.append(result)

    for p in processes:
        p.join()
    return failures


def create_batches(data, batch_size):
    return [data[x : x + batch_size] for x in range(0, len(data), batch_size)]


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


import json
import uuid


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
def concatenate_tables(dest_table, client, transfers, parallelism):
    # create table for first transfer to copy schema from that one
    parq_table = create_ext_table(
        client, transfers[0].srcs, transfers[0].dest_table, None, COLUMN_TYPES
    )

    create_table_stmt = f"create table if not exists {dest_table} as select 'invalid' CDS_ID, * from `{transfers[0].dest_table}` limit 0"
    job = client.query(create_table_stmt)
    job.result()  # wait for completion
    client.delete_table(parq_table)

    # check to make sure the schmea has the right column types. (If we have an existing table which was not created from the line above, the columns might be wrong)
    dest_table_schema = client.get_table(dest_table).schema
    mismatches = []
    all_column_types = dict(COLUMN_TYPES)
    all_column_types["CDS_ID"] = "STRING"
    for field in dest_table_schema:
        expected_type = all_column_types.get(field.name)
        if field.field_type != expected_type:
            mismatches.append(f"{field.name}: {field.field_type} != {expected_type}")
    if len(mismatches) > 0:
        raise Exception(
            f"Looks like {dest_table} was created with an older version of COLUMN_TYPES. Drop table and re-run to recreate table. (Mismatched types: {mismatches})"
        )

    # figure out which cds_ids have already been loaded
    already_loaded = set(
        pd.read_gbq(f"""select distinct cds_id from `{dest_table}` """)["cds_id"]
    )

    # drop transfers already loaded
    remaining_transfers = [x for x in transfers if x.cds_id not in already_loaded]
    # remaining_transfers = remaining_transfers[:2]

    print(
        f"{len(already_loaded)} CDS IDs already loaded. {len(remaining_transfers)} of {len(transfers)} tables need to be loaded"
    )

    failures = copy_in_parallel(
        remaining_transfers, parallism=parallelism, dest_table=dest_table
    )
    print(failures)

    print(f"{len(failures)} failures")


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ids",
        type=str,
        required=True,
        help="CDS IDs to constraint the data query",
    )
    parser.add_argument(
        "--dest",
        default="depmap-omics.maf_staging",
        required=False,
        help="use workspace stored parquet file or local parquet file",
    )
    parser.add_argument(
        "--dumpschema",
        action="store_true",
        help="If set will write out the inferred schema and exit",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    print(args)
    ids = args.ids.split(",")
    print(ids)

    pickled_transfers = "transfers.pickled"
    if os.path.exists(pickled_transfers):
        print(
            f"Reading cached transfers from {pickled_transfers}. Delete this file if you want to get a fresh snapshot of all files from Terra"
        )
        transfers = pickle.load(open(pickled_transfers, "rb"))
    else:
        transfers = []
        for workspace in workspaces:
            print(workspace)
            transfers.extend(get_transfers(workspace, args.dest, ids))
        with open(pickled_transfers, "wb") as fd:
            pickle.dump(transfers, fd)
    print(transfers)

    # Construct a BigQuery client object.
    client = bigquery.Client()

    if args.dumpschema:
        dump_schema(client, transfers)
    else:
        concatenate_tables(
            f"{args.dest}.merged_maf_new", client, transfers, parallelism=8
        )


if __name__ == "__main__":
    main()
