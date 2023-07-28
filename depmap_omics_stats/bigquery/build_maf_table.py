"""bigquery demo from Phil."""
import dalmatian
from dataclasses import dataclass
import pandas as pd
import pickle
import os
import multiprocessing as mp

from google.cloud import bigquery
from tqdm import tqdm

from google.api_core.exceptions import BadRequest
import queue

namespace = "broad-firecloud-ccle"
workspaces = ["DepMap_WES_CN_hg38", "DepMap_WGS_CN"]
dest_dataset = "depmap-omics.maf_staging"


@dataclass
class Transfer:
    """A data class for transfer data state."""

    srcs: str
    dest_table: str
    cds_id: str


def get_transfers(workspace):
    """Get all transfers from Terra workspaces."""
    wm = dalmatian.WorkspaceManager(f"{namespace}/{workspace}")

    sample = wm.get_entities("sample")
    sample = sample.reset_index()

    transfers = []
    for rec in sample.to_dict("records"):
        if isinstance(rec['full_file'], list):
            dest_table = f"{dest_dataset}.stage_maf_{rec['sample_id'].replace('-', '_')}"
            transfers.append(Transfer(rec['full_file'], dest_table, rec["sample_id"]))

    return transfers


def create_ext_table(client,
                     srcs, dest_table, schema=None):
    """Create bigquery table from external parquet files.
    """
    table = bigquery.Table(dest_table)

    external_config = bigquery.ExternalConfig(
        bigquery.external_config.ExternalSourceFormat.PARQUET)
    external_config.source_uris = srcs
    external_config.schema = schema
    table.external_data_configuration = external_config

    client.create_table(table, exists_ok=True)
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
    uris = cleanup_uris(transfer.srcs) # at least one row has a entry which looks like a string containing a comma seperated list instead of a real list
    parq_table = create_ext_table(client, uris, transfer.dest_table, expected_schema)
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

        append_stmt = f"insert into {dest_table} select {selection} from {transfer.dest_table} where hugo_symbol != '' and hugo_symbol is not NULL"
        job = client.query(append_stmt)            
        try:
            job.result() # wait for completion
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
    expected_schema = dest_schema[1:] # everything except the first column

    while True:
        try:
            transfer = q.get(block=False)
        except queue.Empty:
            return
        result = do_transfer(client, transfer, dest_table, column_names, expected_schema)
        results.put(result)


def copy_in_parallel(transfers, parallism, dest_table):
    failures = []
    q = mp.Queue()
    results = mp.Queue()
    processes = [mp.Process(target=start_worker, args=(q,results, dest_table)) for i in range(parallism)]

    # fill up the queue with transfers
    for transfer in transfers:
        q.put(transfer)
    
    for p in processes:
        p.start()

    for _ in tqdm(transfers):
        result = results.get()
        if isinstance(result, Failure):
            print(result)
            print(result.query)
            failures.append(result)
        
    for p in processes:
        p.join()
    return failures


def create_batches(data, batch_size):
    return [data[x:x+batch_size] for x in range(0, len(data), batch_size)]

# Concatenate table adding CDS_ID to the dest table
def concatenate_tables(dest_table, transfers, parallelism):
    # create table for first transfer to copy schema from that one
    parq_table = create_ext_table(client, transfers[0].srcs, transfers[0].dest_table)

    create_table_stmt = f"create table if not exists {dest_table} as select 'invalid' CDS_ID, * from `{transfers[0].dest_table}` limit 0"
    job = client.query(create_table_stmt)
    job.result() # wait for completion
    client.delete_table(parq_table)

    # figure out which cds_ids have already been loaded
    already_loaded = set(pd.read_gbq(f"""select distinct cds_id from `{dest_dataset}.merged_maf` """)["cds_id"])
    
    # drop transfers already loaded
    remaining_transfers = [x for x in transfers if x.cds_id not in already_loaded]
    # remaining_transfers = remaining_transfers[:2]

    print(f"{len(already_loaded)} CDS IDs already loaded. {len(remaining_transfers)} of {len(transfers)} tables need to be loaded")

    failures = copy_in_parallel(remaining_transfers, parallism=parallelism, dest_table=dest_table)

    print(f"{len(failures)} failures")


if __name__ == "__main__":
    pickled_transfers = "transfers.pickled"
    if os.path.exists(pickled_transfers):
        print(f"Reading cached transfers from {pickled_transfers}. Delete this file if you want to get a fresh snapshot of all files from Terra")
        transfers = pickle.load(open(pickled_transfers,"rb"))
    else:
        transfers = []
        for workspace in workspaces:
            transfers.extend(get_transfers(workspace))
        with open(pickled_transfers, "wb") as fd:
            pickle.dump(transfers, fd)

    # Construct a BigQuery client object.
    client = bigquery.Client()

    concatenate_tables(f"{dest_dataset}.merged_maf", transfers, parallelism=8)
