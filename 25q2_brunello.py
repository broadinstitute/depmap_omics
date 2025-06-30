import dalmatian as dm
import numpy as np
import pandas as pd
from depmapomics import constants, mutations
from taigapy import create_taiga_client_v3
from taigapy.client_v3 import UploadedFile, LocalFormat

import depmapomics.patch_firecloud

depmapomics.patch_firecloud.install_patches()

print("creating brunello mutation matrix")
client = create_taiga_client_v3()

# read cds->pr mapping table and construct renaming dictionary
# always read latest version
# print("reading seq -> pr ID mapping table locally")
# omics_id_mapping_table = pd.read_csv("data/25Q2/omics_profile_to_sequencing_id.csv")

print("reading omics ID mapping table from taiga")
omics_id_mapping_table = client.get(
    name="2025-05-01-master-mapping-table-28c2", file="public_release_date.2025-05-01.master_mapping_table"
)

omics_id_mapping_table = omics_id_mapping_table[omics_id_mapping_table.is_default_entry == "True"]
assert len(omics_id_mapping_table) > 0
renaming_dict = dict(
    list(
        zip(
            omics_id_mapping_table["sequencing_id"],
            omics_id_mapping_table["model_id"],
        )
    )
)

wes_wm = dm.WorkspaceManager("broad-firecloud-ccle/DepMap_WES_CN_hg38")
wgs_wm = dm.WorkspaceManager("broad-firecloud-ccle/DepMap_WGS_CN")

folder = constants.WORKING_DIR + "25Q2" + "/merged_"

print("aggregating binary guide mutation matrices")
print("aggregating wes")
wes_seqids = omics_id_mapping_table[omics_id_mapping_table.datatype == "wes"]['sequencing_id'].tolist()
print(len(wes_seqids), "default WES seq ids")
wes_germline_mats = mutations.aggregateGermlineMatrix(
    wes_wm, samples_in_set=wes_seqids, save_output=folder+"wes_"
)
print("aggregating wgs")
wgs_seqids = omics_id_mapping_table[omics_id_mapping_table.datatype == "wgs"]['sequencing_id'].tolist()
print(len(wgs_seqids), "default WGS seq ids")
wgs_germline_mats = mutations.aggregateGermlineMatrix(
    wgs_wm, samples_in_set=wgs_seqids, save_output=folder+"wgs_"
)

for lib, _ in constants.GUIDESBED.items():
    assert lib in wes_germline_mats, "library missing in wes"
    assert lib in wgs_germline_mats, "library missing in wgs"
    # merging wes and wgs
    print("renaming merged wes and wgs germline matrix for library: ", lib)
    germline_mat_merged = pd.concat(
        [wes_germline_mats[lib], wgs_germline_mats[lib].iloc[:, 4:]], axis=1
    )
    guildes_only = germline_mat_merged.iloc[:, :4]
    germline_mat_merged = germline_mat_merged.iloc[:, 4:]
    germline_mat_merged = germline_mat_merged.astype(bool).astype(int)

    # transform from CDSID-level to PR-level
    whitelist_cols = [
        x for x in germline_mat_merged.columns if x in renaming_dict
    ]
    germline_mat_merged = germline_mat_merged[whitelist_cols].rename(columns=renaming_dict)

    germline_mat_merged = guildes_only.join(germline_mat_merged)
    germline_mat_merged["end"] = germline_mat_merged["end"].astype(int)
    print("saving merged binary matrix for library: ", lib)
    germline_mat_merged.to_csv(
        folder + "binary_germline" + "_" + lib + ".csv", index=False
    )
    print("done")

client.update_dataset(
    reason="new " + "25Q2" + " release!",
    permaname="mutations-latest-ed72",
    additions=[
        UploadedFile(
            local_path=folder + "binary_germline_brunello.csv",
            name="binary_mutation_brunello_model_public",
            format=LocalFormat.CSV_TABLE,
            encoding="utf8",
        ),
    ],
)