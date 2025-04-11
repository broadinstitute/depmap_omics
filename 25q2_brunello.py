import dalmatian as dm
import numpy as np
import pandas as pd
from depmapomics import constants, mutations
from taigapy import create_taiga_client_v3
from taigapy.client_v3 import UploadedFile, LocalFormat

print("creating brunello mutation matrix")
client = create_taiga_client_v3()

# read cds->pr mapping table and construct renaming dictionary
# always read latest version
print("reading seq -> pr ID mapping table locally")
omics_id_mapping_table = pd.read_csv("data/25Q2/omics_profile_to_sequencing_id.csv")
renaming_dict = dict(
    list(
        zip(
            omics_id_mapping_table["sequencing_id"],
            omics_id_mapping_table["profile_id"],
        )
    )
)

wes_wm = dm.WorkspaceManager("broad-firecloud-ccle/DepMap_WES_CN_hg38")
wgs_wm = dm.WorkspaceManager("broad-firecloud-ccle/DepMap_WGS_CN")

folder = constants.WORKING_DIR + "25Q2" + "/merged_"

print("aggregating binary guide mutation matrices")
print("aggregating wes")
wes_germline_mats = mutations.aggregateGermlineMatrix(
    wes_wm, "all", save_output=folder+"wes_"
)
print("aggregating wgs")
wgs_germline_mats = mutations.aggregateGermlineMatrix(
    wgs_wm, "all", save_output=folder+"wgs_"
)

for lib, _ in constants.GUIDESBED.items():
    assert lib in wes_germline_mats, "library missing in wes"
    assert lib in wgs_germline_mats, "library missing in wgs"
    # merging wes and wgs
    print("renaming merged wes and wgs germline matrix for library: ", lib)
    germline_mat_merged = pd.concat(
        [wes_germline_mats[lib], wgs_germline_mats[lib].iloc[:, 4:]], axis=1
    )
    germline_mat_merged_noguides = germline_mat_merged.iloc[:, 4:]

    # transform from CDSID-level to PR-level
    whitelist_cols = [
        x for x in germline_mat_merged_noguides.columns if x in renaming_dict
    ]
    whitelist_germline_mat = germline_mat_merged_noguides[whitelist_cols]
    mergedmat = whitelist_germline_mat.rename(columns=renaming_dict)

    mergedmat = mergedmat.astype(bool).astype(int)
    sorted_mat = germline_mat_merged.iloc[:, :4].join(mergedmat)
    sorted_mat["end"] = sorted_mat["end"].astype(int)
    print("saving merged binary matrix for library: ", lib)
    sorted_mat.to_csv(
        folder + "binary_germline" + "_" + lib + ".csv", index=False
    )
    print("done")

client.update_dataset(
    reason="new " + "25Q2" + " release!",
    permaname="mutations-latest-ed72",
    additions=[
        UploadedFile(
            local_path=folder + "binary_germline_avana.csv",
            name="binary_mutation_avana",
            format=LocalFormat.CSV_TABLE,
            encoding="utf8",
        ),
        UploadedFile(
            local_path=folder + "binary_germline_ky.csv",
            name="binary_mutation_ky",
            format=LocalFormat.CSV_TABLE,
            encoding="utf8",
        ),
        UploadedFile(
            local_path=folder + "binary_germline_humagne.csv",
            name="binary_mutation_humagne",
            format=LocalFormat.CSV_TABLE,
            encoding="utf8",
        ),
        UploadedFile(
            local_path=folder + "binary_germline_brunello.csv",
            name="binary_mutation_brunello",
            format=LocalFormat.CSV_TABLE,
            encoding="utf8",
        ),
        UploadedFile(
            local_path=folder + "binary_germline_tkov3.csv",
            name="binary_mutation_tkov3",
            format=LocalFormat.CSV_TABLE,
            encoding="utf8",
        ),
    ],
)