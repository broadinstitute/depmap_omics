from depmapomics import constants
from taigapy import create_taiga_client_v3
from taigapy.client_v3 import UploadedFile, LocalFormat

import depmapomics.patch_firecloud

depmapomics.patch_firecloud.install_patches()

print("creating brunello mutation matrix")
client = create_taiga_client_v3()
folder = constants.WORKING_DIR + "25Q2" + "/merged_"
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
            local_path=folder + "binary_germline_tkov3.csv",
            name="binary_mutation_tkov3",
            format=LocalFormat.CSV_TABLE,
            encoding="utf8",
        ),
    ],
)