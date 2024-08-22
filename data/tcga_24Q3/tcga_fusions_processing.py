from depmapomics.dm_omics import fusionPostProcessing
import pandas as pd
import numpy as np
from taigapy import TaigaClient
import dalmatian as dm
import pytest

tc = TaigaClient()

# pytest.set_trace()
fusionPostProcessing(
    refworkspace = 'broad-firecloud-ccle/tcga_rna_fusion_processing'
)