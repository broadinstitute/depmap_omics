from depmapomics.dm_omics import mutationPostProcessing
import pandas as pd
import numpy as np
from taigapy import TaigaClient
import dalmatian as dm
import pytest

tc = TaigaClient()

# pytest.set_trace()
mutationPostProcessing(
    wgsrefworkspace = 'broad-firecloud-ccle/tcga_mutation_testing'
)