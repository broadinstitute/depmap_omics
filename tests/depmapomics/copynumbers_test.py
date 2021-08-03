import numpy as np
from gsheets import Sheets
from depmapomics import tracker as track
from depmapomics import utils
from depmapomics.qc import cn
from depmapomics.config import *
from IPython.display import Image, display
import dalmatian as dm
import pandas as pd
from taigapy import TaigaClient
tc = TaigaClient()
import os
from genepy import mutations as mut
from genepy.utils import helper as h
from genepy import terra
from genepy.google import gcp
from genepy import rna
import matplotlib.pyplot as plt

def test_loadFromGATKAggregation(tmpdir, monkeypatch):
    loadFromGATKAggregation(refworkspace,  sortby=[SAMPLEID, 'Chromosome', "Start", "End"], 
                            save_output='', doCleanup=True,
                            todrop=[], showPlots=False, colname="combined_seg_file",
                            plotColname="modeled_segments_plot_tumor", tempFolder="temp/",
                            toremove=["readgroup_ubams", ],
                            sampleset="all", colRenaming=COLRENAMING)
