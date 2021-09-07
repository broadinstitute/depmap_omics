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
    # no need to test plotting
    # called by postProcess(). called by CCLEPostProcessing()
    # CCLEPostProcessing only cares about wgs and wes?
    # call params directly from config or define new globals?
    # add dummy columns to ws to be removed later
    dummy_cols = ['a', 'b', 'c']

    

    output = loadFromGATKAggregation(refworkspace,  sortby=[SAMPLEID, 'Chromosome', "Start", "End"], 
                            save_output='', doCleanup=False,
                            todrop=[], colname="combined_seg_file",
                            plotColname="modeled_segments_plot_tumor", tempFolder="temp/",
                            toremove=["readgroup_ubams", ],
                            sampleset="all", colRenaming=colrenaming_sandbox)

    assert(isinstance(output, pd.DataFrame), "Output is not a dataframe")
    # if colremoves are removed
    # if headers are renamed correctly
    # type of each col
    # (how) is output indexed
    

