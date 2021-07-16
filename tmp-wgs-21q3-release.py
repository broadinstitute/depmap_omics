from depmapomics.config import SAMPLESETNAME
from depmapomics import copynumbers as omics_cn
import pandas as pd


# wespriosegs, wgspriosegs = omics_cn.CCLEPostProcessing(samplesetname=SAMPLESETNAME, redoWES=False)

wgspriosegs = pd.read_csv("temp/21Q3/wgs_segments_all_latest.csv")
wespriosegs = pd.read_csv("temp/21Q3/wes_segments_all_latest.csv")


omics_cn.ProcessForAchilles(wespriosegs, wgspriosegs, samplesetname=SAMPLESETNAME)
