from pandas.core.frame import DataFrame
from depmapomics.utils import generateGeneNames
import pandas as pd
import depmapomics.utils
from depmapomics.config import ENSEMBL_SERVER_V

def test_generateGeneNames(tmpdir, monkeypatch):
    def mock(ensemble_server, attributes):
        assert(ENSEMBL_SERVER_V == "http://nov2020.archive.ensembl.org/biomart", "Wrong local biomart version for tests")
        return pd.read_csv("~/.depmapomics/biomart_ensembltohgnc_nov2020.csv")
    monkeypatch.setattr(depmapomics.utils, "_fetchFromServer", mock)
    attr = ['start_position', 'end_position', "chromosome_name"]
    output = generateGeneNames(useCache=True, cache_folder=str(tmpdir) + "/",
      attributes=attr)
    assert(isinstance(output, pd.DataFrame), "Output is not a dataframe")
    assert(len(output.columns) == len(attr) + 5, "Wrong number of columns")

