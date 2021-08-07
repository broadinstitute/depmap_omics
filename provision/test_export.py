import dalmatian
from dump_workspace import export_attributes, export_configs,_fetch_dockstore_wdl

def test_export_configs():
    wm = dalmatian.WorkspaceManager("broad-firecloud-ccle/DepMap_hg38_RNAseq")
    export_configs(wm, "export")


def test_export_attributes():
    wm = dalmatian.WorkspaceManager("broad-firecloud-ccle/DepMap_hg38_RNAseq")
    export_attributes(wm, "export")

def test_dockstore():
    wdl = _fetch_dockstore_wdl('github.com/broadinstitute/depmap_omics/RNA_aggregate', 'master')
    print(wdl)

