import dalmatian as dm
import numpy as np
import pandas as pd


def main():
    wm = dm.WorkspaceManager("broad-firecloud-ccle/DEV_DepMap_hg38_RNAseq")
    attributes_df = wm.get_samples()
    strandness_df = pd.read_csv("/home/ubuntu/depmap_omics/ccle_tasks/rna_strandness_24q2.csv", index_col=0)

    overlapped_ids = np.intersect1d(strandness_df.index, attributes_df.index)
    attributes_df.loc[overlapped_ids, 'strandness'] = strandness_df.loc[overlapped_ids, :].iloc[:, -1].values
    attributes_df = attributes_df.loc[~attributes_df.strandness.isnull(), :]

    print(attributes_df.shape)
    print(attributes_df.loc[attributes_df.loc[:, 'strandness'], :])
    wm.update_sample_set('test_strandness50', attributes_df.loc[attributes_df.loc[:, 'strandness'], :].index[:50])


main()
