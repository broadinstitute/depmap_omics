import seaborn as sns
import matplotlib.pyplot as plt
from taigapy import TaigaClient
import pandas as pd
from depmapomics.test.config import (VIRTUAL_RELEASE, REFERENCE_RELEASE)

tc = TaigaClient()

def create_data_stack(file, number_of_points=500, random_state=0,
                      release1 = REFERENCE_RELEASE, release2 = VIRTUAL_RELEASE):

    data1 = tc.get(file=file, **release1)
    data2 = tc.get(file=file, **release2)

    row = set(data1.index) & set(data2.index)
    col = set(data1.columns) & set(data2.columns)

    data1_stack = data1.loc[row, col].stack()
    data2_stack = data2.loc[row, col].stack()

    if number_of_points > 0:
        data1_stack = data1_stack.sample(number_of_points, random_state=random_state)
    data2_stack = data2_stack.loc[data1_stack.index]

    data_stack = pd.concat([data1_stack, data2_stack], axis=1)

    cols = ['{:s}.{:d}'.format(release1['name'], release1['version']),
            '{:s}.{:d}'.format(release2['name'], release2['version'])]
    data_stack.columns = cols
    data_stack.reset_index(inplace=True)
    data_stack.rename(columns={'level_0': 'DepMap_ID', 'level_1': 'gene'}, inplace=True)
    return data_stack, cols


def CCLE_gene_cn_with_source_change(number_of_points = 500, savefig = False,
                                    release1 = REFERENCE_RELEASE, release2 = VIRTUAL_RELEASE):
    names = [release1['name'], release2['name']]
    CCLE_gene_cn_12, cols = create_data_stack('CCLE_gene_cn', number_of_points=number_of_points, release1 = release1, release2 = release2)

    CCLE_segment_cn_1 = tc.get(name=release1['name'], version=release1['version'], file='CCLE_segment_cn')
    CCLE_segment_cn_2 = tc.get(name=release2['name'], version=release2['version'], file='CCLE_segment_cn')
    sources = pd.merge(CCLE_segment_cn_1[['DepMap_ID', 'Source']].drop_duplicates(),
         CCLE_segment_cn_2[['DepMap_ID', 'Source']].drop_duplicates(),
         on='DepMap_ID', suffixes=['_'+names[0], '_'+names[1]])
    del CCLE_segment_cn_1, CCLE_segment_cn_2
    sources['source_change'] = sources.apply(lambda x: '{:s} -> {:s}'.format(x['Source_'+names[0]], x['Source_'+names[1]]), axis=1)
    sources['source_has_changed'] = (sources['Source_'+names[0]] != sources['Source_'+names[1]])

    CCLE_gene_cn_12 = pd.merge(CCLE_gene_cn_12, sources, on='DepMap_ID')
    return CCLE_gene_cn_12, cols


def plot_matrix_comparison(file, number_of_points = 500,
                           release1 = REFERENCE_RELEASE, release2 = VIRTUAL_RELEASE):
    data_stack, cols = create_data_stack(file, number_of_points = number_of_points,
                                         release1 = release1, release2 = release2)
    data_stack.plot.scatter(*cols)


def plot_gene_cn_comparison(number_of_points = 500, savefig = False,
                            release1 = REFERENCE_RELEASE, release2 = VIRTUAL_RELEASE):

    CCLE_gene_cn_12, cols = CCLE_gene_cn_with_source_change(number_of_points = number_of_points, savefig = False,
                                                            release1 = release1, release2 = release2)
    # CCLE_gene_cn_12, cols = create_data_stack('CCLE_gene_cn', number_of_points=500, release1 = release1, release2 = release2)

    # CCLE_segment_cn_1 = tc.get(name=release1['name'], version=release1['version'], file='CCLE_segment_cn')
    # CCLE_segment_cn_2 = tc.get(name=release2['name'], version=release2['version'], file='CCLE_segment_cn')
    # sources = pd.merge(CCLE_segment_cn_1[['DepMap_ID', 'Source']].drop_duplicates(),
    #      CCLE_segment_cn_2[['DepMap_ID', 'Source']].drop_duplicates(),
    #      on='DepMap_ID', suffixes=['_'+names[0], '_'+names[1]])
    # del CCLE_segment_cn_1, CCLE_segment_cn_2
    # sources['source_change'] = sources.apply(lambda x: '{:s} -> {:s}'.format(x['Source_'+names[0]], x['Source_'+names[1]]), axis=1)
    # sources['source_has_changed'] = (sources['Source_'+names[0]] != sources['Source_'+names[1]])

    # CCLE_gene_cn_12 = pd.merge(CCLE_gene_cn_12, sources, on='DepMap_ID')

    plt.figure(figsize=(20,10))
    sns.scatterplot(data=CCLE_gene_cn_12.sample(number_of_dots, random_state=0), x=cols[0], y=cols[1],
                    hue='source_change', style='source_has_changed', alpha=0.5, cmap='Tab20')

    if savefig:
        filename = ('{}-vs-{}.png'.format(cols[1], cols[0]))#.replace('/', '.')
        print('saving {}'.format(filename))
        plt.savefig(filename, bbox_inches='tight')

    return CCLE_gene_cn_12
