from depmapomics.test.config import (VIRTUAL_RELEASE, REFERENCE_RELEASE, FILE_ATTRIBUTES)
import pytest
import pandas as pd
from taigapy import TaigaClient

tc = TaigaClient()

def tcget_new_old(file):
    data1 = tc.get(name=REFERENCE_RELEASE['name'], file=file, version=REFERENCE_RELEASE['version'])
    data2 = tc.get(name=VIRTUAL_RELEASE['name'], file=file, version=VIRTUAL_RELEASE['version'])
    return data1, data2

@pytest.fixture(scope='module')
def data(request):
    return tcget_new_old(request.param)


PARAMS_nonmatrix_columns_match = [x['file'] for x in FILE_ATTRIBUTES if not x['ismatrix']]
@pytest.mark.parametrize('data', PARAMS_nonmatrix_columns_match, indirect=['data'])
def test_nonmatrix_columns_match(data):
    data1, data2 = data
    assert set(data1.columns) == set(data2.columns), \
        'data columns are not the same.\ncolumns added: {}\ncolumns dropped: {}'.\
            format(set(data2.columns) - set(data1.columns), set(data1.columns) - set(data2.columns))


PARAMS_matrix_correlations = [('CCLE_gene_cn', 0.95)]
PARAMS_matrix_correlations += [(x['file'], 0.99999) for x in FILE_ATTRIBUTES if x['ismatrix'] & (x['omicssource']=='RNA')]
@pytest.mark.parametrize('method', ['spearman', 'pearson'])
@pytest.mark.parametrize('axisname', ['pergene', 'persample'])
@pytest.mark.parametrize('data, threshold', PARAMS_matrix_correlations, indirect=['data'])
def test_matrix_correlations(data, threshold, axisname, method):
    axis = 0 if axisname == 'pergene' else 1 if axisname == 'persample' else None
    data1, data2 = data
    corrs = data1.corrwith(data2, axis=axis, drop=True, method=method)
    # TODO: test for NAs instead of dropping them. It seems like NA can happen if there are all zeros in one of the vectors, so let's dropna for now
    corrs.dropna(inplace=True)

    assert (corrs >= threshold).all(), 'the samples which did not pass the test are:\n{}'.format(corrs[corrs<threshold].sort_values())

# TODO: implement the assert almost equal version of the tests
# from pandas.testing import assert_frame_equal
# PARAMS_matrix_correlations = [('CCLE_expression_full', 0.01), ('CCLE_gene_cn', 0.01)]
# # @pytest.mark.parametrize('axisname', ['pergene', 'persample'])
# @pytest.mark.parametrize('data, rtol', PARAMS_matrix_correlations, indirect=['data'])
# def test_matrix_equality(data, rtol):
#     # axis = 0 if axisname == 'pergene' else 1 if axisname == 'persample' else None
#     data1, data2 = data
#     row = set(data1.index) & set(data2.index)
#     col = set(data1.columns) & set(data2.columns)
#     data1_ = data1.loc[row, col].T
#     data2_ = data2.loc[row, col].T
#     assert_frame_equal(data1_, data2_, rtol=rtol)

@pytest.fixture(scope='module')
def dataframes_merged(request):
    # TODO: figure out how to call the data fixture instead
    data1 = tc.get(name=REFERENCE_RELEASE['name'], file=request.param[0], version=REFERENCE_RELEASE['version'])
    data2 = tc.get(name=VIRTUAL_RELEASE['name'], file=request.param[0], version=VIRTUAL_RELEASE['version'])

    data_merged = pd.merge(data1, data2, on=request.param[1], indicator=True, how='outer')
    return data_merged

PARAMS_fraction_of_unequl_columns_from_merged_file = [(x['file'], x['merge_cols']) for x in FILE_ATTRIBUTES if 'merge_cols' in x]
@pytest.mark.parametrize('dataframes_merged', PARAMS_fraction_of_unequl_columns_from_merged_file,
                         indirect=['dataframes_merged'], ids=[x[0] for x in PARAMS_fraction_of_unequl_columns_from_merged_file])
# @pytest.mark.parametrize('force_dtype_to_str', [False, True], ids=['keep_dtype', 'force_str'])
def test_fraction_of_unequl_columns_from_merged_file(dataframes_merged, force_dtype_to_str = False):
    cols = list(set([x[:-2] for x in dataframes_merged.columns if x.endswith('_x') | x.endswith('_y')]))
    dataframe_merge_both = dataframes_merged[dataframes_merged['_merge'] == 'both']
    if force_dtype_to_str:
        dataframe_merge_both = dataframe_merge_both.astype(str)
    unequal_columns = pd.Series(index=cols, dtype=float)
    for col in cols:
        unequal_columns[col] = (dataframe_merge_both[col+'_x'] != dataframe_merge_both[col+'_y']).mean()
    unequal_columns.sort_values(ascending=False, inplace=True)

    unequal_columns = unequal_columns[unequal_columns > 0]
    assert len(unequal_columns) == 0, 'fraction of unequal columns when subsetted for shared columns {}:\n {}'.format('merge_columns', unequal_columns)
    return unequal_columns

PARAMS_compare_column_dtypes = [x['file'] for x in FILE_ATTRIBUTES if not x['ismatrix']]
@pytest.mark.parametrize('data', PARAMS_compare_column_dtypes, indirect=['data'])
def test_compare_column_dtypes(data):
    data1, data2 = data
    dtypes_compare = pd.concat([data1.dtypes, data2.dtypes], axis=1, keys=['reference', 'virtual'])
    dtypes_compare.dropna(inplace=True)
    dtypes_compare_nonmatching = dtypes_compare.apply(lambda x: x[0]!=x[1], axis=1)
    dtypes_compare = dtypes_compare[dtypes_compare_nonmatching]
    assert len(dtypes_compare) == 0, 'the following columns have changed types between the releases:\n{}'.format(dtypes_compare)
