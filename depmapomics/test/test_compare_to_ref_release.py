from depmapomics.test.config import (FILE_ATTRIBUTES, VIRTUAL_RELEASE, REFERENCE_RELEASE)
import pytest
from taigapy import TaigaClient

tc = TaigaClient()

@pytest.fixture(scope='module')
def data(request):
    data1 = tc.get(name=REFERENCE_RELEASE['name'], file=request.param, version=REFERENCE_RELEASE['version'])
    data2 = tc.get(name=VIRTUAL_RELEASE['name'], file=request.param, version=VIRTUAL_RELEASE['version'])
    return data1, data2


PARAMS_matrix_correlations = [('CCLE_expression_full', 0.95), ('CCLE_gene_cn', 0.95)]
@pytest.mark.parametrize('method', ['spearman', 'pearson'])
@pytest.mark.parametrize('axisname', ['pergene', 'persample'])
@pytest.mark.parametrize('data, threshold', PARAMS_matrix_correlations, indirect=['data'])
def test_matrix_correlations(data, threshold, axisname, method):
    axis = 0 if axisname == 'pergene' else 1 if axisname == 'persample' else None
    data1, data2 = data
    corrs = data1.corrwith(data2, axis=axis, drop=True, method=method)
    # TODO: test for NAs instead of dropping them.
    # it seems like NA can happen if there are all zeros
    # in one of the vectors, so let's dropna for now
    corrs.dropna(inplace=True)

    assert (corrs >= threshold).all(), 'the samples which did not pass the test are:\n{}'.format(corrs[corrs<threshold].sort_values())

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
