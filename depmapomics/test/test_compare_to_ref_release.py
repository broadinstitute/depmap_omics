from depmapomics.test.config import (FILE_ATTRIBUTES, VIRTUAL_RELEASE, REFERENCE_RELEASE)
import pytest
from taigapy import TaigaClient

tc = TaigaClient()

@pytest.fixture(scope='module')
def data(request):
    data1 = tc.get(name=REFERENCE_RELEASE['name'], file=request.param, version=REFERENCE_RELEASE['version'])
    data2 = tc.get(name=VIRTUAL_RELEASE['name'], file=request.param, version=VIRTUAL_RELEASE['version'])
    return data1, data2

PARAMS_matrix_persample_correlations = [('CCLE_expression_full', 0.95), ('CCLE_gene_cn', 0.95)]
@pytest.mark.parametrize('method', ['spearman', 'pearson'])
@pytest.mark.parametrize('axisname', ['pergene', 'persample'])
@pytest.mark.parametrize('data, threshold', PARAMS_matrix_persample_correlations, indirect=['data'])
def test_matrix_correlations(data, threshold, axisname, method):
    if axisname == 'pergene':
        axis = 0
    elif axisname == 'persample':
        axis = 1
    else:
        raise Exception('unknown input provided')
    data1, data2 = data
    corrs = data1.corrwith(data2, axis=axis, drop=True, method=method)
    assert (corrs >= threshold).all(), 'the samples which did not pass the test are:\n{}'.format(corrs[corrs<threshold].sort_values())
