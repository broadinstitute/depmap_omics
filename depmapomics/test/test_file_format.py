import re
from depmapomics.test.config import FILE_ATTRIBUTES, TEMP_VIRTUAL_TAIGA_ID
from taigapy import TaigaClient
import numpy as np
tc = TaigaClient()
import pytest

IS_MATRIX = [x for x in FILE_ATTRIBUTES if x['ismatrix']]
IS_TABLE = [x for x in FILE_ATTRIBUTES if not x['ismatrix']]

def check_symbol_and_entrezid_in_column(matrix):
    matches = matrix.columns.map(lambda x: re.match(r'([a-zA-Z0-9-_/.]+)\s\((\d+)\)$', x))
    assert  matches.notnull().all(), \
        'some columns do not follow the symbol (entrez id) format. The first few are: \n{}'.format(matrix.columns[matches.isnull()][:20])

def check_symbol_and_enst_in_column(matrix):
    p1 = r'(?:[a-zA-Z0-9-_/.]+)' # only gene symbol
    p2 = r'ENST(?:\d{11})' # only ensembl id
    p3 = p1 + r'\s\(' + p2 + r'\)' # gene symbol (ensembl id)
    p4 = r'ERCC-(?:\d{5})' # ERCC
    p  = '|'.join([p2, p3, p4])
    matches = matrix.columns.map(lambda x: re.fullmatch(p, x))
    assert  matches.notnull().all(), \
        'some columns do not follow the symbol (ensembl transcript id) format. The first few are: \n{}'.format(matrix.columns[matches.isnull()][:20])

def check_symbol_and_ensg_in_column(matrix):
    p1 = r'(?:[a-zA-Z0-9-_/.]+)' # only gene symbol
    p2 = r'ENSG(?:\d{11})' # only ensembl id
    p3 = p1 + r'\s\(' + p2 + r'\)' # gene symbol (ensembl id)
    p4 = r'ERCC-(?:\d{5})' # ERCC
    p  = '|'.join([p2, p3, p4])
    matches = matrix.columns.map(lambda x: re.fullmatch(p, x))
    assert  matches.notnull().all(), \
        'some columns do not follow the symbol (ensembl gene id) format. The first few are: \n{}'.format(matrix.columns[matches.isnull()][:20])

def check_symbol_and_id(matrix, gene_id):
    if gene_id == 'ensg':
        check_symbol_and_ensg_in_column(matrix)
    elif gene_id == 'enst':
        check_symbol_and_enst_in_column(matrix)
    elif gene_id == 'entrez':
        check_symbol_and_entrezid_in_column(matrix)
    else:
        raise('unknown id provided')

def check_arxspan_ids(column):
    matches = column.map(lambda x: re.match(r'ACH-[\d]{6}$', x))
    assert  matches.notnull().all(), \
        'some rows do not follow the ACH-xxxxxx format. The first few are: \n{}'.format(column.index[matches.isnull()][:20])

def check_table_format(data):
    assert 'DepMap_ID' in data.columns, 'no DepMap_ID column found'
    check_arxspan_ids(data['DepMap_ID'])

def check_matrix_format(matrix, gene_id=None, hasNA=False):
    # check if column format is symbol (entrez id)
    check_symbol_and_id(matrix, gene_id)

    # check if rows start with arxspan IDs
    check_arxspan_ids(matrix.index)

    # check if there are any NAs in the data
    assert hasNA | matrix.notnull().all().all(), 'null values identified in the matrix'

    # check that data is nothing but float or int
    datatypes = set(matrix.stack().map(type))
    assert datatypes - {float, int} == set(), 'matrix contains data other than int and float. Datatypes are: {}'.format(datatypes)
    assert len(datatypes) == 1, 'there is a mixture of datatypes in the matrix'

@pytest.mark.parametrize('files_attribute', IS_MATRIX)
def test_matrix_format(files_attribute, temp_virtual_taiga_id=TEMP_VIRTUAL_TAIGA_ID, version=None):
    data = tc.get(name=temp_virtual_taiga_id, file=files_attribute['file'], version=version)
    check_matrix_format(data, gene_id=files_attribute['gene_id'], hasNA=files_attribute['hasNA'])

@pytest.mark.parametrize('files_attribute', IS_TABLE)
def test_table_format(files_attribute, temp_virtual_taiga_id=TEMP_VIRTUAL_TAIGA_ID, version=None):
    data = tc.get(name=temp_virtual_taiga_id, file=files_attribute['file'], version=version)
    check_table_format(data)

@pytest.mark.xfail(strict=True)
def test_logtransform(file='CCLE_expression_full', temp_virtual_taiga_id=TEMP_VIRTUAL_TAIGA_ID, version=None):
    CCLE_expression_full = tc.get(name=temp_virtual_taiga_id, version=version, file=file)
    assert CCLE_expression_full.min().min() == 0
    assert not CCLE_expression_full.sum(axis=1).map(lambda x: np.isclose(x, 1e6, rtol=1e-4)).all(), 'expression data is not log-transformed'
