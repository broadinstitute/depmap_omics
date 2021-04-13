import re
from depmapomics.qc.config import FILE_ATTRIBUTES, TEMP_VIRTUAL_TAIGA_ID
from taigapy import TaigaClient

def check_symbol_and_entrezid_in_column(matrix):
    matches = matrix.columns.map(lambda x: re.match(r'([a-zA-Z0-9-_/.]+)\s\((\d+)\)$', x))
    assert  matches.notnull().all(), \
        'some columns do not follow the symbol (entrez id) format. The first few are: \n{}'.format(matrix.columns[matches.isnull()][:20])

def check_symbol_and_enst_in_column(matrix):
    p1 = r'(?:[a-zA-Z0-9-_/.]+)' # only gene symbol
    p2 = r'ENST(?:\d{11})' # only ensembl id
    p3 = p1 + r'\s\(' + p2 + '\)' # gene symbol (ensembl id)
    p4 = r'ERCC-(?:\d{5})' # ERCC
    p  = '|'.join([p2, p3, p4])
    matches = matrix.columns.map(lambda x: re.fullmatch(p, x))
    assert  matches.notnull().all(), \
        'some columns do not follow the symbol (ensembl transcript id) format. The first few are: \n{}'.format(matrix.columns[matches.isnull()][:20])

def check_symbol_and_ensg_in_column(matrix):
    p1 = r'(?:[a-zA-Z0-9-_/.]+)' # only gene symbol
    p2 = r'ENSG(?:\d{11})' # only ensembl id
    p3 = p1 + r'\s\(' + p2 + '\)' # gene symbol (ensembl id)
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

def check_data_format(files_attributes=FILE_ATTRIBUTES,
        temp_virtual_taiga_id=TEMP_VIRTUAL_TAIGA_ID, version=None):
    tc = TaigaClient()

    for files_attribute in files_attributes:
        data = tc.get(name=temp_virtual_taiga_id, file=files_attribute['file'], version=version)
        print('checking ', files_attribute['file'], end=': ')
        if files_attribute['ismatrix']:
            check_matrix_format(data, gene_id=files_attribute['gene_id'], hasNA=files_attribute['hasNA'])
            print('PASS')
        else:
            check_table_format(data)
            print('PASS')
        del data
