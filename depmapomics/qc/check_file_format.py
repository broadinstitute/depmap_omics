import re

def check_matrix_columns_format(matrix):
    # check if column format is symbol (entrez id)
    matches = matrix.columns.map(lambda x: re.match(r'([a-zA-Z0-9-_/.]+)\s\((\d+)\)$', x))
    assert  matches.notnull().all(), \
        'some columns do not follow the symbol (entrez id) format. The first few are: \n{}'.format(matrix.columns[matches.isnull()][:20])

    # check if rows start with arxspan IDs
    matches = matrix.index.map(lambda x: re.match(r'ACH-[\d]{6}$', x))
    assert  matches.notnull().all(), \
        'some rows do not follow the ACH-xxxxxx format. The first few are: \n{}'.format(matrix.index[matches.isnull()][:20])

    assert matrix.notnull().all().all(), 'null values identified in the matrix'

    # check that data is nothing but float or int
    datatypes = set(matrix.stack().map(type))
    assert datatypes - {float, int} == set(), 'matrix contains data other than int and float. Datatypes are: {}'.format(datatypes)
    assert len(datatypes) == 1, 'there is a mixture of datatypes in the matrix'
