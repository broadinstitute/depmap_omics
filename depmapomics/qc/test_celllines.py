import re
from gsheets.api import Sheets
import numpy as np
import pandas as pd
import pytest
from depmapomics.qc.config import (LEGACY_PATCH_FLAGS, LINES_TO_RELEASE,
                                   REFERENCE_RELEASE, VIRTUAL_RELEASE)
from taigapy import TaigaClient

tc = TaigaClient()


def get_arxspan_ids(data, isMatrix):
    if isMatrix:
        arxspans = set(data.index)
    else:
        assert 'DepMap_ID' in data.columns
        arxspans = set(data['DepMap_ID'])

    matches = [re.match(r'ACH-[\d]{6}$', x) for x in arxspans]
    assert all([x is not None for x in matches]), \
        "At least some arxspans do not match the ACH-#### format. Here are a few examples:\n {}".format(list(arxspans)[:5])
    return arxspans


####### FIXTURES ####
def tsv2csv(df):
    df.to_csv('/tmp/data.tsv', index=False)
    df = pd.read_csv('/tmp/data.tsv', sep='\t')
    return df

def get_both_release_lists_from_taiga(file):
    # data1 = tc.get(name='cn-achilles-version-06ca', version=61, file='all_21Q2_segment')
    # data2 = tc.get(name='cn-achilles-version-06ca', version=70, file='achilles_segment')
    data1 = tc.get(name=REFERENCE_RELEASE['name'], file=file, version=REFERENCE_RELEASE['version'])
    data2 = tc.get(name=VIRTUAL_RELEASE['name'], file=file, version=VIRTUAL_RELEASE['version'])
    if LEGACY_PATCH_FLAGS['tsv2csv']:
        # some older taiga data formats (probably 20q1 and older) are tsv and deprecated
        data1 = tsv2csv(data1)

    arxspans1 = get_arxspan_ids(data1, 'DepMap_ID' not in data1.columns)
    arxspans2 = get_arxspan_ids(data2, 'DepMap_ID' not in data2.columns)

    return arxspans1, arxspans2

@pytest.fixture(scope='module')
def arxspans(request):
    return get_both_release_lists_from_taiga(request.param)




# CN
# file = 'CCLE_gene_cn'
# lines_to_drop = {'ACH-001078', 'ACH-002184', 'ACH-001146',
#                  'ACH-002022', 'ACH-001173', 'ACH-001790',
#                  'ACH-002260', 'ACH-001741', 'ACH-000010', 'ACH-002475'}

# expression
file = 'CCLE_expression'
lines_to_drop = {'ACH-001078', 'ACH-002184', 'ACH-001146',
                 'ACH-002022', 'ACH-001173', 'ACH-001790',
                 'ACH-002260', 'ACH-001741', 'ACH-000010',
                 'ACH-001316', 'ACH-001212'}

arxspans1, arxspans2 = get_both_release_lists_from_taiga(file)

lines_to_release_all = set(LINES_TO_RELEASE.stack())
lines_to_add = lines_to_release_all - arxspans1

dropped_lines = arxspans1 - arxspans2
added_lines = arxspans2 - arxspans1

unexpected_added_lines = added_lines - lines_to_release_all
unexpected_dropped_lines = dropped_lines - lines_to_drop
failed_to_drop = lines_to_drop & arxspans2
# missing_lines = lines_to_release_all - arxspans2
print('unexpected_dropped_lines: ', unexpected_dropped_lines)
print('failed_to_drop: ', failed_to_drop)
print('unexpected_added_lines: ', unexpected_added_lines)

# arxspan = 'ACH-001981'
# print(arxspan in lines_to_release_all)
# print(arxspan in unexpected_added_lines)
# print(arxspan in arxspans1, arxspan in arxspans2)
# print('ACH-002707' in arxspans2)
