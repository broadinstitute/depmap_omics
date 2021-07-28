import re
from gsheets.api import Sheets
import numpy as np
import pandas as pd
import pytest
from depmapomics.qc.test_compare_to_ref_release import get_both_release_lists_from_taiga
from depmapomics.qc.config import (LEGACY_PATCH_FLAGS, LINES_TO_RELEASE,
                                   REFERENCE_RELEASE, VIRTUAL_RELEASE)
from taigapy import TaigaClient
tc = TaigaClient()

####### FIXTURES ####
def tsv2csv(df):
    df.to_csv('/tmp/data.tsv', index=False)
    df = pd.read_csv('/tmp/data.tsv', sep='\t')
    return df

# def get_both_releases_from_taiga(file):
#     # data1 = tc.get(name='cn-achilles-version-06ca', version=61, file='all_21Q2_segment')
#     # data2 = tc.get(name='cn-achilles-version-06ca', version=70, file='achilles_segment')
#     data1 = tc.get(name=REFERENCE_RELEASE['name'], file=file, version=REFERENCE_RELEASE['version'])
#     data2 = tc.get(name=VIRTUAL_RELEASE['name'], file=file, version=VIRTUAL_RELEASE['version'])
#     if LEGACY_PATCH_FLAGS['tsv2csv']:
#         # some older taiga data formats (probably 20q1 and older) are tsv and deprecated
#         data1 = tsv2csv(data1)

#     if LEGACY_PATCH_FLAGS['rename_column']:
#         # in 21q1 CCLE_mutations file this column was renamed
#         data1.rename(columns={'Tumor_Allele': 'Tumor_Seq_Allele1'}, inplace=True)
#     return data1, data2


@pytest.fixture(scope='module')
def arxspans(request):
    return get_both_release_lists_from_taiga(request.param)



if __name__ == '__main__':
    # CN + mutations
    file = 'CCLE_mutations'
    lines_to_drop = {'ACH-001078', 'ACH-002184', 'ACH-001146',
                     'ACH-002022', 'ACH-001173', 'ACH-001790',
                     'ACH-002260', 'ACH-001741', 'ACH-000010', 'ACH-002475'}

    # expression
    # file = 'CCLE_expression'
    # lines_to_drop = {'ACH-001078', 'ACH-002184', 'ACH-001146',
    #                 'ACH-002022', 'ACH-001173', 'ACH-001790',
    #                 'ACH-002260', 'ACH-001741', 'ACH-000010',
    #                 'ACH-001316', 'ACH-001212'}

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


    def test_unexpected_dropped_lines(unexpected_dropped_lines):
        assert not unexpected_dropped_lines, \
            'Dropped lines: {}'.format(unexpected_dropped_lines)
