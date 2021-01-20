from gsheets import Sheets
from depmapomics.qc import get_cell_line_set as gcl

def get_expected_lines(sheets_url):
    sheets_obj = Sheets.from_files('~/.client_secret.json', '~/.storage.json')
    sheets = sheets_obj.get(sheets_url).sheets
    release = sheets[0].to_frame(header=0, index_col=None)    
    release.columns = release.columns.str.lower()
    return release


taiga_dict = {'20q3': {'public': ['public-20q3-3d35', 33],
                    'dmc': ['dmc-20q3-deprecated-never-released--5f55', 18],
                    'internal': ['internal-20q3-00d0', 9]},
            '20q4': {'public': ['public-20q4-a4b3', 46],
                    'dmc': ['dmc-20q4-fcf4', 38],
                    'internal': ['internal-20q4-2540', 45]},
            '21q1': {'public': ['public-21q1-4b39', 9],
                    'dmc': ['dmc-21q1-0e11', 10],
                    'ibm': ['ibm-21q1-abd9', 12],
                    'internal': ['internal-21q1-4fc4', 13]}
            }

taiga_dict_expanded = gcl.propagate_taiga_dict_with_filenames(taiga_dict)
arxspan_dict = gcl.get_all_arxspans(taiga_dict_expanded, verbose=True)

import pickle
pickle.dump(arxspan_dict, open('arxspan_dict.pkl', 'bw'))

arxspan_dict['20q3']['ibm'] = arxspan_dict['20q3']['dmc']
arxspan_dict['20q4']['ibm'] = arxspan_dict['20q4']['dmc']

arxspan_dict = gcl.rename_dict(arxspan_dict)

lines_to_release_20q3_20q4_sheet = 'https://docs.google.com/spreadsheets/d/1YuKEgZ1pFKRYzydvncQt9Y_BKToPlHP-oDB-0CAv3gE/edit#gid=0'
lines_to_release = get_expected_lines(lines_to_release_20q3_20q4_sheet)
text = gcl.pretty_print_diff(arxspan_dict, lines_to_release, quarters = ['20q3', '21q1'], savefile=True)

lines_to_release_20q4_sheet = 'https://docs.google.com/spreadsheets/d/14WRGIsdNwJydWHeaofkeH9XvC6_0N5v7zd75-gV10n4/edit#gid=0'
lines_to_release = get_expected_lines(lines_to_release_20q4_sheet)
text = gcl.pretty_print_diff(arxspan_dict, lines_to_release, quarters = ['20q4', '21q1'], savefile=True)