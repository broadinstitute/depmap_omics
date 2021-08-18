from gsheets import Sheets
from depmapomics.qc.archive import get_cell_line_set as gcl
from depmapomics.qc.config import LINES_TO_RELEASE_SHEET
import pickle
import pandas as pd

# PREV_VIRTUAL = {'public': 'public-21q2-110d',
#                 'ibm': 'ibm-21q2-9ed1',
#                 'dmc': 'dmc-21q2-27e1',
#                 'internal': 'internal-21q2-9d16'}

# VIRTUAL = {'internal': 'internal-21q3-fe4c',
#            'ibm': 'ibm-21q3-179f',
#            'dmc': 'dmc-21q3-482c',
#            'public': 'public-21q3-bf1e'}

taiga_dict = {'21q2': {'public': ['public-21q2-110d', 13],
                       'ibm': ['ibm-21q2-9ed1', 15],
                       'dmc': ['dmc-21q2-27e1', 14],
                       'internal': ['internal-21q2-9d16', 17]},
              '21q3': {'internal': ['internal-21q3-fe4c', 12],
                       'ibm': ['ibm-21q3-179f', 8],
                       'dmc': ['dmc-21q3-482c', 7],
                       'public': ['public-21q3-bf1e', 7]}
             }



def main(recompute=False, lines_to_remove={'ACH-001189', 'ACH-002303', 'ACH-002315', 'ACH-002341', 'ACH-002359'}):

  if recompute:
      taiga_dict_expanded = gcl.propagate_taiga_dict_with_filenames(taiga_dict)
      arxspan_dict = gcl.get_all_arxspans(taiga_dict_expanded, verbose=True)
      # arxspan_dict['20q3']['ibm'] = arxspan_dict['20q3']['dmc']
      arxspan_dict['20q4']['ibm'] = arxspan_dict['20q4']['dmc']
      pickle.dump(arxspan_dict, open('arxspan_dict.pkl', 'bw'))
  else:
      arxspan_dict = pickle.load(open('arxspan_dict.pkl', 'br'))


  # lines_to_remove = set()
  lines_to_release_21q1 = gcl.get_expected_lines(lines_to_release_sheet)
  text = gcl.pretty_print_diff(arxspan_dict, LINES_TO_RELEASE_SHEET, lines_to_remove=lines_to_remove, quarters = ['20q4', '21q1'], savefile=True)

  gcl.check_if_fusion_has_expression_released(arxspan_dict, quarter='21q1')
  gcl.check_acciddental_release(arxspan_dict, quarter = '21q1')
#   gcl.plot_diff_heatmap(arxspan_dict, lines_to_release_21q1, lines_to_remove=lines_to_remove, quarters = ['20q4', '21q1'])


  # lines_to_drop_dict = {
  #     'ACH-002217': "no bam found: Sanger",
  #     'ACH-002335': "no bam found: Chordoma",
  #     'ACH-002378': "no bam found: Sanger",
  #     'ACH-001956': "failed qc wes (chordoma)",
  #     'ACH-001955': "failed qc wes (chordoma)",
  #     'ACH-001957': "failed qc wes (chordoma)",
  #     'ACH-001011': 'wrong new (should only have HC)',
  #     'ACH-001108': 'wrong new (should only have HC)',
  #     'ACH-001187': 'wrong new (should only have HC)',
  #     'ACH-001189': 'wrong (from 20Q4)',
  #     'ACH-002303': 'wrong (from 20Q4)',
  #     'ACH-002315': 'wrong (from 20Q4)',
  #     'ACH-001675': 'wrong (from 20Q4)',
  #     'ACH-003000': 'wrong engineered',
  #     'ACH-002875': 'wrong engineered',
  #     'ACH-002874': 'wrong engineered'
  # }
