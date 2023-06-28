#!/usr/bin/env python
import seaborn as sns

import numpy as np
from itertools import combinations
from upsetplot import plot
from upsetplot import UpSet

from google.cloud import bigquery
from google.oauth2 import service_account
import pandas as pd
from matplotlib import pyplot as plt
#from sklearn import datasets
#from sklearn.tree import DecisionTreeClassifier
#from sklearn import tree

credentials = service_account.Credentials.from_service_account_file('broad-qqin.json')
project_id = 'depmap-omics'
client = bigquery.Client(credentials= credentials,project=project_id)


def job_query_to_dataframe(sql):
    query_job = client.query(sql)
    df = query_job.to_dataframe()
    return df


def query_all_variants():
    df_all = job_query_to_dataframe("""SELECT table_id, row_count, size_bytes
    FROM `depmap-omics.maf_staging.__TABLES__`
    WHERE STARTS_WITH(table_id, "merged")
    ORDER BY table_id DESC""")
    return df_all


def filter_fields_query(fields=['PASS', 'multiallelic', 'af', 'dp', 'is_coding', 'germline', 'pon', 'popaf'], 
                        cutoffs=['Y', 'Y', 0.15, 2, 'Y', 'Y', 'Y', 5],
                        groupby=None):
    if groupby is not None:
        group_sql = " GROUP BY CDS_ID "
    else:
        group_sql = ""

    dfs = []
    for field, cutoff in zip(fields, cutoffs):
        print(field,cutoff)
        if isinstance(cutoff, str):
            df = job_query_to_dataframe(f"""SELECT COUNTIF({field}='{cutoff}') AS {field} FROM `depmap-omics.maf_staging.merged_maf` {group_sql}""")
        elif field in ['germline', 'pon']:
            df = job_query_to_dataframe(f"""SELECT COUNTIF({field}!='{cutoff}') AS {field} FROM `depmap-omics.maf_staging.merged_maf` {group_sql}""")
        else:
            df = job_query_to_dataframe(f"""SELECT COUNTIF(SAFE_CAST(SPLIT({field}, ',')[ORDINAL(1)] AS NUMERIC) >= {cutoff}) AS {field} FROM `depmap-omics.maf_staging.merged_maf` {group_sql}""")
        dfs.append([f"{field}:{cutoff}", df.values[0][0]])
    return dfs



def combo_evaluation(fields=['civic_score>=1', 'civic_score>=8', 
                             'SAFE_CAST(SPLIT(brca1_func_score, ",")[ORDINAL(1)] AS NUMERIC)>=-1.328', 
                             "hess_driver='Y'", "likely_driver='Y'", "likely_lof='Y'"], 
                    rename_fields = ['civic>=1', 'civic>=8', 'brca1_func_score>=-1.328', 'hess_driver=Y', 'likely_driver=Y', 'likely_lof=Y'], 
                    figname='test.pdf'):
    n = 0
    references = fields
    ref_bools = []
    metrics = []
    for index,comb in enumerate([combinations(references, ii) for ii in range(1, 6)]): #, combinations(fields, 3)]):
        for co in comb:
            ref_bool = np.zeros(len(references))
            ref_bool[np.isin(np.array(fields), np.array(co))] = 1
            ref_bools.append(tuple(ref_bool.astype(bool)))
            sql_query = ' AND '.join(list(co))
            n += 1
            df = job_query_to_dataframe(f"""SELECT COUNTIF({sql_query}) AS metric{index} FROM `depmap-omics.maf_staging.merged_maf`""")
            metrics.append(df.values.flatten()[0])

    print(np.array(ref_bools))
    print(np.array(ref_bools).shape)
    print(metrics)
    print(len(metrics))
    all_var = query_all_variants()
    fig = plt.figure()
    test = pd.MultiIndex.from_tuples(ref_bools, names=rename_fields)
    res = pd.DataFrame({'count': metrics}, index=test)
    res = res.assign(ratio=metrics / all_var.row_count.values[0])
    print(res)

    upset = UpSet(res, subset_size='sum', sum_over='count', show_counts=True, min_subset_size=15, intersection_plot_elements=3)
    upset.add_catplot(value='ratio', kind='strip', color='blue')
    upset.plot(fig=fig)

    #plot(res, show_counts=True, min_subset_size=15, fig=fig)
    fig.set_size_inches(28, 5)
    fig.savefig(figname)
    print(res)
    print(n)
    return


def whitelisting_fields_query(fields=['civic_score', 'civic_score', 'brca1_func_score', 'hess_driver', 'likely_driver', 'likely_lof'], 
                              cutoffs=[1, 8, -1.328, 'Y', 'Y', 'Y'],
                              groupby=None):
    if groupby is not None:
        group_sql = " GROUP BY CDS_ID "
    else:
        group_sql = ""

    dfs = []
    for field, cutoff in zip(fields, cutoffs):
        print(field,cutoff)
        if isinstance(cutoff, str):
            df = job_query_to_dataframe(f"""SELECT COUNTIF({field}='{cutoff}') AS {field} FROM `depmap-omics.maf_staging.merged_maf` {group_sql}""")
        elif field == 'civic_score':
            df = job_query_to_dataframe(f"""SELECT COUNTIF({field} >= {cutoff}) AS {field} FROM `depmap-omics.maf_staging.merged_maf` {group_sql}""")
        elif field=='brca1_func_score':
            df = job_query_to_dataframe(f"""SELECT COUNTIF(SAFE_CAST(SPLIT({field}, ',')[ORDINAL(1)] AS NUMERIC) <= {cutoff}) AS {field} FROM `depmap-omics.maf_staging.merged_maf` {group_sql}""")
        else:
            df = job_query_to_dataframe(f"""SELECT COUNTIF(SAFE_CAST(SPLIT({field}, ',')[ORDINAL(1)] AS NUMERIC) >= {cutoff}) AS {field} FROM `depmap-omics.maf_staging.merged_maf` {group_sql}""")
        dfs.append([f"{field}:{cutoff}", df.values[0][0]])
    return dfs


def main():
    #combo_evaluation(figname='whitelisting_combo.pdf')
    #combo_evaluation(fields=["SAFE_CAST(SPLIT(af, ',')[ORDINAL(1)] AS NUMERIC)>=0.15", "SAFE_CAST(SPLIT(dp, ',')[ORDINAL(1)] AS NUMERIC)>=2", "is_coding='Y'", "germline!='Y'", "pon!='Y'", "SAFE_CAST(SPLIT(popaf, ',')[ORDINAL(1)] AS NUMERIC) >= 5"], rename_fields=['af>=0.15', 'dp>=2', 'coding=Y', 'germline!=Y', 'pon!=Y', 'gnomad>=5'], figname='filter_combo.pdf')
    df_filter = pd.DataFrame(filter_fields_query())
    df_whitelist = pd.DataFrame(whitelisting_fields_query())
    df_all = query_all_variants()
    print(df_filter)
    df_filter.columns = ['metric', 'count']
    df_whitelist.columns = ['metric', 'count']
    df_filter.loc[:, 'ratio'] = df_filter.iloc[:, 1] / df_all.row_count.values[0]
    df_whitelist.loc[:, 'ratio'] = df_whitelist.iloc[:, 1] / df_all.row_count.values[0]

    fig, ax = plt.subplots(2, 1, figsize=(8, 6))
    for index, data in zip(range(2), [df_filter, df_whitelist]):
        sns.barplot(x='metric', y='count', data=data, capsize=0.2, ax=ax[index])
        
        # show the mean
        for p,t in zip(ax[index].patches, data.ratio.values):
            h, w, x = p.get_height(), p.get_width(), p.get_x()
            xy = (x + w / 2, h / 1.3)
            text = f'ratio:\n{t:0.2E}'
            ax[index].annotate(text=text, xy=xy, ha='center', va='center')
        ax[index].tick_params(axis='x', labelrotation=45)
    fig.subplots_adjust(hspace=0.7)
    fig.savefig("individual_metrics.pdf")


if __name__ == "__main__":
    main()
