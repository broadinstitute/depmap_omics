#!/usr/bin/env python

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
        else:
            df = job_query_to_dataframe(f"""SELECT COUNTIF(SAFE_CAST(SPLIT({field}, ',')[ORDINAL(1)] AS NUMERIC) >= {cutoff}) AS {field} FROM `depmap-omics.maf_staging.merged_maf` {group_sql}""")
        dfs.append([f"{field}:{cutoff}", df.values[0][0]])
    return dfs



def combo_whitelisting():
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
    df_filter = filter_fields_query()
    df_whitelist = whitelisting_fields_query()
    print(df_filter)
    print(df_whitelist)
    df_all = query_all_variants()
    print(df_all)
    df_filter = pd.DataFrame(df_filter)
    df_whitelist = pd.DataFrame(df_whitelist)
    df_filter.iloc[:, 1] /= df_all.row_count.values[0]
    df_whitelist.iloc[:, 1] /= df_all.row_count.values[0]
    print(df_filter)
    print(df_whitelist)
    df_all.to_csv("bigquery_civic.csv")


if __name__ == "__main__":
    main()
