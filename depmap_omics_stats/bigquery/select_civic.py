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


def job_query(sql):
    query_job = client.query(sql)
    df = query_job.to_dataframe()
    return df


def main():
    df_civic = job_query("""SELECT COUNTIF(civic_score >= 8) AS CIVIC_thresh8 FROM `depmap-omics.maf_staging.merged_maf`""")
    #df_civic_samples = job_query("""SELECT COUNTIF(civic_score > 0) AS CIVIC_thresh0, COUNTIF(civic_score > 5) AS CIVIC_thresh5, COUNTIF(civic_score > 8) AS CIVIC_thresh8, COUNTIF(civic_score >20) as CIVIC_thresh20 FROM `depmap-omics.maf_staging.merged_maf` GROUP BY CDS_ID""")

    df_pon = job_query("""SELECT COUNTIF(panel_of_normals!='') AS pon FROM `depmap-omics.maf_staging.merged_maf`""")
    #df_pon_samples = job_query("""SELECT COUNTIF(panel_of_normals='') AS pon FROM `depmap-omics.maf_staging.merged_maf` GROUP BY CDS_ID""")
    print(df_pon)

    df_af = job_query("""SELECT COUNTIF(CAST(SPLIT(af, ',')[ORDINAL(1)] AS NUMERIC) >= 0.15) AS af FROM `depmap-omics.maf_staging.merged_maf`""")

    df_dp = job_query("""SELECT COUNTIF(CAST(SPLIT(dp, ',')[ORDINAL(1)] AS NUMERIC) >= 2) AS dp FROM `depmap-omics.maf_staging.merged_maf`""")

    df_dann = job_query("""SELECT COUNTIF(SAFE_CAST(dann_score AS NUMERIC) >= 0.96) AS dann FROM `depmap-omics.maf_staging.merged_maf`""")

    df_gnomad3 = job_query("""SELECT COUNTIF(SAFE_CAST(popaf AS NUMERIC) >= 3) AS gnomad3 FROM `depmap-omics.maf_staging.merged_maf`""")

    df_gnomad5 = job_query("""SELECT COUNTIF(SAFE_CAST(popaf AS NUMERIC) >= 5) AS gnomad5 FROM `depmap-omics.maf_staging.merged_maf`""")

    df_all = job_query("""SELECT table_id, row_count, size_bytes
    FROM `depmap-omics.maf_staging.__TABLES__`
    WHERE STARTS_WITH(table_id, "merged")
    ORDER BY table_id DESC""")

    print(df_all)
    #df = pd.concat([df_civic, df_civic_samples, af], axis=0)
    df = pd.concat([df_civic, df_pon, df_af, df_dp, df_dann, df_gnomad3, df_gnomad5], axis=0)
    print(df)
    print(df/df_all.values.flatten()[1])
    df.to_csv("bigquery_civic.csv")


if __name__ == "__main__":
    main()

