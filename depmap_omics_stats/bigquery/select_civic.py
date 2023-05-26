from google.cloud import bigquery
from google.oauth2 import service_account

credentials = service_account.Credentials.from_service_account_file('/Users/qqin/.sparkles-cache/service-keys/broad-qqin.json')
project_id = 'depmap-omics'
client = bigquery.Client(credentials= credentials,project=project_id)

def job_query(sql):
    query_job = client.query(sql)
    df = query_job.to_dataframe()
    return df


def main():
    df = job_query("""SELECT COUNTIF(civic_score > 0) AS CIVIC_thresh0, COUNTIF(civic_score >8) AS CIVIC_thresh8, COUNTIF(civic_score >20) as CIVIC_thresh20 FROM `depmap-omics.maf_staging.merged_maf`""")
    print(df)


if __name__ == "__main__":
    main()
    
