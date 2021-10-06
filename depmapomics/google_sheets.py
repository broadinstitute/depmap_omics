import pandas as pd
import numpy as np
import gspread
from google.oauth2 import service_account
import re
import io
import csv
import urllib.parse

SCOPES = [
    "https://www.googleapis.com/auth/spreadsheets",
    "https://www.googleapis.com/auth/drive.metadata.readonly",
]


class DataframeSheets:
    "A service class which is responsible for translating between Google sheets <-> Pandas DataFrames"

    def __init__(self, credentials_file: str):
        credentials = service_account.Credentials.from_service_account_file(
            credentials_file, scopes=SCOPES
        )
        self.client = gspread.authorize(credentials)

    def read(self, document_id: str) -> pd.DataFrame:
        "Read the first worksheet in the given document into a Pandas DataFrame"
        document_id = _url_to_document_id(document_id)
        spreadsheet = self.client.open_by_key(document_id)

        rows = spreadsheet.get_worksheet(0).get_all_values()

        buf = io.StringIO()
        w = csv.writer(buf)
        for row in rows:
            w.writerow(row)

        buf.seek(0)
        return pd.read_csv(buf, index_col=0)

    def replace(self, document_id: str, df: pd.DataFrame):
        "Replace the content of the first worksheet on the specified document with the contents from the dataframe"
        document_id = _url_to_document_id(document_id)
        spreadsheet = self.client.open_by_key(document_id)
        worksheet = spreadsheet.get_worksheet(0)
        worksheet.clear()
        df = df.replace({np.nan: None})
        column_header = [""] + df.columns.values.tolist()
        rows = [[i] + row for i, row in zip(df.index, df.values.tolist())]
        worksheet.update([column_header] + rows)


def _url_to_document_id(url: str) -> str:
    if url.startswith("http"):
        p = urllib.parse.urlparse(url)
        assert p.netloc == "docs.google.com"
        m = re.match("/spreadsheets/d/([^/]+)/edit", p.path)
        assert m, f"path={p.path}"
        return m.group(1)
    return url
