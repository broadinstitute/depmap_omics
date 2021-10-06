from typing import List, Optional, Tuple, Dict
import pandas as pd
from .google_sheets import DataframeSheets
import re


def _validate_tracker_df(df: pd.DataFrame):
    required_fields = [
        ("arxspan_id", str),
        ("version", str),
        ("sm_id", str),
        ("PDO", str),
        ("datatype", str),
        ("size", str),
        ("ccle_name", str),
        ("stripped_cell_line_name", str),
        ("participant_id", str),
        ("cellosaurus_id", str),
        ("bam_public_sra_path", str),
        ("internal_bam_filepath", str),
        ("internal_bai_filepath", str),
        ("legacy_bam_filepath", str),
        ("legacy_bai_filepath", str),
        ("parent_cell_line", str),
        ("sex", str),
        ("matched_normal", str),
        ("age", str),
        ("primary_site", str),
        ("primary_disease", str),
        ("subtype", str),
        ("subsubtype", str),
        ("origin", str),
        ("mediatype", str),
        ("condition", str),
        ("sequencing_type", str),
        ("baits", str),
        ("comments", str),
        ("source", str),
        ("sequencing_date", str),
        ("release_date", str),
        ("bam_qc", str),
        ("processing_qc", str),
        ("crc32c_hash", str),
        ("md5_hash", str),
        ("update_time", str),
        ("legacy_size", str),
        ("legacy_crc32c_hash", str),
        ("low_quality", str),
        ("blacklist", str),
        ("prioritized", str),
    ]
    column_names = set(df.columns)
    for field, type in required_fields:
        assert field in column_names, f"{field} missing from columns: {column_names}"
        # Should we have a check on type? For the moment, I've set them all to str but I bet
        # the -> DataFrame process converts some to floats and I don't know what types the code is assuming
    extra_columns = column_names.difference([field for field, _ in required_fields])
    for column in extra_columns:
        # is it a quarter (ie: 20Q3) column?
        m = re.match("\\d\\d[qQ]\\d", column)
        assert m is not None, f"Did not recognize this column: {column}"

        # buf = io.StringIO()
        # df.to_csv(buf)
        # self.client.import_csv(spreadsheet.id, data=buf.getvalue())


class SampleTracker:
    "A service class can be used to access the sample tracker"

    def __init__(self, sheets_service: DataframeSheets, refsheet_id: str):
        self.sheets = sheets_service
        self.refsheet_id = refsheet_id

    def get_tracker(self):
        df = self.sheets.read(self.refsheet_id)
        _validate_tracker_df(df)
        return df

    def update_tracker(self, df):
        _validate_tracker_df(df)
        assert self.sheets.replace(self.refsheet_id, df)
