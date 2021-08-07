import os
from json import dump
from dump_workspace import dumb_wdl_parse,rewrite_with_relative_imports

def test_wdl_parse():
    wdl="""
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/RNA_pipeline/rsem_aggregate_results.wdl" as rsem_aggregate_results
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/RNA_pipeline/Aggregate_Fusion_Calls.wdl" as Aggregate_Fusion_Calls


workflow RNA_aggregate {

  #rsem_aggregate_results

  call rsem_aggregate_results.rsem_aggregate_results

  call Aggregate_Fusion_Calls.aggregate_set_files

}
"""
    p = dumb_wdl_parse(wdl)
    assert len(p.statements) == 5
    p2 = rewrite_with_relative_imports(p, lambda path: os.path.basename(path))
    assert p2.str() == """
import "rsem_aggregate_results.wdl" as rsem_aggregate_results
import "Aggregate_Fusion_Calls.wdl" as Aggregate_Fusion_Calls


workflow RNA_aggregate {

  #rsem_aggregate_results

  call rsem_aggregate_results.rsem_aggregate_results

  call Aggregate_Fusion_Calls.aggregate_set_files

}
"""

# def dl(path):
#     filename = os.path.join(dest_dir, os.path.basename(stmt.path))
#     with open(filename, "wb") as fd:
#         fd.write(resp.content)


# p = dumb_wdl_parse(wdl)
# p2 = rewrite_with_relative_imports(p, "export")
# print(p2.str())
