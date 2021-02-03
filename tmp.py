#%%

from taigapy import TaigaClient

tc = TaigaClient() # These two steps could be merged in one with `from taigapy import default_tc as tc`

# import default_tc as tc

metadata = tc.get_dataset_metadata("cclf-cn-f083")

# %%
int(max([x['name'] for x in metadata['versions']]))

# %%
tc.get_canonical_id("achilles-v2-4-6.4")
# %%
