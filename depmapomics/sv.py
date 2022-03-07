import pandas as pd 



def load_and_aggregate(names):
    """
    """
    agg = []
    for file in files:
        name = file.split('/')[-1].split('.')[0]
        val = pd.read_csv(file, sep='\t')
        val['sample_id'] = name
        agg.append(val)
    agg = pd.concat(agg)
    return agg
        