import pandas as pd
import numpy as np
import argparse
import os

def get_weighted_avg(dis):
    dis = pd.Series(dis.replace("T:", "").strip().split(" ")).astype(int)
    try:
        return np.average(dis.index+1, weights=dis.values)
    except:
        return np.nan

def get_n_repeats_from_msisensor2_output_dis(path):
    #read the full file
    df = pd.read_csv(path, sep="\t", header=None)
    #split site rows from distribution rows (they alternate)
    df = pd.DataFrame({
        "Site":df.loc[df.index % 2 == 0, 0].values, 
        "Dis":df.loc[df.index % 2 == 1, 0].values
    })
    #split the site into its component info
    site_cols = ["Chromosome", "Location", "LeftFlank", "Repeat", "RightFlank"]
    df = pd.concat([
        pd.DataFrame(df["Site"].str.split(" ").tolist(), columns=site_cols),
        df["Dis"]
    ], axis=1).set_index(site_cols)
    #find weighted average number of repeats
    return df.apply(lambda x: get_weighted_avg(x["Dis"]), axis=1)


if __name__=='__main__':
    #parse args
    parser = argparse.ArgumentParser(description='Aggreate MSISensor Repeats')
    parser.add_argument('msisensor2_output_list', help='File listing MSISensor2 output dis files, ${sample_id}.msisensor2.output_dis')
    parser.add_argument('output_prefix', help='Prefix for output file: ${output_prefix}MicrosatelliteRepeats.csv')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    args = parser.parse_args()

    #get file paths
    with open(args.msisensor2_output_list) as f:
        file_list = f.read().strip().split('\n')

    #find n repeats in each file and combine into a dataframe with rows as microsatellites and columns as models
    ms_repeats = {}
    for path in file_list:
        ms_repeats[os.path.split(path)[1].split('.')[0]] = get_n_repeats_from_msisensor2_output_dis(path)
    ms_repeats = pd.concat(ms_repeats, axis=1)
    ms_repeats.to_csv(os.path.join(args.output_dir, f"{args.output_prefix}.microsatellite_repeats.csv"))