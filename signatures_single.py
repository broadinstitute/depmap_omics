# export GOOGLE_APPLICATION_CREDENTIALS=/localstuff/depmap-omics-9764fbdbe040.json
# fissfc entity_tsv -w DepMap_WGS_CN  -p broad-firecloud-ccle -t sample | sed "s/\[\"//g" | sed 's/\"\]//g' > /localstuff/terra_data_table_wgs_depmapwgs_cn.tsv
import pandas as pd
import numpy as np
import signatureanalyzer as sa
import scipy as sp
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--input_maf", type=str, required=True, help='hgvs_maf_masked column')
parser.add_argument("--hg382bit", type=str, required=True, help='hg38 2bit reference')
parser.add_argument("--sample_id", type=str, required=True, help='hg38 2bit reference')

args = parser.parse_args()
HG_PATH = args.hg382bit
maf_input_file = args.input_maf
sample_id = args.sample_id


sig_list = list()

REF = "pcawg_COMPOSITE"

MAX_ITER = 30000


maf = pd.read_csv(maf_input_file)
ref_df, ref_idx = sa.utils.load_reference_signatures(REF, verbose=False)
spectra_df = sa.spectra.get_spectra_from_maf(maf, hgfile=HG_PATH, reference=REF, real_snps = True)[1]
Wref_df = ref_df.set_index('Somatic Mutation Type').iloc[:,:-2]
# Run supervised NMF
res_supervised = sa.supervised_bnmf.supervised_ardnmf(
	spectra_df,
	Wref_df,
	objective='poisson',
	verbose=True,
	max_iter=MAX_ITER)
sig_list.append(res_supervised['H'].iloc[0])


all_signatures = pd.DataFrame(sig_list)
all_signatures.to_csv(sample_id +".model_id_sig_matrix.csv")
select_cols = ['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']
maf[select_cols].to_csv(sample_id + ".somatic_novel.maf_lite.csv")
