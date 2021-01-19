from depmapomics.qc.compare_cn import plot_gene_cn_comparison

gene_cn_20q3_21q1 = plot_gene_cn_comparison(release1 = {'name': 'internal-20q3-00d0', 'version': 9},
                                            release2 = {'name': 'internal-21q1-4fc4', 'version': 7}, savefig=True)

# %%time
# gene_cn_20q3_21q1 = plot_gene_cn_comparison(release1 = {'name': 'internal-20q3-00d0', 'version': 9},
#                                             release2 = {'name': 'cn-latest-d8d4', 'version': 5},
#                                             filenames2 = {'gene_cn': 'gene_cn_all_merged_21Q1', 
#                                                           'segment_cn': 'segments_all_merged_21Q1'})
