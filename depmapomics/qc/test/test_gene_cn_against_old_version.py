from depmapomics.qc.compare_cn import plot_gene_cn_comparison

gene_cn_20q3_21q1 = plot_gene_cn_comparison(release1 = {'name': 'internal-20q3-00d0', 'version': 9},
                                            release2 = {'name': 'internal-21q1-4fc4', 'version': 13}, savefig=True)

# gene_cn_20q3_21q1 = plot_gene_cn_comparison(release1 = {'name': 'internal-20q3-00d0', 'version': 9},
#                                             release2 = {'name': 'internal-21q1-4fc4', 'version': 7}, savefig=True)

