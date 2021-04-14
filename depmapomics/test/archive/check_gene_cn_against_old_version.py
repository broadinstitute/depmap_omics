from depmapomics.test.archive.compare_cn import plot_gene_cn_comparison, plot_matrix_comparison
from depmapomics.test.config import VIRTUAL_RELEASE, REFERENCE_RELEASE


CCLE_gene_cn_12 = plot_gene_cn_comparison(release1 = REFERENCE_RELEASE,
                                            release2 = VIRTUAL_RELEASE, savefig=True)

plot_matrix_comparison('CCLE_expression_full')
