.PHONY: qc fullqc

qc:
#	pytest -v --pdb -x -s depmapomics/qc/test_plot_comparisons.py::test_plot_per_gene_means
	pytest -v --pdb -x -s depmapomics/qc/test_plot_comparisons.py::test_plot_matrix_comparison
	

fullqc:
	pytest -v --pdb -x -s
