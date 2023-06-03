FROM mambaorg/micromamba:1.4.3

RUN /bin/bash -c "micromamba create -qy -c conda-forge -c bioconda -c defaults -n vep ensembl-vep git gh htslib samtools bcftools ucsc-liftover"
RUN echo "source activate vep" > ~/.bashrc
ENV PATH /opt/conda/envs/vep/bin:${PATH}
RUN git clone https://github.com/qinqian/vcf2maf
