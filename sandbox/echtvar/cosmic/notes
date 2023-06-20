FROM ubuntu


RUN mkdir -p ~/bin/
RUN wget -O ~/bin/echtvar https://github.com/brentp/echtvar/releases/download/v0.1.7/echtvar \
    && chmod +x ~/bin/echtvar \
    && ~/bin/echtvar -h 

RUN sudo apt-get install libz-dev
RUN sudo apt-get install liblzma-dev
RUN sudo apt-get install -y libbz2-dev
RUN sudo apt-get install libcurl4-gnutls-dev


RUN git clone --recurse-submodules https://github.com/samtools/htslib.git
RUN git clone https://github.com/samtools/bcftools.git
RUN cd bcftools && sudo make install
 # The following is optional:
 #   autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters
ENV export BCFTOOLS_PLUGINS=~/bcftools/plugins