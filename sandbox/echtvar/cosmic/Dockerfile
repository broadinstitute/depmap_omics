FROM ubuntu


RUN apt-get update
RUN apt-get install -y wget autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev

RUN apt-get -y install git
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git && git clone https://github.com/samtools/bcftools.git && cd bcftools && make install
 # The following is optional:
 #   autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters
ENV export BCFTOOLS_PLUGINS=~/bcftools/plugins

RUN mkdir -p ~/bin/ && \
    wget -O ~/bin/echtvar https://github.com/brentp/echtvar/releases/download/v0.1.7/echtvar \
    && chmod +x ~/bin/echtvar \
    && ~/bin/echtvar -h 