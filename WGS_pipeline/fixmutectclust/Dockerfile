FROM rocker/rstudio:latest
MAINTAINER Jeremie Kalfon

RUN apt-get update -y && apt-get install -y curl git openssl libxml2-dev libz-dev zlib1g-dev libfontconfig1-dev
RUN R -e "install.packages(c('BiocManager', 'devtools', 'testthat','vcfR', 'fuzzyjoin', 'tidyverse'), dependencies=TRUE, repos='http://cran.rstudio.com/');" &&\
    R -e "BiocManager::install('IRanges')";