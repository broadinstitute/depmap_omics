FROM rocker/rstudio:latest
MAINTAINER Jeremie Kalfon

RUN apt-get update -y && apt-get install -y curl git openssl libxml2-dev libz-dev zlib1g-dev libfontconfig1-dev
RUN R -e "install.packages(c('BiocManager', 'devtools', 'testthat'), dependencies=TRUE, repos='http://cran.rstudio.com/');"
RUN R -e "BiocManager::install(c('data.table', 'GenomicRanges', 'optparse', 'parallel'), ask=FALSE);\
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE);\
  devtools::install_github('mskilab/gUtils');"