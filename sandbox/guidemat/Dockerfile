FROM python:3.9-buster

#RUN pip install pipx && \
#  pipx install poetry==1.2.0

RUN pip install poetry==1.8.3

RUN mkdir /install && mkdir /install/depmapomics && mkdir /install/mgenepy
# first create these as if they are empty modules, 
# just a list of dependencies, so that we can run 
# 'poetry install' and get a layer which has all our
# dependencies. The actual code for these modules will be added
# later, but do it as a different layer, so that the `poetry install`
# step can be avoided if only a py file changed.
COPY pyproject.toml poetry.lock README.md /install/
COPY depmapomics/__init__.py /install/depmapomics/__init__.py
COPY mgenepy/__init__.py /install/mgenepy/__init__.py
RUN cd /install && poetry config virtualenvs.create false
RUN cd /install && poetry install
RUN pip install biomart

# now add on the py code
COPY depmapomics /install/depmapomics
COPY mgenepy /install/mgenepy
RUN cd /install && poetry install


RUN apt-get update && apt-get install -y \
  bzip2 \
  g++ \
  libbz2-dev \
  libcurl4-openssl-dev \
  liblzma-dev \
  make \
  ncurses-dev \
  wget \
  zlib1g-dev

ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
ENV BCFTOOLS_VERSION=1.12

WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
  tar --bzip2 -xf bcftools-$BCFTOOLS_VERSION.tar.bz2

WORKDIR /tmp/bcftools-$BCFTOOLS_VERSION
RUN make prefix=$BCFTOOLS_INSTALL_DIR && \
  make prefix=$BCFTOOLS_INSTALL_DIR install

WORKDIR /
RUN ln -s $BCFTOOLS_INSTALL_DIR/bin/bcftools /usr/bin/bcftools && \
  rm -rf /tmp/bcftools-$BCFTOOLS_VERSION

ENTRYPOINT ["/usr/bin/bcftools"]