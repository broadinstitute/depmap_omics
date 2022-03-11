From r-base:4.1.3

RUN apt-get update && apt-get install -y python3
# install several things that are important for building/installing other packages
RUN apt-get install -y python3-pip python3-dev python3-venv build-essential libxml2-dev zlib1g-dev autoconf libtool bison flex

RUN python3 -m pip install https://github.com/broadinstitute/genepy/archive/refs/heads/master.zip

COPY setup.cfg /install/setup.cfg
COPY pyproject.toml /install/pyproject.toml
COPY depmapomics /install/depmapomics
WORKDIR /install
RUN python3 -m pip install -e .
