FROM python:3.9-buster

#RUN pip install pipx && \
#  pipx install poetry==1.2.0

RUN pip install poetry==1.8.3

# Install poetry package manager
RUN pip install poetry==1.8.3

# Create directory structure for installation
RUN mkdir /install && mkdir /install/depmapomics && mkdir /install/mgenepy
# first create these as if they are empty modules,
# just a list of dependencies, so that we can run
# 'poetry install' and get a layer which has all our
# dependencies. The actual code for these modules will be added
# later, but do it as a different layer, so that the `poetry install`
# step can be avoided if only a py file changed.

# First install dependencies only to create a cacheable layer
# Copy dependency files but not actual code
COPY pyproject.toml poetry.lock README.md /install/
COPY depmapomics/__init__.py /install/depmapomics/__init__.py
COPY depmapomics/__init__.py /install/depmapomics/__init__.py 
COPY mgenepy/__init__.py /install/mgenepy/__init__.py

# Configure poetry and install dependencies
RUN cd /install && poetry config virtualenvs.create false
RUN pip install setuptools==58.1.0 && cd /install && poetry install
RUN pip install biomart

# now add on the py code
# Copy full source code and reinstall with actual code
COPY depmapomics /install/depmapomics
COPY mgenepy /install/mgenepy
RUN cd /install && poetry install
