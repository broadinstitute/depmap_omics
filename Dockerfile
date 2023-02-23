FROM python:3.9-buster

#RUN pip install pipx && \
#  pipx install poetry==1.2.0

RUN pip install poetry==1.2.2

RUN mkdir /install
COPY pyproject.toml poetry.lock README.md /install
COPY depmapomics /install/depmapomics
RUN cd /install && poetry config virtualenvs.create false && poetry install
