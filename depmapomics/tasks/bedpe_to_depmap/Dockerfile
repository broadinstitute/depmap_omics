FROM python:3.8
RUN pip install pandas
RUN pip install fsspec gcsfs
COPY bedpe_to_depmap/bedpe_to_depmap.py /home
COPY oncokb_dm/data/oncogene_oncokb.txt /home
COPY oncokb_dm/data/tumor_suppressor_oncokb.txt /home
WORKDIR /home