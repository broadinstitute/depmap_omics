FROM python:3.9.18-bullseye
RUN pip install pandas==1.5.3
COPY bedpe_to_depmap/bedpe_to_depmap.py /home
COPY oncokb_dm/data/oncogene_oncokb.txt /home
COPY oncokb_dm/data/tumor_suppressor_oncokb.txt /home
WORKDIR /home
