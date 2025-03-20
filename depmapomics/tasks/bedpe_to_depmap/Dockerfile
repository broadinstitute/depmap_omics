FROM python:3.8.19
RUN pip install pandas==2.0.3 numpy==1.24.4
COPY bedpe_to_depmap/bedpe_to_depmap.py /home
COPY oncokb_dm/data/oncogene_oncokb_20250205.txt /home/oncogene_oncokb.txt
COPY oncokb_dm/data/tumor_suppressor_oncokb_20250205.txt /home/tumor_suppressor_oncokb.txt
WORKDIR /home
