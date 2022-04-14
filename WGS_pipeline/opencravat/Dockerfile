FROM python:3.6

ARG PACKAGEDIR=/usr/local/lib/python3.6/site-packages/cravat
ARG BRANCH=master

RUN apt update && apt install -y vim

RUN git clone https://github.com/KarchinLab/open-cravat.git --single-branch --branch $BRANCH && \
    pip install ./open-cravat && \
    pip install gsutil open-cravat-multiuser aiosqlite3 scipy && \


    #    oc version

    RUN mv $PACKAGEDIR/conf /mnt/conf && ln -s /mnt/conf $PACKAGEDIR/conf && \
    mv $PACKAGEDIR/modules /mnt/modules && ln -s /mnt/modules $PACKAGEDIR/modules && \
    mv $PACKAGEDIR/jobs /mnt/jobs && ln -s /mnt/jobs $PACKAGEDIR/jobs && \
    mv $PACKAGEDIR/logs /mnt/logs && ln -s /mnt/logs $PACKAGEDIR/logs
VOLUME /mnt/conf
VOLUME /mnt/modules
VOLUME /mnt/jobs
VOLUME /mnt/logs

VOLUME /tmp/job
WORKDIR /tmp/job