# syntax=docker/dockerfile:1
FROM ubuntu:latest

# install app dependencies
RUN apt-get update && apt-get install -y python3 python3-pip
RUN apt-get -y install git

# install samtools
RUN apt-get -y install gcc &&  apt-get -y install make && apt-get -y install libbz2-dev && apt-get -y install zlib1g-dev && apt-get -y install libncurses5-dev && apt-get -y install libncursesw5-dev && apt-get -y install liblzma-dev
RUN apt-get -y install wget
RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2
RUN tar -vxjf samtools-1.19.2.tar.bz2
RUN cd samtools-1.19.2 && ./configure && make && make install && cd ..
RUN export PATH=$PATH:~/samtools-1.19.2
RUN cp `pwd`/samtools-1.19.2/samtools /usr/local/bin

# install minimap2
RUN git clone https://github.com/lh3/minimap2
RUN cd minimap2 && make
RUN export PATH=$PATH:~/minimap2
RUN cp `pwd`/minimap2/minimap2 /usr/local/bin

# install app
#COPY longReadMapper.wdl /

# final configuration
#ENV FLASK_APP=hello
#EXPOSE 8000
#CMD ["flask", "run", "--host", "0.0.0.0", "--port", "8000"]
