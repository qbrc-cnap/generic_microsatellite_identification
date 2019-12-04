FROM debian:stretch

RUN apt-get update -y && \
    apt-get install -y \
        build-essential \
        wget \
        zlib1g-dev \
        libcurl4-openssl-dev \
        libncurses-dev \
        libbz2-dev \
        liblzma-dev \
        libxml2-dev \
        zip \
        python3-dev \
        python3-pip \
        openjdk-8-jdk \
        pandoc \
        git \
        r-base

# Install all software under /opt/software:
RUN mkdir -p /opt/software 

# Install some Python3 libraries:
ADD requirements.txt /opt/software/
RUN pip3 install -r /opt/software/requirements.txt

# Install cutadapt
RUN pip3 install --user --upgrade cutadapt

# Install bbmap
RUN cd /opt/software && \
    wget https://sourceforge.net/projects/bbmap/files/BBMap_38.73.tar.gz && \
    tar xzvf BBMap_38.73.tar.gz && \
    rm /opt/software/BBMap_38.73.tar.gz

# Install PHOBOS
ADD phobos-v3.3.12-linux.zip /opt/software
RUN cd /opt/software && \
    unzip phobos-v3.3.12-linux.zip && \
    rm /opt/software/phobos-v3.3.12-linux.zip

# Install bowtie2
RUN cd /opt/software && \
    wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip && \
    unzip bowtie2-2.3.5.1-linux-x86_64.zip && \
    rm /opt/software/bowtie2-2.3.5.1-linux-x86_64.zip
ENV PATH="/opt/software/bowtie2-2.3.5.1-linux-x86_64/:${PATH}"

# Install samtools
RUN cd /opt/software && \
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xjf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure --prefix=/opt/software/samtools && \
    make && \
    make install && \
    rm /opt/software/samtools-1.9.tar.bz2
ENV PATH="/opt/software/samtools/bin:${PATH}"

# Install seqtk
RUN cd /opt/software && \
    wget https://github.com/lh3/seqtk/archive/v1.3.tar.gz && \
    tar xzvf v1.3.tar.gz && \
    cd ./seqtk-1.3 && \
    make && \
    rm /opt/software/v1.3.tar.gz
ENV PATH="/opt/software/seqtk-1.3:${PATH}"

# Add the python script for parsing the PHOBOS output file
ADD parse_phobos.py /opt/software