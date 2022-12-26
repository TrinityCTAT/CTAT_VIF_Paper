FROM --platform=linux/amd64 ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive

#~~~~~~~~~~~~~~
# Install Preliminaries
#~~~~~~~~~~~~~~
RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jre \
    g++ \
    gcc \
    git \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libdb-dev \
    libjpeg-dev \
    liblzma-dev \
    libpng-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    make \
    pbzip2 \
    pigz \
    python \
    python-pip \
    python-setuptools \
    python3-pip \
    r-base \
    unzip \
    wget \
    x11-common \
    zlib1g \
    zlib1g-dev \
    zlibc



#~~~~~~~~~~~
# Environmental Variables 
#~~~~~~~~~~~
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

WORKDIR $SRC

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda install -y bioconda::blat=35 

RUN conda install -y bioconda::bedtools=2.28.0 
RUN conda install -y bioconda::bowtie2=2.3.4.1

#RUN conda install -y anaconda::openjdk=8.0.152 && \

RUN conda install -y anaconda::openjdk
RUN conda install -y bioconda::fastqc



RUN conda install -y bioconda::pyfaidx=0.5.5.2 && \
    conda install -y anaconda::python=2.7.15 && \
    conda install -y bioconda::pysam=0.15.2 && \
    conda install -y conda-forge::pandas=0.24.2 

RUN conda install -y bioconda::trim-galore=0.6.2 

#RUN conda install -y bioconda::samtools=1.9
#RUN conda install -c bioconda multiqc=1.7

RUN conda install -y bioconda::deeptools=3.3.0 

#~~~~~~~~~~~
## Bedtools
#~~~~~~~~~~~
WORKDIR $SRC
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
   tar xf bedtools-2.28.0.tar.gz && \
   cd bedtools2 && \
   make && \
   cp bin/* $BIN/ && \
   rm $SRC/bedtools-2.28.0.tar.gz && \
   rm -r $SRC/bedtools2


WORKDIR $SRC
ENV SAMTOOLS_VERSION=1.9
RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    cd $SRC && \
    wget $SAMTOOLS_URL && \
    tar xf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install && \
    rm $SRC/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    rm -r $SRC/samtools-${SAMTOOLS_VERSION}



#~~~~~~~~~~~~~~
# BowTie2
#~~~~~~~~~~~~~~
WORKDIR $SRC
RUN wget --no-check-certificate https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-source.zip && \
    unzip bowtie2-2.3.4.1-source.zip && \
    cd bowtie2-2.3.4.1&& \
    make && \
    cp -r * $BIN
    #rm $SRC/bowtie2-2.3.4.1-linux-x86_64.zip




#~~~~~~~~~~~~~~~~~~~~~~~~~
## Install Nextflow
#~~~~~~~~~~~~~~~~~~~~~~~~~
WORKDIR $SRC
RUN wget -qO- https://get.nextflow.io | bash

## Add Nextflow binary to your PATH:
RUN mv nextflow $BIN/
#export PATH="/home/mbrown/Tools/NextFlow:$PATH"



#~~~~~~~~~~~~~~~~~~~~~~~~~
# multiQC
#~~~~~~~~~~~~~~~~~~~~~~~~~
# install an older version of MultiQC
WORKDIR $SRC
RUN pip install --force-reinstall git+https://github.com/ewels/MultiQC.git@v1.7




#~~~~~~~~~~~~~~~~~~~~~~~~~
# Install the git repo
#~~~~~~~~~~~~~~~~~~~~~~~~~
WORKDIR $SRC
#ARG CACHEBUST=1
RUN git clone https://github.com/bioinfo-pf-curie/nf-VIF.git



RUN apt-get -qq -y remove automake build-essential cmake curl git gcc g++ wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log


ENV USER=trinityctat

