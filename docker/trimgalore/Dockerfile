# Base image for trimgalore/0.6.10,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

# Depedencies of trimgalore/0.6.10:
#  - perl: Larry Wall's Practical Extraction and Report Language
#  - cutadapt: Clean biological sequences from high-throughput sequencing reads
#  - fastqc: quality control for high throughput sequence data
#  - pigz: Parallel Implementation of GZip

LABEL maintainer="Skyler Kuhn <kuhnsa AT nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="trimgalore" \
    about.summary="Docker image for trimgalore/0.6.10" \
    about.home="https://github.com/FelixKrueger/TrimGalore" \
    about.documentation="https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md" \
    about.tags="Genomics"

############### INIT ################
# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with the host's
# filesystem, i.e. /opt and /data
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone

############### SETUP ################
# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        apt-transport-https \
        autotools-dev \
        automake \
        autoconf \
        build-essential \
        ca-certificates \
        cpanminus \
        curl \
        gawk \
        git \
        grep \
        gzip \
        locales \
        pigz \
        unzip \
        wget \
        # Installs cutadapt/3.5-1build1
        cutadapt=3.5-1build1 \
        perl \
        fastqc \
        python3 \
        python3-pip \
        unzip \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install Python Packages 
# and install multiqc/1.14
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip \
	&& pip3 install argparse \
    && pip3 install matplotlib \
    && pip3 install pandas \
    && pip3 install scipy

# Perl fix issue
RUN cpanm FindBin Term::ReadLine LWP::Simple


############### MANUAL ################
# Install trimgalore/0.6.10, for more info see:
# https://github.com/FelixKrueger/TrimGalore
# Requires: perl, cutadapt, fastqc, pigz
RUN mkdir -p /opt2/trimgalore/0.6.10/ \
    && wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -P /opt2/trimgalore/0.6.10/ \
    && tar -vxzf /opt2/trimgalore/0.6.10/0.6.10.tar.gz -C /opt2/trimgalore/0.6.10/ \
    && rm /opt2/trimgalore/0.6.10/0.6.10.tar.gz
ENV PATH="${PATH}:/opt2/trimgalore/0.6.10/TrimGalore-0.6.10"


################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/trimgalore.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH"
WORKDIR /data2