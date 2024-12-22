################################################################################
# SeDNA Dockerfile
################################################################################

################################################################################
# Define and configure the container
################################################################################

# Use the latest Ubuntu LTS as the base image
FROM ubuntu:22.04

# Use the latest Ubuntu LTS as the base image
RUN mkdir -p /bioinfo/software/miniconda3
RUN mkdir -p /bioinfo/resources/
RUN mkdir -p /output/
RUN mkdir -p /input/
RUN mkdir -p /scratch/

# Update the system and install prerequisites
RUN apt-get update && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    bzip2 \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

################################################################################
# Install Miniconda (latest version)
################################################################################
 
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
 -O /bioinfo/software/miniconda3/miniconda.sh \
 --no-check-certificate
RUN chmod +x /bioinfo/software/miniconda3/miniconda.sh 
RUN bash /bioinfo/software/miniconda3/miniconda.sh -b -u -p /bioinfo/software/miniconda3/
RUN rm /bioinfo/software/miniconda3/miniconda.sh;

################################################################################
# Set envs
################################################################################

# Add Conda to the PATH
ARG PATH="/bioinfo/software/miniconda3/bin/:$PATH"
ENV PATH="/bioinfo/software/miniconda3/bin/:$PATH"

################################################################################
# Insall mamba
################################################################################

RUN conda install -n base -c conda-forge mamba

################################################################################
# Verify installations
################################################################################

RUN conda --version && mamba --version && git --version

################################################################################
# Copy SeDNA directory
################################################################################

ADD SeDNA.tar /bioinfo/software/

################################################################################
# Install SeDNA
################################################################################

RUN cd /bioinfo/software/SeDNA; \
./install/install.sh; \
./install/get_databases.sh
