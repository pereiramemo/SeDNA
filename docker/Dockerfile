# Use the latest Ubuntu LTS as the base image
FROM ubuntu:22.04

# Set environment variables to minimize prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Update the system and install prerequisites
RUN apt-get update && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    bzip2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda (latest version)

RUN mkdir -p ~/miniconda3 \
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh \
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3 \
rm ~/miniconda3/miniconda.sh \

# Add Conda to the PATH
ENV PATH="/opt/miniconda/bin:$PATH"

# Install Mamba using Conda
RUN conda install -y mamba -n base -c conda-forge \
    && conda clean -a -y

# Verify installations
RUN conda --version && mamba --version && git --version

# Set default shell to Bash
CMD ["/bin/bash"]
