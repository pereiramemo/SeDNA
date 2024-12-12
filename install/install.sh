#!/bin/bash
# __version__ = "0.1.0"
SCRIPT_PATH=$(realpath $0)
SCRIPT_DIRECTORY=$(realpath $(dirname $0))
SEDNA_REPOSITORY_DIRECTORY=$(realpath ${SCRIPT_DIRECTORY}/../)

LOG_DIRECTORY=${1:-"logs/sedna_installation"}
mkdir -pv ${LOG_DIRECTORY}

CONDA_ENVS_PATH=${2:-"$(conda info --base)/envs/"}

# Update permissions
echo "Updating permissions for scripts in ${SEDNA_REPOSITORY_DIRECTORY}/bin"
chmod 775 ${SEDNA_REPOSITORY_DIRECTORY}/bin/sedna.sh
chmod 775 ${SEDNA_REPOSITORY_DIRECTORY}/bin/*.py
