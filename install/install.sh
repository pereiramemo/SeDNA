#!/bin/bash
SCRIPT_PATH=$(realpath $0)
SCRIPT_DIRECTORY=$(realpath $(dirname $0))
SEDNA_REPOSITORY_DIRECTORY=$(realpath ${SCRIPT_DIRECTORY}/../)

LOG_DIRECTORY=${1:-"logs/sedna_installation"}
mkdir -pv ${LOG_DIRECTORY}

CONDA_ENVS_PATH=${2:-"$(conda info --base)/envs/"}

# Update permissions
echo "Updating permissions for scripts in ${SEDNA_REPOSITORY_DIRECTORY}/bin"
chmod 775 "${SEDNA_REPOSITORY_DIRECTORY}/bin/sedna.sh"
chmod 775 "${SEDNA_REPOSITORY_DIRECTORY}/bin/modules/"*.py

# Environemnts
# Main environment
echo "Creating sedna main environment"

ENV_NAME="sedna_main_env"
(conda create -y -p ${CONDA_ENVS_PATH}/${ENV_NAME} -c conda-forge -c bioconda || echo "Error when creating sedna main environment" ; exit 1) &> ${LOG_DIRECTORY}/sedna.log

if [[ $? -ne 0 ]]; then
    echo "conda create failed. See logs/sedna_installation/sedna.log for details"
    exit $? 
fi

# Copy main executable
echo -e "\t*Copying main sedna executable into ${ENV_NAME} environment path"
cp -r ${SEDNA_REPOSITORY_DIRECTORY}/bin/sedna.sh ${CONDA_ENVS_PATH}/${ENV_NAME}/bin/

# Copy modules to environment bin/
echo -e "\t*Copying sedna modules into ${ENV_NAME} environment path"
cp -r ${SEDNA_REPOSITORY_DIRECTORY}/bin/modules/*.py ${CONDA_ENVS_PATH}/${ENV_NAME}/bin/

# Copy packages to environment bin/
echo -e "\t*Copying sedna modules into ${ENV_NAME} environment path"
cp -r ${SEDNA_REPOSITORY_DIRECTORY}/bin/src ${CONDA_ENVS_PATH}/${ENV_NAME}/bin/

# Version
cp ${SEDNA_REPOSITORY_DIRECTORY}/VERSION.txt ${CONDA_ENVS_PATH}/${ENV_NAME}/bin/VERSION.txt



