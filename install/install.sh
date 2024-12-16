#!/bin/bash
SCRIPT_PATH=$(realpath $0)
SCRIPT_DIRECTORY=$(realpath $(dirname $0))
SEDNA_REPOSITORY_DIRECTORY=$(realpath ${SCRIPT_DIRECTORY}/../)

LOG_DIRECTORY=${1:-"logs/sedna_installation"}
mkdir -pv ${LOG_DIRECTORY}

CONDA_ENVS_PATH=${2:-"$(conda info --base)/envs/"}

###############################################################################
# Update permissions
###############################################################################

echo "Updating permissions for scripts in ${SEDNA_REPOSITORY_DIRECTORY}/bin"
chmod 775 "${SEDNA_REPOSITORY_DIRECTORY}/bin/sedna.sh"
chmod 775 "${SEDNA_REPOSITORY_DIRECTORY}/bin/"*.py

###############################################################################
# Main environemnt
###############################################################################

echo "Creating sedna main environment"

ENV_NAME="sedna_main_env"
conda create -y -p ${CONDA_ENVS_PATH}/${ENV_NAME} -c conda-forge -c bioconda &> ${LOG_DIRECTORY}/sedna.log

EXIT_STATUS=$?
if [[ ${EXIT_STATUS} -ne 0 ]]; then
    echo "conda create failed. See logs/sedna_installation/sedna.log for details."
    exit ${EXIT_STATUS}
fi

# Copy bin directory to sedna main environment
echo -e "\t*Copying sedna bin into ${ENV_NAME} environment path."
cp -r "${SEDNA_REPOSITORY_DIRECTORY}/bin" "${CONDA_ENVS_PATH}/${ENV_NAME}/" &> ${LOG_DIRECTORY}/sedna.log

EXIT_STATUS=$?
if [[ ${EXIT_STATUS} -ne 0 ]]; then
    echo "Copying sedna bin into ${ENV_NAME} environment path. See logs/sedna_installation/sedna.log for details."
    exit ${EXIT_STATUS}
fi

# Version
cp ${SEDNA_REPOSITORY_DIRECTORY}/bin/VERSION.txt ${CONDA_ENVS_PATH}/${ENV_NAME}/bin/VERSION.txt &> ${LOG_DIRECTORY}/sedna.log

EXIT_STATUS=$?
if [[ ${EXIT_STATUS} -ne 0 ]]; then
    echo "Copying sedna version file into ${ENV_NAME} environment path failed. See logs/sedna_installation/sedna.log for details."
    exit ${EXIT_STATUS}
fi

###############################################################################
# Module environemnt
###############################################################################

for ENV_YAML in ${SEDNA_REPOSITORY_DIRECTORY}/install/environments/sedna_*.yml; do
    # Get environment name
    ENV_NAME=$(basename $ENV_YAML .yml)

    # Create conda environment
    echo "Creating ${ENV_NAME} module environment"
    conda create -f "${ENV_YAML}" -p "${CONDA_ENVS_PATH}/${ENV_NAME}" &> "${LOG_DIRECTORY}/${ENV_NAME}.log"
    
    EXIT_STATUS=$?
    if [[ ${EXIT_STATUS} -ne 0 ]]; then
        echo "Creating ${ENV_NAME} module environment. See logs/sedna_installation/sedna.log for details."
        exit ${EXIT_STATUS}
    fi

    # Copy over files to environment bin/
    echo -e "\t*Copying sedna modules into ${ENV_NAME} environment path"
    cp -r "${SEDNA_REPOSITORY_DIRECTORY}/bin/"*.py "${CONDA_ENVS_PATH}/${ENV_NAME}/bin/" &> "${LOG_DIRECTORY}/${ENV_NAME}.log"
    
    EXIT_STATUS=$?
    if [[ ${EXIT_STATUS} -ne 0 ]]; then
        echo "Copying sedna modules into ${ENV_NAME} environment path. See logs/sedna_installation/sedna.log for details."
        exit ${EXIT_STATUS}
    fi

    # Version
    cp ${SEDNA_REPOSITORY_DIRECTORY}/bin/VERSION.txt ${CONDA_ENVS_PATH}/${ENV_NAME}/bin/VERSION.txt &> ${LOG_DIRECTORY}/sedna.log
    
    EXIT_STATUS=$?
    if [[ ${EXIT_STATUS} -ne 0 ]]; then
        echo "Copying sedna version file into ${ENV_NAME} environment path failed. See logs/sedna_installation/sedna.log for details."
        exit ${EXIT_STATUS}
    fi

done

echo -e "..............................."
echo -e "     Installation Complete     "
echo -e "..............................."