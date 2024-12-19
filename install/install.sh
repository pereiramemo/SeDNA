#!/bin/bash

###############################################################################
# 1. Set env
###############################################################################

SCRIPT_PATH=$(realpath $0)
SCRIPT_DIRECTORY=$(realpath $(dirname $0))
SEDNA_REPOSITORY_DIRECTORY=$(realpath ${SCRIPT_DIRECTORY}/../)

LOG_DIRECTORY="logs/sedna_installation"
LOG_FILE="${LOG_DIRECTORY}/sedna_install.log"
CONDA_ENVS_PATH=${1:-"$(conda info --base)/envs/"}

check_exit_status() {
    if [[ $1 -ne 0 ]]; then
        echo "$2"
        echo "See ${LOG_FILE} for details."
        exit $1
    fi
}

check_command() {
    command -v "$1" &> /dev/null || { echo "Error: $1 is not installed. Please install it before running the script."; exit 1; }
}

# Initialize Conda
eval "$(conda shell.bash hook)"

###############################################################################
# 2. Sanity checks
###############################################################################

check_command conda
check_command mamba

if [[ ! -d "${LOG_DIRECTORY}" ]]; then
  mkdir -pv  "${LOG_DIRECTORY}"
  check_exit_status $? "Creating ${LOG_DIRECTORY} failed."
fi

if  [[ -f "${LOG_FILE}" ]]; then
  rm "${LOG_FILE}"
  check_exit_status $? "Removing ${LOG_FILE} failed."
fi

###############################################################################
# 3. Update permissions
###############################################################################

echo "Updating permissions for scripts in ${SEDNA_REPOSITORY_DIRECTORY}/bin"
chmod 775 "${SEDNA_REPOSITORY_DIRECTORY}/bin/sedna.sh" &>> "${LOG_FILE}"
check_exit_status $? "Updating permissions for ${SEDNA_REPOSITORY_DIRECTORY}/bin/sedna.sh failed"
chmod 775 "${SEDNA_REPOSITORY_DIRECTORY}/bin/"*.py &>> "${LOG_FILE}"
check_exit_status $? "Updating permissions for ${SEDNA_REPOSITORY_DIRECTORY}/bin/*.py failed"

###############################################################################
# 4. Create main environemnt
###############################################################################

echo "Creating SeDNA main environment"

ENV_NAME="sedna_main_env"
mamba create -y -p ${CONDA_ENVS_PATH}/${ENV_NAME} -c conda-forge -c bioconda &>> "${LOG_FILE}"
check_exit_status $? "Creating conda envrionment ${CONDA_ENVS_PATH}/${ENV_NAME} failed."

# Copy bin directory to sedna main environment
echo -e "\t*Copying sedna bin into ${ENV_NAME} environment path."
cp -r "${SEDNA_REPOSITORY_DIRECTORY}/bin" "${CONDA_ENVS_PATH}/${ENV_NAME}/" &>> "${LOG_FILE}"
check_exit_status $? "Copying ${SEDNA_REPOSITORY_DIRECTORY}/bin failed."

# Copy version file
cp "${SEDNA_REPOSITORY_DIRECTORY}/bin/VERSION.txt" "${CONDA_ENVS_PATH}/${ENV_NAME}/bin/VERSION.txt" &>> "${LOG_FILE}"
check_exit_status $? "Copying ${SEDNA_REPOSITORY_DIRECTORY}/bin/VERSION.txt failed."

###############################################################################
# 5. Module environemnt
###############################################################################

for ENV_YAML in "${SEDNA_REPOSITORY_DIRECTORY}/install/environments/sedna_"*.yml; do

   if [[ ! -f "${ENV_YAML}" ]]; then
        echo "Error: Environment file ${ENV_YAML} does not exist." &>> "${LOG_FILE}"
        exit 1
   fi

    # Get environment name
    ENV_NAME=$(basename "${ENV_YAML}" .yml)
     
    # Create conda environment
    echo "Creating ${ENV_NAME} module environment"
    mamba env create -y -f "${ENV_YAML}" -p "${CONDA_ENVS_PATH}/${ENV_NAME}" &>> "${LOG_FILE}"
    check_exit_status $? "Creating conda environment ${ENV_YAML} failed."
   
    # Make sure bin dir exists
    if [[ ! -f "${CONDA_ENVS_PATH}/${ENV_NAME}/bin/" ]]; then
        mkdir "${CONDA_ENVS_PATH}/${ENV_NAME}/bin/" &>> "${LOG_FILE}"
        check_exit_status $? "Creating ${CONDA_ENVS_PATH}/${ENV_NAME}/bin/ failed."
    fi

    # Copy over files to environment bin/
    echo -e "\t*Copying sedna modules into ${ENV_NAME} environment path"
    cp -r "${SEDNA_REPOSITORY_DIRECTORY}/bin/"*.py "${CONDA_ENVS_PATH}/${ENV_NAME}/bin/" &>> "${LOG_FILE}"
    check_exit_status $? "Copying sedna modules into ${ENV_NAME} failed."
    cp -r "${SEDNA_REPOSITORY_DIRECTORY}/bin/src" "${CONDA_ENVS_PATH}/${ENV_NAME}/bin/" &>> "${LOG_FILE}"
    check_exit_status $? "Copying sedna src code into ${ENV_NAME} failed."
    
    # Copy version file
    cp "${SEDNA_REPOSITORY_DIRECTORY}/bin/VERSION.txt" "${CONDA_ENVS_PATH}/${ENV_NAME}/bin/VERSION.txt" &>> "${LOG_FILE}"
    check_exit_status $? "Copying ${SEDNA_REPOSITORY_DIRECTORY}/bin/VERSION.txt failed."
    
done

###############################################################################
# 6. Exit installation script
###############################################################################

echo -e "Installation completed successfully."
