#!/bin/bash

###############################################################################
## 1. Set env
###############################################################################

SCRIPT_PATH=$(realpath $0)
SCRIPT_DIRECTORY=$(realpath $(dirname $0))
SEDNA_REPOSITORY_DIRECTORY=$(realpath ${SCRIPT_DIRECTORY}/../)
DATABASE_DIR="${HOME}/.local/share/"

check_exit_status() {
    if [[ $1 -ne 0 ]]; then
        echo "$2"
        exit $1
    fi
}

check_command() {
    command -v "$1" &> /dev/null || { echo "Error: $1 is not installed. Please install it before running the script."; exit 1; }
}

###############################################################################
## 2. Sanity checks
###############################################################################

check_command mamba

###############################################################################
## 3. Remove environments
###############################################################################

for ENV_YAML in "${SEDNA_REPOSITORY_DIRECTORY}/install/environments/sedna_"*.yml; do

    
    ENV_NAME=$(basename "${ENV_YAML}" .yml)
    ENV_PATH=$(conda env list | egrep "${ENV_NAME}" | awk '{print $NF}')
    if [[ -d "${ENV_PATH}" ]]; then
    
        echo "Removing conda environment: ${ENV_NAME} ..."
        mamba env remove -y -n "${ENV_NAME}"
        check_exit_status $? "Removing ${ENV_NAME} failed."

    fi
    
done    

###############################################################################
## 4. Remove databases
###############################################################################

echo "Removing databases ..."
if [[ -d "${DATABASE_DIR}/sedna" ]]; then
    rm -r "${DATABASE_DIR}/sedna"
    check_exit_status $? "Removing ${DATABASE_DIR}/sedna failed."
fi

if [[ -d "${DATABASE_DIR}/antismash" ]]; then 
    rm -r "${DATABASE_DIR}/antismash"
    check_exit_status $? "Removing ${DATABASE_DIR}/antismash failed."
fi

if [[ -d "${DATABASE_DIR}/pfam" ]]; then
    rm -r "${DATABASE_DIR}/pfam"
    check_exit_status $? "Removing ${DATABASE_DIR}/pfam failed."
fi

if [[ -d "${DATABASE_DIR}/deepbgc" ]]; then 
    rm -r "${DATABASE_DIR}/deepbgc"
    check_exit_status $? "Removing ${DATABASE_DIR}/deepbgc failed."
fi    

###############################################################################
# 5. Exit installation script
###############################################################################

echo -e "Uninstallation completed successfully."
