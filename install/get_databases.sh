#!/bin/bash

###############################################################################
# 1. Set env
###############################################################################

SCRIPT_PATH=$(realpath $0)
SCRIPT_DIRECTORY=$(realpath $(dirname $0))
SEDNA_REPOSITORY_DIRECTORY=$(realpath ${SCRIPT_DIRECTORY}/../)

LOG_DIRECTORY="logs/sedna_installation"
LOG_FILE="${LOG_DIRECTORY}/get_databases.log"
DATABASE_DIR="${HOME}/.local/share/"

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

check_command wget
check_command conda

if [[ ! -d "${LOG_DIRECTORY}" ]]; then
  mkdir -pv  "${LOG_DIRECTORY}"
  check_exit_status $? "Creating ${LOG_DIRECTORY} failed."
fi  

if  [[ -f "${LOG_FILE}" ]]; then
  rm "${LOG_FILE}"
  check_exit_status $? "Removing ${LOG_FILE} failed."
fi

###############################################################################
# 3. Create the database directory if it doesn't exist
###############################################################################

if [[ ! -d "${DATABASE_DIR}" ]]; then
     
    mkdir -pv "${DATABASE_DIR}" 2>> "${LOG_FILE}"
    check_exit_status $? "Creating ${DATABASE_DIR} failed."
  
    mkdir -pv "${DATABASE_DIR}/sedna/embeddings" 2>> "${LOG_FILE}"
    check_exit_status $? "Creating ${DATABASE_DIR}/sedna/embeddings failed."

    mkdir -pv "${DATABASE_DIR}/pfam" 2>> "${LOG_FILE}"
    check_exit_status $? "Creating ${DATABASE_DIR}/pfam failed."

    mkdir -pv "${DATABASE_DIR}/antismash" 2>> "${LOG_FILE}"
    check_exit_status $? "Creating ${DATABASE_DIR}/antismash failed."
    
    mkdir -pv "${DATABASE_DIR}/deepbgc" 2>> "${LOG_FILE}"
    check_exit_status $? "Creating ${DATABASE_DIR}/deepbgc failed."

fi

###############################################################################
# 4. Download and unzip the Pfam database
###############################################################################

PFAM_URL="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

if [[ ! -f "${DATABASE_DIR}/Pfam-A.hmm" ]]; then

    echo -e "\t*Downloading Pfam database from ${PFAM_URL}..."
    wget -O "${DATABASE_DIR}/pfam/Pfam-A.hmm.gz" "${PFAM_URL}"  &>> "${LOG_FILE}"
    check_exit_status $? "Downloading Pfam database from ${PFAM_URL} failed."

    echo -e "\t*Unzipping Pfam database..."
    gzip -d "${DATABASE_DIR}/pfam/Pfam-A.hmm.gz"  &>> "${LOG_FILE}"
    check_exit_status $? "Unzipping Pfam database failed."

fi

###############################################################################
# 5. Copy the embeddings and frequency pkl files
###############################################################################

echo -e "\t*Copying sedna embedding files to ${DATABASE_DIR}/sedna/embeddings"
cp -r "${SEDNA_REPOSITORY_DIRECTORY}/resources/embeddings"/*.pkl "${DATABASE_DIR}/sedna/embeddings/" &>> "${LOG_FILE}"
check_exit_status $? "Copying sedna pkl files failed."

###############################################################################
# 6. Download tool specific databases: antismash
###############################################################################

echo -e "\t*Downloading antiSMASH database"
SEDNA_ANTISMASH_ENV=$(conda env list | grep "^sedna_antismash_env " | awk '{print $NF}')
conda activate "${SEDNA_ANTISMASH_ENV}" &>> "${LOG_FILE}"
check_exit_status $? "Activating ${SEDNA_ANTISMASH_ENV} failed."

download-antismash-databases --database-dir "${DATABASE_DIR}/antismash" &>> "${LOG_FILE}"
EXIT_STATUS=$?
conda deactivate 
check_exit_status "${EXIT_STATUS}" "Downloading antismash databases failed."

###############################################################################
# 7. Download tool specific databases: deepbgc
###############################################################################

echo -e "\t*Downloading deepBGC database"
SEDNA_DEEPBGC_ENV=$(conda env list | grep "^sedna_deepbgc_env " | awk '{print $NF}')
conda activate "${SEDNA_DEEPBGC_ENV}" &>> "${LOG_DIRECTORY}/get_databases.log"
check_exit_status $? "Activating ${SEDNA_DEEPBGC_ENV} failed."

deepbgc download &>> "${LOG_DIRECTORY}/get_databases.log"
EXIT_STATUS=$?
conda deactivate 
check_exit_status "${EXIT_STATUS}" "Downloading deepbgc databases failed."

