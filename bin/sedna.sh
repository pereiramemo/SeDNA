#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Initialize variables
CONDA_BASE=$(conda info --base)
CONDA_ENVS_PATH=${CONDA_ENVS_PATH:-"${CONDA_BASE}/envs/"}
BIN_DIRECTORY=$(dirname $0)
PARAMS="--help"

AVAILABLE_MODULES=(
"run_antismash"
"run_deepbgc"
"run_gecco"
"run_all"
"dereplicate"
"annot_cds"
"cluster"
)

usage() {
    echo "Usage: $0 [-m <module>] [-o <options>] [-v|--version] [-h|--help]"
    echo " "
    echo "Options:"
    echo "  -m, --module    Specify the module. Available modules: ${AVAILABLE_MODULES[*]}"
    echo "  -p, --params    Specify parameters to give to each module"
    echo "  -v, --version   Display the version information"
    echo "  -h, --help      Display this help message"
    exit 0
}

# Parse command-line arguments
OPTIONS=$(getopt -o m:p:vVh --long module:,params:,version,help -n "$0" -- "$@")
if [ $? -ne 0 ]; then
  usage
fi

# Set the positional parameters to the parsed options
eval set -- "$OPTIONS"

# Parse options
while true; do
    case "$1" in
        -m|--module)
            MODULE="$2"
            shift 2
            ;;
        -p|--params)
            PARAMS="$2"
            shift 2
            ;;
        -v|--version)
            echo "Version:"
            cat "${BIN_DIRECTORY}/sedna_version.txt"
            exit 0
            ;;
        -h|--help)
            usage
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if module is valid and non empty
TEST_MODULE_FLAG=0
for AVAILABLE_MODULE in "${AVAILABLE_MODULES[@]}"; do
    if [[ "${AVAILABLE_MODULE}" == "${MODULE}" ]]; then
        TEST_MODULE_FLAG=1
        break
    fi        
done    
if [[ "${TEST_MODULE_FLAG}" == 0 ]]; then
    echo "Invalid or empty module name. Must be one of: ${AVAILABLE_MODULES[*]}"
    usage
    exit 1
fi

# Execute specified modules
case "${MODULE}" in
    "run_antismash")
        source "${CONDA_BASE}/bin/activate" ${CONDA_ENVS_PATH}/antismash_env
        "${BIN_DIRECTORY}/modules/run_antismash.py" $PARAMS
        ;;
    "run_deepbgc")
        source "${CONDA_BASE}/bin/activate" ${CONDA_ENVS_PATH}/deepbgc_env
        "${BIN_DIRECTORY}/modules/run_deepbgc.py" $PARAMS
        ;;
    "run_gecco")
        source "${CONDA_BASE}/bin/activate" ${CONDA_ENVS_PATH}/gecco_env
        "${BIN_DIRECTORY}/modules/run_gecco.py" $PARAMS
        ;;
    "run_all")
        source "${CONDA_BASE}/bin/activate" ${CONDA_ENVS_PATH}/deepbgc_env
        "${BIN_DIRECTORY}/modules/run_deepbgc.py" $PARAMS
        source "${CONDA_BASE}/bin/activate" ${CONDA_ENVS_PATH}/antismash_env
        "${BIN_DIRECTORY}/modules/run_antismash.py" $PARAMS
        source "${CONDA_BASE}/bin/activate" ${CONDA_ENVS_PATH}/gecco_env
        "${BIN_DIRECTORY}/modules/run_gecco.py" $PARAMS
        ;;
    "dereplicate")
        source "${CONDA_BASE}/bin/activate" ${CONDA_ENVS_PATH}/sedna_env
        "${BIN_DIRECTORY}/modules/dereplicate.py" $PARAMS
        ;;
    "annot_cds")
        source "${CONDA_BASE}/bin/activate" ${CONDA_ENVS_PATH}/sedna_env
        "${BIN_DIRECTORY}/modules/annot_cds.py" $PARAMS
        ;;    
    "cluster")
        source "${CONDA_BASE}/bin/activate" ${CONDA_ENVS_PATH}/sedna_env
        "${BIN_DIRECTORY}/modules/cluster_mean.py" $PARAMS
        ;;
esac        
