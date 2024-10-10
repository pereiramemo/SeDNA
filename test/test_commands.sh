###############################################################################
### 0. Set env
###############################################################################

REPO_DIR="/home/ec2-user/SageMaker/efs/sandbox/sandbox/development/epereira/dev/bioprospecting_internal/"
DATA_DIR="${REPO_DIR}/test/data"

###############################################################################
### 1. Gather samples
###############################################################################

SAMPLES=$(ls "${DATA_DIR}/OceanDNA"-*.fa)

###############################################################################
### 2. Annotate BGCs
###############################################################################

for s in ${SAMPLES}; do

      SAMPLE_NAME=$(basename "${s}" .fa)
      echo "${SAMPLE_NAME}"
     "${REPO_DIR}/bin/sedna.sh" --module run_all \
     --params "--input_sample ${s} \
               --sample_name ${SAMPLE_NAME} \
               --output_dir ${REPO_DIR}/test/output/sedna_output";
     
done   

###############################################################################
### 3. Dereplicate BGCs
###############################################################################

for s in ${SAMPLES}; do

      SAMPLE_NAME=$(basename "${s}" .fa)
      echo "${SAMPLE_NAME}"
           "${REPO_DIR}/bin/sedna.sh" --module dereplicate \
           --params "--input_dir ${REPO_DIR}/test/output/sedna_output/${SAMPLE_NAME}/bgc_annot/inter"
           
done                     

###############################################################################
### 4. Annotate CDSs
###############################################################################

for s in ${SAMPLES}; do

      SAMPLE_NAME=$(basename "${s}" .fa)
      echo "${SAMPLE_NAME}"
           "${REPO_DIR}/bin/sedna.sh" --module annot_cds \
           --params "--input_dir ${REPO_DIR}/test/output/sedna_output/${SAMPLE_NAME}/bgc_annot/sorted/dereplicated"
done                     

###############################################################################
### 5. Cluster
###############################################################################

"${REPO_DIR}/bin/sedna.sh" --module cluster \
    --params "--input_dir ${REPO_DIR}/test/output/sedna_output/ \
              --threshold 3 \
              --output_tsv ${REPO_DIR}/test/output/sedna_output/clust.tsv"
               




