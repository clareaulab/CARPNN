## Clone the LigandMPNN repo into the public directory
CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
git clone https://github.com/dauparas/LigandMPNN ${CARPNN_DIR}/public/LigandMPNN

## Create the apporpriate ligandmpnn environment
(
    cd ${CARPNN_DIR}/public/LigandMPNN && \
    conda create -n ligandmpnn_env python=3.11 && \
    conda activate ligandmpnn_env && \
    pip3 install -r requirements.txt
)

## download the af2 and mpnn model param weights
(
    cd ${CARPNN_DIR}/public/LigandMPNN && \
    bash ./get_model_params.sh "./model_params"
)

## set up done.
echo "set_up_LigandMPNN_workflow.sh finished. Please update the bash script paths."