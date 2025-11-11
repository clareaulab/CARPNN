## Clone the LigandMPNN repo into the public directory
git clone https://github.com/dauparas/LigandMPNN ../public

## Create the apporpriate ligandmpnn environment
(
    cd ../public/LigandMPNN
    conda create -n ligandmpnn_env python=3.11
    pip3 install -r requirements.txt
)

## download the af2 and mpnn model param weights
(
    cd ../public/LigandMPNN
    bash ./get_model_params.sh "./model_params"
)

## set up done.
echo "set_up_LigandMPNN_workflow.sh finished. Please update the bash script paths."