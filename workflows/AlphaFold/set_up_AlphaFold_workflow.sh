CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
HEMEBINDER_DIR="${CARPNN_DIR}/public/heme_binder_diffusion"

## Clone the repo
git clone https://github.com/ikalvet/heme_binder_diffusion.git ${HEMEBINDER_DIR}

(
    cd "${HEMEBINDER_DIR}" && \
    # 1. CORRECTLY REMOVE the unwanted submodule using git rm 
    # (assuming this was fixed from the previous issue)
    git rm lib/LigandMPNN && \
    
    # 2. COMMIT the removal
    git commit -m "Remove LigandMPNN submodule requirement" && \

    # 3. FIX: Change the submodule URL from SSH to HTTPS
    # Use 'git config' to modify the URL for 'lib/alphafold' in the .gitmodules file.
    git config -f .gitmodules submodule.lib/alphafold.url https://github.com/google-deepmind/alphafold.git && \
    
    # 4. COMMIT the URL change
    git add .gitmodules && \
    git commit -m "Update alphafold submodule URL to HTTPS" && \

    # 5. Now perform Submodule update
    git submodule update --init --recursive 
)

## Down load AF2 model weights
## IF you already downloaded these weights else where you can also skip this step and place them here accordingly
#### AF2 model weights
(
    cd "${HEMEBINDER_DIR}/lib/alphafold" && \
    mkdir -p model_weights/params && cd model_weights/params && \
    wget https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar && \
    tar --extract --verbose --file=alphafold_params_2021-07-14.tar
)

