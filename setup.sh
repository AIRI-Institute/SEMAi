mamba create -y -n sema_env python==3.9
eval "$(conda shell.bash hook)"
conda activate sema_env
mamba install -y pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
# mamba install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia
mamba install -y -c anaconda biopython
mamba install -y tqdm
mamba install -y scipy
mamba install -y scikit-learn
mamba install -y pandas
mamba install -y transformers
mamba install -y biotite
mamba install -y pyg -c pyg

pip install fair-esm
pip install accelerate -U

# mamba install -y torch-scatter
pip install torch-scatter
pip install scikit-plot
pip install biopandas

# # add repository root to env path for saprot_utils usage  
# # export PYTHONPATH=$PYTHONPATH:$(pwd)
# conda config --env --add envs_dirs $(pwd)

mamba install -y ipykernel
ipython kernel install --user --name=sema_env