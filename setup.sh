mamba create -y -n sema_env_test python==3.10
eval "$(conda shell.bash hook)"
conda activate sema_env_test
mamba install -y pytorch pytorch-cuda=11.8 -c pytorch -c nvidia
mamba install -y -c anaconda biopython=1.78
mamba install -y tqdm
mamba install -y scipy
mamba install -y scikit-learn=1.4
mamba install -y pandas=2.2
mamba install -y transformers=4.37
mamba install -y biotite=0.38
mamba install -y pyg -c pyg

pip install fair-esm
pip install accelerate -U

# mamba install -y torch-scatter
# pip install torch-scatter
pip install scikit-plot==0.3
pip install biopandas==0.4

# # add repository root to env path for saprot_utils usage  
# # export PYTHONPATH=$PYTHONPATH:$(pwd)
# conda config --env --add envs_dirs $(pwd)

mamba install -y ipykernel
ipython kernel install --user --name=sema_env_test