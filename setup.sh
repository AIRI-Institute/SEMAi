mamba create -n sema_env python==3.9
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

mamba install -y ipykernel
ipython kernel install --user --name=sema_env

pip install fair-esm
pip install accelerate -U

# mamba install -y torch-scatter
pip install torch-scatter
pip install scikit-plot