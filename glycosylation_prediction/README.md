## Data
Download the train dataset:
```
wget -O data/Nglyco_train.pkl https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_2.0/glycosylation_prediction/data/Nglyco_train.pkl
```
Validation and test sets are already in `./data` \
You can also generate datset from own fasta by `loadTrainingSetFasta` function in [SEMA_PTM_finetuning](https://github.com/AIRI-Institute/SEMAi/blob/main/glycosylation_prediction/SEMA_PTM/PTM_finetuning.ipynb)


## Training the SEMA_PTM model
For trainig the model you can use Jupyter Notebook [SEMA_PTM_finetuning](https://github.com/AIRI-Institute/SEMAi/blob/main/glycosylation_prediction/SEMA_PTM/PTM_finetuning.ipynb).

## Inference
Prepare your model weights or download ours:<br />
1. Create directory *models*: `mkdir models`.
2. Go into the directory: `cd models`
4. Download weights:
```
wget -O sema_1d_ft_cn_atom_r1_8.0_r2_16.0_0.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_2.0/glycosylation_prediction/sema2.0_esm2_ptm.pth 
   ```


Next, run inference using Jupyter Notebook [SEMA_PTM_inference](https://github.com/AIRI-Institute/SEMAi/blob/main/glycosylation_prediction/SEMA_PTM/PTM_inference.ipynb).
