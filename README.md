# SEMA
SEMA (Spatial Epitope Modelling with Artificial intelligence) is a tool for conformational B-cell eptiope prediction from the primary protein sequence or tertiary structure. SEMA involves the use of sequence-based (SEMA-1D) and structure-based (SEMA-3D) approaches. SEMA-1D model is based on an ensemble of [Esm-1v](https://github.com/facebookresearch/esm) transformer deep neural network protein language models. SEMA-3D model is based on an ensemble of inverse folding models, [Esm-IF1](https://github.com/facebookresearch/esm). Both models were fine-tuned to predict the antigen interaction propensity of the amino acid (AA) residue with Fab regions of immunoglobulins. SEMA provides an interpretable score indicating the log-scaled expected number of contacts with antibody residues. 

SEMA is also availble via [web-interface](http://sema.airi.net/).


### Disclaimer:
This code is provided under MIT License

## Dataset
The whole dataset presents as an archived csv-tabel in data folder (dataset.csv.tar.gz). Dataset contains following columns:

* **pdb_id** &#8212; identificator in the [PDB database](https://www.rcsb.org)
* **resi** &#8212; the residue position in the PDB structure
* **resin** &#8212; amino acid 3-letter name
* **res_aa** &#8212; amino acid symbol
* **anigen_chain** &#8212; name of the anigen chain in the PDB structure
* **fab_chains** &#8212; names of the antibody chains in the PDB structure
* **contact_number_R1=i_R2=j** &#8212; contact number values calcualted as the number of antibody residues in contact with any atom of antigen residues within the distance radius R1. Residues between R1 and R2 have a zero contact number. Residues, which located outside R2  distance radius, have a '-100' value. 


## Training the models
For trainig the model you can use Jupyter Notebooks [SEMA-1D_finetuning](https://github.com/AIRI-Institute/SEMAi/blob/main/SEMA_1D/SEMA-1D_finetuning.ipynb) or [SEMA-3D_finetuning](https://github.com/AIRI-Institute/SEMAi/blob/main/SEMA_3D/SEMA-3D_finetuning.ipynb).

## Inference
Prepare your model weights or download ours:<br />
1. Create direcrory "models": `mkdir models`.
2. Go into directory: `cd models`
4. Download weights for SEMA-1D:
```wget -O sema_1d_ft_cn_atom_r1_8.0_r2_16.0_0.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_1d_ft_cn_atom_r1_8.0_r2_16.0_0.pth 
   wget -O sema_1d_ft_cn_atom_r1_8.0_r2_16.0_1.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_1d_ft_cn_atom_r1_8.0_r2_16.0_1.pth
   wget -O sema_1d_ft_cn_atom_r1_8.0_r2_16.0_2.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_1d_ft_cn_atom_r1_8.0_r2_16.0_2.pth
   wget -O sema_1d_ft_cn_atom_r1_8.0_r2_16.0_3.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_1d_ft_cn_atom_r1_8.0_r2_16.0_3.pth
   wget -O sema_1d_ft_cn_atom_r1_8.0_r2_16.0_4.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_1d_ft_cn_atom_r1_8.0_r2_16.0_4.pth
   ```
   
or SEMA-3D:
```
wget -O sema_3d_cn_atom_r1_8.0_r2_18.0_0.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_3d_cn_atom_r1_8.0_r2_18.0_0.pt
wget -O sema_3d_cn_atom_r1_8.0_r2_18.0_1.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_3d_cn_atom_r1_8.0_r2_18.0_1.pt 
wget -O sema_3d_cn_atom_r1_8.0_r2_18.0_2.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_3d_cn_atom_r1_8.0_r2_18.0_2.pt 
wget -O sema_3d_cn_atom_r1_8.0_r2_18.0_3.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_3d_cn_atom_r1_8.0_r2_18.0_3.pt 
wget -O sema_3d_cn_atom_r1_8.0_r2_18.0_4.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_3d_cn_atom_r1_8.0_r2_18.0_4.pt
```
   
6. Unarchive files: `tar -cvzf .tar.gz`

Next, run inference using Jupyter Notebooks [SEMA-1D_inference](https://github.com/AIRI-Institute/SEMAi/blob/main/SEMA_1D/SEMA-1D_inference.ipynb) or [SEMA-3D_inference](https://github.com/AIRI-Institute/SEMAi/blob/main/SEMA_3D/SEMA-3D_inference.ipynb).



