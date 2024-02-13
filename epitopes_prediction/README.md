## Data
<!-- The entire data set with contact number per residue can be downloaded from the [link](https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/dataset.csv.tar.gz).<br /> Dataset contains following columns: -->
The entire data set with contact number per residue is stored at `./data/sema_2.0/`

* **pdb_id** &#8212; identificator in the [PDB database](https://www.rcsb.org)
* **pdb_id_chain** &#8212; identificator in the [PDB database](https://www.rcsb.org) with name of the anigen chain in the PDB structure
* **resi_pos** &#8212; the residue position in the PDB structure
* **resi_name** &#8212; amino acid 3-letter name
* **res_aa** &#8212; amino acid symbol
* **contact_number** &#8212; contact number values calcualted as the number of antibody residues in contact with any atom of antigen residues within the distance radius R1 == 8. Residues between R1 (8) and R2 (16) have a zero contact number. Residues, which located outside R2  distance radius, have a '-100' value. 
* **contact_number_binary** &#8212; epitope mark. binarized contact number values calcualted as the number of antibody residues in contact with any atom of antigen residues within the distance radius R1 == 4.5. If antiden resudue has 1 or more contacts, it is сonsidered as epitope and label is equal to 1. Residues between R1 (4.5) and R2 (16) are considered to be on epitopes and have label equal to 0. Feather residues have masked label equal to -100. In case of masked test we use only 0 and 1 labeled resudues to check quality of the models. For unmasked test we treat -100 labeled residues as non epitopes with 0 label for quality metrics calculation. 

In case you will be training SEMA-3D model, you will need to download additional data:<br/>
An [original set](https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_2.0/data/data_pdbs.zip) of processed pdb-files.

You can generate your own dataset with different R1 and R2 using scripts in the *dataset_generation* directory.

The directory *data* contains example of training and test sets and example of pdb-file for SEMA-3D inference.

## Training the models
For trainig the model you can use Jupyter Notebooks [SEMA-1D_finetuning](?) or [SEMA-3D_finetuning](https://github.com/AIRI-Institute/SEMAi/blob/main/epitopes_prediction/SEMA_3D/train.ipynb).

## Inference
Prepare your model weights or download ours:<br />
1. Create direcrory *models*: `mkdir models`.
2. Go into the directory: `cd models`
4. Download weights for SEMA-1D:
<!-- ```
wget -O sema_1d_ft_cn_atom_r1_8.0_r2_16.0_0.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_1d_ft_cn_atom_r1_8.0_r2_16.0_0.pth 
wget -O sema_1d_ft_cn_atom_r1_8.0_r2_16.0_1.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_1d_ft_cn_atom_r1_8.0_r2_16.0_1.pth
wget -O sema_1d_ft_cn_atom_r1_8.0_r2_16.0_2.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_1d_ft_cn_atom_r1_8.0_r2_16.0_2.pth
wget -O sema_1d_ft_cn_atom_r1_8.0_r2_16.0_3.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_1d_ft_cn_atom_r1_8.0_r2_16.0_3.pth
wget -O sema_1d_ft_cn_atom_r1_8.0_r2_16.0_4.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_weights/sema_1d_ft_cn_atom_r1_8.0_r2_16.0_4.pth
   ``` -->
   
or SEMA-3D:
```
wget -O sema_3d_0.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_2.0/SEMA_3D_weigths/newdata_sema_saprot_continous_noncut_0.pth
wget -O sema_3d_1.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_2.0/SEMA_3D_weigths/newdata_sema_saprot_continous_noncut_1.pth
wget -O sema_3d_2.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_2.0/SEMA_3D_weigths/newdata_sema_saprot_continous_noncut_2.pth
wget -O sema_3d_3.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_2.0/SEMA_3D_weigths/newdata_sema_saprot_continous_noncut_3.pth
wget -O sema_3d_4.pth https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_2.0/SEMA_3D_weigths/newdata_sema_saprot_continous_noncut_4.pth
```
 

Next, run inference using Jupyter Notebooks [SEMA-1D_inference](?) or [SEMA-3D_inference](https://github.com/AIRI-Institute/SEMAi/blob/main/epitopes_prediction/SEMA_3D/evaluate.ipynb).


## dataset generation
To generate entire dataset run `./dataset_generation/prepare_SEMA_dataset.ipynb` 
