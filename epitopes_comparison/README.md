# Epitope-antigen comparison model

## Training
To run training you must first unzip the dataset files in `small_dataset.tar.gz` folder by running:
```
 tar -xvf small_dataset.tar.gz
```
You can run training by running:
```
CUDA_VISIBLE_DEVICES=0 python train.py
```

## Inference
- Download pre-trained model weigths: `wget -O best_model.pt https://bioinformatics-kardymon.obs.ru-moscow-1.hc.sbercloud.ru/SEMA_2.0/epitopes_comparison/best_model.pt` 
- For inference using a pre-trained model run `epi_compare_inference.ipynb`

## Folder overview
* `small_dataset.tar.gz` contains zipped files for the dataset the model was trained on

    * `small_dataset.pkl` contains pairs of antigen and epitope ids and their alignment matrices

    * `small_epitopes_structures.csv`  and `small_full_structures.csv` contain protein sequences, foldseek sequences and SaProt input structural sequences of the epitopes and antigens respectively

* `train.py` contains core code to train the model
* `saprot_epi_comparison_model.py` contains model class
* `saprot_dataset.py` contains dataset class
<!-- * `best_model.pt` contains weights of the best-performing model -->
* `epi_compare_inference.ipynb` is a jupyter notebook with tutorial on how to use a pre-trained model for inference
* `inference_examples/` folder contains examples used in the tutorial
* `dataset_genetration/` code and additinal information to create dataset from set of pdbs
