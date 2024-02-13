# SEMA 2.0
SEMA (Spatial Epitope Modelling with Artificial intelligence) is a set of research tools for sequence- and structure-based conformational B-cell eptiope prediction, accurate identification of N-glycosylation sites, and a distinctive module for comparing the structures of antigen B-cell epitopes enhancing our ability to analyze and understand its immunogenic properties.

SEMA is also availble via [web-interface](http://sema.airi.net/).

### Disclaimer:
This code is provided under MIT License

## Instalation
This script creates new environment `sema_env` 
```
./setup.sh
```

## Conformational B-cell eptiope prediction models:
Involves the use of sequence-based (SEMA-1D) and structure-based (SEMA-3D) approaches. SEMA-1D model is based on an ensemble of [ESM2](https://github.com/facebookresearch/esm) transformer deep neural network protein language models. SEMA-3D model is based on an ensemble of inverse folding models, [SaProt](https://github.com/westlake-repl/SaProt). Both models were fine-tuned to predict the antigen interaction propensity of the amino acid (AA) residue with Fab regions of immunoglobulins. SEMA provides an interpretable score indicating the log-scaled expected number of contacts with antibody residues. \
Code and datset and additional README is stored at `./epitopes_prediction/` folder

## N-glycosylation sites prediction model:
The N-glycosylation prediction model was obtained by adding a fully-connected linear layer on the top layer of the ESM-2 pre-trained model. \
Code and datset and additional README is stored at `./glycosylation_prediction/` folder

## Epitope comparison model
The model is trained to identify local structural similarities within proteins, based on the non-linear transformation of multiplication of the embeddings of PLM with geometric modalities. \
Code and datset and additional README is stored at `./epitopes_comparison/` folder

# Foldseek
We used [SaProt](https://github.com/westlake-repl/SaProt) as base model for SEMA-3D and epitope comparison model. This model utilizes 3Di tokens, prodused by [Foldseek](https://github.com/steineggerlab/foldseek). To create 3Di tokens we used `foldseek_util.py` script from [SaProt](https://github.com/westlake-repl/SaProt?tab=readme-ov-file#convert-protein-structure-into-structure-aware-sequence) foldseek binary file. These files are stored at `./saprot_utils` folder 



