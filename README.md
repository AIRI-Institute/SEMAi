# SEMA 2.0
SEMA (Spatial Epitope Modelling with Artificial intelligence) is a set of research tools for sequence- and structure-based conformational B-cell eptiope prediction, accurate identification of N-glycosylation sites, and a distinctive module for comparing the structures of antigen B-cell epitopes enhancing our ability to analyze and understand its immunogenic properties.

SEMA is also availble via [web-interface](http://sema.airi.net/).

## Instalation
This script creates new environment `sema_env` 
```
./setup.sh
```
## Models
### Conformational B-cell eptiope prediction models:
Involves the use of sequence-based (SEMA-1D) and structure-based (SEMA-3D) approaches. SEMA-1D model is based on an ensemble of [ESM2](https://github.com/facebookresearch/esm) transformer deep neural network protein language models. SEMA-3D model is based on an ensemble of inverse folding models, [SaProt](https://github.com/westlake-repl/SaProt). Both models were fine-tuned to predict the antigen interaction propensity of the amino acid (AA) residue with Fab regions of immunoglobulins. SEMA provides an interpretable score indicating the log-scaled expected number of contacts with antibody residues. \
Code, datset and additional README is stored at `./epitopes_prediction/` folder

### N-glycosylation sites prediction model:
The N-glycosylation prediction model (SEMA_PTM) was obtained by adding a fully-connected linear layer on the top layer of the ESM-2 pre-trained model. \
Code, datset with additional README is stored at `./glycosylation_prediction/` folder

### Epitope comparison model
The model is trained to identify local structural similarities within proteins, based on the non-linear transformation of multiplication of the embeddings of PLM with geometric modalities. \
Code, datset and additional README is stored at `./epitopes_comparison/` folder

## Notes

### Foldseek
We used [SaProt](https://github.com/westlake-repl/SaProt) as base model for SEMA-3D and epitope comparison model. This model utilizes 3Di tokens, prodused by [Foldseek](https://github.com/steineggerlab/foldseek). To create 3Di tokens we used `foldseek_util.py` script from [SaProt](https://github.com/westlake-repl/SaProt?tab=readme-ov-file#convert-protein-structure-into-structure-aware-sequence) foldseek binary file. These files are stored at `./saprot_utils` folder 

### Disclaimer:
This code is provided under MIT License:

The MIT License (MIT) Copyright (c) 2016 AYLIEN Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.




