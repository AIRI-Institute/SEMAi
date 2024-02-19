import os
# set cuda params
# 'TORCH_HOME'directory will be used to save origenal esm-1v weights
os.environ['TORCH_HOME'] = "../torch_hub"
os.environ['CUDA_VISIBLE_DEVICES'] = "0"

import os
import copy
import math
import json
import pandas as pd
import numpy as np
import pickle
import scipy

from pathlib import Path
from biotite.structure.residues import get_residues

import esm
from esm.data import BatchConverter
from esm.inverse_folding.util import CoordBatchConverter

import torch
from torch.utils.data import Dataset
from torch import nn
from tqdm import tqdm

import transformers
from transformers.modeling_outputs import SequenceClassifierOutput
from transformers import Trainer, TrainingArguments, EvalPrediction

import sklearn
import sklearn.metrics as metrics
from sklearn.metrics import r2_score, mean_squared_error, auc, plot_precision_recall_curve


from pathlib import Path
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import random
import os
import copy
import math
import json
os.environ['CUDA_VISIBLE_DEVICES'] = "0"
from esm.data import BatchConverter
from esm.inverse_folding.util import CoordBatchConverter
import torch
from typing import Sequence, Tuple, List

from torch.utils.data import Dataset
from torch import nn
from tqdm import tqdm
import transformers
from transformers.modeling_outputs import SequenceClassifierOutput
from transformers import Trainer, TrainingArguments, EvalPrediction
import scipy
import pandas as pd
import numpy as np
import pickle
from biotite.structure.residues import get_residues
import esm
import sklearn
from sklearn.metrics import r2_score, mean_squared_error, auc, plot_precision_recall_curve
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
_, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
import biotite.structure
from biotite.structure.io import pdbx, pdb
from biotite.structure.residues import get_residues
from biotite.structure import filter_backbone
from biotite.structure import get_chains
from biotite.sequence import ProteinSequence
import gzip


def get_atom_coords_residuewise(atoms: List[str], struct: biotite.structure.AtomArray):
    """
    Example for atoms argument: ["N", "CA", "C"]
    """
    def filterfn(s, axis=None):
        filters = np.stack([s.atom_name == name for name in atoms], axis=1)
        sum = filters.sum(0)
        if not np.all(sum <= np.ones(filters.shape[1])):
            raise RuntimeError("structure has multiple atoms with same name")
        index = filters.argmax(0)
        coords = s[index].coord
        coords[sum == 0] = float("nan")
        return coords

    return biotite.structure.apply_residue_wise(struct, struct, filterfn)



def extract_coords_from_structure(structure: biotite.structure.AtomArray):
    """
    Args:
        structure: An instance of biotite AtomArray
    Returns:
        Tuple (coords, seq)
            - coords is an L x 3 x 3 array for N, CA, C coordinates
            - seq is the extracted sequence
    """
    coords = get_atom_coords_residuewise(["N", "CA", "C"], structure)
    residue_identities = get_residues(structure)[1]
    seq = ''.join([ProteinSequence.convert_letter_3to1(r) for r in residue_identities])
    return coords, seq


def load_structure(pdb_path, chain = None):#"A"):
    if chain is not None:
        chain = chain.upper()
    structure = _load_structure(str(pdb_path), chain)
    resi_index = get_residues(structure)[0]
    resi_aa    = get_residues(structure)[1]

    resi_keys     = []
    cn            = []
    binary        = []
    keep          = []
    
    for resi_index_,resi_aa_ in zip(get_residues(structure)[0],get_residues(structure)[1]):
        key = (str(resi_aa_),chain,resi_index_)
        cn.append(None)
        binary.append(None)
        resi_keys.append(key)
    coords, seq = extract_coords_from_structure(structure)

    out = {"pdb_id":pdb_path,"seq":seq,"chain":chain,"coords":coords,"keys":resi_keys}
    return out

def _load_structure(fpath, chain=None):
    """
    Args:
        fpath: filepath to either pdb or cif file
        chain: the chain id or list of chain ids to load
    Returns:
        biotite.structure.AtomArray
    """
    if fpath.endswith('cif'):
        with open(fpath) as fin:
            pdbxf = pdbx.PDBxFile.read(fin)
        structure = pdbx.get_structure(pdbxf, model=1)
    elif fpath.endswith('pdb'):
        with open(fpath) as fin:
            pdbf = pdb.PDBFile.read(fin)
        structure = pdb.get_structure(pdbf, model=1)
    elif fpath.endswith(".pdb.gz"):
        with gzip.open(fpath,'rt', encoding='utf-8') as fin:
            pdbf = pdb.PDBFile.read(fin)
        structure = pdb.get_structure(pdbf, model=1)

    bbmask = filter_backbone(structure)
    structure = structure[bbmask]
    all_chains = get_chains(structure)
    if len(all_chains) == 0:
        raise ValueError('No chains found in the input file.')
    if chain is None:
        chain_ids = all_chains
    elif isinstance(chain, list):
        chain_ids = chain
    else:
        chain_ids = [chain]
    for chain in chain_ids:
        if chain not in all_chains:
            raise ValueError(f'Chain {chain} not found in input file')
    chain_filter = [a.chain_id in chain_ids for a in structure]
    structure = structure[chain_filter]
    return structure





def esmStructDataset(pdb_path,chain=None):#,chain):
    '''
    Convert PDB-file into dataset format

        Parameters:
            pdb_path (Path): path to pdb-file
            chain (str): antigen chain name
        Returns:
            dict (dict): dictionary, where keys are properties of the protein's tertiary structure
    '''

    return load_structure(str(pdb_path), chain)
    if chain is not None:
        chain = chain.upper()

    structure  = load_structure(str(pdb_path), chain)#chain.upper())
    resi_index = get_residues(structure)[0]
    resi_aa    = get_residues(structure)[1]
    resi_keys     = []
    cn = []
    binary = []
    
    for resi_index_,resi_aa_ in zip(get_residues(structure)[0],get_residues(structure)[1]):
        key = (str(resi_aa_),resi_index_)
        cn.append(None)
        binary.append(None)
        resi_keys.append(key)
    
    coords, seq = esm.inverse_folding.util.extract_coords_from_structure(structure)
    return {"pdb_id":pdb_path,"seq":seq,"chain":chain,"coords":coords,"cn":cn,"binary":binary,"residues":resi_keys}





class ESM1vForTokenClassification(nn.Module):
    def __init__(self, num_labels = 2):
        super().__init__()
        self.num_labels = num_labels    
        self.esm1v, self.esm1v_alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
        self.classifier = nn.Linear(512, self.num_labels)

    def forward(self, coords, padding_mask, confidence, tokens):
        prev_output_tokens = tokens[:, :-1]
        target = tokens[:, 1:]
        target_padding_mask = (target == alphabet.padding_idx)
        feat, x = self.esm1v.forward(coords, padding_mask, confidence, prev_output_tokens, features_only = True)
        f = feat[0,:,:]
        tt = torch.transpose(feat,1,2)
        logits = self.classifier(tt )
        return SequenceClassifierOutput(logits=logits)   


def loadModel(path = "./models/sema_3d_cn_atom_r1_8.0_r2_18.0_0.pth"):
    model = ESM1vForTokenClassification().cuda()
    model.load_state_dict(torch.load(path))
    return model

def get_prediction(sema_models, PDB_path, chain = None):
    struct = esmStructDataset(PDB_path, chain)
    custom_seq = struct["seq"]
    batch = [(struct["coords"], None, custom_seq)]
    batch_converter = CoordBatchConverter(alphabet)
    coords, confidence, strs, tokens, padding_mask = batch_converter(batch)
    with torch.no_grad():
        preds = []
        for sema_model in sema_models:
            pred = sema_model.forward(coords.cuda(),
                                     padding_mask.cuda(),
                                     confidence.cuda(),
                                     tokens.cuda())
            pred = pred.logits[:,:,1].squeeze().detach().cpu().numpy()
            preds.append(pred)
        pred = np.average(np.array(preds),axis=0)
    return pred, struct["keys"]

def screenAFDB():
    struct_path = Path("./afdb_domains/clusters/").glob("*")
    af_db_list  = []
    for p in struct_path:
        struct_path_ = str(p)+"/"+p.name+".pdb.gz" 
        af_db_list.append(struct_path_)
    return af_db_list

def main():
    af_db_list = screenAFDB()
    Path("afdb_sema_predictions").mkdir(exist_ok=True)
    model = loadModel()
    for path in af_db_list:
        out_name = path.split("/")[-1].replace(".pdb.gz",".pkl")
        out_path = f"./afdb_sema_predictions/{out_name}"
        sema_predictions, resi = get_prediction([model], 
                                             path,
                                             chain = "A")
        pickle.dump({"path":path, "score":sema_predictions, "resi":resi}, open(out_path,'wb'))
