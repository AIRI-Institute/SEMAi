import pymol
from pymol import cmd,stored
from pymol.exporting import _resn_to_aa as one_letter
import sys    
import numpy as np
import pickle
import os
import pandas as pd
cmd.do("bg_color white")
import ast
from matplotlib import cm

def makePalette():
    rgb = cm.get_cmap("coolwarm")
    for i in range(21):
        v = float(i)/20.0
        r,g,b,x = rgb(v)
        cmd.do(f"set_color color_{i},[{r},{g},{b}]")


        
def atomSele(r_,obj = "all"):
    resi,resn,chain,name = r_
    return f"({obj} & resi {resi} & resn {resn} & chain {chain} & name {name})"

def loadStruct(path = "./fake_fab_antigen/6ki8_B_A.ent"):
    pdb_id,chain_antigen,chain_fab = path[:-4].split("/")[-1].split("_")
    cmd.do(f"load {path}, antigen_fab")
    cmd.do("remove heta")


    chain_fab_sele = [f"(chain {c})" for c in chain_fab]
    chain_fab_sele = " or ".join(chain_fab_sele)
    
    cmd.do(f"create antigen, antigen_fab & chain {chain_antigen}")
    cmd.do(f"create fab,     antigen_fab & {chain_fab_sele}")
    cmd.do("color white, antigen")
    cmd.do("color gray,  fab")


def getProximityAtoms(atom_sele, target = "fab", R = 6.0):
    stored.r = set()
    cmd.do(f"iterate {target} within {R} of {atom_sele}, stored.r.add((resi,resn,chain,name))")
    return stored.r

    
def findContactingResidues(R = 6.0):
    antigen_sele = f"antigen within {R} of fab"
    stored.contacting_residues = set()
    cmd.do(f"iterate {antigen_sele}, stored.contacting_residues.add((resi,resn,chain,name))")
    cmd.do(f"color red, {antigen_sele}")
    contacts = {}
    for atom in stored.contacting_residues:
        atom_sele = atomSele(atom, obj = "antigen")
        fab_atoms = getProximityAtoms(atom_sele, R=R)
        contacts[atom] = fab_atoms
    return contacts

def atomToResi(r):
    return {(r_[0],r_[1],r_[2]) for r_ in r}

def resiSele(r_,obj = "all"):
    resi,resn,chain = r_
    return f"({obj} & resi {resi} & resn {resn} & chain {chain})"
    
def atomSele(r_,obj = "all"):
    resi,resn,chain,name = r_
    return f"({obj} & resi {resi} & resn {resn} & chain {chain} & name {name})"

def contactNumberColors(cn, cn_type = "CA"):
        
    if cn_type == "cn_raw":
        norm = 200.0
    if cn_type == "cn_non_raw":
        norm = 50.0
    if cn_type == "CA":
        norm = 10.0

    for r in cn:
        resi = r["resi_tuple"]
        v = r[cn_type]
        if v == -100:
            color = "white"
        else:
            vc   = v/norm
            if vc > 1.0:
                vc = 1.0
            vc = int(vc*20)
            color = f"color_{vc}"    
        cmd.do(f"color {color}, {resiSele(resi)}")

        
def contactNumbers(contacts):
    cn_CA      = {}
    cn_non_raw = {}
    cn_raw     = {}
    for antigen_atom,fab_atoms in contacts.items():

        resi,resn,chain,atom_name = antigen_atom
        antigen_resi_key          = (resi,resn,chain)

        cn_CA.setdefault(      antigen_resi_key, set() )
        cn_non_raw.setdefault( antigen_resi_key, set() )
        cn_raw.setdefault(     antigen_resi_key,    [] )

        for fab_atom in fab_atoms:
            cn_non_raw[antigen_resi_key].add(fab_atom)

        for fab_atom in fab_atoms:
            resi_,resn_,chain_,atom_name_ = fab_atom
            cn_CA[antigen_resi_key].add((resi_,resn_,chain_))

        for fab_atom in fab_atoms:
            cn_raw[antigen_resi_key].append(fab_atom)

    cn_CA      = {r:len(v) for r,v in      cn_CA.items()}
    cn_non_raw = {r:len(v) for r,v in cn_non_raw.items()}
    cn_raw     = {r:len(v) for r,v in     cn_raw.items()}

    return {"cn_CA":cn_CA,"cn_non_raw":cn_non_raw,"cn_raw":cn_raw}


def getAntigenData(R0 = 6.0, R1 = 12.0):
    stored.r = []
    cmd.do(f"iterate antigen & name CA, stored.r.append((resi,resn,chain))")
    res  = []
    for resi,resn,chain in stored.r:
        res.append({"aa":one_letter[resn],"resi_tuple":(resi,resn,chain),"cn_CA":-100,"cn_raw":-100,"cn_non_raw":-100,"R0":R0,"R1":R1})

    contacts          = findContactingResidues(R = R0)
    distant_contacts  = findContactingResidues(R = R1)

    distant_contacts  = {d for d in atomToResi(distant_contacts) if d not in atomToResi(contacts)}

    contact_numbers   = contactNumbers(contacts)

    res_out = []
    for i,r in enumerate(res):

        resi_key = r["resi_tuple"]
        if resi_key in contact_numbers["cn_CA"]:
            assert resi_key in contact_numbers["cn_raw"]
            assert resi_key in contact_numbers["cn_non_raw"]    
            r["cn_CA"]      = contact_numbers["cn_CA"][resi_key]
            r["cn_raw"]     = contact_numbers["cn_raw"][resi_key]
            r["cn_non_raw"] = contact_numbers["cn_non_raw"][resi_key]

        if resi_key in distant_contacts:
            assert r["cn_CA"] == -100
            r["cn_CA"]           = 0
            r["cn_raw"]          = 0
            r["cn_non_raw"]      = 0
        
    return res
    
def processFab():
    makePalette()

    pdb_name = sys.argv[-3]
    R0       = sys.argv[-2]
    R1       = sys.argv[-1]
    
    pdb_path = f"./structs_antigen_fab/{pdb_name}.pdb"
    assert os.path.exists(pdb_path)
    loadStruct(pdb_path)
    antigen_data = getAntigenData(R0,R1)

    out_path  = f"./structs_antigen_fab_df"

    if not os.path.exists(out_path):
        os.mkdir(out_path)
        
    out_path += f"/{pdb_name}_{R0}_{R1}.pkl"

    pickle.dump(antigen_data, open(out_path,'wb'))
    
    out_path  = f"./structs_antigen_fab_png"
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    out_path += f"/{pdb_name}.png"
    
    contactNumberColors(antigen_data, "cn_raw")

    cmd.do(f"zoom antigen within {R1} of fab")
    cmd.do("disable antigen_fab")
    cmd.do("show lines")
    cmd.do(f"png {out_path}")
    
    
if __name__ == "__main__":
    processFab()
    
