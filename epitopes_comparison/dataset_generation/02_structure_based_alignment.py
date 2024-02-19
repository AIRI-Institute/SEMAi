import pymol
from   pymol import cmd, stored
from   pathlib import Path
import itertools
import numpy as np
import pickle
import os
from   scipy.spatial import distance
one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', 'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C', "-":"-"}


"""
Get structure-based alignment based on structures aligned in PyMOL
"""
def get_structures_list(path = "afdb_domains/clusters/AF-C6DG98-F1-model_v2_0"):
    paths = Path(path).glob("*.pdb.gz")
    return list(paths)


def ref_vs_all(domain_id = "AF-C6DG98-F1-model_v2_0"):
    path = f"./afdb_domains/clusters/{domain_id}"
    struct_list = get_structures_list(path)
    ref_path      = Path(path+"/"+domain_id+".pdb.gz")
    alignments = []
    out = f"./afdb_domains/clusters_alignments/{domain_id}.pkl"
    for path in struct_list:
        if ref_path.name == path.name:
            continue
        alignment = one_vs_one(path,ref_path)
        alignments.append(alignment)
    if not os.path.exists("./afdb_domains/clusters_alignments/"):
        os.mkdir("./afdb_domains/clusters_alignments/")
    pickle.dump(alignments, open(out,'wb'))
    

def smooth_errors(quality):
    for i in range(5,len(quality)-5):
        c1 = quality[-5+i:i]
        c2 = quality[i:i+5]
        if np.sum(c1)==5:
            continue
        if np.sum(c2)==5:
            continue
        quality[i] = 0
    return quality

def misalignment_detector(alignment):
    """
    Two residues could be within close proximity distance by accident
    This function check if peptides within +- 5 aa are also aligned
    """
    seq_a   = []
    seq_b   = []
    quality = []
    for a in alignment:
        add_a = a["a_id"]
        add_b = a["b_id"]
        n     = len(a["a_id"])-len(a["b_id"])
        rmsd  = a["rmsd"]
        if n < 0:
            add_a   += [[None,"-"] for i in range(-n)]
        if n > 0:
            add_b   += [[None,"-"] for i in range(n)]
        assert len(add_a) == len(add_b)        
        quality     += [rmsd for _ in add_a]
        seq_a += add_a
        seq_b += add_b
    
    seq_a_aa = [one_letter[s[1]] for s in seq_a]
    seq_b_aa = [one_letter[s[1]] for s in seq_b]
    q        = [0 if (q_ is None or q_>=1.5) else 1 for q_ in quality]
    q_before = [_ for _ in q]
    q        = smooth_errors(q)
    
    assert len(seq_b_aa) == len(seq_a_aa)
    assert len(quality) == len(seq_a_aa)

    out = {"sequence_tar"     : seq_a_aa,
           "sequence_ref"     : seq_b_aa,
           "correctly_aligned": q,
           "resi_tar"         : seq_a,
           "resi_ref"         : seq_b,
           "rmsd"             : quality }

    return out
    
    
def one_vs_one(path1, path2):
    """
    Prepare structure-based alignment of two proteins
    Two residue considered to have same alignment positions if distance between CA atoms is lower 2.0 A
    If distance > 2.0, residues are considered to be unaligned    
    """
    
    cmd.do("delete *")

    n1 = f"n1_{path1.name}"
    n2 = f"n2_{path2.name}"
    
    cmd.do(f"load {path1}, {n1}")
    cmd.do(f"load {path2}, {n2}")
    
    stored.xyz1   = []
    stored.xyz2   = []

    stored.resi1  = []
    stored.resi2  = []
    
    cmd.do(f"iterate_state 1, {n1} & name CA, stored.xyz1.append([x,y,z])")
    cmd.do(f"iterate_state 1, {n2} & name CA, stored.xyz2.append([x,y,z])")
            
    stored.xyz1 = np.array(stored.xyz1)
    stored.xyz2 = np.array(stored.xyz2)

    cd = distance.cdist(stored.xyz1,stored.xyz2)
    
    cmd.do(f"iterate_state 1, {n1} & name CA, stored.resi1.append([int(resi),resn])")
    cmd.do(f"iterate_state 1, {n2} & name CA, stored.resi2.append([int(resi),resn])")

    a_list = []
    b_list = []
    for i in range(len(stored.resi1)):
        ids = np.where(cd[i]<2.0)[0]
        a_list.append(stored.resi1[i])

        if len(ids)!=1:
            b_list.append([])
            continue
        j    = ids[0]
        rmsd = cd[i,j]
        b_list.append([stored.resi2[j],rmsd])

    b_resi = {r[0]:r for r in stored.resi2}
    b_range = []

    gap = False

    last_id = min(b_resi)
    alignment = [{"a_id":[],"b_id":[],"rmsd":None}]

    for a_,b_ in zip(a_list,b_list):
        if len(b_)!=0:
            alignment.append({"a_id":[a_],"b_id":[b_[0]],"rmsd":b_[1]})
            continue
        if alignment[-1]["rmsd"] is not None:
            alignment.append({"a_id":[],"b_id":[],"rmsd":None})
        alignment[-1]["a_id"].append(a_)

    for i in range(1,len(alignment)-1):
        a = alignment[i]
        if a["rmsd"] is not None:
            continue
        resi_start = alignment[i-1]["b_id"][0][0]
        resi_end   = alignment[i+1]["b_id"][0][0]
        delta      = resi_end-resi_start
        if delta<=0:
            continue
        alignment[i]["b_id"] = [b_resi[r] for r in range(resi_start+1,resi_end)]

    cmd.do(f"alter {n1}, b=-1.0")
    cmd.do(f"alter {n2}, b=-1.0")
    
    alignment_adj = misalignment_detector(alignment)

    alignment_adj["ref_name"] = path2
    alignment_adj["tar_name"] = path1
    
    return alignment_adj
        
        
def prepare_all_alignments():
    domain_ids = Path("./afdb_domains/clusters/").glob("*")
    for domain_id_ in domain_ids:
        domain_id = domain_id_.name
        ref_vs_all(domain_id)
    
if __name__ == "__main__":
    prepare_all_alignments()
    
