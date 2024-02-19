import pymol
from pymol import cmd,stored
import networkx as nx
from pathlib import Path
import numpy as np
from scipy.spatial import distance
import itertools
import os
import subprocess

def load_protein(name, afdb_path = "./afdb"):
    """
    load protein in PyMOL
    """
    path = f"{afdb_path}/{name}.pdb.gz"
    cmd.do(f"load {path}, {name}")
 
    
def crop_domains(name = "AF-O89110-F1-model_v2", afdb_path = "./databases/alphafold_database/"):
    """
    Function to crop the full-length protein to domains with average high pLDDT value
    Cleavage is done on linkers with length of > 10 aa that are distant from protein regions with pLDDT > 60
    take an input pdb_database
    
    """
    
    cmd.do("delete *")
    load_protein(name, afdb_path)
    stored.xyz  = []
    stored.resi = []
    cmd.do(f"iterate_state 1, {name} & name CA & b>60, stored.xyz.append([x,y,z])")
    cmd.do(f"iterate_state 1, {name} & name CA & b>60, stored.resi.append((int(resi),resn))")

    if len(stored.xyz)<50:
        return []

    xyz = np.array(stored.xyz)
    cd  = distance.cdist(xyz,xyz)

    if len(stored.resi)==0:
        return
    
    fragments = [[stored.resi[0][0]]]
    for r in stored.resi[1:]:
        resi_prev = fragments[-1][-1]
        resi_cur  = r[0]
        if resi_cur == resi_prev:
            continue
        if resi_cur-resi_prev != 1:
            fragments.append([])
        fragments[-1].append(resi_cur)

        
    ids = {r[0]:ii for ii,r in enumerate(stored.resi)}
    Gfrag = nx.Graph()
    for f_id_1 in range(len(fragments)):
        Gfrag.add_edge(f_id_1, f_id_1)
        for f_id_2 in range(f_id_1+1,len(fragments)):
            f1_i = [ids[i] for i in fragments[f_id_1]]
            f2_i = [ids[i] for i in fragments[f_id_2]]
            con = np.where(cd[f1_i][:,f2_i]<8.0)
            n1  = len(set(con[0]))
            n2  = len(set(con[1]))
            if n1>=10:
                Gfrag.add_edge(f_id_1, f_id_2)
            if n2>=10:
                Gfrag.add_edge(f_id_1, f_id_2)

    domain_id = 0

    Path("./afdb_domains/").mkdir(exist_ok=True)
    Path("./afdb_domains/structs/").mkdir(exist_ok=True)
    Path("./afdb_domains/fasta/").mkdir(exist_ok=True)

    for g_ in [Gfrag.subgraph(c).copy() for c in nx.connected_components(Gfrag)]:
        domain_resi = set()
        for f_id in g_.nodes():
            domain_resi|=set(fragments[f_id])
        resi_start,resi_end = min(domain_resi),max(domain_resi)
        if resi_end-resi_start<50:
            continue
        cmd.do(f"create domain_{domain_id}_{name}, {name} & resi {resi_start}-{resi_end}")
        cmd.do(f"save ./afdb_domains/structs/{name}_{domain_id}.pdb.gz, domain_{domain_id}_{name}")
        cmd.do(f"save ./afdb_domains/fasta/{name}_{domain_id}.fasta, domain_{domain_id}_{name}")
        domain_id+=1

    return


def get_struct_list(afdb_path):
    all_path   = Path(afdb_path).glob("*.pdb.gz")
    calculated = Path("./afdb_domains/structs/").glob("*.pdb.gz")
    calculated = set([p.name.split("-")[1] for p in calculated])
    jobs = []
    for p in all_path:
        af_id = p.name.split("-")[1]
        if af_id in calculated:
            continue
        jobs.append(p)
    return jobs

def cluster_sequences():
    """
    Cluster domains with high pLDDT value according to sequence similarity using MMSeqs
    MMSeqs parameters:
    Cut-off 45%
    Coverage 70%
    Cluster-mode 1
    """

    if not os.path.exists("mmseqs_clusters"):
        os.mkdir("mmseqs_clusters")
    
    fastas = list(Path("./afdb_domains/fasta/").glob("*.fasta"))
    print(len(fastas))
    import random
    random.shuffle(fastas)
    fastas_batch = [fastas[i:i+100] for i in range(0,len(fastas),100)]
    
    fo = open("getseq.sh",'w')
    for u,f in enumerate(fastas_batch):
        com  = f"echo {u}\n"
        com += f"cat "+" ".join([str(f_) for f_ in f])+f" > mmseqs_clusters/afdb_{u}.fasta"
        fo.write(com+"\n")
    fo.close()

    fo = open("mmseqs_clusters/afdb.fasta",'w')
    for path in fastas:
        with path.open('r') as f:
            fo.write(f.read())
    fo.close()
    
    with open("mmseqs_clusters/run.sh",'w') as fo:        
        jobs   = f"mmseqs easy-cluster afdb.fasta results temp --min-seq-id 0.45 -c 0.7 --cov-mode 0 -s 8 --cluster-mode 1"
        fo.write(jobs)
    fo.close()

    subprocess.call("sh run.sh",shell=True, cwd="./mmseqs_clusters/")

def align_proteins_on_reference():
    """
    Align all proteins on the cluster center using PyMOL
    Keep only good alignments with RMSD<1.5 A
    """
    clusters = {}
    with open(f"./mmseqs_clusters/results_cluster.tsv") as f:
        for line in f:
            r = line.rstrip().split()
            clusters.setdefault(r[0],set())
            clusters[r[0]].add(r[1])
    print(clusters)

    def _parse_name(name):
        r = name.split("_")
        af_id,domain_id = "_".join(r[2:-1]),r[1]
        return af_id, domain_id

    def _get_path(name):
        #r = name.split("_")
        af_id,domain_id = _parseName(name)
        return f"./afdb_domains/structs/"+af_id+"_"+domain_id+".pdb.gz"

    for k,c in clusters.items():
        if len(c) == 1:
            continue
        cmd.do("delete *")
        ref_path = _getPath(k)
        af_id,domain_id = _parse_name(k)
        ref_name = af_id+"_"+domain_id
        cmd.do(f"load {ref_path}, {ref_name}")
        hits = [ref_name]
        for tar_ in c:
            af_id,domain_id = _parse_name(tar_)
            tar_name = af_id+"_"+domain_id
            tar_path = _get_path(tar_)
            if tar_path == ref_path:
                continue
            cmd.do(f"load {tar_path}, {tar_name}")
            r = cmd.align(f"{tar_name}",f"{ref_name}")
            rmsd = r[0]
            if rmsd>1.5:
                print("fail")
                continue
            hits.append(tar_name)

        if len(hits)<=1:
            continue
        
        if not os.path.exists("./afdb_domains/clusters/"):
            os.mkdir("./afdb_domains/clusters/")
        if not os.path.exists(f"./afdb_domains/clusters/{ref_name}"):
            os.mkdir(f"./afdb_domains/clusters/{ref_name}")
        for hit in hits:
            cmd.do(f"save ./afdb_domains/clusters/{ref_name}/{hit}.pdb.gz, {hit}")

    return clusters

def main():
    """
    Preprocess Alphafold Database
    """
    afdb_path = "./databases/alphafold_database/"
    for protein_path in get_struct_list(afdb_path):
        crop_domains(protein_path.name.split(".")[0], afdb_path)


    """
    Cluster sequences using MMSeqs
    """
    
    cluster_sequences()


    """
    Cluster sequences using MMSeqs
    """
    align_proteins_on_reference()

 
if __name__ == "__main__":
    main()
    
