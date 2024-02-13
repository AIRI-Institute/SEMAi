import pymol
from pymol import cmd,stored
from scipy.spatial import distance
import numpy as np
import networkx as nx
from pathlib import Path
import os
import pickle
import gzip


def findSurfaceAtoms(selection="all", cutoff=2.5, quiet=1):
    """
DESCRIPTION

    Finds those atoms on the surface of a protein
    that have at least 'cutoff' exposed A**2 surface area.

USAGE

    findSurfaceAtoms [ selection, [ cutoff ]]

SEE ALSO

    findSurfaceResidues
    """
    cutoff, quiet = float(cutoff), int(quiet)

    tmpObj = cmd.get_unused_name("_tmp")
    cmd.create(tmpObj, "(" + selection + ") and polymer", zoom=0)

    cmd.set("dot_solvent", 1, tmpObj)
    cmd.get_area(selection=tmpObj, load_b=1)

    # threshold on what one considers an "exposed" atom (in A**2):
    cmd.remove(tmpObj + " and b < " + str(cutoff))

    selName = cmd.get_unused_name("exposed_atm_")
    cmd.select(selName, "(" + selection + ") in " + tmpObj)

    cmd.delete(tmpObj)

    if not quiet:
        print("Exposed atoms are selected in: " + selName)

    return selName


def findSurfaceResidues(selection="all", cutoff=2.5*3, doShow=0, quiet=1):
    """
DESCRIPTION

    Finds those residues on the surface of a protein
    that have at least 'cutoff' exposed A**2 surface area.

USAGE

    findSurfaceResidues [ selection, [ cutoff, [ doShow ]]]

ARGUMENTS

    selection = string: object or selection in which to find exposed
    residues {default: all}

    cutoff = float: cutoff of what is exposed or not {default: 2.5 Ang**2}

RETURNS

    (list: (chain, resv ) )
        A Python list of residue numbers corresponding
        to those residues w/more exposure than the cutoff.

    """
    cutoff, doShow, quiet = float(cutoff), int(doShow), int(quiet)

    selName = findSurfaceAtoms(selection, cutoff, quiet)

    exposed = set()
    cmd.iterate(selName, "exposed.add((chain,resv))", space=locals())

    selNameRes = cmd.get_unused_name("exposed_res_")
    cmd.select(selNameRes, "byres " + selName)

    if not quiet:
        print("Exposed residues are selected in: " + selNameRes)

    if doShow:
        cmd.show_as("spheres", "(" + selection + ") and polymer")
        cmd.color("white", selection)
        cmd.color("yellow", selNameRes)
        cmd.color("red", selName)

    return sorted(exposed)


def genProteinFragment(sel = [140,150], target = "3rdd", surface_residues = None):
    
    
    name = f"{target}_{sel[0]}_{sel[1]}"
    cmd.do(f"create p2p, {target} within 6.0 of ({target} & resi {sel[0]}-{sel[1]})")
    stored.r = []
    cmd.do(f"iterate p2p, stored.r.append((int(resi),resn))")
    
    r_out = [[]]
    for r1,r2 in zip(stored.r[:-1],stored.r[1:]):
        if r1[0] == r2[0]:
            continue
        if r2[0]-r1[0]==1:
            r_out[-1].append([r1,r2])
        else:
            r_out.append([])
        
    r_out = [r for r in r_out if len(r)>3]
    all_resi = set()
    if len(r_out)<2:
        cmd.do("delete p2p")
        return None
    for rr in r_out:
        for r in rr:
            all_resi.add(r[0][0])
            all_resi.add(r[1][0])
    sele_ = " or ".join([f"resi {r}" for r in all_resi])
    cmd.do("delete p2p")
    cmd.do(f"create {name}, ({sele_}) & {target}")
    return name

def keepFirstChain(target = "3rdd"):
    stored.r = []
    cmd.do(f"remove  {target} & not polymer.protein")
    cmd.do(f"iterate {target}, stored.r.append(chain)")
    cmd.do(f"remove  {target} & not chain {stored.r[0]}")

def getResidues(target):
    stored.r = []
    cmd.do(f"iterate {target} & name CA, stored.r.append(int(resi))")
    return stored.r

def splitOnSS(target):
    stored.r = []
    cmd.do(f"iterate {target} & name CA, stored.r.append([int(resi), ss])")
    ss_frags = [[stored.r[0]]]
    for r1, r2 in zip(stored.r[:-1],stored.r[1:]):
        if r1[1]!=r2[1]:
            ss_frags.append([])
        ss_frags[-1].append(r2)
    ss_residues = {}
    for i,s in enumerate(ss_frags):
        s_ = " or ".join([f"resi {q[0]}" for q in s])
        cmd.do(f"create s_{target}_{i}, {s_}")
        
        for resi in getResidues(f"s_{target}_{i}"):
            ss_residues[resi] = f"s_{target}_{i}"
    return ss_residues


def getDistaMatrix(target, test=False):
    stored.r = []
    cmd.do(f"iterate {target} & name CA, stored.r.append([int(resi), ss])")
    ss_frags = [[stored.r[0]]]
    for r1, r2 in zip(stored.r[:-1],stored.r[1:]):
        if r1[1]!=r2[1]:
            ss_frags.append([])
        ss_frags[-1].append(r2)

    G = nx.Graph()
    surface_residues = [r[1] for r in findSurfaceResidues(target)]
    for i,s in enumerate(ss_frags):
        for s_ in s:
            exposed = False
            if s_[0] in surface_residues:
                exposed = True
            G.add_node(s_[0], ss_type = s_[1], ss_id = i, exposed = exposed)
    
    stored.pos = []
    stored.resi = []
    cmd.do(f"iterate_state 1, {target}, stored.pos.append(  (x,y,z) )")
    cmd.do(f"iterate_state 1, {target}, stored.resi.append(int(resi))")
    stored.pos = np.array(stored.pos)
    cd         = distance.cdist(stored.pos,stored.pos)
    resi_ca = []
    resi_cd = []
    ids       = np.where(cd<5.0)
    for i1,i2 in zip(*ids):
        r1 = stored.resi[i1]
        r2 = stored.resi[i2]
        if r1 not in G.nodes():
            continue
        if r2 not in G.nodes():
            continue
        G.add_edge(stored.resi[i1],stored.resi[i2])

    site_id = 0

    site_residues = []
    for resi in G.nodes():

        if not G.nodes()[resi]['exposed']:
            continue

        target_ss   = {G.nodes()[resi]["ss_id"]: [resi]}
        
        for resi_ in list(G.neighbors(resi)):
            ss = G.nodes()[resi_]["ss_id"]
            if ss not in target_ss:
                target_ss[ss] = []
            target_ss[ss].append(resi_)
            
        all_sel = []

        site_ = set()

        if len(target_ss) <= 1:
            continue
        
        for ss in target_ss:
            mi,ma = min(target_ss[ss]),max(target_ss[ss])
            mi-=3
            ma+=3
            if mi<0:
                mi = 0
            sel   = f"resi {mi}-{ma}"
            all_sel.append(sel)
            site_|=set(range(mi,ma+1))
            
        all_sel = " or ".join(all_sel)
        site_id+=1
        
        #if test:

        cmd.do(f"create x{site_id}, {target} & ({all_sel})")
        cmd.do(f"color purple,   resi {resi} & x{site_id}")
        cmd.do(f"show  sticks,   resi {resi} & x{site_id}")
            
        stored.r = []
        cmd.do(f"iterate x{site_id} & name CA, stored.r.append(b)")
        
        b = np.average(stored.r)

        if b < 70:
            cmd.do(f"delete x{site_id}")
            continue

        site_residues.append(site_)        
        
    return site_residues

def processTarget(target_path = "7rwv",test = False):
    """
Script to extract protein regions,
to use it in the dataset of epitope-like regions
The protein is splitted to fragments based on its secondary structural elements
For solvent exposed surface residues  , the fragments within 5.0 A distance from the residue are collected to generate synthetic epitope-like structure.
Regions with pLDDT > 70 are collected
"""
    cmd.do("delete *")
    target = target_path.name.split(".")[0]
    cmd.do(f"load {target_path.resolve()}")
    keepFirstChain(target = target)
    ss_resi = splitOnSS(target = target)
    site_residues = getDistaMatrix(target, test = test)#True)
    return site_residues


def get_structures_list(path = "./afdb_domains/clusters/AF-C6DG98-F1-model_v2_0"):    
    domain_ids =          Path("./afdb_domains/clusters_alignments/").glob("*")
    
    all_path = []
    for domain_id_ in domain_ids:
        domain_id = domain_id_.name.replace(".pkl","")
        ref_path  = f"./afdb_domains/clusters/{domain_id}/{domain_id}.pdb.gz"
        assert os.path.exists(ref_path)
        all_path.append(Path(ref_path))

    return all_path



def epitopeScores(sites,scores):
    """
    extract epitopes based on average epitopes score or based on the size
    best_big_epitopes - epitopes with highest sum of SEMA score    
    best_ave_epitopes - epitopes with highest average SEMA score (currently used)
    """
    values = {i[2]:v for i,v in zip(scores["resi"],scores["score"])}
    
    site_scores_sum = []
    site_scores_ave = []

    for site in sites:
        site_scores_sum.append(np.sum([values[i] for i in site if i in values]))
        site_scores_ave.append(np.average([values[i] for i in site if i in values]))

    ### score by total SEMA epitope score
    sites_sum = sorted(enumerate(sites),key= lambda v:site_scores_sum[v[0]], reverse=True)

    ### score by average SEMA epitope score
    sites_ave = sorted(enumerate(sites),key= lambda v:site_scores_ave[v[0]], reverse=True)
   
    used_big_resi     = set()
    best_big_epitopes = []
    for i,site in sites_sum:
        score  = site_scores_sum[i]
        n_used = len(site & used_big_resi)
        p_used = n_used/len(site)

        ### if 25% of the epitope residues was already used - skip
        if p_used>0.25:
            continue
        
        used_big_resi|=set(site)
        best_big_epitopes.append(site)
        if len(best_big_epitopes) >= 2:
            break

    best_ave_epitopes = []
    for best_big_epitope in best_big_epitopes:
        resis    = list(best_big_epitope)
        scores   = [values[r] for r in resis if r in values]
        resi_max = resis[np.argmax(scores)]
        site_sele = None
        for i,site in sites_ave:
            if resi_max in site:
                site_sele = site
                break
        best_ave_epitopes.append(site_sele)
        
    return best_ave_epitopes


def af_pdb_alignemnt(af_id):
    path = f"/mnt/20tb/nvivanisenko/af_domains/afdb_domains/clusters_alignments/{af_id}"
    data = pickle.load(open(path,'rb'))
    return data

def get_epitope_antigen_map(antigen_residues, af_data):
    epitope_records = []
    for record in af_data:
        record_ids           = [r[0] in antigen_residues for r in record["resi_ref"]]
        cut_data             = {}
        cut_data["tar_name"] = record["tar_name"]
        cut_data["ref_name"] = record["ref_name"]
        
        for r in record:
            if r in {"tar_name","ref_name"}:
                continue
            cut_data[r] = [record[r][i] for i,m in enumerate(record_ids) if m]
            
        epitope_records.append(cut_data)
        
    return epitope_records

def save_sites_map(sites, pdb_alignment, out):
    epitope_map = {}
    for i,site in enumerate(sites):
        site_map = get_epitope_antigen_map(site, pdb_alignment)
        epitope_map[i] = site_map
    pickle.dump(epitope_map, open(out,'wb'))

def generate_epitope_like_structures():
    structs       = get_structures_list()

    for path in structs:
        domain_id      = path.name

        ### load alignment of structures within cluster
        pdb_alignment = af_pdb_alignment(domain_id.replace(".pdb.gz",".pkl"))

        ### extract protein fragments 
        sites         = processTarget(path, test = False)
        scores        = pickle.load(open("./afdb_sema_predictions/"+out.split("/")[-1],'rb'))

        ### crop antigen to domains and select fragment with highest SEMA average score
        small_epitope_sites = epitopeScores(sites,scores)

        ### save epitope-like structures mapped to homologoues full length protein
        save_sites_map(small_epitope_sites,
                       pdb_alignment,
                       out = "./afdb_sema_small_epitopes/" + domain_id.replace(".pdb.gz",".pkl"))
        

def main():
    generate_epitope_like_structures()
    
if __name__ == "__main__":
    main()

