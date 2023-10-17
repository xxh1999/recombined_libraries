# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 14:00:06 2023

@author: yk5
"""

import sys
import os
sys.path.append("../")
sys.path.append("../analysis/")
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
from rdkit.Chem import MolStandardize

import numpy as np

from itertools import product
from joblib import Parallel, delayed
import re
from collections import defaultdict

from IPython.display import clear_output
IPythonConsole.ipython_useSVG = True

from DeLinker_test import DenseGGNNChemModel
import frag_utils
import rdkit_conf_parallel
from data.prepare_data import read_file, preprocess
import example_utils

def get_atom_pair_idx(molA,typ):
    total=[]  
    if typ=='ser':
        amo='OCC(N)C(=O)'
    if typ=='cys':
        amo='SCC(N)C(=O)'        
    mol_amo=Chem.MolFromSmiles(amo)
    match_indices = molA.GetSubstructMatches(mol_amo)
    if match_indices==():
        return None
    print (match_indices)
    for match in match_indices[0]:
        for atom in molA.GetAtomWithIdx(match).GetNeighbors():
            print (atom.GetIdx())
            if atom.GetIdx() not in match_indices[0]:
                idx_S=match
                idx_S_=atom.GetIdx()
                break
    dis_s=Chem.rdmolops.GetDistanceMatrix(molA)[idx_S]
    judge=5
    idx_break=[]
    idx_break_=[]
    while True:
        for i in range(len(dis_s)):
            if int(dis_s[i])==judge:
                if int(i) not in match_indices:
                    if molA.GetAtomWithIdx(int(i)).IsInRing()==False:
                        idx_break.append(int(i))
        if idx_break!=[]:
            frags=[]
            for i in idx_break:
                neighbours = molA.GetAtomWithIdx(i).GetNeighbors()
                list=[]
                if len(neighbours)==1:
                    continue
                for neighbour in neighbours:
                    bond = molA.GetBondBetweenAtoms(neighbour.GetIdx(),i)
                    mol1_f = Chem.FragmentOnBonds(molA, [bond.GetIdx()],addDummies=False)
                    smi_A = Chem.MolToSmiles(mol1_f)
                    print (smi_A)
                    if '.' in smi_A:
                        smi_list=smi_A.split('.')
                        for smi in smi_list:
                            if Chem.MolFromSmiles(smi).GetSubstructMatches(mol_amo)!=():
                                mol_smi=Chem.MolFromSmiles(smi)
                                print (mol_smi.GetNumAtoms())
                                if mol_smi.GetNumAtoms()<20:
                                    list.append(smi)
                list1 = sorted(list, key=lambda i_: len(i_))
                if len(list1)>0:
                    frags.append(list1[0])
                    idx_break_.append(i)        
            if frags!=[]:
                frags_ = sorted(frags, key=lambda i_: len(i_))   
                print (frags_)
                frag=frags_[0]                
                break
        judge-=1
        if judge==1:
            break
    idxA=idx_break_[frags.index(frag)]
    neighbours = molA.GetAtomWithIdx(idxA).GetNeighbors()
    list=[]
    for neighbour in neighbours:
        bond = molA.GetBondBetweenAtoms(neighbour.GetIdx(),idxA)
        mol1_f = Chem.FragmentOnBonds(molA, [bond.GetIdx()])
        smi_B = Chem.MolToSmiles(mol1_f)
        if '.' in smi_B:
            smi_list=smi_B.split('.')
            for smi in smi_list:
                if Chem.MolFromSmiles(smi).GetSubstructMatches(mol_amo)!=():
                    list.append([smi,neighbour.GetIdx()])
    list1 = sorted(list, key=lambda i: len(i[0]))
    idxB=list1[1][1]
    pairs=[[idx_S,idx_S_],[idxA,idxB]]
    return pairs
        


# How many cores for multiprocessing
n_cores = 4
# Whether to use GPU for generating molecules with DeLinker
use_gpu = False

paths=os.listdir(r"D:\xxh\warhead_substitution\res")
for i in paths:
    scaff_1_path = os.path.join(r"D:\xxh\warhead_substitution\res",i)
    print (i)
    try:
        scaff_1_sdf = Chem.SDMolSupplier(scaff_1_path)
        scaff_1_smi = Chem.MolToSmiles(scaff_1_sdf[0])
    except:
        continue
    
    img = Draw.MolsToGridImage([Chem.MolFromSmiles(scaff_1_smi)], molsPerRow=2, subImgSize=(300, 300))
    
    starting_point_2d = Chem.Mol(scaff_1_sdf[0])
    # _ = AllChem.Compute2DCoords(starting_point_2d)
    # example_utils.mol_with_atom_index(starting_point_2d)
    
    pairs=get_atom_pair_idx(starting_point_2d,'cys')
    if pairs!=None:
        atom_pair_idx_1 = pairs[0]
        atom_pair_idx_2 = pairs[1]
    bonds_to_break = [starting_point_2d.GetBondBetweenAtoms(x,y).GetIdx() for x,y in [atom_pair_idx_1, atom_pair_idx_2]]
    
    fragmented_mol = Chem.FragmentOnBonds(starting_point_2d, bonds_to_break)
    _ = AllChem.Compute2DCoords(fragmented_mol)
    
    # Split fragmentation into core and fragments
    fragmentation = Chem.MolToSmiles(fragmented_mol).split('.')
    fragments = []
    for fragment in fragmentation:
        if len([x for x in fragment if x =="*"]) ==2:
            linker=fragment
        else:
            fragments.append(fragment)
    fragments = '.'.join(fragments)
    linker = re.sub('[0-9]+\*', '*', linker)
    fragments = re.sub('[0-9]+\*', '*', fragments)
    
    # Get distance and angle between fragments
    dist, ang = frag_utils.compute_distance_and_angle(scaff_1_sdf[0], linker, fragments)
    
    # Write data to file
    data_path = "./scaffold_hopping_test_data.txt"
    with open(data_path, 'w') as f:
        f.write("%s %s %s" % (fragments, dist, ang))
    raw_data = read_file(data_path)
    preprocess(raw_data, "zinc", "scaffold_hopping_test", True)
    
    import os
    if not use_gpu:
        os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
    # Arguments for DeLinker
    args = defaultdict(None)
    args['--dataset'] = 'zinc'
    args['--config'] = '{"generation": true, \
                          "batch_size": 1, \
                          "number_of_generation_per_valid": 100, \
                          "min_atoms": 9, "max_atoms": 10, \
                          "train_file": "molecules_scaffold_hopping_test.json", \
                          "valid_file": "molecules_scaffold_hopping_test.json", \
                          "output_name": "DeLinker_example_generation_scaffold_hopping.smi"}'
    args['--freeze-graph-model'] = False
    args['--restore'] = '../models/pretrained_DeLinker_model.pickle'
    # Setup model and generate molecules
    model = DenseGGNNChemModel(args)
    model.train()
    # Free up some memory
    model = ''