# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 16:29:31 2023

@author: yk5
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import random

def Mol_Conn(s1,s2,attachment):
    conn = [] 
    m1 = Chem.MolFromSmiles(s1)
    m2 = Chem.MolFromSmiles(s2)    
    m = Chem.CombineMols(m1, m2)
    mw = Chem.RWMol(m)
    for i in range(0,m1.GetNumAtoms()):
        if m1.GetAtomWithIdx(i).GetSymbol() == attachment: 
            neighbor1_idx = m1.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx() 
            for j in range(0,m2.GetNumAtoms()):
                if m2.GetAtomWithIdx(j).GetSymbol() == attachment:
                    neighbor2_idx = m2.GetAtomWithIdx(j).GetNeighbors()[0].GetIdx()
                    mw.AddBond(neighbor1_idx,neighbor2_idx+m1.GetNumAtoms(),Chem.BondType.SINGLE)
                    mw.RemoveAtom(j+m1.GetNumAtoms())
                    mw.RemoveAtom(i)
                    conn.append(Chem.MolToSmiles(mw)) 
                    mw = Chem.RWMol(m) 
    return conn


with open(r"SMIRKS_ser.txt",'r') as h:
    lines_smirks=h.readlines()
    with open(r"delinker_ser_library.txt",'w+') as g:
        with open(r"covint_ser.txt",'r') as f:
            with open(r"delinker_ser_sub.txt",'r') as t:
                lines_t=t.readlines()
                mols=[]
                du = Chem.MolFromSmiles('*')
                for line in lines_t:
                    try:
                        line_=line.split(' ')
                        if line_[0].count('*')>1:
                            mols.append('1')
                            continue
                        mol_frag=Chem.MolFromSmiles(line_[0])
                        c1h=AllChem.DeleteSubstructs(mol_frag,du)
                        mols.append(c1h)
                    except:
                        mols.append('1')
                lines_f=f.readlines()#分子
                for line in lines_f:
                    print (lines_f.index(line))
                    try:
                        mol_=Chem.MolFromSmiles(line)
                    except:
                        continue
                    count=0
                    smiles=[]
                    for mol in mols:
                        print (mols.index(mol))
                        try:
                            line_=line.split(' ')                     
                            if mol_.GetSubstructMatch(mol)!=():
                                count+=1
                                result_mol = Chem.ReplaceSubstructs(mol_, mol, Chem.MolFromSmiles("[*:1]"), replaceAll=True)
                                my_list=lines_t[mols.index(mol)].split(' ')[1:]
                                random_selection = random.sample(my_list, 25)  
                                for i in random_selection:
                                    conn=Mol_Conn(Chem.MolToSmiles(result_mol[0]),i,'*')
                                    if '.' in conn[0]:
                                        break
                                    if '*' in conn[0]:
                                        break
                                    else:
                                        smiles.append(conn[0])
                        except:
                            continue
                    if count==0:
                        for smirks in lines_smirks:        
                            original_mol = mol_
                            reaction = AllChem.ReactionFromSmarts(smirks)
                            transformed_mols = reaction.RunReactants((original_mol,))
                            for transformed_mol in transformed_mols:
                                transformed_smiles = Chem.MolToSmiles(transformed_mol[0])
                                if '*' not in transformed_smiles:
                                    smiles.append(transformed_smiles)
                    if len(smiles)>25:
                        smiles_=random.sample(smiles, 25)  
                    else:
                        smiles_=smiles
                    for smi in smiles_:
                        g.write(smi)
                        g.write('\n')
                        
