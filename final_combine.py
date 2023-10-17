# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 23:17:02 2023

@author: yk5
"""

#final_combine

from rdkit import Chem
from rdkit.Chem import AllChem
import itertools


def flatten_tuple(tuple_3d):
    flattened = [item for sublist in tuple_3d for item in sublist]
    return flattened

def find_element_indices(lst, element):
    indices = [i for i, x in enumerate(lst) if x == element]
    return indices

def get_multi_idx(my_list_):
    my_list=[]
    for i in my_list_:
        if i !=False:
            if i not in my_list:
                my_list.append(i)
        else:
            continue
    multi=[]
    for i in range(2,len(my_list)+1):
        multi_=tuple(multi)
        list_=list(itertools.combinations(my_list, i))
        for j in list_:
            list_c=list(itertools.combinations(j, 2))
            ctrl=0
            for k in list_c:
                list_0=k[0]
                list_1=k[1]
                if len(set(tuple(list_0)) & set(tuple(list_1)))!=0:
                    ctrl=1
                    break
            if ctrl==1:
                continue
            else:
                multi.append(j)
        if tuple(multi_)==tuple(multi):
            break
    remove=[]
    for i in multi:
        for j in multi:
            if i!=j:
                if set(i) & set(j)==set(i):
                    remove.append(i)
    output = list(filter(lambda x: x not in remove, multi))
    for i in my_list:
        if i not in flatten_tuple(output):
            output.append((i,0))
    return output

def delete2dummy(smiles):
    if smiles.count('*')==0:
        return Chem.MolFromSmiles(smiles)
    two=False
    if smiles.count('*')==2:
        two=True
    mol = Chem.MolFromSmiles(smiles)
    idx=[]
    for atom in mol.GetAtoms():
        if atom.GetSymbol()=='*':
            idx=atom.GetIdx()
    edit_mol = Chem.EditableMol(mol)
    edit_mol.RemoveAtom(idx)
    mol_=edit_mol.GetMol()
    if two==True:
        for atom in mol_.GetAtoms():
            if atom.GetSymbol()=='*':
                idx=atom.GetIdx()
        edit_mol = Chem.EditableMol(mol_)
        edit_mol.RemoveAtom(idx)
        mol__=edit_mol.GetMol()
        
        return mol__
    else:
        return mol_

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

def recovery(mol_o):
    symbol_list=['Al','Si','As','Sc','Tl','Sb','Ti','Pb','Te','V','Po','Sb','In','Au','Ga','Ag','Ge','Pt','Bi']
    replace_list=['B','C','N','O','F','P','S','Cl','Se','Br','I','C','N','C','C','C','C','C','C']
    rw_o=Chem.RWMol(mol_o)
    for atom in mol_o.GetAtoms():
        if atom.GetSymbol() in symbol_list:
            idx_symbol=atom.GetSymbol()
            idx_replace=replace_list[symbol_list.index(idx_symbol)]
            rw_o.ReplaceAtom(atom.GetIdx(),Chem.Atom(idx_replace))
    mol=rw_o.GetMol() 
    return mol

def sanitization(smi,replace=False):
    try:
        if replace==True:
            mol=Chem.MolFromSmiles(smi.replace('2*','At'))
        else:
            mol=delete2dummy(smi)
        mol_=Chem.RemoveHs(recovery(mol))
        smi_=Chem.MolToSmiles(mol_,isomericSmiles=False).replace('[C]','C')
        mol__=Chem.MolFromSmiles(smi_)
        return mol__
    except:
        return False

def dummystruct_match(smi,smi_ori,smi_rep):
    try:
        mol=sanitization(smi)
        mol_ori_=sanitization(smi_ori)
        matches=mol.GetSubstructMatches(mol_ori_)[0]
    except:
        return False
    return matches
    
def dummystruct_rep(smi,smi_ori,smi_rep,matches):   
    mol=sanitization(smi)
    mol_ori=sanitization(smi_ori,True)
    break_indices=[]
    for match in matches:
        for nei in mol.GetAtomWithIdx(match).GetNeighbors():
            if nei.GetIdx() not in matches:
                bond_idx=mol.GetBondBetweenAtoms(match,nei.GetIdx()).GetIdx()
                break_indices.append([bond_idx,match,nei.GetIdx()])
    if smi_ori.count('*')==2:
        for idx in break_indices:
            for idx_ in break_indices:
                if idx!=idx_:
                    fragments = Chem.FragmentOnBonds(mol, [idx[0],idx_[0]],dummyLabels=[[1,1],[2,2]])
                    f = Chem.MolToSmiles(fragments).split('.')
                    for f_ in f:
                        if f_.count('*')==2:
                            try:
                                mol_f=Chem.MolFromSmiles(f_.replace('2*','At'))
                                if mol_f.GetSubstructMatches(mol_ori)!=():
                                    for f__ in f:
                                        if f__.count('*')==1:
                                            if '2*' in f__:
                                                f_2=f__.replace('2*','At')
                                            if '1*' in f__:
                                                f_1=f__
                                    m=Mol_Conn(smi_rep.replace('2*','At'),f_1,'*')
                                    p=Mol_Conn(m[0],f_2,'At')
                                    return Chem.MolToSmiles(recovery(Chem.MolFromSmiles(p[0])))
                            except:
                                continue

    else:
        for idx in break_indices:
            fragments = Chem.FragmentOnBonds(mol, [idx[0]],dummyLabels=[[3,3]])
            f = Chem.MolToSmiles(fragments).split('.')
            for f_ in f:
                if sanitization(f_)==False or sanitization(smi_ori)==False:
                    continue
                if sanitization(f_)!=sanitization(smi_ori):        
                    if '2*' in smi_rep:
                        try:
                            p=Mol_Conn(f_,smi_rep.replace('2*','At'),'*')
                            p_=p[0].replace('At','*')
                            p__=Chem.MolToSmiles(delete2dummy(p_))
                            return p__
                        except:
                            continue
                    else:
                        try:
                            p=Mol_Conn(f_,smi_rep,'*')
                            return Chem.MolToSmiles(recovery(Chem.MolFromSmiles(p[0])))
                        except:
                            continue
    return False
    
                            
def multi_replace(l,smi,target_oriented):
    i=0
    p=[smi]  
    
    while i!= len(l):    
        q=[]
        p1=[]
        #l_i一个fragment的一行
        for j in range(len(l[i])):
            if len(p1)<len(l[i]):
                p1=p1+p
            else:
                break
        for l_i in l[i]:
            line__=l_i.split(' ')
            line_=[x.strip() for x in line__ if len(x.strip()) > 0]            
            p_=p1[l[i].index(l_i)]
            smi_ori=line_[1]
            smi_rep=line_[2]
            matches=dummystruct_match(p_,smi_ori,smi_rep)
            if matches==False:
                continue
            else:
                rep=dummystruct_rep(p_,smi_ori,smi_rep,matches)
                if rep !=False:
                    q.append(rep)
        if q==[]:        
            p=p
        else:
            if target_oriented==False:
                p=q
            else:
                p=p+q
        i+=1
    return p


with open(r"BTK_covint_final_frags.txt",'r') as t:
    with open(r"BTK_covint_final.txt",'w+') as f:
        target_oriented=False
        lines=t.readlines()
        threshold=[0]
        for i in range(1,len(lines)):
            line__=lines[i].split(' ')
            line_=[x.strip() for x in line__ if len(x.strip()) > 0]
            line__0=lines[i-1].split(' ')
            line_0=[x.strip() for x in line__0 if len(x.strip()) > 0]
            if line_[0]!=line_0[0]:
                threshold.append(i)
                mol_list=lines[threshold[-2]:threshold[-1]]
                matches_=[]
                matches_idx=[]
                for line in mol_list:
                    line__=line.split(' ')
                    line_=[x.strip() for x in line__ if len(x.strip()) > 0]
                    smi=line_[0]
                    smi_ori=line_[1]
                    smi_rep=line_[2]
                    matches_.append(dummystruct_match(smi,smi_ori,smi_rep))
                    matches_idx.append(line_[3])
                combine_matches=get_multi_idx(matches_)
                for c in combine_matches:
                    l=[]
                    if 0 in c:
                        if matches_idx[matches_.index(c[0])]!='HB\n':
                            continue
                    for c_ in c:
                        if c_==0:
                            continue
                        start=find_element_indices(matches_, c_)[0]
                        end=find_element_indices(matches_, c_)[-1]                
                        l.append(mol_list[start:end+1])
                    print (l)
                    p=multi_replace(l,smi,target_oriented)
                    for p_ in p:
                        if p_!=False:
                            try:
                                recover_m=recovery(Chem.MolFromSmiles(p_))
                                recover_p=Chem.MolToSmiles(recover_m)
                            except:
                                continue
                            f.write(recover_p)
                            f.write('\n')
                
            