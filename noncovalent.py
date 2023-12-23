# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 10:48:23 2023

@author: yk5
"""

import os
from rdkit import Chem
from rdkit.Chem import AllChem
import itertools
import argparse

def replace_atom(smi,symbol=False):
    if symbol:
        mol_=Chem.MolFromSmiles(smi)
        for atom in mol_.GetAtoms():
            if atom.GetSymbol()=='Si':
                idx=atom.GetIdx()
                break
        em = Chem.RWMol(mol_)
        em.ReplaceAtom(idx,Chem.Atom(symbol))
        m=em.GetMol() 
        return m
    else:
        return Chem.MolFromSmiles(smi)

def delete2dummy(smiles,two=False):
    mol = Chem.MolFromSmiles(smiles)
    if mol==None:
        print (mol)
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

def dis_cal(l1,l2):
    dis=(float(l1[0])-float(l2[0]))**2+(float(l1[1])-float(l2[1]))**2+(float(l1[2])-float(l2[2]))**2
    return dis

def sidechain_removal(smi,interaction):
    if interaction=='HI':
        replace_list=['Au','Ga']
    elif interaction=='HB':
        replace_list=['Al','Si','As','Sc','Tl','Sb','Ti','Pb','Te','V','Po','Sb','In']
    elif interaction=='PS':
        replace_list=['Ag','Ge']
    elif interaction=='PC':
        replace_list=['Pt','Bi']
    bond_list=[]
    mol = Chem.MolFromSmiles(smi)
    for atom in mol.GetAtoms():
        for n in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(),n.GetIdx())
            if bond.IsInRing()==False:
                fragments=Chem.FragmentOnBonds(mol,[bond.GetIdx()],addDummies=False)
                f = Chem.MolToSmiles(fragments).split('.')
                for i in f:
                    x=0
                    for sym in replace_list+['*']:
                        if sym in i:
                            x+=1
                    if x==0:
                        bond_list.append(bond.GetIdx())
    
    if bond_list==[]:
        return smi
    else:
        fragments_=Chem.FragmentOnBonds(mol,list(set(bond_list)),addDummies=False)
        f_ = Chem.MolToSmiles(fragments_).split('.')
        for i in f_:
            if '*' in i:
                return i

def get_idx_from_pdb(Atom_coor,lig_pdb_path):
    implicit=[]
    atom_idx=[]
    with open(lig_pdb_path,'r') as t:
        lines=t.readlines()
        for line in lines:
            if line[:4]=='ATOM' or line[:6]=='HETATM':
                line__=line.split(' ')
                line_=[x.strip() for x in line__ if len(x.strip()) > 0]
                if line_[2][0]!='H':
                    implicit.append(line)
        for i in Atom_coor:
            for j in implicit:
                if j[30:38].strip(' ')+','+j[38:46].strip(' ')+','+j[46:54].strip(' ')==i[0]:
                    atom_idx.append(implicit.index(j))   
        return atom_idx
                

def get_coordinate_with_idx(lig_pdb_path,idx):
    lines_=[]
    with open(lig_pdb_path,'r') as t:
        lines=t.readlines()
        for line in lines:
            if line[:4]=='ATOM' or line[:6]=='HETATM':
                line__=line.split(' ')
                line_=[x.strip() for x in line__ if len(x.strip()) > 0]
                if line_[2][0]!='H':
                    lines_.append(line)                    
        l=lines_[int(idx)]
        coor = [l[30:38].strip(' '),l[38:46].strip(' '),l[46:54].strip(' ')]
        return coor
            
def get_ring_systems(lists):
    ls=[]
    while True:
        dup=[]        
        for m in lists:
            for n in lists:
                if m!=n:
                    if len(set(m)&set(n))>0:              
                        ls=tuple(set(m) | set(n))
                        dup=[m,n]
        if dup==[]:           
            list_=list(set(lists))
            return list_
            break
        else:
            for i in dup:
                lists.remove(i)
            lists.append(ls)

def HA_detector(ori_pdb_path,indices):
    with open (ori_pdb_path,'r') as t:
        lines=t.readlines()
        for line in lines:
            if line[:6].strip()=='ATOM':
                if line[6:11].strip() in indices:
                    if line[13:17].strip()=='O' or 'SD' or 'OE1' or 'OE2' or 'OD1' or 'OD2':#只能做氢键受体的氨基酸原子：a
                        indices[indices.index(line[6:11].strip())]='a'
        for i in range(len(indices)):
            if indices[i]!='a':
                indices[i]='c'
        return indices
    
def HA_detector_lig(D_A,idx,lig_pdb_path):
    mol=AllChem.MolFromPDBFile(lig_pdb_path)
    for i in range(len(idx)):
        atom=mol.GetAtomWithIdx(idx[i])
        hasH = atom.GetTotalNumHs()
        if hasH==0:#只能做氢键受体的小分子上原子：b a和b不能配对
            D_A[i]='b'
    return D_A
                


#get coordinates of 2 interaction atoms
def get_interaction_atom_coor(ori_pdb_path,report_path,lig_pdb_path,interaction_type):
    coor_origin=[]    
    with open (lig_pdb_path,'r') as f:
        lines_lig=f.readlines()
        for line in lines_lig:
            line__=line.split(' ')
            line_=[x.strip() for x in line__ if len(x.strip()) > 0]
            if len(line_)>10:
                coor_origin.append(line[30:38].strip(' ')+','+line[38:46].strip(' ')+','+line[46:54].strip(' '))
        Hydrophobic_Interactions=[]
        Hydrogen_Bonds=[]
        Water_Bridges=[]
        Salt_Bridges=[]
        pi_Stacking=[]
        pi_Cation_Interactions=[]
        Halogen_Bonds=[]
        count=0
        atom_idx=[]
        atom_coor=[]
        HA=[]
        with open (report_path,'r') as t:
            lines_=t.readlines()
            idx_UNL=[]             
            for line in lines_:
                if len(line)>3:
                    if ':' == line[3]:
                        idx_UNL.append(lines_.index(line))
            for idx_ in idx_UNL:
                lines=lines_[idx_:]
                for line in lines:
                    if line =='**Hydrophobic Interactions**\n':
                        if interaction_type=='HI':
                            count=count+1                
                            HB = lines.index(line)
                            lines_HB = lines[HB+1:]
                            for line in lines_HB:
                                if line[0]=='*':
                                    break
                                if line =='\n':
                                    continue
                                line__=line.split(' ')
                                line_=[x.strip() for x in line__ if len(x.strip()) > 0]
                                if line_[0]=='|':
                                    if len(line_[9])<4:   
                                        if line_[19]+line_[20]+line_[21] in coor_origin:
                                            atom_coor.append([line_[19]+line_[20]+line_[21]])
                            
                    elif line =='**Hydrogen Bonds**\n':
                        if interaction_type=='HB':
                            count=count+1                
                            HB = lines.index(line)
                            lines_HB = lines[HB+1:]
                            for line in lines_HB:
                                if line[0]=='*':
                                    break
                                if line =='\n':
                                    continue
                                line__=line.split(' ')
                                line_=[x.strip() for x in line__ if len(x.strip()) > 0]
                                if line_[0]=='|':
                                    if len(line_[9])<4:   
                                        if line_[31]+line_[32]+line_[33] in coor_origin:
                                            #The three amino acid groups that can act as hydrogen bond acceptors are the carbonyl group of the peptide bond, the carboxyl group of aspartic acid and glutamic acid, and the sulfur atom of methionine.
                                            if int(line_[23])>int(line_[27]):
                                                HA.append(line_[27])
                                            else:
                                                HA.append(0)                                        
                                            atom_coor.append([line_[31]+line_[32]+line_[33]])                                      
                    elif line =='**Salt Bridges**\n':
                        if interaction_type=='SB':
                            count=count+1                
                            HB = lines.index(line)
                            lines_HB = lines[HB+1:]
                            for line in lines_HB:
                                if line[0]=='*':
                                    break
                                if line =='\n':
                                    continue
                                line__=line.split(' ')
                                line_=[x.strip() for x in line__ if len(x.strip()) > 0]
                                if line_[0]=='|':
                                    if len(line_[11])<4:   
                                        if line_[23]+line_[24]+line_[25] in coor_origin:
                                            atom_coor.append([line_[23]+line_[24]+line_[25]])
                    elif line =='**pi-Stacking**\n':
                        if interaction_type=='PS':                           
                            count=count+1                
                            HB = lines.index(line)
                            lines_HB = lines[HB+1:]
                            for line in lines_HB:
                                if line[0]=='*':
                                    break
                                if line =='\n':
                                    continue
                                line__=line.split(' ')
                                line_=[x.strip() for x in line__ if len(x.strip()) > 0]
                                if line_[0]=='|':
                                    if len(line_[9])<4:   
                                        atom_list=line_[23].split(',')
                            with open (ori_pdb_path,'r') as o:
                                lines_o=o.readlines()
                                for line in lines_o:
                                    if line[:6] in ['HETATM','ATOM']:
                                        if line[77]=='C':
                                            if line[6:11].strip(' ') in atom_list:
                                                atom_coor.append([line[30:38].strip(' ')+','+line[38:46].strip(' ')+','+line[46:54].strip(' ')] )
                    elif line =='**pi-Cation Interactions**\n':
                        if interaction_type=='PC':                           
                            count=count+1                
                            HB = lines.index(line)
                            lines_HB = lines[HB+1:]
                            for line in lines_HB:
                                if line[0]=='*':
                                    break
                                if line =='\n':
                                    continue
                                line__=line.split(' ')
                                line_=[x.strip() for x in line__ if len(x.strip()) > 0]
                                if line_[0]=='|':
                                    if len(line_[11])<4:   
                                        atom_list=line_[23].split(',')
                            with open (ori_pdb_path,'r') as o:
                                lines_o=o.readlines()
                                for line in lines_o:
                                    if line[:6] in ['HETATM','ATOM']:
                                        if line[77]=='C':
                                            if line[6:11].strip(' ') in atom_list:
                                                atom_coor.append([line[30:38].strip(' ')+','+line[38:46].strip(' ')+','+line[46:54].strip(' ')] )
                if atom_idx!=[]:
                    break
            if interaction_type=='HB':
                D_A=HA_detector(ori_pdb_path,HA)
            else:
                D_A=[]
                for i in range(len(atom_coor)):
                    D_A.append('c')
            if count==0:
                return False
            else:
                return [atom_coor,atom_idx,D_A]

#obtained a list of fragments and the coordinates of two bond break positions
def get_pairs(lig_pdb_path,idx_,Atom_coor,mol,interaction,D_A_):
    mol = Chem.MolFromMolFile(r"temporary.mol")
    all_element=['Al','Si','As','Sc','Tl','Sb','Ti','Pb','Te','V','Po','Sb','In','Au','Ga','Ag','Ge','Pt','Bi']
    coor1 = Atom_coor
    # create a molecule from a SMILES string
    try:
        # get the ring information of the molecule
        ring_info = mol.GetRingInfo()
    except:
        return False
    # get the ring systems of the molecule
    ring_systems_ = ring_info.AtomRings()
    ring_systems = get_ring_systems(list(ring_systems_))
    l_=[]
    chain_idx=[]
    for idx in idx_:
        a=idx_.index(idx)
        # loop over the rings in the SSSR
        atoms_idx=[]
        for ring in ring_systems:
            if idx in ring:
                atoms_idx = ring
        if atoms_idx==[]:
            chain_idx.append(idx)
            continue
        break_indices=[]
        for i in atoms_idx:
            neighbours=mol.GetAtomWithIdx(i).GetNeighbors()
            for neighbour in neighbours:
                break_idx=[]
                if neighbour.GetIdx() not in atoms_idx:
                    bond_idx1 = mol.GetBondBetweenAtoms(i,neighbour.GetIdx()).GetIdx()
                    break_idx.append([bond_idx1,i,neighbour.GetIdx()])
                    neighbours1=neighbour.GetNeighbors()
                    for neighbour_ in neighbours1:
                        if neighbour_.IsInRing()==False:
                            bond_idx2 = mol.GetBondBetweenAtoms(neighbour.GetIdx(),neighbour_.GetIdx()).GetIdx()
                            break_idx.append([bond_idx2,neighbour.GetIdx(),neighbour_.GetIdx()])          
                if break_idx !=[]:
                    break_indices.append(break_idx)
        indices_all=break_indices
        break_indices_=tuple(break_indices)
        combs=[]
        if len(break_indices_)>1:
            for i in break_indices:
                indices_all.remove(i)
                for j in break_indices:
                    if i != j:
                        comb = list(itertools.product(i, j))
                        combs.append(comb)
            for c in combs:
                l=[]
                for i in c:
                    f_=[Chem.MolToSmiles(mol)]
                    fragments = Chem.FragmentOnBonds(mol, [i[0][0],i[1][0]],dummyLabels=[[1,1],[2,2]])
                    f = Chem.MolToSmiles(fragments).split('.')
                    x=1
                    for j_ in f:
                        mol_j=Chem.MolFromSmiles(j_)
                        if mol_j==None:
                            x=-1
                            break
                        # delete small fragments
                        # if mol_j.GetNumAtoms()<5:
                        #     x=-1
                        #     break
                    if x==-1:
                        continue
                    matches=[]
                    for j in f:  
                        mol_f=Chem.MolFromSmiles(f_[0])
                        if j.count('*')==2:
                            j_mol=delete2dummy(j,True)
                            try:
                                mol_=AllChem.AssignBondOrdersFromTemplate(mol, mol_f)
                                for match in mol.GetSubstructMatches(j_mol)[0]:
                                    if match in idx_:
                                        matches.append(match)                                     
                                j_1=sidechain_removal(j,interaction)
                            except:
                                j_1=sidechain_removal(j,interaction)
                            f_.append(j_1)    
                        else:
                            f_.append(j)
                    for k in i:
                        coor=get_coordinate_with_idx(lig_pdb_path,k[2])
                        f_.append(coor)  
                    D_A=[]
                    if len(matches)>=2:
                        for match in matches:
                            b=idx_.index(match)
                            f_.extend(coor1[b])
                            D_A.append(D_A_[b])
                    else:
                        f_.extend(coor1[a])
                        D_A.append(D_A_[a])
                    f_.append(D_A)
                    l.append(f_)
                l_.append(l)
        if len(break_indices_)==1:
            print ('aaa')
            matches=[]
            f_=[Chem.MolToSmiles(mol)]
            for only_idx in break_indices_[0]:
                fragments = Chem.FragmentOnBonds(mol, [only_idx[0]],dummyLabels=[[3,3]])
                f = Chem.MolToSmiles(fragments).split('.')
                for j in f:
                    mol_f=Chem.MolFromSmiles(f_[0])
                    mol_=AllChem.AssignBondOrdersFromTemplate(mol, mol_f)
                    j_mol=delete2dummy(j)
                    matches_=mol.GetSubstructMatches(j_mol)
                    if idx in matches_[0]:#j:replace fragments
                        for match in matches_[0]:
                            if match in idx_:
                                matches.append(match)
                        for k in f:
                            if k!=j:
                                coor=get_coordinate_with_idx(lig_pdb_path,only_idx[2])
                                l=f_+[j,k,coor]
                                D_A=[]
                                if len(matches)>=2:                                        
                                    for match in matches:
                                        b=idx_.index(match)
                                        l.extend(coor1[b])
                                        D_A.append(D_A_[b])
                                else:
                                    l.extend(coor1[a])
                                    D_A.append(D_A_[a])
                                l.append(D_A)
                                l_.append([l])
                                print (l_)
    
    for ch in chain_idx:
        mol_=mol
        combs_=()
        a_=idx_.index(ch)
        break_indices=[]
        try:
            neighbours=mol_.GetAtomWithIdx(ch).GetNeighbors()
            ch_symbol=mol.GetAtomWithIdx(ch).GetSymbol()
        except:
            print ([])
        nei_inring=0
        for neighbour in neighbours:
            if neighbour.IsInRing()==True:
                nei_inring=1
                break
        for neighbour in neighbours:
            break_idx=[]
            if neighbour.IsInRing()==True:
                for ring in ring_systems:
                    if neighbour.GetIdx() in ring:
                        ring_indices = ring
                        for j in ring_indices:
                            neighbours_=mol_.GetAtomWithIdx(j).GetNeighbors()
                            for neighbour_ in neighbours_:
                                if neighbour_.GetIdx() not in ring_indices:
                                    if neighbour_.GetIdx()!=ch:
                                        bond_idx1 = mol_.GetBondBetweenAtoms(j,neighbour_.GetIdx()).GetIdx()
                                        break_idx.append([bond_idx1,j,neighbour_.GetIdx()])
                combs_ = combs_+tuple(itertools.product(break_idx, break_idx))
            else:          
                bond_idx1 = mol_.GetBondBetweenAtoms(ch,neighbour.GetIdx()).GetIdx()
                neighbours1=neighbour.GetNeighbors()
                for neighbour_ in neighbours1:
                    if neighbour_.IsInRing()==False:
                        bond_idx2 = mol_.GetBondBetweenAtoms(neighbour.GetIdx(),neighbour_.GetIdx()).GetIdx()
                        if bond_idx2!=bond_idx1:
                            break_idx=[bond_idx2,neighbour.GetIdx(),neighbour_.GetIdx()]
                    if break_idx !=[]:
                        break_indices.append(break_idx)
            indices_all=break_indices
            combs_ = combs_+tuple(itertools.product(break_indices, break_indices))
        combs=[]
        for i in combs_:
            if i[0]!=i[1]:
                combs.append(i)
        l=[]
        for c in combs:
            f_=[Chem.MolToSmiles(mol)]
            fragments = Chem.FragmentOnBonds(mol_, [c[0][0],c[1][0]],dummyLabels=[[1,1],[2,2]])
            f = Chem.MolToSmiles(fragments).split('.')
            x=1
            for j_ in f:
                mol_j=Chem.MolFromSmiles(j_)
                if mol_j==None:
                    x=-1
                    break
                # if mol_j.GetNumAtoms()<5:
                #     x=-1
                #     break
            if x==-1:
                continue
            matches=[]
            for j in f:         
                mol_f=Chem.MolFromSmiles(f_[0])
                if j.count('*')==2:
                    j_mol=delete2dummy(j,True)                
                    try:
                        mol__=AllChem.AssignBondOrdersFromTemplate(mol_, mol_f)
                        for match in mol_.GetSubstructMatches(j_mol)[0]:
                            if match in idx_:
                                matches.append(match)
                        j_1=sidechain_removal(j,interaction)
                    except:
                        j_1=sidechain_removal(j,interaction)
                    f_.append(j_1)
                else:
                    f_.append(j)
            for k in c:
                coor=get_coordinate_with_idx(lig_pdb_path,k[2])
                f_.append(coor)
            D_A=[]
            if len(matches)>=2:
                for match in matches:
                    b_=idx_.index(match)
                    f_.extend(coor1[b_])
                    D_A.append(D_A_[b_])
            else:
                f_.extend(coor1[a_])
                D_A.append(D_A_[a_])
            f_.append(D_A)
            l.append(f_)
        l_.append(l)
    return l_


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', type=str, default = None)
    args = parser.parse_args()
    list_res=['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
    directory_pdb=r"BTK_covint_complex"  
    directory_rep=r"BTK_covint_report"  
    directory_lig=r"BTK_covint_lig"     
    Int=args.mode
    output_path=os.path.join(r"warhead_substitution",'BTK_covint_'+Int+'.txt')          
    with open (output_path,'w+') as t:
        pdbpath=[]
        reppath=[]
        ligpath=[]
        for filename in os.listdir(directory_pdb):
            filepath = os.path.join(directory_pdb, filename)
            if os.path.isfile(filepath):
                pdbpath.append(filepath)
        for filename in os.listdir(directory_rep):
            filepath = os.path.join(directory_rep, filename)
            if os.path.isfile(filepath):
                reppath.append(filepath)
        for filename in os.listdir(directory_lig):
            filepath = os.path.join(directory_lig, filename)
            if os.path.isfile(filepath):
                ligpath.append(filepath)
        for i in range(0,len(pdbpath)):
            t.write(pdbpath[i])
            t.write('\n')
            report_path=reppath[i]
            ori_pdb_path=pdbpath[i]
            lig_pdb_path=ligpath[i]
            gate = False
            with open(lig_pdb_path,'r') as ligs:
                lines=ligs.readlines()
                for line in lines:
                    line__=line.split(' ')
                    line_=[x.strip() for x in line__ if len(x.strip()) > 0]
                #     if line_[0] =='ATOM':
                #         if line_[3] in list_res:
                #             gate = True
                #             break
                # if gate:
                #     continue
                print (report_path)
                int_list=['HI','HB','PS','PC']
                mol_pdb=AllChem.MolFromPDBFile(lig_pdb_path)
                aromatic=0
                try:
                    for atom in mol_pdb.GetAtoms():
                        if atom.GetIsAromatic()==True:
                            aromatic=1                
                            break
                    if aromatic==1:
                        Chem.Kekulize(mol_pdb,clearAromaticFlags=True)
                        mol_oo=mol_pdb
                        mol_o=mol_oo
                    else:
                        mol_oo=mol_pdb
                        mol_o=mol_pdb
                except:
                    continue
                for interaction in int_list:
                    giac=get_interaction_atom_coor(ori_pdb_path,report_path,lig_pdb_path,interaction)
                    if giac==False:
                        indices=[]
                    else:
                        int_coor=giac[0]
                        indices=get_idx_from_pdb(int_coor,lig_pdb_path)
                        if interaction=='HI':
                            symbol_list=['C']
                            replace_list=['Au']
                            replace_list_=['Ga']
                        if interaction=='HB':
                            symbol_list=['B','C','N','O','F','P','S','Cl','Se','Br','I']
                            replace_list=['Al','Si','As','Sc','Tl','Sb','Ti','Pb','Te','V','Po']
                            replace_list_=['Al','Sb','In','Sc','Tl','Sb','Ti','Pb','Te','V','Po']
                        if interaction=='PS':
                            symbol_list=['C']
                            replace_list=['Ag']
                            replace_list_=['Ge']
                        if interaction=='PC':
                            symbol_list=['C']
                            replace_list=['Pt']
                            replace_list_=['Bi']
                        # create a molecule from a SMILES string                
                        rw_o=Chem.RWMol(mol_o)
                        if indices!=[]:
                            for idx in indices:
                                if mol_oo.GetAtomWithIdx(idx).GetIsAromatic()==True:
                                    all_ring=[]
                                    for ring_indices in mol_oo.GetRingInfo().AtomRings():
                                        if idx in ring_indices:
                                            all_ring+=ring_indices
                                    if all_ring.count(idx)>1:
                                        for i in all_ring:
                                            if all_ring.count(i)==1:
                                                if i!=idx:
                                                    idx=i
                                                    break
                                    for ring_indices in mol_oo.GetRingInfo().AtomRings():
                                        if idx in ring_indices:
                                            if len(ring_indices)==5:
                                                idx_symbol=mol_oo.GetAtomWithIdx(idx).GetSymbol()
                                                idx_replace=replace_list[symbol_list.index(idx_symbol)]
                                                rw_o.ReplaceAtom(idx,Chem.Atom(idx_replace))
                                            else:
                                                idx_symbol=mol_oo.GetAtomWithIdx(idx).GetSymbol()
                                                idx_replace=replace_list_[symbol_list.index(idx_symbol)]
                                                rw_o.ReplaceAtom(idx,Chem.Atom(idx_replace))
                                else:
                                    idx_symbol=mol_oo.GetAtomWithIdx(idx).GetSymbol()
                                    idx_replace=replace_list[symbol_list.index(idx_symbol)]
                                    rw_o.ReplaceAtom(idx,Chem.Atom(idx_replace))
                            mol=rw_o.GetMol() 
                            mol_o=mol 
                giac_=get_interaction_atom_coor(ori_pdb_path,report_path,lig_pdb_path,Int)
                if giac_==False:
                    continue
                Atom_coor=giac_[0]
                D_A=giac_[2]
                idx=get_idx_from_pdb(Atom_coor,lig_pdb_path)
                Chem.MolToMolFile(mol_o,r"temporary.mol")
                if Int=='HB':
                    D_A=HA_detector_lig(D_A,idx,lig_pdb_path)
                pairs=get_pairs(lig_pdb_path,idx,Atom_coor,mol,Int,D_A)
                if pairs==False:
                        continue
                for i in range(len(pairs)):
                    for j in pairs[i]:
                        if isinstance(j, list):
                            for k in j:
                                if isinstance(k, list):
                                    for k_i in range(len(k)):
                                        t.write(k[k_i])
                                        if k_i!=len(k)-1:
                                            t.write(',')
                                    t.write(' ')
                                else:
                                    k_=str(k)
                                    t.write(k_)
                                    t.write(' ')                        
                        else:
                            t.write(j)
                            t.write(' ')
                        t.write('\n')
                        
                    
                
