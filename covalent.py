# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 19:17:04 2022

@author: yk5
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd
import os

def get_atom_coor(mol,atom):
    idx=atom.GetIdx()
    coor=list(mol.GetConformer().GetAtomPosition(idx))
    return coor
               
def Mol_Conn(m1_,s2,idxA,attachment):
    atomA = m1_.GetAtomWithIdx(int(idxA))
    neighbours = atomA.GetNeighbors()
    list=[]
    for neighbour in neighbours:
        bond = m1_.GetBondBetweenAtoms(neighbour.GetIdx(),int(idxA))
        mol1_f = Chem.FragmentOnBonds(m1_, [bond.GetIdx()])
        smi_A = Chem.MolToSmiles(mol1_f)
        if '.' in smi_A:
            smi_list=smi_A.split('.')
            for smi in smi_list:
                if 'S' not in smi:
                    list.append(smi)
    list1 = sorted(list, key=lambda i: len(i))
    frag=list1[-1]
    m1=Chem.MolFromSmiles(frag)
    conn = []
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

def amo_filter(pdb,typ):
    for line in pdb:
        list = line.split(' ')
        list_ = [x.strip() for x in list if len(x.strip()) > 0]
        if list_[0]=='TER':
            if typ=='cys':
                if list_[2]=='CYS':
                    return True
            if typ=='ser':
                if list_[2]=='SER':
                    return True
            else:
                return False
    return False

def dis_filter(pdbA,pdbB,molA,typ,noncovalent=False):
    if noncovalent==False:
        total=[]  
        if typ=='ser':
            amo='OCC(N)C(=O)'
        if typ=='cys':
            amo='SCC(N)C(=O)'        
        print (amo)
        mol_amo=Chem.MolFromSmiles(amo)
        match_indices = molA.GetSubstructMatches(mol_amo)
        for match in match_indices[0]:
            for atom in molA.GetAtomWithIdx(match).GetNeighbors():
                if atom.GetIdx() not in match_indices[0]:
                    idx_S=match
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
        idx_list_=[idxA]
        frag_=[frag]
        # for line in pdbA:
        #     list = line.split(' ')
        #     list_ = [x.strip() for x in list if len(x.strip()) > 0]
        #     if len(list_)>10:
        #         total.append(list_)
        for line in pdbA:
            list = line.split(' ')
            list_ = [x.strip() for x in list if len(x.strip()) > 0]
            if len(list_)>10:              
                Xs=float(line[30:38].strip(' '))
                Ys=float(line[38:46].strip(' '))
                Zs=float(line[46:54].strip(' '))
                if line[13:15]=='SG':
                    if typ=='cys':
                        break
                if line[13:15]=='OG':
                    if typ=='ser':
                        break
        for line in pdbA:
            list = line.split(' ')
            list_ = [x.strip() for x in list if len(x.strip()) > 0]
            if len(list_)>10:
                if list_[1]==str(idxA+2):
                    idx = list_[1]
                    Xa = float(line[30:38].strip(' '))
                    Ya = float(line[38:46].strip(' '))
                    Za = float(line[46:54].strip(' '))
                    disA =[((Xs-Xa)**2+(Ys-Ya)**2+(Zs-Za)**2)**0.5]
# =============================================================================
#         else:
#             if molA.GetAtomWithIdx(idxA_).IsInRing()==True:
#                 dis_final = atom_dis[atom_dis.index(disA)-1] 
#                 idx_A = atom_idx[atom_dis_fix.index(dis_final)]
#                 print (idx_A)
#                 break
#                 for i in total:
#                     if i[1]==idx_A:
#                         idxA=total.index(i)
#                         break
# =============================================================================

# noncovalent is a list contains the coor of S atom
    if noncovalent:
        frag=Chem.MolToSmiles(molA)
        coor_list=[]
        dis_list=[]    
        idx_list=[]
        idx_list_=[]
        coor_list_=[]
        for atom in molA.GetAtoms():
            if atom.GetTotalNumHs()>0:
                coor=get_atom_coor(molA,atom)
                dis=((float(coor[0])-noncovalent[0])**2+(float(coor[1])-noncovalent[1])**2+(float(coor[2])-noncovalent[2])**2)**0.5
                idx_list.append(atom.GetIdx())
                coor_list.append(coor)
                dis_list.append(dis)
        for d in range(len(sorted(dis_list))):
            coor_list_.append(coor_list[dis_list.index(dis_list[d])])
            idx_list_.append(idx_list[dis_list.index(dis_list[d])])
            if d==4:
                break
        disA=sorted(dis_list)[:len(coor_list_)]    
    idxB_list=[]
    for disa in disA:
        for line in pdbB:
            list = line.split(' ')
            list_ = [x.strip() for x in list if len(x.strip()) > 0]
            if len(list_)>10:
                Xs_=float(line[30:38].strip(' '))
                Ys_=float(line[38:46].strip(' '))
                Zs_=float(line[46:54].strip(' '))
                if line[13:15]=='SG':
                    if typ=='cys':
                        break
                if line[13:15]=='OG':
                    if typ=='ser':
                        break
        HET_line = []
        dis_list=[]
        for line in pdbB:
            list = line.split(' ')
            list_ = [x.strip() for x in list if len(x.strip()) > 0]
            if len(list_)>10:
                if list_[0]=='HETATM':
                    HET_line.append(list_)
                    X=float(line[30:38].strip(' '))
                    Y=float(line[38:46].strip(' '))
                    Z=float(line[46:54].strip(' '))
                    dis=((Xs_-X)**2+(Ys_-Y)**2+(Zs_-Z)**2)**0.5-disa
                    dis_list.append(abs(dis))
        dis_list_=sorted(dis_list)
        idxB=[]
        for i in dis_list_:
            if i<1.4:
                idx_min=dis_list.index(i)
                idxB.append(HET_line[idx_min][1])
        idxB_list.append(idxB)
    if noncovalent==False:
        idx_list_=idx_list_*len(idxB_list)
        frag=frag_*len(idxB_list)
    return [idx_list_,idxB_list,frag]

def replacement(smi,patt,repl):
    mol=Chem.MolFromSmiles(smi)
    mol_patt=Chem.MolFromSmarts(patt)
    mol_repl = Chem.MolFromSmiles(repl)
    rms = AllChem.ReplaceSubstructs(mol, mol_patt, mol_repl)
    return rms[0]

def warhead_sub(PDB_fileB,idxB_list,typ):
    if typ=='ser':
        amo='OCC(N)C(=O)'
    if typ=='cys':
        amo='SCC(N)C(=O)'
    mol_amo=Chem.MolFromSmiles(amo)
    for idxB in idxB_list:
        molB = AllChem.MolFromPDBFile(PDB_fileB)
        if molB!=None:
            atomB = molB.GetAtomWithIdx(int(idxB)-2)
            if atomB.IsInRing()==True:
                continue
            if atomB.IsInRing()!=True:
                neighbours = atomB.GetNeighbors()
                list=[]
                for neighbour in neighbours:
                    bond = molB.GetBondBetweenAtoms(neighbour.GetIdx(),int(idxB)-2)
                    mol1_f = Chem.FragmentOnBonds(molB, [bond.GetIdx()])
                    smi_B = Chem.MolToSmiles(mol1_f)
                    if '.' in smi_B:
                        smi_list=smi_B.split('.')
                        for smi in smi_list:
                            if Chem.MolFromSmiles(smi).GetSubstructMatches(mol_amo)!=():
                                list.append(smi)
                list1 = sorted(list, key=lambda i: len(i))
                try:
                    frag=list1[1]
                    return frag
                except:
                    continue              
    return None


lig_path=r"covint_lig"
lig_path_=r"covint_lig"
direct_connect=False
noncovalent=False
pattern='cys'
with open(r"C:\Users\yk5\Desktop\precious_f_non_2.txt",'w+') as g:
    res_nameA=os.listdir(lig_path)
    res_nameB=os.listdir(lig_path_)
    pdbs=[]
    pdbsA=[]
    pdb_file=[]
    for i in res_nameA:
        with open (os.path.join(lig_path,i),'r') as t:
            pdbsA.append(t.readlines())  
    for i in res_nameB:
        with open (os.path.join(lig_path,i),'r') as t:
            pdbs.append(t.readlines())      
    for j in pdbsA:
        J=pdbsA.index(j)
        pdbA=j
        for k in pdbs:
            K = pdbs.index(k)
            pdbB = k
            if (amo_filter(pdbA,pattern) or noncovalent) and amo_filter(pdbB,pattern):
                print ('a')
                # try:
                molA = AllChem.MolFromPDBFile(os.path.join(lig_path,res_nameA[J]))
                if molA==None:
                    continue
                try:
                    idxABs = dis_filter(pdbA,pdbB,molA,pattern)
                except:
                    continue
                for idxAB in zip(idxABs[0],idxABs[1],idxABs[2]):                           
                    idxA = idxAB[0]
                    idxB = idxAB[1] 
                    fragA = idxAB[2] 
                    fragB = warhead_sub(os.path.join(lig_path_,res_nameB[K]),idxB,pattern)
                    if fragB != None:                           
                        if direct_connect==True:
                            conn = Mol_Conn(molA,frag,int(idxA),'*')
                            if conn!=[]:
                                g.write(res_name[J])
                                g.write(conn[0])
                                g.write('\n')
                        if direct_connect==False:
                            if fragA and fragB:
                                g.write(Chem.MolToSmiles(molA)+' ')
                                if noncovalent:
                                    g.write(str(idxA)+' ')
                                    g.write(fragB)
                                    g.write('\n')
                                else:
                                    g.write(fragA+' ')
                                    g.write(fragB)
                                    g.write('\n')
