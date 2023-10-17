# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 15:06:55 2023

@author: yk5
"""

import numpy as np
from collections import OrderedDict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
import multiprocessing

def dis_cal(l1,l2):
    dis=(float(l1[0])-float(l2[0]))**2+(float(l1[1])-float(l2[1]))**2+(float(l1[2])-float(l2[2]))**2
    return dis

def rough_cal(l1,l2):
    cal=dis_cal(l1[0],l1[1])-dis_cal(l2[0],l2[1])+dis_cal(l1[0],l1[2])-dis_cal(l2[0],l2[2])+dis_cal(l1[2],l1[1])-dis_cal(l2[2],l2[1])
    return abs(cal)
    
def RMSD_cal(l1,l2):   
    # Define the two sets of coordinates
    coordset1 = np.array(l1)
    coordset2 = np.array(l2)   
    # Calculate the center of mass of each set of coordinates
    com1 = np.mean(coordset1, axis=0)
    com2 = np.mean(coordset2, axis=0)   
    # Translate the sets of coordinates to their respective centers of mass
    coordset1 -= com1
    coordset2 -= com2    
    # Calculate the covariance matrix
    C = np.dot(coordset1.T, coordset2)    
    # Singular Value Decomposition
    U, s, Vt = np.linalg.svd(C)  
    # Calculate the rotation matrix
    R = np.dot(U, Vt)    
    # Apply the rotation and calculate the RMSD
    coordset1_rotated = np.dot(coordset1, R)
    diff = coordset1_rotated - coordset2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))    
    return rmsd

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



def process(lines,lines_h,output):
    choose_list=[70,70,60,40,20,10,7,3,2]
    # choose_list=[1,2,2,2,2,10,7,3,2]#PC
    # choose_list=[1,1,3,4,3,8]#PS
    with open (r"BTK_covint_HB_suggest_frags_new_"+str(output)+".txt",'w+') as f:
        ls_g=[]
        ls=[]
        ls_h=[]
        ls_=[]
        ls1=[]
        ls2=[]
        l_rmsd=[]     
        dis=[]
        dis_h=[]    
        idx_3=[]
        idx3=[]
        pat_=[]
        pat_h_=[]
        for i in range(len(lines_h)):
            coordinate=[]
            if '3*' in lines_h[i]:
                idx_3.append(i)
            line_=lines_h[i].split(' ')
            #delete small frags
            if len(line_)>4:
                try:
                    for coor in line_:
                        if ',' in coor:
                            if coor[0].isalpha()==False:
                                coor_list=coor.split(',')
                                coordinate.append(list(map(float, coor_list)))
                    pat_h=line_[-1].strip('\n')
                except:
                    continue
                pat_h_.append(pat_h)
                ls1.append(coordinate)
                ls_h.append(lines[i])
                dis_h.append(line_[0])#replaced fragments
        lss1=tuple(ls1)
        for i in range(len(lines)):
            coordinate=[]
            if '3*' in lines[i]:
                idx3.append(i)
            line_=lines[i].split(' ')
            #delete small frags
            if len(line_)>4:
                try:
                    for coor in line_:
                        if ',' in coor:
                            if coor[0].isalpha()==False:
                                coor_list=coor.split(',')
                                coordinate.append(list(map(float, coor_list)))
                    pat=line_[-1].strip('\n')
                except:
                    continue
                pat_.append(pat)
                ls.append(coordinate)#coordinate of replaced fragments
                ls_.append(lines[i])
                dis.append(line_[0])#distance of replaced fragments
        lss=tuple(ls)
        #被替换片段
        for i in ls1:              
            idxi=ls1.index(i)
            print (idxi)
            line_listi=lines_h[idxi].split(' ')
            dis_i=dis_h[idxi]
            try:
                for i_ in line_listi:
                    molecule=line_listi[1]
                    if '1*' in i_ and '2*' in i_:
                        frag_ori=i_
                        threshold=choose_list[len(i)-3]
                    if '3*' in i_:
                        frag_ori=i_
                        threshold=choose_list[len(i)-2]
                        break
            except:
                continue
            d={}
            for di in dis:
                if float(dis_i)-1<float(di):
                    idx1=dis.index(di)
                    break
            for di in dis:
                if float(dis_i)+1<float(di):
                    idx2=dis.index(di)
                    break
                else:
                    idx2=len(dis)-1
            #替换片段
            for j in ls[idx1:idx2]:
                idx_1=ls1.index(i)
                idx_=ls.index(j)
                if idx_1 not in idx_3:
                    if idx_ not in idx3:             
                        if len(i)==len(j):
                            judge=1
                            for pat in zip(pat_[idx_].split(','),pat_h_[idx_1].split(',')):
                                if (pat[0]=='a' and pat[1]=='b') or (pat[0]=='b' and pat[1]=='a'):
                                    judge=0
                                    break
                            if judge==0:
                                continue
                            R=RMSD_cal(i,j)
                            d[R]=[j,lss.index(j)]
                if idx_1 in idx_3:
                    if idx_ not in idx3: 
                        if len(i)+1==len(j):
                            judge=1
                            for pat in zip(pat_[idx_].split(','),pat_h_[idx_1].split(',')):
                                if (pat[0]=='a' and pat[1]=='b') or (pat[0]=='b' and pat[1]=='a'):
                                    judge=0
                                    break
                            if judge==0:
                                continue
                            R=RMSD_cal(i,[j[0]]+j[2:])
                            d[R]=[j,lss.index(j)]
                    if idx_ in idx3: 
                        if len(i)==len(j):
                            judge=1
                            for pat in zip(pat_[idx_].split(','),pat_h_[idx_1].split(',')):
                                if (pat[0]=='a' and pat[1]=='b') or (pat[0]=='b' and pat[1]=='a'):
                                    judge=0
                                    break
                            if judge==0:
                                continue
                            R=RMSD_cal(i,j)
                            d[R]=[j,lss.index(j)]
            values = [d[k] for k in sorted(d.keys())]
            bm_scaffolds=[]
            v=0 
            for value in values:
                idx=value[1]
                pat=idx
                line_listj=lines[idx].split(' ')
                for j_ in line_listj:
                    if '1*' in j_ and '2*' in j_:
                        frag_rep=j_
                        mol=Chem.MolFromSmiles(frag_rep)
                        bm_scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
                        if bm_scaffold not in bm_scaffolds:
                            bm_scaffolds.append(bm_scaffold)
                            f.write(molecule+' '+frag_ori+' '+frag_rep)
                            f.write('\n')
                            if len(bm_scaffolds)==threshold:
                                v=1
                    if '3*' in j_:
                        frag_rep=j_
                        mol=Chem.MolFromSmiles(frag_rep)
                        bm_scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
                        if bm_scaffold not in bm_scaffolds:
                            bm_scaffolds.append(bm_scaffold)
                            f.write(molecule+' '+frag_ori+' '+frag_rep)
                            f.write('\n')
                            if len(bm_scaffolds)==threshold:
                                v=1
                        break
                if v==1:
                    break
        return ('finish')
        
if __name__ == '__main__':                 
    with open (r"KLIF_HB_onlydis_ring.txt",'r') as t: #t:frag provider
        with open (r"BTK_covint_HB_onlydis_ring.txt",'r') as h: #h:frag_Accepter       
            lines=t.readlines()
            lines_h=h.readlines()
            sub_lists = [lines_h[i:i + len(lines_h) // 10] for i in range(0, len(lines_h), len(lines_h) // 10)]
            pool = multiprocessing.Pool(processes=10)
            results = []
            for sub_list in sub_lists:
                args = (lines,sub_list, sub_lists.index(sub_list))
                result = pool.apply_async(process, args=args)
                results.append(result)
            
            results = [result.get() for result in results]







