# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 20:20:53 2023

@author: yk5
"""

import os
from rdkit import Chem
from rdkit.Chem import AllChem

path=r'ser_delinker'
path_sdf=r'covint_lig_sdf'
files=os.listdir(path)
du = Chem.MolFromSmiles('*')
with open (r'delinker_generate_ser.txt','w+') as f:
    for file in files:
        print (file[-8:-4])
        with open (os.path.join(path,file),'r') as t:
            suppl=Chem.SDMolSupplier(os.path.join(path_sdf,file[-8:-4]+'.sdf'))
            mol=suppl[0]
            lines=t.readlines()
            for line in lines:
                line_=line.split(' ')
                smi_=line_[2]
                for i in line_[0].split('.'):
                    if i!='[*]SCC([NH3+])C=O':
                        frag=i
                mol_frag=Chem.MolFromSmiles(frag)
                c1h=AllChem.ReplaceSubstructs(mol_frag,du,Chem.MolFromSmiles('[H]'),True)[0]
                c1h = Chem.RemoveHs(c1h)
                mol_rep=Chem.MolFromSmiles(smi_)
                frag_ori=Chem.ReplaceSubstructs(mol, c1h, Chem.MolFromSmiles("[*:1]"), replaceAll=True)
                frag_rep=Chem.ReplaceSubstructs(mol_rep, c1h, Chem.MolFromSmiles("[*:1]"), replaceAll=True)
                f.write(Chem.MolToSmiles(mol)+' ')
                f.write(Chem.MolToSmiles(frag_ori[0])+' ')
                f.write(Chem.MolToSmiles(frag_rep[0])+'\n')
            
        