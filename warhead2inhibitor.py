# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 09:14:53 2023

@author: yk5
"""

from rdkit import Chem
from rdkit.Chem import AllChem

with open(r"ser_sub.txt",'w+') as g:
    with open(r"SMIRKS_ser.txt",'r') as f:
        lines_smirks=f.readlines()
        with open (r"covint_ser.txt",'r') as t:
            lines=t.readlines()
            l=[]
            for i in range(1,len(lines)):
                print (i)
                line_=lines[i].split(' ')
                if lines[i].split(' ')[1]==lines[i-1].split(' ')[1]:
                    original_smiles=lines[i].split(' ')[2]
                    for smirks in lines_smirks:        
                        original_mol = Chem.MolFromSmiles(original_smiles)
                        reaction = AllChem.ReactionFromSmarts(smirks)
                        try:
                            transformed_mols = reaction.RunReactants((original_mol,))
                        except:
                            continue
                        for transformed_mol in transformed_mols:
                            transformed_smiles = Chem.MolToSmiles(transformed_mol[0])
                            l.append(transformed_smiles)
                if lines[i].split(' ')[1]!=lines[i-1].split(' ')[1]:
                    l_set=set(l)
                    g.write(lines[i-1].split(' ')[1]+' ')
                    cc=0
                    for j in l_set:
                        if '*' in j:
                            g.write(j+' ')
                    g.write('\n')

    
        
