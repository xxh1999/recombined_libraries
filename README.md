![Re_F1](https://github.com/xxh1999/reconstructed_libraries/assets/94356070/86ffa86b-29d9-4cdb-89de-8974da1d3892)
# A split-and-mix computational approach for the design of novel (non)covalent virtual screening libraries.
## requirements
python>=3.6  
biopandas  
multiprocessing  
numpy  
pandas=1.1.2  
rdkit=2020.03.3  

## noncovalent recombination
#### 1.Get complex files in pdb format and rectify inappropriate structures in it.
#### 2.Use PLIP to analysis noncovalent interactions in PDB files and generate reports. 
``git clone https://github.com/pharmai/plip.git``  
``call report_generation.bat``
#### 3.Extract noncovalent replacement fragments with report and PDB files.
``python noncovalent.py --mode 'HB/HI/PC/PS'``  
#### 4.Select appropriate replacement fragments for each pair of retained fragments.  
``python RMSD_calculation.py --provider path_to_frag_provider --acceptor path_to_frag_acceptor``  
#### 5.Fill in every vacancy in the original molecules.  
``python final_combine.py``  

## covalent recombination
#### 1.Get complex files in pdb format and rectify inappropriate structures in it.
#### 2.Use PLIP to analysis noncovalent interactions in PDB files and generate reports. 
``git clone https://github.com/pharmai/plip.git``  
``call report_generation.bat``  
#### 3.Extract covalent adducts from complex PDB files.
``python pdb_covalent_extract.py``
#### 4.Extract covalent fragments with report and PDB files.
##### 1) Extract covalent fragments from existing adducts.
``python covalent.py``   
##### 2) generate covalent fragments with delinker.
``git clone https://github.com/oxpig/DeLinker.git``  
``python delinker_generate.py``   
``python delinker_combine.py`` 
#### 5.Extract noncovalent replacement fragments with report and PDB files.
``python noncovalent.py --mode 'HB/HI/PC/PS'``  
#### 6.Select appropriate replacement fragments for each pair of retained fragments.  
``python RMSD_calculation.py --provider path_to_frag_provider --acceptor path_to_frag_acceptor``  
#### 7.Fill in every vacancy in the original molecules.  
``python final_combine.py``  
#### 8.Transform adducts to covalent inhibitors.  
``python warhead2inhibitor.py``  
``python warhead_recover.py``  
