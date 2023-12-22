![Re_F1](https://github.com/xxh1999/reconstructed_libraries/assets/94356070/86ffa86b-29d9-4cdb-89de-8974da1d3892)
# A split-and-mix computational approach for the design of novel (non)covalent virtual screening libraries.
## noncovalent recombination
### data preprocessing
#### 1.Get complex files in pdb format and rectify inappropriate structures in it.
#### 2.Use PLIP to analysis noncovalent interactions in PDB files and generate reports. 
``git clone https://github.com/pharmai/plip.git``  
``call report_generation.bat``
#### 3.Extract replacement fragments with report and PDB files.
``python noncovalent.py``  
#### 4.Select appropriate replacement fragments for each pair of retained fragments.  
``python RMSD_calculation.py``  
#### 5.Fill in every vacancy in the original molecules.  
``python final_combine.py``  
