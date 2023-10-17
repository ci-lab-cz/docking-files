## Complexes prepared for docking

This is a repository of protein-ligand complexes prepared for docking in PDB nad PDBQT formats.

The detailed pipeline is described in the document `Target-prepare_desgn.docx`. Please follow it if you contribute a new structure.    

The repository is composed of directories for individual proteins, where every directory has the following structure:
```
UNIPROT_ID                 (protein)
|-PDB_ID_1                 (particular complex)
|  |-CHAIN_ID              (particular chain)
|    |-.._protein.pdb      (prepared protein)
|    |-.._protein.pdbqt    (prepared protein)
|    |-.._ligand.pdb       (prepared ligand)
|    |-.._liagnd.mol(sdf)  (prepared ligand)
|    |-..                  (other files can be presented)  
|-PDB_ID_2
|  |-CHAIN_ID
|    |-.._protein.pdb
|    |-.._protein.pdbqt
|    |-.._ligand.pdb
|    |-.._liagnd.mol(sdf)
|-box.txt                  (grid box, may have different names)
|-info.txt                 (auxillary information)
```
The directory `scripts` contains some useful scripts for protein preparation.
 