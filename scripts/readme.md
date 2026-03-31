#### Note: 
* Create *~/docking directory*  
``
mkdir ~/docking
``
* Save all dependencies and scripts to the *~/docking* directory
## convert2pdb2pdbqt.pbs
#### script for the preparation of ligands
__Dependency__  
* RDKit  
```
conda create -n my-rdkit-env
source activate my-rdkit-env
conda install -c conda-forge rdkit

```
* Autodock vina  
http://mgltools.scripps.edu/downloads  
download and install to the *~/docking/ directory*

* Chemaxon  
https://chemaxon.com/

* molchemaxon2pdb.py  
https://github.com/DrrDom/rdkit-scripts  
download to the *~/docking/ directory*

__Input:__ *.smi or *.sdf 
 
To keep things simple, the name should be *ligand.smi / ligand.sdf* or, if you don't want it, you can make small changes to the script files.  

__Output:__ new folder (folder name is equal to smi/sdf file name - default: ./ligand) which contains ligand pdbqt files  

__Run__: 
```
qsub ~/docking/convert2pdb2pdbqt.pbs -v path=$(pwd)
```
## run-vina.pbs
__Dependency__  
* Vina   
http://vina.scripps.edu/download.html  
download, extract and copy  *vina* file to the *~/docking/ directory*

__Input:__ 
* __./target__ folder which contains *pdbqts of targets* and *config.txt*
* __./ligand__ folder which contains *pdbqts of ligands*  

This script uses 32 ncpus. You can change it manually.

__Output:__ ./out_vina folder which contains out docked pdbqt files.
Out files are saved as ./out_vina/*target-fname*_*ligand-fname*.pdbqt

__Run__: 
```
qsub ~/docking/run-vina.pbs -v path=$(pwd)
```

## pdb2mol.py
__Dependency__  
* RDKit (required)
* ChemAxon `cxcalc` (for `--protonation chemaxon` mode)
* Easydock with apptainer support + molgpka/unipka `.sif` image (for `--protonation /path/to/image.sif`)   
https://github.com/ci-lab-cz/easydock
- unipka SIF image can be obtained from https://zenodo.org/records/17854824.
- molgpka SIF image will be available on zenodo soon   

__Input:__ one or more *.pdb files and a reference SMILES string or *.smi file (names in the .smi must match PDB basenames)  
__Output:__ *.mol files with restored bond orders (hydrogens and coordinates are preserved); skipped molecules are listed in `mol2pdb_error.smi`

__Run examples__: 
```
# Protonate reference via ChemAxon (you will require ChemAxon licence) and convert PDB to MOL
python scripts/pdb2mol.py -i ligand.pdb -s ligands.smi --protonation chemaxon

# Use molgpka/unipka SIF protonation; download .sif from Zenodo
python scripts/pdb2mol.py -i ligand.pdb -s ligands.smi --protonation /path/to/unipka.sif

# Do not externally protonate reference and convert PDB to MOL with original protonation state   
python scripts/pdb2mol.py -i ligand.pdb -s ligands.smi

```
