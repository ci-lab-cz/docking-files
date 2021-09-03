## convert2pdb2pdbqt.pbs
#### script for the preparation of ligands
__Run__: 
```
qsub ./convert2pdb2pdbqt.pbs -v path=$(pwd)
```
__Input:__ ligand.smi or ligand.sdf  
To keep things simple, the name should be *ligand* or, if you don't want it, you can make small changes to the script files.  

__Output:__ new folder (folder name is equal to smi/sdf file name) which contains ligand pdbqt files  
## run-vina.pbs
__Run__: 
```
qsub ./run-vina_1.pbs  -v path=$(pwd)
```
__Input:__ 
* __./target__ folder which contains *pdbqts of targets* and *config.txt*
* __./ligand__ folder which contains *pdbqts of ligands*  
