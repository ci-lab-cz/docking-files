#!/usr/bin/env python3

import argparse
import sys
import os
import re
import logging
import datetime
import tempfile
import subprocess

from rdkit import Chem
from rdkit.Chem import *

def fix_valence_charge(m, mol_name):
    m.UpdatePropertyCache(strict=False)
    ps = Chem.DetectChemistryProblems(m)
    for p in ps:
        if p.GetType() != 'KekulizeException':
            at = m.GetAtomWithIdx(p.GetAtomIdx())

            print(f'Warning!!!!!!! Your reference tautomeric form doesnot match with your pdb. {mol_name}. '
                  f'Atom - ID:{at.GetIdx()} Atomic Num:{at.GetAtomicNum()}. {p.GetType()}')

            if at.GetAtomicNum() == 7 and at.GetTotalValence() == 4 and at.GetFormalCharge() != 1:
                at.SetFormalCharge(1)

            if at.GetAtomicNum() == 5 and at.GetTotalValence() == 4 and at.GetFormalCharge() != -1:
                at.SetFormalCharge(-1)
    return m

def update_charge(mol, mol_name):
    for a in mol.GetAtoms():
        if a.GetNumImplicitHs() != 0 :
            print(f'Warning!!!!!!! Your reference tautomeric form doesnot match with your pdb. {mol_name}. '
                  f'Atom - ID:{a.GetIdx()} Atomic Num:{a.GetAtomicNum()}')
            a.SetFormalCharge(a.GetFormalCharge()-a.GetNumImplicitHs())
            a.SetNoImplicit(True)

    return mol

def AssignBondOrdersFromTemplate_fix(refmol, mol, mol_name):
  refmol2 = rdchem.Mol(refmol)
  mol2 = rdchem.Mol(mol)
  # do the molecules match already?
  matching = mol2.GetSubstructMatch(refmol2)
  if not matching:  # no, they don't match
    # check if bonds of mol are SINGLE
    for b in mol2.GetBonds():
      if b.GetBondType() != BondType.SINGLE:
        b.SetBondType(BondType.SINGLE)
        b.SetIsAromatic(False)
    # set the bonds of mol to SINGLE
    for b in refmol2.GetBonds():
      b.SetBondType(BondType.SINGLE)
      b.SetIsAromatic(False)
    # set atom charges to zero;
    for a in refmol2.GetAtoms():
      a.SetFormalCharge(0)
    for a in mol2.GetAtoms():
      a.SetFormalCharge(0)

    matching = mol2.GetSubstructMatches(refmol2, uniquify=False)
    # do the molecules match now?
    if matching:
      if len(matching) > 1:
        print("More than one matching pattern found - picking one")
      matching = matching[0]
      # apply matching: set bond properties
      for b in refmol.GetBonds():
        atom1 = matching[b.GetBeginAtomIdx()]
        atom2 = matching[b.GetEndAtomIdx()]
        b2 = mol2.GetBondBetweenAtoms(atom1, atom2)
        b2.SetBondType(b.GetBondType())
        b2.SetIsAromatic(b.GetIsAromatic())
      # apply matching: set atom properties
      for a in refmol.GetAtoms():
        a2 = mol2.GetAtomWithIdx(matching[a.GetIdx()])
        a2.SetHybridization(a.GetHybridization())
        a2.SetIsAromatic(a.GetIsAromatic())
        a2.SetNumExplicitHs(0)
        a2.SetFormalCharge(a.GetFormalCharge())

      mol2 = Chem.RemoveHs(mol2, implicitOnly=True)
      Chem.KekulizeIfPossible(mol2)
      mol2 = fix_valence_charge(mol2, mol_name)

      mol2 = update_charge(mol2, mol_name)

      status = Chem.SanitizeMol(mol2, catchErrors=True)
      if status is not Chem.rdmolops.SanitizeFlags.SANITIZE_NONE:
          mol2.UpdatePropertyCache(strict=False)
          ps = Chem.DetectChemistryProblems(mol2)
          return mol2, ps
    else:
      return None, 'no matching found'
  return mol2, []


def get_major_microspecies(smi_fname, h='7.4', tautomerize=True):
    canon_smi_dict = {}
    output = os.path.join(os.path.dirname(smi_fname), f'{os.path.basename(smi_fname).split(".")[0]}'
                                                       f'_protonation{"_majortautomer" if tautomerize else ""}'
                                                       f'_pH{h}.smi')
    cmd_run = ['cxcalc', '-S', '--ignore-error', 'majormicrospecies', '-H', h,
               f'{"-M" if tautomerize else ""}', smi_fname]
    with open(output, 'w') as file:
        subprocess.run(cmd_run, stdout=file, text=True)
    for mol in Chem.SDMolSupplier(output, sanitize=False):
        if mol:
            smi = mol.GetPropsAsDict().get('MAJORMS', None)
            mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else smi
            if smi is not None:
                try:
                    cansmi = Chem.CanonSmiles(smi)
                    canon_smi_dict[mol_name.lower()] = cansmi
                except:
                    sys.stderr.write(f'ERROR: {mol_name}, Reference smiles {smi} obtained after protonation '
                                     f'could not be read by RDKit. The molecule was skipped.\n')

    return canon_smi_dict

def convertpdb2mol(input_fnames, input_smi, regex, no_protonation, no_tautomerization):

    smis = dict()
    if (input_smi is not None) and (os.path.isfile(input_smi)) and (input_smi.lower().endswith('.smi') or input_smi.lower().endswith('.smiles')):
        if no_protonation:
            with open(input_smi) as f:
                for line in f:
                    values = line.strip().split()
                    if len(values) >= 2:
                        smis[values[1].lower()] = values[0]
                    else:
                        sys.stderr.write(
                            f'Line "{line}" in input smiles does not have two fields - SMILES and mol name. Skipped.\n')
        else:
            smis = get_major_microspecies(input_smi, tautomerize=not no_tautomerization)


    for in_fname in input_fnames:
        print(in_fname)
        mol_name = None
        if regex:
            mol_name = re.search(regex, os.path.basename(in_fname).replace('.pdb',''))
            if mol_name:
                mol_name = mol_name.group()
        if not mol_name:
            mol_name = os.path.basename(in_fname).replace('.pdb','')

        if smis:
            if mol_name.lower() in smis:
                smi = smis[mol_name.lower()]
            else:
                logging.warning(f'Molecule {in_fname}: {mol_name}-smiles pair was not found. Molecule will be skipped')
                continue
        else:
            smi = input_smi

        mol = Chem.MolFromPDBFile(in_fname, sanitize=False, removeHs=False)
        template = Chem.MolFromSmiles(smi)
        mol_new, ps = AssignBondOrdersFromTemplate_fix(template, mol, mol_name)
        if not ps:
            # to make all hydrogens explicit
            mol_new = Chem.AddHs(mol_new, addCoords=True)
            mol_new.SetProp('_Name', mol_name)
            Chem.MolToMolFile(mol_new, in_fname.replace('.pdb', '.mol'))
            print(in_fname+' finished successfully')
        else:
            if isinstance(ps,tuple):
                logging.warning(f'Molecule {in_fname} fails to be transformed. {"; ".join([i.Message() for i in ps])}')
            else:
                logging.warning(f'Molecule {in_fname} fails to be transformed. {ps}')



def main():

    parser = argparse.ArgumentParser(description='Conversion of input PDB file to MOL with RDKit. '
                                                 'Save hydrogens and input charges from pdb.')
    parser.add_argument('-i', '--input', metavar='input.pdb', required=True, nargs='+',
                        help='input PDB file.')
    parser.add_argument('-s', '--smi', metavar='STRING or FILE', required=True,
                        help='SMILES or File of SMILES of a molecule in PDB to restore bond orders.')
    parser.add_argument('--regex', metavar='REGEX', required=False, default=None,
                        help='Use it if there are complex names of pdbqt files. '
                             'Use regex search to establish a relationship between reference smiles name and pdbqt '
                             'filename. If None filename of pdbqt file will be taken to find reference smiles '
                             'Examples: MOLID[0-9]* or .*(?=_out.pdbqt)')
    parser.add_argument('--no_protonation', action='store_true', default=False,
                        help='disable protonation of molecules before docking. Protonation requires installed '
                             'cxcalc chemaxon utility.')
    parser.add_argument('--no_tautomerization', action='store_true', default=False,
                        help='disable tautomerization of molecules during protonation.')

    args = parser.parse_args()
    convertpdb2mol(input_fnames=args.input, input_smi=args.smi, regex=args.regex, no_protonation=args.no_protonation, no_tautomerization=args.no_tautomerization)



if __name__ == '__main__':
    out_time = f'{datetime.datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}'
    log_file = f'pdb2mol_savechargesfrompdb_{out_time}.log'
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])
    main()
