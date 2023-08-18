#!/usr/bin/env python3

import argparse
import sys
import os
import re
import logging
import datetime

from rdkit import Chem
from rdkit.Chem import *

def fix_valence_charge(m):
    m.UpdatePropertyCache(strict=False)
    ps = Chem.DetectChemistryProblems(m)
    for p in ps:
        if p.GetType() != 'KekulizeException':
            at = m.GetAtomWithIdx(p.GetAtomIdx())
            if at.GetAtomicNum() == 7 and at.GetTotalValence() == 4 and at.GetFormalCharge() != 1:
                at.SetFormalCharge(1)

            if at.GetAtomicNum() == 5 and at.GetTotalValence() == 4 and at.GetFormalCharge() != -1:
                at.SetFormalCharge(-1)

            if at.GetNumImplicitHs() > 0:
                at.SetFormalCharge(at.GetFormalCharge() - at.GetNumImplicitHs())
                at.SetNoImplicit(True)
    return m


def AssignBondOrdersFromTemplate_fix(refmol, mol):
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
      mol2 = fix_valence_charge(mol2)

      status = Chem.SanitizeMol(mol2, catchErrors=True)
      if status is not Chem.rdmolops.SanitizeFlags.SANITIZE_NONE:
          mol2.UpdatePropertyCache(strict=False)
          ps = Chem.DetectChemistryProblems(mol2)
          return mol2, ps
    else:
      return None, 'no matching found'
  return mol2, []

def convertpdb2mol(input_fnames, input_smi, regex):

    smis = dict()
    if (input_smi is not None) and (os.path.isfile(input_smi)) and (input_smi.lower().endswith('.smi') or input_smi.lower().endswith('.smiles')):
        with open(input_smi) as f:
            for line in f:
                values = line.strip().split()
                if len(values) >= 2:
                    smis[values[1].lower()] = values[0]
                else:
                    sys.stderr.write(
                        f'Line "{line}" in input smiles does not have two fields - SMILES and mol name. Skipped.\n')

    for in_fname in input_fnames:
        print(in_fname)
        mol_name = None
        if regex:
            mol_name = re.search(regex, os.path.basename(in_fname))
            if mol_name:
                mol_name = mol_name.group()
        if not mol_name:
            mol_name = os.path.basename(in_fname)

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
        mol_new, ps = AssignBondOrdersFromTemplate_fix(template, mol)
        if not ps:
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

    args = parser.parse_args()
    convertpdb2mol(input_fnames=args.input, input_smi=args.smi, regex=args.regex)



if __name__ == '__main__':
    out_time = f'{datetime.datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}'
    log_file = f'pdb2mol_savechargesfrompdb_{out_time}.log'
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])
    main()
