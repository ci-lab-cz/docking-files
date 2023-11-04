#!/usr/bin/env python3

import argparse
import os
import re
import logging
import datetime
import subprocess

from rdkit import Chem
from rdkit.Chem import AllChem


def get_major_microspecies(smi_fname, h='7.4', tautomerize=False):
    print('Protonation is turned on')
    canon_smi_dict = {}
    output = os.path.join(os.path.dirname(smi_fname), f'{os.path.basename(smi_fname).split(".")[0]}'                                         f'_protonation{"_majortautomer" if tautomerize else " "}'
                                                       f'_pH{h}')
    if tautomerize:
        print('Tautomerization is turned on')
        cmd_run = ['cxcalc', '-S', '--ignore-error', 'majormicrospecies', '-H', h, "-M", smi_fname]
    else:
        cmd_run = ['cxcalc', '-S', '--ignore-error', 'majormicrospecies', '-H', h, smi_fname]
    print(cmd_run)
    with open(output+'.sdf', 'w') as file:
        subprocess.run(cmd_run, stdout=file, text=True)
    for mol in Chem.SDMolSupplier(output+'.sdf', sanitize=False):
        if mol:
            smi = mol.GetPropsAsDict().get('MAJORMS', None)
            mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else smi
            if smi is not None:
                try:
                    cansmi = Chem.CanonSmiles(smi)
                    canon_smi_dict[mol_name.lower()] = cansmi
                except:
                    logging.error(f'ERROR: {mol_name}, Reference smiles {smi} obtained after protonation '
                                     f'could not be read by RDKit. The molecule was skipped.\n')
    with open(output+'.smi', 'w') as file:
        file.write('\n'.join([f'{j}\t{i}' for i, j in canon_smi_dict.items()]))

    return canon_smi_dict

def get_smi_dict(input_smi):
    smis = dict()
    with open(input_smi) as f:
        for line in f:
            values = line.strip().split()
            if len(values) >= 2:
                smis[values[1].lower()] = values[0]
            else:
                logging.warning(
                    f'Line "{line}" in input smiles does not have two fields - SMILES and mol name. Skipped.\n')
    return smis

def convertpdb2mol(input_fnames, input_smi, regex, no_protonation, tautomerization):

    smis = None
    if (input_smi is not None) and (os.path.isfile(input_smi)) and (input_smi.lower().endswith('.smi') or input_smi.lower().endswith('.smiles')):
        if no_protonation:
            smis = get_smi_dict(input_smi)
        else:
            smis = get_major_microspecies(input_smi, tautomerize=tautomerization)

    error_smi = {}

    for in_fname in input_fnames:
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
        ref_smi = Chem.MolToSmiles(template)
        template = Chem.AddHs(template)

        mol_new = None

        try:
            mol_new = AllChem.AssignBondOrdersFromTemplate(template, mol)
            mol_new.SetProp('_Name', mol_name)
        except Exception as e:
            logging.error(f'Fail to convert. Your PDB and smiles have different protonation. Problem: {e}. Mol: {ref_smi}\t{in_fname}')

        if mol_new:
            Chem.MolToMolFile(mol_new, in_fname.replace('.pdb', '.mol'))
        else:
            error_smi[mol_name] = ref_smi

    if error_smi:
        with open('mol2pdb_error.smi', 'w') as out:
            out.write('\n'.join([f'{j}\t{i}' for i, j in error_smi.items()]))


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
    parser.add_argument('--tautomerization', action='store_true', default=False,
                        help='Enable tautomerization of molecules during protonation.')

    args = parser.parse_args()
    logging.info(args)
    convertpdb2mol(input_fnames=args.input, input_smi=args.smi, regex=args.regex, no_protonation=args.no_protonation, tautomerization=args.tautomerization)



if __name__ == '__main__':
    out_time = f'{datetime.datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}'
    log_file = f'pdb2mol_savechargesfrompdb_{out_time}.log'
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])
    main()
