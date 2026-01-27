#!/usr/bin/env python3

import argparse
import os
import re
import logging
import datetime
import subprocess
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem


def protonate_sif(smiles, sif_path):
    # protonate molecules with uni-pka
    try:
        from easydock.containers import apptainer_exec
        from easydock.auxiliary import expand_path
    except ImportError as exc:
        raise RuntimeError('Protonation via SIF requires easydock to be installed.') from exc

    sif_path = expand_path(sif_path)
    if not os.path.isfile(sif_path):
        raise FileNotFoundError(f'Provided SIF file does not exist: {sif_path}')

    nmols = len(smiles)
    protonated_mols = dict()
    with tempfile.NamedTemporaryFile(suffix='.smi', mode='w', encoding='utf-8') as tmp:
        for i, smi in enumerate(smiles):
            tmp.write(f'{smi}\t{i}\n')
        tmp.flush()

        tmpdir = tempfile.mkdtemp()
        output = os.path.join(tmpdir, "output.smi")
        try:
            bind_path = set()
            bind_path.add(os.path.dirname(expand_path(tmp.name)))
            bind_path.add(os.path.dirname(expand_path(output)))
            apptainer_exec(sif_path,
                           ['protonate', '-i', tmp.name, '-o', output],
                           list(bind_path))
            # protonate_apptainer(tmp.name, output, '/home/polishcp/imgs/unipka.sif')
            with open(output) as f:
                for line in f:
                    tmp = line.strip().split()
                    mol_prot = Chem.MolFromSmiles(tmp[0])
                    if mol_prot:
                        mol_prot.SetProp('_Name', tmp[1])
                        protonated_mols[int(tmp[1])] = mol_prot
        finally:
            if os.path.exists(output):
                os.remove(output)
            os.rmdir(tmpdir)
    # keep the order of input mols and replace with None missing mols
    output = [protonated_mols.get(i, None) for i in range(nmols)]
    return output


def protonate_smi_file_with_sif(smi_fname, sif_path):
    """Protonate SMILES file using SIF and return dict name->canonical SMILES."""
    smiles_list = []
    names = []
    with open(smi_fname) as f:
        for line in f:
            values = line.strip().split()
            if len(values) >= 2:
                smiles_list.append(values[0])
                names.append(values[1])
            else:
                logging.warning(
                    f'Line "{line}" in input smiles does not have two fields - SMILES and mol name. Skipped.\n')

    if not smiles_list:
        return {}

    protonated_mols = protonate_sif(smiles_list, sif_path)
    canon_smi_dict = {}
    for name, mol in zip(names, protonated_mols):
        if mol:
            try:
                canon_smi_dict[name.lower()] = Chem.CanonSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))
            except Exception:
                logging.error(f'ERROR: {name}, Reference smiles obtained after protonation '
                              f'could not be read by RDKit. The molecule was skipped.\n')
    return canon_smi_dict

def get_major_microspecies(smi_fname, h='7.4', tautomerize=False):
    print('Protonation is turned on')
    canon_smi_dict = {}
    output = os.path.join(os.path.dirname(smi_fname), f'{os.path.basename(smi_fname).split(".")[0]}'                                     
                          f'_protonation{"_majortautomer" if tautomerize else ""}'
                                                       f'_pH{h}')
    if tautomerize:
        print('Tautomerization is turned on')
        cmd_run = ['cxcalc', '-S', '--ignore-error', 'majormicrospecies', '-H', h, "-M", smi_fname]
    else:
        cmd_run = ['cxcalc', '-S', '--ignore-error', 'majormicrospecies', '-H', h, smi_fname]
    print(cmd_run)
    with open(output+'.sdf', 'w') as file:
        res = subprocess.run(cmd_run, stdout=file, text=True, check=False)
    if res.returncode != 0:
        logging.error("cxcalc failed: %s", res.stderr.strip())
        raise RuntimeError("cxcalc protonation failed")
    
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


def resolve_protonation(protonation_arg):
    """Resolve user input to protonation mode and sif path."""
    if protonation_arg:
        arg_lower = protonation_arg.lower()
        if arg_lower == 'chemaxon':
            return 'chemaxon', None
        if arg_lower == 'none':
            return 'none', None
        if protonation_arg.endswith('.sif') or os.path.isfile(protonation_arg):
            return 'sif', protonation_arg
        raise ValueError('Protonation must be "chemaxon", "none", or a path to an existing .sif file.')

    # default: no protonation
    return 'none', None


def convertpdb2mol(input_fnames, input_smi, regex, protonation_mode, tautomerization, protonation_path=None,
                   add_hs_reference=True):

    smis = None
    if tautomerization and protonation_mode != 'chemaxon':
        logging.warning('Tautomerization flag is ignored when protonation is not chemaxon.')

    if (input_smi is not None) and (os.path.isfile(input_smi)) and (input_smi.lower().endswith('.smi') or input_smi.lower().endswith('.smiles')):
        if protonation_mode == 'none':
            smis = get_smi_dict(input_smi)
        elif protonation_mode == 'sif':
            if not protonation_path:
                raise ValueError('Protonation mode "sif" selected but no SIF file was provided.')
            if not os.path.isfile(protonation_path):
                raise FileNotFoundError(f'Provided SIF file does not exist: {protonation_path}')
            smis = protonate_smi_file_with_sif(input_smi, protonation_path)
        else:  # chemaxon
            smis = get_major_microspecies(input_smi, tautomerize=tautomerization)
    elif protonation_mode == 'sif':
        if not protonation_path:
            raise ValueError('Protonation mode "sif" selected but no SIF file was provided.')
        if not os.path.isfile(protonation_path):
            raise FileNotFoundError(f'Provided SIF file does not exist: {protonation_path}')
        protonated = protonate_sif([input_smi], protonation_path)
        if protonated and protonated[0]:
            try:
                input_smi = Chem.MolToSmiles(protonated[0], isomericSmiles=True)
            except Exception:
                logging.error('ERROR: reference smiles obtained after protonation could not be read by RDKit. '
                              'Original SMILES will be used.\n')

    error_smi = {}

    for in_fname in input_fnames:
        base_name = os.path.splitext(os.path.basename(in_fname))[0]
        mol_name = None
        if regex:
            mol_name = re.search(regex, base_name)
            if mol_name:
                mol_name = mol_name.group()
        if not mol_name:
            mol_name = base_name

        if smis:
            if mol_name.lower() in smis:
                smi = smis[mol_name.lower()]
            else:
                logging.warning(f'Molecule {in_fname}: {mol_name}-smiles pair was not found. Molecule will be skipped')
                continue
        else:
            smi = input_smi

        mol = Chem.MolFromPDBFile(in_fname, sanitize=False, removeHs=False)
        reference = Chem.MolFromSmiles(smi)
        ref_smi = Chem.MolToSmiles(reference)
        if add_hs_reference:
            reference = Chem.AddHs(reference)

        mol_new = None

        try:
            mol_new = AllChem.AssignBondOrdersFromTemplate(reference, mol)
            mol_new.SetProp('_Name', mol_name)
        except Exception as e:
            logging.error(f'Fail to convert. Your PDB and smiles have different protonation. Problem: {e}. Mol: {ref_smi}\t{in_fname}')
        try:
            Chem.SanitizeMol(mol_new)
        except Exception as e:
            logging.warning(f'SanitizeMol failed for {mol_name}: {e}. Continuing with unsanitized molecule.')

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
    parser.add_argument('--protonation', default=None,
                        help='Protonation selection: "chemaxon", or path to a .sif image. '
                             'If omitted, no protonation is performed.')
    parser.add_argument('--tautomerization', action='store_true', default=False,
                        help='Enable tautomerization of molecules during protonation.')
    parser.add_argument('--noaddHs_reference', action='store_true', default=False,
                        help='Do not add hydrogens to the reference molecule before assigning bond orders.')

    args = parser.parse_args()
    logging.info(args)
    protonation_mode, sif_path = resolve_protonation(args.protonation)
    convertpdb2mol(input_fnames=args.input, input_smi=args.smi, regex=args.regex,
                   protonation_mode=protonation_mode, tautomerization=args.tautomerization,
                   protonation_path=sif_path, add_hs_reference=not args.noaddHs_reference)


if __name__ == '__main__':
    out_time = f'{datetime.datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}'
    log_file = f'pdb2mol_savechargesfrompdb_{out_time}.log'
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])
    main()
