import argparse
import datetime
import logging
import os
from multiprocessing import Pool

from rdkit import Chem
from requests import get

url_pdb_download = 'https://files.rcsb.org/download/{0}.pdb'
url_lig_download = 'https://models.rcsb.org/v1/{pdbid}/ligand?auth_asym_id={chain}&auth_seq_id={resNum}&encoding={format}'
url_fasta_download = 'https://www.rcsb.org/fasta/entry/{0}/download'

#pdb_get_inform_pdb = '''http://www.rcsb.org/pdb/rest/customReport.csv?pdbids={}&customReportColumns=structureTitle,experimentalTechnique,releaseDate,resolution,classification,macromoleculeType,residueCount,atomSiteCount,pdbDoi,ligandId,ligandName,ligandMolecularWeight,ligandSmiles,Ki,Kd,EC50,IC50,deltaG,compound,source,expressionHost,crystallizationMethod,crystallizationTempK,phValue,densityPercentSol,pdbxDetails,rFree,refinementResolution,publicationYear,title,journalName,doi&service=wsdisplay&format=json&ssa=null&primaryOnly=1&service=wsdisplay&format=csv'''
#pdb_get_inform_lig = 'http://rest.rcsb.org/rest/ligands/{0}'


def get_sdffrompdb(pdbID, ligID, chain, resNum, format):
    res = get(url_lig_download.format(pdbid=pdbID.lower(), chain=chain,
                                      resNum=resNum, format=format))
    if res.ok:
        return res.text
    else:
        print('Error. SDF Ligand Download.: {0} {1} {2} {3}.\n***{4}***'.format(pdbID, ligID,
                                                                                chain, resNum, res))
        return None


def get_pdb(pdbId):
    res = get(url_pdb_download.format(pdbId.upper()))
    if res.ok:
        return {pdbId: res.text}
    else:
        print('Error. PDB Download. PDBID: {0}.\n***{1}***'.format(pdbId, res))
        return None

def get_fasta(pdbId):
    res = get(url_fasta_download.format(pdbId.upper()))
    if res.ok:
        return {pdbId: res.text}
    else:
        print('Error. PDB FASTA Download. PDBID: {0}.\n***{1}***'.format(pdbId, res))
        return None


def get_ligs_field(pdb_data):
    lig_info = [line.strip() for line in pdb_data.split('\n') if line.startswith('HET ')]
    return lig_info


def get_clean_protein(pdb_data):
    clean_protein_block = [line for line in pdb_data.split('\n') if not line.startswith('HETATM')]
    return '\n'.join(clean_protein_block)


def get_water_only(pdb_data):
    water = [line for line in pdb_data.split('\n') if line.startswith('HETATM') and line[17:20] == 'HOH']
    return '\n'.join(water)


def main(pdbid_list, out_path, split, ncpu):
    p = Pool(ncpu)
    pdb_all_structrs = {}
    for i in p.map(get_pdb, pdbid_list):
        if i is not None:
            pdb_all_structrs.update(i)
    p.close()
    p.join()

    # save pdb
    all_smiles = {}
    for pdbid, pdbblock in pdb_all_structrs.items():
        outdir = os.path.join(out_path, pdbid)
        os.makedirs(outdir, exist_ok=True)

        with open(os.path.join(outdir, pdbid + '.pdb'), 'w') as out:
            out.write(pdbblock)
        if split:
            clean_protein = get_clean_protein(pdbblock)
            water = get_water_only(pdbblock)
            with open(os.path.join(outdir, pdbid + '_clean.pdb'), 'w') as out:
                out.write(clean_protein)
            with open(os.path.join(outdir, pdbid + '_water.pdb'), 'w') as out:
                out.write(water)

    # get fasta
    p = Pool(ncpu)
    fasta_all_dict = {}
    for i in p.map(get_fasta, pdbid_list):
        if i is not None:
            fasta_all_dict.update(i)
    p.close()
    p.join()

    for pdbid, fasta in fasta_all_dict.items():
        outdir = os.path.join(out_path, pdbid)
        with open(os.path.join(outdir, pdbid + '.fasta'), 'w') as out:
            out.write(fasta)

    # get ligands_info
    ligands_info_dict = {}

    for pdbid, struct in pdb_all_structrs.items():
        # HET field in pdb file
        lig_info = get_ligs_field(struct)
        lig_info_dict = {}
        for i in lig_info:
            ligid, chain, resnum = i[7:10].strip(), i[12], i[13:22].strip()
            if ligid in lig_info_dict:
                if chain in lig_info_dict[ligid]:
                    lig_info_dict[ligid][chain].append(resnum)
                else:
                    lig_info_dict[ligid].update({chain: [resnum]})
            else:
                lig_info_dict[ligid] = {chain: [resnum]}

        outdir = os.path.join(out_path, pdbid)
        with open(os.path.join(outdir, 'ligands_list.log'), 'w') as output:
            output.write('\n'.join(lig_info))

        if lig_info_dict:
            ligands_info_dict[pdbid] = lig_info_dict

    for pdbid, ligs_dict in ligands_info_dict.items():
        sdf_out = Chem.SDWriter(os.path.join(out_path, pdbid, f'{pdbid}_ligands_frompdb.sdf'))
        smi_dict = {}
        for ligid, ligdata in ligs_dict.items():
            for chain, resnum_list in ligdata.items():
                for resnum in resnum_list:
                    lig_molblock = get_sdffrompdb(pdbID=pdbid, ligID=ligid, chain=chain, resNum=resnum,
                                                  format='mol')
                    if lig_molblock is not None:
                        mol = Chem.MolFromMolBlock(lig_molblock, removeHs=False)
                        if not mol:
                            logging.error(f'Error to transform mol: {pdbid} {ligid} {chain} {resnum}')
                            mol = Chem.MolFromMolBlock(lig_molblock, removeHs=False, sanitize=False)

                        if mol:
                            mol.SetProp('_Name', f'{pdbid}_{ligid}_{chain}_{resnum}')
                            sdf_out.write(mol)
                            if ligid not in smi_dict:
                                smi = Chem.MolToSmiles(mol)
                                smi_dict[ligid] = smi
                        else:
                            logging.error(f'Error to save mol: {pdbid} {ligid} {chain} {resnum}')


        sdf_out.close()
        with open(os.path.join(out_path, pdbid, f'{pdbid}_ligands_frompdb.smi'), 'w') as output:
            for ligid, smi in smi_dict.items():
                output.write(f'{smi}\t{ligid}\n')

        all_smiles[pdbid] = smi_dict

        print('Ready: ', pdbid)

    with open(os.path.join(out_path, f'ligands_frompdb.smi'), 'w') as output:
        for pdbcode, smi_dict_pdbid in all_smiles.items():
            for ligid, smi in smi_dict_pdbid.items():
                molid = f'{pdbcode}_{ligid}_ligand'
                output.write(f'{smi}\t{molid}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Conversion of input PDB file to MOL with RDKit. '
                                                 'Save hydrogens and input charges from pdb.')
    parser.add_argument('-i', '--input', metavar='pdbid', required=True, nargs='+',
                        help='pdbid(s)')
    parser.add_argument('-o', '--outdir', metavar='UniprotID', required=True,
                        help='dirname to save')
    parser.add_argument('-n', '--ncpu', metavar='int', default=1,
                        help='ncpu')
    parser.add_argument('--split', action='store_true', default=False,
                        help='save clean protein and water pdb files')


    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    out_time = f'{datetime.datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}'
    log_file = os.path.join(args.outdir, f'GetfromPDB_{out_time}.log')
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])
    logging.info(args)

    main(pdbid_list=args.input, out_path=args.outdir, split=args.split, ncpu=args.ncpu)
