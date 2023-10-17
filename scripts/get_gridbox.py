#!/usr/bin/env python3

import argparse
import numpy as np
import os
from glob import glob

from rdkit import Chem


def get_ligbox_coords(pdb_fname):
    mol = Chem.MolFromPDBFile(pdb_fname, removeHs=True, sanitize=False)
    xyz = mol.GetConformer().GetPositions()
    min_values = np.min(xyz, axis=0)
    max_values = np.max(xyz, axis=0)
    return min_values, max_values


def main():

    parser = argparse.ArgumentParser(description='Get coordinates for a gridbox based on a position of a ligand(s).')
    parser.add_argument('-i', '--input', metavar='FILENAME(S) or DIRNAME', required=True, nargs='+',
                        help='PDB file(s) with ligand(s) or a directory name which will be searched for all entries '
                             'ending with _ligand.pdb. If multiple ligands are supplied a single gridbox will be '
                             'generated to fit all ligands.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='text file where the center and the size of a gridbox will be stored.')
    parser.add_argument('-d', '--dist',  metavar='NUMERIC', required=False, default=5, type=float,
                        help='the distance which will be added to minimum and maximum values of ligand coordinates. '
                             'Default: 5.')

    args = parser.parse_args()
    values = []   # (min_values, max_values)
    for name in args.input:
        if os.path.isdir(name):
            for fname in glob(name + '/**/*_ligand.pdb', recursive=True):
                values.append(get_ligbox_coords(fname))
        if os.path.isfile(name):
            values.append(get_ligbox_coords(name))
    if values:
        min_values = np.array([values[0][0]])  # reshape 1d to 2d array with ine row
        max_values = np.array([values[0][1]])  # reshape 1d to 2d array with ine row
        for v in values[1:]:
            min_values = np.vstack((min_values, v[0]))
            max_values = np.vstack((max_values, v[1]))
        min_values = np.min(min_values, axis=0)
        max_values = np.min(max_values, axis=0)
        center = np.round((min_values + max_values) / 2, 2)
        size = np.round(max_values - min_values + 2 * args.dist, 2)
        with open(args.output, 'wt') as f:
            f.write(f'size_x = {size[0]}\n')
            f.write(f'size_y = {size[1]}\n')
            f.write(f'size_z = {size[2]}\n')
            f.write(f'center_x = {center[0]}\n')
            f.write(f'center_y = {center[1]}\n')
            f.write(f'center_z = {center[2]}\n')


if __name__ == '__main__':
    main()
