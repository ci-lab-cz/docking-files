import argparse

from chimera import openModels
from chimera import replyobj
from chimera import runCommand as rc

def main(complex, out, resn):
    replyobj.status("Processing " + complex)
    rc("open " + complex)
    rc("split ligand")
    complex_model_list = openModels.list()
    if len(complex_model_list) > 2:
        print('Warning! {0}. More than 1 ligand found. Check which one do you need to use.'.format(complex))
        with open(complex.split('.pdb')[0]+'_split.log', 'w') as output_log:
            output_log.write('Warning! {0}. More than 1 ligand found. Check which one do you need to use.'.format(complex))

    if resn:
        ligand_model = None
        for n, model in enumerate(complex_model_list, 1):
            # hierarchy model0.1, model0.2_ligid1
            model_name = model.name.lower().split(' ')
            # ligand model
            if resn.lower() in model_name:
                ligand_model = n
            # protein model
            if len(model_name) == 1:
                if len(complex_model_list) > 2:
                    # there are protein, ligand and cofactor
                    out_curr = out + '_protein_MD.pdb'
                else:
                    out_curr = out + '_protein.pdb'
            else:
                # other ligand model
                out_curr = out + '_' + model.name.split(' ')[1] + '_ligand.pdb'
            rc("write format pdb #0.{n} {out}".format(n=n, out=out_curr))

        if ligand_model:
            if len(complex_model_list) > 2:
                # save protein-cofactor model
                rc("close #0.{n}".format(n=ligand_model))
                rc("combine #0 newchainids false")
                rc("write format pdb #1 {out}".format(out=out + '_protein.pdb'))
    else:
        for n, model in enumerate(complex_model_list, 1):
            out_curr = '_'
            if n == 1:
                out_curr = out + '_protein.pdb'
            if n > 1:
                out_curr = out + '_' + model.name.split(' ')[1] + '_ligand.pdb'
            rc("write format pdb #0.{n} {out}".format(n=n, out=out_curr))

    rc("stop now")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', metavar='complex.pdb', help='input pdb')
    parser.add_argument('--out', metavar='pdbid', help='PDBID prefix for output files')
    parser.add_argument('--resn', metavar='UNL', default=None)
    args = parser.parse_args()

    main(args.input, args.out, args.resn)
