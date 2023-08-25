import argparse
import shutil

from chimera import replyobj  # for emitting status messages
from chimera import runCommand as rc  # use 'rc' as shorthand for runCommand


def main(ref, mobile_list):
    replyobj.status("Processing " + ref)
    rc("open " + ref)
    for n, model in enumerate(mobile_list, 1):
        shutil.copy(model, model.split('.pdb')[0] + '_orig.pdb')
        replyobj.status("Processing " + model)
        rc("open " + model)
        if ref != model:
            rc("matchmaker #0 #1")

        rc("write format pdb #1 {model}".format(model=model))
        rc("close #1")

    rc("stop now")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref', metavar='ref.pdb', help='reference pdb model')
    parser.add_argument('--mobile', metavar='input_mobile.pdb',
                        nargs='+', help='mobile pdb model')
    args = parser.parse_args()
    print(args)
    main(args.ref, args.mobile)
