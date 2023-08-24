from chimera import runCommand as rc
from chimera import replyobj
import argparse

def main(fname_list):
	for fn in fname_list:
		pdb = fn.replace('.pdb','')+'_HIS.pdb'
		replyobj.status("Processing " + fn) # show what file we're working on
		rc("open " + fn)
		rc("setattr r type HIS :HIZ")
		rc("setattr r type HID :HIS@HD1,DD1,TD1,HND,H01")
		rc("setattr r type HIP :HID@HE2,DE2,TE2")
		rc("setattr r type HIE :HIS@HNE,HE2,H02")
		rc("write format pdb #0 %s" % pdb)

	rc("stop now")

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',  required=False, nargs='*', metavar='fname.pdb', help='input pdb')
	args = parser.parse_args()

	main(fname_list=args.input)
