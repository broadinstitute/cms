
olderfile = "/idi/sabeti-scratch/jvitti/1kg_components/popihh/unmasked_aa_v2/chr4_CEU_ihh.txt.gz"
newerfile = "/idi/sabeti-scratch/jvitti/1kg_components/popihh/unmasked/chr4_CEU_ihh.txt.gz"#truncated

import gzip
from scipy.stats.stats import pearsonr

def read_file(infilename):
	positions, values = [], []
	openfile = gzip.open(infilename, 'rb')
	for line in openfile:
		entries = line.split()
		positions.append(int(entries[0]))
		values.append(float(entries[-1]))
	openfile.close()
	return positions, values 

def main():
	old_pos, old_vals = read_file(olderfile)
	new_pos, new_vals = read_file(newerfile)	

	shared = [item for item in old_pos if item in new_pos]

	shared_old_vals = [old_vals[old_pos.index(item)] for item in shared]
	shared_new_vals = [new_vals[new_pos.index(item)] for item in shared]

	infostring = str(pearsonr(shared_old_vals, shared_new_vals))
	print(infostring)

main()
