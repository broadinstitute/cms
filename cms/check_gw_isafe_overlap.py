##	can I stitch together values like this? how well do overlap SNP scores correlate

basedir = "/idi/sabeti-scratch/jvitti/iSAFE/1kg_gw/"
filename1, filename2 = "chr22_CEU_22_31050115_36050115.iSAFE.out", "chr22_CEU_22_33550115_38550115.iSAFE.out"

from scipy.stats import pearsonr

def read_isafe_file(infilename):
	positions, values = [], []
	openfile = open(infilename, 'r')
	openfile.readline()
	for line in openfile:
		entries = line.split()
		thisPos, thisVal = int(entries[0]), float(entries[1])
		positions.append(thisPos)
		values.append(thisVal)
	openfile.close()
	return positions, values

def main():
	pos1, val1 = read_isafe_file(basedir + filename1)
	pos2, val2 = read_isafe_file(basedir + filename2)
	overlap_pos = [item for item in pos1 if item in pos2]
	overlap_vals1 = [val1[pos1.index(item)] for item in overlap_pos]
	overlap_vals2 = [val2[pos2.index(item)] for item in overlap_pos]
	r2 = pearsonr(overlap_vals1, overlap_vals2)
	print(r2)

main()

