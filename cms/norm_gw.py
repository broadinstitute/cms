import numpy as np
import sys
import os
selpop, model, likessuffix = sys.argv[1], sys.argv[2], sys.argv[3]
if len(sys.argv) != 4:
	print('call program as: python norm_gw.py <selpop> <model> <likessuffix>')
	sys.exit()

chroms = range(1,23)


def normalize(rawscore, mean, sd):
	rawscore, mean, sd = float(rawscore), float(mean), float(sd)
	normalizedvalue = (rawscore - mean) / sd
	return normalizedvalue

scores = []
for chrom in chroms:
	unnormedfile = "/idi/sabeti-scratch/jvitti/synth/cms_composite/"+ selpop  +"_" + str(model) + "_" + likessuffix + ".chr" + str(chrom) + ".txt"
	assert os.path.isfile(unnormedfile)
	openfile = open(unnormedfile, 'r')
	for line in openfile:
		entries=line.split()
		scores.append(float(entries[-1]))
	openfile.close()

print('loaded ' + str(len(scores)) + " scores")
print(str(max(scores)))
print(str(min(scores)))
print(str(np.mean(scores)))
print(str(np.var(scores)))
mean = np.mean(scores)
var = np.var(scores)
sd = np.sqrt(var)

for chrom in chroms:
	unnormedfile = "/idi/sabeti-scratch/jvitti/synth/cms_composite/"+ selpop  +"_" + str(model) + "_" + likessuffix + ".chr" + str(chrom) + ".txt"
	assert os.path.isfile(unnormedfile)
	normedfile = "/idi/sabeti-scratch/jvitti/synth/cms_composite/"+ selpop  +"_" + str(model) + "_" + likessuffix + ".chr" + str(chrom) + ".txt.norm"


	readfile = open(unnormedfile, 'r')
	writefile = open(normedfile, 'w')
	for line in readfile:
		line = line.strip('\n')
		entries=line.split()
		rawscore = float(entries[-1])
		normedscore = normalize(rawscore, mean, sd)
		writeline = line + "\t" + str(normedscore) + '\n'
		writefile.write(writeline)
	readfile.close()
	writefile.close
	print('wrote to '  + normedfile)
