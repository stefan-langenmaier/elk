#!/usr/bin/python

# DOS transformer v0.1 by Markus Meinert
# Usage: dos_transformer.py PDOS_Sxx_Ayyyy.OUT
# Or: dos_transformer.py TDOS.OUT

import sys

print "\n DOS transformer v0.1 \n"

# Read the file.
filename = sys.argv[1]
f = open(filename, 'r')
data = f.readlines()
f.close()

# Analyze the file.
# Number of lines.
nlines = len(data)
print " Number of lines: %i " % nlines

# Count blank lines to determine number of datasets.
ndatasets = 0
for line in data:
	if line.split() == []:
		ndatasets += 1
print " Number of datasets: %i " % ndatasets

# Number of datasets is either 1, 2, 16 or 32: [TDOS or (1 x s, 3 x p, 5 x d, 7 x f)] * nspins.
if ndatasets == 1:
	nspins = 1
elif ndatasets == 2:
	nspins = 2
elif ndatasets == 16:
	nspins = 1
elif ndatasets == 32:
	nspins = 2
else:
	print " Strange number of datasets. Exiting."
	sys.exit()
print " Number of spins: %i " % nspins

# Number of energies is:
nenergies = (nlines - ndatasets)/ndatasets
print " Number of energies: %i " % nenergies

# Collect the datasets into a list of lists with a double-loop over datasets and energies.
datasets = []
for i in range(0,ndatasets):
        currentset = []
        for j in range(i*nenergies + i, (i*nenergies + i) + nenergies):
                currentset.append(data[j].split()) # Split each line by empty spaces.
        datasets.append(currentset)

# Merge the datasets line-wise. Output format is modified and converted to eV.
output = ""
for i in range(0,nenergies):
	line = '%.6e\t' % (float(datasets[0][i][0]) * 27.211) # First column is energy.
	# Append the DOS values as columns.
	for j in range(0, ndatasets):
		line += '%.6e\t' % (float(datasets[j][i][1]) / 27.211)
	line += "\n"
	output += line

filename = filename + ".transformed"
f = open(filename, 'w')
f.write(output)
f.close()

print
print " Output filename: %s" % filename
print " Output unit is eV."
if nspins == 1 and ndatasets == 1:
	print " Output format is: E 1"
elif nspins == 1 and ndatasets == 16:
	print " Output format is: E s p p p d d d d d f f f f f f f"
elif nspins == 2 and ndatasets == 2:
	print " Output format is: E 1 2"
elif nspins == 2 and ndatasets == 32:
	print " Output format is: E (s p p p d d d d d f f f f f f f)_1 (s p p p d d d d d f f f f f f f)_2"
print
print " Done."
print
