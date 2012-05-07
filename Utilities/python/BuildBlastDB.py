# This script build a BLAST database for all
# completely sequenced genomes.

import sys, os
import fileinput
import platform

if len(sys.argv) != 3:
	print 'BuildBlastDb v1.0 by Donovan Parks, Norm MacDonald, and Rob Beiko'
	print ''
	print 'Usage: python BuildBlastDb.py <makeblastdb-path> <sequence-file>'
	print ''
	print 'Required parameters:'
	print '  <makeblastdb-path>  Full path to makeblastdb.'
	print '  <sequence-file>     File indicating all sequences to contain in database'
	print ''
	print 'Typical usage: '
	print '  python BuildBlastDb.py /path/to/executable/makeblastdb ./training/sequences.txt'
	print ''
	exit()
	
dbExeFile = sys.argv[1]
seqFile = sys.argv[2]

# concatenate all sequences into a single file
if not os.path.exists('./blast_data'):
	os.makedirs('./blast_data')

fout = open('./blast_data/_temp_.txt', 'w')
	
print 'Concatenating genome files...'
for seq in fileinput.input([seqFile]): 
	fin = open(seq.strip())
	data = fin.readlines()
	fin.close()
	
	for line in data:
		fout.write(line)
fout.close()

# create blast database
print 'Creating BLAST database...'

os.chdir('./blast_data')
os.system(dbExeFile + ' -max_file_sz 100GB -in _temp_.txt -dbtype nucl -title BacteriaAndArchaeaGenomes -out BacteriaAndArchaeaGenomesDB')

if platform.system() == 'Windows':
	os.system('del _temp_.txt')
else:
	os.system('rm _temp_.txt')
	
print 'Done.'