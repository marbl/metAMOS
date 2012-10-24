# This script classifies query fragments using BLASTN

import sys, os
import fileinput

if len(sys.argv) != 4:
	print 'BLASTN v1.0 by Donovan Parks, Norm MacDonald, and Rob Beiko'
	print 'Usage: python BLASTN.py <blastn-path> <query-file> <results-file>'
	print ''
	print 'Required parameters:'
	print '  <blastn-path>   Full path to blastn.'
	print '  <query-file>    Multi-FASTA file containing query fragments to classify.'
	print '  <results-file>  File to write classification results to.'
	print ''
	print 'Typical usage:'
	print '  python BLASTN.py /path/to/blastn test.fasta blastn_results.txt'
	print ''
	exit()

blastnEXE = sys.argv[1]
queryFile = sys.argv[2]
resultsFile = sys.argv[3]

print 'Running BLASTN...'

os.system(blastnEXE + ' -query ' + queryFile + ' -db ./blast_data/BacteriaAndArchaeaGenomesDB -evalue 10 -outfmt 7 -task blastn -out ' + resultsFile)

print 'Done.'
