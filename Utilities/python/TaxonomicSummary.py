# Produce flat file summarizing taxonomic assignments

import sys
import fileinput

if len(sys.argv) != 4:
	print 'TaxonomicSummary v1.0 by Donovan Parks, Norm MacDonald, and Rob Beiko'
	print ''
	print 'Usage: python TaxonomicSummary.py <query-file> <results-file> <summary-file>'
	print ''
	print 'Required parameters:'
	print '  <query-file>    Multi-FASTA file containing query fragments to classify.'
	print '  <results-file>  File indicating taxonomic classification of each fragment.'
	print '  <summary-file>  File where taxonomic summary information to.'
	print ''
	print 'Typical usage:'
	print '  python TaxonomicSummary.py test.fasta nb_topModels.txt nb_taxonomicSummary.txt'
	print ''
	exit()

queryFile = sys.argv[1]
classificationFile = sys.argv[2]
resultsFile = sys.argv[3]

# check for valid classification file
fin = open(classificationFile)
header = fin.readline()
fin.close()

if 'Fragment Id	Taxonomic classification' not in header:
	print ''
	print 'Error: Count not process classification file.'
	print ''
	print 'Results files produced by nb-classify or BLASTN.py must'
	print 'first be run through Epsilon-NB (with epsilon = 0) or' 
	print 'LCA (with percentage = 0), respectively, to be compatible'
	print 'with this script.'
	print ''
	exit()

# get length of each query fragment
print 'Determining length of each query fragment...'

fragmentLens = {}
fragmentId = ''
totalBasePairs = 0.0
totalFragments = 0.0
for line in fileinput.input([queryFile]):
	if line[0] == '>':
		if fragmentId != '':
			fragmentLens[fragmentId] = fragmentLen
			totalBasePairs += fragmentLen
		fragmentId = line[1:line.find(' ')].strip()
		fragmentLen = 0
		totalFragments += 1.0
	else:
		fragmentLen += len(line.strip())
fragmentLens[fragmentId] = fragmentLen
totalBasePairs += fragmentLen

# get classifications at each rank
print 'Determining classification at each rank...'
classification = [{}, {}, {}, {}, {}, {}, {}, {}]

bHeader = True
for line in fileinput.input([classificationFile]):
	if bHeader == True:
		bHeader = False
		continue
		
	lineSplit = line.split('\t')
	
	fragmentId = lineSplit[0]
	taxonomy = lineSplit[1].split(';')
	
	for r in xrange(0, 8):
		category = taxonomy[r]
		counts = classification[r].get(category, [0, 0])
		counts[0] += 1
		counts[1] += fragmentLens[fragmentId]
		classification[r][category] = counts

# write out results
print 'Writing out taxonomic summary...'
fout = open(resultsFile, 'w')
fout.write('Category\tFragments\tPercentage of fragments\tBase pairs\tPercentage of base pairs\n\n')
ranks = ['DOMAIN', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES', 'STRAIN']
for r in xrange(0, 8):
	fout.write(ranks[r] + '\n')
	
	categories = classification[r].keys()
	categories.sort()
	for category in categories:
		counts = classification[r][category]
		fout.write(category + '\t' + str(counts[0]) + '\t' + str(counts[0] / totalFragments) + '\t' + str(counts[1]) + '\t' + str(counts[1] / totalBasePairs) + '\n')
	fout.write('\n')
fout.close()

print 'Done.'