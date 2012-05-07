# Classify fragments using LCA.

import sys
import fileinput
import math

if len(sys.argv) != 5:
	print 'LCA v1.0.2 by Donovan Parks, Norm MacDonald, and Rob Beiko'
	print ''
	print 'Usage: python LCA.py <blastn-results> <E-value> <percentage> <results-file>'
	print ''
	print 'Required parameters:'
	print '  <blastn-results>  Results of BLASTN classifier.'
	print '  <E-value>         Ignore hits with an E-value above this threshold.'
	print '  <percentage>      Use all hits with a bit score with this percentage of the '
	print '                      hit with the highest bit score to classify a fragment.'
	print '  <results-file>    File to write classification results to.'
	print ''
	print 'Typical usage:'
	print '  python LCA.py blastn_results.txt 1E-5 15 lca_results.txt'
	print ''
	exit()
	
blastnResults = sys.argv[1]
eValueThreshold = float(sys.argv[2])
percentageThreshold = float(sys.argv[3]) / 100.0
resultsFile = sys.argv[4]

# read in complete taxonomy of each accession number
accessionToTaxonomy = {}
for line in  fileinput.input(['taxonomy.txt']):
	lineSplit = line.split('\t')
	
	accession = lineSplit[0]
	taxonomy = lineSplit[1].split(';')
	
	accessionToTaxonomy[accession] = taxonomy
		
# determine all hits above the specified E-value threshold for each fragment
print 'Determing hits for each query fragment...'
fragmentHits = {}
for line in  fileinput.input([blastnResults]):
	if '# Query:' in line:
		fragmentId = line[line.rfind(' '):]
		fragmentHits[fragmentId.strip()] = [[0, ['unclassified']*8]]
		continue
		
	if line[0] == '#':
		continue
		
	lineSplit = line.split('\t')
	fragmentId = lineSplit[0]
	accession = lineSplit[1]
	evalue = float(lineSplit[10])
	bitscore = float(lineSplit[11])
	
	if evalue > eValueThreshold:
		continue
	
	fragmentHits[fragmentId].append([bitscore, list(accessionToTaxonomy[accession])])

# Determine consensus classification amongst all hits
print 'Classifying fragments with LCA...'
topTaxonomies = {}
for fragmentId in fragmentHits:
	hits = fragmentHits[fragmentId]
	
	# determine top bitscore
	topBitScore = -1
	for bitScore, taxonomy in hits:
		if bitScore > topBitScore:
			topTaxonomy = taxonomy
			topBitScore = bitScore
			
	# determine concensus classification
	for bitScore, taxonomy in hits:
		if bitScore >= (1.0 - percentageThreshold)*topBitScore:
			for r in xrange(0, 8):
				if topTaxonomy[r] != taxonomy[r]:
					topTaxonomy[r] = 'unclassified'
				
	topTaxonomies[fragmentId] = topTaxonomy

# write out classification results
print 'Writing out results...'
fout = open(resultsFile, 'w')
fout.write('Fragment Id\tTaxonomic classification\n')
for fragmentId in topTaxonomies:
	fout.write(fragmentId + '\t')
	
	topTaxonomy = topTaxonomies[fragmentId]
	for r in xrange(0, 8):
		fout.write(topTaxonomy[r] + ';')
	fout.write('\n')
		
fout.close()

print 'Done.'