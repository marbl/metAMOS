# Classify fragments using LCA + Epsilon-NB.

import sys
import fileinput
import math

if len(sys.argv) != 7:
	print 'LCA+Epsilon-NB v1.0 by Donovan Parks, Norm MacDonald, and Rob Beiko'
	print ''
	print 'Usage: python LCA+Epsilon-NB.py <blastn-results> <nb-results> <E-value> <percentage> <epsilon> <results-file>'
	print ''
	print 'Required parameters:'
	print '  <blastn-results>  Results of BLASTN classifier.'
	print '  <nb-results>      Results of NB classifier with T=0.'
	print '  <E-value>         Ignore hits with an E-value above this threshold.'
	print '  <percentage>      Use all hits with a bit score with this percentage of the '
	print '                       hit with the highest bit score to classify a fragment.'
	print '  <epsilon>         Use all models with a likelihood at most epsilon times smaller'
	print '                       than the maximum likelihood model to classify a fragment.'
	print '  <results-file>    File to write classification results to.'
	print ''
	print 'Typical usage:'
	print 'python LCA+Epsilon-NB.py blastn_results.txt nb_results.txt 1E-5 15 1E10 lca+epsilon-nb_results.txt'
	print ''
	exit()

blastnResults = sys.argv[1]
nbResults = sys.argv[2]
eValueThreshold = float(sys.argv[3])
percentageThreshold = float(sys.argv[4]) / 100.0
epsilon = float(sys.argv[5])
resultsFile = sys.argv[6]

# read in complete taxonomy of each accession number
strainToTaxonomy = {}
accessionToTaxonomy = {}
for line in  fileinput.input(['taxonomy.txt']):
	lineSplit = line.split('\t')
	
	accession = lineSplit[0]
	taxonomy = lineSplit[1].split(';')
	
	accessionToTaxonomy[accession] = taxonomy
	
	strain = taxonomy[7].replace(':', '_').replace(' ', '_')
	strainToTaxonomy[strain] = taxonomy
		
# determine all hits above the specified E-value threshold for each fragment
print 'Determing hits for each query fragment...'
fragmentHits = {}
for line in fileinput.input([blastnResults]):
	if '# Query:' in line:
		fragmentId = line[line.rfind(' '):]
		fragmentHits[fragmentId.strip()] = [[0, ['u']*8]]
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
blastnTopTaxonomies = {}
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
				
	blastnTopTaxonomies[fragmentId] = topTaxonomy
	
# calculate epislon-NB classification
print 'Calculating epsilon-NB classification for each fragment...'
nbTopTaxonomies = {}
bHeaderLine = True
strains = []
for line in fileinput.input([nbResults]):
	lineSplit = line.split('\t')
	
	if bHeaderLine:
		for i in xrange(3, len(lineSplit)):
			strains.append(lineSplit[i].strip())
		bHeaderLine = False
		continue
	
	fragmentId = lineSplit[0]

	# find maximum likelihood model
	topModel = ''
	topLogLikelihood = -1e100
	for i in xrange(3, len(lineSplit)):
		if float(lineSplit[i]) > topLogLikelihood:
			topModel = strains[i-3]
			topLogLikelihood = float(lineSplit[i])
	
	# find all NB models within a given distance of the maximum likelihood model
	topTaxonomy = list(strainToTaxonomy[topModel])
	for i in xrange(3, len(lineSplit)):
		logLikelihood = float(lineSplit[i])

		if epsilon + logLikelihood >= topLogLikelihood:
			t = strainToTaxonomy[strains[i-3]]
			
			for r in xrange(0, 8):
				if topTaxonomy[r] != t[r]:
					topTaxonomy[r] = 'unclassified'
								
	nbTopTaxonomies[fragmentId] = topTaxonomy
	
# Calculate consensus classification between LCA and epsilon-NB
print 'Calculating LCA + Epsilon-NB classification...'
topTaxonomies = {}
for fragmentId in nbTopTaxonomies:
	topTaxonomy = nbTopTaxonomies[fragmentId]
	blastnTaxonomy = blastnTopTaxonomies[fragmentId]
	
	for r in xrange(0, 8):
		if topTaxonomy[r] != blastnTaxonomy[r]:
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