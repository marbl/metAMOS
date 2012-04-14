# Classify fragments using NB-BL

import sys
import fileinput
import math

if len(sys.argv) != 4:
	print 'NB-BL v1.0 by Donovan Parks, Norm MacDonald, and Rob Beiko'
	print ''
	print 'Usage: python NB-BL.py <nb-results> <blastn-results> <results-file>'
	print ''
	print 'Required parameters:'
	print '  <blastn-path>   Results of NB classifier with T=0.'
	print '  <query-file>    Results of BLASTN classifier.'
	print '  <results-file>  File to write classification results to.'
	print ''
	print 'Typical usage:'
	print '  python NB-BL.py nb_results.txt blastn_results.txt nb-bl_results.txt'
	print ''
	exit()
	
nbResults = sys.argv[1]
blastnResults = sys.argv[2]
resultsFile = sys.argv[3]

# read in complete taxonomy of each strain/genome
strainToTaxonomy = {}
accessionToTaxonomy = {}
for line in  fileinput.input(['taxonomy.txt']):
	lineSplit = line.split('\t')
	
	accession = lineSplit[0]
	taxonomy = lineSplit[1].split(';')
	
	accessionToTaxonomy[accession] = taxonomy
	
	strain = taxonomy[7].replace(':', '_').replace(' ', '_')
	strainToTaxonomy[strain] = taxonomy
		
# determine smallest blastn E-value for each strain (genome) for each fragment
print 'Processing blastn results...'
blastnEvalues = {}
for line in  fileinput.input([blastnResults]):
	if line[0] == '#':
		continue
		
	lineSplit = line.split('\t')
	queryId = lineSplit[0]
	accession = lineSplit[1]
	evalue = float(lineSplit[10])
	
	taxonomy = accessionToTaxonomy[accession]
	strain = accessionToTaxonomy[accession][7]
	
	evalues = blastnEvalues.get(queryId, {})
	if strain in evalues:
		if evalue < evalues[strain]:
			evalues[strain] = evalue
	else:
		evalues[strain] = evalue
		
	blastnEvalues[queryId] = evalues

# get NB likelihood for each fragment
print 'Processing NB results...'
nbLogLikelihoods = {}
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
	
	logLikelihoods = {}
	for i in xrange(3, len(lineSplit)):
		strain = strains[i-3]
		logLikelihoods[strain] = float(lineSplit[i])
		
	nbLogLikelihoods[fragmentId] = logLikelihoods

# calculate NB-BL score for each fragment against each model
print 'Calculating NB-BL scores...'
fragmentScores = {}
for fragmentId in nbLogLikelihoods:
	scores = {}
	for strain in nbLogLikelihoods[fragmentId]:
		logLikelihood = nbLogLikelihoods[fragmentId][strain]
		
		if (fragmentId not in blastnEvalues) or (strain not in blastnEvalues[fragmentId]):
			# just take NB score
			scores[strain] = logLikelihood
		else:
			evalue = blastnEvalues[fragmentId][strain]
			
			if evalue == 0:
				scores[strain] = logLikelihood + 10000	# big bonus to score
			else:
				scores[strain] = logLikelihood + 12*(4 - math.log(evalue)) # seems like this should be log10, but PhymmBL uses the natural log
		
	fragmentScores[fragmentId] = scores

# write out classification results
print 'Writing out results...'
fout = open(resultsFile, 'w')
fout.write('Fragment Id\tTaxonomic classification\n')
for fragmentId in fragmentScores:
	# find top score
	topScore = -1e100	
	scores = fragmentScores[fragmentId]
	for strain in scores:
		if scores[strain] > topScore:
			topTaxonomy = strainToTaxonomy[strain]
			topScore = scores[strain]
			
	# print results
	fout.write(fragmentId + '\t')
	for r in xrange(0, 8):
		fout.write(topTaxonomy[r] + ';')
	fout.write('\n')
		
fout.close()

print 'Done.'