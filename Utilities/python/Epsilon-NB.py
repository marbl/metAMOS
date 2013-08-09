# Classify fragments using Epsilon-NB.

import sys
import fileinput
import math

if len(sys.argv) != 4:
	print 'Epsilon-NB v1.0 by Donovan Parks, Norm MacDonald, and Rob Beiko'
	print ''
	print 'Usage: python Epsilon-NB.py <nb-results> <epsilon> <results-file>'
	print ''
	print 'Required parameters:'
	print '  <nb-results>    Results of NB classifier with T=0.'
	print '  <epsilon>       Use all models with a likelihood at most epsilon times smaller'
	print '                    than the maximum likelihood model to classify a fragment.'
	print '  <results-file>  File to write classification results to.'
	print ''
	print 'Typical usage:'
	print '  python Epsilon-NB.py nb_results.txt 1E10 epsilon-nb_results.txt'
	print ''
	sys.exit()
	
nbResults = sys.argv[1]
if sys.argv[2] == '0' or sys.argv[2] == '0.0' or float(sys.argv[2]) == 0.0:
	epsilon = 0.0
else:
	epsilon = math.log(float(sys.argv[2]))
resultsFile = sys.argv[3]

# read in complete taxonomy of each strain/genome
strainToTaxonomy = {}
for line in  fileinput.input(['taxonomy.txt']):
	lineSplit = line.split('\t')
	taxonomy = lineSplit[1].split(';')
	
	strain = taxonomy[7].replace(':', '_').replace(' ', '_')
	strainToTaxonomy[strain] = taxonomy
		
# calculate epsilon-NB classification
print 'Calculating epsilon-NB classification for each fragment...'
topTaxonomies = {}
# SK - begin changes for metAMOS
topLogTaxonomies = {}
# SK - end changes for metAMOS
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

        # SK - begin changes for metAMOS - we create an array initialized to the likelyhood for each taxonomy level
        topLogTaxonomy = list(strainToTaxonomy[topModel])
        for r in xrange(0, 8):
           topLogTaxonomy[r] = topLogLikelihood
        # SK - end changes for metAMOS

	for i in xrange(3, len(lineSplit)):
		logLikelihood = float(lineSplit[i])

		if epsilon + logLikelihood >= topLogLikelihood:
			t = strainToTaxonomy[strains[i-3]]
			
			for r in xrange(0, 8):
				if topTaxonomy[r] != t[r]:
					topTaxonomy[r] = 'unclassified'
                                        # SK - begin metAMOS changes
                                        if logLikelihood < topLogTaxonomy[r]:
                                           topLogTaxonomy[r] = logLikelihood
                                        # SK - end metAMOS changes
								
	topTaxonomies[fragmentId] = topTaxonomy
        # SK - begin changes for metAMOS
        topLogTaxonomies[fragmentId] = topLogTaxonomy
        # SK - end changes for metAMOS

# write out classification results
print 'Writing out results...'
fout = open(resultsFile, 'w')
fout.write('Fragment Id\tTaxonomic classification\n')
for fragmentId in topTaxonomies:
	fout.write(fragmentId + '\t')
	
	topTaxonomy = topTaxonomies[fragmentId]
	for r in xrange(0, 8):
		fout.write(topTaxonomy[r] + ';')

        # SK - begin changes for metAMOS
        topLogTaxonomy = topLogTaxonomies[fragmentId]
        fout.write('\t')
        for r in xrange(0, 8):
                fout.write("%f;"%(topLogTaxonomy[r]))
        # SK - end changes for metAMOS
	fout.write('\n')
		
fout.close()

print 'Done.'
