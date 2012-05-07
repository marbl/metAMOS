# This scrip downloads all complete sequenced bacterial and
# archaeal genomes in the NCBI RefSeq database. Due to the size
# of this file (> 3.5GB) it can take several hours to download.
# Taxonomic information is also obtained by downloading a portion
# of the NCBI Taxonomy database.

import os
import sys
import platform
import tarfile
import urllib
import fileinput
from ftplib import FTP

# Write a file to disk obtained via FTP
class FtpWriter:
	def __init__(self, file):
		self.f = open(file, 'wb')
		self.count = 0

	def __call__(self, block):
		self.f.write(block)
		
		if self.count % 100 == 0:
			print '.',
			
		self.count += 1
		
	def close(self):
		self.f.close()

# Download bacterial and archaeal genomes in NCBI RefSeq
def DownloadGenomes(genomeFile):
	bDownload = True
	if os.path.exists('./' + genomeFile):
		bValidResponse = False
		while not bValidResponse:
			response = raw_input('NCBI genome file ' + genomeFile + ' already exists. Would you like to download the latest version [Y/N]? ')
			if response[0] == 'Y' or response[0] == 'y':
				bDownload = True
				bValidResponse = True
			elif response[0] == 'N' or response[0] == 'n':
				bDownload = False
				bValidResponse = True
			
	if bDownload:
		ncbiFTP = 'ftp.ncbi.nih.gov'
		genomeDir = '/genomes/Bacteria'
		
		# Connect to NBCI's FTP site using an anonymous account
		print 'Connecting to NCBI FTP site (' + ncbiFTP + ')...'
		ftp = FTP(ncbiFTP)
		print ftp.login()
		print '\n'

		# Change to directory containing bacterial and archaeal genomes
		print 'Changing to directory ' + genomeDir
		print ftp.cwd(genomeDir)
		print '\n'

		# Download bacterial and archaeal genomes
		print 'Downloading bacterial and archaeal genomes (' + genomeFile + ')...'
		print '  This file is over 3GB and may take awhile to download.'
		ftpWriter = FtpWriter(genomeFile)
		msg = ftp.retrbinary('RETR ' + genomeFile, ftpWriter, 32*1024*1024)
		print '\n'
		print msg
		ftpWriter.close()

		ftp.quit()
	
# Download NCBI taxonomy database
def DownloadTaxonomy(taxonomyDump):
	ncbiFTP = 'ftp.ncbi.nih.gov'
	taxonomyDir = '/pub/taxonomy'
	
	# Connect to NBCI's FTP site using an anonymous account
	print 'Connecting to NCBI FTP site (' + ncbiFTP + ')...'
	ftp = FTP(ncbiFTP)
	print ftp.login()
	print '\n'

	# Change to directory containing taxonomy files
	print 'Changing to directory ' + taxonomyDir
	print ftp.cwd(taxonomyDir)
	print '\n'

	# Download taxonomy files
	print 'Downloading taxonomy database files...'
	print '  It may take a few minutes to download these files.'
	
	ftpWriter = FtpWriter(taxonomyDump)
	msg = ftp.retrbinary('RETR ' + taxonomyDump, ftpWriter, 32*1024*1024)
	print '\n'
	print msg
	ftpWriter.close()
	
	ftp.quit()

# Decompress genome file
def DecompressGenomes(genomeFile):
	tar = tarfile.open(genomeFile, 'r:gz')
	tar.extractall('./ncbi_genomes/')
	tar.close()
	
# Decompress taxonomy files
def DecompressTaxonomy(taxonomyDump):
	tar = tarfile.open(taxonomyDump, 'r:gz')
	tar.extractall('./taxonomy/')
	tar.close()
	
# Get full taxonomy of all prokaryotes
def BuildTaxonomyFile():
	# read taxon Id number of all contigs
	print 'Extracting taxon Id from each contig...'
	assessionToTaxonId = {}
	accessionToSource = {}
	genomeDirs = os.listdir('./ncbi_genomes/')
	for dir in genomeDirs:
		for filename in os.listdir('./ncbi_genomes/' + dir):
			accession = filename.split('.')[0]
			for line in fileinput.input(['./ncbi_genomes/' + dir + '/' + filename]): 
				if 'SOURCE' in line:
					source = line[len('SOURCE'):].strip()
					accessionToSource[accession] = source.replace('/', '_')
				if '/db_xref="taxon:' in line:
					taxonId = line.split(':')[1]
					taxonId = int(taxonId[0:taxonId.rfind('\"')])
					assessionToTaxonId[accession] = taxonId
					fileinput.close()
					break
	
	print 'Number of contigs: ' + str(len(assessionToTaxonId))

	# extract taxonomy of each contig
	print 'Extracting taxonomy of each contig...'

	nodeIdToName = {}
	for line in fileinput.input(['./taxonomy/names.dmp']): 
		lineSplit = line.split('|')
		id = int(lineSplit[0])
		name = lineSplit[1].strip()
		type = lineSplit[3].strip()
		
		if type == 'scientific name':
			nodeIdToName[id] = name
		
	taxonIdToNode = {}
	for line in fileinput.input(['./taxonomy/nodes.dmp']): 
		lineSplit = line.split('|')
		taxonId = int(lineSplit[0])
		parentId = int(lineSplit[1])
		rank = lineSplit[2].strip()
		
		taxonIdToNode[taxonId] = [rank, parentId]
				
	ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
	fout = open('taxonomy.txt', 'w')
	for assession in assessionToTaxonId:
		taxonId = assessionToTaxonId[assession]
		source = accessionToSource[assession]

		fout.write(assession + '\t')
		
		taxonomy = ['','','','','','','','']
		rankIndex = 0
		while nodeIdToName[taxonId] != 'root':
			node = taxonIdToNode[taxonId]
			if node[0] in ranks:
				while rankIndex < ranks.index(node[0]):
					if rankIndex != 0:
						taxonomy[rankIndex] = nodeIdToName[taxonId] + ' (' + ranks[rankIndex] + ')'
					else:
						taxonomy[rankIndex] = source
					rankIndex += 1
					
				taxonomy[ranks.index(node[0])] = nodeIdToName[taxonId]
				rankIndex += 1
			
			taxonId = node[1]
		
		for r in xrange(7, -1, -1):
			fout.write(taxonomy[r] + ';')
		fout.write('\n')
		
	fout.close()

# create genome-level input files
def CreateStrainSeqFiles():
	# determine genome of each sequence
	assessionToGenome = {}
	for line in fileinput.input(['taxonomy.txt']): 
		lineSplit = line.split('\t')
		
		seqId = lineSplit[0]
		category = lineSplit[1].split(';')[7]

		assessionToGenome[seqId] = category
		
	# creat required directories
	if not os.path.exists('./training'):
		os.makedirs('./training')
		os.makedirs('./training/sequences')
		os.makedirs('./training/custom')
	
	# remove any previously created models
	for assession in assessionToGenome:
		genome = assessionToGenome[assession]
		genomeFile = genome.replace(' ', '_')
		genomeFile = genomeFile.replace(':', '_')
		genomeFile += '.fasta'
		
		if os.path.exists('./training/sequences/' + genomeFile):
			os.remove('./training/sequences/' + genomeFile)
		
	# convert genbank files to fasta files
	genomeDirs = os.listdir('./ncbi_genomes/')
	for dir in genomeDirs:
		for filename in os.listdir('./ncbi_genomes/' + dir):
			fullFilename = './ncbi_genomes/' + dir + '/' + filename
			
			# read sequence data from genbank file
			data = open(fullFilename).read()
			origin = data.rfind('ORIGIN')
			start = data.find('1', origin)
			end = data.find('//', origin)
			
			seqLines = data[start:end].split('\n')
			
			seq = ''
			for line in seqLines:
				subseq = line.split()
				seq += ''.join(subseq[1:])
			
			# write fasta file
			assession = filename.split('.')[0]
			genome = assessionToGenome[assession]
			
			print assession
			
			genomeFile = genome.replace(' ', '_')
			genomeFile = genomeFile.replace(':', '_')
			
			fout = open('./training/sequences/' + genomeFile + '.fasta', 'a')
			fout.write('>' + assession + '\n')
			
			index = 0
			while index+60 < len(seq):
				fout.write(seq[index:index+60] + '\n')
				index += 60
			fout.write(seq[index:] + '\n')
			fout.close()
			
	# create training file for genome models
	trainingSet = open('./training/sequences.txt', 'w')
	for filename in os.listdir('./training/sequences/'):
		trainingSet.write('./training/sequences/' + filename + '\n')
	trainingSet.close()
	
# Build Naive Bayes models
def BuildNaiveBayesModels():
	# creat required directories
	if not os.path.exists('./models'):
		os.makedirs('./models')
		os.makedirs('./models/genomes')
		
	if not os.path.exists('./nb-temp-results'):
		os.makedirs('./nb-temp-results')
		
	# build stain-level models
	if platform.system() == 'Windows':
		print 'Building genome-level models...'
		os.system('nb-train-windows.exe -n 10 -t ./taxonomy.txt -s ./training/sequences.txt -m ./models/genomes/')
	else: # assume the system can build the executable from source
		print 'Compiling nb-train...'
		os.chdir('./nb-train-src')
		os.system('make')
		os.chdir('..')
		os.system('cp ./nb-train-src/nb-train .')
		
		print 'Compiling nb-classify...'
		os.chdir('./nb-classify-src')
		os.system('make')
		os.chdir('..')
		os.system('cp ./nb-classify-src/nb-classify .')
		
		print 'Building genome-level models...'
		os.system('./nb-train -n 10 -t ./taxonomy.txt -s ./training/sequences.txt -m ./models/genomes/')

	# create model file for classifying query fragments
	modelFile = open('./models/models.txt', 'w')
	for line in fileinput.input(['./training/sequences.txt']): 
		genome = line[line.rfind('/')+1:line.rfind('.')]
		modelFile.write('./models/genomes/' + genome + '.txt' + '\n')
	modelFile.close()

genomeFile = 'all.gbk.tar.gz'
taxonomyDump = 'taxdump.tar.gz'

print 'This script is maintained by Donovan Parks, Norm MacDonald, and Rob Beiko (beiko@cs.dal.ca).'
print ''
print 'Changes to the NCBI FTP site or NCBI file formats may break this script.'
print 'Please contact us if this script is broken and we will try to resolve the issue.'

print '\n'
print 'Downloading bacterial and archaeal genomes from NCBI:'
DownloadGenomes(genomeFile)
print '\n----------------------------------------------------------\n'

print 'Decompressing genomes:'
DecompressGenomes(genomeFile)
print '\n----------------------------------------------------------\n'

print 'Downloading NCBI taxonomy database:'
DownloadTaxonomy(taxonomyDump)
print '\n----------------------------------------------------------\n'

print 'Decompressing taxonomy files:'
DecompressTaxonomy(taxonomyDump)
print '\n----------------------------------------------------------\n'

print 'Building taxonomy file for genomes:'
BuildTaxonomyFile()
print '\n----------------------------------------------------------\n'

print 'Creating input sequence file for each genome:'
CreateStrainSeqFiles()
print '\n----------------------------------------------------------\n'

print 'Building Naive Bayes models for each genomes:'
BuildNaiveBayesModels()
print '\n----------------------------------------------------------\n'

print 'Installation complete. '
