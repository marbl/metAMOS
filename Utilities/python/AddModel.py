import sys
import os
import fileinput
import platform

if len(sys.argv) != 10:
	print 'AddModel v1.0 by Donovan Parks, Norm MacDonald, and Rob Beiko'
	print ''
	print 'Usage: python AddModel.py <N> <sequence-file> <domain> <phylum> <class> <order> <family> <genus> <species> <strain>'
	print ''
	print 'Required parameters:'
	print '  <N>              Desired n-mer length (must be the same for all models).'
	print '  <sequence-file>  Multi-FASTA file containing sequence data to build model from.'
	print '  <domain>         Domain of model.'
	print '  ...'
	print '  <stain>          Strain of model.'
	print ''
	print 'Typical usage:'
	print '  python AddModel.py ./training/custom/MySeqData.fasta Bacteria Proteobacteria'
	print '            Betaproteobacteria Burkholderiales Burkholderiaceae Ralstonia '
	print '            "Ralstonia pickettii" "Ralstonia pickettii 12D"'
	print ''
	exit()

fastaFile = sys.argv[1]
domain = sys.argv[2]
phylum = sys.argv[3]
classRank = sys.argv[4]
order = sys.argv[5]
family = sys.argv[6]
genus = sys.argv[7]
species = sys.argv[8]
strain = sys.argv[9]

taxonomy = [domain, phylum, classRank, order, family, genus, species, strain]

# get id for each sequence in new sequence file
seqIds = []
for line in fileinput.input([fastaFile]): 
	if line[0] == '>':
		id = line[1:].strip()
		seqIds.append(id)

# add new sequence to taxonomy file
print 'Adding sequence to taxonomy file (taxonomy.txt)...'
fout = open('taxonomy.txt', 'a')
for id in seqIds:
	fout.write(id + '\t')
	
	for r in xrange(0, len(taxonomy)):
		fout.write(taxonomy[r] + ';')
	fout.write('\n')
fout.close()

# add new sequence to sequence file
print 'Adding sequence to sequence file (sequence.txt)...'
fout = open('./training/sequences.txt', 'a')
fout.write('.' + fastaFile + '\n')
fout.close()

# write tempory sequence file
fout = open('./training/__temp__.txt', 'w')
fout.write('.' + fastaFile + '\n')
fout.close()

# build new model
print 'Building model for new sequence...'
if platform.system() == 'Windows':	
	os.system('nb-train-windows.exe 10 ./taxonomy.txt ./training/__temp__.txt ./models/genomes/')
else: # assume the system can build the executable from source
	os.system('./nb-train 10 ./taxonomy.txt ./training/__temp__.txt ./models/genomes/')
	
# removing temporary sequence file
os.chdir('..')
os.remove('./training/__temp__.txt')

# add new model to model file
print 'Adding new model to model file (models.txt)...'
fout = open('./models/models.txt', 'a')

filename = fastaFile[fastaFile.rfind('/')+1:fastaFile.rfind('.')]
fout.write('./models/genomes/' + filename + '.txt' + '\n')
fout.close()

print 'Done.'

