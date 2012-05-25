#!/fs/szdevel/core-cbcb-software/Linux-x86_64/bin/python

# requirements for plots
#  biopython
#  matplotlib
#  numpy

from Bio import SeqIO, Seq
import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np

# represent info related to each sample
class Sample:
    def __init__(self, id, directory, line_type, args):
        self.id = id
        self.directory = directory
        self.args = args
        self.ctg_sizes = {}
        self.sf_sizes = {}
        self.ctg_features = {}
        self.cum_sizes = []
        self.sfcum_sizes = []
        self.list_sizes = []
        self.sfsort_sizes = []
        self.line_type = line_type

# contig size cutoff for contig included in plots
SIZE_CUTOFF = 100

# Takes a fasta file from an assembler and produces three plots
#
# 1. Histogram of sizes
# 2. cummulative size versus contig size
# 3. cummulative size versus contig number
# 4. if # of features are provided, FR curve
#
# INPUT - tab delimited file with sampleID, directory, clinical parameters
#

#############
#
#  Values that may need to be changed based on MetaAMOS RUN
#  Currently they have been changing
#
##############
METAMOS_prefix = "proba"
METAMOS_PATH = "/Scaffold/out/"
Contig_file = ".contigs"
Scaffold_file = ".linearize.scaffolds.final"
METAPHYLER = "/Abundance/out/"

handle = open(sys.argv[1], "rU")
TITLE = sys.argv[2]
samples = [ ]

# process assembly files
for line in handle:

    line = line.strip()

    if(line != ""):
        sample_id, dir, METAMOS_prefix, line_type, args = line.split('\t')
        sample = Sample(sample_id, dir, line_type, args)

        for record in SeqIO.parse(dir + METAMOS_PATH + METAMOS_prefix + Contig_file, "fasta") :
            l = len(record.seq)
            if(l > SIZE_CUTOFF):
                sample.ctg_sizes[record.id] = l


        print " ctg_size " + str(len(sample.ctg_sizes))
        sample.list_sizes = sorted(sample.ctg_sizes.values(), reverse=True)


        for record in SeqIO.parse(dir + METAMOS_PATH + METAMOS_prefix + Scaffold_file, "fasta") :
            l = len(record.seq)
            if(l > SIZE_CUTOFF):
                sample.sf_sizes[record.id] = l


        print " sf_size " + str(len(sample.sf_sizes))
        sample.sfsort_sizes = sorted(sample.sf_sizes.values(), reverse=True)
        samples.append(sample)


### histograms
c = 0
h = []
for s in samples:
    h.append(s.list_sizes)

n, bins, patches = plt.hist(h,  100)
plt.ylabel('Contigs Count')
plt.xlabel('Contigs Size')
plt.title(r'Male and Female Contig Size histogram')
plt.grid(True)
plt.savefig('hist_contigs.png')
plt.close()
c += 1

c = 0
h = []
for s in samples:
    h.append(s.sfsort_sizes)

n, bins, patches = plt.hist(h,  100)
plt.ylabel('Scaffold Count')
plt.xlabel('Scaffold Size')
plt.title(r'Male and Female Scaffold Size histogram')
plt.grid(True)
plt.savefig('hist_scaffold.png')
plt.close()
c += 1
######


# contig size plot
c = 0
for s in samples:
    total = 0
    for x in s.list_sizes:
        total += x
        s.cum_sizes.append(total)

    total = 0
    for x in s.sfsort_sizes:
        total += x
        s.sfcum_sizes.append(total)

    print " id " + s.id
    print " total " + str(total)
    print " total contig " + str(len(s.list_sizes))

    plt.plot(s.list_sizes, s.cum_sizes, s.line_type, label="Contig " + str(s.id))
    c += 1

ax = plt.gca()
ax.set_xlim(ax.get_xlim()[::-1])
plt.ylabel('Total Size')
plt.xlabel('Contig Size')
plt.title(TITLE + 'N50 plot')
plt.legend()
plt.savefig('ContigSizes.png')
plt.close()


# scaffold size plot
c = 0
for s in samples:

    plt.plot(s.sfsort_sizes, s.sfcum_sizes, s.line_type, label="Scaffold " + str(s.id))
    c += 1

ax = plt.gca()
ax.set_xlim(ax.get_xlim()[::-1])
plt.ylabel('Total Size')
plt.xlabel('Scaffold Size')
plt.title(TITLE + ' N50 plot')
plt.legend()
plt.savefig('ScaffoldSizes.png')
plt.close()


c = 0
for s in samples:
    plt.plot(range(1, len(s.cum_sizes) + 1), s.cum_sizes, s.line_type, label="Contig " + str(s.id))
    c += 1

plt.ylabel('Total Size')
plt.xlabel('Cumulative Contig')
plt.title(TITLE + ' Contig (Cumulative) ')
plt.legend(loc=4)
plt.savefig('ContigSizes2.png')
plt.close()



c = 0
for s in samples:
    plt.plot(range(1, len(s.sfcum_sizes) + 1), s.sfcum_sizes, s.line_type,  label="Scaffold " + str(s.id))
    c += 1

plt.ylabel('Total Size')
plt.xlabel('Cumulative Scaffold')
plt.title(TITLE + ' Scaffold (Cumulative) ')
plt.legend(loc=4)
plt.savefig('ScaffoldSizes2.png')
plt.close()



handle.close()

                  
