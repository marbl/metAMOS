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
        self.orf_sizes = {}
        self.scforf_sizes = {}
        self.cum_sizes = []
        self.sfcum_sizes = []
        self.orfcum_sizes = []
        self.scforfcum_sizes = []
        self.list_sizes = []
        self.orflist_sizes = []
        self.scforflist_sizes = []
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
ORF_PATH = "/FindORFS/out/"
ORF_file = ".orfs.ffn"
SCFORF_PATH = "/FindScaffoldORFS/out/"
SCFORF_file = ".orfs.ffn"

handle = open(sys.argv[1], "r")
TITLE = sys.argv[2]
samples = [ ]

# process assembly files
for line in handle.xreadlines():

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

        for record in SeqIO.parse(dir + ORF_PATH + METAMOS_prefix + ORF_file, "fasta") :
            l = len(record.seq)
            if(l > SIZE_CUTOFF):
                sample.orf_sizes[record.id] = l
        sample.orflist_sizes = sorted(sample.orf_sizes.values(), reverse=True)
        print " orf_size " + str(len(sample.orf_sizes))

        for record in SeqIO.parse(dir + SCFORF_PATH + METAMOS_prefix + SCFORF_file, "fasta") :
            l = len(record.seq)
            if(l > SIZE_CUTOFF):
                sample.scforf_sizes[record.id] = l
        sample.scforflist_sizes = sorted(sample.scforf_sizes.values(), reverse=True)
        print " sforf_size " + str(len(sample.scforf_sizes))


        samples.append(sample)


### histograms
c = 0
h = []
for s in samples:
    h.append(s.list_sizes)

n, bins, patches = plt.hist(h,  100)
plt.ylabel('Contigs Count')
plt.xlabel('Contigs Size')
plt.title(r'Contig Size histogram')
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
plt.title(r'Scaffold Size histogram')
plt.grid(True)
plt.savefig('hist_scaffold.png')
plt.close()
c += 1

c = 0
h = []
for s in samples:
    h.append(s.orflist_sizes)

n, bins, patches = plt.hist(h,  100)
plt.ylabel('ORF Count')
plt.xlabel('ORF Size')
plt.title(r'ORF Size histogram')
plt.grid(True)
plt.savefig('hist_orf.png')
plt.close()
c += 1

c = 0
h = []
for s in samples:
    h.append(s.scforflist_sizes)

n, bins, patches = plt.hist(h,  100)
plt.ylabel('Scaffold ORF Count')
plt.xlabel('Scaffold ORF Size')
plt.title(r'Scaffold ORF Size histogram')
plt.grid(True)
plt.savefig('hist_scforf.png')
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

    total = 0
    for x in s.orflist_sizes:
        total += x
        s.orfcum_sizes.append(total)

    total = 0
    for x in s.scforflist_sizes:
        total += x
        s.scforfcum_sizes.append(total)

    print " id " + s.id
    print " total " + str(total)
    print " total contig " + str(len(s.list_sizes))
    print " total ORF " + str(len(s.orflist_sizes))
    print " total Scaffold ORF " + str(len(s.scforflist_sizes))

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

# ORF size plot
c = 0
for s in samples:

    plt.plot(s.orflist_sizes, s.orfcum_sizes, s.line_type, label="ORF " + str(s.id))
    c += 1

ax = plt.gca()
ax.set_xlim(ax.get_xlim()[::-1])
plt.ylabel('Total Size')
plt.xlabel('ORF Size')
plt.title(TITLE + ' N50 plot')
plt.legend()
plt.savefig('ORFSizes.png')
plt.close()

# Scaffold ORF size plot
c = 0
for s in samples:

    plt.plot(s.scforflist_sizes, s.scforfcum_sizes, s.line_type, label="Scaffold ORF " + str(s.id))
    c += 1

ax = plt.gca()
ax.set_xlim(ax.get_xlim()[::-1])
plt.ylabel('Total Size')
plt.xlabel('Scaffold ORF Size')
plt.title(TITLE + ' N50 plot')
plt.legend()
plt.savefig('SCFORFSizes.png')
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

                  
