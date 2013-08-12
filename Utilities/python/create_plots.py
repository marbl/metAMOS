import sys, os, string
ROOT = os.path.dirname(os.path.abspath(__file__))

# requirements for plots
#  matplotlib
#  numpy

import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np

from itertools import groupby

#fasta parser from http://biostar.stackexchange.com/questions/711/correct-way-to-parse-a-fasta-file-in-python
def fasta_iter(fasta_name):
   """
   given a fasta file. yield tuples of header, sequence
   """
   try:
       fh = open(fasta_name, 'r')
       # ditch the boolean (x[0]) and just keep the header or sequence since
       # we know they alternate.
       faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
       for header in faiter:
           # drop the ">"
           header = header.next()[1:].strip()
           # join all sequence lines to one.
           seq = "".join(s.strip() for s in faiter.next())
           yield header, seq
   except Exception, e:
       print "error: %s"%(e)
       pass

#use:

# represent info related to each sample
class Sample:
    def __init__(self, id, directory, line_type, args):
        self.id = id
        self.directory = directory
        self.args = args
        self.ctg_sizes = {}
        self.ctgcvg_values = {}
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

def create_plots(handlef,TITLE):
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
    ORF_file = ".fna"
    SCFORF_PATH = "/FindScaffoldORFS/out/"
    SCFORF_file = ".fna"
    ASSEMBLY = "/Assemble/out/"
    Cvg_file = ".contig.cvg"
    
    handle = open(handlef, "r")
    samples = [ ]
    
    METDIR =""
    # process assembly files
    for line in handle.xreadlines():
    
        line = line.strip()
    
        if(line != ""):
            sample_id, dir, METAMOS_prefix, line_type, args = line.split('\t')
            METDIR = dir
            sample = Sample(sample_id, dir, line_type, args)
    
            for header, seq in (fasta_iter(dir + METAMOS_PATH + METAMOS_prefix + Contig_file)):
                l = len(seq)
                if(l > SIZE_CUTOFF):
                    sample.ctg_sizes[header] = l
    
    
            #print " ctg_size " + str(len(sample.ctg_sizes))
            sample.list_sizes = sorted(sample.ctg_sizes.values(), reverse=True)
    
    
            for header, seq in (fasta_iter(dir + METAMOS_PATH + METAMOS_prefix + Scaffold_file)):
                l = len(seq)
                if(l > SIZE_CUTOFF):
                    sample.sf_sizes[header] = l
    
            #print " sf_size " + str(len(sample.sf_sizes))
            sample.sfsort_sizes = sorted(sample.sf_sizes.values(), reverse=True)
    
            for header, seq in (fasta_iter(dir + ORF_PATH + METAMOS_prefix + ORF_file)):
                l = len(seq)
                if(l > SIZE_CUTOFF):
                    sample.orf_sizes[header] = l
            sample.orflist_sizes = sorted(sample.orf_sizes.values(), reverse=True)
            #print " orf_size " + str(len(sample.orf_sizes))
    
            for header, seq in (fasta_iter(dir + SCFORF_PATH + METAMOS_prefix + SCFORF_file)):
                l = len(seq)
                if(l > SIZE_CUTOFF):
                    sample.scforf_sizes[header] = l
            sample.scforflist_sizes = sorted(sample.scforf_sizes.values(), reverse=True)
            #print " sforf_size " + str(len(sample.scforf_sizes))
    
            cvgs = open(dir + ASSEMBLY + METAMOS_prefix + Cvg_file, 'r')
            counter = 1
            for line in cvgs:
               splitLine = line.replace("\n","").split("\t")
               if len(splitLine) < 2:
                  continue
               sample.ctgcvg_values[counter] = float(splitLine[1])
               counter += 1
            sample.ctgcvg_values = sorted(sample.ctgcvg_values.values(), reverse=True)
            #print " ctgcvg_size " + str(len(sample.ctgcvg_values))
            cvgs.close() 
    
            samples.append(sample)
    
    
    ### histograms
    c = 0
    h = []
    for s in samples:
        c += len(s.list_sizes)
        h.append(s.list_sizes)
    
    if c > 0:
       n, bins, patches = plt.hist(h,  100)
       plt.ylabel('Contigs Count')
       plt.xlabel('Contigs Size')
       plt.title(r'Contig Size histogram')
       plt.grid(True)
       plt.savefig('hist_contigs.png')
       plt.savefig('hist_contigs.pdf', format='pdf')
       plt.close()
    c += 1
    
    c = 0
    h = []
    for s in samples:
        c += len(s.sfsort_sizes)
        h.append(s.sfsort_sizes)
    
    if c > 0:
       n, bins, patches = plt.hist(h,  100)
       plt.ylabel('Scaffold Count')
       plt.xlabel('Scaffold Size')
       plt.title(r'Scaffold Size histogram')
       plt.grid(True)
       plt.savefig('hist_scaffold.png')
       plt.savefig('hist_scaffold.pdf', format='pdf')
       plt.close()
    c += 1
    
    c = 0
    h = []
    for s in samples:
        c += len(s.orflist_sizes)
        h.append(s.orflist_sizes)
    
    if c > 0:
       n, bins, patches = plt.hist(h,  100)
       plt.ylabel('ORF Count')
       plt.xlabel('ORF Size')
       plt.title(r'ORF Size histogram')
       plt.grid(True)
       plt.savefig('hist_orf.png')
       plt.savefig('hist_orf.pdf', format='pdf')
       plt.close()
    c += 1
    
    c = 0
    h = []
    for s in samples:
        c + len(s.scforflist_sizes)
        h.append(s.scforflist_sizes)
    
    if c > 0:
       n, bins, patches = plt.hist(h,  100)
       plt.ylabel('Scaffold ORF Count')
       plt.xlabel('Scaffold ORF Size')
       plt.title(r'Scaffold ORF Size histogram')
       plt.grid(True)
       plt.savefig('hist_scforf.png')
       plt.savefig('hist_scforf.pdf', format='pdf')
       plt.close()
    c += 1
    
    c = 0
    h = []
    for s in samples:
        c += len(s.ctgcvg_values)
        h.append(s.ctgcvg_values)
    
    if c > 0:
       n, bins, patches = plt.hist(h,  20)
       ff2 = open("mapreads.tsv","w")
       ff2.write("cov\tfrequency\n")
       i = 0
       while i < len(n):
           try:
               ff2.write("%.2f\t%f\n"%(bins[i],n[i]))
           except TypeError:
               #incorrect matplotlib version? skip..
               ff2.write("%.2f\t%f\n"%("0.0","0.0"))
               #continue
           i+=1
       ff2.close()
       plt.ylabel('Contig Count')
       plt.xlabel('Contig Coverage')
       plt.title(r'Contig Coverage histogram')
       plt.grid(True)
       plt.savefig('hist_ctgcvg.png')
       plt.savefig('hist_ctgcvg.pdf', format='pdf')
       plt.close()
    c += 1
    ######
    
    # contig size plot
    c = 0
    tsvdata = {}
    ref_cuml = []
    sample_ids = []
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
    
        #print " id " + s.id
        #print " total " + str(total)
        #print " total contig " + str(len(s.list_sizes))
        #print " total ORF " + str(len(s.orflist_sizes))
        #print " total Scaffold ORF " + str(len(s.scforflist_sizes))
        plt.plot(s.list_sizes, s.cum_sizes, s.line_type, label="Contig " + str(s.id))
        sample_ids.append(str(s.id))
        i = 0
        while i < len(s.list_sizes):
            try:
                if c == 0:
                    tsvdata[s.cum_sizes[i]].append(s.list_sizes[i])
                    ref_cuml.append(s.cum_sizes[i])
                else:
                    tsvdata[ref_cuml[i]].append(s.list_sizes[i])
            except KeyError:
                if c == 0:
                    tsvdata[s.cum_sizes[i]] = [s.list_sizes[i]]
                    ref_cuml.append(s.cum_sizes[i])
                else:
                    tsvdata[ref_cuml[i]] = [s.list_sizes[i]]
            i+=1
        c += 1
    ff = open("assemble.tsv", 'w')
    ff.write("total\tcontig\n")#%(sample_ids[0]))
    #if len(sample_ids) > 1:
    #    for id in sample_ids[1:]:
    #        ff.write("\t%s"%(id))
    #ff.write("\n")
    
    tkeys = tsvdata.keys()
    tkeys.sort()
    for row in tkeys:
        
        i = 0
        ff.write("%d\t"%(int(row)))
        while i < len(tsvdata[row])-1:
            ff.write("%d\t"%(tsvdata[row][i+1]))
            i+=1
        ff.write("%d\n"%(tsvdata[row][-1]))
    
    ff.close()
    ax = plt.gca()
    ax.set_xlim(ax.get_xlim()[::-1])
    plt.ylabel('Total Size')
    plt.xlabel('Contig Size')
    plt.title(TITLE + ' N50 plot')
    plt.legend()
    plt.savefig('ContigSizes.png')
    plt.savefig('ContigSizes.pdf', format='pdf')
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
    plt.savefig('ORFSizes.pdf', format='pdf')
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
    plt.savefig('SCFORFSizes.pdf', format='pdf')
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
    plt.savefig('ScaffoldSizes.pdf', format='pdf')
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
    plt.savefig('ContigSizes2.pdf', format='pdf')
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
    plt.savefig('ScaffoldSizes2.pdf', format='pdf')
    plt.close()
    
    
    
    handle.close()
    
                      
