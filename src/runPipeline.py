#!python

## $Id$

#########################
## runPipeline.py - main pipeline driver for metAMOS
#########################


## The usual library dependencies
import os
import sys
import string
import time
import BaseHTTPServer
import getopt
import re
import subprocess
import webbrowser
import multiprocessing
from operator import itemgetter
from ruffus import *


## Setting up paths
INITIAL_SRC   = "%s%ssrc"%(sys.path[0], os.sep)

sys.path.append(INITIAL_SRC)
import utils
sys.path.append(utils.INITIAL_UTILS)

## Get start time
t1 = time.time()

## Hardcode a k-mer size
DEFAULT_KMER  = 31

def usage():
    print "usage: runPipeline.py [options] -d projectdir"
    print "   -h = <bool>:   print help [this message]"
    print "   -v = <bool>:   verbose output? (default = NO)"
    print "   -d = <string>: directory created by initPipeline (REQUIRED)"

    print "\n[options]: [pipeline_opts] [misc_opts]"
    print "\n[pipeline_opts]: options that affect the pipeline execution"

    print "Pipeline consists of the following steps:"
    print "  Preprocess, Assemble, FindORFS, MapReads, Abundance, Annotate,"
    print "  Scaffold, Propagate, Classify, Postprocess"

    print "Each of these steps can be referred to by the following options:" 
    print "   -f = <string>: force this step to be run (default = NONE)"
    print "   -s = <string>: start at this step in the pipeline (default = Preprocess)"
    print "   -e = <string>: end at this step in the pipeline (default = Postprocess)"
    print "   -n = <string>: step to skip in pipeline (default=NONE)"
    print "   -j = <bool>:   just output all of the programs and citations then exit (default = NO)"

    print "\nFor each step you can fine-tune the execution as follows"
    print "[Preprocess]"
    print "   -t = <bool>:   filter input reads? (default = NO)"
    print "   -q = <bool>:   produce FastQC quality report for reads with quality information (fastq or sff)? (default = NO)"
    print "[Assemble]"
    print "   -a = <string>: genome assembler to use (default = SOAPdenovo)"
    print "   -k = <kmer size>: k-mer size to be used for assembly (default = " + str(DEFAULT_KMER) +  ")"
    print "   -o = <int>>:   min overlap length"
    print "[MapReads]"
    print "   -m = <string>: read mapper to use? (default = bowtie)"
    print "   -i = <bool>:   save bowtie (i)ndex? (default = NO)"
    print "   -b = <bool>:   create library specific per bp coverage of assembled contigs (default = NO)"
    print "[FindORFS]"
    print "   -g = <string>: gene caller to use (default=FragGeneScan)"
    print "   -l = <int>:    min contig length to use for ORF call (default = 300)"
    print "   -x = <int>>:   min contig coverage to use for ORF call (default = 3X)"
    print "[Classify]"
    print "   -c = <string>: classifier to use for annotation (default = FCP)"
    print "   -u = <bool>:   annotate unassembled reads? (default = NO)"

    print "\n[misc_opts]: Miscellaneous options"
    print "   -r = <bool>:   retain the AMOS bank?  (default = NO)"
    print "   -p = <int>:    number of threads to use (be greedy!) (default=1)"
    print "   -4 = <bool>:   454 data? (default = NO)"    

try:
    opts, args = getopt.getopt(sys.argv[1:], "hrbd:s:e:o:k:c:a:n:p:qtf:vm:4g:iul:x:",\
                                   ["help", \
                                        "retainBank", \
                                        "libspeccov",\
                                        "projectdir",\
                                        "startat",\
                                        "endat", \
                                        "minoverlap",\
                                        "kmersize",\
                                        "classifier",\
                                        "assembler",\
                                        "skipsteps",\
                                        "threads",\
                                        "fastqreport",\
                                        "filter",\
                                        "forcesteps",\
                                        "verbose",\
                                        "mapper",\
                                        "454",\
                                        "genecaller",\
                                        "bowtieindex",\
                                        "unassembledreads",\
                                        "minlen",\
                                        "mincov"])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

supported_programs = {}
supported_genecallers = ["FragGeneScan", "fraggenescan","metagenemark","glimmermg"]
supported_assemblers = ["soap","soapdenovo","newbler","ca","velvet","metavelvet",\
                            "metaidba","sparse","sparseassembler","minimus"]
supported_mappers = ["bowtie"]
supported_abundance = ["metaphyler"]
supported_classifiers = ["FCP","fcp","PhyloSift","phylosift","phmmer","blast",\
                             "metaphyler"]
supported_scaffolders = ["bambus2"]
supported_programs["findorfs"] = supported_genecallers
supported_programs["assemble"] = supported_assemblers
supported_programs["mapreads"] = supported_mappers
supported_programs["abundance"] = supported_abundance
supported_programs["classify"] = supported_classifiers
supported_programs["scaffold"] = supported_scaffolders

# Why 4 times????
pub_dict = {}
pub_dict["fraggenescan"] = "Li et al. "
pub_dict["fraggenescan"] = "Li et al. "
pub_dict["fraggenescan"] = "Li et al. "
pub_dict["fraggenescan"] = "Li et al. "

allsteps = ["Preprocess","Assemble","FindORFS","Abundance","Annotate",\
                "Scaffold","Propagate","Classify","Postprocess"]

## Need comments here and further down

output = None
reads = None
quals = None
format = None
savebtidx = False
verbose = False
bowtie_mapping = 1
startat = None
stopat = None
filter = False
forcesteps = []
skipsteps = []
run_fastqc = False
run_metaphyler = False
runfast = False
cls = None
retainBank = False
asm = "soap"
orf = "metagenemark"
fff = ""
readlen = 75
fqlibs = {}
fqfrags = []
rlibs = []
mapper = "bowtie"
ctgbpcov = False
min_ctg_len = 300
min_ctg_cvg = 3
annotate_unassembled = False
output_programs = 0
settings = utils.Settings(DEFAULT_KMER, multiprocessing.cpu_count() - 1, "")

for o, a in opts:
    if o in ("-v","--verbose"):
        utils.Settings.VERBOSE = True
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-i","--indexbowtie"):
        bowtie_mapping = 1
    elif o in ("-j","--justprogs"):
        output_programs = 1
    elif o in ("-u","--unassembledreads"):
        annotate_unassembled = 1
    elif o in ("-x","--xcov"):
        min_ctg_cvg = int(a)
    elif o in ("-l","--lencontigorf"):
        min_ctg_len = int(a)
    elif o in ("-b","--bpctgcov"):
        bpctgcov = True
    elif o in ("-s","--startat"):
        startat = a
        if startat not in allsteps:
            print "cannot start at %s, step does not exist in pipeline"%(startat)
            print allsteps 
    elif o in ("-e","--endat"):
        pass
    elif o in ("-o", "--minoverlap"):
        pass
    elif o in ("-k", "--kmersize"):
        utils.Settings.kmer = int(a)
    elif o in ("-4", "--454"):
        fff = "-454"
    elif o in ("-f", "--forcesteps"):
        forcesteps = a.split(",")
    elif o in ("-n", "--skipsteps"):
        skipsteps = a.split(",")
    elif o in ("-p", "--threads"):
        utils.Settings.threads = int(a)
    elif o in ("-d", "--projectdir"):
        utils.Settings.rundir = os.path.abspath(a)
        if not os.path.exists(a):
          print "project dir %s does not exist!"%(settings.rundir)
          usage()
          sys.exit(1)
    elif o in ("-q", "--fastqc"):
        run_fastqc = True
    elif o in ("-t", "--filter"):
        filter = True

    elif o in ("-m", "--mapper"):
        mapper = a
    elif o in ("-r", "--retainBank"):
        retainBank = True
    elif o in ("-c", "--classifier"):
        cls = a
        if cls == "phylosift" or cls == "PhyloSift":
            cls = "phylosift"
        if cls not in supported_classifiers:
            print "!!Sorry, %s is not a supported classification method. Using FCP instead"%(cls)
            cls = "fcp"
    elif o in ("-a","--assembler"):
        asm = a.lower()
        if asm == "metaidba":
            bowtie_mapping = 1
        if asm not in supported_assemblers:
            print "!!Sorry, %s is not a supported assembler. Using SOAPdenovo instead"%(asm)
            asm = "soap"
        
    elif o in ("-g","--genecaller"):
        orf = a
        if orf not in supported_genecallers:
            print "!!Sorry, %s is not a supported gene caller. Using FragGeneScan instead"%(orf)
            orf = "fraggenescan"
    elif o in ("-f","--fastest"):
        runfast = True
    elif o in ("-b","--savebowtieidx"):
        savebtidx = True
    else:
        assert False, "unhandled option"

if not os.path.exists(settings.rundir) or settings.rundir == "":
    print "project dir %s does not exist!"%(settings.rundir)
    usage()
    sys.exit(1)

#parse frag/libs out of pipeline.ini out of rundir
inifile = settings.rundir+os.sep+"pipeline.ini"
inf = open(inifile,'r')
libs = []
readlibs = []
readobjs = []
frgs = []
format = ""
mean = 0
stdev = 0
mmin = 0
mmax = 0
mated = True
interleaved = False
innie = True

# This should be an option somewhere and probably belongs to initPipeline
linkerType = "titanium"
frg = ""
f1 = ""
f2 = ""
currlibno = 0
newlib = ""
usecontigs = False
libadded = False

for line in inf:
    line = line.replace("\n","")
    if "#" in line:
        continue
    elif "asmcontigs:" in line:
        asmc = line.replace("\n","").split("\t")[-1]
        if len(asmc) <= 2:
            continue
        utils.run_process(settings, "cp %s %s/Assemble/out/%s"%(asmc,settings.rundir,\
                          "proba.asm.contig"),"RunPipeline")
        usecontigs = True
        asm = "none"
        bowtie_mapping = 1
    elif "format:" in line:

        if f1 and not libadded:
            nread1 = utils.Read(format,f1,mated,interleaved)
            readobjs.append(nread1)
            nread2 = ""
            nlib = utils.readLib(format,mmin,mmax,nread1,nread2,mated,interleaved,innie,linkerType)
            readlibs.append(nlib)
        libadded = False
        format = line.replace("\n","").split("\t")[-1]
    elif "mated:" in line:
        mated = utils.str2bool(line.replace("\n","").split("\t")[-1])
    elif "interleaved:" in line:
        interleaved = utils.str2bool(line.replace("\n","").split("\t")[-1])
    elif "innie:" in line:
        innie = utils.str2bool(line.replace("\n","").split("\t")[-1])
    elif "linker:" in line:
        linkerType = line.replace("\n","").split("\t")[-1]
    elif "f1:" in line:
        data = line.split("\t")

        fqlibs[data[0]] = data[1]
 
        f1 = "%s/Preprocess/in/%s"%(settings.rundir,data[1].split(",")[0])
        inf = data[1].split(",")
        mean = int(inf[3])
        stdev = int(inf[4])
        mmin = int(inf[1])
        mmax = int(inf[2])
        libs.append(f1)

    elif "f2:" in line:
        data = line.split("\t")

        fqlibs[data[0]] = data[1]
        f2 = "%s/Preprocess/in/%s"%(settings.rundir,data[1].split(",")[0])
        inf = data[1].split(",")
        mean = int(inf[3])
        stdev = int(inf[4])
        mmin = int(inf[1])
        mmax = int(inf[2])
        libs.append(f2)
        
        nread1 = utils.Read(format,f1,mated,interleaved)
        readobjs.append(nread1)
        nread2 = utils.Read(format,f2,mated,interleaved)
        readobjs.append(nread2)
        nlib = utils.readLib(format,mmin,mmax,nread1,nread2,mated,interleaved,\
                             innie,linkerType)
        readlibs.append(nlib)
        libadded = True
    elif "frg" in line:

        data = line.split("\t")
        frg = "%s/Preprocess/in/%s"%(settings.rundir,data[1].split(",")[0])
        mated = False
        f1 = frg
        libs.append(frg)
if f1 and not libadded:
    nread1 = utils.Read(format,f1,mated,interleaved)
    readobjs.append(nread1)
    nread2 = ""
    nlib = utils.readLib(format,mmin,mmax,nread1,nread2,mated,interleaved,innie,\
                         linkerType)
    readlibs.append(nlib)

if len(readlibs) > 1 and asm == "metaidba":
    print "ERROR: meta-IDBA only supports 1 library, please select different assembler or reduce libraries"
    sys.exit(1)

infile = ""

if (format == "fastq" and mated):
    infile = f1
elif (format == "fastq" and not mated):
    infile = frg
elif (format == "fasta" and not mated):
    infile = frg
elif format == "sff":
    if frg == "":
       infile = f1
    else:
       infile = frg

readpaths = []
filtreadpaths = []
for lib in readlibs:
   for read in lib.reads:
      readpaths.append("%s/Preprocess/in/"%(settings.rundir)+read.fname)
      filtreadpaths.append("%s/Preprocess/out/"%(settings.rundir)+read.fname)

if "Preprocess" in forcesteps:
   for path in readpaths:
      utils.run_process(settings, "touch %s"%(path),"RunPipeline")
utils.Settings.readpaths = readpaths

asmfiles = []

for lib in readlibs:
    if "MapReads" in forcesteps:
        utils.run_process(settings, \
           "touch %s/Assemble/out/%s.asm.contig"%(settings.rundir,settings.PREFIX),\
           "RunPipeline")
    if "Assemble" in forcesteps:
        utils.run_process(settings, \
           "touch %s/Preprocess/out/lib%d.seq"%(settings.rundir,lib.id),\
           "RunPipeline")
    asmfiles.append("%s/Preprocess/out/lib%d.seq"%(settings.rundir,lib.id))


utils.Settings.asmfiles = asmfiles

if "Assemble" not in skipsteps and "Assemble" in forcesteps:
    utils.run_process(settings, \
          "rm %s/Assemble/out/%s.asm.contig"%(settings.rundir,settings.PREFIX),\
          "RunPipeline")

if "FINDORFS" in forcesteps or "findorfs" in forcesteps or "FindORFS" in forcesteps:
   utils.run_process(settings, \
          "rm %s/FindORFS/out/%s.faa"%(settings.rundir,settings.PREFIX),"RunPipeline")
   utils.run_process(settings, \
          "rm %s/FindORFS/out/%s.fna"%(settings.rundir,settings.PREFIX),"RunPipeline")
   utils.run_process(settings, \
          "touch %s/Assemble/out/%s.asm.contig"%(settings.rundir,settings.PREFIX),\
          "RunPipeline")

if "Annotate" in forcesteps:
   utils.run_process(settings, \
          "rm %s/Annotate/out/%s.hits"%(settings.rundir,settings.PREFIX),"RunPipeline")

if "Abundance" in forcesteps:
   utils.run_process(settings, \
          "touch %s/FindORFS/out/%s.faa"%(settings.rundir,settings.PREFIX),\
          "RunPipeline")
   utils.run_process(settings, \
          "rm %s/Abundance/out/%s.taxprof.pct.txt"%(settings.rundir,settings.PREFIX),\
          "RunPipeline")

if "Scaffold" in forcesteps:
    utils.run_process(settings, \
          "rm %s/Scaffold/out/%s.scaffolds.final"%(settings.rundir,settings.PREFIX),\
          "RunPipeline")

if "FindScaffoldORFS" in forcesteps:
    utils.run_process(settings, \
          "touch %s/Scaffold/out/%s.linearize.scaffolds.final"%(settings.rundir,settings.PREFIX),\
          "RunPipeline")

if "Propagate" in forcesteps:
    utils.run_process(settings, "touch %s/DB/class_key.tab"%(settings.METAMOS_UTILS),\
                      "RunPipeline")

if __name__ == "__main__":
    print "Starting metAMOS pipeline"
    settings = utils.initConfig(settings.kmer, settings.threads, settings.rundir)

    import preprocess
    import assemble
    import mapreads
    import findorfs
    import findreps
    import abundance
    import annotate
    import scaffold
    import findscforfs
    import propagate
    import classify
    import postprocess

    # initialize submodules
    preprocess.init(readlibs, skipsteps, asm, run_fastqc,filter)
    assemble.init(readlibs, skipsteps, asm, usecontigs)
    mapreads.init(readlibs, skipsteps, asm, mapper, savebtidx,ctgbpcov)
    findorfs.init(readlibs, skipsteps, asm, orf, min_ctg_len, min_ctg_cvg)
    findreps.init(readlibs, skipsteps)
    annotate.init(readlibs, skipsteps, cls)
    abundance.init(readlibs, skipsteps, forcesteps, cls)
    scaffold.init(readlibs, skipsteps, retainBank, asm)
    findscforfs.init(readlibs, skipsteps, orf)
    propagate.init(readlibs, skipsteps, cls)
    classify.init(readlibs, skipsteps, cls)
    postprocess.init(readlibs, skipsteps, cls)

    try:
       dlist = []
       pipeline_printout(sys.stdout,[preprocess.Preprocess,assemble.Assemble, \
                         findorfs.FindORFS, findreps.FindRepeats, annotate.Annotate, \
                         abundance.Abundance, scaffold.Scaffold, \
                         findscforfs.FindScaffoldORFS, propagate.Propagate, \
                         classify.Classify, postprocess.Postprocess], verbose=1)
       pipeline_printout_graph (   'flowchart.svg',
                            'svg',
                            [postprocess.Postprocess],
                            no_key_legend = True)
       pipeline_run([preprocess.Preprocess, assemble.Assemble,findorfs.FindORFS, \
                    findreps.FindRepeats, annotate.Annotate, abundance.Abundance, \
                    scaffold.Scaffold, findscforfs.FindScaffoldORFS, \
                    propagate.Propagate, classify.Classify, postprocess.Postprocess],\
                    verbose = 1) 

       #multiprocess threads
       t2 = time.time()
       elapsed = float(t2)-float(t1)

       #print elapsed
       print "done! pipeline took %.2f minutes"%(float(elapsed)/float(60.0))
    except JobSignalledBreak:
       print "Done with errors\n"
