#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser, multiprocessing
from operator import itemgetter

INITIAL_SRC   = "%s%ssrc"%(sys.path[0], os.sep)
DEFAULT_KMER  = 31

sys.path.append(INITIAL_SRC)
import utils

t1 = time.time()
sys.path.append(utils.INITIAL_UTILS)
from ruffus import *

def usage():
    print "usage: runPipeline.py [options] -d projectdir (required)"
    print "options:  -a <assembler> -k <kmer size> -c <classification method> -m <enable metaphyler?> -p <num threads>  "
    print "-h: help?"
    print "-r: retain the AMOS bank?  (default = NO)"
    print "-m: read mapper to use? (default = bowtie)"
    print "-d = <project dir>: directory created by initPipeline"
    print "-s = <runPipeline step>: start at this step in the pipeline"
    print "-e = <runPipeline step>: end at this step in the pipeline"
    print "-o = <int>>: min overlap length"
    print "-k = <int>: kmer size for assembly"
    print "-c = <classifier>: classifier to use for annotation"
    print "-a = <assembler>: genome assembler to use"
    print "-n = <runPipeline step>: step to skip in pipeline"
    print "-p = <int>: number of threads to use (be greedy!)"
    print "-q: produce FastQC quality report for reads with quality information (fastq or sff)? (default = NO)"
    print "-t: filter input reads? (default = NO)"
    print "-f = <runPipeline step>: force this step to be run"
    print "-v: verbose output? (default = NO)"
    print "-4: 454 data? (default = NO)"
    print "-b: save (bowtie) index? (default = NO)"
    
    #print "options: annotate, stopafter, startafter, fq, fa"

try:
    opts, args = getopt.getopt(sys.argv[1:], "hrbd:s:e:o:k:c:a:n:p:qtf:vm:4g:", ["help", "retainBank""bowtie","projectdir","startat","endat", "minoverlap","kmersize","classifier","assembler","skipsteps","threads","filter","forcesteps","verbose","mapper","454","genecaller"])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

supported_genecallers = ["fraggenescan","metagenemark"]
supported_assemblers = ["soap","soapdenovo","newbler","ca","velvet","metavelvet","metaidba","sparse","sparseassembler","minimus"]
supported_mappers = ["bowtie"]
supported_abundance = ["metaphyler"]
supported_classifiers = ["FCP","fcp","PhyloSift","phylosift","phmmer","blast","metaphyler"]
supported_scaffolders = ["bambus2"]

allsteps = ["Preprocess","Assemble","FindORFS","Abundance","Annotate","Scaffold","Propagate","Classify","Postprocess"]
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
settings = utils.Settings(DEFAULT_KMER, multiprocessing.cpu_count() - 1, "")

for o, a in opts:
    if o in ("-v","--verbose"):
        utils.Settings.VERBOSE = True
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-b","--bowtie"):
        bowtie_mapping = 1
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
        #print o,a
        forcesteps = a.split(",")
        #print forcesteps
    elif o in ("-n", "--skipsteps"):
        #print o, a
        skipsteps = a.split(",")
        #print skipsteps
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
        #blast,fcp,etc 
        #default: fcp?
        cls = a#"phmmer"
        if cls == "phylosift" or cls == "PhyloSift":
            cls = "phylosift"
        if cls not in supported_classifiers:
            print "!!Sorry, %s is not a supported classification method. Using FCP instead"%(cls)
            cls = "fcp"
    elif o in ("-a","--assembler"):
        #maximus,CA,soap
        #default: maximus?
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
        #tweak all parameters to run fast
        #bambus2, use SOAP, etc
        runfast = True
    elif o in ("-b","--savebowtieidx"):
        savebtidx = True
    else:
        assert False, "unhandled option"

    #sys.exit(2)

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
        #move to proba.asm.contigs
        #skip Assembly
        asmc = line.replace("\n","").split("\t")[-1]
        if len(asmc) <= 2:
            continue
        utils.run_process(settings, "cp %s %s/Assemble/out/%s"%(asmc,settings.rundir,"proba.asm.contig"))
        #skipsteps.append("Assemble")
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
    elif "f1:" in line:# or "f2:" in line:
        data = line.split("\t")

        fqlibs[data[0]] = data[1]
        #f1 = data[1].split(",")[0]
        f1 = "%s/Preprocess/in/%s"%(settings.rundir,data[1].split(",")[0])
        inf = data[1].split(",")
        mean = int(inf[3])
        stdev = int(inf[4])
        mmin = int(inf[1])
        mmax = int(inf[2])
        libs.append(f1)

    elif "f2:" in line:# or "f2:" in line:
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
        nlib = utils.readLib(format,mmin,mmax,nread1,nread2,mated,interleaved,innie,linkerType)
        readlibs.append(nlib)
        libadded = True
    elif "frg" in line:

        data = line.split("\t")
        #frg = data[1]
        frg = "%s/Preprocess/in/%s"%(settings.rundir,data[1].split(",")[0])
        mated = False
        f1 = frg
        #fqfrags[data[0]] = data[1]
        #frgs.append(data[1])
        libs.append(frg)
if f1 and not libadded:
    nread1 = utils.Read(format,f1,mated,interleaved)
    readobjs.append(nread1)
    nread2 = ""
    nlib = utils.readLib(format,mmin,mmax,nread1,nread2,mated,interleaved,innie,linkerType)
    readlibs.append(nlib)
    #libadded = True


if len(readlibs) > 1 and asm == "metaidba":
    print "ERROR: meta-IDBA only supports 1 library, please select different assembler or reduce libraries"
    sys.exit(1)

infile = ""

#for lib in libs:
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

#if asm == "soap":
if "Preprocess" in forcesteps:
   for path in readpaths:
      utils.run_process(settings, "touch %s"%(path))
utils.Settings.readpaths = readpaths

asmfiles = []
#if asm == "soap"

for lib in readlibs:
    #print "touch"
    if "MapReads" in forcesteps:
        utils.run_process(settings, "touch %s/Assemble/out/%s.asm.contig"%(settings.rundir,settings.PREFIX))
    if "Assemble" in forcesteps:
        #print lib.id
        utils.run_process(settings, "touch %s/Preprocess/out/lib%d.seq"%(settings.rundir,lib.id))

    asmfiles.append("%s/Preprocess/out/lib%d.seq"%(settings.rundir,lib.id))
utils.Settings.asmfiles = asmfiles

if "Assemble" not in skipsteps and "Assemble" in forcesteps:
    utils.run_process(settings, "rm %s/Assemble/out/%s.asm.contig"%(settings.rundir,settings.PREFIX))

if "FindORFS" in forcesteps:
   utils.run_process(settings, "rm %s/FindORFS/out/%s.faa"%(settings.rundir,settings.PREFIX))

if "Annotate" in forcesteps:
   utils.run_process(settings, "rm %s/Annotate/out/%s.hits"%(settings.rundir,settings.PREFIX))

if "Abundance" in forcesteps:
   utils.run_process(settings, "touch %s/FindORFS/out/%s.faa"%(settings.rundir,settings.PREFIX))
   utils.run_process(settings, "rm %s/Abundance/out/%s.taxprof.pct.txt"%(settings.rundir,settings.PREFIX))

if "Scaffold" in forcesteps:
    #utils.run_process(settings, "touch %s/Assemble/out/%s.asm.contig"%(settings.rundir,settings.PREFIX))
    utils.run_process(settings, "rm %s/Scaffold/out/%s.scaffolds.final"%(settings.rundir,settings.PREFIX))

if "FindScaffoldORFS" in forcesteps:
    utils.run_process(settings, "touch %s/Scaffold/out/%s.linearize.scaffolds.final"%(settings.rundir,settings.PREFIX))

if "Propagate" in forcesteps:
    utils.run_process(settings, "touch %s/DB/class_key.tab"%(settings.METAMOS_UTILS))

if __name__ == "__main__":
    #pid = start_http()
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
    mapreads.init(readlibs, skipsteps, asm, mapper, savebtidx)
    findorfs.init(readlibs, skipsteps, asm, orf)
    findreps.init(readlibs, skipsteps)
    annotate.init(readlibs, skipsteps, cls)
    abundance.init(readlibs, skipsteps, forcesteps, cls)
    scaffold.init(readlibs, skipsteps, retainBank, asm)
    findscforfs.init(readlibs, skipsteps, orf)
    propagate.init(readlibs, skipsteps, cls)
    classify.init(readlibs, skipsteps, cls)
    postprocess.init(readlibs, skipsteps, cls)

    try:
       #files = os.listdir(".")
       dlist = []
       pipeline_printout(sys.stdout,[preprocess.Preprocess,assemble.Assemble, findorfs.FindORFS, findreps.FindRepeats, annotate.Annotate, abundance.Abundance, scaffold.Scaffold, findscforfs.FindScaffoldORFS, propagate.Propagate, classify.Classify, postprocess.Postprocess], verbose=1)
       pipeline_printout_graph (   'flowchart.svg',
                            'svg',
                            [postprocess.Postprocess],
                            no_key_legend = True)
       pipeline_run([preprocess.Preprocess,assemble.Assemble,findorfs.FindORFS, findreps.FindRepeats, annotate.Annotate, abundance.Abundance, scaffold.Scaffold, findscforfs.FindScaffoldORFS, propagate.Propagate, classify.Classify, postprocess.Postprocess], verbose = 1) 
       #multiprocess threads
       t2 = time.time()#clock()
       elapsed = float(t2)-float(t1)
       #print elapsed
       print "done! pipeline took %.2f minutes"%(float(elapsed)/float(60.0))
    except JobSignalledBreak:
       print "Done with errors\n"
