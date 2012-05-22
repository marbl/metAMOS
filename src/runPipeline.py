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
    opts, args = getopt.getopt(sys.argv[1:], "hrjbd:s:e:o:k:c:a:n:p:qtf:vm:4g:iul:x:",\
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

## always use long names, search will auto-detect abbreviations

supported_programs = {}
supported_genecallers = ["fraggenescan","metagenemark","glimmermg"]
supported_assemblers = ["soapdenovo","newbler","ca","velvet","metavelvet",\
                            "metaidba","sparseassembler","minimus"]
supported_mappers = ["bowtie"]
supported_abundance = ["metaphyler"]
supported_classifiers = ["fcp","phylosift","phmmer","blast",\
                             "metaphyler"]
supported_scaffolders = ["bambus2"]
supported_programs["findorfs"] = supported_genecallers
supported_programs["assemble"] = supported_assemblers
supported_programs["mapreads"] = supported_mappers
supported_programs["abundance"] = supported_abundance
supported_programs["classify"] = supported_classifiers
supported_programs["scaffold"] = supported_scaffolders

## a list of citations for each program used in the pipeline
progname_dict = {}
progname_dict["fraggenescan"] = "FragGeneScan"
progname_dict["metagenemark"] = "MetaGeneMark"
progname_dict["soapdenovo"] = "SOAPdenovo"
progname_dict["metaidba"] = "Meta-IDBA"
progname_dict["velvet"] = "Velvet"
progname_dict["metavelvet"] = "MetaVelvet"
progname_dict["ca"] = "Celera Assembler"
progname_dict["bambus2"] = "Bambus 2"
progname_dict["fcp"] = "FCP,Naive Bayesian Classifier"
progname_dict["bowtie"] = "Bowtie"
progname_dict["blast"] = "BLAST"
progname_dict["phmmer"] = "PHMMER"
progname_dict["phylosift"] = "PhyloSift"
progname_dict["minimus"] = "Minimus"
progname_dict["glimmermg"] = "Glimmer-MG"
progname_dict["sparseassembler"] = "Sparse Assembler"
progname_dict["metaphyler"] = "MetaPhyler"
progname_dict["newbler"] = "Newbler"

pub_dict = {}
pub_dict["fraggenescan"] = "Rho M, Tang H, Ye Y: FragGeneScan: predicting genes in short and error-prone reads. Nucleic Acids Research 2010, 38:e191-e191."
pub_dict["metagenemark"] = "Borodovsky M, Mills R, Besemer J, Lomsadze A: Prokaryotic gene prediction using GeneMark and GeneMark.hmm.Current protocols in bioinformatics editoral board Andreas D Baxevanis et al 2003, Chapter 4:Unit4.6-Unit4.6."
pub_dict["soapdenovo"] = "Li Y, Hu Y, Bolund L, Wang J: State of the art de novo assembly of human genomes from massively parallel sequencing data.Human genomics 2010, 4:271-277."
pub_dict["metaidba"] = "Peng Y, Leung HCM, Yiu SM, Chin FYL: Meta-IDBA: a de Novo assembler for metagenomic data. Bioinformatics 2011, 27:i94-i101."
pub_dict["metavelvet"] = "Namiki T, Hachiya T, Tanaka H, Sakakibara Y: MetaVelvet : An extension of Velvet assembler to de novo metagenome assembly from short sequence reads. In; 2011."
pub_dict["bowtie"] = "Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol. 2009;10(3):R25. Epub 2009 Mar 4."
pub_dict["bambus2"] = "Koren S, Treangen TJ, Pop M. Bambus 2: scaffolding metagenomes. Bioinformatics. 2011 Nov 1;27(21):2964-71. Epub 2011 Sep 16."
pub_dict["fcp"] = "Macdonald NJ, Parks DH, Beiko RG. Rapid identification of high-confidence taxonomic assignments for metagenomic data. Nucleic Acids Res. 2012 Apr 24."
pub_dict["metaphyler"] = "Liu B, Gibbons T, Ghodsi M, Treangen T, Pop M. Accurate and fast estimation of taxonomic profiles from metagenomic shotgun sequences. BMC Genomics. 2011;12 Suppl 2:S4. Epub 2011 Jul 27."
pub_dict["glimmermg"] = "Kelley DR, Liu B, Delcher AL, Pop M, Salzberg SL. Gene prediction with Glimmer for metagenomic sequences augmented by classification and clustering. Nucleic Acids Res. 2012 Jan;40(1):e9. Epub 2011 Nov 18."
pub_dict["sparseassembler"] = "Ye C, Ma ZS, Cannon CH, Pop M, Yu DW. Exploiting sparseness in de novo genome assembly. BMC Bioinformatics. 2012 Apr 19;13 Suppl 6:S1."
pub_dict["minimus"] = "Sommer DD, Delcher AL, Salzberg SL, Pop M. Minimus: a fast, lightweight genome assembler. BMC Bioinformatics. 2007 Feb 26;8:64."
pub_dict["velvet"] = "Zerbino DR, Birney E. Velvet: algorithms for de novo short read assembly using de Bruijn graphs. Genome Res. 2008 May;18(5):821-9. Epub 2008 Mar 18."
pub_dict["phylosift"] = "http://phylosift.wordpress.com/"
pub_dict["ca"] = "Miller JR, Delcher AL, Koren S, Venter E, Walenz BP, Brownley A, Johnson J, Li K, Mobarry C, Sutton G. Aggressive assembly of pyrosequencing reads with mates.Bioinformatics. 2008 Dec 15;24(24):2818-24. Epub 2008 Oct 24."
pub_dict["phmmer"] = "Eddy SR. Accelerated Profile HMM Searches. PLoS Comput Biol. 2011 Oct;7(10):e1002195. Epub 2011 Oct 20."
pub_dict["blast"] = "Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. J Mol Biol. 1990 Oct 5;215(3):403-10." 
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
        print "\n======Supported programs and citations (if available)=======\n"
        for type in supported_programs.keys():
            print "[" + type + "]"
            ccnt = 1
            for prog in supported_programs[type]:
                citation = "NA"
                try: 
                    citation = pub_dict[prog]
                except KeyError:
                    citation = "NA"
                print "  %d)"%(ccnt)+" "+progname_dict[prog]
                print "    "+citation+"\n"
                ccnt +=1
        sys.exit(0)
                
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
        mapper = a.lower()
        foundit = False
        for sm in supported_mappers:
            if mapper not in sm:
                continue
            else:
                mapper = sm
                foundit = True
                break
        if not foundit:
            print "!!Sorry, %s is not a supported read alignment method. Using bowtie instead"%(mapper)
            mapper = "bowtie"
    elif o in ("-r", "--retainBank"):
        retainBank = True
    elif o in ("-c", "--classifier"):
        cls = a.lower()
        foundit = False
        for sc in supported_classifiers:
            if cls not in sc:
                continue
            else:
                cls = sc
                foundit = True
                break
        if not foundit:
            print "!!Sorry, %s is not a supported classification method. Using FCP instead"%(fcp)
            orf = "fcp"

    elif o in ("-a","--assembler"):
        asm = a.lower()
        if asm == "metaidba":
            bowtie_mapping = 1
            
        foundit = False
        
        for sa in supported_assemblers:
            if asm not in sa:
                continue
            else:
                if asm != "velvet":
                    #some special cases required, velvet would trigger MetaVelvet, not velvet, etc
                    asm = sa
                foundit = True
                break
        
        if not foundit:
            print "!!Sorry, %s is not a supported assembler. Using SOAPdenovo instead"%(asm)
            asm = "soap"

        
    elif o in ("-g","--genecaller"):
        orf = a.lower()
        foundit = False
        for sg in supported_assemblers:
            if orf not in sg:
                continue
            else:
                orf = sg
                foundit = True
                break
        if not foundit:
            print "!!Sorry, %s is not a supported gene caller. Using FragGeneScan instead"%(orf)
            orf = "soap"

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
