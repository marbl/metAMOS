#!python

## $Id$

#########################
## runPipeline.py - main pipeline driver for metAMOS
#########################
##first imports
import os,sys
sys.tracebacklimit = 0
shellv = os.environ["SHELL"]
## Setting up paths
INITIAL_SRC   = "%s%ssrc"%(sys.path[0], os.sep)
## Hardcode a k-mer size
DEFAULT_KMER  = 31
## Hardcode a default taxonomic classification level
DEFAULT_TAXA_LEVEL = "class"
sys.path.append(INITIAL_SRC)
import check_install
validate_install = 0
if validate_install:
    rt = check_install.validate_dir(sys.path[0].strip(),sys.path[0]+os.sep+'required_file_list.txt')
    if rt == -1:
        print "MetAMOS not properly installed, please reinstall or contact development team for assistance"
        sys.exit(1)

import utils
import workflow
ppath = ""
if "PYTHONPATH" not in os.environ:
   os.environ["PYTHONPATH"] = ""
else:
   ppath = os.environ["PYTHONPATH"] 
   #os.environ["PYTHONPATH"] = ""
os.environ["PYTHONPATH"]+=utils.INITIAL_UTILS+os.sep+"python"+os.pathsep
os.environ["PYTHONPATH"]+=utils.INITIAL_UTILS+os.sep+"ruffus"+os.pathsep
os.environ["PYTHONPATH"] += utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.pathsep
os.environ["PYTHONPATH"] += utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"+os.pathsep
os.environ["PYTHONPATH"] += utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.pathsep
os.environ["PYTHONPATH"] += utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.sep+"python"+os.pathsep
os.environ["PYTHONPATH"] += utils.INITIAL_UTILS+os.pathsep

if "PERL5LIB" not in os.environ:
    os.environ["PERL5LIB"] =  INITIAL_SRC+os.sep+"phylosift"+os.sep+"lib"+os.sep
else:
    os.environ["PERL5LIB"] =  INITIAL_SRC+os.sep+"phylosift"+os.sep+"lib"+os.sep + os.pathsep + os.environ["PERL5LIB"]
try:
    os.environ["PYTHONPATH"] += sys._MEIPASS + os.pathsep
    os.environ["PYTHONHOME"] = sys._MEIPASS + os.pathsep
except Exception:
    pass

try:
    sys._MEIPASS
    #if we are here, frozen binary
except Exception:
    #else normal mode, add site dir
    import site
    site.addsitedir(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python")
    site.addsitedir(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.sep+"python")

sys.path.append(utils.INITIAL_UTILS)
sys.path.append(utils.INITIAL_UTILS+os.sep+"python")
sys.path.append(utils.INITIAL_UTILS+os.sep+"ruffus")
sys.path.append(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python")
sys.path.append(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.sep+"python")
try:
    sys.path.append(sys._MEIPASS)
except Exception:
    pass
sys.path.append("/usr/lib/python")

#remove imports from pth file, if exists
nf = []
if 'bash' in shellv or utils.cmdExists('export'):
   os.system("export PYTHONPATH=%s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"))
   os.system("export PYTHONPATH=%s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"))
elif utils.cmdExists('setenv'):
   os.system("setenv PYTHONPATH %s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"))
   os.system("setenv PYTHONPATH %s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"))
else:
   print "Warning: could not set PYTHONPATH. Unknown shell %s, some functionality may not work\n"%(shellv)
## The usual library dependencies
import string
import time
import datetime
import BaseHTTPServer
import getopt
import re
import subprocess
import webbrowser
import multiprocessing
from operator import itemgetter
from ruffus import *
from task import JobSignalledBreak
skipsteps = ["FindRepeats"]
isolate_genome = False
asmScore = utils.SCORE_TYPE.LAP

## Get start time
t1 = time.time()

def usage():
    print "usage: runPipeline [options] -d projectdir"
    print "   -h = <bool>:   print help [this message]"
    print "   -j = <bool>:   just output all of the programs and citations then exit (default = NO)"
    print "   -v = <bool>:   verbose output? (default = NO)"
    print "   -d = <string>: directory created by initPipeline (REQUIRED)"

    print "\n[options]: [pipeline_opts] [misc_opts]"
    print "\n[pipeline_opts]: options that affect the pipeline execution"

    print "Pipeline consists of the following steps:"
    print "  Preprocess, Assemble, FindORFS, MapReads, Abundance, Annotate,"
    print "  FunctionalAnnotation, Scaffold, Propagate, Classify, Postprocess"

    print "Each of these steps can be referred to by the following options:" 
    print "   -f = <string>: force this step to be run (default = NONE)"
    print "   -s = <string>: start at this step in the pipeline (default = Preprocess)"
    print "   -e = <string>: end at this step in the pipeline (default = Postprocess)"
    print "   -n = <string>: step to skip in pipeline (default=NONE)"

    print "\nFor each step you can fine-tune the execution as follows"
    print "[Preprocess]"
    print "   -t = <bool>:   filter input reads? (default = NO)"
    print "   -q = <bool>:   produce FastQC quality report for reads with quality information (fastq or sff)? (default = NO)"
    print "[Assemble]"
    print "   -a = <string>: genome assembler to use (default = %s, supported = %s)"%(selected_programs["assemble"], ",".join(supported_programs["assemble"]))
    print "   -k = <kmer size>: k-mer size to be used for assembly (default = " + str(DEFAULT_KMER) +  ")"
    print "   -o = <int>>:   min overlap length"
    print "[MapReads]"
    print "   -m = <string>: read mapper to use? (default = %s, supported = %s)"%(selected_programs["mapreads"], ",".join(supported_programs["mapreads"]))
    print "   -i = <bool>:   save bowtie (i)ndex? (default = NO)"
    print "   -b = <bool>:   create library specific per bp coverage of assembled contigs (default = NO)"
    print "[FindORFS]"
    print "   -g = <string>: gene caller to use (default = %s, supported = %s)"%(selected_programs["findorfs"], ",".join(supported_programs["findorfs"]))
    print "   -l = <int>:    min contig length to use for ORF call (default = 300)"
    print "   -x = <int>>:   min contig coverage to use for ORF call (default = 3X)"
    print "[Annotate]"
    print "   -c = <string>: classifier to use for annotation (default = %s, supported = %s"%(selected_programs["annotate"], ",".join(supported_programs["annotate"]))
    print "   -u = <bool>:   annotate unassembled reads? (default = NO)"

    print "[Classify]"
    print "   -z = <string>: taxonomic level to categorize at (default = %s)"%(DEFAULT_TAXA_LEVEL)

    print "\n[misc_opts]: Miscellaneous options"
    print "   -B = <bool>:   blast DBs not available (default = NO)"
    print "   -r = <bool>:   retain the AMOS bank?  (default = NO)"
    print "   -p = <int>:    number of threads to use (be greedy!) (default=1)"
    print "   -4 = <bool>:   454 data? (default = NO)"    
    print "   -L = <bool>:   generate local Krona plots. Local Krona plots can only be viewed on the machine they are generated on but will work on a system with no internet connection (default = NO)"

def updateCounter():
    ##if user says its ok, create MetAMOS counter git repo and update counter each time its run!
    if os.path.exists("%s/config/usage.ok"%(utils.Settings.METAMOS_UTILS)):
        FORMAT = '%Y%m%d%H%M%S'
        filestamp = datetime.datetime.now().strftime(FORMAT)
        #print filestamp, utils.Settings.METAMOS_UTILS
        #"%s/pipeline.run"%(settings.rundir)) 
        os.system("cat %s/pipeline.run > %s/%s.txt"%(utils.Settings.rundir,utils.Settings.rundir,filestamp))
        #print "%s/config/updateftpcounter.sh %s/%s.txt %s.txt"%(utils.Settings.METAMOS_UTILS,utils.Settings.rundir,filestamp,filestamp)
        os.system("%s/config/updateftpcounter.sh %s/%s.txt %s.txt >& ftp.out"%(utils.Settings.METAMOS_UTILS,utils.Settings.rundir,filestamp,filestamp))

def printConfiguration(fileName=None):
    configurationText = []
    configurationText.append("metAMOS configuration summary:\n")
    configurationText.append("metAMOS Version:\t%s\n"%(utils.getVersion()))
    configurationText.append("Time and Date:\t\t%s\n"%(str(datetime.date.today())))
    configurationText.append("Working directory:\t%s\n"%(utils.Settings.rundir))
    configurationText.append("Prefix:\t\t\t%s\n"%(utils.Settings.PREFIX))
    configurationText.append("K-Mer:\t\t\t%d\n"%(utils.Settings.kmer))
    configurationText.append("Threads:\t\t%d\n"%(utils.Settings.threads)) 
    configurationText.append("Taxonomic level:\t%s\n"%(utils.Settings.taxa_level))
    configurationText.append("Verbose:\t\t%s\n"%(utils.Settings.VERBOSE))
    configurationText.append("Steps to skip:\t\t%s\n"%(", ".join(skipsteps)))
    configurationText.append("Steps to force:\t\t%s\n"%(", ".join(forcesteps)))

    # get metamos citations
    configurationText.append("\n")
    configurationText.append("[citation]\n")
    (progName, citation) = utils.getProgramCitations(settings, "metamos")
    configurationText.append(progName + "\n")
    configurationText.append("\t" + citation + "\n\n")
    if isolate_genome == True:
       (progName, citation) = utils.getProgramCitations(settings, "metamos_isolate")
       configurationText.append(progName + "\n")
       configurationText.append("\t" + citation + "\n\n") 

    configurationText.append("\n")
    configurationText.append("Step-specific configuration:\n")
    for type in selected_programs.keys():
        configurationText.append("[" + type + "]\n")
        progs = set(selected_programs[type].split(","))
        for prog in progs:
           if prog == None or prog == "none":
              configurationText.append("None\n\n")
           else:
              (progName, citation) = utils.getProgramCitations(settings, prog)
              if progName == "":
                 configurationText.append(prog + "\n")
              else:
                 configurationText.append(progName + "\n")
              try:
                 configurationText.append("\t" + eval("utils.Settings.%s"%(prog.replace("-", "_").upper()))+"\n")
              except AttributeError:
                 try:
                    configurationText.append("\t" + generic.getLocation(utils.STEP_NAMES.mapping[type.upper()], prog) + "\n")
                 except KeyError:
                     configurationText.append("UNKNOWN")
              if citation != "":
                 configurationText.append("\t" + citation + "\n\n")

    # add step-indepent citations that are always run
    configurationText.append("[other]\n")
    for prog in always_run_programs:
       if prog == None or prog == "none":
          configurationText.append("None\n\n")
       else:
          (progName, citation) = utils.getProgramCitations(settings, prog)
          configurationText.append(progName + "\n")
          try:
             configurationText.append("\t" + eval("utils.Settings.%s"%(prog.upper()))+"\n")
          except AttributeError:
             configurationText.append("\tUNKNOWN\n")
          configurationText.append("\t" + citation + "\n\n")

    if fileName == "" or fileName == None:
        print ''.join(configurationText)
    else:
        conf = open(fileName, 'w')
        conf.write(''.join(configurationText))
        conf.close()

shortOptions = "hM:IR:rjwbd:s:e:o:k:c:a:n:p:qtf:vm:4g:iu1l:x:yz:LBVX:S:"
longOptions = ["help", \
                                        "multialigner",\
                                        "isolate",\
                                        "refgenomes",\
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
                                        "lowmem",\
                                        "minlen",\
                                        "mincov", \
                                        "justprogs", \
                                        "what", \
                                        "lowcpu",\
                                        "taxalevel",\
                                        "localKrona",\
                                        "noblastdb",\
                                        "version",\
                                        "validator",\
                                        "asmscore"]
try:
    opts, args = getopt.getopt(sys.argv[1:], shortOptions, longOptions)
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

## always use long names, search will auto-detect abbreviations
settings = utils.Settings(DEFAULT_KMER, multiprocessing.cpu_count() - 1, "", DEFAULT_TAXA_LEVEL)
import generic

supported_programs = {}
supported_genecallers = ["fraggenescan","metagenemark","glimmermg"]
supported_assemblers = ["newbler", "soapdenovo","soapdenovo2","ca","velvet","velvet-sc","metavelvet",\
                            "metaidba","sparseassembler","minimus"]
supported_assemblers.extend(generic.getSupportedList(utils.INITIAL_UTILS, utils.STEP_NAMES.ASSEMBLE))

supported_mappers = ["bowtie","bowtie2"]
supported_abundance = ["metaphyler"]
supported_aligners = ["mgcat"]
supported_classifiers = ["fcp","phylosift","phmmer","blast",\
                             "metaphyler", "phymm"]
supported_classifiers.extend(generic.getSupportedList(utils.INITIAL_UTILS, utils.STEP_NAMES.ANNOTATE))
supported_validators = ["lap", "ale", "quast", "frcbam", "freebayes", "cgal"]
supported_fannotate = ["blast"]
supported_scaffolders = ["bambus2"]
supported_programs["findorfs"] = supported_genecallers
supported_programs["assemble"] = supported_assemblers
supported_programs["mapreads"] = supported_mappers
supported_programs["abundance"] = supported_abundance
supported_programs["annotate"] = supported_classifiers
supported_programs["fannotate"] = supported_fannotate
supported_programs["scaffold"] = supported_scaffolders
supported_programs["multialign"] = supported_aligners
supported_programs["validate"] = supported_validators

supported_taxonomic = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

selected_programs = {}
selected_programs["assemble"] = "soapdenovo"
selected_programs["findorfs"] = "fraggenescan"
selected_programs["mapreads"] = "bowtie"
selected_programs["abundance"] = "metaphyler"
selected_programs["annotate"] = "kraken"
selected_programs["fannotate"] = "blast"
selected_programs["scaffold"] = "bambus2"
selected_programs["multialign"] = "mgcat"
selected_programs["validate"] = "lap"

always_run_programs = ["krona"]


allsteps = ["Preprocess","Assemble","MapReads","MultiAlign","FindORFS","FindRepeats","Abundance","Annotate",\
                "FunctionalAnnotation","Scaffold","FindScaffoldORFS","Propagate","Classify","Postprocess"]

## Need comments here and further down

output = None
reads = None
quals = None
format = None
savebtidx = False
verbose = False
bowtie_mapping = 1
startat = None
endat = None
filter = False
forcesteps = []

run_fastqc = False
runfast = False
retainBank = False

fff = ""
readlen = 75
fqlibs = {}
fqfrags = []
rlibs = []
ctgbpcov = False
min_ctg_len = 300
read_orfs = False
min_ctg_cvg = 3
lowmem= False
output_programs = 0
settings = utils.Settings(DEFAULT_KMER, multiprocessing.cpu_count() - 1, "", DEFAULT_TAXA_LEVEL)
nofcpblast = False
noblastdb = False
asmSpecified = False
refgenomes = ""

# get the required rundirectory so we can get further options from the workflow
for o, a in opts:
    if o in ("-d", "--projectdir"):
        utils.Settings.rundir = os.path.abspath(a)
        if not os.path.exists(a):
          print "project dir %s does not exist!"%(settings.rundir)
          usage()
          sys.exit(1)
if not os.path.exists(settings.rundir) or settings.rundir == "":
    print "project dir %s does not exist!"%(settings.rundir)
    usage()
    sys.exit(1)

inifile = settings.rundir+os.sep+"pipeline.ini"
inf = open(inifile,'r')
(asmcontigs, readlibs, wfName) = utils.readConfigInfo(inf, "%s/Preprocess/in/"%(settings.rundir))

if wfName != "":
   #parse frag/libs out of pipeline.ini out of rundir
   availableWf = workflow.getSupportedWorkflows("%s/workflows"%(utils.INITIAL_UTILS), True)
   availableWf.extend(workflow.getSupportedWorkflows(os.getcwd(), True))
   availableWorkflows = dict()
   for wf in availableWf:
      availableWorkflows[wf.name] = wf

   if wfName.lower() not in availableWorkflows.keys():
      print "Error: unknown wofkflow %s specified. Please choose one of %s."%(wfName, ",".join(availableWorkflows.keys()))
      sys.exit(1)
   wf = availableWorkflows[wfName]
   try:
      wfopts, wfargs = getopt.getopt(wf.commandList.strip().split(), shortOptions, longOptions)
      if wf.canModify():
         wfopts.extend(opts)
         wfargs.extend(args)
         opts = wfopts
         args = wfargs
      else:
         for o, a in opts:
            if o in ("-V", "--version") or o in ("-h", "--help"):
               wfopts.append([o, a])
         opts = wfopts
         args = wfargs
   except getopt.GetoptError, err:
      # print help information and exit:
       print str(err) # will print something like "option -a not recognized"
       usage()
       sys.exit(2)

for o, a in opts:
    if o in ("-V", "--version"):
       print "metAMOS Version %s"%(utils.getVersion())
       sys.exit()
    elif o in ("-v","--verbose"):
        utils.Settings.VERBOSE = True
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-i","--indexbowtie"):
        bowtie_mapping = 1
    elif o in ("-B","--noblastdb"):
        noblastdb = True
        #skip Metaphyler
        #skipsteps.append("Abundance")
        skipsteps.append("FunctionalAnnotation")
        #skip
        nofcpblast = True
    elif o in ("-y","--lowcpu"):
        nofcpblast = True
    elif o in ("-w","--what"):
        utils.Settings.OUTPUT_ONLY = True
    elif o in ("-j","--justprogs"):
        output_programs = 1
        print "\n======Supported programs and citations (if available)=======\n"
        for type in supported_programs.keys():
            print "[" + type + "]"
            ccnt = 1
            for prog in supported_programs[type]:
                (progName, citation) = utils.getProgramCitations(settings, prog)
                #citation = "NA"
                try: 
                    citation = pub_dict[prog]
                except KeyError:
                    citation = "NA"
                print "  %d)"%(ccnt)+" "+progName
                print "    "+citation+"\n"
                ccnt +=1
        sys.exit(0)
                
    elif o in ("-u","--unassembledreads"):
        utils.Settings.annotate_unmapped = 1
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
        skipsteps.extend(allsteps[:allsteps.index(startat)]) 
    elif o in ("-e","--endat"):
        endat = a
        if endat not in allsteps:
            print "cannot end at %s, step does not exist in pipeline"%(endat)
            print allsteps 
        skipsteps.extend(allsteps[allsteps.index(endat)+1:])
    elif o in ("-o", "--minoverlap"):
        pass
    elif o in ("-k", "--kmersize"):
        utils.Settings.kmer = int(a)
    elif o in ("-4", "--454"):
       selected_programs["assemble"] = "newbler,%s"%(selected_programs["assemble"])
       selected_programs["mapreads"] = "bowtie2"
    elif o in ("-f", "--forcesteps"):
        forcesteps = a.split(",")
    elif o in ("-n", "--skipsteps"):
        skipsteps.extend(a.split(","))
    elif o in ("-p", "--threads"):
        if int(a) > 0:
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
    elif o in ("-I", "--isolate"):
        isolate_genome = True
    elif o in ("-R", "--refgenomes"):
        if not os.path.exists(a):
            print "ref genome dir %s does not exist!"%(a)
            usage()
            sys.exit(1)
        refgenomes = a
    elif o in ("-M", "--multialigner"):
        selected_programs["multialign"] = a.lower()
        foundit = False
        for sm in supported_aligners:
            if selected_programs["multialign"] not in sm:
                continue
            else:
                selected_programs["multialign"] = sm
                foundit = True
                break
        if not foundit:
            print "!!Sorry, %s is not a supported multi alignment method. Using mgcat instead"%(selected_programs["multialign"])
            selected_programs["multialign"] = "mgcat"
    elif o in ("-m", "--mapper"):
        selected_programs["mapreads"] = a.lower()
        foundit = False
        for sm in supported_mappers:
            if selected_programs["mapreads"] not in sm:
                continue
            else:
                selected_programs["mapreads"] = sm
                foundit = True
                break
        if not foundit:
            print "!!Sorry, %s is not a supported read alignment method. Using bowtie instead"%(selected_programs["mapreads"])
            selected_programs["mapreads"] = "bowtie"
        #mapper = a
    elif o in ("-r", "--retainBank"):
        retainBank = True
    elif o in ("-c", "--classifier"):
        selected_programs["annotate"] = a.lower()
        foundit = False
        for sc in supported_classifiers:
            if selected_programs["annotate"] not in sc:
                continue
            else:
                selected_programs["annotate"] = sc
                foundit = True
                break
        if sc == "metaphyler":
            #not quite ready for primetime, need krona import script and annots file
            skipsteps.append("Propagate")
            skipsteps.append("Classify")
        if not foundit:
            print "!!Sorry, %s is not a supported classification method. Using FCP instead"%(selected_programs["annotate"])
            selected_programs["annotate"] = "fcp"
    elif o in ("-z", "--taxalevel"):
        utils.Settings.taxa_level = a.lower()

        if utils.Settings.taxa_level not in supported_taxonomic:
           print "!!Sorry, %s is not a valid taxonomic level. Using class instead"%(utils.Settings.taxa_level)

    elif o in ("-1","--lowmem"):
        lowmem = True
    elif o in ("-a","--assembler"):
        if "metaidba" in a.lower():
            bowtie_mapping = 1
            
        assemblers = a.lower().strip().split(",")
        selected_programs["assemble"] = None 

        for assembler in assemblers:
           if assembler == "soap2":
              assembler="soapdenovo2"
           elif assembler == "soap":
              assembler="soapdenovo"

           foundit = False
           for sa in supported_assemblers:
              if assembler not in sa:
                 continue
              else:
                 if (assembler != "velvet" and assembler != "soapdenovo") or assembler == sa:
                    #some special cases required, velvet would trigger MetaVelvet, not velvet, etc
                    if selected_programs["assemble"] != None:
                       selected_programs["assemble"] += ","
                    else:
                       selected_programs["assemble"] = ""
                    selected_programs["assemble"] += sa
                 foundit = True
                 break
        
           if not foundit:
               print "!!Sorry, %s is not a supported assembler."%(assembler)

        if selected_programs["assemble"] == None:
           print "!!Sorry, no valid assembler specified. Using SOAPdenovo instead"
           selected_programs["assemble"] = "soapdenovo"

        asmSpecified = True
        
    elif o in ("-g","--genecaller"):
        selected_programs["findorfs"] = a.lower()
        foundit = False
        for sg in supported_genecallers:
            if selected_programs["findorfs"] not in sg:
                continue
            else:
                selected_programs["findorfs"] = sg
                foundit = True
                break
        if selected_programs["findorfs"] == "metagenemark" and "Darwin" in settings.OSTYPE:
            print "!!Sorry, %s is not a supported gene caller for Mac OSX. Using FragGeneScan instead"%(selected_programs["findorfs"])
            selected_programs["findorfs"] = "fraggenescan"
        if not foundit:
            print "!!Sorry, %s is not a supported gene caller. Using FragGeneScan instead"%(selected_programs["findorfs"])
            selected_programs["findorfs"] = "fraggenescan"

    elif o in ("-S", "--asmscore"):
       found = False
       for sa in utils.SCORE_TYPE.reverse_mapping.keys():
          print "Comparing %s to %s"%(a.upper(), utils.SCORE_TYPE.reverse_mapping[sa])
          if a.upper() in utils.SCORE_TYPE.reverse_mapping[sa]:
             asmScore = sa
             found = True
             break
       if not found:   
          print "Warning: invalid score type %s specified, supported: %s"%(a, ",".join(utils.SCORE_TYPE.reverse_mapping.values()))
          asmScore = utils.SCORE_TYPE.LAP
    elif o in ("-X", "--validator"):
        validators = a.lower().split(",")
        selected_programs["validate"] = None

        for validator in validators:
           if "snp" in validator:
              validator = "freebayes"
           foundit = False
           for sa in supported_validators:
              if validator not in sa:
                 continue
              else:
                 if selected_programs["validate"] != None:
                    selected_programs["validate"] += ","
                 else:
                    selected_programs["validate"] = ""
                 selected_programs["validate"] += sa
                 foundit = True
                 break

           if not foundit:
               print "!!Sorry, %s is not a supported validator."%(validator)

        if selected_programs["validate"] == None:
           print "!!Sorry, no valid assembly validator specified. Using LAP instead"
           selected_programs["validate"] = "lap"

    elif o in ("-f","--fastest"):
        runfast = True
    elif o in ("-b","--savebowtieidx"):
        savebtidx = True
    elif o in ("-L", "--localKrona"):
        utils.Settings.local_krona = True
    else:
        assert False, "unhandled option"

if (settings.noblastdb or noblastdb) and (selected_programs["annotate"] == "blast" or selected_programs["annotate"] == "fcp"):
    print "**no DB directory available, cannot run blast or FCP for classification (model files in DB dir). replacing with phylosift!"
    selected_programs["annotate"] = "phylosift"

print "[Steps to be skipped]: ", skipsteps
#remove started & ok flags in Logs
if os.path.exists("%s%sLogs%s*.started"%(settings.rundir,os.sep,os.sep)):
    os.system("rm %s%sLogs%s*.started"%(settings.rundir,os.sep,os.sep))

if not isolate_genome:
  skipsteps.append("MultiAlign")
else:
  try:
      settings.doscaffolding = True
      selected_programs["multialign"]
      selected_programs["validate"] = ",".join(supported_programs["validate"])
      selected_programs["assemble"] = selected_programs["assemble"] + ",velvet"
      selected_programs["findorfs"] = "prokka"
      asmScore = utils.SCORE_TYPE.ALL
      skipsteps.append("Scaffold")
      skipsteps.append("Propagate")
  except KeyError:
      skipsteps.append("MultiAlign")

if len(asmcontigs) != 0 and not asmSpecified:
   selected_programs["assemble"] = "none"

if len(readlibs) > 1 and "metaidba" in selected_programs["assemble"]:
    print "ERROR: meta-IDBA only supports 1 library, please select different assembler or reduce libraries"
    sys.exit(1)

if "scaffold" in skipsteps or "Scaffold" in skipsteps:
   skipsteps.append("Propagate")

inf.close()

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
      print "Processing lib %s read %s"%(lib.id, read.fname)

if "Preprocess" in forcesteps:
   for path in readpaths:
      utils.run_process(settings, "touch %s"%(path),"RunPipeline")
   os.system("rm %s%sLogs%s*.ok"%(settings.rundir,os.sep,os.sep))
utils.Settings.readpaths = readpaths

asmfiles = []

for lib in readlibs:
    if "MapReads" in forcesteps:
     if "MapReads" in forcesteps:
        for a in selected_programs["assemble"].strip().split(","):
           if os.path.exists("%s/Assemble/out/%s.asm.contig"%(settings.rundir,a)):
              utils.run_process(settings, \
              "touch %s/Assemble/out/%s.asm.contig"%(settings.rundir,a),\
              "RunPipeline")
        for a in asmcontigs:
           if os.path.exists("%s/Assemble/out/%s.asm.contig"%(settings.rundir,os.path.splitext(a)[0])):
                utils.run_process(settings, \
                   "touch %s/Assemble/out/%s.asm.contig"%(settings.rundir,os.path.splitext(a)[0]),\
                   "RunPipeline")
    if "Assemble" in forcesteps:
        utils.run_process(settings, \
           "touch %s/Preprocess/out/lib%d.seq"%(settings.rundir,lib.id),\
           "RunPipeline")
    asmfiles.append("%s/Preprocess/out/lib%d.seq"%(settings.rundir,lib.id))


utils.Settings.asmfiles = asmfiles

if "Propagate" in forcesteps:
    utils.run_process(settings, "rm %s/Logs/propagate.ok"%(settings.rundir), "RunPipeline")
    utils.run_process(settings, "touch %s/Annotate/out/%s.annots"%(settings.rundir, settings.PREFIX),\
                      "RunPipeline")

if "FindScaffoldORFS" in forcesteps:
    utils.run_process(settings, \
          "touch %s/Scaffold/out/%s.linearize.scaffolds.final"%(settings.rundir,settings.PREFIX),\
          "RunPipeline")

if "Scaffold" in forcesteps:
    utils.run_process(settings, \
          "rm %s/Scaffold/out/%s.linearize.scaffolds.final"%(settings.rundir,settings.PREFIX),\
          "RunPipeline")

if "Abundance" in forcesteps:
   utils.run_process(settings, \
          "touch %s/FindORFS/out/%s.faa"%(settings.rundir,settings.PREFIX),\
          "RunPipeline")
   utils.run_process(settings, \
          "rm %s/Abundance/out/%s.taxprof.pct.txt"%(settings.rundir,settings.PREFIX),\
          "RunPipeline")

if "Annotate" in forcesteps:
   utils.run_process(settings, \
          "rm %s/Annotate/out/%s.hits"%(settings.rundir,settings.PREFIX),"RunPipeline")

if "FINDORFS" in forcesteps or "findorfs" in forcesteps or "FindORFS" in forcesteps:
   utils.run_process(settings, \
          "rm %s/FindORFS/out/%s.faa"%(settings.rundir,settings.PREFIX),"RunPipeline")
   utils.run_process(settings, \
          "rm %s/FindORFS/out/%s.fna"%(settings.rundir,settings.PREFIX),"RunPipeline")

if "Assemble" not in skipsteps and "Assemble" in forcesteps:
    utils.run_process(settings, \
          "rm %s/Logs/assemble.ok"%(settings.rundir),\
          "RunPipeline")

if __name__ == "__main__":
    print "Starting metAMOS pipeline"
    if settings.threads < 1:
        settings.threads = 1
    settings = utils.initConfig(settings.kmer, settings.threads, settings.rundir, settings.taxa_level, settings.local_krona, settings.annotate_unmapped, settings.doscaffolding, settings.VERBOSE, settings.OUTPUT_ONLY)
    # add krona to system path
    currPath = os.environ["PATH"]
    if utils.Settings.KRONA not in currPath:
       os.environ["PATH"]="%s:%s"%(utils.Settings.KRONA, currPath)

    # check for memory/cpu
    if not settings.nopsutil and lowmem == False:
        import psutil

        avram = utils.getAvailableMemory(settings)
        print "[Available RAM: %d GB]"%(avram)
        lowmem= False
        if avram <= 32:
            print utils.WARNING_YELLOW+"\tThere is *%d GB of RAM currently available on this machine, suggested minimum of 32 GB"%(avram)+utils.ENDC
            print utils.WARNING_YELLOW+"\t*Enabling low MEM mode, might slow down some steps in pipeline"+utils.ENDC
            print utils.WARNING_YELLOW+"\t*If more RAM is available than what is listed above, please close down other programs and restart runPipeline"+utils.ENDC
            lowmem= True
        else:
            print utils.OK_GREEN+"\t*ok"+utils.ENDC
    if not settings.nopsutil and nofcpblast == False:
        numcpus = psutil.NUM_CPUS
        print "[Available CPUs: %d]"%(numcpus)
        if numcpus < 8:
            print utils.WARNING_YELLOW+"\t*Only %d CPU cores available, some modules might take awhile to complete"%(numcpus)+utils.ENDC
            print utils.WARNING_YELLOW+"\t*Disabling all BLAST (where possible)"+utils.ENDC
            nofcpblast = True
            skipsteps.append("FunctionalAnnotation")
        else:
            print utils.OK_GREEN+"\t*ok"+utils.ENDC

    if settings.nopysam:
       #need pysam for bowtie2 support
       supported_mappers = ["bowtie"]
       if "bowtie2" == selected_programs["mapreads"]:
          selected_programs["mapreads"] = "bowtie"

    if asmScore == utils.SCORE_TYPE.ALL:
       selected_programs["validate"] = ",".join(supported_programs["validate"])
    elif utils.SCORE_TYPE.reverse_mapping[asmScore].lower() not in selected_programs["validate"]:
       selected_programs["validate"] = "%s,%s"%(selected_programs["validate"], asmScore.lower())

    import preprocess
    import assemble
    import mapreads
    import validate
    import multialign
    import findorfs
    import findreps
    import abundance
    import annotate
    import fannotate
    import scaffold
    import findscforfs
    import propagate
    import classify
    import postprocess

    # initialize submodules
    preprocess.init(readlibs, asmcontigs, skipsteps, selected_programs["assemble"], run_fastqc,filter)
    assemble.init(readlibs, skipsteps, selected_programs["assemble"], asmcontigs)
    mapreads.init(readlibs, skipsteps, selected_programs["mapreads"], savebtidx,ctgbpcov,lowmem)
    validate.init(readlibs, skipsteps, selected_programs["validate"], asmScore)
    findorfs.init(readlibs, skipsteps, selected_programs["findorfs"], min_ctg_len, min_ctg_cvg,read_orfs)
    findreps.init(readlibs, skipsteps)
    multialign.init(readlibs, skipsteps, forcesteps, selected_programs["multialign"],refgenomes)
    annotate.init(readlibs, skipsteps, selected_programs["annotate"], nofcpblast)
    fannotate.init(skipsteps)
    abundance.init(readlibs, skipsteps, forcesteps, selected_programs["annotate"])
    scaffold.init(readlibs, skipsteps, retainBank)
    findscforfs.init(readlibs, skipsteps, selected_programs["findorfs"])
    propagate.init(readlibs, skipsteps, selected_programs["annotate"])
    classify.init(readlibs, skipsteps, selected_programs["annotate"], lowmem)
    postprocess.init(readlibs, skipsteps, selected_programs["annotate"])
    generic.init(skipsteps, readlibs)

    try:

       dlist = []
       #pipeline_printout(sys.stdout,[preprocess.Preprocess],verbose=1)                                                                                                                           
       tasks_to_run = ["preprocess.Preprocess"]

       if "ASSEMBLE" in skipsteps or "Assemble" in skipsteps or "assemble" in skipsteps or "asm" in skipsteps:
           pass
       else:
           tasks_to_run.append("assemble.Assemble")

       if "FINDORFS" in skipsteps or "FindORFS" in skipsteps or "findorfs" in skipsteps:
           pass
       else:
           tasks_to_run.append("findorfs.FindORFS")
       tasks_to_run.append("postprocess.Postprocess")
       #pipeline_printout(sys.stdout,tasks_to_run,verbose=2)                                                                                                                                      

       pipeline_printout(sys.stdout,[preprocess.Preprocess,assemble.Assemble, \
                         mapreads.MapReads, \
                         findorfs.FindORFS, findreps.FindRepeats, annotate.Annotate, \
                         abundance.Abundance, fannotate.FunctionalAnnotation, scaffold.Scaffold, \
                         findscforfs.FindScaffoldORFS, propagate.Propagate, \
                         classify.Classify, postprocess.Postprocess], verbose=1)

       if not utils.getFromPath("dot", "Graphviz") == "":
          pipeline_printout_graph (   'flowchart.svg',
                               'svg',
                               [postprocess.Postprocess],
                               no_key_legend = True)

       if 1 or not utils.Settings.BINARY_DIST:
           printConfiguration()
           printConfiguration("%s/pipeline.run"%(settings.rundir))                                                                                                                                   
       updateCounter()
       forcetasks = []
       for item in forcesteps:
           if "ASSEMBLE" in string.upper(item):
               forcetasks.append("assemble.Assemble")

           elif "FINDORFS" in string.upper(item):
               forcetasks.append("findorfs.FindORFS")

       #pipeline_run(tasks_to_run,forcedtorun_tasks=forcesteps,verbose=2)                                                                                                                         

       pipeline_run([preprocess.Preprocess, assemble.Assemble,findorfs.FindORFS, \
                    mapreads.MapReads, \
                    findreps.FindRepeats, annotate.Annotate, abundance.Abundance, \
                    fannotate.FunctionalAnnotation, scaffold.Scaffold, findscforfs.FindScaffoldORFS, \
                    propagate.Propagate, classify.Classify, postprocess.Postprocess],\
                    verbose = 1)

       #multiprocess threads
       t2 = time.time()
       elapsed = float(t2)-float(t1)

       #print elapsed
       print "done! pipeline took %.2f minutes"%(float(elapsed)/float(60.0))
    except JobSignalledBreak:
       print "Done with errors\n"
