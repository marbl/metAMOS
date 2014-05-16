#!python

## $Id$

#########################
## runPipeline.py - main pipeline driver for metAMOS
#########################
##first imports
import os,sys,traceback
sys.tracebacklimit = 0
shellv = os.environ["SHELL"]
## Setting up paths
INITIAL_SRC   = "%s%ssrc"%(sys.path[0], os.sep)
## Hardcode a k-mer size
DEFAULT_KMER  = "31"
## Hardcode a default taxonomic classification level
DEFAULT_TAXA_LEVEL = "class"
sys.path.insert(1, INITIAL_SRC)
validate_install = 0
if validate_install:
    import check_install
    rt = check_install.validate_dir(sys.path[0].strip(),sys.path[0]+os.sep+'required_file_list.txt')
    if rt == -1:
        print "MetAMOS not properly installed, please reinstall or contact development team for assistance"
        sys.exit(1)

import utils
import workflow
utils.configureEnvironment(utils.INITIAL_UTILS)

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
skipsteps = set()
skipsteps.add("FindRepeats")
isolate_genome = False
userKmerSupplied = False
asmScores = "%d"%(utils.SCORE_TYPE.LAP)
userScoresSupplied = False
asmScoreWeights = dict()

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
    print "   -t = <string>:   filter input reads? (default = %s, supported = %s)"%(selected_programs["preprocess"], ",".join(supported_programs["preprocess"]))
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
    print "[Validate]"
    print "   -X = <string>: comma-separated list of validators to run on the assembly. (default = %s, supported = %s)"%(selected_programs["validate"], ",".join(supported_programs["validate"]))
    print "   -S = <string>: comma-separated list of scores to use to select the winning assembly. By default, all validation tools specified by -X will be run. For each score, an optional weight can be specified as SCORE:WEIGHT. For example, LAP:1,CGAL:2 (supported = %s)"%(",".join(utils.SCORE_TYPE.reverse_mapping.values()).lower())
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
    configurationText.append("K-Mer:\t\t\t%s\n"%(utils.Settings.kmer))
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
        progs = set(selected_programs[type].split(","))
        configurationText.append("[" + type + "]\n")
        for prog in progs:
           # special case for validation, orf is a reference to selected orf finder, n50 is nothing
           if type == "validate" and prog == "n50":
              continue
           if type == "validate" and prog == "orf":
              prog = selected_programs["findorfs"]
           if type == "preprocess" and prog == "none":
              configurationText.append(prog.upper() + "\n\tN/A\n\n")
              continue
           if type == "preprocess" and prog == "metamos":
              configurationText.append("metAMOS built-in filtering\n\tN/A\n\n") 
              continue

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

shortOptions = "hM:IR:rjwbd:s:e:o:k:c:a:n:p:qt:f:vm:4g:iu1l:x:yz:LBVX:S:"
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
supported_preprocessors = ["none", "metamos", "eautils", "pbcr"]
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
supported_validators = ["reapr", "orf", "lap", "ale", "quast", "frcbam", "freebayes", "cgal", "n50"]
supported_fannotate = ["blast"]
supported_scaffolders = ["bambus2"]
supported_programs["preprocess"] = supported_preprocessors
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
selected_programs["preprocess"] = "metamos"
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
#turn on by default
forcesteps = set()

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
    if o in ("-V", "--version"):
       print "metAMOS Version %s"%(utils.getVersion())
       sys.exit()
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-d", "--projectdir"):
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
      print "Error: unknown workflow %s specified. Please choose one of %s."%(wfName, ",".join(availableWorkflows.keys()))
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
         opts = wfopts
         args = wfargs
   except getopt.GetoptError, err:
      # print help information and exit:
       print str(err) # will print something like "option -a not recognized"
       usage()
       sys.exit(2)

# finally reload any commands we had
pip = workflow.Workflow("pipeline", settings.rundir + os.sep)
pip.read()
if len(pip.commandList.strip()) > 0:
   try:
      wfopts, wfargs = getopt.getopt(pip.commandList.strip().split(), shortOptions, longOptions)
      wfopts.extend(opts)
      wfargs.extend(args)
      opts = wfopts
      args = wfargs
   except getopt.GetoptError, err:
      # print help information and exit:
       print str(err) # will print something like "option -a not recognized"
       usage()
       sys.exit(2)

# update the list of options we have
utils.updateConfigCommands(inifile, opts)

for o, a in opts:
    if o in ("-v","--verbose"):
        utils.Settings.VERBOSE = True
    elif o in ("-i","--indexbowtie"):
        bowtie_mapping = 1
    elif o in ("-B","--noblastdb"):
        noblastdb = True
        #skip Metaphyler
        #skipsteps.append("Abundance")
        skipsteps.add("FunctionalAnnotation")
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
        skipsteps.update(allsteps[:allsteps.index(startat)]) 
    elif o in ("-e","--endat"):
        endat = a
        if endat not in allsteps:
            print "cannot end at %s, step does not exist in pipeline"%(endat)
            print allsteps 
        skipsteps.update(allsteps[allsteps.index(endat)+1:])
    elif o in ("-o", "--minoverlap"):
        pass
    elif o in ("-k", "--kmersize"):
        utils.Settings.kmer = a
        userKmerSupplied=True
    elif o in ("-4", "--454"):
       selected_programs["assemble"] = "newbler,%s"%(selected_programs["assemble"])
       selected_programs["mapreads"] = "bowtie2"
    elif o in ("-f", "--forcesteps"):
        for step in a.split(","):
           forcesteps.add(step)
           skipsteps.discard(step)
    elif o in ("-n", "--skipsteps"):
        for step in a.split(","):
           if (step.lower() in "mapreads"):
              print "Warning: MapReads cannot currently be skipped, re-enabling step!"
           elif (step.lower() in "assemble"):
              print "Warning: Assemble step cannot be skipped when no contigs are specified, re-enabling step! Please specify contigs with -c to initPipeline to skip assembly." 
           else:
              skipsteps.add(step)
              forcesteps.discard(step)
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
        selected_programs["preprocess"] = a.lower()
        found = False
        for sa in supported_preprocessors:
           if selected_programs["preprocess"] in sa:
              found = True
              selected_programs["preprocess"] = sa
              break
        if not found:
           print "Warning: invalid preprocessor %s specified, supported: %s"%(a.lower(), ",".join(supported_preprocessors))
           selected_programs["preprocess"] = "metamos"
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
            skipsteps.add("Propagate")
            skipsteps.add("Classify")
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
        noAsm = False

        for assembler in assemblers:
           if assembler == "soap2":
              assembler="soapdenovo2"
           elif assembler == "soap":
              assembler="soapdenovo"
           elif assembler == "idba":
              assembler="idba-ud"

           foundit = False
           for sa in supported_assemblers:
              if assembler == "none":
                 selected_programs["assemble"] = "none"
                 asmSpecified = False
                 foundit = True
                 noAsm = True
                 break
              elif assembler not in sa:
                 continue
              else:
                 if (assembler != "velvet" and assembler != "soapdenovo" and assembler != "idba") or assembler == sa:
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
              if not noAsm or len(asmcontigs) == 0:
                 print "!!Sorry, no valid assembler specified. Using SOAPdenovo instead"
                 selected_programs["assemble"] = "soapdenovo"
           if not noAsm:
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
       asmScores = ""
       scores = a.lower().strip().split(",")
       for s in scores:
          sSplit = s.upper().strip().split(":")
          score = sSplit[0]
          found = False
          for sa in utils.SCORE_TYPE.reverse_mapping.keys():
             if score.upper() in utils.SCORE_TYPE.reverse_mapping[sa]:
                if asmScores != "":
                    asmScores += ","
                else:
                    asmScores = ""
                asmScores += "%d"%(sa)
                weight = utils.getDefaultWeight(sa)
                if len(sSplit) > 1:
                    weight = float(sSplit[1])
                asmScoreWeights[sa] = weight
                found = True
                break
          if not found:
             print "Warning: invalid score type %s specified, supported: %s"%(score, ",".join(utils.SCORE_TYPE.reverse_mapping.values()))

       if asmScores == "":
          asmScores = "%d"%(utils.SCORE_TYPE.LAP)
       else:
          userScoresSupplied = True

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

if (settings.noblastdb or noblastdb):
    print "**no blast DB directory available, disabling steps requiring BLAST DB"
    if (selected_programs["annotate"] == "blast"):
        print "Cannot run blast for classification (model files in DB dir). replacing with kraken!"
        selected_programs["annotate"] = "kraken"
    skipsteps.add("FunctionalAnnotation")
    nofcpblast = True

print "[Steps to be skipped]: ", skipsteps
#remove started & ok flags in Logs
if os.path.exists("%s%sLogs%s*.started"%(settings.rundir,os.sep,os.sep)):
    os.system("rm %s%sLogs%s*.started"%(settings.rundir,os.sep,os.sep))

if not isolate_genome:
  skipsteps.add("MultiAlign")
else:
  try:
      settings.doscaffolding = True
      selected_programs["multialign"]
      selected_programs["validate"] = ",".join(supported_programs["validate"])
      selected_programs["assemble"] = selected_programs["assemble"] + ",velvet"
      selected_programs["findorfs"] = "prokka"
      asmScores = "%d"%(utils.SCORE_TYPE.ALL)
      skipsteps.add("Scaffold")
      skipsteps.add("Propagate")
  except KeyError:
      skipsteps.add("MultiAlign")

# by default don't do functional annotate or scaffold orf finding
if "FunctionalAnnotation" not in forcesteps:
   skipsteps.add("FunctionalAnnotation")
if "FindScaffoldORFS" not in forcesteps:
   skipsteps.add("FindScaffoldORFS")

if len(asmcontigs) != 0 and not asmSpecified:
   selected_programs["assemble"] = "none"

if len(readlibs) > 1 and "metaidba" in selected_programs["assemble"]:
    print "ERROR: meta-IDBA only supports 1 library, please select different assembler or reduce libraries"
    sys.exit(1)

if "scaffold" in skipsteps or "Scaffold" in skipsteps:
   skipsteps.add("Propagate")

#if we have pacbio reads use bowtie 2 (since bowtie 1 expects 1024bp sequences) and run CA
if selected_programs["preprocess"] == "pbcr":
   selected_programs["mapreads"] = "bowtie2"
   selected_programs["assemble"] = selected_programs["assemble"] + ",ca"

if userKmerSupplied == False and isolate_genome:
   always_run_programs.append("kmergenie")

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

if "Preprocess" in forcesteps:
   for path in readpaths:
      utils.run_process(settings, "touch %s"%(path),"RunPipeline")
   os.system("rm %s%sLogs%s*.ok"%(settings.rundir,os.sep,os.sep))
utils.Settings.readpaths = readpaths

asmfiles = []

for lib in readlibs:
    if "Assemble" in forcesteps:
        utils.run_process(settings, \
           "touch %s/Preprocess/out/lib%d.seq"%(settings.rundir,lib.id),\
           "RunPipeline")
    asmfiles.append("%s/Preprocess/out/lib%d.seq"%(settings.rundir,lib.id))

if "MapReads" in forcesteps:
    utils.run_process(settings, "rm %s/Assemble/out/*.contig.cvg"%(settings.rundir), "RunPipeline")
    utils.run_process(settings, "rm %s/Logs/mapreads.ok"%(settings.rundir), "RunPipeline")

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
          "rm %s/Logs/scaffold.ok"%(settings.rundir),\
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

if "Validate" in forcesteps:
   utils.run_process(settings, \
          "rm %s/Logs/validate.ok"%(settings.rundir), "RunPipeline")
   utils.run_process(settings, \
          "rm %s/Validate/out/%s.asm.selected"%(settings.rundir, settings.PREFIX), "RunPipeline")

if "FINDORFS" in forcesteps or "findorfs" in forcesteps or "FindORFS" in forcesteps:
   utils.run_process(settings, \
          "rm %s/Assemble/out/*.faa"%(settings.rundir),"RunPipeline")
   utils.run_process(settings, \
          "rm %s/Assemble/out/*.fna"%(settings.rundir),"RunPipeline")
   utils.run_process(settings, \
          "rm %s/FindORFS/out/*.faa"%(settings.rundir),"RunPipeline")
   utils.run_process(settings, \
          "rm %s/FindORFS/out/*.fna"%(settings.rundir),"RunPipeline")

if "Assemble" not in skipsteps and "Assemble" in forcesteps:
    utils.run_process(settings, \
          "rm %s/Logs/assemble.ok"%(settings.rundir),\
          "RunPipeline")
    utils.run_process(settings, \
          "rm %s/Assemble/out/*.asm.contig"%(settings.rundir),\
           "RunPipeline")

if "Classify" not in skipsteps and "Classify" in forcesteps:
    utils.run_process(settings, \
          "rm %s/Logs/classify.ok"%(settings.rundir), \
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
            skipsteps.add("FunctionalAnnotation")
        else:
            print utils.OK_GREEN+"\t*ok"+utils.ENDC

    if settings.nopysam:
       #need pysam for bowtie2 support
       supported_mappers = ["bowtie"]
       if "bowtie2" == selected_programs["mapreads"]:
          selected_programs["mapreads"] = "bowtie"

    if userScoresSupplied == False:
       for validator in selected_programs["validate"].split(","):
          validator = validator.upper()
          if (validator not in utils.SCORE_TYPE.mapping):
             continue
          if str(utils.SCORE_TYPE.mapping[validator]) not in asmScores:
             asmScores = "%s,%s"%(asmScores, utils.SCORE_TYPE.mapping[validator])

    for asmScore in asmScores.split(","):
       if int(asmScore) == utils.SCORE_TYPE.ALL:
          selected_programs["validate"] = ",".join(supported_programs["validate"])
       elif utils.SCORE_TYPE.reverse_mapping[int(asmScore)].lower() not in selected_programs["validate"]:
          selected_programs["validate"] = "%s,%s"%(selected_programs["validate"], utils.SCORE_TYPE.reverse_mapping[int(asmScore)].lower())

    # intialize weights
    utils.initValidationScores(asmScoreWeights)

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
    preprocess.init(readlibs, asmcontigs, skipsteps, selected_programs["assemble"], run_fastqc,selected_programs["preprocess"])
    assemble.init(readlibs, skipsteps, selected_programs["assemble"], asmcontigs, (userKmerSupplied == False and isolate_genome))
    mapreads.init(readlibs, skipsteps, selected_programs["mapreads"], savebtidx,ctgbpcov,lowmem)
    validate.init(readlibs, skipsteps, selected_programs["validate"], asmScores)
    findorfs.init(readlibs, skipsteps, selected_programs["findorfs"], min_ctg_len, min_ctg_cvg,read_orfs)
    findreps.init(readlibs, skipsteps)
    multialign.init(readlibs, skipsteps, forcesteps, selected_programs["multialign"],refgenomes)
    annotate.init(readlibs, skipsteps, selected_programs["annotate"], nofcpblast)
    fannotate.init(skipsteps)
    abundance.init(readlibs, skipsteps, forcesteps, selected_programs["annotate"])
    scaffold.init(readlibs, skipsteps, retainBank)
    findscforfs.init(readlibs, skipsteps, selected_programs["findorfs"])
    propagate.init(readlibs, skipsteps, selected_programs["annotate"])
    classify.init(readlibs, skipsteps, selected_programs["annotate"], lowmem, 0 if not isolate_genome else 100)
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
       print "Please see results in %s/Postprocess/out"%(settings.rundir)
       html_results = "%s/Postprocess/out/html/summary.html"%(settings.rundir)
       if os.path.exists(html_results):
          print "metAMOS summary report is file://%s"%(html_results)

    except Exception, e:
       print "Oops, MetAMOS finished with errors! see text in red above for details."
       if utils.Settings.VERBOSE:
          exc_type, exc_value, exc_traceback = sys.exc_info()
          traceback.print_exception(exc_type, exc_value, exc_traceback,
                              limit=100, file=sys.stderr)
       sys.exit()
