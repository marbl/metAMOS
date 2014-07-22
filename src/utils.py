#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from datetime import date
from datetime import time
from datetime import datetime
from operator import itemgetter
import multiprocessing
import hashlib
shellv = os.environ["SHELL"]
_BINARY_DIST = False
def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS                                                                                                                           
        base_path = sys._MEIPASS
        _BINARY_DIST = True
        #print sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

application_path = ""
if getattr(sys, 'frozen', False):
    application_path = os.path.dirname(sys.executable)
elif __file__:
    application_path = os.path.dirname(__file__)

#print application_path

#NEED INSTALL DIR
#CWD print os.path.abspath(".")
#internal DIR print sys.path[0]

CSI="\x1B["
reset=CSI+"m"
OK_GREEN = CSI+'32m'
WARNING_YELLOW = CSI+'\033[93m'
ERROR_RED = CSI+'\033[91m'
ENDC = CSI+'0m'

_METAMOSDIR    = resource_path(sys.path[0])
INITIAL_UTILS = "%s%sUtilities"%(_METAMOSDIR, os.sep)
INITIAL_SRC   = "%s%ssrc"%(_METAMOSDIR, os.sep)
_NUM_LINES    = 10
_PROG_NAME_DICT = {}
_PUB_DICT = {}

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    reverse = dict((value, key) for key, value in enums.iteritems())
    mapping = dict((key, value) for key, value in enums.iteritems())
    enums['reverse_mapping'] = reverse
    enums['mapping'] = mapping
    return type('Enum', (), enums)

STEP_NAMES = enum("ASSEMBLE", "ANNOTATE", "SCAFFOLD")
STEP_OUTPUTS = enum(".asm.contig", ".hits", ".linearize.scaffolds.final")
INPUT_TYPE = enum("FASTQ", "FASTA", "CONTIGS", "SCAFFOLDS", "ORF_FA", "ORF_AA")

SCORE_TYPE = enum("ALL", "LAP", "ALE", "CGAL", "SNP", "FRCBAM", "ORF", "REAPR", "N50")
SCORE_WEIGHTS = dict()

_failFast = True

class AtomicCounter(object):
  def __init__(self, initval=0):
     self.val = multiprocessing.RawValue('i', initval)
     self.lock = multiprocessing.Lock()

  def increment(self):
     with self.lock:
         origVal = self.val.value
         self.val.value += 1
         return origVal

_atomicCounter = AtomicCounter(0)
_envCounter = AtomicCounter(0)

class Settings:
   asmfiles = []
   runfiles = []

   kmer = "55"
   threads = 16
   rundir = ""
   taxa_level = "class"
   local_krona = False
   annotate_unmapped = False
   task_dict = []
   noblastdb = False
   doscaffolding = False
   VERBOSE = False
   OUTPUT_ONLY = False

   PREFIX = ""

   OSTYPE = ""
   OSVERSION = ""
   MACHINETYPE = ""

   METAMOSDIR = ""
   METAMOS_UTILS = ""
   METAMOS_JAVA = ""

   FASTQC = ""
   AMOS = ""
   BAMBUS2 = ""

   SOAPDENOVO = ""
   SOAPDENOVO2 = ""
   METAIDBA = ""
   CA = ""
   BLASR = ""
   NEWBLER = ""
   VELVET = ""
   VELVET_SC = ""
   METAVELVET = ""
   SPARSEASSEMBLER = ""

   EAUTILS   = ""
   KMERGENIE = ""
   R         = ""

   MGCAT = ""

   METAPHYLER = ""
   BOWTIE = ""
   BOWTIE2 = ""
   SAMTOOLS = ""

   METAGENEMARK = ""
   FRAGGENESCAN = ""
   PROKKA = ""
   SIGNALP = ""

   FCP = ""
   PHMMER = ""
   PHYMM = ""
   BLAST = ""
   PHYLOSIFT = ""
   DB_DIR = ""
   BLASTDB_DIR = ""
   KRONA = ""
   REPEATOIRE = ""

   LAP = ""
   ALE = ""
   FRCBAM = ""
   FREEBAYES = ""
   CGAL = ""
   REAPR = ""
   QUAST = ""

   MPI = ""

   BINARY_DIST = 0

   nopsutil = False
   nopysam = False

   def __init__(self, kmer = None, threads = None, rundir = None, taxa_level = "", localKrona = False, annotateUnmapped = False, doScaffolding = False, verbose = False, outputOnly = False, update = False):

      configureEnvironment(INITIAL_UTILS)

      if (Settings.rundir != "" and update == False):
         return

      if (kmer == None or threads == None or rundir == None):
         print "Error settings is uninitialized and no intialization provided\n"
         raise(Exception)

      _BINARY_DIST = False
      try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS                                                                                                                           
        base_path = sys._MEIPASS
        _BINARY_DIST = True
        #print sys._MEIPASS
      except Exception:
        pass

      try:
         import pysam
         if verbose:
            print "Found pysam in %s"%(pysam.__file__)
      except ImportError:
         Settings.nopysam = True
         print "Could not import pysam, disabling."
      try:
         import psutil
         if verbose:
            print "Found psutil in %s"%(psutil.__file__)
      except ImportError:
         Settings.nopsutil = True
         print "Could not import psutil, disabling."

      Settings.rundir = rundir
      Settings.kmer = kmer
      Settings.threads = threads 
      Settings.rundir = rundir
      Settings.taxa_level = taxa_level
      Settings.local_krona = localKrona
      Settings.doscaffolding = doScaffolding
      Settings.annotate_unmapped = annotateUnmapped
      Settings.task_dict = []

      Settings.PREFIX = "proba"
      Settings.VERBOSE = verbose
      Settings.OUTPUT_ONLY = outputOnly

      Settings.OSTYPE        = "Linux"
      Settings.OSVERSION     = "0.0"
      Settings.MACHINETYPE   = "x86_64"
      getMachineType()

      Settings.METAMOSDIR    = sys.path[0]
      Settings.METAMOS_DOC   = "%s%sdoc"%(Settings.METAMOSDIR, os.sep)
      Settings.METAMOS_UTILS = "%s%sUtilities"%(Settings.METAMOSDIR, os.sep) 
      Settings.METAMOS_JAVA  = "%s%sjava:%s"%(Settings.METAMOS_UTILS,os.sep,os.curdir)

      Settings.noblastdb = False
      _DB_PATH = "%s/DB/"%(Settings.METAMOS_UTILS)
      _BLASTDB_PATH = _DB_PATH
      
      if _BINARY_DIST:
          #need to change KronaTools.pm to external Taxonomy directory
          try:
              _DB_PATH = "%s/DB/"%(application_path)
              _BLASTDB_PATH = _DB_PATH + os.sep + "blastdbs"+os.sep
              if "BLASTDB" in os.environ and len(os.environ["BLASTDB"]) != 0:
                  _BLASTDB_PATH == os.environ["BLASTDB"]
                  if not os.path.exists(_BLASTDB_PATH):
                      print "Error: cannot find BLAST DB directory, yet path set via $BLASTDB: %s. Disabling blastdb dependent programs"%(os.environ["BLASTDB"])
                      Settings.noblastdb = True
              elif not os.path.exists(_BLASTDB_PATH):
                 print "Error: cannot find BLAST DB directory, expected it in %s. Disabling blastdb dependent programs"%(_BLASTDB_PATH)
                 Settings.noblastdb = True
          except KeyError:
              #_DB_PATH = "./DB/"
              Settings.noblastdb = True
              pass

          if not os.path.exists(_DB_PATH):
              print "Error: cannot find DB directory in %s, was it deleted? oops, it is required to run MetAMOS!"%(_DB_PATH)
              sys.exit(1)
      elif Settings.rundir != "":
         if "BLASTDB" in os.environ and len(os.environ["BLASTDB"]) != 0:
             _BLASTDB_PATH == os.environ["BLASTDB"]
             if not os.path.exists(_BLASTDB_PATH):
                print "Error: cannot find BLAST DB directory, yet path set via $BLASTDB: %s. Disabling blastdb dependent programs"%(os.environ["BLASTDB"])
                Settings.noblastdb = True
         elif not os.path.exists("%s%srefseq_protein.00.pin"%(_BLASTDB_PATH, os.sep)):
            print "Error: cannot find BLAST DB directory, expected it in %s. Disabling blastdb dependent programs"%(_BLASTDB_PATH)
            Settings.noblastdb = True
       
      Settings.DB_DIR        = _DB_PATH 
      Settings.BLASTDB_DIR   = _BLASTDB_PATH 
      Settings.BINARY_DIST   = _BINARY_DIST
      Settings.AMOS          = "%s%sAMOS%sbin"%(Settings.METAMOSDIR, os.sep, os.sep)
      Settings.BAMBUS2       = Settings.AMOS

      Settings.SOAPDENOVO    = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.SOAPDENOVO2   = "%s%scpp%s%s-%ssoap2/bin"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.METAIDBA      = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.CA            = "%s%sCA%s%s-%s%sbin"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE.replace("x86_64", "amd64"), os.sep)
      Settings.NEWBLER       = "%s%snewbler%s%s-%s"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.VELVET        = "%s%scpp%s%s-%s%svelvet"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.VELVET_SC     = "%s%scpp%s%s-%s%svelvet-sc"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.METAVELVET    = "%s%scpp%s%s-%s%sMetaVelvet"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.SPARSEASSEMBLER = "%s%scpp%s%s-%s%sSparseAssembler"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.EAUTILS       = "%s%scpp%s%s-%s%seautils"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.KMERGENIE     = "%s%scpp%s%s-%s%skmergenie"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.R             = "%s%scpp%s%s-%s%sR"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.PHYMM = "%s%sperl%sphymm%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, os.sep)

      Settings.METAPHYLER    = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)

      Settings.BOWTIE        = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.BOWTIE2       = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.SAMTOOLS      = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)

      Settings.METAGENEMARK  = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.FRAGGENESCAN  = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.PROKKA        = "%s%scpp%s%s-%s/prokka/bin"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      runp = True
      if 1:
           try:
               kronalibf = open("%s%scpp%s%s-%s/prokka/bin/prokka"%(Settings.METAMOS_UTILS,os.sep,os.sep, Settings.OSTYPE, Settings.MACHINETYPE))
           except IOError:
               #this is initPipeline, skip                                                                                                                                                                                                    
               runp = False

      if _BINARY_DIST and runp:
           #need to change PROKKA to external db directory
           kronalibf = open("%s%scpp%s%s-%s/prokka/bin/prokka"%(Settings.METAMOS_UTILS,os.sep,os.sep, Settings.OSTYPE, Settings.MACHINETYPE))
           data = kronalibf.read()
           if "my $DBDIR = \"$FindBin::RealBin/../db\";" not in data:
               kronalibf.close()
           else:
               dd = data.replace("my $DBDIR = \"$FindBin::RealBin/../db\";","my $DBDIR = \"%s/prokka/db\";"%(Settings.DB_DIR))
               kronalibf.close()
               kronalibf = open("%s%scpp%s%s-%s/prokka/bin/prokka"%(Settings.METAMOS_UTILS,os.sep,os.sep, Settings.OSTYPE, Settings.MACHINETYPE), 'w')
               kronalibf.write(dd)
               kronalibf.close()


           # also need to change phylosift to external DB
           os.system("cp %s%sphylosift%sphylosiftrc %s%sphylosift%sphylosiftrc.orig"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.METAMOSDIR, os.sep, os.sep))
           testIn = open("%s%sphylosift%sphylosiftrc.orig"%(Settings.METAMOSDIR, os.sep, os.sep), 'r')
           testOut = open("%s%sphylosift%sphylosiftrc"%(Settings.METAMOSDIR, os.sep, os.sep), 'w')
           for line in testIn.xreadlines():
              if "marker_path" in line:
                 testOut.write("$marker_path=\"%s%sshare%sphylosift\";\n"%(Settings.DB_DIR, os.sep, os.sep))
              elif "ncbi_path" in line:
                 testOut.write("$ncbi_path=\"%s%sshare%sphylosift\";\n"%(Settings.DB_DIR, os.sep, os.sep))
              else:
                 testOut.write(line.strip() + "\n")
           testIn.close()
           testOut.close()

      Settings.SIGNALP       = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)

      Settings.FCP           = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.PHMMER        = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.MGCAT         = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.BLAST         = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.PHYLOSIFT     = "%s%sPhyloSift"%(Settings.METAMOSDIR, os.sep)

      Settings.KRONA         = "%s%sKronaTools%sbin"%(Settings.METAMOSDIR,os.sep,os.sep)
      if _BINARY_DIST and runp:
          #need to change KronaTools.pm to external Taxonomy directory
           kronalibf = open("%s%sKronaTools%slib%sKronaTools.pm"%(Settings.METAMOSDIR,os.sep,os.sep,os.sep))
           data = kronalibf.read()
           if "my $taxonomyDir = \"$libPath/../taxonomy\";" not in data:
               kronalibf.close()
           else:
               dd = data.replace("my $taxonomyDir = \"$libPath/../taxonomy\";","my $taxonomyDir = \"%s/taxonomy\";"%(Settings.DB_DIR))
               kronalibf.close()
               kronalibf = open("%s%sKronaTools%slib%sKronaTools.pm"%(Settings.METAMOSDIR,os.sep,os.sep,os.sep),'w')
               kronalibf.write(dd)
               kronalibf.close()
               os.system("ln -s %s/taxonomy %s%sKronaTools%staxonomy"%(Settings.DB_DIR,Settings.METAMOSDIR,os.sep,os.sep))
      Settings.REPEATOIRE    = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.LAP	     = "%s%sLAP"%(Settings.METAMOSDIR, os.sep)
      Settings.ALE           = "%s%scpp%s%s-%s/ALE/src"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.CGAL          = "%s%scpp%s%s-%s/cgal"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.REAPR          = "%s%scpp%s%s-%s/REAPR"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.FRCBAM        = "%s%scpp%s%s-%s/FRCbam/bin"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.FREEBAYES     = "%s%scpp%s%s-%s/freebayes/bin"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.QUAST         = "%s%squast"%(Settings.METAMOSDIR, os.sep)

      Settings.MPI	     = "%s%smpiexec"%(Settings.METAMOSDIR, os.sep)


libcounter = 1
readcounter = 1
class Read:
    format = ""
    maxlen = 150
    qformat = "Sanger"
    filtered = False
    mated = True
    path = ""
    fname = ""
    id = 0
    sid = ""
    def __init__(self,format,path,mated=True,interleaved=False):
        global readcounter 
        self.id = readcounter
        readcounter +=1
        self.format = format
        self.path = path
        self.fname = os.path.basename(self.path)
        self.mated = mated
        self.interleaved = interleaved
        #self.init()
        #self.validate()

class readLib:
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
    f12 = ""
    reads = []
    readDict = {}
    pairDict = {}
    def __init__(self,format,mmin,mmax,f1,f2="",mated=True,interleaved=False,innie=True,linkerType="titanium"):
        global libcounter
        self.id = libcounter
        self.sid = "lib"+str(libcounter)
        libcounter +=1
        self.format = format
        self.mated=mated
        self.interleaved=interleaved
        self.innie=innie
        self.linkerType=linkerType
        self.mmin = mmin
        self.mmax = mmax
        self.f1 = f1
        self.f2 = f2
        self.f1.sid = self.sid
        self.readDict[f1.id] = self.f1
        if f2 != "":
            self.readDict[f2.id] = self.f2
            self.pairDict[f1.id] = f2.id
            self.pairDict[f2.id] = f1.id
            self.f2.sid = self.sid
        self.reads = []
        self.reads.append(f1)
        if self.f2 != "":
            self.reads.append(f2)
        self.initLib()
        self.validateLib()
    def getPair(self,readId):
        try:
            return self.readDict[self.pairDict[readId]]
        except KeyError:
            #no pair for read
            return -1
    def initLib(self):
        self.mean = (self.mmin+self.mmax)/2.0
        self.stdev = 0.1*self.mean
        #count num reads
        #check pairs
        #if self.interleaved:
        #    f12 = self.f1.path
        #else:
        #need to shuffle em
        #    if self.f1.format == "fasta" and self.f2.format == "fasta":
        #        pass
        #    elif self.f2.format = "fastq" and self.f1.format == "fastq":
        #        pass
        #    f12 = ""
    def validateLib(self):
        pass

    def concatLib(self):
        #combine two libs of same format & w/ same insert size into one 
        pass
   
    def toggleInterleaved(self):
        #if interleaved switch to two files, else vice versa
        pass

    def filterReads(self):
        #remove all Reads with N, etc
        pass

    def __str__(self):
        pass

def getDefaultWeight(sa):
    if sa == SCORE_TYPE.LAP or sa == SCORE_TYPE.ALE or sa == SCORE_TYPE.CGAL:
       return 0.333333333
    elif sa == SCORE_TYPE.ORF:
       return 0
    else:
       return 1

def nearlyEqual(a, b, epsilon = 0.0001):
    absA = abs(float(a))
    absB = abs(float(b))
    diff = abs(float(a) - float(b))

    if a == b:
        return True
    elif (a * b == 0): # a or b or both are zero
        # relative error is not meaningful here
        return diff < (epsilon * epsilon)
    else: # use relative error
        return diff / (absA + absB) < epsilon

def initValidationScores(weights = dict()):
   for score in SCORE_TYPE.reverse_mapping.keys():
      if score in weights:
         SCORE_WEIGHTS[score] = weights[score]
      elif len(weights) == 0:
         SCORE_WEIGHTS[score] = getDefaultWeight(score)
      else:
         SCORE_WEIGHTS[score] = 0

def updateConfigCommands(infileName, opts):
   # build the list of commands
   commands = ""
   for o, a in opts:
      if o == "-f" or o == "--force":
         continue
      if o == "-d" or o == "--projectdir":
         continue
      if "--" in o:
         commands = "%s %s=%s"%(commands, o, a)
      else:
         commands = "%s %s %s"%(commands, o, a)

   tempFileName = "%s.tmp"%(infileName)
   tempFile = open(tempFileName, 'w')
   infile = open(infileName, 'r')
   for line in infile.xreadlines():
      if "command:" in line:
          tempFile.write("command:\t%s\n"%(commands.strip()))
      else:
          tempFile.write(line)
   infile.close()
   tempFile.close()
   os.system("mv %s %s"%(tempFileName, infileName))

def updateLibInfo(infileName, lib):
   tempFileName = "%s.tmp"%(infileName)
   tempFile = open(tempFileName, 'w')
   infile = open(infileName, 'r')
   written = False
   for line in infile.xreadlines():
      if "lib%d"%(lib.id) in line:
         if written == False:
            written = True
            tempFile.write("lib%dformat:\t%s\n"%(lib.id, lib.format))
            tempFile.write("lib%dmated:\t%s\n"%(lib.id, lib.mated))
            tempFile.write("lib%dinnie:\t%s\n"%(lib.id, lib.innie))
            tempFile.write("lib%dinterleaved\t%s\n"%(lib.id, lib.interleaved))
            if lib.mated:
               if lib.interleaved:
                  tempFile.write("lib%df1:\t%s,%d,%d,%d,%d\n"%(lib.id, lib.f1.fname, lib.mmin, lib.mmax, lib.mean, lib.stdev))
               else:
                  tempFile.write("lib%df1:\t%s,%d,%d,%d,%d\n"%(lib.id, lib.f1.fname, lib.mmin, lib.mmax, lib.mean, lib.stdev))
                  tempFile.write("lib%df2:\t%s,%d,%d,%d,%d\n"%(lib.id, lib.f2.fname, lib.mmin, lib.mmax, lib.mean, lib.stdev))
            else:
               tempFile.write("lib%dfrg:\t%s\n"%(lib.id, lib.f1.fname))
      else:
          tempFile.write(line)
   infile.close()
   tempFile.close()
   os.system("mv %s %s"%(tempFileName, infileName))

def readConfigInfo(infile, filePrefix=""):
   readlibs = []
   asmcontigs = []
   workflow = ""

   libs = []
   readobjs = []
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
   libadded = False
   nlib = None
   lib = None

   for line in infile.xreadlines():
      line = line.replace("\n","")

      if "#" in line:
         continue
      elif "inherit:" in line:
         wfc = line.replace("\n", "").split(":")
         if len(wfc) < 2:
            continue
         workflow = wfc[1].strip()
      elif "asmcontigs:" in line:
         asmc = line.replace("\n","").split("asmcontigs:")
         if len(asmc) < 2 or len(asmc[1].strip()) == 0:
            continue
         contigs = asmc[1].strip().split(",")
         for contig in contigs:
            if (len(contig.strip()) > 0):
               asmcontigs.append(contig)
      elif "format:" in line:
         if f1 and not libadded:
            nread1 = Read(format,f1,mated,interleaved)
            readobjs.append(nread1)
            nread2 = ""
            nlib = readLib(format,mmin,mmax,nread1,nread2,mated,interleaved,innie,linkerType)
            readlibs.append(nlib)
         libadded = False
         format = line.replace("\n","").split(":")[-1].strip()
      elif "mated:" in line:
         mated = str2bool(line.replace("\n","").split(":")[-1].strip())
      elif "interleaved:" in line:
         interleaved = str2bool(line.replace("\n","").split(":")[-1].strip())
      elif "innie:" in line:
         innie = str2bool(line.replace("\n","").split(":")[-1].strip())
      elif "linker:" in line:
          linkerType = line.replace("\n","").split(":")[-1].strip()
      elif "f1:" in line:
         data = line.split("f1:")
         f1 = "%s%s"%(filePrefix, data[1].strip().split(",")[0])
         inf = data[1].strip().split(",")
         mean = int(inf[3])
         stdev = int(inf[4])
         mmin = int(inf[1])
         mmax = int(inf[2])
         libs.append(f1)

      elif "f2:" in line:
         data = line.split("f2:")
         f2 = "%s%s"%(filePrefix,data[1].strip().split(",")[0])
         inf = data[1].split(",")
         mean = int(inf[3])
         stdev = int(inf[4])
         mmin = int(inf[1])
         mmax = int(inf[2])
         libs.append(f2)
        
         nread1 = Read(format,f1,mated,interleaved)
         readobjs.append(nread1)
         nread2 = Read(format,f2,mated,interleaved)
         readobjs.append(nread2)
         nlib = readLib(format,mmin,mmax,nread1,nread2,mated,interleaved,\
                              innie,linkerType)
         readlibs.append(nlib)
         libadded = True
      elif "frg" in line:
         data = line.split("frg:")
         frg = "%s%s"%(filePrefix,data[1].strip().split(",")[0])
         mated = False
         f1 = frg
         libs.append(frg)

   if f1 and not libadded:
      nread1 = Read(format,f1,mated,interleaved)
      readobjs.append(nread1)
      nread2 = ""
      nlib = readLib(format,mmin,mmax,nread1,nread2,mated,interleaved,innie,\
                          linkerType)
      readlibs.append(nlib)

   return (asmcontigs, readlibs, workflow)

def concatContig(ctgfile):
    if len(sys.argv) < 3:
        print "usage: contig_file out_file"
    contig_file = open(ctgfile,'r')
    out_file = open(ctgfile+".merged",'w')
    out_data = ""
    for line in contig_file.xreadlines():
        if ">" not in line:
             out_data += line.replace("\n","")
    width = 60
    pp = 0
    out_file.write(">seq\n")
    while pp+60 < len(out_data):
        out_file.write(out_data[pp:pp+60]+"\n")
        pp +=60
    out_file.write(out_data[pp:]+"\n")
    out_file.close()
    contig_file.close()

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def sizeFastaFile(fileName):
   if not os.path.exists(fileName):
      return 0

   p = subprocess.Popen("java -cp %s/java:. SizeFasta -t %s"%(Settings.METAMOS_UTILS, fileName), shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr != "":
      print "Warning: cannot size file, return size 0\n"
      return 0
   else:
      return int(checkStdout)

def getMD5Sum(fileName):
   if not os.path.exists(fileName):
      return ""

   md5 = hashlib.md5()
   with open(fileName,'rb') as f: 
      for chunk in iter(lambda: f.read(128*md5.block_size), ''): 
         md5.update(chunk)

   return md5.hexdigest()

def getMachineType():
   p = subprocess.Popen("echo `uname`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr != "":
      print "Warning: Cannot determine OS, defaulting to %s"%(Settings.OSTYPE)
   else:
      Settings.OSTYPE = checkStdout.strip()

   p = subprocess.Popen("echo `uname -r`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr != "":
      print "Warning: Cannot determine OS version, defaulting to %s"%(Settings.OSVERSION)
   else:
      Settings.OSVERSION = checkStdout.strip()

   p = subprocess.Popen("echo `uname -m`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr != "":
      print "Warning: Cannot determine system type, defaulting to %s"%(Settings.MACHINETYPE)
   else:
      Settings.MACHINETYPE = checkStdout.strip()

def getCommandOutput(theCommand, checkForStderr):
    p = subprocess.Popen(theCommand, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (checkStdout, checkStderr) = p.communicate()
    if checkForStderr and checkStderr != "":
       return ""
    else:
       return checkStdout.strip()

def getFromPath(theCommand, theName, printWarning = True):
    deprecated_list = ["METAIDBA","PHYMM"]
    result = getCommandOutput("which %s"%(theCommand), True)
    if theName.upper() not in deprecated_list and result == "" and printWarning:
       print "Warning: %s is not found, some functionality will not be available"%(theName)
       return ""
    else:
       return os.path.dirname(result.strip())

def cmdExists(cmd):
    result = False
    try:
       result = subprocess.call(cmd,
           shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
    except OSError:
       result = False

    return result

def initConfig(kmer, threads, theRundir, taxaLevel, localKrona, annotateUnmapped, verbose, outputOnly, doScaffolding):
    Settings(kmer, threads, theRundir, taxaLevel, localKrona, annotateUnmapped, verbose, outputOnly, doScaffolding, True)

    getMachineType()

    if not os.path.exists(Settings.METAMOS_UTILS):
       Settings.METAMOSDIR = os.getcwd()
       print "Script is running from: %s"%(Settings.METAMOSDIR)
   
       Settings.METAMOS_UTILS = "%s%sUtilities"%(Settings.METAMOSDIR, os.sep) 
       if not os.path.exists(Settings.METAMOS_UTILS):
          print "Error: cannot find metAMOS utilities. Will not run pipeline"
          sys.exit(1)   

       Settings.METAMOS_JAVA  = "%s%sjava:%s"%(Settings.METAMOS_UTILS, os.sep, os.curdir)
       Settings.METAMOS_DOC   = "%s%sdoc"%(Settings.METAMOS_UTILS, os.sep)

    # FastQC
    Settings.FASTQC = "%s%sFastQC"%(Settings.METAMOSDIR, os.sep)
    if not os.path.exists(Settings.FASTQC + os.sep + "fastqc"):
       Settings.FASTQC = getFromPath("fastqc", "FastQC")
    fastqcMD5 = getMD5Sum(Settings.FASTQC + os.sep + "fastqc")
    
    # now check for assemblers
    # 1. AMOS
    Settings.AMOS = "%s%sAMOS%s%s-%s%sbin"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
    if not os.path.exists(Settings.AMOS + os.sep + "toAmos_new"):
       Settings.AMOS = getFromPath("toAmos_new", "AMOS") 
       if not os.path.exists(Settings.AMOS + os.sep + "toAmos_new"):
          print "Error: cannot find AMOS in %s or %s. Please check your path and try again."%(Settings.METAMOSDIR + os.sep + "AMOS", Settings.AMOS)
          sys.exit(1)
    amosMD5 = getMD5Sum(Settings.AMOS + os.sep + "toAmos_new")
    Settings.BAMBUS2 = Settings.AMOS
    bambusMD5 = getMD5Sum(Settings.BAMBUS2 + os.sep + "OrientContigs")

    # 2. Soap
    Settings.SOAPDENOVO = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE) 
    if not os.path.exists(Settings.SOAPDENOVO + os.sep + "SOAPdenovo-63mer"):
       Settings.SOAPDENOVO = ""
    soapMD5 = getMD5Sum(Settings.SOAPDENOVO + os.sep + "SOAPdenovo-63mer")

    Settings.SOAPDENOVO2 = "%s%scpp%s%s-%s/soap2/bin"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.SOAPDENOVO2 + os.sep + "SOAPdenovo-63mer"):
       Settings.SOAPDENOVO2 = ""
    soapMD5 = getMD5Sum(Settings.SOAPDENOVO2 + os.sep + "SOAPdenovo-63mer")

    # 3. CA
    Settings.CA = "%s%sCA%s%s-%s%sbin"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE.replace("x86_64","amd64"), os.sep)
    if not os.path.exists(Settings.CA + os.sep + "gatekeeper"):
       Settings.CA = getFromPath("gatekeeper", "Celera Assembler") 
    CAMD5 = getMD5Sum(Settings.CA + os.sep + "gatekeeper")

    # BLASR goes with CA
    Settings.BLASR = "%s/../../../smrtanalysis/current/analysis/bin"%(Settings.CA)
    if not os.path.exists(Settings.BLASR + os.sep + "blasr"):
       Settings.BLASR = getFromPath("blasr", "BLASR")
    blasrMD5 = getMD5Sum(Settings.BLASR + os.sep + "blasr")

    # 4. Newbler
    Settings.NEWBLER = "%s%snewbler%s%s-%s"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.NEWBLER + os.sep + "runProject"):
       Settings.NEWBLER = getFromPath("runProject", "Newbler")
    newblerMD5 = getMD5Sum(Settings.NEWBLER + os.sep + "runProject")

    # 5. meta-IDBA
    Settings.METAIDBA = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE) 
    if not os.path.exists(Settings.METAIDBA + os.sep + "metaidba"):
       Settings.METAIDBA = getFromPath("metaidba", "METAIDBA")
    metaidbaMD5 = getMD5Sum(Settings.METAIDBA + os.sep + "metaidba")

    # when searching for velvet, we ignore paths because there are so many variations of velvet (velvet, velvet-sc, meta-velvet that all have a velveth/g and we have no way to tell if we got the right one
    #6. velvet
    Settings.VELVET = "%s%scpp%s%s-%s%svelvet"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
    if not os.path.exists(Settings.VELVET + os.sep + "velvetg"):
       Settings.VELVET = ""
    velvetMD5 = getMD5Sum(Settings.VELVET + os.sep + "velvetg")

    #7. velvet-sc
    Settings.VELVET_SC = "%s%scpp%s%s-%s%svelvet-sc"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
    if not os.path.exists(Settings.VELVET_SC + os.sep + "velvetg"):
       Settings.VELVET_SC = ""
    velvetSCMD5 = getMD5Sum(Settings.VELVET_SC + os.sep + "velvetg")

    #8. metavelvet
    Settings.METAVELVET = "%s%scpp%s%s-%s%sMetaVelvet"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
    if not os.path.exists(Settings.METAVELVET + os.sep + "meta-velvetg"):
       Settings.METAVELVET = getFromPath("meta-velvetg", "METAVELVET")
    metaVelvetMD5 = getMD5Sum(Settings.METAVELVET + os.sep + "meta-velvetg")

    # 8. SparseAssembler
    Settings.SPARSEASSEMBLER = "%s%scpp%s%s-%s%sSparseAssembler"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
    if not os.path.exists(Settings.SPARSEASSEMBLER + os.sep + "SparseAssembler"):
       Settings.SPARSEASSEMBLER = getFromPath("SparseAssembler", "SparseAssembler")
    sparseAssemblerMD5 = getMD5Sum(Settings.SPARSEASSEMBLER + os.sep + "SparseAssembler")

    Settings.KRONA         = "%s%sKronaTools%sbin"%(Settings.METAMOSDIR,os.sep,os.sep)
    if not os.path.exists(Settings.KRONA + os.sep + "ktImportTaxonomy"):
       Settings.KRONA = getFromPath("Krona", "ktImportTaxonomy")
    kronaMD5 = getMD5Sum(Settings.KRONA + os.sep + "ktImportTaxonomy")

    # now for repeatoire
    Settings.REPEATOIRE = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.REPEATOIRE + os.sep + "repeatoire"):
       Settings.REPEATOIRE = getFromPath("repeatoire", "Repeatoire")
    else:
       Settings.REPEATOIRE += os.sep + "repeatoire"
    repeatoireMD5 = getMD5Sum(Settings.REPEATOIRE + os.sep + "repeatoire")

    # now for the mappers
    Settings.BOWTIE = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.BOWTIE + os.sep + "bowtie"):
       Settings.BOWTIE = getFromPath("bowtie", "Bowtie")
    bowtieMD5 = getMD5Sum(Settings.BOWTIE + os.sep + "bowtie")

    Settings.BOWTIE2 = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.BOWTIE2 + os.sep + "bowtie2"):
       Settings.BOWTIE2 = getFromPath("bowtie2", "Bowtie2")
    bowtie2MD5 = getMD5Sum(Settings.BOWTIE + os.sep + "bowtie2")

    Settings.SAMTOOLS = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.SAMTOOLS + os.sep + "samtools"):
       Settings.SAMTOOLS = getFromPath("samtools", "samtools")
    samtoolsMD5 = getMD5Sum(Settings.SAMTOOLS + os.sep + "samtools")

    # now the gene callers
    Settings.METAGENEMARK = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.METAGENEMARK + os.sep + "gmhmmp"):
       Settings.METAGENEMARK = getFromPath("gmhmmp", "MetaGeneMark")
    gmhmmpMD5 = getMD5Sum(Settings.METAGENEMARK + os.sep + "gmhmmp")

    Settings.PROKKA = "%s%scpp%s%s-%s/prokka/bin"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.PROKKA + os.sep + "prokka"):
       Settings.PROKKA = getFromPath("prokka", "Prokka")
    prokkaMD5 = getMD5Sum(Settings.PROKKA + os.sep + "prokka")

    Settings.SIGNALP = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.SIGNALP + os.sep + "signalp"):
       Settings.SIGNALP = getFromPath("signalp", "SignalP+")
    signalpMD5 = getMD5Sum(Settings.SIGNALP + os.sep + "signalp")

    Settings.FRAGGENESCAN = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.FRAGGENESCAN + os.sep + "FragGeneScan"):
       Settings.FRAGGENESCAN = getFromPath("FragGeneScan","FragGeneScan")
    fraggenescanMD5 = getMD5Sum(Settings.FRAGGENESCAN + os.sep + "FragGeneScan")

    # now for the annotation
    Settings.METAPHYLER = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.METAPHYLER + os.sep + "metaphylerClassify"):
       Settings.METAPHYLER = getFromPath("metaphylerClassify", "metaphylerClassify")
    metaphylerMD5 = getMD5Sum(Settings.METAPHYLER + os.sep + "metaphylerClassify")

    Settings.FCP = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.FCP + os.sep + "nb-classify"):
       Settings.FCP = getFromPath("nb-classify", "FCP")
    fcpMD5 = getMD5Sum(Settings.FCP + os.sep + "nb-classify")

    Settings.PHMMER = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.PHMMER + os.sep + "phmmer"):
       Settings.PHMMER = getFromPath("phmmer", "PHmmer")
    phmmerMD5 = getMD5Sum(Settings.PHMMER + os.sep + "phmmer")

    Settings.MGCAT = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.MGCAT + os.sep + "mgcat"):
       Settings.MGCAT = getFromPath("mgcat", "mgcat")
    mgcatMD5 = getMD5Sum(Settings.MGCAT + os.sep + "mgcat")

    Settings.PHYMM = "%s%sperl%sphymm%s"%(Settings.METAMOS_UTILS, os.sep, os.sep,os.sep)
    if not os.path.exists(Settings.PHYMM + os.sep + "scoreReads.pl"):
       Settings.PHYMM = getFromPath("phymm", "Phymm")
    phymmMD5 = getMD5Sum(Settings.PHYMM + os.sep + "scoreReads.pl")

    Settings.BLAST = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.BLAST + os.sep + "blastall"):
       Settings.BLAST = getFromPath("blastall", "blast")
    blastMD5 = getMD5Sum(Settings.BLAST + os.sep + "blastall")

    # currently only supported on Linux 64-bit and only from one location
    Settings.PHYLOSIFT = "%s%sphylosift"%(Settings.METAMOSDIR, os.sep)
    if not os.path.exists(Settings.PHYLOSIFT + os.sep + "bin" + os.sep + "phylosift"):
       print "Warning: PhyloSift was not found, will not be available\n"
       Settings.PHYLOSIFT = ""
    phylosiftMD5 = getMD5Sum(Settings.PHYLOSIFT + os.sep + "bin" + os.sep + "phylosift")

    Settings.EAUTILS = "%s%scpp%s%s-%s%seautils"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
    if not os.path.exists(Settings.EAUTILS + os.sep + "fastq-mcf"):
       Settings.EAUTILS = getFromPath("fastq-mcf", "EA-UTILS")
    eautilsMD5 = getMD5Sum(Settings.EAUTILS + os.sep + "fastq-mcf")

    Settings.KMERGENIE = "%s%scpp%s%s-%s%skmergenie"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
    if not os.path.exists(Settings.KMERGENIE + os.sep + "kmergenie"):
       Settings.KMERGENIE = getFromPath("kmergenie", "KmerGenie")
    kmergenieMD5 = getMD5Sum(Settings.KMERGENIE + os.sep + "kmergenie")

    Settings.R = "%s%scpp%s%s-%s%sR"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
    if not os.path.exists(Settings.R + os.sep + "R"):
       Settings.R = getFromPath("R", "R package")
    rMD5 = getMD5Sum(Settings.R + os.sep + "R")

    # now for the validators
    Settings.LAP = "%s%sLAP"%(Settings.METAMOSDIR, os.sep)
    if not os.path.exists(Settings.LAP + os.sep + "aligner" + os.sep + "calc_prob.py"):
       Settings.LAP = getFromPath("calc_prop.py", "LAP")
    lapMD5 = getMD5Sum(Settings.LAP + os.sep + "aligner" + os.sep + "calc_prob.py")

    Settings.ALE = "%s%scpp%s%s-%s/ALE/src"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.ALE + os.sep + "ALE"):
       Settings.ALE = getFromPath("ALE", "ALE")
    aleMD5 = getMD5Sum(Settings.ALE + os.sep + "ALE")

    Settings.CGAL = "%s%scpp%s%s-%s/cgal"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.CGAL + os.sep + "cgal"):
       Settings.CGAL = getFromPath("cgal", "CGAL")
    cgalMD5 = getMD5Sum(Settings.CGAL + os.sep + "cgal")

    Settings.REAPR = "%s%scpp%s%s-%s/REAPR"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.REAPR + os.sep + "reapr"):
       Settings.REAPR = getFromPath("reapr", "REAPR")
    reaprMD5 = getMD5Sum(Settings.REAPR + os.sep + "reapr")
    
    Settings.FRCBAM = "%s%scpp%s%s-%s/FRCbam/bin"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.FRCBAM + os.sep + "FRC"):
       Settings.FRCBAM = getFromPath("FRC", "FRCbam")
    frcMD5 = getMD5Sum(Settings.FRCBAM + os.sep + "FRC")

    Settings.FREEBAYES = "%s%scpp%s%s-%s/freebayes/bin"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.FREEBAYES + os.sep + "freebayes"):
       Settings.FREEBAYES = getFromPath("freebayes", "FreeBayes")
    freebayesMD5 = getMD5Sum(Settings.FREEBAYES + os.sep + "freebayes")

    Settings.QUAST = "%s%squast"%(Settings.METAMOSDIR, os.sep)
    if not os.path.exists(Settings.QUAST + os.sep + "quast.py"):
       Settings.QUAST = getFromPath("quast.py", "QUAST")
    quastMD5 = getMD5Sum(Settings.QUAST + os.sep + "quast.py")

    Settings.MPI = "%s%smpiexec"%(Settings.METAMOSDIR, os.sep)
    if not os.path.exists(Settings.MPI):
       Settings.MPI = getFromPath("mpiexec", "MPI", False)
       if Settings.MPI == "":
          Settings.MPI = getFromPath("openmpiexec", "OPENMPI", False)
          if Settings.MPI != "":
             Settings.MPI = "%s%s%s"%(Settings.MPI, os.sep, "openmpiexec")
       else:
          Settings.MPI = "%s%s%s"%(Settings.MPI, os.sep, "mpiexec")
    if not os.path.exists(Settings.MPI):
       print "Warning: MPI is not available, some functionality may not be available"
    mpiMD5 = getMD5Sum(Settings.MPI)

    # finally store the configuration 
    conf = open("%s/pipeline.conf"%(Settings.rundir),'w')
    if Settings.BINARY_DIST and 1: 
          prevtmpdirs = []
          try:
              bdf = open("%s/prevruns.tmp"%(application_path),'r')
              for line in bdf.xreadlines():
                  prevtmpdirs.append(line.replace("\n",""))
              for pdir in prevtmpdirs:
                  if os.path.exists("%s"%(pdir)):
                      os.system("rm -rf %s"%(pdir))
              bdf.close()
              bdf = open("%s/prevruns.tmp"%(application_path),'w')
              bdf.close()

          except IOError:
              #do not have permissions to write to install dir, store in tmp?
              #tf, tf_path = tempfile.mkstemp("prevruns.tmp",'w')
              bdf = open("%s/prevruns.tmp"%(application_path),'w')
              bdf.write("%s\n"%(sys._MEIPASS))
              bdf.close()

          except TypeError:
              bdf = open("%s/prevruns.tmp"%(application_path),'w')
              bdf.write("%s\n"%(sys._MEIPASS))
              bdf.close()

    conf.write("#Configuration summary\n")
    conf.write("OS:\t\t\t%s\nOS Version:\t\t%s\nMachine:\t\t%s\n"%(Settings.OSTYPE, Settings.OSVERSION, Settings.MACHINETYPE))
    conf.write("metAMOS main dir:\t%s\nmetAMOS Utilities:\t%s\nmetAMOS Java:\t\t%s\n"%(Settings.METAMOSDIR, Settings.METAMOS_UTILS, Settings.METAMOS_JAVA))
    conf.write("AMOS:\t\t\t%s\t%s\n"%(Settings.AMOS, amosMD5))
    conf.write("BAMBUS2:\t\t%s\t%s\n"%(Settings.BAMBUS2, bambusMD5))
    conf.write("SOAPDENOVO:\t\t\t%s\t%s\n"%(Settings.SOAPDENOVO, soapMD5))
    conf.write("SOAPDENOVO2:\t\t\t%s\t%s\n"%(Settings.SOAPDENOVO2, soapMD5))
    conf.write("METAIDBA:\t\t%s\t%s\n"%(Settings.METAIDBA, metaidbaMD5))
    conf.write("Celera Assembler:\t%s\t%s\n"%(Settings.CA, CAMD5))
    conf.write("NEWBLER:\t\t%s\t%s\n"%(Settings.NEWBLER, newblerMD5))
    conf.write("Velvet:\t\t\t%s\t%s\nVelvet-SC:\t\t%s\t%s\n"%(Settings.VELVET, velvetMD5, Settings.VELVET_SC, velvetSCMD5))
    conf.write("MetaVelvet:\t\t%s\t%s\n"%(Settings.METAVELVET, metaVelvetMD5))
    conf.write("SparseAssembler:\t%s\t%s\n"%(Settings.SPARSEASSEMBLER, sparseAssemblerMD5))
    conf.write("metaphylerClassify:\t\t\t%s\t%s\n"%(Settings.METAPHYLER, metaphylerMD5))
    conf.write("Bowtie:\t\t\t%s\t%s\n"%(Settings.BOWTIE, bowtieMD5))
    conf.write("Bowtie2:\t\t\t%s\t%s\n"%(Settings.BOWTIE2, bowtie2MD5))
    conf.write("samtools:\t\t\t%s\t%s\n"%(Settings.SAMTOOLS, samtoolsMD5))
    conf.write("M-GCAT:\t\t\t%s\t%s\n"%(Settings.MGCAT, mgcatMD5))
    conf.write("METAGENEMARK:\t\t\t%s\t%s\n"%(Settings.METAGENEMARK, gmhmmpMD5))
    conf.write("FRAGGENESCAN:\t\t%s\t%s\n"%(Settings.FRAGGENESCAN, fraggenescanMD5))
    conf.write("PROKKA:\t\t\t%s\t%s\n"%(Settings.PROKKA, prokkaMD5))
    conf.write("SIGNALP:\t\t\t%s\t%s\n"%(Settings.SIGNALP, signalpMD5))
    conf.write("FCP:\t\t\t%s\t%s\n"%(Settings.FCP, fcpMD5))
    conf.write("PHMMER:\t\t\t%s\t%s\n"%(Settings.PHMMER, phmmerMD5))
    conf.write("PHYMM:\t\t\t%s\t%s\n"%(Settings.PHYMM, phymmMD5))
    conf.write("BLAST:\t\t\t%s\t%s\n"%(Settings.BLAST, blastMD5))
    conf.write("PHYLOSIFT:\t\t%s\t%s\n"%(Settings.PHYLOSIFT, phylosiftMD5))
    conf.write("FASTQC:\t\t\t%s\t%s\n"%(Settings.FASTQC, fastqcMD5))
    conf.write("EAUTILS:\t\t%s\t%s\n"%(Settings.EAUTILS, eautilsMD5))
    conf.write("KMERGENIE:\t\t%s\t%s\n"%(Settings.KMERGENIE, kmergenieMD5))
    conf.write("REPEATOIRE:\t\t%s\t%s\n"%(Settings.REPEATOIRE, repeatoireMD5))
    conf.write("KRONA:\t\t\t%s\t%s\n"%(Settings.KRONA, kronaMD5))
    conf.write("LAP:\t\t\t%s\t%s\n"%(Settings.LAP, lapMD5))
    conf.write("ALE:\t\t\t%s\t%s\n"%(Settings.ALE, aleMD5))
    conf.write("CGAL:\t\t\t%s\t%s\n"%(Settings.CGAL, cgalMD5))
    conf.write("REAPR:\t\t\t%s\t%s\n"%(Settings.REAPR, reaprMD5))
    conf.write("FRCBAM:\t\t\t%s\t%s\n"%(Settings.FRCBAM, frcMD5))
    conf.write("FREEBAYES:\t\t\t%s\t%s\n"%(Settings.FREEBAYES, freebayesMD5))
    conf.write("QUAST:\t\t\t%s\t%s\n"%(Settings.QUAST, quastMD5))
    conf.close()

    return Settings

def setFailFast(fail):
   global _failFast

   _failFast = fail

def run_process(settings,command,step=""):
       outf = ""
       workingDir = ""
       if step != "":
           workingDir = "%s/%s/out"%(settings.rundir, step)
           if not os.path.exists(workingDir):
              workingDir = ""
           step = string.upper(step)
           if not os.path.exists(settings.rundir+os.sep+"Logs"):
               # create Log directory
               os.system("mkdir %s/Logs"%(settings.rundir))

               # create the log of commands
               commandf = open(settings.rundir + os.sep + "Logs" + os.sep + "COMMANDS.log", 'w')
               commandf.close()

           # open command log file for appending (it should have been created above)
           commandf = open(settings.rundir + os.sep + "Logs" + os.sep + "COMMANDS.log", 'a')

           if step not in settings.task_dict:
              print "Starting Task = %s.%s"%(step.lower(), step)
              dt = datetime.now().isoformat(' ')[:-7]
              commandf.write("|%s|# [%s]\n"%(dt,step))
              outf = open(settings.rundir+os.sep+"Logs"+os.sep+step+".log",'w')
              settings.task_dict.append(step)

              # create started file
              startedf = open(settings.rundir + os.sep + "Logs" + os.sep + step.lower() + ".started", 'w')
              startedf.close()
           else:
              outf = open(settings.rundir+os.sep+"Logs"+os.sep+step+".log",'a')

       if settings.VERBOSE or settings.OUTPUT_ONLY:
           print "*** metAMOS running command: %s\n"%(command)

       if settings.OUTPUT_ONLY == False:
          stdout = ""
          stderr = ""
          if workingDir == "":
              p = subprocess.Popen(command, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,close_fds=True,executable="/bin/bash")
          else:
              p = subprocess.Popen(command, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,close_fds=True,executable="/bin/bash", cwd=workingDir)
          fstdout,fstderr = p.communicate()
          rc = p.returncode
          if rc != 0 and _failFast and "rm " not in command and "ls " not in command and "unlink " not in command and "ln " not in command and "mkdir " not in command and "mv " not in command and "cat" not in command:
              # flush all error/output streams
              outf.flush()
              outf.write(fstdout+fstderr)
              outf.close()
              commandf.flush()
              dt = datetime.now().isoformat(' ')[:-7]
              commandf.write("|%s| "%(dt)+command+"\n")
              commandf.close()

              global _atomicCounter
              if _atomicCounter.increment() == 0: 
                 print ERROR_RED+"*****************************************************************"
                 print "*************************ERROR***********************************"
                 print "During %s, the following command failed with return code %d:"%(step.lower(), rc)
                 print ">>",command
                 print ""
                 print "*************************DETAILS***********************************"
                 print "Last %d commands run before the error (%s/Logs/COMMANDS.log)"%(_NUM_LINES, settings.rundir)
                 p = subprocess.Popen("tail -n %d %s/Logs/COMMANDS.log"%(_NUM_LINES, settings.rundir), shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,close_fds=True, executable="/bin/bash")
                 (checkStdout, checkStderr) = p.communicate()
                 val = p.returncode
                 print "%s"%(checkStdout)
                 print "Last %d lines of output (%s/Logs/%s.log)"%(_NUM_LINES, settings.rundir, step)
                 p = subprocess.Popen("tail -n %d %s/Logs/%s.log"%(_NUM_LINES, settings.rundir, step), shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,close_fds=True, executable="/bin/bash")
                 (checkStdout, checkStderr) = p.communicate()
                 val = p.returncode
                 print "%s"%(checkStdout)
                 print "Please veryify input data and restart MetAMOS. If the problem persists please contact the MetAMOS development team."
                 print "*************************ERROR***********************************"
                 print "*****************************************************************"+ENDC

              # also make sure this step will be re-run on restart
              os.system("rm %s%sLogs%s%s.ok"%(settings.rundir, os.sep, os.sep, step.lower())) 
              #sys.exit(rc)
              raise 
          if step == "":
              print fstdout,fstderr
          else:
              outf.write(fstdout+fstderr)
              outf.close()
              dt = datetime.now().isoformat(' ')[:-7]
              commandf.write("|%s| "%(dt)+command+"\n")
              commandf.close()

def recruitGenomes(settings,query,genomeDir,outDir,stepName, top=1):
   print "recruiting genomes.."
   setFailFast(False)
   run_process(settings, "%s/mgcat -M -r %s -d %s -o %s -p %d"%(settings.MGCAT,query,genomeDir,outDir,settings.threads), stepName.title())
   setFailFast(True)

   gtr = []
   if os.path.exists("%s/recruited_genomes.lst"%(outDir)):
      rg = open("%s/recruited_genomes.lst"%(outDir),'r')
      rglist = []
      cnt = 0
      for genome in rg.xreadlines():
         genome = genome.replace("\n","")
         seq,mumi = genome.split(",")
         if os.path.exists(seq):
            rglist.append([float(mumi),seq])
            cnt +=1
      print "done! recruited %d genomes!"%(cnt)
      rglist.sort()
      i = 0
      while i < len(rglist) and i < top:
         gtr.append(rglist[i][1])
         i+=1
   else:
      print "Error: recruiting references failed"

   return gtr

      
def getProgramCitations(settings, programName, comment="#"):
   global _PUB_DICT
   global _PROG_NAME_DICT
   cite = ""
   if len(_PUB_DICT) == 0:
      try:
         cite = open("%s/%s"%(settings.METAMOS_DOC, "citations.rst"), 'r')
      except IOError:
         #print "no citations file! cannot print!"
         return ("","")

      for line in cite.xreadlines():
         (line, sep, commentLine) = line.partition(comment)
         splitLine = line.strip().split("\t")
         if len(splitLine) >= 3:
            name = splitLine[0]
            commonName = splitLine[1]
            citation = splitLine[2]
         elif len(splitLine) >= 2:
           name = splitLine[0]
           commonName = splitLine[1]
           citation = "NA" 
         else:
            continue
         _PROG_NAME_DICT[name] = commonName
         _PUB_DICT[name] = citation

   try:
      return (_PROG_NAME_DICT[programName], _PUB_DICT[programName]) 
   except KeyError:
      return(programName, "UNKNOWN")

def getProgramParams(configDir, fileName, module="", prefix="", comment="#", separator=""):
    # we process parameters in the following priority:
    # first: current directory
    # second: user home directory
    # third: metAMOS directory
    # a parameter specifeid in the current directory takes priority over all others, and so on down the line
    dirs = [configDir + os.sep + "config", os.path.expanduser('~') + os.sep + ".metAMOS", os.getcwd()]
    optDict = {} 

    cmdOptions = ""

    for curDir in dirs:
       spec = ""
       curFile = curDir + os.sep + fileName
       try:
          spec = open(curFile, 'r')
       except IOError as e:
          continue

       read = False
       if module == "":
          read = True

       for line in spec.xreadlines():
          (line, sep, commentLine) = line.partition(comment)
          line = line.strip()

          if line == "[" + module + "]":
             read = True
             continue
          elif read == True and line.startswith("["):
             break

          if read:
             if (line != ""):
                if (line.endswith("\\")):
                   for next in spec:
                      next = next.strip()
                      line = line.replace("\\", "") + next.replace("\\", "")
                      if not next.endswith("\\"):
                         break
                splitLine = line.split();
                optDict[splitLine[0]] = separator.join(splitLine[1:]).strip() 
       spec.close()

    for option in optDict:
       cmdOptions += prefix + option + " " + optDict[option] + " "

    return cmdOptions

def getAvailableMemory(settings):
   if settings.nopsutil:
      return 0

   import psutil

   cacheusage=0
   if 'linux' in settings.OSTYPE.lower():
      cacheusage = psutil.cached_phymem()
   memusage =  `psutil.phymem_usage()`.split(",")
   freemem = long(memusage[2].split("free=")[-1])+long(cacheusage)
   percentfree = float(memusage[3].split("percent=")[-1].split(")")[0])
   avram = (freemem/1000000000)

   return avram

def getSelectedAssembler(settings):
   if settings.rundir == "":
      print "Error: attempted to get selected assembler before initialization"
      raise (JobSignalledBreak)
   elif not os.path.exists("%s/Validate/out/%s.asm.selected"%(settings.rundir, settings.PREFIX)):
      print "Error: attempted to get selected assembler before validation"
      raise (JobSignalledBreak)
   else:
      return getCommandOutput("cat %s/Validate/out/%s.asm.selected"%(settings.rundir, settings.PREFIX), False)

def getSelectedKmer(settings):
   kmer = ""
   if os.path.exists("%s/Assemble/out/%s.kmer"%(settings.rundir, settings.PREFIX)):
      stats = open("%s/Assemble/out/%s.kmer"%(settings.rundir, settings.PREFIX), 'r')
      kmer = stats.read().strip()
      stats.close()
   return kmer

def getEstimatedGenomeSize(settings):
   genomeSize = 0
   if os.path.exists("%s/Assemble/out/%s.genomesize"%(settings.rundir, settings.PREFIX)):
      stats = open("%s/Assemble/out/%s.genomesize"%(settings.rundir, settings.PREFIX), 'r')
      genomeSize  = int(stats.read().strip())
      stats.close()
   return genomeSize 

def getVersion():
   #look for pattern like: MetAMOS [VERSION] README
   version = "UNKNOWN"
   filePath = "%s%sREADME.md"%(sys.path[0], os.sep)
   try:
      sys._MEIPASS
      filePath = "%s%sREADME.md"%(sys._MEIPASS, os.sep)
   except Exception:
      filePath = "%s%sREADME.md"%(sys.path[0], os.sep)
   
   if os.path.exists(filePath):
      readme_file = open(filePath, 'r')
      for line in readme_file.xreadlines():
         if "# MetAMOS" in line:
            version = line.strip().split("# MetAMOS")[1]
            version = version.strip().split("README")[0]
            break
      readme_file.close()

   import workflow
   wfs = workflow.getSupportedWorkflowNames("%s/Utilities/workflows"%(sys.path[0]), False)

   return version + " workflows: " + ",".join(wfs)

def configureEnvironment(utilPath):
   global _envCounter
   if _envCounter.increment() == 0:
      if "PYTHONPATH" not in os.environ:
         os.environ["PYTHONPATH"] = ""
      else:
         ppath = os.environ["PYTHONPATH"]
         #os.environ["PYTHONPATH"] = ""
      os.environ["PYTHONPATH"]+=utilPath+os.sep+"python"+os.pathsep
      os.environ["PYTHONPATH"]+=utilPath+os.sep+"ruffus"+os.pathsep
      os.environ["PYTHONPATH"] += utilPath+os.sep+"python"+os.sep+"lib"+os.pathsep
      os.environ["PYTHONPATH"] += utilPath+os.sep+"python"+os.sep+"lib"+os.sep+"python"+os.pathsep
      os.environ["PYTHONPATH"] += utilPath+os.sep+"python"+os.sep+"lib64"+os.pathsep
      os.environ["PYTHONPATH"] += utilPath+os.sep+"python"+os.sep+"lib64"+os.sep+"python"+os.pathsep
      os.environ["PYTHONPATH"] += utilPath+os.pathsep

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
         site.addsitedir(utilPath+os.sep+"python"+os.sep+"lib"+os.sep+"python")
         site.addsitedir(utilPath+os.sep+"python"+os.sep+"lib64"+os.sep+"python")

         sys.path.append(utilPath)
         sys.path.append(utilPath+os.sep+"python")
         sys.path.append(utilPath+os.sep+"ruffus")
         sys.path.append(utilPath+os.sep+"python"+os.sep+"lib"+os.sep+"python")
         sys.path.append(utilPath+os.sep+"python"+os.sep+"lib64"+os.sep+"python")
         try:
            sys.path.append(sys._MEIPASS)
         except Exception:
            pass
         sys.path.append("/usr/lib/python")

         #remove imports from pth file, if exists
         nf = []
         if 'bash' in shellv or cmdExists('export'):
            os.system("export PYTHONPATH=%s:$PYTHONPATH"%(utilPath+os.sep+"python"))
            os.system("export PYTHONPATH=%s:$PYTHONPATH"%(utilPath+os.sep+"python"+os.sep+"lib"+os.sep+"python"))
         elif cmdExists('setenv'):
            os.system("setenv PYTHONPATH %s:$PYTHONPATH"%(utilPath+os.sep+"python"))
            os.system("setenv PYTHONPATH %s:$PYTHONPATH"%(utilPath+os.sep+"python"+os.sep+"lib"+os.sep+"python"))
         else:
            print "Warning: could not set PYTHONPATH. Unknown shell %s, some functionality may not work\n"%(shellv)

   # finally set LD path
   libPath = os.path.abspath(utilPath + os.sep + ".." + os.sep + "lib")
   if os.path.exists(libPath):
      oldLDPath = ""
      needToAdd = True
      if "LD_LIBRARY_PATH" in os.environ:
        oldLDPath = os.environ["LD_LIBRARY_PATH"]
        if libPath in oldLDPath:
           needToAdd = False
      elif "DYLD_FALLBACK_LIBRARY_PATH" in os.environ:
         oldLDPath = os.environ["DYLD_FALLBACK_LIBRARY_PATH"]
         if libPath in oldLDPath:
            needToAdd = False
      if needToAdd:
         os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = libPath + os.pathsep + oldLDPath
         os.environ["LD_LIBRARY_PATH"] = libPath + os.pathsep + oldLDPath

def translateToSRAURL(settings, name):
   oldDyLD = ""
   if "DYLD_FALLBACK_LIBRARY_PATH" in os.environ:
      oldDyLD = os.environ["DYLD_FALLBACK_LIBRARY_PATH"]
      del os.environ["DYLD_FALLBACK_LIBRARY_PATH"]
   command = "%s/cpp%s%s-%s%s/sra/bin"%(settings.METAMOS_UTILS, os.sep, settings.OSTYPE, settings.MACHINETYPE, os.sep)
   result = getCommandOutput("%s/srapath %s"%(command, name), True)
   if result == name:
      result = ""

   if oldDyLD != "":
      os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = oldDyLD
   return result
