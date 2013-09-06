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

STEP_NAMES = enum("ASSEMBLE", "ANNOTATE")
STEP_OUTPUTS = enum(".asm.contig", ".hits")
INPUT_TYPE = enum("FASTQ", "FASTA", "CONTIGS", "SCAFFOLDS", "ORF_FA", "ORF_AA")

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

class Settings:
   asmfiles = []
   runfiles = []

   kmer = 55
   threads = 16
   rundir = ""
   taxa_level = "class"
   local_krona = False
   annotate_unmapped = False
   task_dict = []
   noblastdb = False
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
   METAIDBA = ""
   CA = ""
   NEWBLER = ""
   VELVET = ""
   VELVET_SC = ""
   METAVELVET = ""
   SPARSEASSEMBLER = ""

   MGCAT = ""

   METAPHYLER = ""
   BOWTIE = ""
   BOWTIE2 = ""
   SAMTOOLS = ""

   METAGENEMARK = ""
   FRAGGENESCAN = ""
   FCP = ""
   PHMMER = ""
   PHYMM = ""
   BLAST = ""
   PHYLOSIFT = ""
   DB_DIR = ""
   BLASTDB_DIR = ""
   KRONA = ""
   REPEATOIRE = ""

   BINARY_DIST = 0

   nopsutil = False
   nopysam = False

   def __init__(self, kmer = None, threads = None, rundir = None, taxa_level = "", localKrona = False, annotateUnmapped = False, verbose = False, outputOnly = False, update = False):

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
      Settings.annotate_unmapped = annotateUnmapped
      Settings.task_dict = []

      Settings.PREFIX = "proba"
      Settings.VERBOSE = verbose
      Settings.OUTPUT_ONLY = outputOnly

      Settings.OSTYPE        = "Linux"
      Settings.OSVERSION     = "0.0"
      Settings.MACHINETYPE   = "x86_64"


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
              if len(os.environ["BLASTDB"]) != 0:
                  _BLASTDB_PATH == os.environ["BLASTDB"]
                  if not os.path.exists(_BLASTDB_PATH):
                      print "Error: cannot find BLAST DB directory, yet path set via $BLASTDB: %s. Disabling blastdb dependent programs"%(os.environ["BLASTDB"])
                      Settings.noblastdb = True
                      #sys.exit(1)
              #print "BINARY DIST", _DB_PATH
          except KeyError:
              #_DB_PATH = "./DB/"
              pass

          if not os.path.exists(_DB_PATH):
              print "Error: cannot find DB directory in %s, was it deleted? oops, it is required to run MetAMOS!"%(_DB_PATH)
              sys.exit(1)
             
      Settings.DB_DIR        = _DB_PATH 
      Settings.BLASTDB_DIR   = _BLASTDB_PATH 
      Settings.BINARY_DIST   = _BINARY_DIST
      Settings.AMOS          = "%s%sAMOS%sbin"%(Settings.METAMOSDIR, os.sep, os.sep)
      Settings.BAMBUS2       = Settings.AMOS

      Settings.SOAPDENOVO    = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.METAIDBA      = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.CA            = "%s%sCA%s%s-%s%sbin"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE.replace("x86_64", "amd64"), os.sep)
      Settings.NEWBLER       = "%s%snewbler%s%s-%s"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.VELVET        = "%s%scpp%s%s-%s%svelvet"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.VELVET_SC     = "%s%scpp%s%s-%s%svelvet-sc"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.METAVELVET    = "%s%scpp%s%s-%s%sMetaVelvet"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.SPARSEASSEMBLER = "%s%scpp%s%s-%s%sSparseAssembler"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.PHYMM = "%s%sperl%sphymm%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, os.sep)

      Settings.METAPHYLER        = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)

      Settings.BOWTIE        = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.BOWTIE2        = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.SAMTOOLS        = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)

      Settings.METAGENEMARK  = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.FRAGGENESCAN  = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)

      Settings.FCP           = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.PHMMER        = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.MGCAT         = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.BLAST         = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.PHYLOSIFT     = "%s%sPhyloSift"%(Settings.METAMOSDIR, os.sep)

      Settings.KRONA         = "%s%sKronaTools%sbin"%(Settings.METAMOSDIR,os.sep,os.sep)
      if _BINARY_DIST:
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

def getFromPath(theCommand, theName):
    result = getCommandOutput("which %s"%(theCommand), True)
    if result == "":
       print "Warning: %s is not found, some functionality will not be available"%(theName)
       return ""
    else:
       return result.replace(theCommand, "").strip()

def cmdExists(cmd):
    result = False
    try:
       result = subprocess.call(cmd,
           shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
    except OSError:
       result = False

    return result

def initConfig(kmer, threads, theRundir, taxaLevel, localKrona, annotateUnmapped, verbose, outputOnly):
    Settings(kmer, threads, theRundir, taxaLevel, localKrona, annotateUnmapped, verbose, outputOnly, True)

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
    if not os.path.exists(Settings.SOAPDENOVO + os.sep + "soap63"):
       Settings.SOAPDENOVO = getFromPath("soap63", "SOAPDENOVO")
    soapMD5 = getMD5Sum(Settings.SOAPDENOVO + os.sep + "soap63")

    # 3. CA
    Settings.CA = "%s%sCA%s%s-%s%sbin"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE.replace("x86_64","amd64"), os.sep)
    if not os.path.exists(Settings.CA + os.sep + "gatekeeper"):
       Settings.CA = getFromPath("gatekeeper", "Celera Assembler") 
    CAMD5 = getMD5Sum(Settings.CA + os.sep + "gatekeeper")

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
    metaVelvetMD5 = getMD5Sum(Settings.SOAPDENOVO + os.sep + "meta-velvetg")

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

    # now for the annotation
    Settings.METAGENEMARK = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.METAGENEMARK + os.sep + "gmhmmp"):
       Settings.METAGENEMARK = getFromPath("gmhmmp", "MetaGeneMark")
    gmhmmpMD5 = getMD5Sum(Settings.METAGENEMARK + os.sep + "gmhmmp")
    
    Settings.METAPHYLER = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.METAPHYLER + os.sep + "metaphylerClassify"):
       Settings.METAPHYLER = getFromPath("metaphylerClassify", "metaphylerClassify")
    metaphylerMD5 = getMD5Sum(Settings.METAPHYLER + os.sep + "metaphylerClassify")

    Settings.FRAGGENESCAN = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.FRAGGENESCAN + os.sep + "FragGeneScan"):
       Settings.FRAGGENESCAN = getFromPath("FragGeneScan","FragGeneScan")
    fraggenescanMD5 = getMD5Sum(Settings.FRAGGENESCAN + os.sep + "FragGeneScan")

    Settings.FCP = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.FCP + os.sep + "nb-classify"):
       Settings.FCP = getFromPath("nb-classify", "FCP")
    fcpMD5 = getMD5Sum(Settings.FCP + os.sep + "nb-classify")

    Settings.PHMMER = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.PHMMER + os.sep + "phmmer"):
       Settings.PHMMER = getFromPath("phmmer", "PHmmer")
    phmmerMD5 = getMD5Sum(Settings.PHMMER + os.sep + "phmmer")

    Settings.MGCAT = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.PHMMER + os.sep + "mgcat"):
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
    conf.write("FCP:\t\t\t%s\t%s\n"%(Settings.FCP, fcpMD5))
    conf.write("PHMMER:\t\t\t%s\t%s\n"%(Settings.PHMMER, phmmerMD5))
    conf.write("PHYMM:\t\t\t%s\t%s\n"%(Settings.PHYMM, phymmMD5))
    conf.write("BLAST:\t\t\t%s\t%s\n"%(Settings.BLAST, blastMD5))
    conf.write("PHYLOSIFT:\t\t%s\t%s\n"%(Settings.PHYLOSIFT, phylosiftMD5))
    conf.write("FASTQC:\t\t\t%s\t%s\n"%(Settings.FASTQC, fastqcMD5))

    conf.write("REPEATOIRE:\t\t%s\t%s\n"%(Settings.REPEATOIRE, repeatoireMD5))
    conf.write("KRONA:\t\t\t%s\t%s\n"%(Settings.KRONA, kronaMD5))
    conf.close()

    return Settings

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
          if rc != 0 and "rm " not in command and "ls " not in command and "unlink " not in command and "ln " not in command and "mkdir " not in command and "mv " not in command:
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
              raise (JobSignalledBreak)
          if step == "":
              print fstdout,fstderr
          else:
              outf.write(fstdout+fstderr)
              outf.close()
              dt = datetime.now().isoformat(' ')[:-7]
              commandf.write("|%s| "%(dt)+command+"\n")
              commandf.close()

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
