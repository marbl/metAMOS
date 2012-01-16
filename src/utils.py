#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

import hashlib

_METAMOSDIR    = sys.path[0]
INITIAL_UTILS = "%s%sUtilities"%(_METAMOSDIR, os.sep)

class Settings:
   asmfiles = []
   runfiles = []

   kmer = 55
   threads = 16
   rundir = ""
   PREFIX = ""
   VERBOSE = False
   OSTYPE = ""
   OSVERSION = ""
   MACHINETYPE = ""

   METAMOSDIR = ""
   METAMOS_UTILS = ""
   METAMOS_JAVA = ""

   FASTQC = ""
   AMOS = ""

   SOAP = ""
   METAIDBA = ""
   CA = ""
   NEWBLER = ""
   VELVET = ""
   VELVET_SC = ""

   BOWTIE = ""

   GMHMMP = ""

   FCP = ""
   PHMMER = ""
   BLAST = ""
   AMPHORA = ""

   KRONA = ""
   REPEATOIRE = ""

   def __init__(self, kmer = None, threads = None, rundir = None, update = False):

      if (Settings.rundir != "" and update == False):
         return

      if (kmer == None or threads == None or rundir == None):
         print "Error settings is uninitialized and no intialization provided\n"
         raise(Exception)

      Settings.rundir = rundir
      Settings.kmer = kmer
      Settings.threads = threads 
      Settings.rundir = rundir

      Settings.PREFIX = "proba"
      Settings.VERBOSE = False
      Settings.OSTYPE        = "Linux"
      Settings.OSVERSION     = "0.0"
      Settings.MACHINETYPE   = "x86_64"

      Settings.METAMOSDIR    = sys.path[0]
      Settings.METAMOS_UTILS = "%s%sUtilities"%(Settings.METAMOSDIR, os.sep) 
      Settings.METAMOS_JAVA  = "%s%sjava:%s"%(Settings.METAMOS_UTILS,os.sep,os.curdir)

      Settings.AMOS          = "%s%sAMOS%sbin"%(Settings.METAMOSDIR, os.sep, os.sep)

      Settings.SOAP          = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.METAIDBA      = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.CA            = "%s%sCA%s%s-%s%sbin"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE.replace("x86_64", "amd64"), os.sep)
      Settings.NEWBLER       = "%s%snewbler%s%s-%s"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.VELVET        = "%s%scpp%s%s-%s%svelvet"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)
      Settings.VELVET_SC     = "%s%scpp%s%s-%s%svelvet-sc"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep)

      Settings.BOWTIE        = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)

      Settings.GMHMMP        = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)

      Settings.FCP           = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.PHMMER        = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.BLAST         = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
      Settings.AMPHORA       = "%s%sAmphora-2"%(Settings.METAMOSDIR, os.sep)

      Settings.KRONA         = "%s%skrona"%(Settings.METAMOS_UTILS,os.sep)
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

def getFromPath(theCommand, theName):
    p = subprocess.Popen("which %s"%(theCommand), shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (checkStdout, checkStderr) = p.communicate()
    if checkStderr != "":
       print "Warning: %s is not found, some functionality will not be available"%(theName)
       return ""
    else:
       return checkStdout.replace(theCommand, "").strip()

def initConfig(kmer, threads, theRundir):
    Settings(kmer, threads, theRundir, True)

    getMachineType()

    if not os.path.exists(Settings.METAMOS_UTILS):
       Settings.METAMOSDIR = sys.path[0]
       print "Script is running from: %s"%(Settings.METAMOSDIR)
   
       Settings.METAMOS_UTILS = "%s%sUtilities"%(Settings.METAMOSDIR, os.sep) 
       if not os.path.exists(Settings.METAMOS_UTILS):
          print "Error: cannot find metAMOS utilities. Will not run pipeline"
          sys.exit(1);   

       Settings.METAMOS_JAVA  = "%s%sjava:%s"%(Settings.METAMOS_UTILS, os.sep, os.curdir)

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

    # 2. Soap
    Settings.SOAP = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE) 
    if not os.path.exists(Settings.SOAP + os.sep + "soap63"):
       Settings.SOAP = getFromPath("soap63", "SOAP")
    soapMD5 = getMD5Sum(Settings.SOAP + os.sep + "soap63")

    # 3. CA
    Settings.CA = "%s%sCA%s%s-%s%sbin"%(Settings.METAMOSDIR, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE.replace("x86_64","amd64"), os.sep)
    if not os.path.exists(Settings.CA + os.sep + "gatekeeper"):
       Settings.CA = getFromPath("gatekeeper", "Celera Assembler") 
    CAMD5 = getMD5Sum(Settings.CA + os.sep + "gatekeeper")

    # 4. Newbler
    Settings.NEWBLER = "%s%snewbler"%(Settings.METAMOSDIR, os.sep);
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
    Settings.VELVET = "%s%scpp%s%s-%s%svelvet"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep);
    if not os.path.exists(Settings.VELVET + os.sep + "velvetg"):
       Settings.VELVET = ""
    velvetMD5 = getMD5Sum(Settings.VELVET + os.sep + "velvetg")

    #7. velvet-sc
    Settings.VELVET_SC = "%s%scpp%s%s-%s%svelvet-sc"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE, os.sep);
    if not os.path.exists(Settings.VELVET_SC + os.sep + "velvetg"):
       Settings.VELVET_SC = ""
    velvetSCMD5 = getMD5Sum(Settings.VELVET_SC + os.sep + "velvetg")

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

    # now for the annotation
    Settings.GMHMMP = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.GMHMMP + os.sep + "gmhmmp"):
       Settings.GMHMMP = getFromPath("gmhmmp", "GeneMark.hmm")
    gmhmmpMD5 = getMD5Sum(Settings.GMHMMP + os.sep + "gmhmmp")

    Settings.FCP = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.FCP + os.sep + "fcp"):
       Settings.FCP = getFromPath("fcp", "FCP")
    fcpMD5 = getMD5Sum(Settings.FCP + os.sep + "fcp")

    Settings.PHMMER = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.PHMMER + os.sep + "phmmer"):
       Settings.PHMMER = getFromPath("phmmer", "PHmmer")
    phmmerMD5 = getMD5Sum(Settings.PHMMER + os.sep + "phmmer")

    Settings.BLAST = "%s%scpp%s%s-%s"%(Settings.METAMOS_UTILS, os.sep, os.sep, Settings.OSTYPE, Settings.MACHINETYPE)
    if not os.path.exists(Settings.BLAST + os.sep + "blastall"):
       Settings.BLAST = getFromPath("blastall", "blast")
    blastMD5 = getMD5Sum(Settings.BLAST + os.sep + "blastall")

    # currently only supported on Linux 64-bit and only from one location
    Settings.AMPHORA = "%s%sAmphora-2"%(Settings.METAMOSDIR, os.sep)
    if not os.path.exists(Settings.AMPHORA + os.sep + "amphora2"):
       print "Warning: Amphora 2 was not found, will not be available\n";
       Settings.AMPHORA = ""
    if Settings.AMPHORA != "" and (Settings.OSTYPE != "Linux" or Settings.MACHINETYPE != "x86_64"):
       print "Warning: Amphora 2 not compatible with %s-%s. It requires Linux-x86_64\n"%(Settings.OSTYPE, Settings.MACHINETYPE)
       Settings.AMPHORA = "" 
    amphoraMD5 = getMD5Sum(Settings.AMPHORA + os.sep + "amphora2")

    # finally store the configuration 
    conf = open("%s/pipeline.conf"%(Settings.rundir),'w')

    conf.write("#Configuration summary\n")
    conf.write("THREADS:\t\t\t%d\n"%(Settings.threads))
    conf.write("KMER:\t\t\t%d\n"%(Settings.kmer))
    conf.write("PREFIX:\t\t\t%s\n"%(Settings.PREFIX))
    conf.write("VERBOSE:\t\t%s\n"%(Settings.VERBOSE)) 
    conf.write("OS:\t\t\t%s\nOS Version:\t\t%s\nMachine:\t\t%s\n"%(Settings.OSTYPE, Settings.OSVERSION, Settings.MACHINETYPE))
    conf.write("metAMOS main dir:\t%s\nmetAMOS Utilities:\t%s\nmetAMOS Java:\t\t%s\n"%(Settings.METAMOSDIR, Settings.METAMOS_UTILS, Settings.METAMOS_JAVA))
    conf.write("AMOS:\t\t\t%s\t%s\n"%(Settings.AMOS, amosMD5))
    conf.write("SOAP:\t\t\t%s\t%s\n"%(Settings.SOAP, soapMD5))
    conf.write("METAIDBA:\t\t%s\t%s\n"%(Settings.METAIDBA, metaidbaMD5))
    conf.write("Celera Assembler:\t%s\t%s\n"%(Settings.CA, CAMD5))
    conf.write("NEWBLER:\t\t%s\t%s\n"%(Settings.NEWBLER, newblerMD5))
    conf.write("Velvet:\t\t\t%s\t%s\nVelvet-SC:\t\t%s\t%s\n"%(Settings.VELVET, velvetMD5, Settings.VELVET_SC, velvetSCMD5))
    conf.write("Bowtie:\t\t\t%s\t%s\n"%(Settings.BOWTIE, bowtieMD5))
    conf.write("GMHMMP:\t\t\t%s\t%s\n"%(Settings.GMHMMP, gmhmmpMD5))
    conf.write("FCP:\t\t\t%s\t%s\n"%(Settings.FCP, fcpMD5))
    conf.write("PHMMER:\t\t\t%s\t%s\n"%(Settings.PHMMER, phmmerMD5))
    conf.write("BLAST:\t\t\t%s\t%s\n"%(Settings.BLAST, blastMD5))
    conf.write("AMPHORA:\t\t%s\t%s\n"%(Settings.AMPHORA, amphoraMD5))
    conf.write("FASTQC:\t\t%s\t%s\n"%(Settings.FASTQC, fastqcMD5))

    conf.write("REPEATOIRE:\t\t%s\t%s\n"%(Settings.REPEATOIRE, repeatoireMD5))
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
           if not os.path.exists(settings.rundir+"/Logs"):
               os.system("mkdir %s/Logs"%(settings.rundir))
           outf = open(settings.rundir+"/Logs/"+step+".log",'a')
       if settings.VERBOSE:
           print command
       stdout = ""
       stderr = ""
       if workingDir == "":
           p = subprocess.Popen(command, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,close_fds=True,executable="/bin/bash")
       else:
           p = subprocess.Popen(command, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,close_fds=True,executable="/bin/bash", cwd=workingDir)
       fstdout,fstderr = p.communicate()

       if step == "":
           print fstdout,fstderr
       else:
           outf.write(fstdout+fstderr)
           outf.close()


def getProgramParams(configDir, fileName, module="", prefix="", comment="#"):
    # we process parameters in the following priority:
    # first: current directory
    # second: user home directory
    # third: metAMOS directory
    # a parameter specifeid in the current directory takes priority over all others, and so on down the line
    dirs = [configDir + os.sep + "config", os.path.expanduser('~') + os.sep + ".metAMOS", os.getcwd()]
    optDict = {} 

    cmdOptions = ""

    for curDir in dirs:
       curFile = curDir + os.sep + fileName;
       try:
          spec = open(curFile, 'r')
       except IOError as e:
          continue

       read = False
       if module == "":
          read = True

       for line in spec:
          (line, sep, commentLine) = line.partition(comment)
          line = line.strip()

          if line == "[" + module + "]":
             read = True
             continue;
          elif read == True and line.startswith("["):
             break;

          if read:
             if (line != ""):
                splitLine = line.split();
                optDict[splitLine[0]] = "".join(splitLine[1:]).strip() 
       spec.close()

    for option in optDict:
       cmdOptions += prefix + option + " " + optDict[option] + " ";

    return cmdOptions
