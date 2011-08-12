import os, sys, string, time, BaseHTTPServer, getopt, subprocess#
from operator import itemgetter

PREFIX = "proba"
VERBOSE = False
OSTYPE        = "Linux"
OSVERSION     = "0.0"
MACHINETYPE   = "x86_64"

METAMOSDIR    = sys.path[0]
METAMOS_UTILS = "%s%sUtilities"%(METAMOSDIR, os.sep) 
METAMOS_JAVA  = "%s%sjava:%s"%(METAMOS_UTILS,os.sep,os.curdir)
AMOS          = "%s%sAMOS%sbin"%(METAMOSDIR, os.sep, os.sep)
SOAP          = "%s%scpp"%(METAMOS_UTILS, os.sep)
CA            = "%s%sCA%s%s-%s%sbin"%(METAMOSDIR, os.sep, os.sep, OSTYPE, MACHINETYPE.replace("x86_64", "amd64"), os.sep)
NEWBLER       = "%s%snewbler"%(METAMOSDIR, os.sep)
BOWTIE        = "%s%scpp"%(METAMOS_UTILS, os.sep)
GMHMMP        = "%s%scpp"%(METAMOS_UTILS, os.sep)
libcounter = 1
readcounter = 1
t1 = time.time()#clock()
sys.path.append(METAMOS_UTILS)
from ruffus import *


class Read:
    format = ""
    maxlen = 150
    qformat = "Sanger"
    filtered = False
    mated = True
    path = ""
    fname = ""
    id = 0
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
    linkerType = "titanium"
    frg = ""
    f1 = ""
    f2 = ""
    f12 = ""
    reads = []
    readDict = {}
    pairDict = {}
    def __init__(self,format,mmin,mmax,f1,f2="",mated=True,interleaved=False):
        global libcounter
        self.id = libcounter
        self.sid = "lib"+str(libcounter)
        libcounter +=1
        self.format = format
        self.mated=mated
        self.interleaved=interleaved
        self.mmin = mmin
        self.mmax = mmax
        self.f1 = f1
        self.f2 = f2
        self.readDict[f1.id] = self.f1
        if f2 != "":
            self.readDict[f2.id] = self.f2
            self.pairDict[f1.id] = f2.id
            self.pairDict[f2.id] = f1.id

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

     
def parseInterleaved(rf,wf,fastq=True):
    if 1:
        if 1:
            if 1:
                   #this means we have this entire lib in one file
                   #parse out paired record (8 lines), rename header to be filename + "/1" or "/2", and remove reads with N
                   rf = open(read.path,'r')
                   wf = open(read.path.replace("/in/","/out/"),'w')
                   start = 1
                   rcnt = 0
                   recordcnt = 0
                   record = []
                   shdr = ""
                   for line in rf.xreadlines():
                       if start:
                           s1hdr = line
                           record.append(line)
                           start = 0
                           rcnt =1

                       else:
                           if rcnt == 7:
                               #end of record
                               record.append(line)
                               rcnt +=1
                               if len(record) != 8:
                                   #something went wrong
                                   continue
                               rseq = string.upper(record[0]+record[5])                               
                               if "N" in rseq:
                                   #skip both, dont' want Ns
                                   continue
                               #update hdrs to be filename /1 or /2
                               recordcnt +=1
                               hdr = lib.sid+"r"+str(recordcnt)+"/"
                               if fastq == True:
                                   wf.writelines("@"+hdr+"1\n")
                                   wf.writelines(record[1])
                                   wf.writelines("+"+hdr+"1\n")
                                   wf.writelines(record[3])
                                   wf.writelines("@"+hdr+"2\n")
                                   wf.writelines(record[5])
                                   wf.writelines("+"+hdr+"2\n")
                                   wf.writelines(record[7])
                               else:
                                   wf.writelines(">"+hdr+"1\n")
                                   wf.writelines(record[1])
                                   wf.writelines(">"+hdr+"1\n")
                                   wf.writelines(record[3])
                                   wf.writelines(">"+hdr+"2\n")
                                   wf.writelines(record[5])
                                   wf.writelines(">"+hdr+"2\n")
                                   wf.writelines(record[7])
                           elif rcnt % 4 == 0:
                               s2hdr = line
                               rlcs = LCS(s1hdr,s2hdr)
                               #these should almost identical
                               if float(len(rlcs))/float(len(s1hdr)) < 0.9:
                                   #missing record somewhere, start over with this one
                                   s1hdr = line
                                   record = [line]
                                   start = 0
                                   rcnt = 1
                               else:
                                   record.append(line)
                                   rcnt +=1
                           elif rcnt % 2 == 0:
                               #quality hdr
                               record .append(line)
                               rcnt +=1
                           else:
                               record .append(line)
                               rcnt +=1
                   #update to new path
                   read.path = read.path.replace("/in/","/out/")            


def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def getMachineType():
   global OSTYPE
   global OSVERSION
   global MACHINETYPE

   p = subprocess.Popen("echo `uname`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr != "":
      print "Warning: Cannot determine OS, defaulting to %s"%(OSTYPE)
   else:
      OSTYPE = checkStdout.strip()

   p = subprocess.Popen("echo `uname -r`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr != "":
      print "Warning: Cannot determine OS version, defaulting to %s"%(OSVERSION)
   else:
      OSVERSION = checkStdout.strip()

   p = subprocess.Popen("echo `uname -m`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr != "":
      print "Warning: Cannot determine system type, defaulting to %s"%(MACHINETYPE)
   else:
      MACHINETYPE = checkStdout.strip()

def getFromPath(theCommand, theName):
    p = subprocess.Popen("which %s"%(theCommand), shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (checkStdout, checkStderr) = p.communicate()
    if checkStderr != "":
       print "Warning: %s is not found, some functionality will not be available"%(theName)
       return ""
    else:
       return checkStdout.replace(theCommand, "").strip()

def guessPaths():
    global METAMOSDIR
    global METAMOS_UTILS
    global METAMOS_JAVA
    global AMOS
    global SOAP
    global CA
    global NEWBLER
    global BOWTIE
    global GMHMMP

    getMachineType()

    if not os.path.exists(METAMOS_UTILS):
       METAMOSDIR = sys.path[0]
       print "Script is running from: %s"%(METAMOSDIR)
   
       METAMOS_UTILS = "%s%sUtilities"%(METAMOSDIR, os.sep) 
       if not os.path.exists(METAMOS_UTILS):
          print "Error: cannot find metAMOS utilities. Will not run pipeline"
          sys.exit(1);   

       METAMOS_JAVA  = "%s%sjava:%s"%(METAMOS_UTILS, os.sep, os.curdir)

    # now check for assemblers
    # 1. AMOS
    AMOS = "%s%sAMOS%sbin"%(METAMOSDIR, os.sep, os.sep)
    #if not os.path.exists(AMOS + os.sep + "bank-transact"):
    if not os.path.exists(AMOS + os.sep + "toAmos_new"):
       AMOS = getFromPath("bank-transact", "AMOS") 
    # 2. Soap
    SOAP = "%s%scpp"%(METAMOS_UTILS, os.sep) 
    if not os.path.exists(SOAP + os.sep + "SOAPdenovo-63mer"):
       SOAP = getFromPath("SOAPdenovo-63mer", "SOAP")
       #SOAP = getFromPath("SOAPdenovo-63mer", "SOAP")
    # 3. CA
    CA = "%s%sCA%s%s-%s%sbin"%(METAMOSDIR, os.sep, os.sep, OSTYPE, MACHINETYPE.replace("x86_64","amd64"), os.sep)
    if not os.path.exists(CA + os.sep + "gatekeeper"):
       CA = getFromPath("gatekeeper", "Celera Assembler") 
    # 4. Newbler
    NEWBLER = "%s%snewbler"%(METAMOSDIR, os.sep);
    #print NEWBLER
    #sys.exit(0)
    if not os.path.exists(NEWBLER + os.sep + "runProject"):
       NEWBLER = getFromPath("runProject", "Newbler")

    # now for the mappers
    BOWTIE = "%s%scpp"%(METAMOS_UTILS, os.sep)
    if not os.path.exists(BOWTIE + os.sep + "bowtie"):
       BOWTIE = getFromPath("bowtie", "Bowtie")

    # now for the annotation
    GMHMMP = "%s%scpp"%(METAMOS_UTILS, os.sep)
    if not os.path.exists(GMHMMP + os.sep + "gmhmmp"):
       GMHMMP = getFromPath("gmhmmp", "GeneMark.hmm")

    # finally add the utilities to our path
    print "Configuration summary:"
    print "OS:\t\t\t%s\nOS Version:\t\t%s\nMachine:\t\t%s\n"%(OSTYPE, OSVERSION, MACHINETYPE)
    print "metAMOS main dir:\t%s\nmetAMOS Utilities:\t%s\nmetAMOS Java:\t\t%s\n"%(METAMOSDIR, METAMOS_UTILS, METAMOS_JAVA)
    print "AMOS:\t\t\t%s\nSOAP:\t\t\t%s\nCelera Assembler:\t%s\nNEWBLER:\t\t%s\n"%(AMOS, SOAP, CA, NEWBLER)
    print "Bowtie:\t\t\t%s"%(BOWTIE)
    print "GMHMMP:\t\t\t%s"%(GMHMMP)

def getProgramParams(fileName, module="", prefix="", comment="#"):
    cmdOptions = ""
    read = False
    if module == "":
       read = True

    spec = open("%s/config/%s"%(METAMOS_UTILS, fileName),'r')
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
             cmdOptions += " " + prefix + line
    spec.close()
    return cmdOptions

def usage():
    print "usage: runPipeline.py [options] -d projectdir (required)"
    print "options:  -a <assembler> -k <kmer size> -c <classification method> -m <enable metaphyler?> -p <num threads>  "
    #print "options: annotate, stopafter, startafter, fq, fa"

try:
    opts, args = getopt.getopt(sys.argv[1:], "hb:d:s:e:o:k:c:a:n:p:tf:vm4", ["help", "bowtie","projectdir","startat","endat", "minoverlap","kmersize","classifier","assembler","skipsteps","threads","filter","forcesteps","verbose","metaphyler","454"])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

allsteps = ["Preprocess","Assemble","FindORFS","Metaphyler","Annotate","Scaffold","Propagate","Classify","Postprocess"]
output = None
reads = None
quals = None
format = None
verbose = False
bowtie_mapping = False
startat = None
stopat = None
filter = False
forcesteps = []
skipsteps = []
run_metaphyler = False
runfast = False
cls = "phmmer"
asm = "soap"
rundir = ""
fff = ""
threads = 16
readlen = 75
kmer = 31
fqlibs = {}
fqfrags = []
rlibs = []
for o, a in opts:
    if o in ("-v","--verbose"):
        VERBOSE = True
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-b","--bowtie"):
        bowtie_mapping = True
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
        kmer = int(a)
    elif o in ("-4", "--454"):
        fff = "-454"
    elif o in ("-f", "--forcesteps"):
        print o,a
        forcesteps = a.split(",")
        print forcesteps
    elif o in ("-n", "--skipsteps"):
        print o, a
        skipsteps = a.split(",")
        print skipsteps
    elif o in ("-p", "--threads"):
        threads = int(a)
    elif o in ("-d", "--projectdir"):
        rundir = a
        if not os.path.exists(a):
          print "project dir %s does not exist!"%(rundir)
          usage()
          sys.exit(1)
    elif o in ("-t", "--filter"):
        filter = True

    elif o in ("-m", "--metaphyler"):
        run_metaphyler = True
    elif o in ("-c", "--classifier"):
        #blast,fcp,etc 
        #default: fcp?
        cls = a#"phmmer"
    elif o in ("-a","--assembler"):
        #maximus,CA,soap
        #default: maximus?
        asm = a
    elif o in ("-f","--fastest"):
        #tweak all parameters to run fast
        #bambus2, use SOAP, etc
        runfast = True
    
    else:
        assert False, "unhandled option"

    #sys.exit(2)

if not os.path.exists(rundir) or rundir == "":
    print "project dir %s does not exist!"%(rundir)
    usage()
    sys.exit(1)

#parse frag/libs out of pipeline.ini out of rundir
inifile = os.curdir+os.sep+rundir+os.sep+"pipeline.ini"
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
linkerType = "titanium"
frg = ""
f1 = ""
f2 = ""
currlibno = 1
newlib = ""
for line in inf:
    line = line.replace("\n","")
    if "#" in line:
        continue
    elif "format:" in line:
        newlibno = int(line.split("format")[0].split("lib")[1])


        if newlibno != currlibno:
            nread1 = Read(format,f1,mated,interleaved)
            readobjs.append(nread1)
            nread2 = ""
            if f2 != "":
                nread2 = Read(format,f2,mated,interleaved)
                readobjs.append(nread2)
                #nread12 = Read(format,f1+".interleaved",mated,1)
            nlib = readLib(format,mmin,mmax,nread1,nread2,mated,interleaved)
            readlibs.append(nlib)
            
        format = line.replace("\n","").split("\t")[-1]
        currlibno = newlibno
    elif "mated:" in line:
        mated = str2bool(line.replace("\n","").split("\t")[-1])
    elif "interleaved:" in line:
        interleaved = str2bool(line.replace("\n","").split("\t")[-1])
    elif "linker:" in line:
        linkerType = line.replace("\n","").split("\t")[-1]
    elif "f1:" in line:# or "f2:" in line:
        data = line.split("\t")

        fqlibs[data[0]] = data[1]
        #f1 = data[1].split(",")[0]
        f1 = "%s/Preprocess/in/%s"%(rundir,data[1].split(",")[0])
        inf = data[1].split(",")
        mean = int(inf[3])
        stdev = int(inf[4])
        mmin = int(inf[1])
        mmax = int(inf[2])
        libs.append(f1)

    elif "f2:" in line:# or "f2:" in line:
        data = line.split("\t")

        fqlibs[data[0]] = data[1]
        f2 = "%s/Preprocess/in/%s"%(rundir,data[1].split(",")[0])
        inf = data[1].split(",")
        mean = int(inf[3])
        stdev = int(inf[4])
        mmin = int(inf[1])
        mmax = int(inf[2])
        libs.append(f2)
    elif "frg" in line:

        data = line.split("\t")
        frg = data[1]
        mated = False
        f1 = frg
        #fqfrags[data[0]] = data[1]
        #frgs.append(data[1])
        libs.append(frg)
nread1 = Read(format,f1,mated,interleaved)
nread2 = ""
readobjs.append(nread1)
if f2 != "":
    nread2 = Read(format,f2,mated,interleaved)
    readobjs.append(nread2)
nlib = readLib(format,mmin,mmax,nread1,nread2,mated,interleaved)
readlibs.append(nlib)
print len(readlibs)
def run_process(command):
       if VERBOSE:
           print command
       stdout = ""
       stderr = ""
       p = subprocess.Popen(command, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
       #p = subprocess.Popen(command, shell=True, stdin=None, stdout=stdout, stderr=stderr)
       if "clk" not in command and "Bundler" not in command and "MarkRepeats" not in command and VERBOSE:
           print p.stdout.read(), p.stderr.read()
       (checkStdout, checkStderr) = p.communicate()

def map2contig(fasta):
    bowtie_mapping = 1
    
    readDir = ""
    asmDir = ""
    threads = 0
    tigr_file = open("%s/Assemble/out/%s.asm.tigr"%(rundir,PREFIX),'w')
    contigfile = open("%s/Assemble/out/%s.asm.contig"%(rundir,PREFIX),'r')

    seqdict = {}
    hdr = ""
    cnt = 0
    contigdict = {}
    contigdict2 = {}
    readdict = {}
    matedict = {}
    ctgmates = 0
    matectgdict = {}
    mateotdict = {}
    read_lookup = {}
    readcnt = 1

    for lib in readlibs:
         
        seqfile = open("%s/Preprocess/out/lib%d.seq.btfilt"%(rundir,lib.id),'w')
        matefile = open("%s/Preprocess/out/lib%d.seq.mates"%(rundir,lib.id),'r')
        new_matefile = open("%s/Assemble/out/%s.lib%d.mappedmates"%(rundir,PREFIX,lib.id),'w')
        matedict[lib.id] = {}
        for line in matefile.xreadlines():
            line = line.replace("\n","")
            mate1, mate2 = line.split("\t")
            mate1 = mate1.replace("@","").replace(">","")
            mate2 = mate2.replace("@","").replace(">","")
            matedict[lib.id][mate2] = mate1
            #matedict[lib.id][mate1] = mate2
            read_lookup[readcnt] = mate1
            read_lookup[readcnt+1] = mate2
            readcnt += 2

        if bowtie_mapping == 1:
            #trim to 25bp
            trim = 0
            if trim:
                f1 = open("%s/Preprocess/out/lib%d.seq"%(rundir,lib.id))
                f2 = open("%s/Preprocess/out/lib%d.seq.trim"%(rundir,lib.id),'w')
                linecnt = 1
                for line in f1.xreadlines():
                    if linecnt % 2 == 0:
                        f2.write(line[0:25]+"\n")
                    else:
                        f2.write(line)
                    linecnt +=1
                f1.close()
                f2.close()
            if not os.path.exists("%s/Assemble/out/IDX.1.ebwt"%(rundir)):
                run_process("%s/bowtie-build %s/Assemble/out/%s.asm.contig %s/Assemble/out/IDX"%(BOWTIE, rundir,PREFIX,rundir))
            #run_process("%s/bowtie-build %s/Assemble/out/%s.asm.contig %s/Assemble/out/IDX"%(BOWTIE, rundir,PREFIX,rundir))
            if "bowtie" not in skipsteps and fasta:
                if trim:
                    run_process("%s/bowtie -p %d -f -v 1 -M 2 %s/Assemble/out/IDX %s/Preprocess/out/lib%d.seq.trim >& %s/Assemble/out/%s.bout"%(BOWTIE,threads,rundir,rundir,lib.id,rundir,PREFIX))
                else:
                    run_process("%s/bowtie -p %d -f -l 28 -M 2 %s/Assemble/out/IDX %s/Preprocess/out/lib%d.seq >& %s/Assemble/out/%s.bout"%(BOWTIE,threads,rundir,rundir,lib.id,rundir,PREFIX))
            elif "bowtie" not in skipsteps:
                if trim:
                    run_process("%s/bowtie  -p %d -v 1 -M 2 %s/Assemble/out/IDX %s/Preprocess/out/lib%d.seq.trim >& %s/Assemble/out/%s.bout"%(BOWTIE,threads,rundir,rundir,lib.id,rundir,PREFIX))
                else:
                    run_process("%s/bowtie  -p %d -l 28 -M 2 %s/Assemble/out/IDX %s/Preprocess/out/lib%d.seq >& %s/Assemble/out/%s.bout"%(BOWTIE,threads,rundir,rundir,lib.id,rundir,PREFIX))
            infile = open("%s/Assemble/out/%s.bout"%(rundir,PREFIX),'r')
            for line1 in infile.xreadlines():
                line1 = line1.replace("\n","")
                ldata = line1.split("\t")
                if len(ldata) < 6:
                    continue
                read, strand, contig, spos,read_seq, read_qual  = ldata[:6]
                read = read.split(" ")[0]
                epos = int(spos)+len(read_seq)
                try:
                    contigdict[contig].append([int(spos), int(spos)+epos, strand, read])
                except KeyError:
                    contigdict[contig] = [[int(spos),int(spos)+epos,strand,read]]
            
                seqdict[read] = read_seq
                seqfile.write(">%s\n%s\n"%(read,read_seq))
                seqfile.flush()

        else:
 
            #open soap ReadOnContig
            #some contigs are missing!
            infile = open("%s/Assemble/out/%s.asm.readOnContig"%(rundir,PREFIX),'r')
            #readID, ContigID, startpos, strand
            hdr = infile.readline()
            linecnt = 1
            for line in infile.xreadlines():
                if linecnt % 100000 == 0:
                    #print linecnt,
                    sys.stdout.flush()
                data = line.replace("\n","").split("\t")
                #print data
                if len(data) < 4:
                    continue
                contig = data[1]
                spos = int(data[2])
                if spos < 0:
                    spos = 0
                epos = spos+readlen
                strand = data[3]
                read = int(data[0])

                try:
                    contigdict[contig].append([int(spos), int(spos)+epos, strand, read_lookup[read]])
                except KeyError:
                    contigdict[contig] = [[int(spos),int(spos)+epos,strand,read_lookup[read]]]
                read_seq = "TEST"
            
                seqdict[read_lookup[read]] = read_seq
                linecnt +=1
        
    contig_data = contigfile.read()
    contig_data = contig_data.split(">")
    errfile = open("%s/Assemble/out/contigs_wo_location_info.txt"%(rundir),'w')
    new_ctgfile = open("%s/Assemble/out/%s.seq100.contig"%(rundir,PREFIX),'w')
    ctgcnt = 1
    ctgseq = 0
    ctgsizes = []
    n50_size = 0
    n50_mid = 955,000
    for item in contig_data:
        if item == '':
            continue

        item = item.split("\n",1)
        ref = item[0].split(" ")[0]
        ref = ref.replace("\n","")
        cseq = item[1].replace("\n","")
        ctgseq+=len(cseq)
        ctgsizes.append(len(cseq))
        i = 0
        cpos = 0
        width = 70
        cseq_fmt = ""
        while i+width < len(cseq):
            cseq_fmt += cseq[i:i+width]+"\n"
            i+= width
        cseq_fmt += cseq[i:]+"\n"
        ctgslen = len(item[1])
        #contigdict2[ref] = item[1]
        try:
            tigr_file.write("##%d %d %d bases, 00000000 checksum.\n"%(ctgcnt,len(contigdict[ref]), len(item[1])))
            tigr_file.flush()
        except KeyError:
            #print "oops, not in mapping file\n"
            errfile.write("%s\n"%ref)
            continue
        new_ctgfile.write(">%d\n%s"%(ctgcnt,cseq_fmt))#item[1]))
        ctgcnt +=1
        tigr_file.write(cseq_fmt)#item[1])
        contigdict[ref].sort()
        #print contigdict[ref]
        for read in contigdict[ref]:
            
            try:
                #if read[0] <= 500 and ctgslen - (int(read[1])) <= 500:
                matectgdict[read[-1]] = ref
                mateotdict[read[-1]] = read[2]
            except KeyError:
                pass
            if read[2] == "-":
                tigr_file.write("#%s(%d) [RC] %d bases, 00000000 checksum. {%d 1} <%d %s>\n"%(read[-1],read[0]-1, len(seqdict[read[-1]]), len(read[-1]), read[0], read[1]))
            else:
                tigr_file.write("#%s(%d) [] %d bases, 00000000 checksum. {1 %d} <%d %s>\n"%(read[-1],read[0]-1, len(seqdict[read[-1]]), len(read[-1]), read[0], read[1]))
            tigr_file.write(seqdict[read[-1]]+"\n")

   

    for lib in readlibs:
        new_matefile.write("library\t%d\t%d\t%d\n"%(lib.id,lib.mmin,lib.mmax))
    for lib in readlibs:
        linked_contigs = {}
        for mate in matedict[lib.id].keys():
            new_matefile.write("%s\t%s\t%d\n"%(mate,matedict[lib.id][mate],lib.id))
            new_matefile.flush()
            continue
    
    tigr_file.close()
        
    
        
def LCS(S1, S2):
    M = [[0]*(1+len(S2)) for i in xrange(1+len(S1))]
    longest, x_longest = 0, 0
    for x in xrange(1,1+len(S1)):
        for y in xrange(1,1+len(S2)):
            if S1[x-1] == S2[y-1]:
                M[x][y] = M[x-1][y-1] + 1
                if M[x][y]>longest:
                    longest = M[x][y]
                    x_longest  = x
            else:
                M[x][y] = 0
    return S1[x_longest-longest: x_longest]


def start_http(server_class=BaseHTTPServer.HTTPServer,
        handler_class=BaseHTTPServer.BaseHTTPRequestHandler):
    #pid = os.fork()
    server_address = ('localhost', 8111)
    httpd = server_class(server_address, handler_class)
    httpd.serve_forever()
    #return pid
def validate_run(dir):
    run_process("./%s/run.sh"%(dir))
    #check to see if all output files (listed in README) were generated
    readme = open("./%s/README"%(dir),'r')
    outf = 0
    for line in readme:
        if "[OUTPUT]" in line:
            outf = 1
        elif "[RUN]" in line:
            print "all output files successfully generated!"
            outf = 0
        if outf:
            if "." in line:
                ff = line.split(".")[-1].split(",")
                for file in ff:
                    if not os.path.exists("./%s/out/%s"%(dir,file)):
                        print "%s failed"%(dir)
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

#if asm == "soap":
readpaths = []
filtreadpaths = []
for lib in readlibs:
    for read in lib.reads:
        readpaths.append("%s/Preprocess/in/"%(rundir)+read.fname)
        filtreadpaths.append("%s/Preprocess/out/"%(rundir)+read.fname)
   
if "Preprocess" in forcesteps:
    for path in readpaths:
        run_process("touch %s"%path)

#@transform(readpaths,["%s/Preprocess/out/all.seq"%(rundir),"%s/Preprocess/out/all.seq.mates"%(rundir)])
@files(readpaths,filtreadpaths)
def Preprocess(input,output):
   #move input files into Preprocess ./in dir
   #output will either be split fastq files in out, or AMOS bank
   if "Preprocess" in skipsteps or "preprocess" in skipsteps:
       for lib in readlibs:
           for read in lib.reads:
               run_process("ln -t ./%s/Preprocess/out/ -s ../../Preprocess/in/%s"%(rundir,read.fname))
       return 0
   if filter == True:
       #print "filtering.."
     
       #for reads+libs
       cnt = 1
       for lib in readlibs:
           print lib.interleaved, lib.mated, lib.format
           for read in lib.reads:
               print read.interleaved, read.mated, read.format
               if not read.filtered and read.format == "fastq" and read.mated and read.interleaved:
                   #this means we have this entire lib in one file
                   #parse out paired record (8 lines), rename header to be filename + "/1" or "/2", and remove reads with N
                   rf = open(read.path,'r')
                   npath = read.path.replace("/in/","/out/")
                   #readpath,base = os.path.split(npath)
                   #newpath = readpath+"lib%d"%(lib.id)
                   wf = open(npath,'w')
                   #wf = open(read.path.replace("/in/","/out/"),'w')
                   #wf = open(readpath+"lib%d"%(lib.id),'w')
                   start = 1
                   rcnt = 0
                   recordcnt = 0
                   record = []
                   shdr = ""
                   for line in rf.xreadlines():
                       if start:
                           s1hdr = line
                           record.append(line)
                           start = 0
                           rcnt =1

                       else:
                           if rcnt == 7:
                               #end of record
                               record.append(line)
                               rcnt +=1
                               if len(record) != 8:
                                   #something went wrong
                                   continue
                               rseq = string.upper(record[0]+record[5])                               
                               if "N" in rseq:
                                   #skip both, dont' want Ns
                                   continue
                               #update hdrs to be filename /1 or /2
                               recordcnt +=1
                               hdr = lib.sid+"r"+str(recordcnt)+"/"
                               wf.writelines("@"+hdr+"1\n")
                               wf.writelines(record[1])
                               wf.writelines("+"+hdr+"1\n")
                               wf.writelines(record[3])
                               wf.writelines("@"+hdr+"2\n")
                               wf.writelines(record[5])
                               wf.writelines("+"+hdr+"2\n")
                               wf.writelines(record[7])
                           elif rcnt % 4 == 0:
                               s2hdr = line
                               rlcs = LCS(s1hdr,s2hdr)
                               #these should almost identical
                               if float(len(rlcs))/float(len(s1hdr)) < 0.9:
                                   #missing record somewhere, start over with this one
                                   s1hdr = line
                                   record = [line]
                                   start = 0
                                   rcnt = 1
                               else:
                                   record.append(line)
                                   rcnt +=1
                           elif rcnt % 2 == 0:
                               #quality hdr
                               record .append(line)
                               rcnt +=1
                           else:
                               record .append(line)
                               rcnt +=1
                   #update to new path
                   read.path = read.path.replace("/in/","/out/")            
                   #read.fname = "lib%d"%(lib.id)
                   read.filtered = True
               elif not read.filtered and read.format == "fastq" and read.mated and not read.interleaved:
                   readpair = lib.getPair(read.id)
                   if readpair == -1:
                       #not interleaved and mated, yet do not have 2nd file..
                       continue
                   rf1 = open(read.path,'r')
                   wf1 = open(read.path.replace("/in/","/out/"),'w')
                   rf2 = open(readpair.path,'r')
                   wf2 = open(readpair.path.replace("/in/","/out/"),'w')
                   recordcnt = 0
                   while 1:                   
                       rs1 = rf1.readline()
                       rs2 = rf1.readline()
                       rs3 = rf1.readline()
                       rs4 = rf1.readline()
                       if rs1 == "" or rs2 == "" or rs3 == "" or rs4 == "":
                           #EOF or something went wrong, break
                           break 
                       rseq = string.upper(rs2)                               
                       if "N" in rseq:
                           continue
                       rp1 = rf2.readline()
                       rp2 = rf2.readline()
                       rp3 = rf2.readline()
                       rp4 = rf2.readline()
                       if rp1 == "" or rp2 == "" or rp3 == "" or rp4 == "":
                           #EOF or something went wrong, break
                           break 
                       rseq = string.upper(rp2)                               
                       if "N" in rseq:
                           continue

                       #rlcs = LCS(rs1,rp1)
                       if 0:#float(len(rlcs))/float(len(rs1)) < 0.9:
                           #not aligned! something needs to be removed?
                           #go on to next
                           continue
                       else:
                           #record.append(line)
                           #rcnt +=1
                           recordcnt +=1
                           hdr = lib.sid+"r"+str(recordcnt)+"/"
                           wf1.writelines("@"+hdr+"1\n")
                           wf1.writelines(rs2)
                           wf1.writelines("+"+hdr+"1\n")
                           wf1.writelines(rs4)
                           wf2.writelines("@"+hdr+"2\n")
                           wf2.writelines(rp2)
                           wf2.writelines("+"+hdr+"2\n")
                           wf2.writelines(rp4)
                           wf1.flush()
                           wf2.flush()
                   readpair.filtered = True
                   read.filtered = True
                   read.path = read.path.replace("/in/","/out/")
                   readpair.path = readpair.path.replace("/in/","/out/")
               elif not read.filtered and read.format == "fastq" and not read.mated:
                   #this is easy, just throw out reads with Ns
                   rf = open(read.path,'r')
                   wf = open(read.path.replace("/in/","/out/"),'w')
                   while 1:
                       rs1 = rf.readline()
                       rs2 = rf.readline()
                       rs3 = rf.readline()
                       rs4 = rf.readline()
                       if rs1 == "" or rs2 == "" or rs3 == "" or rs4 == "":
                           #EOF or something went wrong, break
                           break 
                       rseq = string.upper(rs2)                               
                       if "N" in rseq:
                           continue
                       wf.writelines(rs1)
                       wf.writelines(rs2)
                       wf.writelines(rs3)
                       wf.writelines(rs4)
                   read.path = read.path.replace("/in/","/out/")
               elif not read.filtered and read.format == "fasta" and read.mated and read.interleaved:
                   #this means we have this entire lib in one file
                   #parse out paired record (4 lines), rename header to be filename + "/1" or "/2", and remove reads with N
                   rf = open(read.path,'r')
                   npath = read.path.replace("/in/","/out/")
                   print npath
                   #readpath,base = os.path.split(npath)
                   #newpath = readpath+"lib%d"%(lib.id)
                   wf = open(npath,'w')
                   #wf = open(read.path.replace("/in/","/out/"),'w')
                   start = 1
                   rcnt = 0
                   recordcnt = 0
                   record = []
                   shdr = ""
                   reads = rf.read().split(">")[1:]
                   if len(reads) % 4 != 0:
                       print "Read file corrupted, please fix and restart!"
                       sys.exit(1)

                   prevok = False
                   first = True
                   second = False
                   prevseq = ""
                   readcnt = 1
                   for rd in reads:
                       if first:
                           hdr,seq = rd.split("\n",1)
                           if "N" in string.upper(seq) or len(seq) < 2:
                               prevok = False
                           else: 
                               prevok = True
                               prevseq = seq

                           second = True
                           first = False
                       elif second:
                           hdr,seq = rd.split("\n",1)
                           if "N" in string.upper(seq) or len(seq) < 2:
                               pass
                           elif prevok:
                               hdr = lib.sid+"r"+str(readcnt)+"/"
                               wf.writelines(">"+hdr+"1\n")
                               wf.writelines(prevseq)
                               wf.writelines(">"+hdr+"2\n")
                               wf.writelines(seq)
                               readcnt +=1
                           second = False
                           first = True

                   #update to new path
                   read.path = read.path.replace("/in/","/out/")            

                   #read.path = newpath#read.path.replace("/in/","/out/")            
                   #read.fname = "lib%d"%(lib.id)
               elif not read.filtered and read.format == "fasta" and read.mated and not read.interleaved:
                   readpair = lib.getPair(read.id)
                   if readpair == -1:
                       #not interleaved and mated, yet do not have 2nd file..
                       continue
                   rf1 = open(read.path,'r')
                   wf1 = open(read.path.replace("/in/","/out/"),'w')
                   rf2 = open(readpair.path,'r')
                   wf2 = open(readpair.path.replace("/in/","/out/"),'w')
                   recordcnt = 0
                   while 1:                   
                       rs1 = rf1.readline()
                       rs2 = rf1.readline()

                       if rs1 == "" or rs2 == "":
                           #EOF or something went wrong, break
                           break 
                       rseq = string.upper(rs2)                               
                       if "N" in rseq:
                           continue
                       rp1 = rf2.readline()
                       rp2 = rf2.readline()

                       if rp1 == "" or rp2 == "":
                           #EOF or something went wrong, break
                           break 
                       rseq = string.upper(rp2)                               
                       if "N" in rseq:
                           continue

                       rlcs = LCS(rs1,rp1)
                       if 0:#float(len(rlcs))/float(len(rs1)) < 0.9:
                           #not aligned! something needs to be removed?
                           #go on to next
                           continue
                       else:
                           #record.append(line)
                           #rcnt +=1
                           recordcnt +=1
                           hdr = lib.sid+"r"+str(recordcnt)+"/"
                           wf1.writeline(">"+hdr+"1\n")
                           wf1.writeline(rs2)
                           wf2.writeline(">"+hdr+"2\n")
                           wf2.writeline(rp2)

                   readpair.filtered = True
                   read.filtered = True
                   read.path = read.path.replace("/in/","/out/")
                   readpair.path = readpair.path.replace("/in/","/out/")
               elif not read.filtered and read.format == "fasta" and not read.mated:
                   #easiest case, check for Ns
                   rf = open(read.path,'r')
                   wf = open(read.path.replace("/in/","/out/"),'w')
                   while 1:
                       rs1 = rf.readline()
                       rs2 = rf.readline()
                       if rs1 == "" or rs2 == "":
                           #EOF or something went wrong, break
                           break
                       rseq = string.upper(rs2)                               
                       if "N" in rseq:
                           continue
                       wf.writelines(rs1)
                       wf.writelines(rs2)
                   read.path = read.path.replace("/in/","/out/")
           cnt +=1
   else:
       for lib in readlibs:
           for read in lib.reads:
               run_process("ln -t ./%s/Preprocess/out/ -s ../../Preprocess/in/%s"%(rundir,read.fname))
   #PUNT HERE
   for lib in readlibs:
      if 1:
           #this means interleaved, single file
           if lib.format == "sff":
               # generate the fasta files from the sff file
               sffToCACmd = "%s/sffToCA -clear 454 -clear discard-n -trim chop -libraryname lib%d -output %s/Preprocess/out/%s"%(CA, lib.id,rundir, PREFIX)
               if (read.mated == True):
                   run_process("%s -linker %s -insertsize %d %d %s"%(sffToCACmd, read.linkerType, lib.mean, lib.stdev, read.path))
               else:
                   run_process("%s %s"%(sffToCACmd, read.path))
               run_process("%s/gatekeeper -T -F -o %s/Preprocess/out/%s.gkpStore %s/Preprocess/out/%s.frg"%(CA, rundir, PREFIX, rundir, PREFIX))
               run_process("%s/gatekeeper -dumpnewbler %s/Preprocess/out/%s %s/Preprocess/out/%s.gkpStore"%(CA, rundir, PREFIX, rundir, PREFIX))
               run_process("%s/gatekeeper -dumplibraries -tabular %s/Preprocess/out/%s.gkpStore |awk '{if (match($3, \"U\") == 0 && match($1, \"UID\") == 0) print \"library\t\"$1\"\t\"$4-$5*3\"\t\"$4+$5*3}' > %s/Preprocess/out/lib%d.seq.mates"%(CA, rundir, PREFIX, rundir,lib.id))
               run_process("%s/gatekeeper -dumpfragments -tabular %s/Preprocess/out/%s.gkpStore|awk '{if ($3 != 0 && match($1, \"UID\")==0 && $1 < $3) print $1\"\t\"$3\"\t\"$5}' >> %s/Preprocess/out/lib%d.seq.mates"%(CA, rundir, PREFIX, rundir,lib.id))
               run_process("unlink %s/Preprocess/out/lib%d.seq"%(rundir,lib.id))
               run_process("ln -s  ../../Preprocess/out/%s.fna %s/Preprocess/out/lib%d.seq"%(PREFIX, rundir,lib.id))
               run_process("ln -s ../../Preprocess/out/%s.fna.qual %s/Preprocess/out/lib%d.seq.qual"%(PREFIX, rundir,lib.id))
               run_process("rm -rf %s/Preproces/out/%s.gkpStore"%(rundir, PREFIX))
               run_process("unlink %s/Preprocess/out/%s.frg"%(rundir, PREFIX))
           elif lib.format == "fasta" and not lib.mated:
               run_process("ln -s  ../../Preprocess/in/%s %s/Preprocess/out/lib%d.seq"%(lib.id,lib.f1.fname, rundir))
               run_process("ln -s ../../Preprocess/in/%s.qual %s/Preprocess/out/lib%d.seq.qual"%(lib.id,lib.f1.fname, rundir))
               run_process("touch %s/Preprocess/out/lib%d.seq.mates"%(lib.id,rundir))
           elif format == "fasta" and mated and not interleaved:
               #FIXME, make me faster!
               run_process("perl %s/perl/shuffleSequences_fasta.pl  %s/Preprocess/out/%s %s/Preprocess/out/%s %s/Preprocess/out/lib%d.seq"%(METAMOS_UTILS,rundir,lib.f1.fname, rundir,lib.f2.fname,rundir,lib.id))
               run_process("python %s/python/extract_mates_from_fasta.py %s/Preprocess/out/lib%d.seq"%(METAMOS_UTILS,rundir,lib.id))
               run_process("unlink ./%s/Preprocess/out/lib%d.seq.mates"%(rundir, lib.id))
               run_process("ln -t ./%s/Preprocess/out/ -s ../../Preprocess/in/lib%d.seq.mates"%(rundir,lib.id))
           elif format == "fastq" and mated and not interleaved:
               #extract mates from fastq
               run_process("perl %s/perl/shuffleSequences_fastq.pl  %s/Preprocess/out/%s %s/Preprocess/out/%s %s/Preprocess/out/lib%d.seq"%(METAMOS_UTILS,rundir,lib.f1.fname, rundir,lib.f2.fname,rundir,lib.id))
               run_process("python %s/python/extract_mates_from_fastq.py %s/Preprocess/out/lib%d.seq"%(METAMOS_UTILS,rundir,lib.id))
           elif mated and interleaved:
               os.system("cp %s/Preprocess/out/%s %s/Preprocess/out/lib%d.seq"%(rundir,lib.f1.fname,rundir,lib.id))
               if format == "fastq":
                   run_process("python %s/python/extract_mates_from_fastq.py %s/Preprocess/out/lib%d.seq"%(METAMOS_UTILS,rundir,lib.id))
               else:
                   run_process("python %s/python/extract_mates_from_fasta.py %s/Preprocess/out/lib%d.seq"%(METAMOS_UTILS,rundir,lib.id))
           #update_soap_config()
           elif asm == "ca":
               #useful for 454, need to get SFF to FRG?
               #/fs/wrenhomes/sergek/wgs-assembler/Linux-amd64/bin/sffToCA
               pass
           elif asm == "amos":
               #call toAmos_new              
               pass

asmfiles = []
#if asm == "soap"

for lib in readlibs:
    if "Assemble" in forcesteps:
        run_process("touch %s/Preprocess/out/lib%d.seq"%(rundir,lib.id))
    asmfiles.append("%s/Preprocess/out/lib%d.seq"%(rundir,lib.id))

@files(asmfiles,["%s/Assemble/out/%s.asm.contig"%(rundir,PREFIX)])
#@posttask(create_symlink,touch_file("completed.flag"))
@follows(Preprocess)
def Assemble(input,output):
   #pick assembler
   if "Assemble" in skipsteps or "assemble" in skipsteps:
      return 0
   if asm == "soap":
      #open & update config
      soapf = open("%s/config.txt"%(rundir),'r')
      soapd = soapf.read()
      soapf.close()
      cnt = 1
      libno = 1
      #print libs
      for lib in readlibs:
          if (lib.format == "fastq" or lib.format == "fasta")  and lib.mated and not lib.interleaved:
              soapd = soapd.replace("LIB%dQ1REPLACE"%(lib.id),"%s/Preprocess/out/%s"%(rundir,lib.f1.fname))
              soapd = soapd.replace("LIB%dQ2REPLACE"%(lib.id),"%s/Preprocess/out/%s"%(rundir,lib.f2.fname))

          elif lib.format == "fastq"  and lib.mated and lib.interleaved:
              #this is NOT supported by SOAP, make sure files are split into two..
              #need to update lib.f2 path
              run_process("perl spilt_fastq.pl %s/Preprocess/out/%s %s/Preprocess/out/%s %s/Preprocess/out/%s"%(lib.f1.fname,lib.f1.fname,lib.f2.fname))
              soapd = soapd.replace("LIB%dQ1REPLACE"%(lib.id),"%s/Preprocess/out/%s"%(rundir,lib.f1.fname))
              soapd = soapd.replace("LIB%dQ2REPLACE"%(lib.id),"%s/Preprocess/out/%s"%(rundir,lib.f2.fname))

          elif format == "fasta"  and mated and interleaved:
              soapd = soapd.replace("LIB%dQ1REPLACE"%(lib.id),"%s/Preprocess/out/%s"%(rundir,lib.f1.fname))
          else:
              soapd = soapd.replace("LIB%dQ1REPLACE"%(lib.id),"%s"%(lib.f1))#frg))
      #cnt +=1
      soapw = open("%s/soapconfig.txt"%(rundir),'w')
      soapw.write(soapd)
      soapw.close()
      print "Running SOAPdenovo on input reads..."
      #start stopwatch
      run_process("%s/soap63 all  -D -d -R -p %d -K %d -M 3 -L 300 -s %s/soapconfig.txt -o %s/Assemble/out/%s.asm"%(SOAP, threads, kmer, rundir,rundir,PREFIX))#SOAPdenovo config.txt

      #run_process("%s/SOAPdenovo-63mer all -D -d -R -p %d -K %d -M 3 -s %s/soapconfig.txt -o %s/Assemble/out/%s.asm"%(SOAP, threads, kmer, rundir,rundir,PREFIX))#SOAPdenovo config.txt
      #run_process("%s/SOAPdenovo-31mer all  -D 3 -d 2 -R -p %d -M 3 -K %d -s %s/soapconfig.txt -o %s/Assemble/out/%s.asm"%(SOAP, threads, kmer, rundir,rundir,PREFIX))#SOAPdenovo config.txt

      #run_process("ln -s %s/Assemble/out/%s.asm.contig ./%s/FindORFS/in/%s.asm.contig"%(rundir,PREFIX, rundir, PREFIX)) 

      #if OK, convert output to AMOS

   elif asm == "newbler":
      run_process("%s/newAssembly -force %s/Assemble/out"%(NEWBLER, rundir));
      for lib in readlibs:
          if lib.format == "fasta"  and lib.interleaved:
              run_process("%s/addRun %s/Assemble/out %s/Preprocess/out/lib%d.seq"%(NEWBLER, rundir, rundir,lib.id));
          elif lib.format == "fastq":
              print "WARNING!! FASTQ + Newbler only supported in version 2.6+, Newbler may fail"
              sys.exit(1)
      newblerCmd = "%s%srunProject"%(NEWBLER, os.sep)
      # read spec file to input to newbler parameters
      newblerCmd += getProgramParams("newbler.spec", "", "-")
      run_process("%s -cpu %d %s/Assemble/out"%(newblerCmd,threads,rundir));

      # convert to AMOS
      run_process("%s/toAmos -o %s/Assemble/out/%s.mates.afg -m %s/Preprocess/out/all.seq.mates -ace %s/Assemble/out/assembly/454Contigs.ace"%(AMOS,rundir, PREFIX, rundir, rundir));
      # get info on EID/IIDs for contigs
      run_process("cat %s/Assemble/out/%s.mates.afg | grep -A 3 \"{CTG\" |awk '{if (match($1, \"iid\") != 0) {IID = $1} else if (match($1, \"eid\") != 0) {print $1\" \"IID; } }'|sed s/eid://g |sed s/iid://g > %s/Assemble/out/454eidToIID"%(rundir, PREFIX, rundir))
      run_process("java -cp %s convert454GraphToCTL %s/Assemble/out/454eidToIID %s/Assemble/out/assembly/454ContigGraph.txt > %s/Assemble/out/%s.graph.cte"%(METAMOS_JAVA, rundir, rundir, rundir, PREFIX));
      run_process("cat %s/Assemble/out/%s.mates.afg %s/Assemble/out/%s.graph.cte > %s/Assemble/out/%s.afg"%(rundir, PREFIX, rundir, PREFIX, rundir, PREFIX))
    
      # make symlink for subsequent steps
      run_process("rm %s/Assemble/out/%s.asm.contig"%(rundir, PREFIX));
      run_process("ln -s ../../Assemble/out/assembly/454AllContigs.fna %s/Assemble/out/%s.asm.contig"%(rundir, PREFIX))
      if mated == True:
         run_process("ln -s ../../Assemble/out/assembly/454Scaffolds.fna %s/Assemble/out/%s.asm.scafSeq"%(rundir, PREFIX))
      else:
         run_process("ln -s ../../Assemble/out/assembly/454AllContigs.fna %s/Assemble/out/%s.asm.scafSeq"%(rundir, PREFIX))

   elif asm == "amos":
      run_process("Minimus ./Preprocess/out/bank")
   elif asm == "CA":
      #runCA script
      frglist = ""
      for lib in readlibs:
          for read in lib.reads:
              if read.format == "fastq":
                  run_process("%s/fastqToCA -insertsize %d %d -libraryname %s -t illumina -innie -fastq %s/Preprocess/in/%s"%(CA, lib.mean,lib.stdev, read.path, rundir,PREFIX))
              elif read.format == "fasta":
                  run_process("%s/convert-fasta-to-v2.pl -l %s -mean %d -stddev %d -s %s/Preprocess/in/%s -q %s/Preprocess/in/%s.qual -m matepairids %s/Preprocess/out/%s.mateids %s"%(CA,read.path, lib.mean, lib.stdev, rundir, read.fname, rundir, read.fname, rundir, read.fname,fff))
              frglist += "%s.frg"%(lib)
      run_process("%s/runCA -p asm -d %s/Assemble/out/ -s %/config/asm.spec %s"%(CA,rundir,METAMOS_UTILSPREFIX,frglist))
      #convert CA to AMOS
      run_process("%s/gatekeeper -dumpfrg -allreads -format2 asm.gkpStore > asm.frg bzip2 asm.frg"%(CA))
      run_process("%s/terminator -g asm.gkpStore -t asm.tigStore/ 2 -o asm bzip2 asm.asm"%(CA))
      run_process("%s/toAmos_new -a asm.asm.bz2 -f asm.frg.bz2 -b asm.bnk -U "%(AMOS))
   #stop here, for now
   #sys.exit(0)
   #check if sucessfully completed   

if "FindORFS" in forcesteps:
   run_process("touch %s/Assemble/out/%s.asm.contig"%(rundir,PREFIX))
@follows(Assemble)
@files("%s/Assemble/out/%s.asm.contig"%(rundir,PREFIX),"%s/FindORFS/out/%s.faa"%(rundir,PREFIX))
def FindORFS(input,output):
   if "FindORFS" in skipsteps:
      run_process("touch %s/FindRepeats/in/%s.fna"%(rundir, PREFIX))
      run_process("touch %s/FindORFS/out/%s.faa"%(rundir, PREFIX))
      return 0

   if asm == "soap":
         
       #if not os.path.exists("%s/Assemble/out/%s.asm.scafSeq.contigs"%(rundir,PREFIX)):
       #    run_process("python %s/python/extract_soap_contigs.py %s/Assemble/out/%s.asm.scafSeq"%(METAMOS_UTILS,rundir,PREFIX))
       #run_process("unlink %s/FindORFS/in/%s.asm.scafSeq.contigs"%(rundir,PREFIX))
       #run_process("unlink %s/FindORFS/in/%s.asm.contig"%(rundir,PREFIX))
       #run_process("ln -t ./%s/FindORFS/in/ -s ../../Assemble/out/%s.asm.scafSeq.contigs"%(rundir,PREFIX))
       #run_process("mv %s/FindORFS/in/%s.asm.scafSeq.contigs  %s/FindORFS/in/%s.asm.contig"%(rundir,PREFIX,rundir,PREFIX))
       #try using contigs instead of contigs extracted from scaffolds
       run_process("cp %s/Assemble/out/%s.asm.contig  %s/FindORFS/in/%s.asm.contig"%(rundir,PREFIX,rundir,PREFIX))
   else:

       run_process("unlink %s/FindORFS/in/%s.asm.contig"%(rundir,PREFIX))
       run_process("ln -t ./%s/FindORFS/in/ -s ../../Assemble/out/%s.asm.contig"%(rundir,PREFIX))


   #run_process("ln -t ./%s/FindORFS/in/ -s ../../Assemble/out/%s.asm.scafSeq.contigs"%(rundir,PREFIX))
   run_process("%s/gmhmmp -o %s/FindORFS/out/%s.orfs -m %s/config/MetaGeneMark_v1.mod -d -a %s/FindORFS/in/%s.asm.contig"%(GMHMMP,rundir,PREFIX,METAMOS_UTILS,rundir,PREFIX))
   parse_genemarkout("%s/FindORFS/out/%s.orfs"%(rundir,PREFIX))
   run_process("unlink ./%s/Annotate/in/%s.faa"%(rundir,PREFIX))
   run_process("unlink ./%s/Annotate/in/%s.fna"%(rundir,PREFIX))
   run_process("unlink ./%s/FindRepeats/in/%s.fna"%(rundir,PREFIX))
   run_process("ln -t ./%s/Annotate/in/ -s ../../FindORFS/out/%s.faa"%(rundir,PREFIX))
   run_process("ln -t ./%s/FindRepeats/in/ -s ../../FindORFS/out/%s.fna"%(rundir,PREFIX))

@follows(FindORFS)
@files("%s/FindRepeats/in/%s.fna"%(rundir,PREFIX),"%s/FindRepeats/out/%s.repeats"%(rundir,PREFIX))
def FindRepeats(input,output):
   if "FindORFS" in skipsteps or "FindRepeats" in skipsteps:
     return 0

   run_process("python %s/python/getContigRepeats.py  %s/FindRepeats/in/%s.fna %s/FindRepeats/out/%s.repeats"%(METAMOS_UTILS,rundir,PREFIX,rundir,PREFIX))


#@follows(FindRepeats)
@files("%s/Annotate/in/%s.faa"%(rundir,PREFIX),"%s/Annotate/out/%s.hits"%(rundir,PREFIX))
def Annotate(input,output):
   #annotate contigs > 1000bp with FCP
   #lets start by annotating ORFs with phmmer
   if cls == "phmmer":

       run_process("phmmer --cpu %d --F1 0.01 --F2 0.0001 --F3 0.000001 -E 0.01 -o %s/Annotate/out/%s.phm.out --tblout %s/Annotate/out/%s.phm.tbl --notextw %s/Annotate/in/%s.faa %s/DB/allprots.faa"%(threads,rundir,PREFIX,rundir,PREFIX,rundir,PREFIX,METAMOS_UTILS))
       parse_phmmerout("%s/Annotate/out/%s.phm.tbl"%(rundir,PREFIX))
       run_process("mv %s/Annotate/out/%s.phm.tbl  %s/Annotate/out/%s.hits"%(rundir,PREFIX,rundir,PREFIX))
       #run_process("mv %s/Annotate/out/%s.phm.tbl  %s/Annotate/out/%s.annotate"%(rundir,PREFIX,rundir,PREFIX))
   elif cls == "blast":
       run_process("blastall -v 1 -b 1 -a %d -p blastp -m 8 -e 0.00001 -i %s/Annotate/in/%s.faa -d %s/DB/new_all_complete_bacteria.faa -o %s/Annotate/out/%s.blastout"%(threads, rundir,PREFIX,METAMOS_UTILS,rundir,PREFIX))
       run_process("mv %s/Annotate/out/%s.blastout  %s/Annotate/out/%s.hits"%(rundir,PREFIX,rundir,PREFIX))
   elif cls == "fcp":
       print "FCP not yet supported.. stay tuned!"



if "Metaphyler" in forcesteps:
   run_process("touch %s/FindORFS/out/%s.faa"%(rundir,PREFIX))
   #run_process("rm %s/Metaphyler/out/%s.classify.txt"%(rundir,PREFIX))

@follows(FindORFS)
@files("%s/FindORFS/out/%s.faa"%(rundir,PREFIX),"%s/Metaphyler/out/%s.classify.txt"%(rundir,PREFIX))
def Metaphyler(input,output):
   if "FindORFS" in skipsteps or "Metaphyler" in skipsteps:
      return 0

   run_process("unlink ./%s/Metaphyler/in/%s.contig.cvg"%(rundir,PREFIX))
   run_process("unlink ./%s/Metaphyler/in/%s.faa"%(rundir,PREFIX))
   run_process("ln -t ./%s/Metaphyler/in/ -s ../../FindORFS/out/%s.contig.cvg"%(rundir,PREFIX))
   run_process("ln -t ./%s/Metaphyler/in/ -s ../../FindORFS/out/%s.faa"%(rundir,PREFIX))
   blastfile = PREFIX+".blastx"
   run_process("formatdb  -p T -i %s/DB/markers.pfasta"%(METAMOS_UTILS))
   #run_process("perl %s/perl/runblast.pl  %s/Metaphyler/in/%s.faa %s/Metaphyler/out/%s.blastx %s/DB/markers.fna"%(METAMOS_UTILS,rundir,PREFIX, rundir,PREFIX,METAMOS_UTILS))

   run_process("%s/cpp/blastall -p blastp -i %s/Metaphyler/in/%s.faa -d %s/DB/markers.pfasta -m8 -b10 -v10 -a %s -o %s/Metaphyler/out/%s.blastp"%(METAMOS_UTILS,rundir,PREFIX,METAMOS_UTILS,threads,rundir,PREFIX))

   run_process("perl %s/perl/metaphyler_contigs.pl %s/Metaphyler/out/%s.blastp %s %s/Metaphyler/in/%s.contig.cvg %s/Metaphyler/out %s"%(METAMOS_UTILS,rundir,PREFIX,PREFIX,rundir,PREFIX,rundir,METAMOS_UTILS))

   #run_process("perl %s/perl/metaphyler_contigs.pl %s/Metaphyler/out/%s.blastx %s/Metaphyler/out/%s %s/Metaphyler/in/%s.contig.cvg"%(METAMOS_UTILS,rundir,PREFIX, rundir, PREFIX, rundir, PREFIX))

   
if "Scaffold" in forcesteps:
    run_process("touch %s/Assemble/out/%s.asm.contig"%(rundir,PREFIX))
@follows(Metaphyler)
@files(["%s/Assemble/out/%s.asm.contig"%(rundir,PREFIX)],"%s/Scaffold/out/%s.scaffolds.final"%(rundir,PREFIX))
def Scaffold(input,output):
   # check if we need to do scaffolding
   numMates = 0
   if asm == "newbler":
      p = subprocess.Popen("cat %s/Assemble/out/%s.graph.cte |grep \"{CTL\" |wc -l"%(rundir, PREFIX), stdin=None, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      (checkStdout, checkStderr) = p.communicate()
      numMates = int(checkStdout.strip())

   if mated == False and numMates == 0:
      print "No mate pair info available for scaffolding, skipping"
      return 0

   if asm == "soap":
       for lib in readlibs:
           #run_process("ln -t ./%s/Metaphyler/in/ -s ../../FindORFS/out/%s.faa"%(rundir,PREFIX))
           run_process("rm -rf %s/Scaffold/in/%s.bnk"%(rundir,PREFIX))
           if lib.format == "fasta":
               map2contig(1)
               #PUNT: here, add toAmos_new call for each lib
               run_process("%s/toAmos_new -s %s/Preprocess/out/lib%d.seq -m %s/Assemble/out/%s.lib%d.mappedmates -b %s/Scaffold/in/%s.bnk "%(AMOS,rundir,lib.id,rundir, PREFIX,lib.id,rundir,PREFIX))

           elif format == "fastq":
               #if "bowtie" not in skipsteps:
               map2contig(0)
               run_process("%s/toAmos_new -Q %s/Preprocess/out/lib%d.seq -m %s/Assemble/out/%s.lib%d.mappedmates -b %s/Scaffold/in/%s.bnk "%(AMOS,rundir,lib.id,rundir,PREFIX, lib.id,rundir,PREFIX))

       run_process("%s/toAmos_new -c %s/Assemble/out/%s.asm.tigr -b %s/Scaffold/in/%s.bnk "%(AMOS,rundir,PREFIX,rundir,PREFIX))


   elif asm == "newbler":
      run_process("rm -rf %s/Scaffold/in/%s.bnk"%(rundir, PREFIX))
      # build the bank for amos
      run_process("%s/bank-transact -b %s/Scaffold/in/%s.bnk -c -m %s/Assemble/out/%s.afg"%(AMOS,rundir, PREFIX, rundir, PREFIX));

   #calls to Bambus2, goBambus2 script
   # first, parse the parameters
   markRepeatParams = getProgramParams("bambus.spec", "MarkRepeats", "-")
   orientContigParams = getProgramParams("bambus.spec", "OrientContigs", "-")

   run_process("%s/clk -b %s/Scaffold/in/%s.bnk"%(AMOS,rundir,PREFIX))
   run_process("%s/Bundler -b %s/Scaffold/in/%s.bnk"%(AMOS,rundir,PREFIX))
   run_process("%s/MarkRepeats %s -b %s/Scaffold/in/%s.bnk > %s/Scaffold/in/%s.reps"%(AMOS,markRepeatParams,rundir,PREFIX,rundir,PREFIX))
   run_process("%s/OrientContigs %s -b %s/Scaffold/in/%s.bnk -repeats %s/Scaffold/in/%s.reps "%(AMOS,orientContigParams,rundir,PREFIX, rundir, PREFIX))

   # output results
   run_process("%s/bank2fasta  -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.contigs"%(AMOS,rundir,PREFIX,rundir,PREFIX))
   run_process("%s/OutputMotifs -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.motifs"%(AMOS,rundir,PREFIX,rundir,PREFIX))
   run_process("%s/OutputResults -b %s/Scaffold/in/%s.bnk -p %s/Scaffold/out/%s "%(AMOS,rundir,PREFIX,rundir,PREFIX))
   run_process("%s/OutputScaffolds -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.scaffolds.final"%(AMOS,rundir,PREFIX,rundir,PREFIX))

   # generate linearize results
   run_process("%s/Linearize -b %s/Scaffold/in/%s.bnk"%(AMOS,rundir,PREFIX))
   run_process("%s/OutputResults -b %s/Scaffold/in/%s.bnk -p %s/Scaffold/out/%s.linearize "%(AMOS,rundir,PREFIX,rundir,PREFIX))
   run_process("%s/OutputScaffolds -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.linearize.scaffolds.final"%(AMOS,rundir,PREFIX,rundir,PREFIX))


if "FindScaffoldORFS" in forcesteps:
    run_process("touch %s/Scaffold/out/%s.linearize.scaffolds.final"%(rundir,PREFIX))
@follows(Scaffold)
@files("%s/Scaffold/out/%s.linearize.scaffolds.final"%(rundir,PREFIX),"%s/FindORFS/out/%s.scaffolds.faa"%(rundir,PREFIX))
def FindScaffoldORFS(input,output):
   if "FindScaffoldORFS" in skipsteps:
      run_process("touch %s/FindScaffoldORFS/out/%s.scaffolds.faa"%(rundir, PREFIX))
      return 0

   run_process("%s/gmhmmp -o %s/FindScaffoldORFS/out/%s.scaffolds.orfs -m %s/config/MetaGeneMark_v1.mod -d -a %s/%s/Scaffold/out/%s.linearize.scaffolds.final"%(GMHMMP,rundir,PREFIX,METAMOS_UTILS,METAMOSDIR,rundir,PREFIX))
   #print"%s/cpp/gmhmmp -o %s/FindORFS/out/%s.scaffolds.orfs -m %s/config/MetaGeneMark_v1.mod -d -a %s/Scafffold/out/%s.linearize.scaffolds.final"%(METAMOS_UTILS,rundir,PREFIX,METAMOS_UTILS,rundir,PREFIX)
   parse_genemarkout("%s/FindScaffoldORFS/out/%s.scaffolds.orfs"%(rundir,PREFIX),1)
   #run_process("unlink ./%s/FindORFS/in/%s.scaffolds.faa"%(rundir,PREFIX))
   #run_process("ln -t ./%s/Annotate/in/ -s ../../FindORFS/out/%s.scaffolds.faa"%(rundir,PREFIX))

if "Propagate" in forcesteps:
    run_process("touch %s/Metaphyler/out/%s.classify.txt"%(rundir,PREFIX))
@follows(FindScaffoldORFS)
@files("%s/Metaphyler/out/%s.classify.txt"%(rundir,PREFIX),"%s/Propagate/out/%s.clusters"%(rundir,PREFIX))
def Propagate(input,output):
   #run propogate java script
   # create s12.annots from Metaphyler output
   run_process("python %s/python/create_mapping.py %s/DB/class_key.tab %s/Metaphyler/out/%s.classify.txt %s/Propagate/in/%s.annots"%(METAMOS_UTILS,METAMOS_UTILS,rundir,PREFIX,rundir,PREFIX))
   # strip headers from file and contig name prefix

   run_process("cat %s/Propagate/in/%s.annots |sed s/contig_//g |grep -v contigID > %s/Propagate/in/%s.clusters"%(rundir,PREFIX,rundir,PREFIX))
   run_process("%s/cpp/FilterEdgesByCluster -b %s/Scaffold/in/%s.bnk -clusters in/s12.clusters -noRemoveEdges > %s/Propagate/out/%s.clusters"%(METAMOS_UTILS,rundir,PREFIX,rundir,PREFIX))

@follows(Propagate)
@files("%s/Propagate/out/%s.clusters"%(rundir,PREFIX),"%s/Classify/out/sorted.txt"%(rundir))
def Classify(input,output):
   run_process("python %s/python/sort_contigs.py %s/Propagate/in/%s.clusters %s/DB/class_key.tab %s/Classify/out %s/Scaffold/in/%s.bnk"%(METAMOS_UTILS, rundir, PREFIX, METAMOS_UTILS,rundir, rundir, PREFIX))

@follows(Classify)
@files("%s/Classify/out/sorted.txt"%(rundir),"%s/Postprocess/%.scf.fa"%(rundir,PREFIX))
def Postprocess():
#create_report.py <metaphyler tab file> <AMOS bnk> <output prefix> <ref_asm>
   #copy files into output for createReport   
   #generate reports
   #linearize
   run_process("cp %s/Metaphyler/out/%s.classify.txt %s/Postprocess/out/. "%(rundir,PREFIX,rundir))
   run_process("cp %s/Scaffold/out/%s.linearize.scaffolds.final %s/Postprocess/out/%s.scf.fa"%(rundir,PREFIX,rundir,PREFIX))
   run_process("ln -t %s/Postprocess/out/ -s ../../Scaffold/in/%s.bnk "%(rundir,PREFIX))
   run_process("python %s/python/create_report.py %s/Postprocess/out/%s.taxprof.pct.txt  %s/Postprocess/out/%s.bnk %s %s/Postprocess/out/%s.scf.fa"%(METAMOS_UTILS,rundir,PREFIX,rundir,PREFIX,PREFIX,rundir,PREFIX))   
   


def parse_genemarkout(orf_file,is_scaff=False):
    coords = open(orf_file,'r')
    coords.readline()
#    outf = open("proba.orfs",'w')
    prevhdr = 0
    prevhdraa = 0
    prevhdrnt = 0

    curcontig = ""
    curseqaa = ""
    curseqnt = ""
    reads = {}
    gene_dict = {}
    fna_dict = {}
    cvg_dict = {}
    for line in coords:
        if ">gene" in line[0:10]:
            if "_nt|" in line:
                #print prevhdraa, prevhdrnt#, curseqaa, curseqnt
                if prevhdraa and curseqaa != "":
                    try:
                        gene_dict[curcontig].append(curseqaa)
                    except KeyError:
                        gene_dict[curcontig] = []
                        gene_dict[curcontig].append(curseqaa)
                    curseqaa = ""

                elif prevhdrnt and curseqnt != "":
                    try:
                        fna_dict[curcontig].append(curseqnt)
                    except KeyError:
                        fna_dict[curcontig] = []
                        fna_dict[curcontig].append(curseqnt)
                    curseqnt = ""

                prevhdrnt = 1
                prevhdraa = 0

            elif "_aa|" in line:

                if prevhdrnt and curseqnt != "":
                    try:
                        fna_dict[curcontig].append(curseqnt)
                    except KeyError:
                        fna_dict[curcontig] = []
                        fna_dict[curcontig].append(curseqnt)
                    curseqnt = ""
                elif prevhdraa and curseqaa != "":
                    try:
                        gene_dict[curcontig].append(curseqaa)
                    except KeyError:
                        gene_dict[curcontig] = []
                        gene_dict[curcontig].append(curseqaa)
                    curseqaa = ""
                prevhdraa = 1
                prevhdrnt = 0

            prevhdr = 1
            lined = line.replace("\n","")
            data = line[1:].split(">",1)[1]
            
            curcontig = data.split(" ")[0]
            if len(data.split(" ")) == 1:
                curcontig = data.split("\t")[0]
            curcontig = curcontig.strip()
            #print curcontig, len(curcontig)
            cvg = data.split(" ")[-1]
            if asm == "soap":
                try:
                    cvg = float(cvg.split("_")[1])
                    cvg_dict[curcontig] = cvg
                except IndexError:
                    #print "Coverage not found, skip?"
                    cvg = 1.0
                    cvg_dict[curcontig] = cvg
            elif asm == "newbler":
                try:
                    #run_process("./sergeTest/Assemble/out/assembly/454ContigGraph.txt")
                    run_process("cat %s/%s/Assemble/out/assembly/454ContigGraph.txt | grep %s | awk \'{print $4}\' > %s/%s/cvg1.out"%(METAMOSDIR,rundir, curcontig, METAMOSDIR,rundir))
                    #print "cat %s/%s/Assemble/out/assembly/454ContigGraph.txt | grep %s | awk \'{print $4}\' > %s/%s/cvg1.out"%(METAMOSDIR,rundir, curcontig, METAMOSDIR,rundir)
                    #run_process("cat %s/%s/Assemble/out/assembly/454ContigGraph.txt | grep %s  > %s/%s/cvg1.out"%(METAMOSDIR,rundir, curcontig, METAMOSDIR,rundir))
                    fin = open("%s/%s/cvg1.out"%(METAMOSDIR,rundir),'r')
                    cvg = fin.readline().replace("\n","")
                    cvg = float(cvg)
                    #print cvg
                    cvg_dict[curcontig] = cvg
                except ValueError:
                    #print "Coverage not found, skip?"
                    cvg = 1.0
                    cvg_dict[curcontig] = cvg                
            prevhdr = 1

        elif len(line) > 2 and prevhdraa == 1 and prevhdr:
            curseqaa += line
        elif len(line) > 2 and prevhdrnt == 1 and prevhdr:
            curseqnt += line
        elif len(line) <= 2 or "Nucleotide" in line: #and prevhdr == 1:
            prevhdr = 0
            #prevhdraa = 0
            #prevhdrnt = 0

        else:
            continue
    if prevhdraa and curseqaa != "":
        try:
          gene_dict[curcontig].append(curseqaa)
        except KeyError:
          gene_dict[curcontig] = []
          gene_dict[curcontig].append(curseqaa)
          curseqaa = ""

    elif prevhdrnt and curseqnt != "":
        try:
          fna_dict[curcontig].append(curseqnt)
        except KeyError:
          fna_dict[curcontig] = []
          fna_dict[curcontig].append(curseqnt)
    if is_scaff:
        outf = open("%s/FindScaffoldORFS/out/%s.faa"%(rundir,PREFIX),'w')
        outf2 = open("%s/FindScaffoldORFS/out/%s.fna"%(rundir,PREFIX),'w')
        cvgf = open("%s/FindScaffoldORFS/out/%s.contig.cvg"%(rundir,PREFIX),'w')
    else:
        outf = open("%s/FindORFS/out/%s.faa"%(rundir,PREFIX),'w')
        outf2 = open("%s/FindORFS/out/%s.fna"%(rundir,PREFIX),'w')
        cvgf = open("%s/FindORFS/out/%s.contig.cvg"%(rundir,PREFIX),'w')
    #print len(gene_dict.keys())
    orfs = {}
    for key in gene_dict.keys():
        genecnt = 1
        if not is_scaff:
            cvgf.write("%s_gene%d\t%s\n"%(key,genecnt,cvg_dict[key])) 
        for gene in gene_dict[key]:
            #min aa length, read depth
            if len(gene) < 100:# or cvg_dict[key] < 5:
                continue
            try:
                #print "contig"+key
                orfs["%s"%(key)] +=1
            except KeyError:
                orfs["%s"%(key)] =1
            outf.write(">%s_gene%d\n%s"%(key,genecnt,gene))

            genecnt +=1
    for key in fna_dict.keys():
        for gene in fna_dict[key]:
            if len(gene) < 300:# or cvg_dict[key] < 5:
                continue
            outf2.write(">%s_gene%d\n%s"%(key,genecnt,gene))
#        print gene_dict[key][0]
    outf.close()
    cvgf.close()
def parse_phmmerout(phmmerout):

    hit_dict = {}
    #phmout = open("%s.phm.tbl"%(prefix),'r')
    phmout = open(phmmerout,'r')
    phmmer_hits = {}
    ctghits = {}
    annot = {}
    for line in phmout:
        line = line.replace("\n","")

        if "gene" in line:
            tts = line.split("[",1)
            if len(tts) < 2:
                 phage_annot = "NA"
            else:
                 line,phage_annot = line.split("[",1)
            phage_annot = phage_annot.replace("]","")
            data = line.split(" ")
            data2 = []
            for item in data:
                if item == "" or item == "-" or item == "\n":
                    continue
                else:
                    data2.append(item)
            try:
                data2[16]
                for git in data2[16:]:
                    phage_annot += " "+git + " "

            except IndexError:
                pass

            data2 = data2[:15]
            #print phage_annot
            #print data2
            #print data2[1].split("_",1)[0]
            try:
                ctghits[data2[1]]
                continue
            except KeyError:
                ctghits[data2[1]] = 1
                pass
            phage_annot = phage_annot.replace(",","")
            try:
                annot[data2[1].split("_",1)[0]] += phage_annot
            except KeyError:
                annot[data2[1].split("_",1)[0]] = phage_annot
            try:
                phmmer_hits[data2[1].split("_",1)[0]] +=1
            except KeyError:
                phmmer_hits[data2[1].split("_",1)[0]] = 1
            try:
                hit_dict[data2[1]]
            except KeyError:
                hit_dict[data2[1]] = [float(data2[2]),int(float(data2[3])),phage_annot]
    #print len(hit_dict.keys())
    #for key in hit_dict.keys():
    #    print hit_dict[key]

if __name__ == "__main__":
    #pid = start_http()
    print "Starting metAMOS pipeline"
    guessPaths()


    
    files = os.listdir(".")
    dlist = []
    pipeline_printout(sys.stdout,[Preprocess,Assemble, FindORFS, FindRepeats, Metaphyler, Scaffold, FindScaffoldORFS, Propagate, Classify, Postprocess], verbose=5)
    pipeline_printout_graph (   'flowchart.svg',
                            'svg',
                            [Postprocess],
                            no_key_legend = True)
    pipeline_run([Preprocess,Assemble, FindORFS, FindRepeats, Metaphyler, Scaffold, FindScaffoldORFS, Propagate, Classify, Postprocess], multiprocess=threads,verbose = 3) 
   
    t2 = time.time()#clock()
    elapsed = float(t2)-float(t1)
    #print elapsed
    print "done! pipeline took %.2f minutes"%(float(elapsed)/float(60.0))
    #os.kill(pid)
