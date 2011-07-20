import os, sys, string, time, BaseHTTPServer, getopt, subprocess#
from operator import itemgetter

PREFIX = "proba"

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

sys.path.append(METAMOS_UTILS)
from ruffus import *

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
    # 3. CA
    CA = "%s%sCA%s%s-%s%sbin"%(METAMOSDIR, os.sep, os.sep, OSTYPE, MACHINETYPE.replace("x86_64","amd64"), os.sep)
    if not os.path.exists(CA + os.sep + "gatekeeper"):
       CA = getFromPath("gatekeeper", "Celera Assembler") 
    # 4. Newbler
    NEWBLER = "%s%snewbler"%(METAMOSDIR, os.sep);
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
    opts, args = getopt.getopt(sys.argv[1:], "hd:s:e:o:k:c:a:n:p:tf:vm4", ["help", "projectdir","startat","endat", "minoverlap","kmersize","classifier","assembler","skipsteps","threads","filter","forcesteps","verbose","metaphyler","454"])
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
for o, a in opts:
    if o in ("-v","--verbose"):
        verbose = True
    elif o in ("-h", "--help"):
        usage()
        sys.exit()

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
frgs = []
format = ""
mean = 0
stdev = 0
mmin = 0
mmax = 0
mated = True
linkerType = "titanium"
frg = ""
f1 = ""
f2 = ""

for line in inf:
    line = line.replace("\n","")
    if "#" in line:
        continue
    elif "format:" in line:
        format = line.replace("\n","").split("\t")[-1]
    elif "mated:" in line:
        mated = str2bool(line.replace("\n","").split("\t")[-1])
    elif "linker:" in line:
        linkerType = line.replace("\n","").split("\t")[-1]
    elif "f1:" in line:# or "f2:" in line:
        data = line.split("\t")

        fqlibs[data[0]] = data[1]
        f1 = data[1].split(",")[0]
        inf = data[1].split(",")
        mean = int(inf[3])
        stdev = int(inf[4])
        mmin = int(inf[1])
        mmax = int(inf[2])
        libs.append(f1)

    elif "f2:" in line:# or "f2:" in line:
        data = line.split("\t")

        fqlibs[data[0]] = data[1]
        f2 = data[1].split(",")[0]
        inf = data[1].split(",")
        mean = int(inf[3])
        stdev = int(inf[4])
        mmin = int(inf[1])
        mmax = int(inf[2])
        libs.append(f2)
    elif "frg" in line:

        data = line.split("\t")
        frg = data[1]
        #fqfrags[data[0]] = data[1]
        #frgs.append(data[1])
        libs.append(frg)

def map2contig(min,max,fasta=True):
    bowtie_mapping = 1
    
    readDir = ""
    asmDir = ""
    threads = 0

    contigfile = open("%s/Assemble/out/%s.asm.contig"%(rundir,PREFIX),'r')
    #seqfile = open("%s/Preprocess/out/all.seq.btfilt"%(rundir),'w')
    matefile = open("%s/Preprocess/out/all.seq.mates"%(rundir),'r')
    tigr_file = open("%s/Assemble/out/%s.asm.tigr"%(rundir,PREFIX),'w')
    new_matefile = open("%s/Assemble/out/%s.mappedmates"%(rundir,PREFIX),'w')
    seqdict = {}
    hdr = ""
    cnt = 0
    matedict = {}
    contigdict = {}
    contigdict2 = {}
    readdict = {}
    matedict = {}
    ctgmates = 0
    matectgdict = {}
    mateotdict = {}
    read_lookup = {}
    readcnt = 1
    for line in matefile.xreadlines():
        line = line.replace("\n","")
        mate1, mate2 = line.split("\t")
        mate1 = mate1.replace("@","").replace(">","")
        mate2 = mate2.replace("@","").replace(">","")
        matedict[mate2] = mate1
        read_lookup[readcnt] = mate1
        read_lookup[readcnt+1] = mate2
        readcnt += 2

            
        #matedict[mate2] = mate1
        #print mate1, mate2

    if bowtie_mapping == 1:
        #trim to 25bp
        trim = 0
        if trim:
            f1 = open("%s/Preprocess/out/all.seq"%(rundir))
            f2 = open("%s/Preprocess/out/all.seq.trim"%(rundir),'w')
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
            os.system("%s/bowtie-build %s/Assemble/out/%s.asm.contig %s/Assemble/out/IDX"%(BOWTIE, rundir,PREFIX,rundir))
        if "bowtie" not in skipsteps and fasta:
            if trim:
                os.system("%s/bowtie -p 20 -f -v 1 -M 2 %s/Assemble/out/IDX %s/Preprocess/out/all.seq.trim >& %s/Assemble/out/%s.bout"%(BOWTIE,rundir,rundir,rundir,PREFIX))
            else:
                os.system("%s/bowtie -p 20 -f -l 28 -M 2 %s/Assemble/out/IDX %s/Preprocess/out/all.seq >& %s/Assemble/out/%s.bout"%(BOWTIE,rundir,rundir,rundir,PREFIX))
        elif "bowtie" not in skipsteps:
            if trim:
                os.system("%s/bowtie  -p 20 -v 1 -M 2 %s/Assemble/out/IDX %s/Preprocess/out/all.seq.trim >& %s/Assemble/out/%s.bout"%(BOWTIE,rundir,rundir,rundir,PREFIX))
            else:
                os.system("%s/bowtie  -p 20 -l 28 -M 2 %s/Assemble/out/IDX %s/Preprocess/out/all.seq >& %s/Assemble/out/%s.bout"%(BOWTIE,rundir,rundir,rundir,PREFIX))
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
            #seqfile.write(">%s\n%s\n"%(read,read_seq))
            #seqfile.flush()

    else:
 
        #open soap ReadOnContig
        #some contigs are missing!
        infile = open("%s/Assemble/out/%s.asm.readOnContig"%(rundir,PREFIX),'r')
        #readID, ContigID, startpos, strand
        hdr = infile.readline()
        linecnt = 1
        for line in infile.xreadlines():
            if linecnt % 100000 == 0:
                print linecnt,
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
            #print read, read_lookup[read]
            try:
                contigdict[contig].append([int(spos), int(spos)+epos, strand, read_lookup[read]])
            except KeyError:
                contigdict[contig] = [[int(spos),int(spos)+epos,strand,read_lookup[read]]]
            read_seq = "TEST"
            seqdict[read] = read_seq
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
            print "oops, not in mapping file\n"
            errfile.write("%s\n"%ref)
            continue
        new_ctgfile.write(">%d\n%s"%(ctgcnt,cseq_fmt))#item[1]))
        ctgcnt +=1
        tigr_file.write(cseq_fmt)#item[1])
        contigdict[ref].sort()
        #print contigdict[ref]
        for read in contigdict[ref]:
            
            try:
                if read[0] <= 500 and ctgslen - (int(read[1])) <= 500:
                    matectgdict[read[-1]] = ref
                    mateotdict[read[-1]] = read[2]
            except KeyError:
                pass
            if read[2] == "-":
#                matedict[read[-1]]
                tigr_file.write("#%s(%d) [RC] %d bases, 00000000 checksum. {%d 1} <%d %s>\n"%(read[-1],read[0]-1, len(seqdict[read[-1]]), len(seqdict[read[-1]]), read[0], read[1]))
            else:
                tigr_file.write("#%s(%d) [] %d bases, 00000000 checksum. {1 %d} <%d %s>\n"%(read[-1],read[0]-1, len(seqdict[read[-1]]), len(seqdict[read[-1]]), read[0], read[1]))
            tigr_file.write(seqdict[read[-1]]+"\n")

   
    nonctgmates = 0
    sffcnt = 0
    sfrcnt = 0
    srfcnt = 0
    srrcnt = 0
    nffcnt = 0
    nfrcnt = 0
    nrfcnt = 0
    nrrcnt = 0
    new_matefile.write("library\t110110\t%d\t%d\n"%(mmin,mmax))
    linked_contigs = {}
    for mate in matedict.keys():
        try:
            #matectgdict[mate]
            #matectgdict[matedict[mate]]
            #mateotdict[mate]
            pass
        except KeyError:
            continue
        new_matefile.write("%s\t%s\t110110\n"%(mate,matedict[mate]))
        new_matefile.flush()
        continue
        if matectgdict[mate] == matectgdict[matedict[mate]]:
            ctgmates +=1
            if mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "+":
                sffcnt +=1
            elif mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "-":
                sfrcnt +=1
            elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "-":
                srrcnt +=1
            elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "+":
                srfcnt +=1
        else:
            nonctgmates +=1
            ctgid = [matectgdict[mate],matectgdict[matedict[mate]]]
            ctgid.sort()
            try:
                linked_contigs[ctgid[0]+"*"+ctgid[1]] += 1
            except KeyError:
                linked_contigs[ctgid[0]+"*"+ctgid[1]] = 1

            if mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "+":
                nffcnt +=1
            elif mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "-":
                nfrcnt +=1
            elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "-":
                nrrcnt +=1
            elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "+":
                nrfcnt +=1

    ffcnt = 0
    frcnt = 0
    rfcnt = 0
    rrcnt = 0


    
    for mate in matedict.keys():
        try:
            matectgdict[mate]
            matectgdict[matedict[mate]]
        except KeyError:
            continue
        if mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "+":
            ffcnt +=1
        elif mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "-":
            frcnt +=1
        elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "-":
            rrcnt +=1
        elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "+":
            rfcnt +=1

    #ctgcnt = 1
    #ctgseq = 0
    #ctgsizes = []
    #n50_size = 0
    #n50_mid = 955,000
    print "total contigs: ", ctgcnt
    #print "max contig size: ", max(ctgsizes)
    #print "min contig size: ", min(ctgsizes)
    print "average contig size: ",float(sum(ctgsizes))/float(ctgcnt)
    i = 0
    cp = 0
    ctgsizes.sort()
    ctgsizes = ctgsizes[::-1]
    while i < len(ctgsizes):
        cp += ctgsizes[i]
        if cp >= 955000:
            n50_size = ctgsizes[i]
            break
        i+=1
    print "N50 size: ", n50_size
    print "total pairs of contigs with links:", len(linked_contigs.keys())
    items = linked_contigs.items()
    items.sort(key=itemgetter(1),reverse=True)
    print items[0:10]
    print "---->    ---->", ffcnt
    print "---->    <----", frcnt
    print "<----    ---->", rfcnt
    print "<----    <----", rrcnt
    print "mates located in same contig: ", ctgmates
    print "---->    ---->", sffcnt
    print "---->    <----", sfrcnt
    print "<----    ---->", srfcnt
    print "<----    <----", srrcnt
    print "mates located in diff contig: ", nonctgmates
    print "---->    ---->", nffcnt
    print "---->    <----", nfrcnt
    print "<----    ---->", nrfcnt
    print "<----    <----", nrrcnt

    
    tigr_file.close()
        
    
        
        

def start_http(server_class=BaseHTTPServer.HTTPServer,
        handler_class=BaseHTTPServer.BaseHTTPRequestHandler):
    #pid = os.fork()
    server_address = ('localhost', 8111)
    httpd = server_class(server_address, handler_class)
    httpd.serve_forever()
    #return pid
def validate_run(dir):
    os.system("./%s/run.sh"%(dir))
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
#ASSUME SOAP + 1 LIB for now..
if "Preprocess" in forcesteps:
    os.system("touch %s/Preprocess/in/"%(rundir)+infile)

@files("%s/Preprocess/in/"%(rundir)+infile,["%s/Preprocess/out/all.seq"%(rundir),"%s/Preprocess/out/all.seq.mates"%(rundir)])
def Preprocess(input,output):
   #move input files into Preprocess ./in dir
   #output will either be split fastq files in out, or AMOS bank
   if "Preprocess" in skipsteps or "preprocess" in skipsteps:
       return 0
   if filter == True:
       #for reads+libs
       cnt = 1
       for lib in libs:
           if (format == "fastq"):
               #os.system(" perl %s/perl/prinseq.pl -stats_all -verbose -fastq %s/Preprocess/in/%s"%(METAMOS_UTILS,rundir,lib))
               #print "perl %s/perl/prinseq.pl -fastq %s -min_qual_score 10 -ns_max_p 5 -seq_id f%d_ -out_good %s/Preprocess/out/%s.clean"%(METAMOS_UTILS,lib,cnt,rundir,lib)
               os.system("perl %s/perl/prinseq.pl -fastq %s -seq_id %s -out_good %s/Preprocess/out/%s"%(METAMOS_UTILS,lib,os.path.basename(lib).rsplit(".",1)[0],rundir,lib.rsplit(".",1)[0]))
           else:
               #os.system(" perl %s/perl/prinseq.pl -stats_all -verbose -fasta %s/Preprocess/in/%s"%(METAMOS_UTILS,rundir,lib))
               #os.system("perl %s/perl/prinseq.pl -ns_max_p 5 -fasta %s -min_qual_score 10 -ns_max_p 5 -seq_id %s -out_good %s/Preprocess/out/%s"%(METAMOS_UTILS,lib,os.path.basename(lib).rsplit(".",1)[0],rundir,lib.rsplit(".",1)[0]))
               os.system("perl %s/perl/prinseq.pl -fasta %s -seq_id %s -out_good %s/Preprocess/out/%s"%(METAMOS_UTILS,lib,os.path.basename(lib).rsplit(".",1)[0],rundir,lib.rsplit(".",1)[0]))
       cnt +=1
   else:
       for lib in libs:
           os.system("ln -t ./%s/Preprocess/out/ -s ../../Preprocess/in/%s"%(rundir,lib))
   if format == "sff":
      # generate the fasta files from the sff file
      sffToCACmd = "%s/sffToCA -clear 454 -clear discard-n -trim chop -libraryname sff -output %s/Preprocess/out/%s"%(CA, rundir, PREFIX)
      if (mated == True):
         os.system("%s -linker %s -insertsize %d %d %s"%(sffToCACmd, linkerType, mean, stdev, infile))
      else:
         os.system("%s %s"%(sffToCACmd, infile))
      os.system("%s/gatekeeper -T -F -o %s/Preprocess/out/%s.gkpStore %s/Preprocess/out/%s.frg"%(CA, rundir, PREFIX, rundir, PREFIX))
      os.system("%s/gatekeeper -dumpnewbler %s/Preprocess/out/%s %s/Preprocess/out/%s.gkpStore"%(CA, rundir, PREFIX, rundir, PREFIX))
      os.system("%s/gatekeeper -dumplibraries -tabular %s/Preprocess/out/%s.gkpStore |awk '{if (match($3, \"U\") == 0 && match($1, \"UID\") == 0) print \"library\t\"$1\"\t\"$4-$5*3\"\t\"$4+$5*3}' > %s/Preprocess/out/all.seq.mates"%(CA, rundir, PREFIX, rundir))
      os.system("%s/gatekeeper -dumpfragments -tabular %s/Preprocess/out/%s.gkpStore|awk '{if ($3 != 0 && match($1, \"UID\")==0 && $1 < $3) print $1\"\t\"$3\"\t\"$5}' >> %s/Preprocess/out/all.seq.mates"%(CA, rundir, PREFIX, rundir))
      os.system("unlink %s/Preprocess/out/all.seq"%(rundir))
      os.system("ln -s  ../../Preprocess/out/%s.fna %s/Preprocess/out/all.seq"%(PREFIX, rundir))
      os.system("ln -s ../../Preprocess/out/%s.fna.qual %s/Preprocess/out/all.seq.qual"%(PREFIX, rundir))
      os.system("rm -rf %s/Preproces/out/%s.gkpStore"%(rundir, PREFIX))
      os.system("unlink %s/Preprocess/out/%s.frg"%(rundir, PREFIX))
   elif format == "fasta" and not mated:
      os.system("ln -s  ../../Preprocess/in/%s %s/Preprocess/out/all.seq"%(infile, rundir))
      os.system("ln -s ../../Preprocess/in/%s.qual %s/Preprocess/out/all.seq.qual"%(infile, rundir))
      os.system("touch %s/Preprocess/out/all.seq.mates"%(rundir))
   elif format == "fasta" and mated:
       #FIXME, make me faster!
       os.system("perl %s/perl/shuffleSequences_fasta.pl  %s/Preprocess/out/%s %s/Preprocess/out/%s %s/Preprocess/out/all.seq"%(METAMOS_UTILS,rundir,f1, rundir,f2,rundir))
       os.system("python %s/python/extract_mates_from_fasta.py %s/Preprocess/out/all.seq"%(METAMOS_UTILS,rundir))
       os.system("unlink ./%s/Preprocess/out/%s.mates"%(rundir, lib))
       os.system("ln -t ./%s/Preprocess/out/ -s ../../Preprocess/in/%s.mates"%(rundir,lib))
   elif format == "fastq" and mated:
       #extract mates from fastq
       os.system("perl %s/perl/shuffleSequences_fastq.pl  %s/Preprocess/out/%s %s/Preprocess/out/%s %s/Preprocess/out/all.seq"%(METAMOS_UTILS,rundir,f1, rundir,f2,rundir))
       os.system("python %s/python/extract_mates_from_fastq.py %s/Preprocess/out/all.seq"%(METAMOS_UTILS,rundir))

   #update_soap_config()
   elif asm == "ca":
       #useful for 454, need to get SFF to FRG?
       #/fs/wrenhomes/sergek/wgs-assembler/Linux-amd64/bin/sffToCA
       pass
   elif asm == "amos":
       #call toAmos_new              
       pass
   #params
   #toAmos (-m mates|-f frg)
   #      (-c contig|-s fasta|-q qual|-Q fastq)
   #      -t fastqQualtyType [SCUFL]
   #no mates, fastq, lib info?
   #toAmos_new -Q file.fastq -t fastq_type -b bank.bnk
   

#if asm == "soap"
if "Assemble" in forcesteps:
    os.system("touch %s/Preprocess/out/all.seq"%(rundir))

@files("%s/Preprocess/out/all.seq"%(rundir),["%s/Assemble/out/%s.asm.contig"%(rundir,PREFIX),"%s/Assemble/out/%s.asm.scafSeq"%(rundir,PREFIX)])
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
      #print libs
      #for lib in libs:
      if (format == "fastq" or format == "fasta")  and mated:
          soapd = soapd.replace("LIBQ1REPLACE","%s/Preprocess/out/%s"%(rundir,f1))
          soapd = soapd.replace("LIBQ2REPLACE","%s/Preprocess/out/%s"%(rundir,f2))
      else:
          soapd = soapd.replace("LIBQ1REPLACE","%s"%frg)
      #cnt +=1
      soapw = open("%s/soapconfig.txt"%(rundir),'w')
      soapw.write(soapd)
      soapw.close()
      print "Running SOAPdenovo on input reads..."
      #start stopwatch
      os.system("%s/SOAPdenovo-63mer all  -D -d -R -p %d -K %d -M 3 -L 300 -s %s/soapconfig.txt -o %s/Assemble/out/%s.asm"%(SOAP, threads, kmer, rundir,rundir,PREFIX))#SOAPdenovo config.txt
      #os.system("%s/SOAPdenovo-63mer all -D -d -R -p %d -K %d -M 3 -s %s/soapconfig.txt -o %s/Assemble/out/%s.asm"%(SOAP, threads, kmer, rundir,rundir,PREFIX))#SOAPdenovo config.txt
      #os.system("%s/SOAPdenovo-31mer all  -D 3 -d 2 -R -p %d -M 3 -K %d -s %s/soapconfig.txt -o %s/Assemble/out/%s.asm"%(SOAP, threads, kmer, rundir,rundir,PREFIX))#SOAPdenovo config.txt

      #os.system("ln -s %s/Assemble/out/%s.asm.contig ./%s/FindORFS/in/%s.asm.contig"%(rundir,PREFIX, rundir, PREFIX)) 

      #if OK, convert output to AMOS

   elif asm == "newbler":
      os.system("%s/newAssembly -force %s/Assemble/out"%(NEWBLER, rundir));
      os.system("%s/addRun %s/Assemble/out %s/Preprocess/out/all.seq"%(NEWBLER, rundir, rundir));
      newblerCmd = "%s%srunProject"%(NEWBLER, os.sep)
      # read spec file to input to newbler parameters
      newblerCmd += getProgramParams("newbler.spec", "", "-")
      os.system("%s -cpu %d %s/Assemble/out"%(newblerCmd,threads,rundir));

      # convert to AMOS
      os.system("%s/toAmos -o %s/Assemble/out/%s.mates.afg -m %s/Preprocess/out/all.seq.mates -ace %s/Assemble/out/assembly/454Contigs.ace"%(AMOS,rundir, PREFIX, rundir, rundir));
      # get info on EID/IIDs for contigs
      os.system("cat %s/Assemble/out/%s.mates.afg | grep -A 3 \"{CTG\" |awk '{if (match($1, \"iid\") != 0) {IID = $1} else if (match($1, \"eid\") != 0) {print $1\" \"IID; } }'|sed s/eid://g |sed s/iid://g > %s/Assemble/out/454eidToIID"%(rundir, PREFIX, rundir))
      os.system("java -cp %s convert454GraphToCTL %s/Assemble/out/454eidToIID %s/Assemble/out/assembly/454ContigGraph.txt > %s/Assemble/out/%s.graph.cte"%(METAMOS_JAVA, rundir, rundir, rundir, PREFIX));
      os.system("cat %s/Assemble/out/%s.mates.afg %s/Assemble/out/%s.graph.cte > %s/Assemble/out/%s.afg"%(rundir, PREFIX, rundir, PREFIX, rundir, PREFIX))
    
      # make symlink for subsequent steps
      os.system("rm %s/Assemble/out/%s.asm.contig"%(rundir, PREFIX));
      os.system("ln -s ../../Assemble/out/assembly/454AllContigs.fna %s/Assemble/out/%s.asm.contig"%(rundir, PREFIX))
      if mated == True:
         os.system("ln -s ../../Assemble/out/assembly/454Scaffolds.fna %s/Assemble/out/%s.asm.scafSeq"%(rundir, PREFIX))
      else:
         os.system("ln -s ../../Assemble/out/assembly/454AllContigs.fna %s/Assemble/out/%s.asm.scafSeq"%(rundir, PREFIX))

   elif asm == "amos":
      os.system("Minimus ./Preprocess/out/bank")
   elif asm == "CA":
      #runCA script
      frglist = ""
      for lib in libs:
          if format == "fastq":
              os.system("%s/fastqToCA -insertsize %d %d -libraryname %s -t illumina -innie -fastq %s/Preprocess/in/%s"%(CA, mean[lib],stdev[lib], lib, rundir,PREFIX))
          elif format == "fasta":
              os.system("%s/convert-fasta-to-v2.pl -l %s -mean %d -stddev %d -s %s/Preprocess/in/%s -q %s/Preprocess/in/%s.qual -m matepairids %s/Preprocess/out/%s.mateids %s"%(CA, lib, mean[lib], stdev[lib], rundir, lib, rundir, lib, lib,fff))
          frglist += "%s.frg"%(lib)
      os.system("%s/runCA -p asm -d %s/Assemble/out/ -s %/config/asm.spec %s"%(CA,rundir,METAMOS_UTILSPREFIX,frglist))
      #convert CA to AMOS
      os.system("%s/gatekeeper -dumpfrg -allreads -format2 asm.gkpStore > asm.frg bzip2 asm.frg"%(CA))
      os.system("%s/terminator -g asm.gkpStore -t asm.tigStore/ 2 -o asm bzip2 asm.asm"%(CA))
      os.system("%s/toAmos_new -a asm.asm.bz2 -f asm.frg.bz2 -b asm.bnk -U "%(AMOS))
   #stop here, for now
   #sys.exit(0)
   #check if sucessfully completed   

@follows(Assemble)
@files("%s/Assemble/out/%s.asm.scafSeq"%(rundir,PREFIX),"%s/FindORFS/out/%s.faa"%(rundir,PREFIX))
def FindORFS(input,output):
   if "FindORFS" in skipsteps:
      os.system("touch %s/FindRepeats/in/%s.fna"%(rundir, PREFIX))
      os.system("touch %s/FindORFS/out/%s.faa"%(rundir, PREFIX))
      return 0

   if asm == "soap":
       if not os.path.exists("%s/Assemble/out/%s.asm.scafSeq.contigs"%(rundir,PREFIX)):
           os.system("python %s/python/extract_soap_contigs.py %s/Assemble/out/%s.asm.scafSeq"%(METAMOS_UTILS,rundir,PREFIX))
       os.system("unlink %s/FindORFS/in/%s.asm.scafSeq.contigs"%(rundir,PREFIX))
       os.system("unlink %s/FindORFS/in/%s.asm.contig"%(rundir,PREFIX))
       os.system("ln -t ./%s/FindORFS/in/ -s ../../Assemble/out/%s.asm.scafSeq.contigs"%(rundir,PREFIX))
       os.system("mv %s/FindORFS/in/%s.asm.scafSeq.contigs  %s/FindORFS/in/%s.asm.contig"%(rundir,PREFIX,rundir,PREFIX))
   else:

       os.system("unlink %s/FindORFS/in/%s.asm.contig"%(rundir,PREFIX))
       os.system("ln -t ./%s/FindORFS/in/ -s ../../Assemble/out/%s.asm.contig"%(rundir,PREFIX))


   #os.system("ln -t ./%s/FindORFS/in/ -s ../../Assemble/out/%s.asm.scafSeq.contigs"%(rundir,PREFIX))
   os.system("%s/gmhmmp -o %s/FindORFS/out/%s.orfs -m %s/config/MetaGeneMark_v1.mod -d -a %s/FindORFS/in/%s.asm.contig"%(GMHMMP,rundir,PREFIX,METAMOS_UTILS,rundir,PREFIX))
   parse_genemarkout("%s/FindORFS/out/%s.orfs"%(rundir,PREFIX))
   os.system("unlink ./%s/Annotate/in/%s.faa"%(rundir,PREFIX))
   #os.system("unlink ./%s/Annotate/in/%s.fna"%(rundir,PREFIX))
   os.system("unlink ./%s/FindRepeats/in/%s.fna"%(rundir,PREFIX))
   os.system("ln -t ./%s/Annotate/in/ -s ../../FindORFS/out/%s.faa"%(rundir,PREFIX))
   os.system("ln -t ./%s/FindRepeats/in/ -s ../../FindORFS/out/%s.fna"%(rundir,PREFIX))

@follows(FindORFS)
@files("%s/FindRepeats/in/%s.fna"%(rundir,PREFIX),"%s/FindRepeats/out/%s.repeats"%(rundir,PREFIX))
def FindRepeats(input,output):
   if "FindORFS" in skipsteps or "FindRepeats" in skipsteps:
     return 0

   os.system("python %s/python/getContigRepeats.py  %s/FindRepeats/in/%s.fna %s/FindRepeats/out/%s.repeats"%(METAMOS_UTILS,rundir,PREFIX,rundir,PREFIX))


#@follows(FindRepeats)
@files("%s/Annotate/in/%s.faa"%(rundir,PREFIX),"%s/Annotate/out/%s.hits"%(rundir,PREFIX))
def Annotate(input,output):
   #annotate contigs > 1000bp with FCP
   #lets start by annotating ORFs with phmmer
   if cls == "phmmer":

       os.system("phmmer --cpu %d --F1 0.01 --F2 0.0001 --F3 0.000001 -E 0.01 -o %s/Annotate/out/%s.phm.out --tblout %s/Annotate/out/%s.phm.tbl --notextw %s/Annotate/in/%s.faa %s/DB/allprots.faa"%(threads,rundir,PREFIX,rundir,PREFIX,rundir,PREFIX,METAMOS_UTILS))
       parse_phmmerout("%s/Annotate/out/%s.phm.tbl"%(rundir,PREFIX))
       os.system("mv %s/Annotate/out/%s.phm.tbl  %s/Annotate/out/%s.hits"%(rundir,PREFIX,rundir,PREFIX))
       #os.system("mv %s/Annotate/out/%s.phm.tbl  %s/Annotate/out/%s.annotate"%(rundir,PREFIX,rundir,PREFIX))
   elif cls == "blast":
       os.system("blastall -v 1 -b 1 -a %d -p blastp -m 8 -e 0.00001 -i %s/Annotate/in/%s.faa -d %s/DB/new_all_complete_bacteria.faa -o %s/Annotate/out/%s.blastout"%(threads, rundir,PREFIX,METAMOS_UTILS,rundir,PREFIX))
       os.system("mv %s/Annotate/out/%s.blastout  %s/Annotate/out/%s.hits"%(rundir,PREFIX,rundir,PREFIX))
   elif cls == "fcp":
       print "FCP not yet supported.. stay tuned!"



if "Metaphyler" in forcesteps:
   os.system("touch %s/FindORFS/out/%s.faa"%(rundir,PREFIX))
   #os.system("rm %s/Metaphyler/out/%s.phylum.tab"%(rundir,PREFIX))

@follows(FindORFS)
@files("%s/FindORFS/out/%s.faa"%(rundir,PREFIX),"%s/Metaphyler/out/%s.phylum.tab"%(rundir,PREFIX))
def Metaphyler(input,output):
   if "FindORFS" in skipsteps or "Metaphyler" in skipsteps:
      return 0;

   os.system("unlink ./%s/Metaphyler/in/%s.contig.cvg"%(rundir,PREFIX))
   os.system("unlink ./%s/Metaphyler/in/%s.faa"%(rundir,PREFIX))
   os.system("ln -t ./%s/Metaphyler/in/ -s ../../FindORFS/out/%s.contig.cvg"%(rundir,PREFIX))
   os.system("ln -t ./%s/Metaphyler/in/ -s ../../FindORFS/out/%s.faa"%(rundir,PREFIX))
   blastfile = PREFIX+".blastx"
   os.system("formatdb  -p F -i %s/DB/markers.fna"%(METAMOS_UTILS))
   os.system("perl %s/perl/runblast.pl  %s/Metaphyler/in/%s.faa %s/Metaphyler/out/%s.blastx %s/DB/markers.fna"%(METAMOS_UTILS,rundir,PREFIX, rundir,PREFIX,METAMOS_UTILS))
   print "perl %s/perl/metaphyler_contigs.pl %s/Metaphyler/out/%s.blastx %s/Metaphyler/out/%s %s/Metaphyler/in/%s.contig.cvg"%(METAMOS_UTILS,rundir,PREFIX, rundir, PREFIX, rundir, PREFIX)
   os.system("perl %s/perl/metaphyler_contigs.pl %s/Metaphyler/out/%s.blastx %s/Metaphyler/out/%s %s/Metaphyler/in/%s.contig.cvg"%(METAMOS_UTILS,rundir,PREFIX, rundir, PREFIX, rundir, PREFIX))

   

@follows(Metaphyler)
@files(["%s/Assemble/out/%s.asm.contig"%(rundir,PREFIX),"%s/Preprocess/out/all.seq.mates"%(rundir)],"%s/Scaffold/out/%s.genus"%(rundir,PREFIX))
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
       if 1:
           #os.system("ln -t ./%s/Metaphyler/in/ -s ../../FindORFS/out/%s.faa"%(rundir,PREFIX))
           os.system("rm -rf %s/Scaffold/in/%s.bnk"%(rundir,PREFIX))
           if format == "fasta":
               #if "bowtie" not in skipsteps:
               map2contig(min,max,1)
               os.system("%s/toAmos_new -c %s/Assemble/out/%s.asm.tigr -s %s/Preprocess/out/all.seq -m %s/Preprocess/out//all.seq.mates -b %s/Scaffold/in/%s.bnk "%(AMOS,rundir,PREFIX,rundir, rundir,rundir,PREFIX))
           elif format == "fastq":
               #if "bowtie" not in skipsteps:
               map2contig(min,max,0)
               print "%s/toAmos_new -c %s/Assemble/out/%s.asm.tigr -Q %s/Preprocess/out/all.seq -m %s/Assemble/out/%s.mappedmates -b %s/Scaffold/in/%s.bnk "%(AMOS,rundir,PREFIX,rundir, rundir,PREFIX,rundir,PREFIX)
               os.system("%s/toAmos_new -c %s/Assemble/out/%s.asm.tigr -Q %s/Preprocess/out/all.seq -m %s/Assemble/out/%s.mappedmates -b %s/Scaffold/in/%s.bnk "%(AMOS,rundir,PREFIX,rundir, rundir,PREFIX,rundir,PREFIX))
   elif asm == "newbler":
      os.system("rm -rf %s/Scaffold/in/%s.bnk"%(rundir, PREFIX))
      # build the bank for amos
      os.system("%s/bank-transact -b %s/Scaffold/in/%s.bnk -c -m %s/Assemble/out/%s.afg"%(AMOS,rundir, PREFIX, rundir, PREFIX));

   #calls to Bambus2, goBambus2 script
   # first, parse the parameters
   markRepeatParams = getProgramParams("bambus.spec", "MarkRepeats", "-")
   orientContigParams = getProgramParams("bambus.spec", "OrientContigs", "-")

   os.system("%s/clk -b %s/Scaffold/in/%s.bnk"%(AMOS,rundir,PREFIX))
   os.system("%s/Bundler -b %s/Scaffold/in/%s.bnk"%(AMOS,rundir,PREFIX))
   os.system("%s/MarkRepeats %s -b %s/Scaffold/in/%s.bnk > %s/Scaffold/in/%s.reps"%(AMOS,markRepeatParams,rundir,PREFIX,rundir,PREFIX))
   os.system("%s/OrientContigs %s -b %s/Scaffold/in/%s.bnk -repeats %s/Scaffold/in/%s.reps "%(AMOS,orientContigParams,rundir,PREFIX, rundir, PREFIX))

   # output results
   os.system("%s/bank2fasta -d -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.contigs"%(AMOS,rundir,PREFIX,rundir,PREFIX))
   os.system("%s/OutputMotifs -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.motifs"%(AMOS,rundir,PREFIX,rundir,PREFIX))
   os.system("%s/OutputResults -b %s/Scaffold/in/%s.bnk -p %s/Scaffold/out/%s "%(AMOS,rundir,PREFIX,rundir,PREFIX))
   os.system("%s/OutputScaffolds -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.scaffolds.final"%(AMOS,rundir,PREFIX,rundir,PREFIX))

   # generate linearize results
   os.system("%s/Linearize -b %s/Scaffold/in/%s.bnk"%(AMOS,rundir,PREFIX))
   os.system("%s/OutputResults -b %s/Scaffold/in/%s.bnk -p %s/Scaffold/out/%s.linearize "%(AMOS,rundir,PREFIX,rundir,PREFIX))
   os.system("%s/OutputScaffolds -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.linearize.scaffolds.final"%(AMOS,rundir,PREFIX,rundir,PREFIX))

@follows(Scaffold)
@files("%s/Scaffold/out/%s.linearize.scaffolds.final"%(rundir,PREFIX),"%s/FindORFS/out/%s.scaffolds.faa"%(rundir,PREFIX))
def FindScaffoldORFS(input,output):
   if "FindScaffoldORFS" in skipsteps:
      os.system("touch %s/FindScaffoldORFS/out/%s.scaffolds.faa"%(rundir, PREFIX))
      return 0

   os.system("%s/cpp/gmhmmp -o %s/FindORFS/out/%s.scaffolds.orfs -m %s/config/MetaGeneMark_v1.mod -d -a %s/FindScaffoldORFS/out/%s.linearize.scaffolds.final"%(METAMOS_UTILS,rundir,PREFIX,METAMOS_UTILS,rundir,PREFIX))
   parse_genemarkout("%s/FindORFS/out/%s.scaffolds.orfs"%(rundir,PREFIX))
   os.system("unlink ./%s/Annotate/in/%s.scaffolds.faa"%(rundir,PREFIX))
   os.system("ln -t ./%s/Annotate/in/ -s ../../FindORFS/out/%s.scaffolds.faa"%(rundir,PREFIX))

@follows(Scaffold)
@files("%s/Scaffold/in/%s.bnk"%(rundir,PREFIX),"%s/Scaffold/out/%s.genus"%(rundir,PREFIX))
def Propagate(input,output):
   #run propogate java script
   # create s12.annots from Metaphyler output
   os.system("python %s/python/create_mapping.py %s/DB/class_key.tab %s/Metaphyler/out/%s.blastx %s/Propagate/in/%s.annots"%(METAMOS_UTILS,METAMOS_UTILS,rundir,PREFIX,rundir,PREFIX))
   # strip headers from file and contig name prefix

   os.system("cat %s/Propagate/in/%s.annots |sed s/contig_//g |grep -v contigID > %s/Propagate/in/%s.clusters"%(rundir,PREFIX,rundir,PREFIX))
   os.system("%s/FilterEdgesByCluster -noRemoveEdges -b %s/Scaffold/in/%s.bnk -clusters in/s12.clusters -noRemoveEdges > %s/Propagate/out/%s.clusters"%(AMOS,rundir,PREFIX,rundir,PREFIX))

@follows(Propagate)
@files("%s/Propagate/out/%s.clusters"%(rundir,PREFIX),"%s/Classify/out/sorted.txt"%(rundir))
def Classify(input,output):
   #run Dan's classify script
   os.system("python %s/python/sort_contigs.py %s/Propagate/in/%s.clusters %s/DB/class_key.tab %s/Classify/out %s/Scaffold/in/%s.bnk"%(rundir, PREFIX, METAMOS_UTILS, METAMOS_UTILS,rundir, rundir, PREFIX))

@follows(Scaffold)
def Postprocess():
#create_report.py <metaphyler tab file> <AMOS bnk> <output prefix> <ref_asm>
   #copy files into output for createReport   
   #generate reports
   #linearize
   os.system("cp %s/Metaphyler/out/%s.phylum.tab %s/Postprocess/out/. "%(rundir,PREFIX,rundir))
   os.system("cp %s/Scaffold/out/%s.linearize.scaffolds.final %s/Postproces/out/%s.scf.fa"%(rundir,PREFIX,rundir,PREFIX))
   os.system("ln -t %s/Postprocess/out/ -s ../../Scaffold/in/%s.bnk "%(rundir,PREFIX))
   os.system("python %s/python/create_report.py %s/Postprocess/out/%s.phylum.tab  %s/Postprocess/out/%s.bnk %s %s/Postprocess/out/%s.scf.fa"%(METAMOS_UTILS,rundir,PREFIX,rundir,PREFIX,PREFIX,rundir,PREFIX))   
   


def parse_genemarkout(orf_file):
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
        if ">gene" in line:
            if "_nt" in line:
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

            elif "_aa" in line:

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
            cvg = data.split(" ")[-1]
            cvg = float(cvg.split("_")[1])
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
    outf = open("%s/FindORFS/out/%s.faa"%(rundir,PREFIX),'w')
    outf2 = open("%s/FindORFS/out/%s.fna"%(rundir,PREFIX),'w')
    cvgf = open("%s/FindORFS/out/%s.contig.cvg"%(rundir,PREFIX),'w')
    print len(gene_dict.keys())
    orfs = {}
    for key in gene_dict.keys():
        genecnt = 1
        cvgf.write("%s_gene%d\t%s\n"%(key,genecnt,cvg_dict[key])) 
        for gene in gene_dict[key]:
            #min aa length, read depth
            if len(gene) < 100 or cvg_dict[key] < 5:
                continue
            try:
                #print "contig"+key
                orfs["%s"%(key)] +=1
            except KeyError:
                orfs["%s"%(key)] =1
            outf.write(">%s_gene%d\n%s"%(key,genecnt,gene))

            genecnt +=1
        for gene in fna_dict[key]:
            if len(gene) < 300 or cvg_dict[key] < 5:
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
    print len(hit_dict.keys())
    #for key in hit_dict.keys():
    #    print hit_dict[key]

if __name__ == "__main__":
    #pid = start_http()
    print "Starting metAMOS pipeline"
    guessPaths()

    t1 = time.clock()
    
    files = os.listdir(".")
    dlist = []
    pipeline_printout(sys.stdout,[Preprocess,Assemble, FindORFS, FindRepeats, Metaphyler, Scaffold, Propagate, FindScaffoldORFS, Classify, Postprocess], verbose=5)
    pipeline_printout_graph (   'flowchart.svg',
                            'svg',
                            [Postprocess],
                            no_key_legend = True)
    pipeline_run([Preprocess,Assemble, FindORFS, FindRepeats, Metaphyler, Scaffold, Propagate, FindScaffoldORFS, Classify, Postprocess], verbose = 5) 
   
    t2 = time.clock()
    elapsed = t2-t1
    print "done! pipeline took %.2f minutes"%(float(elapsed)/float(60))
    #os.kill(pid)
