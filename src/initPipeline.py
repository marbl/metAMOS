#!python

import os, sys, string, time, BaseHTTPServer, getopt, time, datetime
#from datetime import date
#from ruffus import *

#code to discover frozen binary location
application_path = ""
if getattr(sys, 'frozen', False):
    application_path = os.path.dirname(sys.executable)
elif __file__:
    application_path = os.path.dirname(__file__)

#code to update frozen binary runPipeline if newer one exists
try:

   sys._MEIPASS
   #print "Download/install latest version of MetAMOS frozen binary, if available? (note, this might require root/sudo access!)"
   dl = 'y'
   #dl = raw_input("Enter Y/N: ")

   #add autoupdate.txt file in install dir if you would like MetAMOS to auto-update
   try:
       f = open("%s/autoupdate.ok","r")
       f.close()
   except IOError:
       dl = 'n'
   if dl == 'y' or dl == 'Y':
       os.system("wget -P %s -N http://www.cbcb.umd.edu/confcour/temp/runPipeline"%(application_path))
       os.system("wget -P %s -N http://www.cbcb.umd.edu/confcour/temp/initPipeline"%(application_path))
except Exception:
   pass

#check remote filestamp, if newer ask to download & replace

libcounter = 1
class readLib:
    format = ""
    innie = True
    mean = 0
    stdev = 0
    mmin = 0
    mmax = 0
    mated = True
    interleaved = False
    linkerType = "titanium"
    f1 = ""
    f2 = ""
    f12 = ""
    def __init__(self,format,f1,f2="",innie=True,linkertype="titanium",mated=False,interleaved=False):
        global libcounter
        self.id = libcounter
        self.sid = "lib"+str(libcounter)
        libcounter +=1
        self.format = format
        self.mated=mated
        self.interleaved=interleaved
        self.mmin = 0
        self.mmax = 0 
        self.f1 = f1
        self.f2 = f2
        self.linkerType = linkertype
        self.innie = innie

        if interleaved:
           self.f1 = self.f2 = ""
           self.f12 = f1

    def initLib(self):
        self.mean = (self.mmin+self.mmax)/2.0
        self.stdev = ((self.mean-self.mmin)+(self.mmax-self.mean))/2#0.1*self.mean
        

    def validateLib(self):
        pass

    def setMinMax(self,mmin,mmax):
       self.mmin = mmin
       self.mmax = mmax
       self.initLib()

    def setSingleEnd(self,f1):
       self.f1 = f1
       mated = false
       interleaved = false

    def setPairs(self,f2):
       self.f2 = f2
       self.mated = True
       self.interleaved=False

    def __str__(self):
        pass

def usage():
    print "usage: initPipeline -f/-q -1 file.fastq.1 -2 file.fastq.2 -d projectDir -i 300:500 "
    print "options: -s -c -q, -f, -1, -2, -d, -m, -i"
    print "-1: either non-paired file of reads or first file in pair, can be list of multiple separated by a comma"
    print "-2: second paired read file, can be list of multiple separated by a comma"
    print "-c:  fasta file containing contigs"
    print "-d: output project directory (required)"
    print "-f: boolean, reads are in fasta format (default is fastq)"
    print "-h: display help message"
    print "-i: insert size of library, can be list separated by commas for multiple libraries"
    print "-l: SFF linker type"
    print "-m: interleaved file of paired reads"
    print "-o: reads are in outtie orientation (default innie)"
    print "-q: boolean, reads are in fastq format (default is fastq)"
    print "-s/--sff: boolean, reads are in SFF format (default is fastq)"
    
if len(sys.argv) < 2:
    usage()
    sys.exit(1)
allsteps = ["Preprocess","Assemble","MultiAlign","FindORFS","FindRepeats","Abundance","Annotate","FunctionalAnnotation","Scaffold","Propagate","FindScaffoldORFS","Classify","Postprocess"]

today = datetime.datetime.now()
#todaytime = date.fromtimestamp(time.time())
timestamp = "P_"+today.isoformat().replace("-","_").replace(".","").replace(":","").replace("T","_")
#print timestamp
try:
    opts, args = getopt.getopt(sys.argv[1:], "hfsq1:2:m:c:i:d:or:l:", ["help", "fasta","fastq","sff","f1=","f2=","matelib=","asmcontig=","insertlen=","dir=","outtie=","readlen="])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

id = os.path.abspath(timestamp)
libs = {}
frags = []
cf = ""
format ="fastq"
inserts = []
maxreadlen = 150
innie = True
readlibs = []
SFFLinkerType = "titanium"
contigs = ""
lastLib = 0
libs = ""
for o, a in opts:
    if o == "-v":
        verbose = True
    elif o == "-o":
        innie = False
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-c"):
        contigs = a
    elif o in ("-q", "--fastq"):
        #reads = a
        format = "fastq"
    elif o in ("-r"):
        maxreadlen = int(a)
    elif o in ("-f", "--fasta"):
        #reads = a
        format = "fasta"
    # 454 specfic options
    elif o in ("-s", "--sff"):
        format = "sff"
    elif o in ("-l"):
        SFFLinkerType = a

    elif o in ("-1"):
        libs = a.split(",")
        i = 0
        lastLib = len(readlibs)
        while i < len(libs):
           nlib = readLib(format, libs[i], "", innie, SFFLinkerType)
           readlibs.append(nlib)
           i+= 1 
    elif o in ("-2"):
        #lib1, min, max  = a.split(",")
        #libs[lib1] = [int(min),int(max)]
        libs = a.split(",")
        i = 0
        if (lastLib + len(libs) > len(readlibs)):
           print "Error: mismatch in the number of files supplied to the -1 and -2 options\n"
           sys.exit(2) 
        while i < len(libs):
           readlibs[i+lastLib].setPairs(libs[i])
           i+= 1
    elif o in ("-m"):
        #lib1, min, max  = a.split(",")
        #libs[lib1] = [int(min),int(max)]
        libs = a.split(",")
        i = 0
        while i < len(libs):
           nlib = readLib(format, libs[i], "", innie, SFFLinkerType, True, True)
           readlibs.append(nlib)
           i+= 1
        
    elif o in ("-i"):
        ins = a.split(",")
        for insert in ins:
            data = insert.split(":")
            if len(data) < 2:
                print "For library %d, need to provide both min & max insert size using min:max syntax!"%(len(inserts)+1)
                sys.exit(1)
            min,max = data[0],data[1]#insert.split(",")
            inserts.append([min,max])
           
    elif o in ("-d"):
        id = os.path.abspath(a)
        #print a

if len(readlibs) == 0:
   print "no reads specified!"
   usage()
   sys.exit(2)

# now go through all our libs and set their insert sizes
i = 0
j = 0
filesOK = True
errorMessage = ""

while i < len(readlibs):
    mylib = readlibs[i]
    if mylib.interleaved:
        if not os.path.exists(mylib.f12):
            filesOK = False
            errorMessage += "File %s from library %d does not exist\n"%(mylib.f12, i+1)
    elif mylib.mated:
        if not os.path.exists(mylib.f1) or not os.path.exists(mylib.f2):
            filesOK = False
            if not os.path.exists(mylib.f1):
                errorMessage += "File %s from library %d does not exist\n"%(mylib.f1, i+1)
            if not os.path.exists(mylib.f2):
                errorMessage += "File %s from library %d does not exist\n"%(mylib.f2, i+1)
    else:
        if not os.path.exists(mylib.f1):
            filesOK = False  
            errorMessage += "File %s from library %d does not exist\n"%(mylib.f1, i+1)
     
    if readlibs[i].mated:
        if (len(inserts) <= j):
           print "Error: no insert size specified for library %d\n"%(i+1)
           sys.exit(2)
        readlibs[i].setMinMax(int(inserts[j][0]), int(inserts[j][1]))
        j += 1
    i += 1
if filesOK == False:
    print "Error, provided files do not exist: \n%s"%(errorMessage)
    sys.exit(1)

if len(contigs) > 1 and not os.path.exists(contigs):
    print "Error, provided contig file does not exist: ", contigs
    sys.exit(1)
if os.path.exists(id):
    print "Project directory already exists, please specify another"
    print "Alternatively, use runPipeline to run an existing project"
    sys.exit(1)
else:
    os.system("mkdir " + id)
    #create config file
    for dir in allsteps:
        os.system("mkdir %s/"%(id)+dir)
        os.system("mkdir %s/%s/in"%(id,dir))   
        os.system("mkdir %s/%s/out"%(id,dir))   
    #os.system("cp soapconfig.txt %s/config.txt"%(id))
    soapf = open("%s/config.txt"%(id),'w')
    soapf.write("max_rd_len=%d\n"%(maxreadlen))
    soapf.close()



cf = open(id+"/pipeline.ini",'w')
cf.write("#metAMOS pipeline configuration file\n")
#if len(contigs) > 0:
#   #user specified a contig file
cf.write("asmcontigs:\t%s\n"%(contigs))
#cnt = 1
i = 0
while i < len(readlibs):
    mylib = readlibs[i]

    f1 = f2 = ""
    if mylib.interleaved:
       f1 = mylib.f12
    elif mylib.mated:
       f1 = mylib.f1
       f2 = mylib.f2
    else:
       f1 = mylib.f1

    if mylib.format == "fastq":
         cf.write("lib%dformat:\tfastq\n"%(i+1))

    elif mylib.format == "sff":
        cf.write("lib%dformat:\tsff\n"%(i+1))
        cf.write("lib%dlinker:\t%s\n"%(i+1,mylib.linkerType))
    else:
        cf.write("lib%dformat:\tfasta\n"%(i+1))
        if mylib.interleaved or not mylib.mated:
            filen = os.path.basename(f1)
            if os.path.exists("%s.qual"%(f1)):
                os.system("cp %s.qual %s/Preprocess/in/. "%(f1,id))
        else:
            filen1 = os.path.basename(f1)
            filen2 = os.path.basename(f2)
            if os.path.exists("%s.qual"%(f1)):
                os.system("cp %s.qual %s/Preprocess/in/."%(f1,id))
            if os.path.exists("%s.qual"%(f2)):
                os.system("cp %s.qual %s/Preprocess/in/. "%(f2,id))
    if not mylib.mated:
        filen = os.path.basename(f1)
        cf.write("lib%dmated:\tFalse\n"%(i+1))
        cf.write("lib%dinterleaved:\tFalse\n"%(i+1))
        cf.write("lib%dfrg:\t%s\n"%(i+1,filen))
        if os.path.exists("%s"%(f1)):
            os.system("cp %s %s/Preprocess/in/. "%(f1,id))

    #os.system("ln -t %s -s %s/Preprocess/in/%s"%(frg,id,filen))
    elif mylib.mated:
        cf.write("lib%dmated:\tTrue\n"%(i+1))
        if mylib.innie:
            cf.write("lib%dinnie:\tTrue\n"%(i+1))
        else:
            cf.write("lib%dinnie:\tFalse\n"%(i+1))

        if not mylib.interleaved:
            cf.write("lib%dinterleaved:\tFalse\n"%(i+1))
            filen1 =  os.path.basename(f1)
            filen2 =  os.path.basename(f2)
            min = mylib.mmin
            max = mylib.mmax
            mean = mylib.mean
            stdev = mylib.stdev
            cf.write("lib%df1:\t%s,%d,%d,%d,%d\n"%(i+1,filen1,min,max,mean,stdev))
            cf.write("lib%df2:\t%s,%d,%d,%d,%d\n"%(i+1,filen2,min,max,mean,stdev))
            os.system("cp %s %s/Preprocess/in/. "%(f1,id))
            os.system("cp %s  %s/Preprocess/in/. "%(f2,id))
        elif mylib.interleaved:
            cf.write("lib%dinterleaved:\tTrue\n"%(i+1))
            filen1 =  os.path.basename(f1)
            min = mylib.mmin
            max = mylib.mmax
            mean = mylib.mean
            stdev = mylib.stdev
            cf.write("lib%df1:\t%s,%d,%d,%d,%d\n"%(i+1,filen1,min,max,mean,stdev))
            os.system("cp %s %s/Preprocess/in/. "%(f1,id))

    if 1:
        soapf = open("%s/config.txt"%(id),'a')
        soaplib = "[LIB]\n"
        if mylib.mated:
            min = mylib.mmin
            max = mylib.mmax
            soaplib += "avg_ins="+str(mylib.mean)+"\n"
        else:
            soaplib += "avg_ins=0\n"
        if mylib.innie:
            soaplib += "reverse_seq=0\n"
        else:
            soaplib += "reverse_seq=1\n"
        soaplib += "asm_flags=3\n"
        soaplib += "rank=1\n"
        if mylib.format == "fastq" and mylib.mated:
            soaplib += "q1=LIB%dQ1REPLACE\n"%(i+1)
            soaplib += "q2=LIB%dQ2REPLACE\n"%(i+1)
        elif mylib.format == "fasta" and mylib.mated and not mylib.interleaved:
            soaplib += "f1=LIB%dQ1REPLACE\n"%(i+1)
            soaplib += "f2=LIB%dQ2REPLACE\n"%(i+1)
        elif mylib.format == "fasta" and not mylib.mated:
            soaplib += "f=LIB%dQ1REPLACE\n"%(i+1)
        elif mylib.format == "fastq" and not mylib.mated:
            soaplib += "q=LIB%dQ1REPLACE\n"%(i+1)
        elif mylib.format == "fasta" and mylib.mated and mylib.interleaved:
            soaplib += "p=LIB%dQ1REPLACE\n"%(i+1)
        elif mylib.format == "sff" and mylib.mated:
            soaplib += "q1=LIB%dQ1REPLACE\n"%(i+1)
            soaplib += "q2=LIB%dQ2REPLACE\n"%(i+1)
        elif mylib.format == "sff":
            soaplib += "q=LIB%dQ1REPLACE\n"%(i+1)
        soapf.write(soaplib)
        soapf.close()
    i+=1

print "Project dir %s successfully created!"%(id)
print "Use runPipeline.py to start Pipeline"
