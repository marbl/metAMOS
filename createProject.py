import os, sys, string, time, BaseHTTPServer, getopt, time, datetime
#from datetime import date

from ruffus import *

def usage():
    print "usage: createProject.py -1 file.fastq.1 -2 file.fastq.2 -d projectDir -i 300,500 -f/-q"
    print "options: -s -q, -f, -1, -2, -d, -i"

allsteps = ["Preprocess","Assemble","FindORFS","FindRepeats","Metaphyler","Annotate","Scaffold","FindScaffoldORFS","Propagate","Classify","Postprocess"]

today = datetime.datetime.now()
#todaytime = date.fromtimestamp(time.time())
timestamp = "P_"+today.isoformat().replace("-","_").replace(".","").replace(":","").replace("T","_")
#print timestamp
try:
    opts, args = getopt.getopt(sys.argv[1:], "hfsq1:2:i:d:or", ["help", "fasta=","fastq=","sff=","f1=","f2=","insertlen=","dir=","outtie=","readlen="])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)


id = timestamp
libs = {}
frags = []
cf = ""
format ="fastq"
mated = False
min = ""
max = ""
maxreadlen = 150
f1 = ""
f2 = ""
innie = True
SFFLinkerType = "titanium"

for o, a in opts:
    if o == "-v":
        verbose = True
    elif o == "-o":
        innie = False
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-q"):
        #reads = a
        format = "fastq"
    elif o in ("-r"):
        maxreadlen = int(a)
    elif o in ("-f"):
        #reads = a
        format = "fasta"
    # 454 specfic options
    elif o in ("-s"):
        format = "sff"
    elif o in ("-l"):
        SFFLinkerType = a

    elif o in ("-1"):
        #lib1,min,max = a.split(",")
        #libs[lib1] = [int(min),int(max)]
        f1 = a

    elif o in ("-2"):
        #lib1, min, max  = a.split(",")
        #libs[lib1] = [int(min),int(max)]
        mated = True
        f2 = a
        
    elif o in ("-i"):
        data = a.split(",")
        if len(data) < 2:
            print "Need to provide both min & max!"
            sys.exit(1)
        min,max = a.split(",")
           
    elif o in ("-d"):
        id = a
        #print a

if os.path.exists(id):
    print "Project directory already exists, please specify another"
    print "Alternatively, use runPipeline to run an existing project"
    sys.exit(1)
else:
    os.system("mkdir " + id)
    #create config file
    for dir in allsteps:
        os.system("mkdir ./%s/"%(id)+dir)
        os.system("mkdir ./%s/%s/in"%(id,dir))   
        os.system("mkdir ./%s/%s/out"%(id,dir))   
    #os.system("cp soapconfig.txt ./%s/config.txt"%(id))
    soapf = open("%s/config.txt"%(id),'w')
    soapf.write("max_rd_len=%d\n"%(maxreadlen))
    soapf.close()
#if len(frags) == 0 and len(libs.keys()) == 0:
if f1 == "" and f2 == "":
   print "no reads specified!"
   usage()
   sys.exit(2)




cf = open("./"+id+"/pipeline.ini",'w')
cf.write("#metAMOS pipeline configuration file\n")
#cnt = 1
if format == "fastq":
    cf.write("format:\tfastq\n")

elif format == "sff":
    cf.write("format:\tsff\n")
    cf.write("linker:\t%s\n"%(SFFLinkerType))
else:

    cf.write("format:\tfasta\n")
    if not mated:
        filen = os.path.basename(f1)
        os.system("cp %s.qual %s/Preprocess/in/. "%(f1,id))
    else:
        filen1 = os.path.basename(f1)
        filen2 = os.path.basename(f2)
        os.system("cp %s.qual %s/Preprocess/in/."%(f1,id))
        os.system("cp %.qual %s/Preprocess/in/. "%(f2,id))
if not mated:
    filen = os.path.basename(f1)
    if format == "sff" and min != "":
       cf.write("mated:\tTrue\n")
       min = int(min)
       max = int(max)
       mean = (min+max)/2
       stdev = mean * 0.2
       cf.write("f1:\t%s,%d,%d,%d,%d\n"%(filen,min,max,mean,stdev))
    else:
       cf.write("mated:\tFalse\n")
       cf.write("frg:\t%s\n"%(filen))
    os.system("cp %s %s/Preprocess/in/. "%(f1,id))

    #os.system("ln -t %s -s ./%s/Preprocess/in/%s"%(frg,id,filen))
elif mated:
    cf.write("mated:\tTrue\n")
    filen1 =  os.path.basename(f1)
    filen2 =  os.path.basename(f2)
    min = int(min)
    max = int(max)
    mean = (min+max)/2
    stdev = mean * 0.2
    cf.write("f1:\t%s,%d,%d,%d,%d\n"%(filen1,min,max,mean,stdev))
    cf.write("f2:\t%s,%d,%d,%d,%d\n"%(filen2,min,max,mean,stdev))
    #print "ln -t %s/Preprocess/in/ -s %s"%(id,f1)
    #print "ln -t %s/Preprocess/in/ -s %s"%(id,f2)

    os.system("cp %s %s/Preprocess/in/. "%(f1,id))
    os.system("cp %s  %s/Preprocess/in/. "%(f2,id))
    #open config.txt, edit LIBs
if mated:
    soapf = open("%s/config.txt"%(id),'a')
    soaplib = "[LIB]\n"
    soaplib += "avg_ins="+str((min+max)/2)+"\n"
    if innie:
        soaplib += "reverse_seq=0\n"
    else:
        soaplib += "reverse_seq=1\n"
    soaplib += "asm_flags=3\n"
    soaplib += "rank=1\n"
    if format == "fastq" and mated:
        soaplib += "q1=LIBQ1REPLACE\n"
        soaplib += "q2=LIBQ2REPLACE\n"
    elif format == "fasta" and mated:
        soaplib += "f1=LIBQ1REPLACE\n"
        soaplib += "f2=LIBQ2REPLACE\n"
    elif format == "fastq" and not mated:
        soaplib += "q=LIBQ1REPLACE\n"
    elif format == "fasta" and not mated:
        soaplib += "p=LIBQ1REPLACE\n"
    soapf.write(soaplib)
    soapf.close()
     
    #os.system("ln -t %s/Preprocess/in/%s.1 -s %s.1"%(lib,id,filen))
    #os.system("ln -t %s/Preprocess/in/%s.2 -s %s.2"%(lib,id,filen))


print "Project dir %s successfully created!"%(id)
print "Use runPipeline.py to start Pipeline"
