import os, sys, string, time, BaseHTTPServer, getopt, time, datetime
#from datetime import date


#from ruffus import *


def usage():
    print "usage: createProject.py -1 file.fastq.1 -2 file.fastq.2 -d projectDir -i 300:500 -f/-q"
    print "options: -s -q, -f, -1, -2, -d, -m, -i"

if len(sys.argv) < 2:
    usage()
    sys.exit(1)
allsteps = ["Preprocess","Assemble","FindORFS","FindRepeats","Metaphyler","Annotate","Scaffold","Propagate","FindScaffoldORFS","Classify","Postprocess"]

today = datetime.datetime.now()
#todaytime = date.fromtimestamp(time.time())
timestamp = "P_"+today.isoformat().replace("-","_").replace(".","").replace(":","").replace("T","_")
#print timestamp
try:
    opts, args = getopt.getopt(sys.argv[1:], "hfsq1:2:m:i:d:or", ["help", "fasta=","fastq=","sff=","f1=","f2=","matelib=","insertlen=","dir=","outtie=","readlen="])
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
interleaved = False
min = ""
max = ""
inserts = []
maxreadlen = 150
f1 = ""
f2 = ""
libs1 = []
libs2 = []
innie = True
innies = []
numlibs = 1
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
        libs1 = a.split(",")
        f1 = libs1[0]#a
        numlibs = len(libs1)
    elif o in ("-2"):
        #lib1, min, max  = a.split(",")
        #libs[lib1] = [int(min),int(max)]
        libs2 = a.split(",")
        mated = True
        interleaved = False
        f2 = libs2[0]#a

    elif o in ("-m"):
        #lib1, min, max  = a.split(",")
        #libs[lib1] = [int(min),int(max)]
        mated = True
        interleaved = True
        libs1 = a.split(",")
        f1 = libs1[0]
        
    elif o in ("-i"):

        ins = a.split(",")
        for insert in ins:
            data = insert.split(":")
            if len(data) < 2:
                print "Need to provide both min & max!"
                sys.exit(1)
            min,max = data[0],data[1]#insert.split(",")
            inserts.append([min,max])
           
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
i = 0
while i < numlibs:
    f1 = libs1[i]
    f2 = ""
    if not interleaved or mated:
        f2 = libs2[i]
    if format == "fastq":
         cf.write("lib%dformat:\tfastq\n"%(i+1))

    elif format == "sff":
        cf.write("lib%dformat:\tsff\n"%(i+1))
        cf.write("lib%dlinker:\t%s\n"%(i+1,SFFLinkerType))
    else:
 
        cf.write("lib%dformat:\tfasta\n"%(i+1))
        if interleaved or not mated:
            filen = os.path.basename(f1)
            os.system("cp %s.qual %s/Preprocess/in/. "%(f1,id))
        else:
            filen1 = os.path.basename(f1)
            filen2 = os.path.basename(f2)
            os.system("cp %s.qual %s/Preprocess/in/."%(f1,id))
            os.system("cp %s.qual %s/Preprocess/in/. "%(f2,id))
    if not mated:
        filen = os.path.basename(f1)
        if format == "sff" and min != "":
            cf.write("lib%dmated:\tTrue\n"%(i+1))
            cf.write("lib%dinterleaved:\tTrue\n"%(i+1))
            min = int(inserts[i][0])
            max = int(inserts[i][1])
            mean = (min+max)/2
            stdev = mean * 0.2
            cf.write("lib%df1:\t%s,%d,%d,%d,%d\n"%(i+1,filen,min,max,mean,stdev))
        else:
            cf.write("lib%dmated:\tFalse\n"%(i+1))
            cf.write("lib%dinterleaved:\tFalse\n"%(i+1))
            cf.write("lib%dfrg:\t%s\n"%(i+1,filen))
        os.system("cp %s %s/Preprocess/in/. "%(f1,id))

    #os.system("ln -t %s -s ./%s/Preprocess/in/%s"%(frg,id,filen))
    elif mated and not interleaved:
        cf.write("lib%dmated:\tTrue\n"%(i+1))
        cf.write("lib%dinterleaved:\tFalse\n"%(i+1))
        filen1 =  os.path.basename(f1)
        filen2 =  os.path.basename(f2)

        min = int(inserts[i][0])
        max = int(inserts[i][1])
        mean = (min+max)/2
        stdev = mean * 0.2
        cf.write("lib%df1:\t%s,%d,%d,%d,%d\n"%(i+1,filen1,min,max,mean,stdev))
        cf.write("lib%df2:\t%s,%d,%d,%d,%d\n"%(i+1,filen2,min,max,mean,stdev))
        os.system("cp %s %s/Preprocess/in/. "%(f1,id))
        os.system("cp %s  %s/Preprocess/in/. "%(f2,id))


    elif mated and interleaved:
        cf.write("lib%dmated:\tTrue\n"%(i+1))
        cf.write("lib%dinterleaved:\tTrue\n"%(i+1))
        filen1 =  os.path.basename(f1)
        #min = int(min)
        #max = int(max)
        min = int(inserts[i][0])
        max = int(inserts[i][1])
        mean = (min+max)/2
        stdev = mean * 0.2
        cf.write("lib%df1:\t%s,%d,%d,%d,%d\n"%(i+1,filen1,min,max,mean,stdev))
        os.system("cp %s %s/Preprocess/in/. "%(f1,id))

    if 1:
        soapf = open("%s/config.txt"%(id),'a')
        soaplib = "[LIB]\n"
        if min != "" and max != "":
            min = int(inserts[i][0])
            max = int(inserts[i][1])
            soaplib += "avg_ins="+str((min+max)/2)+"\n"
        else:
            soaplib += "avg_ins=0\n"
        if innie:
            soaplib += "reverse_seq=0\n"
        else:
            soaplib += "reverse_seq=1\n"
        soaplib += "asm_flags=3\n"
        soaplib += "rank=1\n"
        if format == "fastq" and mated and not interleaved:
            soaplib += "q1=LIB%dQ1REPLACE\n"%(i+1)
            soaplib += "q2=LIB%dQ2REPLACE\n"%(i+1)
        elif format == "fasta" and mated and not interleaved:
            soaplib += "f1=LIB%dQ1REPLACE\n"%(i+1)
            soaplib += "f2=LIB%dQ2REPLACE\n"%(i+1)
        elif format == "fasta" and not mated:
            soaplib += "f=LIB%dQ1REPLACE\n"%(i+1)
        elif format == "fastq" and not mated:
            soaplib += "q=LIB%dQ1REPLACE\n"%(i+1)
        elif format == "fasta" and mated and interleaved:
            soaplib += "p=LiB%dQ1REPLACE\n"%(i+1)
        soapf.write(soaplib)
        soapf.close()
    i+=1

print "Project dir %s successfully created!"%(id)
print "Use runPipeline.py to start Pipeline"
