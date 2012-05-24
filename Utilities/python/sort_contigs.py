#!/usr/bin/python

import sys
import os

contigs_by_class = { }
id_class = { }
id_class["0"] = "UNKNOWN"

orig_class_file = open(sys.argv[1])
class_file = open(sys.argv[2])
#pass class_key as argument
class_key = open(sys.argv[3])
#pass outdir as argument
out_dir = sys.argv[4]
#pass AMOS bank as argument
amos_bnk = sys.argv[5]
# parse in key file
amos_dir = sys.argv[6]
for line in class_key:
    line = line.strip()
    fields = line.split("\t")
    # f1 is id, f2 is class name

    if len(fields) != 2:
        print "Error in file format\n"
    else:
        id_class[fields[0]] = fields[1]

# parse original file to identity ambiguous assignment (which is one more than max previous ID)
maxClassID = 0;
for line in orig_class_file:
    line = line.strip()
    fields = line.split()
    # f1 is contig, f2 is class

    if len(fields) != 2:
        print "Error in file format\n"

    elif maxClassID < int(fields[1]):
       maxClassID = int(fields[1])

print "The max class id is %d\n"%(maxClassID)
id_class[str(maxClassID+1)] = "AMBIGUOUS"

# parse contig class file
for line in class_file:
    line = line.strip()
    fields = line.split()
    # f1 is contig, f2 is class

    if len(fields) != 2:
        print "Error in file format\n"

    elif contigs_by_class.has_key(fields[1]):
        contigs_by_class[fields[1]].append(fields[0])
        
    else:
        contigs_by_class[fields[1]] = [fields[0]]


for key in contigs_by_class:
    class_name = id_class[key]
    path = out_dir + os.sep + class_name + os.sep
    if not os.path.exists(path):
        os.mkdir(path)

    f = open(path + class_name + ".eid", 'w')
    f.write("\n".join(contigs_by_class[key]) + "\n")
    f.close()
    ret = os.system("%s/bank2fasta -b %s -eid -E '%s%s%s.eid' > '%s%s%s.fasta'"%(amos_dir,amos_bnk,path,os.sep,class_name,path,os.sep,class_name))
    
    





    
