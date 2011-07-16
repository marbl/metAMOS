#!/usr/bin/python

import sys
import os

contigs_by_class = { }
id_class = { }
id_class["0"] = "UNKNOWN"

class_key = open("in/class_key.tab")

# parse in key file
for line in class_key:
    line = line.strip()
    fields = line.split()
    # f1 is id, f2 is class name

    if len(fields) != 2:
        print "Error in file format\n"
    else:
        id_class[fields[0]] = fields[1]



# parse contig class file
for line in sys.stdin:
    line = line.strip()
    fields = line.split()
    # f1 is contig, f2 is class

    if len(fields) != 2:
        print "Error in file format\n"

    elif contigs_by_class.has_key(fields[1]):
        contigs_by_class[fields[1]].append(fields[0])
        
    else:
        contigs_by_class[fields[1]] = [fields[0]]
        contigs_by_class[fields[1]].append(fields[0])


for key in contigs_by_class:
    class_name = id_class[key]
    path = "out/" + class_name + "/"
    if not os.path.exists(path):
        os.mkdir(path)

    f = open(path + class_name + ".iid", 'w')
    f.write("\n".join(contigs_by_class[key]))
    f.close()
    ret = os.system("bank2fasta -b in/s12.bnk -i " + path + class_name + ".iid > " + path + class_name + ".fasta")
    
    





    
