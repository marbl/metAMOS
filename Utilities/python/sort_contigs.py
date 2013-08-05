#!/usr/bin/python

import sys
import os
from utils import *
_settings = Settings()
def sort_contigs(ocf,cf,rcf,ck,orf_fasta,orf_protein,orf_mapf,out_dir,amos_bnk,amos_dir):
    contigs_by_class = { }
    reads_by_class = { }
    orf_to_src = { }
    orf_by_class = { }
    id_class = { }
    id_class["0"] = "UNKNOWN"
    
    orig_class_file = open(ocf)
    class_file = open(cf)
    read_class_file = open(rcf)
    #pass class_key as argument
    class_key = open(ck)
    # pass orf indo
    #orf_fasta = sys.argv[5]
    #orf_protein = sys.argv[6]
    orf_map = open(orf_mapf)
    
    #pass outdir as argument
    #out_dir = sys.argv[8]
    #pass AMOS bank as argument
    #amos_bnk = sys.argv[9]
    # parse in key file
    #amos_dir = sys.argv[10]
    for line in class_key:
        line = line.strip()
        fields = line.split("\t")
        # f1 is id, f2 is class name
    
        if len(fields) != 2:
            print "Error in file format\n"
        else:
            id_class[fields[0]] = fields[1]
    
    # parse orf map to contigs/reads
    for line in orf_map:
        line = line.strip()
        fields = line.split()
    
        if len(fields) != 2:
           print "Error in file format\n"
        elif orf_to_src.has_key(fields[0]):
            orf_to_src[fields[0]].append(fields[1])
        else:
            orf_to_src[fields[0]] = [fields[1]]
    
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
    
    #print "The max class id is %d\n"%(maxClassID)
    id_class[str(maxClassID+1)] = "AMBIGUOUS"
    orig_class_file.close()
    
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
    
        if orf_to_src.has_key(fields[0]):
           if orf_by_class.has_key(fields[1]):
              orf_by_class[fields[1]].extend(orf_to_src[fields[0]])
           else: 
              orf_by_class[fields[1]] = orf_to_src[fields[0]]
    class_file.close()
    
    for line in read_class_file:
        line = line.strip()
        fields = line.split()
        # f1 is read, f2 is class
    
        if len(fields) != 2:
            print "Error in file format\n"
    
        elif reads_by_class.has_key(fields[1]):
            reads_by_class[fields[1]].append(fields[0])
    
        else:
            reads_by_class[fields[1]] = [fields[0]]
    
        if orf_to_src.has_key(fields[0]):
           if orf_by_class.has_key(fields[1]):
              orf_by_class[fields[1]].extend(orf_to_src[fields[0]])
           else:
              orf_by_class[fields[1]] = orf_to_src[fields[0]]
    read_class_file.close()
    
    # output contigs
    for key in contigs_by_class:
        if key not in id_class:
           continue
        class_name = id_class[key]
        path = out_dir + os.sep + class_name + os.sep
        if not os.path.exists(path):
            os.mkdir(path)
    
        f = open(path + class_name + ".eid", 'w')
        f.write("\n".join(contigs_by_class[key]) + "\n")
        f.close()
        run_process(_settings,"%s/bank2fasta -b %s -eid -E '%s%s%s.eid' > '%s%s%s.ctg.fasta'"%(amos_dir,amos_bnk,path,os.sep,class_name,path,os.sep,class_name),"Classify")
    
    # output reads
    for key in reads_by_class:
        if key not in id_class:
           continue
        class_name = id_class[key]
        path = out_dir + os.sep + class_name + os.sep
        if not os.path.exists(path):
            os.mkdir(path)
    
        f = open(path + class_name + ".read.eid", 'w')
        f.write("\n".join(reads_by_class[key]) + "\n")
        f.close()
        run_process(_settings,"%s/dumpreads %s -e -E '%s%s%s.read.eid' > '%s%s%s.read.fasta'"%(amos_dir,amos_bnk,path,os.sep,class_name,path,os.sep,class_name),"Classify")
        run_process(_settings,"%s/dumpreads %s -q -e -E '%s%s%s.read.eid' > '%s%s%s.read.qual'"%(amos_dir,amos_bnk,path,os.sep,class_name,path,os.sep,class_name),"Classify")
        run_process(_settings,"%s/dumpreads %s -f -e -E '%s%s%s.read.eid' > '%s%s%s.read.fastq'"%(amos_dir,amos_bnk,path,os.sep,class_name,path,os.sep,class_name),"Classify")
    
    # finally output the orfs
    for key in orf_by_class:
       if key not in id_class:
          continue
       class_name = id_class[key]
       path = out_dir + os.sep + class_name + os.sep
       if not os.path.exists(path):
          os.mkdir(path)
    
       f = open(path + class_name + ".orf.eid", "w")
       f.write("\n".join(orf_by_class[key]) + "\n")
       f.close()
    
       run_process(_settings,"%s/dumpreads %s -e -E '%s%s%s.orf.eid' > '%s%s%s.orf.fna'"%(amos_dir,orf_fasta,path,os.sep,class_name,path,os.sep,class_name),"Classify")
       run_process(_settings,"%s/dumpreads %s -e -E '%s%s%s.orf.eid' > '%s%s%s.orf.faa'"%(amos_dir,orf_protein,path,os.sep,class_name,path,os.sep,class_name),"Classify")
