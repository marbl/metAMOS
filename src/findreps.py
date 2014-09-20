#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from validate import Validate

sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_settings = Settings()

def init(reads, skipsteps):
   global _readlibs
   global _skipsteps

   _readlibs = reads
   _skipsteps = skipsteps


def getContigRepeats(contigFile,outFile):

    contig_repeats = ""
    contig_file = ""
    try:
        contig_repeats = open(outFile,'w')
    except IOError, errmsg:
        print "Error creating output file %s "%(sys.argv[2]), errmsg
        sys.exit(1)

    try:
        contig_file = open(contigFile,'r')
    except IOError, errmsg:
        print "Error opening input file %s "%(sys.argv[1]), errmsg
        sys.exit(1)
    contig_file.close()
    contig_file = open(contigFile,'r')
    concatContig(contigFile)
    run_process(_settings, "%s --minreplen=200 --z=17 --sequence=%s.merged --xmfa=%s.xmfa"%(_settings.REPEATOIRE,contigFile,contigFile),"FindRepeats")
    try:
        repeat_file = open(contigFile+".xmfa",'r')
    except IOError:
        print "Repeatoire failed, check log for details, skipping.."
        return -1
    ctg_dict = {}
    seq_map = {}
    contig_data = contig_file.read()
    num_contigs = contig_data.count(">")
    contig_data = contig_data.split(">")[1:]
    prev_pos = 0
    eid = ""
    iid = 1
    for contig in contig_data:
        hdr,seq = contig.split("\n",1)
        id = hdr.split(" ")[0]
        hdr = hdr.replace(">","").replace("\n","")
        start = prev_pos
        clen = len(seq.replace("\n",""))
        end = prev_pos+clen
        ctg_dict[iid] = [start, end, seq]
        i = 0
        while i < clen:
            seq_map[prev_pos+i] = hdr#iid
            i+=1
        prev_pos = end+1
        iid +=1

    repeat_data = repeat_file.readlines()
    repfam = 1
    reppos = []
    clc = 1
    for line in repeat_data:
        if "=" in line:
          repfam +=1
          ctg_list = []
          for copy in reppos:
             try:
                #print seq_map[int(copy[0])]
                if seq_map[int(copy[0])] == seq_map[int(copy[1])]:
                    ctg_list.append(seq_map[int(copy[0])])
                    #ctg_list.append(seq_map[copy[1]])
             except KeyError:
                 continue
          #print ctg_list

          if len(ctg_list) > 1 and ctg_list.count(ctg_list[0]) != len(ctg_list):
              for item in ctg_list:
                   contig_repeats.write("%d:"%repfam+str(item)+"\n")
          clc +=1
          reppos = []
        if ">" not in line:
            continue
        gg, info = line.split(":",1)
        spos,info = info.split("-",1)
        epos,info = info.split(" ",1)
        orient, info = info.split(" ",1)
#        print spos, epos, orient
        reppos.append([spos,epos])


@follows(Assemble)
@posttask(touch_file("%s/Logs/findrepeats.ok"%(_settings.rundir)))
@files("%s/FindRepeats/in/%s.fna"%(_settings.rundir,_settings.PREFIX),"%s/FindRepeats/out/%s.repeats"%(_settings.rundir,_settings.PREFIX))
def FindRepeats(input,output):
   if "FindRepeats" in _skipsteps or "findrepeats" in _skipsteps or "FindORFS" in _skipsteps or "findorfs" in _skipsteps:
     run_process(_settings, "touch %s/Logs/findrepeats.skip"%(_settings.rundir), "FindRepeats")
     return 0
   #run_process(_settings, "python %s/python/getContigRepeats.py %s/FindRepeats/in/%s.fna %s/FindRepeats/out/%s.repeats"%(_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Findrepeats")
   getContigRepeats("%s/FindRepeats/in/%s.fna"%(_settings.rundir,_settings.PREFIX), "%s/FindRepeats/out/%s.repeats"%(_settings.rundir,_settings.PREFIX))

