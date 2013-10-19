#!python                                                                                                                                                                                                                                                  

import os, sys, math, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from preprocess import Preprocess
from annotate import Annotate
sys.path.append(INITIAL_UTILS)
from ruffus import *
import pysam
_skipsteps = []
_settings = Settings()
def init(skipsteps):
   global _skipsteps
   _skipsteps = skipsteps

@follows(Annotate)
@posttask(touch_file("%s/Logs/functionalannotation.ok"%(_settings.rundir)))
@files("%s/FindORFS/out/%s.faa"%(_settings.rundir,_settings.PREFIX),["%s/FunctionalAnnotation/out/blast.out"%(_settings.rundir),"%s/FunctionalAnnotation/out/krona.ec.input"%(_settings.rundir)])
def FunctionalAnnotation(input,output):
   if "FunctionalAnnotation" in _skipsteps:
      run_process(_settings, "touch %s/Logs/functionalannotation.skip"%(_settings.rundir), "FunctionalAnnotation")
      run_process(_settings, "touch %s/FunctionalAnnotation/out/blast.out"%(_settings.rundir), "FunctionalAnnotation")
      return 0
   # uniprot_sprot_enz_set
   
   if os.path.exists("%s/uniprot_sprot.fasta"%(_settings.BLASTDB_DIR)):
       
       run_process(_settings,"%s/blastall -p blastp -i %s/FindORFS/out/proba.faa -d %s/uniprot_sprot.fasta -a %s -e 0.001 -m 8 -b 1 > %s/FunctionalAnnotation/out/blast.out"%(_settings.BLAST,_settings.rundir,_settings.BLASTDB_DIR,_settings.threads,_settings.rundir),"FunctionalAnnotation")
   #run_process(_settings,"%s/blastall -p blastx -a %d -m 8 -b 1 -e 1e-2 -i %s -d %s/perl/metaphyler/test/test.ref.protein > %s/Annotate/out/%s.query.blastx"%(_settings.BLAST,_settings.threads,orfFA,_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX))
   #create index of EC codes
   eclines = []
   if os.path.exists("%s/uniprot_sprot_enz_set"%(_settings.BLASTDB_DIR)):
       ecdata = open("%s/uniprot_sprot_enz_set"%(_settings.BLASTDB_DIR),'r')
       eclines = ecdata.readlines()
   ecdict = {}
   for line in eclines:
       line = line.replace("\n","")
       data = line.split(" ")
       #print data
       data2 = []
       for item in data:
           if len(item) <= 1:
               continue
           else:
               data2.append(item)
       seqid = data2[0]
       ecid = data2[-1]
       ecdict[seqid] = ecid
   
   blastout = ""
   blastdict = {}
   #process blast output
   if os.path.exists("%s/FunctionalAnnotation/out/blast.out"%(_settings.rundir)):
       blastout = open("%s/FunctionalAnnotation/out/blast.out"%(_settings.rundir),'r')
   else:
      print "blastall in FunctionalAnnotation failed.."
      sys.exit(1)
   blastdata = blastout.readlines()
   foutput = open("%s/FunctionalAnnotation/out/krona.ec.input"%(_settings.rundir),'w')
   for line in blastdata:
       line = line.replace("\n","")
       items = line.split("\t")
       if len(items) < 10:
           continue
       seqid = items[1]
       pid = float(items[2])
       hitlen = float(items[3])
       evalue = items[10]
       if pid > 80 and hitlen > 50:
           try:
               foutput.write("%s\t%s\t%s\n"%(seqid,ecdict[seqid],evalue))
           except KeyError:
               continue
   foutput.close() 
   #for top hit for each seq, report id, e-vlue and EC value
   #create krona plot
   run_process(_settings,"%s/KronaTools/bin/ktImportEC %s %s/FunctionalAnnotation/out/krona.ec.input"%(_settings.METAMOSDIR,"-l" if _settings.local_krona else "",_settings.rundir), "FunctionalAnnotation")
   
 
