#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from findscforfs import FindScaffoldORFS

sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_forcesteps = []
_cls = None
_settings = Settings()

def init(reads, skipsteps, forcesteps, cls):
   global _readlibs
   global _skipsteps
   global _forcesteps
   global _cls
   _readlibs = reads
   _skipsteps = skipsteps
   _forcesteps = forcesteps
   _cls = cls

def parse_metaphyler(giMapping, toTranslate, output):
   giDictionary = {};
   try:
      GIs = open(giMapping, 'r')
   except IOError as e:
      return
   for line in GIs:
      line = line.replace("\n","")
      splitLine = line.split("\t")
      giDictionary[splitLine[0]] = splitLine[1]
   GIs.close();
   try:
      GIs = open(toTranslate, 'r')
   except IOError as e:
      print "Exception opening file %s"%(e)
      return
   outf = open(output, 'w')
   for line in GIs:
      line = line.replace("\n","")
      splitLine = line.split("\t")
      if splitLine[1] in giDictionary:
         outf.write(line.replace(splitLine[1], giDictionary[splitLine[1]]) + "\n")
   GIs.close()
   outf.close()

@follows(FindScaffoldORFS)
@posttask(touch_file("%s/Logs/abundance.ok"%(_settings.rundir)))
@files("%s/Assemble/out/%s.asm.contig"%(_settings.rundir,_settings.PREFIX),"%s/Abundance/out/%s.taxprof.pct.txt"%(_settings.rundir,_settings.PREFIX))
def Abundance(input,output):
   if "Abundance" not in _forcesteps and ("FindORFS" in _skipsteps or "FindScaffoldORFS" in _skipsteps or "Abundance" in _skipsteps):
      # can this be done automatically by ruffus pipeline?
      run_process(_settings, "touch %s/Logs/abundance.skip"%(_settings.rundir), "Abundance")
      run_process(_settings, "touch %s/Abundance/out/%s.taxprof.pct.txt"%(_settings.rundir, _settings.PREFIX), "Abundance")
      run_process(_settings, "touch %s/Abundance/out/%s.classify.txt"%(_settings.rundir, _settings.PREFIX), "Abundance"); 
      return 0

   blastfile = _settings.PREFIX+".blastx"
   blastc = _settings.BLAST + os.sep + "blastall"
   formatc = _settings.BLAST + os.sep + "formatdb"
   run_process(_settings, "ln %s/markers.pfasta %s/Abundance/out/markers.pfasta"%(_settings.DB_DIR, _settings.rundir), "Abundance")
   run_process(_settings, "%s  -p T -i %s/Abundance/out/markers.pfasta"%(formatc,_settings.rundir),"Abundance")
   #update to MetaPhyler 1.25
   run_process(_settings, "%s -p blastp -i %s/FindORFS/out/%s.faa -d %s/Abundance/out/markers.pfasta -m 8 -b 10 -v 10 -a %s -o %s/Abundance/out/%s.blastp"%(blastc, _settings.rundir,_settings.PREFIX,_settings.rundir,_settings.threads,_settings.rundir,_settings.PREFIX),"Abundance")

   run_process(_settings, "perl %s/perl/metaphyler_contigs.pl %s/Abundance/out/%s.blastp %s %s/FindORFS/out/%s.gene.cvg %s/Abundance/out %s"%(_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX,_settings.PREFIX,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.METAMOS_UTILS),"Abundance")
   #run_process(_settings, "./installMetaphyler.pl")
   #run_process(_settings, "$Bin/metaphylerClassify $Bin/markers/markers.$blast.classifier $Bin/markers/markers.taxonomy $prefix.$blast > $prefix.classification")

   # finally add the GI numbers to the results where we can
   parse_metaphyler("%s/markers.toGI.txt"%(_settings.DB_DIR), "%s/Abundance/out/%s.blastp"%(_settings.rundir, _settings.PREFIX), "%s/Abundance/out/%s.gi.blastp"%(_settings.rundir, _settings.PREFIX))
   # generate Krona output
   importMetaPhyler = "%s%sperl%sImportMetaPhyler.pl"%(_settings.METAMOS_UTILS, os.sep, os.sep)
   if not os.path.exists(importMetaPhyler):
      print "Error: Krona importer for MetaPhyler not found in %s. Please check your path and try again.\n"%(importMetaPhyler)
      raise(JobSignalledBreak)
   run_process(_settings, "perl %s %s -c -v -i -p %s/Abundance/out/%s.taxprof.pct.txt:%s"%(importMetaPhyler,"-l" if _settings.local_krona else "",_settings.rundir,_settings.PREFIX, _settings.taxa_level),"Abundance")

   if _cls == 'metaphyler' or _cls == None:
      print "!!No classification selected, using MetaPhyler for taxonomic composition classification\n"
      run_process(_settings, "cp %s/Abundance/out/%s.gi.blastp %s/Postprocess/in/%s.hits"%(_settings.rundir, _settings.PREFIX,_settings.rundir,_settings.PREFIX),"Abundance")

