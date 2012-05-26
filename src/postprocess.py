#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from classify import Classify

sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_cls = None
_settings = Settings()

def init(reads, skipsteps, cls):
   global _readlibs
   global _skipsteps
   global _cls

   _readlibs = reads
   _skipsteps = skipsteps
   _cls = cls

openbrowser = False
if os.environ.get('DISPLAY') != None:
    openbrowser = True

def start_http(server_class=BaseHTTPServer.HTTPServer,
        handler_class=BaseHTTPServer.BaseHTTPRequestHandler):
    #pid = os.fork()
    server_address = ('localhost', 8111)
    httpd = server_class(server_address, handler_class)
    httpd.serve_forever()
    #return pid

def validate_run(dir):
    run_process(_settings, "./%s/run.sh"%(dir))
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

@follows(Classify)
@posttask(touch_file("%s/Logs/postprocess.ok"%(_settings.rundir)))
@files("%s/Assemble/out/%s.asm.contig"%(_settings.rundir,_settings.PREFIX),"%s/Postprocess/%s.scf.fa"%(_settings.rundir,_settings.PREFIX))
def Postprocess(input,output):
#create_report.py <metaphyler tab file> <AMOS bnk> <output prefix> <ref_asm>
   #copy files into output for createReport   
   #generate reports
   #linearize
   #call KronaReports
   if _cls == 'phmmer':
       if not os.path.exists(_settings.KRONA + os.sep + "ImportPHMMER.pl"):
          print "Error: Krona importer for PHMMER not found in %s. Please check your path and try again.\n"%()
          raise(JobSignalledBreak)
       run_process(_settings, "perl %s/ImportPHMMER.pl -c -v -i %s/Postprocess/in/%s.hits"%(_settings.KRONA,_settings.rundir,_settings.PREFIX),"Postprocess")
   elif _cls == 'blast' or _cls == 'metaphyler' or _cls == None:
       if not os.path.exists(_settings.KRONA + os.sep + "ImportBLAST.pl"):
          print "Error: Krona importer for BLAST not found in %s. Please check your path and try again.\n"%()
          raise(JobSignalledBreak)
       run_process(_settings, "perl %s/ImportBLAST.pl -c -v -i %s/Postprocess/in/%s.hits"%(_settings.KRONA,_settings.rundir,_settings.PREFIX),"Postprocess")
   elif _cls == 'fcp':
       # now ran in Annotate step
       pass
       #if not os.path.exists(_settings.KRONA + os.sep + "ImportFCP.pl"):
       #   print "Error: Krona importer for FCP not found in %s. Please check your path and try again.\n"%()
       #   raise(JobSignalledBreak)
       #run_process(_settings, "perl %s/ImportFCP.pl -c -v -i -p %s/Postprocess/in/%s.epsilon-nb_results.txt"%(_settings.KRONA,_settings.rundir,_settings.PREFIX),"Postprocess")
       run_process(_settings, "unlink %s/Postprocess/out/report.krona.html"%(_settings.rundir), "Postprocess")
       run_process(_settings, "ln -s %s/Annotate/out/report.krona.html %s/Postprocess/out/report.krona.html"%(_settings.rundir, _settings.rundir), "Postprocess")
   elif _cls == 'phymm':
       if not os.path.exists(_settings.KRONA + os.sep + "ImportPHYMM.pl"):
          print "Error: Krona importer for PHYMM not found in %s. Please check your path and try again.\n"%()
          raise(JobSignalledBreak)
       run_process(_settings, "perl %s/ImportPhymmBL.pl -c -v -i %s/Postprocess/in/%s.hits"%(_settings.KRONA,_settings.rundir,_settings.PREFIX),"Postprocess")
   elif _cls == 'phylosift':
       #now ran in Annotate step to generate file for Propogate/Classsify
       pass
       #if not os.path.exists(_settings.KRONA + os.sep + "ImportPhyloSift.pl"):
       #    print "Error: Krona importer for PhyloSift not found in %s. Please check your path and try again.\n"%()
       #    raise(JobSignalledBreak)
       #run_process(_settings, "perl %s/ImportPhyloSift.pl -c -v -i %s/Postprocess/in/%s.hits:%s/Assemble/out/%s.contig.cvg"%(_settings.KRONA,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX), "Postprocess") 
       run_process(_settings, "unlink %s/Postprocess/out/report.krona.html"%(_settings.rundir), "Postprocess")
       run_process(_settings, "ln -s %s/Annotate/out/report.krona.html %s/Postprocess/out/report.krona.html"%(_settings.rundir, _settings.rundir), "Postprocess")

   #command to open webbrowser?
   #try to open Krona output
   if openbrowser:
       if os.path.exists(_settings.rundir + os.sep + "Postprocess" + os.sep + "out" + os.sep + "report.krona.html"):
           webbrowser.open_new("%s%sPostprocess%sout%sreport.krona.html"%(_settings.rundir, os.sep, os.sep, os.sep))
       else:
           print "ERROR: No Krona html file available! skipping"
   #webbrowser.open_new(output.html)
   #webbrowser.open_new_tab(output.html)
   run_process(_settings, "cp -r %s/Preprocess/out/*.fastqc %s/Postprocess/out"%(_settings.rundir, _settings.rundir), "Postprocess")
   run_process(_settings, "cp %s/Abundance/out/%s.classify.txt %s/Postprocess/out/. "%(_settings.rundir,_settings.PREFIX,_settings.rundir),"Postprocess")
   run_process(_settings, "cp %s/Scaffold/out/%s.linearize.scaffolds.final %s/Postprocess/out/%s.scf.fa"%(_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Postprocess")
   run_process(_settings, "ln -t %s/Postprocess/out/ -s %s/Scaffold/in/%s.bnk "%(_settings.rundir,_settings.rundir,_settings.PREFIX),"Postprocess")
#   print "python %s/python/create_report.py %s/Abundance/out/%s.taxprof.pct.txt  %s/Postprocess/out/%s.bnk %s/Postprocess/out/ %s/Postprocess/out/%s.scf.fa %s %s %d"%(_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.rundir,_settings.PREFIX,_settings.METAMOS_UTILS,_settings.AMOS, len(_readlibs))
   run_process(_settings, "python %s/python/create_summary.py %s/Abundance/out/%s.taxprof.pct.txt  %s/Postprocess/out/%s.bnk %s/Postprocess/out/ %s/Postprocess/out/%s.scf.fa %s %s %d"%(_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.rundir,_settings.PREFIX,_settings.METAMOS_UTILS,_settings.AMOS, len(_readlibs)),"Postprocess")
   #webbrowser.open_new_tab(createreport.html)
   if openbrowser:
       if os.path.exists("%s/Postprocess/out/summary.html"%(_settings.rundir)):
           webbrowser.open_new_tab("%s/Postprocess/out/summary.html"%(_settings.rundir))
       else:
           print "ERROR: No Summary html file available! skipping"

