#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter
from math import ceil

from utils import *
from findreps import FindRepeats
#from findorfs import FindORFS

sys.path.append(INITIAL_UTILS)
from ruffus import *
from splitfasta import *

from multiprocessing import *

import generic

_MIN_SEQ_LENGTH = 10000000
_USE_GRID = 0

_readlibs = []
_skipsteps = []
_settings = Settings()
_cls = None
_noblast = False

_FCP_MODELS = "%s"%(_settings.METAMOS_UTILS)
if _settings.BINARY_DIST:
   _FCP_MODELS = "%s"%(_settings.DB_DIR)
def init(reads, skipsteps, cls, noblast):
   global _readlibs
   global _skipsteps
   global _cls
   global _noblast
   _readlibs = reads
   _skipsteps = skipsteps
   _cls = cls
   _noblast = noblast

def parse_phmmerout(phmmerout):

    hit_dict = {}
    #phmout = open("%s.phm.tbl"%(prefix),'r')
    phmout = open(phmmerout,'r')
    phmmer_hits = {}
    ctghits = {}
    annot = {}
    for line in phmout:
        line = line.replace("\n","")

        if "gene" in line:
            tts = line.split("[",1)
            if len(tts) < 2:
                 phage_annot = "NA"
            else:
                 line,phage_annot = line.split("[",1)
            phage_annot = phage_annot.replace("]","")
            data = line.split(" ")
            data2 = []
            for item in data:
                if item == "" or item == "-" or item == "\n":
                    continue
                else:
                    data2.append(item)
            try:
                data2[16]
                for git in data2[16:]:
                    phage_annot += " "+git + " "

            except IndexError:
                pass

            data2 = data2[:15]
            #print phage_annot
            #print data2
            #print data2[1].split("_",1)[0]
            try:
                ctghits[data2[1]]
                continue
            except KeyError:
                ctghits[data2[1]] = 1
                pass
            phage_annot = phage_annot.replace(",","")
            try:
                annot[data2[1].split("_",1)[0]] += phage_annot
            except KeyError:
                annot[data2[1].split("_",1)[0]] = phage_annot
            try:
                phmmer_hits[data2[1].split("_",1)[0]] +=1
            except KeyError:
                phmmer_hits[data2[1].split("_",1)[0]] = 1
            try:
                hit_dict[data2[1]]
            except KeyError:
                hit_dict[data2[1]] = [float(data2[2]),int(float(data2[3])),phage_annot]
    #print len(hit_dict.keys())
    #for key in hit_dict.keys():
    #    print hit_dict[key]


def annotateSeq(cls, contigs, orfAA, orfFA, output):
   #annotate contigs > 1000bp with FCP
   #lets start by annotating ORFs with phmmer

   if cls == "phmmer":
       if not os.path.exists(_settings.PHMMER + os.sep + "phmmer"):
          print "Error: PHMMER not found in %s. Please check your path and try again.\n"%(_settings.PHMMER)
          raise(JobSignalledBreak)

       if not os.path.exists("%s/allprots.faa"%(_settings.BLASTDB_DIR)):
          print "Error: You indicated you would like to run phmmer but DB allprots.faa not found in %s. Please check your path and try again.\n"%(_settings.BLASTDB_DIR)
          raise(JobSignalledBreak)

       run_process(_settings, "%s/phmmer --cpu %d -E 0.0000000000000001 -o %s/Annotate/out/%s.phm.out --tblout %s/Annotate/out/%s.phm.tbl --notextw %s %s/allprots.faa"%(_settings.PHMMER, _settings.threads,_settings.rundir,output,_settings.rundir,output,orfAA,_settings.BLASTDB_DIR),"Annotate")
       parse_phmmerout("%s/Annotate/out/%s.phm.tbl"%(_settings.rundir,output))
       run_process(_settings, "mv %s/Annotate/out/%s.phm.tbl  %s/Annotate/out/%s.intermediate.hits"%(_settings.rundir,output,_settings.rundir,output),"Annotate")

   elif cls == "metaphyler":
       if not os.path.exists(_settings.BLAST + os.sep + "blastall"):
           print "Error: BLAST not found in %s. Please check your path and try again.\n"%(_settings.BLAST)
           raise(JobSignalledBreak)

       #run_process(_settings, "perl %s/perl/installMetaphyler.pl"%(_settings.METAMOS_UTILS)
       run_process(_settings,"%s/blastall -p blastx -a %d -m 8 -b 1 -e 1e-2 -i %s -d %s/perl/metaphyler/test/test.ref.protein > %s/Annotate/out/%s.query.blastx"%(_settings.BLAST,_settings.threads,orfFA,_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX))
       run_process(_settings, "%s/metaphylerClassify %s/perl/metaphyler/markers/markers.blastx.classifier %s/perl/metaphyler/markers/markers.taxonomy %s/Annotate/out/%s.query.blastx > %s/Annotate/out/%s.classification"%(_settings.METAPHYLER,_settings.METAMOS_UTILS,_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX,_settings.rundir,output) )
       #Krona import Metaphyler script!
       run_process(_settings, "mv %s/Annotate/out/%s.classification  %s/Annotate/out/%s.intermediate.hits"%(_settings.rundir,output,_settings.rundir,output),"Annotate")

   elif cls == "blast":
       if not os.path.exists(_settings.BLAST + os.sep + "blastall"):
          print "Error: BLAST not found in %s. Please check your path and try again.\n"%(_settings.BLAST)
          raise(JobSignalledBreak)

       run_process(_settings, "%s/blastall -v 1 -b 1 -a %d -p blastp -m 8 -e 0.00001 -i %s -d %s/refseq_protein -o %s/Annotate/out/%s.blastout"%(_settings.BLAST, _settings.threads,orfAA,_settings.DB_DIR,_settings.rundir,output),"Annotate")
       run_process(_settings, "mv %s/Annotate/out/%s.blastout  %s/Annotate/out/%s.intermediate.hits"%(_settings.rundir,output,_settings.rundir,output),"Annotate")
   elif cls == "phylosift":
       if _settings.PHYLOSIFT == "" or not os.path.exists(_settings.PHYLOSIFT + os.sep + "bin" + os.sep + "phylosift"):
          print "Error: PhyloSift not found in %s. Please check your path and try again.\n"%(_settings.PHYLOSIFT)
          raise(JobSignalledBreak)

       phylosiftCmd =  "%s/bin/phylosift all --threads=%d"%(_settings.PHYLOSIFT, _settings.threads)
       phylosiftCmd += " %s"%(getProgramParams("phylosift.spec", "", "--"))
       # run on contigs for now
       #for lib in readlibs:
       #   if lib.mated:
       #       if not lib.innie or lib.interleaved:
       #          print "Warning: PhyloSift only supports innie non-interleaved libraries now, skipping library %d"%(lib.id)
       #       else:
       #          run_process(_settings, "%s -paired %s/Preprocess/in/%s %s/Preprocess/in/%s"%(phylosiftCmd,_settings.rundir,lib.f1.fname,_settings.rundir,lib.f2.fname), "Annotate")
       #   else:
       #      run_process(_settings, "%s %s/Preprocess/out/lib%d.seq"%(phylosiftCmd,_settings.rundir,lib.id), "Annotate")

       run_process(_settings, "rm -rf %s/Annotate/out/PS_temp"%(_settings.rundir), "Annotate")
       run_process(_settings, "%s %s --coverage=%s/Assemble/out/%s.contig.cnt "%(phylosiftCmd, contigs, _settings.rundir,_settings.PREFIX), "Annotate")

       # save the results
       run_process(_settings, "unlink %s/Annotate/out/%s.intermediate.hits"%(_settings.rundir, output), "Annotate")
       run_process(_settings, "ln %s/Annotate/out/PS_temp/%s/sequence_taxa_summary.txt %s/Annotate/out/%s.intermediate.hits"%(_settings.rundir, os.path.basename(contigs), _settings.rundir, output), "Annotate") 
       
   elif cls == "fcp":
       run_process(_settings, "ln -s %s/models"%(_FCP_MODELS), "Annotate")
       run_process(_settings, "ln -s %s/models/taxonomy.txt"%(_FCP_MODELS), "Annotate")

       # normal nb classify
       run_process(_settings, "%s/nb-classify -q %s -m %s/models/models.txt -r %s/Annotate/out/%s.nb_results.txt -e %s"%(_settings.FCP,contigs,_FCP_MODELS,_settings.rundir,output,output), "Annotate")

       # for blast options
       if not _noblast and os.path.exists("%s/blast_data/BacteriaAndArchaeaGenomesDB.nin"%(_settings.BLASTDB_DIR)):
          run_process(_settings, "ln -s %s/blast_data blast_data"%(_settings.BLASTDB_DIR), "Annotate")
          #run_process(_settings, "python %s/python/BLASTN.py %s/blastn %s %s/Annotate/out/%s.bl_results.txt %d"%(_settings.METAMOS_UTILS, _settings.BLAST, contigs, _settings.rundir, output, _settings.threads), "Annotate") 
          run_process(_settings, "%s/blastn -query %s -out  %s/Annotate/out/%s.bl_results.txt -db %s/blast_data/BacteriaAndArchaeaGenomesDB -evalue 10 -outfmt 7 -task blastn -num_threads %d"%(_settings.BLAST,contigs,_settings.rundir,output,_settings.BLASTDB_DIR,_settings.threads),"Annotate")
          #os.system(blastnEXE + ' -query ' + queryFile + ' -db ./blast_data/BacteriaAndArchaeaGenomesDB -evalue 10 -outfmt 7 -task blastn -num_threads %s -out '%(threads) + resultsFile)
          #combine the results
          if _settings.BINARY_DIST:
              run_process(_settings, "%s/NB-BL %s/Annotate/out/%s.nb_results.txt %s/Annotate/out/%s.bl_results.txt %s/Annotate/out/%s.intermediate.epsilon-nb_results.txt"%(_settings.FCP, _settings.rundir, output, _settings.rundir, output, _settings.rundir, output), "Annotate")
          else:
              run_process(_settings, "python %s/python/NB-BL.py %s/Annotate/out/%s.nb_results.txt %s/Annotate/out/%s.bl_results.txt %s/Annotate/out/%s.intermediate.epsilon-nb_results.txt"%(_settings.METAMOS_UTILS, _settings.rundir, output, _settings.rundir, output, _settings.rundir, output), "Annotate")

       else:
          if _settings.BINARY_DIST:
              run_process(_settings, "%s/Epsilon-NB %s/Annotate/out/%s.nb_results.txt 1E5 %s/Annotate/out/%s.intermediate.epsilon-nb_results.txt"%(_settings.FCP,_settings.rundir,output,_settings.rundir,output),"Annotate")
          else:
              run_process(_settings, "python %s/python/Epsilon-NB.py %s/Annotate/out/%s.nb_results.txt 1E5 %s/Annotate/out/%s.intermediate.epsilon-nb_results.txt"%(_settings.METAMOS_UTILS,_settings.rundir,output,_settings.rundir,output),"Annotate")

       #need python TaxonomicSummary.py test.fasta nb_topModels.txt nb_taxonomicSummary.txt
       #run_process(_settings, "python %s/python/TaxonomicSummary.py %s/Annotate/in/%s.fna %s/Annotate/out/%s.nb_results.txt %s/Annotate/out/%s.epsilon-nb_results.txt"%(_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Annotate")

       run_process(_settings, "unlink %s/Annotate/out/%s.intermediate.hits"%(_settings.rundir, output), "Annotate")
       run_process(_settings, "unlink %s/Annotate/out/%s.nb_results.txt"%(_settings.rundir, output), "Annotate")
       run_process(_settings, "unlink %s/Annotate/out/%s.bl_results.txt"%(_settings.rundir, output), "Annotate")
       run_process(_settings, "ln %s/Annotate/out/%s.intermediate.epsilon-nb_results.txt %s/Annotate/out/%s.intermediate.hits"%(_settings.rundir, output, _settings.rundir, output), "Annotate")

   elif cls == "phymm":
       if not os.path.exists("%s"%(_settings.PHYMM)):
           print "Error: Phymm not found in %s but selected as classifier. Please check your path and try again.\n"%(_settings.PHYMM)
           raise(JobSignalledBreak)
       # link to the files phymm expects locally
       run_process(_settings, "ln -s %s/.blastData"%(_settings.PHYMM), "Annotate")
       run_process(_settings, "ln -s %s/.genomeData"%(_settings.PHYMM), "Annotate")
       run_process(_settings, "ln -s %s/.scripts"%(_settings.PHYMM), "Annotate")
       run_process(_settings, "ln -s %s/.taxonomyData"%(_settings.PHYMM), "Annotate")
       run_process(_settings, "mkdir -f %sAnnotate/out/.logs"%(_settings.rundir), "Annotate")

       run_process(_settings, "perl %s/scoreReads.pl %s > %s/Annotate/out/%s.phymm.err"%(_settings.PHYMM, contigs,_settings.rundir,output),"Annotate")
       run_process(_settings, "cat results.03.phymmBL_%s.txt | grep -v \"QUERY_ID\" > %s/Annotate/out/%s.intermediate.phymm.out"%(contigs.replace(os.sep, "_").replace(".", "_"), _settings.rundir, output), "Annotate")
       run_process(_settings, "ln %s/Annotate/out/%s.intermediate.phymm.out %s/Annotate/out/%s.intermediate.hits"%(_settings.rundir, output, _settings.rundir, output), "Annotate")
       run_process(_settings, "rm %s/Annotate/out/*_%s.txt "%(_settings.rundir, contigs.replace(os.sep, "_").replace(".","_")),"Annotate")
       
   elif cls == None:
       print "No method specified, skipping"

def parallelWrapper(params):
   try:
      jobID = params["jobID"]
      result = {}
      result["jobID"] = jobID
      result["status"] = 0
 
      annotateSeq(params["cls"], params["contigs"], params["orfAA"], params["orfFA"], params["out"])

      result["status"] = 1
      return result
   except KeyboardInterrupt:
      result["status"] = 0
      print "Keyboard error in thread %d, quitting\n"%(jobID)
      return result
   except Exception:
      result["status"] = 0
      print "Other error in thread %d, quitting\n"%(jobID)
      return result

@follows(FindRepeats)
@posttask(touch_file("%s/Logs/annotate.ok"%(_settings.rundir)))
@files("%s/Annotate/in/%s.faa"%(_settings.rundir,_settings.PREFIX),"%s/Annotate/out/%s.hits"%(_settings.rundir,_settings.PREFIX))
def Annotate(input,output):
   if "Annotate" in _skipsteps or _cls == None:
      run_process(_settings, "touch %s/Logs/annotate.skip"%(_settings.rundir), "Annotate")
      run_process(_settings, "touch %s/Annotate/out/%s.hits"%(_settings.rundir, _settings.PREFIX), "Annotate")
      return 0

   listOfFiles = "%s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX)

   # clean up any existing files
   run_process(_settings, "touch %s/Annotate/out/%s.annots"%(_settings.rundir, _settings.PREFIX), "Annotate")
   run_process(_settings, "unlink %s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX), "Annotate")
   run_process(_settings, "ln -s %s/Assemble/out/%s.asm.contig %s/Annotate/in/"%(_settings.rundir, _settings.PREFIX, _settings.rundir), "Annotate")
   run_process(_settings, "unlink %s/Annotate/out/%s.hits"%(_settings.rundir, _settings.PREFIX), "Annotate")
   run_process(_settings, "rm -f %s/Annotate/out/*.hits"%(_settings.rundir), "Annotate")
   run_process(_settings, "rm -f %s/Annotate/out/*.epsilon-nb_results.txt"%(_settings.rundir), "Annotate")
   run_process(_settings, "rm -f %s/Annotate/out/*.phymm.out"%(_settings.rundir), "Annotate")

   pool = Pool(processes=_settings.threads)
   tasks = []

   if "fcp" in _cls or "phymm" in _cls:
      # hack to use gridX
      if _USE_GRID:
         size = sizeFastaFile("%s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX))
         perThread = ceil(float(size) / 200)
         #print "The size of the contigs is %d per thread %d\n"%(size, perThread)
         #run_process(_settings, "python %s/python/splitfasta.py %s/Annotate/in/%s.asm.contig %d %s/Annotate/in/%s %d"%(_settings.METAMOS_UTILS, _settings.rundir, _settings.PREFIX, perThread, _settings.rundir, _settings.PREFIX, 1), "Annotate")
         splitfasta("%s/Annotate/in/%s.asm.contig"%(_settings.rundir,_settings.PREFIX),"%d"%(perThread),"%s/Annotate/in/%s"%(_settings.rundir,_settings.PREFIX),"%d"%(1))
         totalJobs = 0
         for partFile in os.listdir("%s/Annotate/in/"%(_settings.rundir)):
            if "_part" in partFile and "%s_part"%(_settings.PREFIX) in partFile:
               print "A file I have to process is %s\n"%(partFile)
               totalJobs += 1

         for lib in _readlibs:
            listOfFiles += ":%s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, lib.id)
            run_process(_settings, "ln %s/Assemble/out/lib%d.unaligned.fasta %s/Annotate/in/lib%d.unaligned.fasta"%(_settings.rundir, lib.id, _settings.rundir, lib.id), "Annotate")
            size = sizeFastaFile("%s/Annotate/in/lib%d.unaligned.fasta"%(_settings.rundir, lib.id))
            perThread = ceil(float(size) / 200)
            #print "The size of the lib %d is %d per one %d\n"%(lib.id, size, perThread)
            #run_process(_settings, "python %s/python/splitfasta.py %s/Annotate/in/lib%d.unaligned.fasta %d %s/Annotate/in/%s %d"%(_settings.METAMOS_UTILS, _settings.rundir, lib.id, perThread, _settings.rundir, _settings.PREFIX, totalJobs+1), "Annotate")
            #splitfasta("%s/Annotate/in/%s.asm.contig,%d,%s/Annotate/in/%s,%d"%(_settings.rundir,_settings.PREFIX,perThread,_settings.rundir,_settings.PREFIX,totalJobs+1))
            splitfasta("%s/Annotate/in/lib%d.unaligned.fasta"%(_settings.rundir,lib.id),"%d"%(perThread),"%s/Annotate/in/%s"%(_settings.rundir,_settings.PREFIX),"%d"%(totalJobs+1))

         totalJobs = 0
         for partFile in os.listdir("%s/Annotate/in/"%(_settings.rundir)):
            if "_part" in partFile and "%s_part"%(_settings.PREFIX) in partFile:
               #print "A file I have to process is %s\n"%(partFile)
               totalJobs += 1

         cmdfile = open("%s/Annotate/out/runAnnot.sh"%(_settings.rundir), "w")
         cmdfile.write("#!/bin/sh\n")
         cmdfile.write("\n")
         cmdfile.write("jobid=$GRID_TASK\n")
         cmdfile.write("if [ x$jobid = x -o x$jobid = xundefined -o x$jobid = 0 ]; then\n")
         cmdfile.write("   jobid=$1\n")
         cmdfile.write("fi\n")
         cmdfile.write("if test x$jobid = x; then\n")
         cmdfile.write("  echo Error: I need a job index on the command line\n")
         cmdfile.write("  exit 1\n")
         cmdfile.write("fi\n")
         cmdfile.write("if [ $jobid -gt %d ]; then\n"%(totalJobs))
         cmdfile.write("   echo Job id $jobid is out of range %d\n"%(totalJobs))
         cmdfile.write("   exit 0\n")
         cmdfile.write("fi\n")
         cmdfile.write("if test -e %s/Annotate/out/$jobid.success ; then\n"%(_settings.rundir))
         cmdfile.write("   echo Job previously completed successfully.\n")
         cmdfile.write("else\n")
         cmdfile.write("ln -s %s/.blastData\n"%(_settings.PHYMM))
         cmdfile.write("ln -s %s/.genomeData\n"%(_settings.PHYMM))
         cmdfile.write("ln -s %s/.scripts\n"%(_settings.PHYMM))
         cmdfile.write("ln -s %s/.taxonomyData\n"%(_settings.PHYMM))
         cmdfile.write("mkdir .logs\n")
         cmdfile.write("perl %s/scoreReads.pl %s/Annotate/in/%s_part$jobid.fa"%(_settings.PHYMM,_settings.rundir,_settings.PREFIX))
         cmdfile.write(" && touch %s/Annotate/out/$jobid.success\n"%(_settings.rundir))
         cmdfile.write("fi\n")
         cmdfile.close()
         run_process(_settings, "chmod u+x %s/Annotate/out/runAnnot.sh"%(_settings.rundir), "Annotate")

         #run_process(_settings, "gridx -p %d -r %d -T -c %s/Annotate/out/runAnnot.sh"%(min(totalJobs+1, 200), totalJobs+1, _settings.rundir), "Annotate")
         run_process(_settings, "cat %s/Annotate/out/gridx-ibissub00*/wrk_*/results.03.phymmBL_%s_Annotate_in_*%s* | grep -v \"QUERY_ID\" > %s/Annotate/out/%s.phymm.out"%(_settings.rundir, _settings.rundir.replace(os.sep, "_").replace(".", "_"), _settings.PREFIX, _settings.rundir, _settings.PREFIX))
         run_process(_settings, "ln %s/Annotate/out/%s.phymm.out %s/Annotate/out/%s.hits"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX))

         # for now we only work as phymm
         # generate Krona output ImportPhymmBL.pl
         importPhymm = "%s%sperl%sImportPhymmBL.pl"%(_settings.METAMOS_UTILS, os.sep, os.sep)
         if not os.path.exists(importPhymm):
            print "Error: Krona importer for Phymm not found in %s. Please check your path and try again.\n"%(importPhymm)
            raise(JobSignalledBreak)
         run_process(_settings, "perl %s %s -f %s %s/Annotate/out/%s.phymm.out:%s/Assemble/out/%s.contig.cnt:%s"%(importPhymm, "-l" if _settings.local_krona else "", listOfFiles,_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.taxa_level),"Annotate") # TODO: local url (after next KronaTools release)

         # generate taxonomic-level annots
         readctg_dict = {}
         for lib in _readlibs:
            ctgfile = open("%s/Assemble/out/%s.lib%dcontig.reads"%(_settings.rundir, _settings.PREFIX, lib.id), 'r')
            for line in ctgfile.xreadlines():
               line = line.replace("\n","")
               read, ctg = line.split()
               if ctg in readctg_dict:
                  readctg_dict[ctg].append(read)
               else:
                  readctg_dict[ctg] = [read,]
            ctgfile.close()

         annotsfile = open("%s/Annotate/out/%s.annots"%(_settings.rundir, _settings.PREFIX), 'r')
         annotreads = open("%s/Annotate/out/%s.reads.annots"%(_settings.rundir, _settings.PREFIX), 'w')
         for line in annotsfile.xreadlines():
            line = line.replace("\n", "")
            ctg, annot = line.split()
            if ctg in readctg_dict:
               for x in readctg_dict[ctg]:
                  annotreads.write("%s\t%s\n"%(x, annot))
            else:
               annotreads.write("%s\t%s\n"%(ctg, annot))
         annotsfile.close()
         annotreads.close()
         readctg_dict.clear()

         return 
      # we should also split the fna and faa file but for now this is good enough
      size = sizeFastaFile("%s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX))
      perThread = max(ceil(float(size) / _settings.threads), _MIN_SEQ_LENGTH)
      #print "The size of the contigs is %d per thread %d\n"%(size, perThread)
      #run_process(_settings, "python %s/python/splitfasta.py %s/Annotate/in/%s.asm.contig %d"%(_settings.METAMOS_UTILS, _settings.rundir, _settings.PREFIX, perThread), "Annotate")
      #splitfasta("%s/Annotate/in/%s.asm.contig,%d,%s/Annotate/in/%s,%d"%(_settings.rundir,_settings.PREFIX,perThread,_settings.rundir,_settings.PREFIX,1))
      splitfasta("%s/Annotate/in/%s.asm.contig"%(_settings.rundir,_settings.PREFIX),"%d"%(perThread))
      for partFile in os.listdir("%s/Annotate/in/"%(_settings.rundir)):
         if "_part" in partFile and "%s.asm.contig"%(_settings.PREFIX) in partFile:
            partStart = partFile.find("_part")+5
            partEnd = partFile.find(".fa", partStart, len(partFile))
            partNumber = int(partFile[partStart:partEnd])
            params = {}
            params["jobID"] = len(tasks) 
            params["cls"] = _cls
            params["contigs"] = "%s/Annotate/in/%s"%(_settings.rundir, partFile)
            params["orfAA"] = ""
            params["orfFA"] = ""
            params["out"] = "%s.ctg_%d"%(_settings.PREFIX, partNumber)
            tasks.append(params) 
   else:
      annotateSeq(_cls, "%s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX), "%s/Annotate/in/%s.faa"%(_settings.rundir, _settings.PREFIX), "%s/Annotate/in/%s.fna"%(_settings.rundir, _settings.PREFIX), "%s.ctg"%(_settings.PREFIX))

   # annotate all the unmapped sequences using FCP
   if _cls == "blast" or _cls == "phmmer" or not _settings.annotate_unmapped:
      #print "Warning: blast and PHMMER is not supported for annotating unmapped sequences"
      #print "Warning: unmapped/unaligned sequences will not be annotated!"
      pass
   else:
      for lib in _readlibs:
         listOfFiles += ":%s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, lib.id)
         run_process(_settings, "ln %s/Assemble/out/lib%d.unaligned.fasta %s/Annotate/in/lib%d.unaligned.fasta"%(_settings.rundir, lib.id, _settings.rundir, lib.id), "Annotate")

         if "fcp" in _cls or "phymm" in _cls:
            size = sizeFastaFile("%s/Annotate/in/lib%d.unaligned.fasta"%(_settings.rundir, lib.id))
            perThread = max(ceil(float(size) / _settings.threads), _MIN_SEQ_LENGTH)
            #run_process(_settings, "python %s/python/splitfasta.py %s/Annotate/in/lib%d.unaligned.fasta %d"%(_settings.METAMOS_UTILS, _settings.rundir, lib.id, perThread), "Annotate")
            splitfasta("%s/Annotate/in/lib%d.unaligned.fasta"%(_settings.rundir,lib.id),"%d"%(perThread))
            for partFile in os.listdir("%s/Annotate/in/"%(_settings.rundir)):
               if "_part" in partFile and "lib%d.unaligned.fasta"%(lib.id) in partFile:
                  partStart = partFile.find("_part")+5
                  partEnd = partFile.find(".fa", partStart, len(partFile))
                  partNumber = int(partFile[partStart:partEnd])
                  params = {}
                  params["jobID"] = len(tasks)
                  params["cls"] = _cls
                  params["contigs"] = "%s/Annotate/in/%s"%(_settings.rundir, partFile)
                  params["orfAA"] = ""
                  params["orfFA"] = ""
                  params["out"] = "%s.lib%d_%d"%(_settings.PREFIX, lib.id, partNumber)
                  tasks.append(params)
         else:
            annotateSeq(_cls, "%s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, lib.id), "", "", "%s.lib%d"%(_settings.PREFIX, lib.id))
   if "fcp" in _cls or "phymm" in _cls:
         result = pool.map_async(parallelWrapper, tasks).get(sys.maxint)
         for i in result:
            if (i["status"] == 1):
               run_process(_settings, "rm %s"%(tasks[i["jobID"]]["contigs"]), "Annotate")
            else:
               print "Error: parallel annotation job %d failed\n"%(i["jobID"])
               raise(JobSignalledBreak)
   pool.close()
   pool.join()

   if generic.checkIfExists(STEP_NAMES.ANNOTATE, _cls.lower()):
      generic.execute(STEP_NAMES.ANNOTATE, _cls.lower())
   else:
      #  merge results
      run_process(_settings, "cat %s/Annotate/out/*.intermediate.hits > %s/Annotate/out/%s.hits"%(_settings.rundir, _settings.rundir, _settings.PREFIX), "Annotate")
 
   if _cls == "phylosift":
       importPS = "%s%sperl%sImportPhyloSift.pl"%(_settings.METAMOS_UTILS, os.sep, os.sep)
       if not os.path.exists(importPS):
           print "Error: Krona importer for PhyloSift not found in %s. Please check your path and try again.\n"%(_settings.KRONA)
           raise(JobSignalledBreak)
       run_process(_settings, "perl %s %s -c -i -f %s %s/Annotate/out/%s.hits:%s/Assemble/out/%s.contig.cnt:%s"%(importPS,"-l" if _settings.local_krona else "",listOfFiles,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX, _settings.taxa_level), "Annotate")

   elif _cls == "fcp":
       # generate Krona output
       importFCP = "%s%sperl%sImportFCP.pl"%(_settings.METAMOS_UTILS, os.sep, os.sep)
       if not os.path.exists(importFCP):
          print "Error: Krona importer for FCP not found in %s. Please check your path and try again.\n"%(importFCP)
          raise(JobSignalledBreak)
       run_process(_settings, "cat %s/Annotate/out/*.intermediate.epsilon-nb_results.txt | grep -v 'Fragment Id' > %s/Annotate/out/%s.epsilon-nb_results.txt"%(_settings.rundir, _settings.rundir, _settings.PREFIX), "Annotate")

       run_process(_settings, "perl %s %s -c -i -f %s %s/Annotate/out/%s.epsilon-nb_results.txt:%s/Assemble/out/%s.contig.cnt:%s"%(importFCP, "-l" if _settings.local_krona else "", listOfFiles, _settings.rundir,_settings.PREFIX,_settings.rundir, _settings.PREFIX, _settings.taxa_level),"Annotate") # TODO: local url (after next KronaTools release)

   elif _cls == "phymm":
       # generate Krona output ImportPhymmBL.pl
       importPhymm = "%s%sperl%sImportPhymmBL.pl"%(_settings.METAMOS_UTILS, os.sep, os.sep)
       if not os.path.exists(importPhymm):
          print "Error: Krona importer for Phymm not found in %s. Please check your path and try again.\n"%(importPhymm)
          raise(JobSignalledBreak)
       run_process(_settings, "cat %s/Annotate/out/*.intermediate.phymm.out > %s/Annotate/out/%s.phymm.out"%(_settings.rundir, _settings.rundir, _settings.PREFIX), "Annotate")
       run_process(_settings, "perl %s %s -f %s %s/Annotate/out/%s.phymm.out:%s/Assemble/out/%s.contig.cnt:%s"%(importPhymm, "-l" if _settings.local_krona else "", listOfFiles,_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.taxa_level),"Annotate") # TODO: local url (after next KronaTools release)
   elif generic.checkIfExists(STEP_NAMES.ANNOTATE, _cls.lower()):
      genericImport = "%s%sperl%sImport%s.pl"%(_settings.METAMOS_UTILS, os.sep, os.sep, _cls.title())
      if os.path.exists(genericImport):
         run_process(_settings, "perl %s %s -c -i -f %s %s/Annotate/out/%s.hits:%s/Assemble/out/%s.contig.cnt:%s"%(genericImport, "-l" if _settings.local_krona else "", listOfFiles, _settings.rundir,_settings.PREFIX,_settings.rundir, _settings.PREFIX, _settings.taxa_level),"Annotate") # TODO: local url (after next KronaTools release)
      else:
         genericImport = "%s%sperl%sImportGeneric.pl"%(_settings.METAMOS_UTILS, os.sep, os.sep)
         if not os.path.exists(genericImport):
            print "Error: Krona importer for generic classifier not found in %s. Please check your path and try again.\n"%(genericImport)
            raise(JobSignalledBreak)
         run_process(_settings, "perl %s %s -c -i -f %s %s/Annotate/out/%s.hits:%s/Assemble/out/%s.contig.cnt:%s"%(genericImport, "-l" if _settings.local_krona else "", listOfFiles, _settings.rundir,_settings.PREFIX,_settings.rundir, _settings.PREFIX, _settings.taxa_level),"Annotate") # TODO: local url (after next KronaTools release)


   run_process(_settings, "unlink %s/Postprocess/in/%s.hits"%(_settings.rundir, _settings.PREFIX), "Annotate")
   run_process(_settings, "unlink %s/Postprocess/out/%s.hits"%(_settings.rundir, _settings.PREFIX), "Annotate")
   run_process(_settings, "ln %s/Annotate/out/%s.hits %s/Postprocess/in/%s.hits"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Annotate")
   run_process(_settings, "ln %s/Annotate/out/%s.hits %s/Postprocess/out/%s.hits"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Annotate")

   # generate taxonomic-level annots
   readctg_dict = {}
   for lib in _readlibs:
      ctgfile = open("%s/Assemble/out/%s.lib%dcontig.reads"%(_settings.rundir, _settings.PREFIX, lib.id), 'r')
      for line in ctgfile.xreadlines():
         line = line.replace("\n","")
         read, ctg = line.split()
         if ctg in readctg_dict:
            readctg_dict[ctg].append(read)
         else:
            readctg_dict[ctg] = [read,]
      ctgfile.close()

   annotsfile = open("%s/Annotate/out/%s.annots"%(_settings.rundir, _settings.PREFIX), 'r')
   annotreads = open("%s/Annotate/out/%s.reads.annots"%(_settings.rundir, _settings.PREFIX), 'w')
   for line in annotsfile.xreadlines():
     line = line.replace("\n", "")
     ctg, annot = line.split()
     if ctg in readctg_dict:
        for x in readctg_dict[ctg]:
           annotreads.write("%s\t%s\n"%(x, annot))
     else:
        annotreads.write("%s\t%s\n"%(ctg, annot))
   annotsfile.close()
   annotreads.close()
   readctg_dict.clear()
