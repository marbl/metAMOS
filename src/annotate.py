#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from findreps import FindRepeats
#from findorfs import FindORFS

sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_settings = Settings()
_cls = None

def init(reads, skipsteps, cls):
   global _readlibs
   global _skipsteps
   global _cls

   _readlibs = reads
   _skipsteps = skipsteps
   _cls = cls

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

       if not os.path.exists("%s/DB/allprots.faa"%(_settings.METAMOS_UTILS)):
          print "Error: You indicated you would like to run phmmer but DB allprots.faa not found in %s/DB. Please check your path and try again.\n"%(_settings.METAMOS_UTILS)
          raise(JobSignalledBreak)

       run_process(_settings, "%s/phmmer --cpu %d -E 0.0000000000000001 -o %s/Annotate/out/%s.phm.out --tblout %s/Annotate/out/%s.phm.tbl --notextw %s %s/DB/allprots.faa"%(_settings.PHMMER, _settings.threads,_settings.rundir,output,_settings.rundir,output,orfAA,_settings.METAMOS_UTILS),"Annotate")
       parse_phmmerout("%s/Annotate/out/%s.phm.tbl"%(_settings.rundir,output))
       run_process(_settings, "mv %s/Annotate/out/%s.phm.tbl  %s/Annotate/out/%s.hits"%(_settings.rundir,output,_settings.rundir,output),"Annotate")

   elif cls == "metaphyler":
       if not os.path.exists(_settings.BLAST + os.sep + "blastall"):
           print "Error: BLAST not found in %s. Please check your path and try again.\n"%(_settings.BLAST)
           raise(JobSignalledBreak)

       if not os.path.exists("%s/DB/allprots.faa"%(_settings.METAMOS_UTILS)):
          print "Error: You indicated you would like to run MetaPhyler but it depends on BLAST and DB allprots.faa not found in %s/DB. Please check your path and try again.\n"%(_settings.METAMOS_UTILS)
          raise(JobSignalledBreak)

       #run_process(_settings, "perl %s/perl/installMetaphyler.pl"%(_settings.METAMOS_UTILS)
       run_process(_settings,"%s/blastall -p blastx -a %d -m 8 -b 1 -e 1e-2 -i %s -d %s/perl/metaphyler/test/test.ref.protein > %s/Annotate/out/%s.query.blastx"%(_settings.BLAST,_settings.threads,orfFA,_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX))
       run_process(_settings, "%s/metaphylerClassify %s/perl/metaphyler/markers/markers.blastx.classifier %s/perl/metaphyler/markers/markers.taxonomy %s/Annotate/out/%s.query.blastx > %s/Annotate/out/%s.classification"%(_settings.METAPHYLER,_settings.METAMOS_UTILS,_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX,_settings.rundir,output) )
       #Krona import Metaphyler script!
       run_process(_settings, "mv %s/Annotate/out/%s.classification  %s/Annotate/out/%s.hits"%(_settings.rundir,output,_settings.rundir,output),"Annotate")

   elif cls == "blast":
       if not os.path.exists(_settings.BLAST + os.sep + "blastall"):
          print "Error: BLAST not found in %s. Please check your path and try again.\n"%(_settings.BLAST)
          raise(JobSignalledBreak)

       if not os.path.exists("%s/DB/allprots.faa"%(_settings.METAMOS_UTILS)):
          print "Error: You indicated you would like to run BLAST but DB allprots.faa not found in %s/DB. Please check your path and try again.\n"%(_settings.METAMOS_UTILS)
          raise(JobSignalledBreak)
       run_process(_settings, "%s/blastall -v 1 -b 1 -a %d -p blastp -m 8 -e 0.00001 -i %s -d %s/DB/refseq_protein -o %s/Annotate/out/%s.blastout"%(_settings.BLAST, _settings.threads,orfAA,_settings.METAMOS_UTILS,_settings.rundir,output),"Annotate")
       run_process(_settings, "mv %s/Annotate/out/%s.blastout  %s/Annotate/out/%s.hits"%(_settings.rundir,output,_settings.rundir,output),"Annotate")
   elif cls == "phylosift":
       if _settings.PHYLOSIFT == "" or not os.path.exists(_settings.PHYLOSIFT + os.sep + "bin" + os.sep + "phylosift"):
          print "Error: PhyloSift not found in %s. Please check your path and try again.\n"%(_settings.PHYLOSIFT)
          raise(JobSignalledBreak)

       phylosiftCmd =  "%s/bin/phylosift all --threaded=%d"%(_settings.PHYLOSIFT, _settings.threads)
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
       run_process(_settings, "unlink %s/Annotate/out/%s.hits"%(_settings.rundir, output), "Annotate")
       run_process(_settings, "ln %s/Annotate/out/PS_temp/%s/sequence_taxa_summary.txt %s/Annotate/out/%s.hits"%(_settings.rundir, os.path.basename(contigs), _settings.rundir, output), "Annotate") 
       
   elif cls == "fcp":
       run_process(_settings, "ln -s %s/models"%(_settings.METAMOS_UTILS), "Annotate")
       run_process(_settings, "ln -s %s/models/taxonomy.txt"%(_settings.METAMOS_UTILS), "Annotate")

       # normal nb classify
       run_process(_settings, "%s/nb-classify -q %s -m %s/models/models.txt -r %s/Annotate/out/%s.nb_results.txt"%(_settings.FCP,contigs,_settings.METAMOS_UTILS,_settings.rundir,output), "Annotate")

       # for blast options
       if os.path.exists("%s/DB/blast_data/BacteriaAndArchaeaGenomesDB.nin"%(_settings.METAMOS_UTILS)):
          run_process(_settings, "ln -s %s/DB/blast_data blast_data"%(_settings.METAMOS_UTILS), "Annotate")
          run_process(_settings, "python %s/python/BLASTN.py %s/blastn %s %s/Annotate/out/%s.bl_results.txt"%(_settings.METAMOS_UTILS, _settings.BLAST, contigs, _settings.rundir, output), "Annotate") 

          #combine the results
          run_process(_settings, "python %s/python/NB-BL.py %s/Annotate/out/%s.nb_results.txt %s/Annotate/out/%s.bl_results.txt %s/Annotate/out/%s.epsilon-nb_results.txt"%(_settings.METAMOS_UTILS, _settings.rundir, output, _settings.rundir, output, _settings.rundir, output), "Annotate")
       else:
          run_process(_settings, "python %s/python/Epsilon-NB.py %s/Annotate/out/%s.nb_results.txt 1E5 %s/Annotate/out/%s.epsilon-nb_results.txt"%(_settings.METAMOS_UTILS,_settings.rundir,output,_settings.rundir,output),"Annotate")

       #need python TaxonomicSummary.py test.fasta nb_topModels.txt nb_taxonomicSummary.txt
       #run_process(_settings, "python %s/python/TaxonomicSummary.py %s/Annotate/in/%s.fna %s/Annotate/out/%s.nb_results.txt %s/Annotate/out/%s.epsilon-nb_results.txt"%(_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Annotate")

       run_process(_settings, "unlink %s/Annotate/out/%s.hits"%(_settings.rundir, output), "Annotate")
       run_process(_settings, "ln %s/Annotate/out/%s.epsilon-nb_results.txt %s/Annotate/out/%s.hits"%(_settings.rundir, output, _settings.rundir, output), "Annotate")

   elif cls == "phymm":
       if not os.path.exists("%s"%(_settings.PHYMM)):
           print "Error: Phymm not found in %s but selected as classifier. Please check your path and try again.\n"%(_settings.PHYMM)
           raise(JobSignalledBreak)
       #run_process(_settings, "ln -s %s.genomeData"%(_settings.PHYMM), "Annotate")
       #run_process(_settings, "ln -s \"%s.blastData\""%(_settings.PHYMM), "Annotate")
       #run_process(_settings, "ln -s \"%s.taxonomyData\""%(_settings.PHYMM), "Annotate")
       #run_process(_settings, "ln -s \"%s.scripts\""%(_settings.PHYMM), "Annotate")
       #run_process(_settings, "ln -s \"%s.logs\""%(_settings.PHYMM), "Annotate")
       #run_process(_settings, "ln -s %s .logs"%(_settings.PHYMM), "Annotate")
       os.chdir("%s"%(_settings.PHYMM))
       #print os.getcwd()
       run_process(_settings, "perl scoreReads.pl %s > %s/Annotate/out/%s.phymm.out"%(contigs,_settings.rundir,output),"")
       os.chdir("%s"%(_settings.rundir))
       #perl scoreReads.pl mammoth17.fna > mammoth17.fna.out &
       
   elif cls == None:
       print "No method specified, skipping"


@follows(FindRepeats)
@posttask(touch_file("%s/Logs/annotate.ok"%(_settings.rundir)))
@files("%s/Annotate/in/%s.faa"%(_settings.rundir,_settings.PREFIX),"%s/Annotate/out/%s.hits"%(_settings.rundir,_settings.PREFIX))
def Annotate(input,output):
   if "Annotate" in _skipsteps or _cls == None:
      run_process(_settings, "touch %s/Logs/annotate.skip"%(_settings.rundir), "Annotate")
      run_process(_settings, "touch %s/Annotate/out/%s.hits"%(_settings.rundir, _settings.PREFIX), "Annotate")
      return 0

   listOfFiles = "%s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX)

   run_process(_settings, "touch %s/Annotate/out/%s.annots"%(_settings.rundir, _settings.PREFIX), "Annotate")
   run_process(_settings, "unlink %s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX), "Annotate")
   run_process(_settings, "ln -s %s/Assemble/out/%s.asm.contig %s/Annotate/in/"%(_settings.rundir, _settings.PREFIX, _settings.rundir), "Annotate")
   run_process(_settings, "unlink %s/Annotate/out/%s.hits"%(_settings.rundir, _settings.PREFIX), "Annotate")

   annotateSeq(_cls, "%s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX), "%s/Annotate/in/%s.faa"%(_settings.rundir, _settings.PREFIX), "%s/Annotate/in/%s.fna"%(_settings.rundir, _settings.PREFIX), "%s.ctg"%(_settings.PREFIX))

   # annotate all the unmapped sequences using FCP
   if _cls == "blast" or _cls == "phmmer":
      print "Warning: blast and PHMMER is not supported for annotating unmapped sequences"
   else:
      for lib in _readlibs:
         listOfFiles += ":%s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, lib.id)
         annotateSeq(_cls, "%s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, lib.id), "", "", "%s.lib%d"%(_settings.PREFIX, lib.id))

   # merge results
   run_process(_settings, "cat %s/Annotate/out/*.hits > %s/Annotate/out/%s.hits"%(_settings.rundir, _settings.rundir, _settings.PREFIX), "Annotate")
 
   if _cls == "phylosift":
       importPS = "%s%sperl%sImportPhyloSift.pl"%(_settings.METAMOS_UTILS, os.sep, os.sep)
       if not os.path.exists(importPS):
           print "Error: Krona importer for PhyloSift not found in %s. Please check your path and try again.\n"%(_settings.KRONA)
           raise(JobSignalledBreak)
       run_process(_settings, "perl %s -c -v -i -f %s %s/Annotate/out/%s.hits:%s/Assemble/out/%s.contig.cnt:%s"%(importPS,listOfFiles,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX, _settings.taxa_level), "Annotate")

   elif _cls == "fcp":
       # generate Krona output
       importFCP = "%s%sperl%sImportFCP.pl"%(_settings.METAMOS_UTILS, os.sep, os.sep)
       if not os.path.exists(importFCP):
          print "Error: Krona importer for FCP not found in %s. Please check your path and try again.\n"%(importFCP)
          raise(JobSignalledBreak)
       run_process(_settings, "cat %s/Annotate/out/*.epsilon-nb_results.txt | grep -v 'Fragment Id' > %s/Annotate/out/%s.epsilon-nb_results.txt"%(_settings.rundir, _settings.rundir, _settings.PREFIX), "Annotate")

       run_process(_settings, "perl %s -c -v -i -f %s %s/Annotate/out/%s.epsilon-nb_results.txt:%s/Assemble/out/%s.contig.cnt:%s"%(importFCP, listOfFiles, _settings.rundir,_settings.PREFIX,_settings.rundir, _settings.PREFIX, _settings.taxa_level),"Annotate") # TODO: local url (after next KronaTools release)

   elif _cls == "phymm":
       # generate Krona output ImportPhymmBL.pl
       importPhymm = "%s%skrona%sImportPhymmBL.pl"%(_settings.METAMOS_UTILS, os.sep, os.sep)
       if not os.path.exists(importFCP):
          print "Error: Krona importer for Phymm not found in %s. Please check your path and try again.\n"%(importPHYMM)
          raise(JobSignalledBreak)
       #run_process(_settings, "cat %s/Annotate/out/*.epsilon-nb_results.txt | grep -v 'Fragment Id' > %s/Annotate/out/%s.epsilon-nb_results.txt"%(_settings.rundir, _settings.rundir, _settings.PREFIX), "Annotate")

       run_process(_settings, "perl %s -c -v -i %s/Annotate/out/%s.phymm.out:%s/Assemble/out/%s.contig.cnt:%s"%(importPHYMM, _settings.rundir,_settings.PREFIX,_settings.rundir, _settings.PREFIX, _settings.taxa_level),"Annotate") # TODO: local url (after next KronaTools release)

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
