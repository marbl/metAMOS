#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from mapreads import SplitForORFs

sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_settings = Settings()
_orf = None 
_min_ctg_len = 300
_min_ctg_cvg = 3
_read_orfs = False
_run_fast = True

def setRunFast(fast):
   global _run_fast
   _run_fast = fast

def init(reads, skipsteps, orf, min_ctg_len, min_ctg_cvg,read_orfs):
   global _readlibs
   global _skipsteps
   global _orf
   global _min_ctg_cvg
   global _min_ctg_len
   global _read_orfs
   _readlibs = reads
   _skipsteps = skipsteps
   _orf = orf
   _min_ctg_cvg = min_ctg_cvg
   _min_ctg_len = min_ctg_len
   _read_orfs = read_orfs

def parse_genemarkout(orf_file,is_scaff=False, error_stream="FindORFS",min_len=_min_ctg_len,min_cvg=_min_ctg_cvg):
    coverageFile = open("%s/Assemble/out/%s.contig.cvg"%(_settings.rundir, _settings.PREFIX), 'r')
    cvg_dict = {} 
    for line in coverageFile:
        data = line.split()
        cvg_dict[data[0]] = float(data[1])
    coverageFile.close()

    coords = open(orf_file,'r')
    coords.readline()
#    outf = open("proba.orfs",'w')
    prevhdr = 0
    prevhdraa = 0
    prevhdrnt = 0

    curcontig = ""
    curseqaa = ""
    curseqnt = ""
    reads = {}
    gene_dict = {}
    fna_dict = {}
    for line in coords:
        if ">gene" in line[0:10]:
            if "_nt|" in line:
                #print prevhdraa, prevhdrnt#, curseqaa, curseqnt
                if prevhdraa and curseqaa != "":
                    try:
                        gene_dict[curcontig][curseqaa] = 1
                    except KeyError:
                        gene_dict[curcontig] = {}
                        gene_dict[curcontig][curseqaa] =1
                    curseqaa = ""

                elif prevhdrnt and curseqnt != "":
                    try:
                        fna_dict[curcontig][curseqnt] = 1
                    except KeyError:
                        fna_dict[curcontig] = {}
                        fna_dict[curcontig][curseqnt] = 1
                    curseqnt = ""

                prevhdrnt = 1
                prevhdraa = 0

            elif "_aa|" in line:

                if prevhdrnt and curseqnt != "":
                    try:
                        fna_dict[curcontig][curseqnt] = 1
                    except KeyError:
                        fna_dict[curcontig] = {}
                        fna_dict[curcontig][curseqnt] = 1
                    curseqnt = ""
                elif prevhdraa and curseqaa != "":
                    try:
                        gene_dict[curcontig][curseqaa] = 1
                    except KeyError:
                        gene_dict[curcontig] = {}
                        gene_dict[curcontig][curseqaa] = 1
                    curseqaa = ""
                prevhdraa = 1
                prevhdrnt = 0

            prevhdr = 1
            lined = line.replace("\n","")
            data = line[1:].split(">",1)[1]
            
            curcontig = data.split(" ")[0]
            if len(data.split(" ")) == 1:
                curcontig = data.split("\t")[0]
            curcontig = curcontig.strip()
            #print curcontig, len(curcontig)
            prevhdr = 1

        elif len(line) > 2 and prevhdraa == 1 and prevhdr:
            curseqaa += line
        elif len(line) > 2 and prevhdrnt == 1 and prevhdr:
            curseqnt += line
        elif len(line) <= 2 or "Nucleotide" in line: #and prevhdr == 1:
            prevhdr = 0
            #prevhdraa = 0
            #prevhdrnt = 0

        else:
            continue
    if prevhdraa and curseqaa != "":
        try:
          gene_dict[curcontig][curseqaa] = 1
        except KeyError:
          gene_dict[curcontig] = {}
          gene_dict[curcontig][curseqaa] = 1
          curseqaa = ""

    elif prevhdrnt and curseqnt != "":
        try:
          fna_dict[curcontig][curseqnt] = 1
        except KeyError:
          fna_dict[curcontig] = {}
          fna_dict[curcontig][curseqnt] = 1
    if is_scaff:
        outf = open("%s/FindScaffoldORFS/out/%s.faa"%(_settings.rundir,_settings.PREFIX),'w')
        outf2 = open("%s/FindScaffoldORFS/out/%s.fna"%(_settings.rundir,_settings.PREFIX),'w')
        #cvgf = open("%s/FindScaffoldORFS/out/%s.contig.cvg"%(_settings.rundir,_settings.PREFIX),'w')
        cvgg = open("%s/FindScaffoldORFS/out/%s.gene.cvg"%(_settings.rundir,_settings.PREFIX),'w')
        genectg = open("%s/FindScaffoldORFS/out/%s.gene.map"%(_settings.rundir, _settings.PREFIX), 'w')
    else:
        outf = open("%s/FindORFS/out/%s.faa"%(_settings.rundir,_settings.PREFIX),'w')
        outf2 = open("%s/FindORFS/out/%s.fna"%(_settings.rundir,_settings.PREFIX),'w')
        #cvgf = open("%s/FindORFS/out/%s.contig.cvg"%(_settings.rundir,_settings.PREFIX),'w')
        cvgg = open("%s/FindORFS/out/%s.gene.cvg"%(_settings.rundir,_settings.PREFIX),'w')
        genectg = open("%s/FindORFS/out/%s.gene.map"%(_settings.rundir, _settings.PREFIX), 'w')
    #print len(gene_dict.keys())
    orfs = {}

    genecnt = 1
    for key in gene_dict.keys():
        genecnt = 1
        for gene in gene_dict[key].keys():
            if not is_scaff:
                if key in cvg_dict.keys():
                    cvgg.write("%s_gene%d\t%.2f\t%.2f\t%.2f\n"%(key,genecnt,cvg_dict[key]*len(gene),len(gene),cvg_dict[key]))                         
                else:
                    cvgg.write("%s_gene%d\t%.2f\t%.2f\t%.2f\n"%(key,genecnt, 1.0,len(gene),1.0))

            #min aa length, read depth
            if key in cvg_dict:
                if len(gene) < min_len/3 or cvg_dict[key] < min_cvg:# or cvg_dict[key] < 5:
                   continue
            else:
                if len(gene) < min_len/3 or 1 > min_cvg: 
                   continue
            try:
                #print "contig"+key
                orfs["%s"%(key)] +=1
            except KeyError:
                orfs["%s"%(key)] =1
            outf.write(">%s_gene%d\n%s"%(key,genecnt,gene))
            genectg.write("%s\t%s_gene%d\n"%(key, key, genecnt))

            genecnt +=1
    for key in fna_dict.keys():
        genecnt = 1
        for gene in fna_dict[key].keys():
            #gene = fna_dict[key][gkey]
            if key in cvg_dict:
                if len(gene) < min_len or cvg_dict[key] < min_cvg:# or cvg_dict[key] < 5:
                    continue
            else:
                if len(gene) < min_len or 1 > min_cvg:
                   continue;
            outf2.write(">%s_gene%d\n%s"%(key,genecnt,gene))
            genecnt +=1

#        print gene_dict[key][0]
    outf.close()
    outf2.close()
    genectg.close()
    cvgg.close()

def parse_fraggenescanout(orf_file,is_scaff=False, error_stream="FindORFS",min_len=_min_ctg_len,min_cvg=_min_ctg_cvg):
    coverageFile = open("%s/Assemble/out/%s.contig.cvg"%(_settings.rundir, _settings.PREFIX), 'r')
    cvg_dict = {} 
    len_dict = {}

    for line in coverageFile:
        data = line.split()
        cvg_dict[data[0]] = float(data[1])
    coverageFile.close()
    genefile = ""
    if is_scaff:
        genefile = open("%s/FindScaffoldORFS/out/%s.orfs.ffn"%(_settings.rundir,_settings.PREFIX),'r')
        cvgg = open("%s/FindScaffoldORFS/out/%s.gene.cvg"%(_settings.rundir,_settings.PREFIX),'w')
        genectg = open("%s/FindScaffoldORFS/out/%s.gene.map"%(_settings.rundir, _settings.PREFIX), 'w')
    else:
        genefile = open("%s/FindORFS/out/%s.orfs.ffn"%(_settings.rundir,_settings.PREFIX),'r')
        genectg = open("%s/FindORFS/out/%s.gene.map"%(_settings.rundir, _settings.PREFIX), 'w')
        cvgg = open("%s/FindORFS/out/%s.gene.cvg"%(_settings.rundir,_settings.PREFIX),'w')
    orfs = {}
  
    data = genefile.read()
    seqs = data.split(">")[1:]
    gene_ids = []
    orfhdrs = {}
    for seq in seqs:
        hdr,gene = seq.split("\n",1)
        hdr = hdr.split("\n")[0]
        gene_ids.append(hdr)
        len_dict[hdr] = len(seq)

    for seq in seqs:
       hdr,gene = seq.split("\n",1)
       #hdr = hdr.split("\n")[0]
       hdr = hdr.rstrip("\n")
       #gene_ids.append(hdr)
       #split the header in two
       orfkey = '_'.join(hdr.split('_')[:1])
       #orfval = '_'.join(hdr.split('_')[2:])
       orfhdrs[orfkey]=hdr

       genectg.write("%s\t%s\n"%(orfkey, hdr))
       if orfkey in cvg_dict:
           if len_dict[hdr] > min_len and cvg_dict[orfkey] >=  min_cvg:
               cvgg.write("%s\t%.2f\n"%(hdr,len_dict[hdr]*cvg_dict[orfkey]))
       else:
           if len_dict[hdr] > min_len and min_cvg >= 1:
               cvgg.write("%s\t%s\n"%(hdr,str(1.0)))

    cvgg.close()
    genectg.close()
    #for key in gene_ids:
    #    genecnt = 1
    #    gkey = ""
    #    if not is_scaff:
    #        for ckey in cvg_dict.keys():
    #            if ckey in key:
    #                gkey = ckey
    #        
    #        if gkey != "":
    #            cvgg.write("%s\t%s\n"%(key,cvg_dict[gkey])) 
    #        else:
    #            cvgg.write("%s\t%s\n"%(key,1.0))
    #    genecnt +=1
    #cvgg.close()

def parse_prokka(orf_file,is_scaff=False, error_stream="FindORFS",min_len=_min_ctg_len,min_cvg=_min_ctg_cvg):
    coverageFile = open("%s/Assemble/out/%s.contig.cvg"%(_settings.rundir, _settings.PREFIX), 'r')
    cvg_dict = {}
    len_dict = {}

    for line in coverageFile:
        data = line.split()
        cvg_dict[data[0]] = float(data[1])
    coverageFile.close()
    genefile = ""
    if is_scaff:
        genefile = open("%s/FindScaffoldORFS/out/prokka/%s.gff"%(_settings.rundir,_settings.PREFIX),'r')
        genectg = open("%s/FindScaffoldORFS/out/%s.gene.map"%(_settings.rundir, _settings.PREFIX), 'w')
        cvgg = open("%s/FindScaffoldORFS/out/%s.gene.cvg"%(_settings.rundir,_settings.PREFIX),'w')
    else:
        genefile = open("%s/FindORFS/out/%s.prokka/%s.gff"%(_settings.rundir,_settings.PREFIX,_settings.PREFIX),'r')
        genectg = open("%s/FindORFS/out/%s.gene.map"%(_settings.rundir, _settings.PREFIX), 'w')
        cvgg = open("%s/FindORFS/out/%s.gene.cvg"%(_settings.rundir,_settings.PREFIX),'w')

    for line in genefile.xreadlines():
       if "#" in line:
          continue
       elif line.startswith(">"):
          break
       elif len(line.strip().split("\t")) < 9:
          continue
       (ctgID, geneA, geneB, start, end, index, geneori, geneC, id) = line.strip().split("\t", 9)
       id = id.split(";")[0].replace("ID=", "")
       genectg.write("%s\t%s\n"%(ctgID, id))
       try:
          cvgg.write("%s\t%s\n"%(id, cvg_dict[ctgID]))
       except KeyError:
          cvgg.write("%s\t%s\n"%(id, 1))
    genectg.close()
    cvgg.close()
    genefile.close()

def findFastaORFs(orf, contigs, outputFNA, outputFAA, outputCVG, outputMAP, min_len, min_cvg):
   if orf == "metagenemark":
       if not os.path.exists(_settings.METAGENEMARK + os.sep + "gmhmmp"):
          print "Error: MetaGeneMark not found in %s. Please check your path and try again.\n"%(_settings.METAGENEMARK)
          raise(JobSignalledBreak)
       run_process(_settings, "%s/gmhmmp -o %s/FindORFS/out/%s.orfs -m %s/config/MetaGeneMark_v1.mod -d -a %s"%(_settings.METAGENEMARK,_settings.rundir,_settings.PREFIX,_settings.METAMOS_UTILS,contigs),"FindORFS")
       parse_genemarkout("%s/FindORFS/out/%s.orfs"%(_settings.rundir,_settings.PREFIX), 0, "FindORFS", min_len, min_cvg)
       run_process(_settings, "mv %s/FindORFS/out/%s.fna %s/FindORFS/out/%s"%(_settings.rundir, _settings.PREFIX, _settings.rundir, outputFNA), "FindORFS")
       run_process(_settings, "mv %s/FindORFS/out/%s.faa %s/FindORFS/out/%s"%(_settings.rundir, _settings.PREFIX, _settings.rundir, outputFAA), "FindORFS")
       run_process(_settings,"mv %s/FindORFS/out/%s.gene.cvg %s/FindORFS/out/%s"%(_settings.rundir, _settings.PREFIX, _settings.rundir, outputCVG), "FindORFS")
       run_process(_settings,"mv %s/FindORFS/out/%s.gene.map %s/FindORFS/out/%s"%(_settings.rundir, _settings.PREFIX, _settings.rundir, outputMAP), "FindORFS")
   elif orf == "fraggenescan":
       if not os.path.exists(_settings.FRAGGENESCAN + os.sep + "FragGeneScan"):
          print "Error: FragGeneScan not found in %s. Please check your path and try again.\n"%(_settings.FRAGGENESCAN)
          raise(JobSignalledBreak)
       run_process(_settings,"%s/FragGeneScan -s %s -o %s/FindORFS/out/%s.orfs -w 0 -t complete"%(_settings.FRAGGENESCAN,contigs,_settings.rundir,_settings.PREFIX), "FindORFS")
       
       parse_fraggenescanout("%s/FindORFS/out/%s.orfs"%(_settings.rundir,_settings.PREFIX),0, "FindORFS", min_len, min_cvg)
       run_process(_settings,"mv %s/FindORFS/out/%s.orfs.ffn %s/FindORFS/out/%s"%(_settings.rundir,_settings.PREFIX,_settings.rundir,outputFNA), "FindORFS")
       run_process(_settings,"mv %s/FindORFS/out/%s.orfs.faa %s/FindORFS/out/%s"%(_settings.rundir,_settings.PREFIX,_settings.rundir,outputFAA), "FindORFS")
       run_process(_settings,"mv %s/FindORFS/out/%s.gene.cvg %s/FindORFS/out/%s"%(_settings.rundir, _settings.PREFIX, _settings.rundir, outputCVG), "FindORFS")
       run_process(_settings,"mv %s/FindORFS/out/%s.gene.map %s/FindORFS/out/%s"%(_settings.rundir, _settings.PREFIX, _settings.rundir, outputMAP), "FindORFS")
   elif orf == "prokka":
      if not os.path.exists(_settings.PROKKA + os.sep + "prokka"):
         print "Error: Prokka not found in %s. Please check your paths and try again"%(_settings.PROKKA)
         raise(JobSignalledBreak)

      libUpdate = "%s/lib/"%(os.path.dirname(_settings.PROKKA))
      if "PERL5LIB" in os.environ:
         libUpdate = "%s%s%s"%(os.environ["PERL5LIB"], os.pathsep, libUpdate)
      os.environ["PERL5LIB"]=libUpdate
      prokkaOptions = getProgramParams(_settings.METAMOS_UTILS, "prokka.spec", "", "-")
      if "--gram" in prokkaOptions and not os.path.exists(_settings.SIGNALP + os.sep + "signalp"):
         print "Warning: Prokka option --gram requires SignalP which is not found. Disabling"
         prokkaOptions = prokkaOptions.replace("--gram", "")

      run_process(_settings, "rm -rf %s/FindORFS/out/%s.prokka"%(_settings.rundir, _settings.PREFIX))
      run_process(_settings, "%s/prokka %s --outdir %s/FindORFS/out/%s.prokka --prefix %s --force %s"%(_settings.PROKKA,"--fast" if _run_fast else "",  _settings.rundir, _settings.PREFIX, _settings.PREFIX, contigs), "FindORFS")
      parse_prokka("%s/FindORFS/out/prokka/%s.gff"%(_settings.rundir,_settings.PREFIX), 0, "FindORFS", min_len, min_cvg)
      run_process(_settings, "ln %s/FindORFS/out/%s.prokka/%s.ffn %s/FindORFS/out/%s"%(_settings.rundir, _settings.PREFIX, _settings.PREFIX, _settings.rundir, outputFNA), "FindORFS")
      run_process(_settings, "ln %s/FindORFS/out/%s.prokka/%s.faa %s/FindORFS/out/%s"%(_settings.rundir, _settings.PREFIX, _settings.PREFIX, _settings.rundir, outputFAA), "FindORFS")
      run_process(_settings,"mv %s/FindORFS/out/%s.gene.cvg %s/FindORFS/out/%s"%(_settings.rundir, _settings.PREFIX, _settings.rundir, outputCVG), "FindORFS")
      run_process(_settings,"mv %s/FindORFS/out/%s.gene.map %s/FindORFS/out/%s"%(_settings.rundir, _settings.PREFIX, _settings.rundir, outputMAP), "FindORFS")

   else:
       #not recognized
       return 1

# warning: this is not thread safe so cannot be run in parallel
@posttask(touch_file("%s/Logs/findorfs.ok"%(_settings.rundir)))
@transform(SplitForORFs, suffix(".contig.cvg"), ".faa")
def FindORFS(input,output):
   originalPrefix = _settings.PREFIX
   _settings.PREFIX = output.replace("%s/Assemble/out/"%(_settings.rundir), "")
   _settings.PREFIX = _settings.PREFIX.replace(".faa", "")

   if "FindORFS" in _skipsteps:
      run_process(_settings, "touch %s/Logs/findorfs.skip"%(_settings.rundir), "FindORFS")
      run_process(_settings, "touch %s/FindORFS/out/%s.faa"%(_settings.rundir, _settings.PREFIX), "FindORFS")
      run_process(_settings, "touch %s/FindORFS/out/%s.fna"%(_settings.rundir, _settings.PREFIX), "FindORFS")
      run_process(_settings, "touch %s/Assemble/out/%s.faa"%(_settings.rundir, _settings.PREFIX), "FindORFS")
      run_process(_settings, "touch %s/Assemble/out/%s.fna"%(_settings.rundir, _settings.PREFIX), "FindORFS")
      run_process(_settings, "touch %s/FindORFS/out/%s.gene.cvg"%(_settings.rundir, _settings.PREFIX), "FindORFS")
      run_process(_settings, "touch %s/FindORFS/out/%s.gene.map"%(_settings.rundir, _settings.PREFIX), "FindORFS")
      _settings.PREFIX = originalPrefix
      return 0

   #if _asm == "soapdenovo":
       #if not os.path.exists("%s/Assemble/out/%s.asm.scafSeq.contigs"%(_settings.rundir,_settings.PREFIX)):
       #    run_process(_settings, "python %s/python/extract_soap_contigs.py %s/Assemble/out/%s.asm.scafSeq"%(_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX))
       #run_process(_settings, "unlink %s/FindORFS/in/%s.asm.scafSeq.contigs"%(_settings.rundir,_settings.PREFIX))
       #run_process(_settings, "unlink %s/FindORFS/in/%s.asm.contig"%(_settings.rundir,_settings.PREFIX))
       #run_process(_settings, "ln -t %s/FindORFS/in/ -s %s/Assemble/out/%s.asm.scafSeq.contigs"%(_settings.rundir, _settings.rundir,_settings.PREFIX))
       #run_process(_settings, "cp %s/FindORFS/in/%s.asm.scafSeq.contigs  %s/FindORFS/in/%s.asm.contig"%(_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX))
       #try using contigs instead of contigs extracted from scaffolds
       #run_process(_settings, "cp %s/Assemble/out/%s.asm.contig  %s/FindORFS/in/%s.asm.contig"%(_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"FindORFS")
   #else:
   run_process(_settings, "unlink %s/FindORFS/in/%s.asm.contig"%(_settings.rundir,_settings.PREFIX),"FindORFS")
   run_process(_settings, "ln %s/Assemble/out/%s.asm.contig %s/FindORFS/in/%s.asm.contig"%(_settings.rundir,_settings.PREFIX,_settings.rundir, _settings.PREFIX),"FindORFS")

   findFastaORFs(_orf, "%s/FindORFS/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX), "%s.ctg.fna"%(_settings.PREFIX), "%s.ctg.faa"%(_settings.PREFIX), "%s.ctg.gene.cvg"%(_settings.PREFIX), "%s.ctg.gene.map"%(_settings.PREFIX), _min_ctg_len, _min_ctg_cvg)
   
   #don't call ORFs on unassembled reads? flag?
   if _read_orfs:
       for lib in _readlibs:
           run_process(_settings, "ln %s/Assemble/out/%s.lib%d.unaligned.fasta %s/FindORFS/in/%s.lib%d.unaligned.fasta"%(_settings.rundir,_settings.PREFIX,lib.id,_settings.rundir, _settings.PREFIX, lib.id),"FindORFS")
           findFastaORFs(_orf, "%s/FindORFS/in/%s.lib%d.unaligned.fasta"%(_settings.rundir, _settings.PREFIX, lib.id), "%s.lib%d.fna"%(_settings.PREFIX, lib.id), "%s.lib%d.faa"%(_settings.PREFIX, lib.id), "%s.lib%d.gene.cvg"%(_settings.PREFIX, lib.id), "%s.lib%d.gene.map"%(_settings.PREFIX, lib.id), 0, 1)

   # merge results
   run_process(_settings, "rm -r %s/FindORFS/out/%s.fna"%(_settings.rundir, _settings.PREFIX), "FindORFS")
   run_process(_settings, "cat %s/FindORFS/out/%s*.fna > %s/FindORFS/out/%s.fna"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "FindORFS")
   run_process(_settings, "rm -r %s/FindORFS/out/%s.fna.bnk"%(_settings.rundir, _settings.PREFIX), "FindORFS")
   run_process(_settings, "%s/toAmos_new -s %s/FindORFS/out/%s.fna -b %s/FindORFS/out/%s.fna.bnk"%(_settings.AMOS, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "FindORFS")

   run_process(_settings, "rm -r %s/FindORFS/out/%s.faa"%(_settings.rundir, _settings.PREFIX), "FindORFS")
   run_process(_settings, "cat %s/FindORFS/out/%s*.faa > %s/FindORFS/out/%s.faa"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "FindORFS")
   run_process(_settings, "rm -r %s/FindORFS/out/%s.faa.bnk"%(_settings.rundir, _settings.PREFIX), "FindORFS")
   run_process(_settings, "%s/toAmos_new -s %s/FindORFS/out/%s.faa -b %s/FindORFS/out/%s.faa.bnk"%(_settings.AMOS, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "FindORFS")

   run_process(_settings, "rm -r %s/FindORFS/out/%s.gene.cvg"%(_settings.rundir, _settings.PREFIX), "FindORFS")
   run_process(_settings, "cat %s/FindORFS/out/%s*.gene.cvg > %s/FindORFS/out/%s.gene.cvg"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "FindORFS")

   run_process(_settings, "rm -r %s/FindORFS/out/%s.gene.map"%(_settings.rundir, _settings.PREFIX), "FindORFS")
   run_process(_settings, "cat %s/FindORFS/out/%s*.gene.map > %s/FindORFS/out/%s.gene.map"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "FindORFS")

   run_process(_settings, "ln %s/FindORFS/out/%s.faa %s/Assemble/out/%s.faa"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "FindORFS")
   run_process(_settings, "ln %s/FindORFS/out/%s.fna %s/Assemble/out/%s.fna"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "FindORFS")

   _settings.PREFIX = originalPrefix

