#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from mapreads import MapReads 
sys.path.append(INITIAL_UTILS)
from ruffus import *

import generic

CGAL_NUM_ALIGN = 100

_readlibs = []
_skipsteps = []
_settings = Settings()
_scores = None
_validators = []

def init(reads, skipsteps, validators, scoreType):
   global _readlibs
   global _skipsteps
   global _scores
   global _validators

   _readlibs = reads
   _skipsteps = skipsteps
   _scores = "%s"%(scoreType)
   _validators = set(validators.strip().split(","))

def minScore():
   return -1 * (sys.maxint - 1)


def getAsmName(input_file_name):
   asmName = input_file_name.replace("%s/Assemble/out/"%(_settings.rundir), "")
   return asmName.replace(".contig.cvg", "").replace(".asm.contig", "")

def getAssemblerName(assembler):
   (asmName, citation) = getProgramCitations(_settings, assembler.lower())
   if citation == "UNKNOWN":
      asmName = assembler.upper()

   return asmName

def getBAMMapped(inputsam):
   numMapped = int(getCommandOutput("%s/samtools flagstat %s |awk '{if (match($4, \"mapped\")) {print $1+$3; }}'"%(_settings.SAMTOOLS, inputsam), False))
   return numMapped

def convertSamToBAM(inputsam):
   run_process(_settings, "%s/samtools view -bS %s > %s.bam"%(_settings.SAMTOOLS, inputsam, inputsam), "Validate")
   run_process(_settings, "%s/samtools sort %s.bam %s.sorted"%(_settings.SAMTOOLS, inputsam, inputsam), "Validate")
   run_process(_settings, "%s/samtools index %s.sorted.bam"%(_settings.SAMTOOLS, inputsam), "Validate")
   return getBAMMapped("%s.sorted.bam"%(inputsam))

def runQUAST(asmNames, asmFiles, min, max, genomeSize, reference):
   if not os.path.exists("%s/metaquast.py"%(_settings.QUAST)):
      return None
   if not os.path.exists(os.path.abspath(reference)):
      return minScore()

   # quast cannot handle dots in names of the files
   asmNames = asmNames.replace(".", "_")
   run_process(_settings, "cat %s |awk -F \"|\" '{if (match($1, \">\")) { if (NF >= 2) { line=$1\"_\"$2; for (i=3; i<= NF; i++) { line=line\"\t\"$i; } } else { line = $1; } print line; } else {print $0; }}' > %s/Validate/out/%s.ref.fastas"%(os.path.baspath(reference), _settings.rundir, _settings.PREFIX), "Validate")
   #run_process(_settings, "ln %s %s/Validate/out/%s.ref.fasta"%(os.path.abspath(reference), _settings.rundir, _settings.PREFIX), "Validate")
   run_process(_settings, "%s/metaquast.py --est-ref-size %s -o %s/Validate/out/quast -R %s/Validate/out/%s.ref.fasta -T %d -l \"%s\" %s"%(_settings.QUAST, genomeSize, _settings.rundir, _settings.rundir, _settings.PREFIX, _settings.threads, asmNames, asmFiles), "Validate")

   return 0

def runSNP(inputsam, prefix, assembly, min, max, genomeSize):
   if not os.path.exists("%s/freebayes"%(_settings.FREEBAYES)) or not os.path.exists("%s.sorted.bam"%(inputsam)):
      return None

   if getBAMMapped("%s.sorted.bam"%(inputsam)) == 0:
      return minScore()

   freebayesOptions = getProgramParams(_settings.METAMOS_UTILS, "freebayes.spec", "", "-")

   run_process(_settings, "%s/freebayes %s -E 0 -X -u -p 1 -b %s.sorted.bam -v %s/Validate/out/%s.vcf -f %s"%(_settings.FREEBAYES, freebayesOptions, inputsam, _settings.rundir, prefix, assembly), "Validate")
   num_snps = getCommandOutput("cat %s/Validate/out/%s.vcf |grep -v \"#\" |wc -l"%(_settings.rundir, prefix), False)
   return -1 * int(num_snps)

def runLAP(prefix, assembly, pairedReads, unpairedReads, abundanceFile = ""):
   if abundanceFile != "":
      abundanceFile = "-n %s"%(abundanceFile)

   run_process(_settings, "python %s/aligner/calc_prob.py -k --output_sam_file %s.sam -p %d -a %s %s %s > %s/Validate/out/%s.prob"%(_settings.LAP, prefix, _settings.threads, assembly, abundanceFile, pairedReads if pairedReads != "" else unpairedReads, _settings.rundir, prefix), "Validate")
   lapScore = getCommandOutput("python %s/aligner/sum_prob.py -i %s/Validate/out/%s.prob"%(_settings.LAP, _settings.rundir, prefix), True).split()[0]
   return lapScore

def runFRCBAM(inputsam, prefix, assembly, min, max, genomeSize):
   if not os.path.exists("%s/FRC"%(_settings.FRCBAM)) or not os.path.exists("%s.sorted.bam"%(inputsam)):
      return None

   if getBAMMapped("%s.sorted.bam"%(inputsam)) == 0:
      return minScore()

   if "LD_LIBRARY_PATH" in os.environ:
      os.environ["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"] + os.pathsep + _settings.FRCBAM
   else:
      os.environ["LD_LIBRARY_PATH"] = _settings.FRCBAM
   run_process(_settings, "%s/FRC --pe-sam %s.sorted.bam --pe-min-insert %d --pe-max-insert %d --genome-size %s --output frc_%s"%(_settings.FRCBAM, inputsam, min, max, genomeSize, prefix), "Validate")
   num_features = getCommandOutput("cat %s/Validate/out/frc_%s_FRC.txt |tail -n 1 |awk '{print $1}'"%(_settings.rundir, prefix), False)

   return -1 * int(num_features)

def runALE(inputsam, prefix, assembly, min, max, genomeSize):
   if not os.path.exists("%s/ALE"%(_settings.ALE)):
      return None

   run_process(_settings, "%s/ALE %s %s %s.ale"%(_settings.ALE, inputsam, assembly, prefix), "Validate")
   ale_score = getCommandOutput("cat %s/Validate/out/%s.ale |head -n 1 |awk '{print $NF}'"%(_settings.rundir, prefix), False)
   return ale_score 

def runCGAL(inputsam, prefix, assembly, min, max, genomeSize):
   if not os.path.exists("%s/cgal"%(_settings.CGAL)):
      return None

   if getBAMMapped("%s.sorted.bam"%(inputsam)) == 0:
      return minScore()

   run_process(_settings, "%s/bowtie2convert %s"%(_settings.CGAL, inputsam), "Validate")
   run_process(_settings, "%s/align %s %d %d"%(_settings.CGAL, assembly, CGAL_NUM_ALIGN, _settings.threads), "Validate")
   run_process(_settings, "%s/cgal %s > %s/Validate/out/%s.cgal"%(_settings.CGAL, assembly, _settings.rundir, prefix), "Validate")
   cgal_score = getCommandOutput("cat %s/Validate/out/%s.cgal"%(_settings.rundir, prefix), False)

   return cgal_score

@posttask(touch_file("%s/Logs/validate.ok"%(_settings.rundir)))
@merge(MapReads, ["%s/Logs/validate.ok"%(_settings.rundir)])
def Validate (input_file_names, output_file_name):
   bestScores = dict()
   bestAssemblers = dict()
   bestAssemblies = dict()
   bestAssembler = ""
   bestAssembly = ""

   selectedAsm = open("%s/Validate/out/%s.asm.selected"%(_settings.rundir, _settings.PREFIX), 'w')
   selectedReferences = open("%s/Validate/out/%s.ref.selected"%(_settings.rundir, _settings.PREFIX), 'w')
   lapfile = open("%s/Validate/out/%s.lap"%(_settings.rundir,_settings.PREFIX),'w')

   failedOutput = ""

   if "Validate" in _skipsteps or "validate" in _skipsteps:
      run_process(_settings, "touch %s/Logs/validate.skip"%(_settings.rundir), "Validate")
      bestAssembler = getAsmName(input_file_names[0])
      bestAssembly = input_file_names[0].replace(".contig.cvg", ".asm.contig") 
   elif len(input_file_names) == 1:
      run_process(_settings, "touch %s/Logs/validate.skip"%(_settings.rundir), "Validate")
      bestAssembler = getAsmName(input_file_names[0])
      bestAssembly = input_file_names[0].replace(".contig.cvg", ".asm.contig") 
   elif not os.path.exists("%s/bowtie2"%(_settings.BOWTIE2)) or not os.path.exists("%s/aligner/calc_prob.py"%(_settings.LAP)):
      run_process(_settings, "touch %s/Logs/validate.skip"%(_settings.rundir), "Validate")
      print "Warning! LAP is not available, cannot select best assembly, chosing first available: %s!"%(getAsmName(input_file_names[0]))
      bestAssembler = getAsmName(input_file_names[0])
      bestAssembly = input_file_names[0].replace(".contig.cvg", ".asm.contig") 
   else:
      # build string of files to use for validation
      os.environ["BT2_HOME"]=_settings.BOWTIE2

      pairedReads = ""
      pairedMin = 0
      pairedMax = 0
      unpairedReads = ""
      for lib in _readlibs:
         if lib.mated and pairedReads == "":
            pairedReads="-1 %s/Preprocess/out/lib%d.1.fastq -2 %s/Preprocess/out/lib%d.2.fastq -m %d -t %d -I %d -X %d"%(_settings.rundir, lib.id, _settings.rundir, lib.id, lib.mean, lib.stdev, (lib.mean-5*lib.stdev), (lib.mean+5*lib.stdev))
            pairedMin = lib.mmin
            pairedMax = lib.mmax
         elif not lib.paired and unpairedReads == "":
            unpaired=" -i %s/Preprocess/out/lib%d.fastq"%(_settings.rundir, lib.id)

      asmNames = ""
      asmFiles = ""
      firstLine = True
      for input_file_name in input_file_names:
         assembler = getAsmName(input_file_name)
         assembly = input_file_name.replace(".contig.cvg", ".asm.contig")
         abundanceFile = ""
         if os.path.exists("%s/Assemble/out/%s.contig.cvg"%(_settings.rundir, assembler)):
            abundanceFile = "%s/Assemble/out/%s.contig.cvg"%(_settings.rundir, assembler)
         scores = dict()
         scores[SCORE_TYPE.LAP] = runLAP(assembler, assembly, pairedReads, unpairedReads, abundanceFile)

         scoreOutput = ""
         inputSam = "%s/Validate/out/%s.sam_0"%(_settings.rundir, assembler)
         numMapped = 0
         if os.path.exists("%s/samtools"%(_settings.SAMTOOLS)):
            numMapped = convertSamToBAM(inputSam)
         genomeSize = 0

         if numMapped != 0:
            if asmNames == "":
               asmNames = assembler
               asmFiles = assembly
            else:
               asmNames = "%s,%s"%(asmNames, assembler)
               asmFiles = "%s %s"%(asmFiles, assembly)

         for validator in _validators:
            if validator.lower() == "lap" or validator.lower() == "quast":
               continue
            elif validator.lower() == "frcbam":
               if pairedReads != "":
                  scores[SCORE_TYPE.FRCBAM] = runFRCBAM(inputSam, assembler, assembly, pairedMin, pairedMax, genomeSize)
               else:
		  scores[SCORE_TYPE.FRCBAM] = scores[SCORE_TYPE.LAP]
            elif validator.lower() == "ale":
               scores[SCORE_TYPE.ALE] = runALE(inputSam, assembler, assembly, pairedMin, pairedMax, genomeSize)
            elif validator.lower() == "snp" or validator.lower() == "freebayes":
               scores[SCORE_TYPE.SNP] = runSNP(inputSam, assembler, assembly, pairedMin, pairedMax, genomeSize)
            elif validator.lower() == "cgal":
               scores[SCORE_TYPE.CGAL] = runCGAL(inputSam, assembler, assembly, pairedMin, pairedMax, genomeSize)
            else:
               print "Unknown validator %s requested, skipping"%(validator.lower())

         debugString = "*** metAMOS assembler %s scores: "%(assembler)
         for type in scores.keys():
            outputScore = scores[type]
            if scores[type] != None and scores[type] == minScore():
               outputScore = "None"
            elif scores[type] != None and (type == SCORE_TYPE.FRCBAM or type == SCORE_TYPE.SNP):
               outputScore = -1 * scores[type]

            if scoreOutput == "":
               scoreOutput = "%s"%(outputScore)
               failedOutput = "None"
               header = "%s"%(SCORE_TYPE.reverse_mapping[type])
            else:
               scoreOutput = "%s\t%s"%(scoreOutput, outputScore)
               failedOutput = "%s\tNone"%(failedOutput)
               header = "%s\t%s"%(header, SCORE_TYPE.reverse_mapping[type])

            if scores[type] == None:
               continue

            debugString = "%s %s:%s"%(debugString, SCORE_TYPE.reverse_mapping[type], scores[type])
            if type not in bestScores or float(scores[type]) > bestScores[type]:
               bestScores[type] = float(scores[type])
               bestAssemblers[type] = assembler
               bestAssemblies[type] = assembly
         if _settings.VERBOSE:
            print "%s"%(debugString)

         if firstLine:
            lapfile.write("Assembly\t%s\n"%(header))
            firstLine = False
         asmName = getAssemblerName(assembler)
         lapfile.write("%s\t%s\n"%(asmName, scoreOutput))
         lapfile.flush()

   # output failed assemblers
   for file in os.listdir("%s/Assemble/out/"%(_settings.rundir)): 
      if (file.endswith(".failed")):
          assembler = os.path.splitext(os.path.basename(file))[0]
          asmName = getAssemblerName(assembler)
          lapfile.write("%s\t%s\n"%(asmName, failedOutput))
          lapfile.flush()

   # select assembly
   global _scores
   if bestAssembler == "": 
      if "%s"%(SCORE_TYPE.ALL) in _scores:
         _scores = SCORE_TYPE.reverse_mapping.keys()

      assemblerVotes = dict()
      assemblyVotes = dict()
      currMaxVote = 0
      processedType = dict()

      for type in _scores:
         if type == "":
            continue

         type = int(type)
         if type == SCORE_TYPE.ALL:
            continue
         if type not in bestAssemblers.keys():
            type = SCORE_TYPE.LAP
            print "Warning: selected score %s was not available, defaulting to using LAP"%(SCORE_TYPE.reverse_mapping[type].upper())
         if type in processedType.keys():
            continue
         processedType[type] = True

         if bestAssemblers[type] not in assemblerVotes:
            assemblerVotes[bestAssemblers[type]] = 0
            assemblyVotes[bestAssemblies[type]] = 0
         weight = 1
         try:
            weight = SCORE_WEIGHTS[type]
         except KeyError:
            print "Uknown Type %s using weight %d"%(SCORE_TYPE.reverse_mapping[type], weight)
            weight = 1
         assemblerVotes[bestAssemblers[type]] += weight
         assemblyVotes[bestAssemblies[type]] += weight
         print "Votes for assembler %s is %s (type %s just voted)"%(bestAssemblers[type], assemblerVotes[bestAssemblers[type]], SCORE_TYPE.reverse_mapping[type])
     
      for assembler in assemblerVotes:
         if assemblerVotes[assembler] > currMaxVote:
            bestAssembler = assembler
            currMaxVote = assemblerVotes[assembler]

      currMaxVote = 0
      for assembly in assemblyVotes:
         if assemblyVotes[assembly] > currMaxVote:
            bestAssembly = assembly
            currMaxVote = assemblyVotes[assembly]
      if getAsmName(bestAssembly) != bestAssembler:
         print "Error: inconsistent assembly and assembler chosen %s %s"%(bestAssembly, bestAssembler)
         raise(JobSignalledBreak)

   if "quast" in _validators:
      # recruit a reference
      references = recruitGenomes(_settings,bestAssembly,"%s/refseq"%(_settings.DB_DIR), "%s/Validate/out/recruit"%(_settings.rundir), "Validate", 1)
      if len(references) > 0:
         for g in references:
            selectedReferences.write(os.path.splitext(os.path.basename(g))[0] + "\n")
         runQUAST(asmNames, asmFiles, pairedMin, pairedMax, genomeSize, references[0])
      else:
         print "Warning: could not recruit any references"

   if _settings.VERBOSE:
      print "*** metAMOS assembler %s selected."%(bestAssembler)  
   run_process(_settings, "ln %s %s/Assemble/out/%s.asm.contig"%(bestAssembly, _settings.rundir, _settings.PREFIX), "Validate")
   if os.path.exists("%s/Assemble/out/%s.linearize.scaffolds.final"%(_settings.rundir, bestAssembler)):
      run_process(_settings, "ln %s/Assemble/out/%s.linearize.scaffolds.final %s/Assemble/out/%s.linearize.scaffolds.final"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")
   else:
      run_process(_settings, "ln %s %s/Assemble/out/%s.linearize.scaffolds.final"%(bestAssembly, _settings.rundir, _settings.PREFIX), "Validate")
   run_process(_settings, "ln %s/Assemble/out/%s.asm.tigr %s/Assemble/out/%s.asm.tigr"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")
   run_process(_settings, "ln %s/Assemble/out/%s.contig.cnt %s/Assemble/out/%s.contig.cnt"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")
   run_process(_settings, "ln %s/Assemble/out/%s.contig.cvg %s/Assemble/out/%s.contig.cvg"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")
   if os.path.exists("%s/Assemble/out/%s.afg"%(_settings.rundir, bestAssembler)):
      run_process(_settings, "ln %s/Assemble/out/%s.afg %s/Assemble/out/%s.afg"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")

   for lib in _readlibs: 
      run_process(_settings, "ln %s/Assemble/out/%s.lib%d.badmades %s/Assemble/out/%s.lib%d.badmates"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.lib%d.hdr %s/Assemble/out/%s.lib%d.hdr"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.lib%d.mappedmates %s/Assemble/out/%s.lib%d.mappedmates"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.lib%d.mates_in_diff_contigs %s/Assemble/out/%s.lib%d.mates_in_diff_contigs"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.lib%dcontig.reads %s/Assemble/out/%s.lib%dcontig.reads"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.lib%d.unaligned.fasta %s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, lib.id), "Validate")
      if os.path.exists("%s/Assemble/out/%s.lib%d.unaligned.fastq"%(_settings.rundir, bestAssembler, lib.id)):
         run_process(_settings, "ln %s/Assemble/out/%s.lib%d.unaligned.fastq %s/Assemble/out/lib%d.unaligned.fastq"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, lib.id), "Validate")

   selectedAsm.write("%s"%(bestAssembler))
   selectedAsm.close()
   selectedReferences.close()
   lapfile.close()

