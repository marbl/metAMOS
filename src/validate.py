#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from findorfs import FindORFS,setRunFast
from mapreads import getMeanSD
sys.path.append(INITIAL_UTILS)
from ruffus import *

import generic

CGAL_NUM_ALIGN = 500
TOP_PERCENTILE = 0.10
QUALIFY_RATIO = 2

_readlibs = []
_skipsteps = []
_settings = Settings()
_scores = None
_validators = []
_tieBreakerType = SCORE_TYPE.LAP
_skipTypesForCandidates = [ SCORE_TYPE.ORF ]

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
   return asmName.replace(".faa", "").replace(".asm.contig", "")

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
   setFailFast(False)

   if not os.path.exists("%s/metaquast.py"%(_settings.QUAST)):
      return None
   if not os.path.exists(os.path.abspath(reference)):
      return minScore()

   # quast cannot handle dots in names of the files
   asmNames = asmNames.replace(".", "_")
   run_process(_settings, "cat %s |awk -F \"|\" '{if (match($1, \">\")) { if (NF >= 2) { line=$1\"_\"$2; for (i=3; i<= NF; i++) { line=line\"\t\"$i; } } else { line = $1; } print line; } else {print $0; }}' > %s/Validate/out/%s.ref.fasta"%(reference, _settings.rundir, _settings.PREFIX), "Validate")
   #run_process(_settings, "ln %s %s/Validate/out/%s.ref.fasta"%(os.path.abspath(reference), _settings.rundir, _settings.PREFIX), "Validate")
   run_process(_settings, "%s/metaquast.py --est-ref-size %s -o %s/Validate/out/quast -R %s/Validate/out/%s.ref.fasta -T %d -l \"%s\" %s"%(_settings.QUAST, genomeSize, _settings.rundir, _settings.rundir, _settings.PREFIX, _settings.threads, asmNames, asmFiles), "Validate")

   setFailFast(True)
   return 0

def runSNP(inputsam, prefix, assembly, min, max, genomeSize):
   if not os.path.exists("%s/freebayes"%(_settings.FREEBAYES)) or not os.path.exists("%s.sorted.bam"%(inputsam)):
      return None

   if getBAMMapped("%s.sorted.bam"%(inputsam)) == 0:
      return minScore()

   setFailFast(False)
   freebayesOptions = getProgramParams(_settings.METAMOS_UTILS, "freebayes.spec", "", "-")

   run_process(_settings, "%s/freebayes %s -E 0 -X -u -p 1 -b %s.sorted.bam -v %s/Validate/out/%s.vcf -f %s"%(_settings.FREEBAYES, freebayesOptions, inputsam, _settings.rundir, prefix, assembly), "Validate")
   num_snps = getCommandOutput("cat %s/Validate/out/%s.vcf |grep -v \"#\" |wc -l"%(_settings.rundir, prefix), False)
   setFailFast(True)
   if num_snps == "":
      return minScore()

   return -1 * int(num_snps)

def runLAP(prefix, assembly, pairedReads, unpairedReads, abundanceFile = ""):
   if abundanceFile != "":
      abundanceFile = "-n %s"%(abundanceFile)

   run_process(_settings, "python %s/aligner/calc_prob.py -q -k --output_sam_file %s.sam -p %d -a %s %s %s > %s/Validate/out/%s.prob"%(_settings.LAP, prefix, _settings.threads, assembly, abundanceFile, pairedReads if pairedReads != "" else unpairedReads, _settings.rundir, prefix), "Validate")
   lapScore = getCommandOutput("python %s/aligner/sum_prob.py -i %s/Validate/out/%s.prob"%(_settings.LAP, _settings.rundir, prefix), True).split()[0]

   # paired-end lap produces more than one sam file, while unpaired does not
   # since we always use a single paired-end library for validation, we can link to the pair
   if pairedReads != "":
      run_process(_settings, "ln %s/Validate/out/%s.sam_0 %s/Validate/out/%s.sam"%(_settings.rundir, prefix, _settings.rundir, prefix), "Validate")

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

   setFailFast(False)
   run_process(_settings, "%s/ALE %s %s %s.ale"%(_settings.ALE, inputsam, assembly, prefix), "Validate")
   ale_score = getCommandOutput("cat %s/Validate/out/%s.ale |head -n 1 |awk '{print $NF}'"%(_settings.rundir, prefix), False)
   setFailFast(True)

   if ale_score == "":
      ale_score = minScore()

   return ale_score 

def runCGAL(inputsam, prefix, assembly, min, max, genomeSize):
   if not os.path.exists("%s/cgal"%(_settings.CGAL)):
      return None

   if getBAMMapped("%s.sorted.bam"%(inputsam)) == 0:
      return minScore()

   setFailFast(False)
   run_process(_settings, "%s/bowtie2convert %s"%(_settings.CGAL, inputsam), "Validate")
   run_process(_settings, "%s/align %s %d %d"%(_settings.CGAL, assembly, CGAL_NUM_ALIGN, _settings.threads), "Validate")
   run_process(_settings, "%s/cgal %s > %s/Validate/out/%s.cgal"%(_settings.CGAL, assembly, _settings.rundir, prefix), "Validate")
   cgal_score = getCommandOutput("cat %s/Validate/out/%s.cgal"%(_settings.rundir, prefix), False)
   setFailFast(True)
   if cgal_score == "":
      cgal_score = minScore()

   return cgal_score

def runORF(inputsam, prefix, assembly, min, max, genomeSize):
   orf_score = getCommandOutput("grep -c '>' %s |awk -F \":\" '{print $NF}'"%(assembly.replace(".asm.contig", ".fna")), False)
   return orf_score

def runN50(inputsam, prefix, assembly, min, max, genomeSize):
   n50 = getCommandOutput("java -cp %s GetFastaStats -n50 -genomeSize %s %s |awk '{print $NF}'"%(_settings.METAMOS_JAVA, genomeSize, assembly), False)
   if n50 == "":
      return minScore()
   return n50

def runREAPR(pairedFiles, prefix, assembly, min, max, genomeSize):
   if not os.path.exists("%s/reapr"%(_settings.REAPR)):
      return None

   # reapr has many reasons for failing, including inconsistent line lengths in fasta file or insufficient sequences mapped to estimate insert sizes
   # thus, if it fails, return min score
   setFailFast(False)
   libUpdate = "%s/lib/"%(_settings.REAPR)
   if "PERL5LIB" in os.environ:
      libUpdate = "%s%s%s"%(os.environ["PERL5LIB"], os.pathsep, libUpdate)
   os.environ["PERL5LIB"]=libUpdate
   run_process(_settings, "%s/reapr facheck %s %s/Validate/out/%s.reapr"%(_settings.REAPR, assembly, _settings.rundir, prefix), "Validate")
   # does not support multiple threads so runs slowly even for bacterial genomes, run in parallel manually
   run_process(_settings, "%s/reapr smaltmap -x %s/Validate/out/%s.reapr.fa %s %s/Validate/out/%s.reapr.bam > %s/Validate/out/%s.reapr.cmds"%(_settings.REAPR, _settings.rundir, prefix, pairedFiles, _settings.rundir, prefix, _settings.rundir, prefix), "Validate")
   commandList = open("%s/Validate/out/%s.reapr.cmds"%(_settings.rundir, prefix), 'r') 
   for cmd in commandList.xreadlines():
      if "smalt map" in cmd and _settings.OSTYPE != "Darwin":
         cmd = cmd.replace("smalt map", "smalt map -n %d"%(_settings.threads))
      run_process(_settings, "%s"%(cmd), "Validate")
   commandList.close()

   run_process(_settings, "%s/reapr perfectmap %s/Validate/out/%s.reapr.fa %s %s %s/Validate/out/%s.reapr.perfect"%(_settings.REAPR, _settings.rundir, prefix, pairedFiles, (min+max)/2, _settings.rundir, prefix), "Validate")

   # finally run reapr
   perfectMapCount = getCommandOutput("cat %s/Validate/out/%s.reapr.perfect.hist | awk -v SUM=0 '{if ($1 == 0) { ZERO = $NF; }  SUM+=$NF; } END {print ZERO/SUM*100}'"%(_settings.rundir, prefix), False)

   if perfectMapCount != "" and float(perfectMapCount) < 50:
      run_process(_settings, "%s/reapr pipeline %s/Validate/out/%s.reapr.fa %s/Validate/out/%s.reapr.bam %s/Validate/out/%s.reapr %s/Validate/out/%s.reapr.perfect"%(_settings.REAPR, _settings.rundir, prefix, _settings.rundir, prefix, _settings.rundir, prefix, _settings.rundir, prefix), "Validate")
   else:
      run_process(_settings, "%s/reapr pipeline %s/Validate/out/%s.reapr.fa %s/Validate/out/%s.reapr.bam %s/Validate/out/%s.reapr"%(_settings.REAPR, _settings.rundir, prefix, _settings.rundir, prefix, _settings.rundir, prefix), "Validate")

   setFailFast(True)

   if os.path.exists("%s/Validate/out/%s.reapr/05.summary.report.tsv"%(_settings.rundir, prefix)):
      reapr_score = getCommandOutput("cat %s/Validate/out/%s.reapr/05.summary.report.tsv |grep -v \"#\" | awk -F \"\t\" '{printf(\"%%.2f\",($38/$2)*100)}'"%(_settings.rundir, prefix), False)
   else:
      reapr_score = minScore()
   return reapr_score

@posttask(touch_file("%s/Logs/validate.ok"%(_settings.rundir)))
@merge(FindORFS, ["%s/Logs/validate.ok"%(_settings.rundir)])
def Validate (input_file_names, output_file_name):
   validationScores = dict()

   bestAssembler = ""
   bestAssembly = ""

   tieBreakingScores = dict()

   validatedAsms = dict()
   scoreOrder = dict()
   genomeSize = getEstimatedGenomeSize(_settings)
   needToOutput = False

   if os.path.exists("%s/Validate/out/%s.lap"%(_settings.rundir,_settings.PREFIX)) and os.path.getsize("%s/Validate/out/%s.lap"%(_settings.rundir, _settings.PREFIX)):
      lapfile = open("%s/Validate/out/%s.lap"%(_settings.rundir,_settings.PREFIX),'r')
      for line in lapfile.xreadlines():
         scores = line.split()
         if len(scoreOrder) == 0:
            i = 0
            for score in scores:
               if score.upper() in SCORE_TYPE.mapping:
                  scoreOrder[i] = SCORE_TYPE.mapping[score.upper()]
               i += 1
         else:
            i = 0
            asmScores = dict()
            for score in scores:
               if i in scoreOrder:
                   asmScores[scoreOrder[i]] = minScore() if score.lower() == "none" else float(score)
                   if (scoreOrder[i] == SCORE_TYPE.FRCBAM or scoreOrder[i] == SCORE_TYPE.SNP):
                      asmScores[scoreOrder[i]] = -1 * asmScores[scoreOrder[i]]
               i += 1 
            validatedAsms[scores[0].lower()] = asmScores
      lapfile.close()
      lapfile = open("%s/Validate/out/%s.lap"%(_settings.rundir,_settings.PREFIX),'a')
   else:
      lapfile = open("%s/Validate/out/%s.lap"%(_settings.rundir,_settings.PREFIX),'w')

   failedOutput = ""
   totalRun = 0

   if "Validate" in _skipsteps or "validate" in _skipsteps:
      run_process(_settings, "touch %s/Logs/validate.skip"%(_settings.rundir), "Validate")
      bestAssembler = getAsmName(input_file_names[0])
      bestAssembly = "%s/Assemble/out/%s.asm.contig"%(_settings.rundir, bestAssembler)
   elif len(input_file_names) == 1 and len(_validators) == 1:
      run_process(_settings, "touch %s/Logs/validate.skip"%(_settings.rundir), "Validate")
      bestAssembler = getAsmName(input_file_names[0])
      bestAssembly = "%s/Assemble/out/%s.asm.contig"%(_settings.rundir, bestAssembler)
   elif not os.path.exists("%s/bowtie2"%(_settings.BOWTIE2)) or not os.path.exists("%s/aligner/calc_prob.py"%(_settings.LAP)):
      run_process(_settings, "touch %s/Logs/validate.skip"%(_settings.rundir), "Validate")
      print "Warning! LAP is not available, cannot select best assembly, chosing first available: %s!"%(getAsmName(input_file_names[0]))
      bestAssembler = getAsmName(input_file_names[0])
      bestAssembly = "%s/Assemble/out/%s.asm.contig"%(_settings.rundir, bestAssembler)
   elif os.path.exists("%s/Validate/out/%s.asm.selected"%(_settings.rundir, _settings.PREFIX)):
      selectedAsm = open("%s/Validate/out/%s.asm.selected"%(_settings.rundir, _settings.PREFIX), 'r')
      bestAssembler = selectedAsm.read().strip()
      bestAssembly = "%s/Assemble/out/%s.asm.contig"%(_settings.rundir, bestAssembler)
      selectedAsm.close()
   else:
      # build string of files to use for validation
      os.environ["BT2_HOME"]=_settings.BOWTIE2

      pairedReads = ""
      pairedFiles = ""
      pairedMin = 0
      pairedMax = 0
      unpairedReads = ""
      for lib in _readlibs:
         if lib.mated and pairedReads == "" and lib.innie:
            (pairedMin, pairedMax, pairedMean, pairedSD) = getMeanSD(lib.id) 
            pairedReads="-1 %s/Preprocess/out/lib%d.1.fastq -2 %s/Preprocess/out/lib%d.2.fastq -m %d -t %d -I %d -X %d"%(_settings.rundir, lib.id, _settings.rundir, lib.id, pairedMean, pairedSD, 0 if pairedMean < 5*pairedSD else (pairedMean-5*pairedSD), (pairedMean+5*pairedSD))
            pairedFiles="%s/Preprocess/out/lib%d.1.fastq %s/Preprocess/out/lib%d.2.fastq"%(_settings.rundir, lib.id, _settings.rundir, lib.id)
         elif not lib.mated and unpairedReads == "":
            unpairedReads = "-i %s/Preprocess/out/lib%d.fastq"%(_settings.rundir, lib.id)

      asmNames = ""
      asmFiles = ""
      firstLine = True

      for input_file_name in input_file_names:
         assembler = getAsmName(input_file_name)
         assembly = input_file_name.replace(".faa", ".asm.contig")
         abundanceFile = ""
         if os.path.exists("%s/Assemble/out/%s.contig.cvg"%(_settings.rundir, assembler)) and os.path.getsize("%s/Assemble/out/%s.contig.cvg"%(_settings.rundir, assembler)) > 0:
            abundanceFile = "%s/Assemble/out/%s.contig.cvg"%(_settings.rundir, assembler)
         scores = dict()

         if assembler == _settings.PREFIX:
            continue

         scoreOutput = ""
         needToRun = True
         # if we already had a score for an assembler,we dont need to revalidate it
         if assembler.lower() in validatedAsms:
            asmScores = validatedAsms[assembler.lower()]
            needToRun = False
            firstLine = False
            missingScores = ""

            inputSam = "%s/Validate/out/%s.sam.sorted.bam"%(_settings.rundir, assembler)
            numMapped = 0
            if os.path.exists("%s"%(inputSam)):
               numMapped = getBAMMapped(inputSam)
               if numMapped != 0:
                  if asmNames == "":
                     asmNames = assembler
                     asmFiles = assembly
                  else:
                     asmNames = "%s,%s"%(asmNames, assembler)
                     asmFiles = "%s %s"%(asmFiles, assembly)
            else:
               needToRun = True

            # check for missing scores, if the validator list was changed on the restart, we will need to rerun 
            #    or if any failed for a given assembler
            # this would be better if it did not duplicate code below for running validators
            for validator in _validators:
               if validator.lower() == "quast":
                  continue
               elif validator.lower() == "lap":
                  if SCORE_TYPE.LAP not in asmScores:
                     needToRun = True
                     missingScores = missingScores + " %s"%(validator.lower())
                  else:
                     scores[SCORE_TYPE.LAP] = asmScores[SCORE_TYPE.LAP]
               elif validator.lower() == "frcbam":
                  if SCORE_TYPE.FRCBAM not in asmScores:
                     needToRun = True
                     missingScores = missingScores + " %s"%(validator.lower())
                  else:
                     scores[SCORE_TYPE.FRCBAM] = asmScores[SCORE_TYPE.FRCBAM]
               elif validator.lower() == "ale":
                  if SCORE_TYPE.ALE not in asmScores:
                     needToRun = True
                     missingScores = missingScores + " %s"%(validator.lower())
                  else:
                     scores[SCORE_TYPE.ALE] = asmScores[SCORE_TYPE.ALE]
               elif validator.lower() == "snp" or validator.lower() == "freebayes":
                  if SCORE_TYPE.SNP not in asmScores:
                     needToRun = True
                     missingScores = missingScores + " %s"%(validator.lower())
                  else:
                     scores[SCORE_TYPE.SNP] = asmScores[SCORE_TYPE.SNP]
               elif validator.lower() == "cgal":
                  if SCORE_TYPE.CGAL not in asmScores:
                     needToRun = True
                     missingScores = missingScores + " %s"%(validator.lower())
                  else:
                     scores[SCORE_TYPE.CGAL] = asmScores[SCORE_TYPE.CGAL]
               elif validator.lower() == "orf":
                  if SCORE_TYPE.ORF not in asmScores:
                     needToRun = True
                     missingScores = missingScores + " %s"%(validator.lower())
                  else:
                     scores[SCORE_TYPE.ORF] = asmScores[SCORE_TYPE.ORF]
               elif validator.lower() == "n50":
                  if SCORE_TYPE.N50 not in asmScores:
                     needToRun = True
                     missingScores = missingScores + " %s"%(validator.lower())
                  else:
                     scores[SCORE_TYPE.N50] = asmScores[SCORE_TYPE.N50]
               elif validator.lower() == "reapr":
                  if SCORE_TYPE.REAPR not in asmScores:
                     needToRun = True
                     missingScores = missingScores + " %s"%(validator.lower())
                  else:
                     scores[SCORE_TYPE.REAPR] = asmScores[SCORE_TYPE.REAPR]

            for type in asmScores:
               if scoreOutput == "":
                  scoreOutput = "%s"%(asmScores[type])
               else:
                  scoreOutput = "%s\t%s"%(scoreOutput, asmScores[type])

            if needToRun == True:
               # as soon as one assembler fails, we truncate the file so we need to output all those that could be resumed as well
               needToOutput = True
               scoreOutput = ""

               del validatedAsms[assembler.lower()]
               print "*** metAMOS Warning: validation resume for assembler %s not possible, missing scores %s"%(assembler.upper(), missingScores.strip())
               asmNames = ""
               firstLine = True
               lapfile.seek(0)
               lapfile.truncate()

         header = ""
         if needToRun:
            totalRun += 1
            scores[SCORE_TYPE.LAP] = runLAP(assembler, assembly, pairedReads, unpairedReads, abundanceFile)

            inputSam = "%s/Validate/out/%s.sam"%(_settings.rundir, assembler)
            numMapped = 0
            if os.path.exists("%s/samtools"%(_settings.SAMTOOLS)):
               numMapped = convertSamToBAM(inputSam)

            if numMapped != 0:
               if asmNames == "":
                  asmNames = assembler
                  asmFiles = assembly
               else:
                  asmNames = "%s,%s"%(asmNames, assembler)
                  asmFiles = "%s %s"%(asmFiles, assembly)

            # run each selected validator, lap always run so skip it. quast runs at the end to compare asms
            for validator in _validators:
               if validator.lower() == "lap" or validator.lower() == "quast":
                  continue
               elif validator.lower() == "frcbam":
                  if pairedReads != "":
                     scores[SCORE_TYPE.FRCBAM] = runFRCBAM(inputSam, assembler, assembly, pairedMin, pairedMax, genomeSize)
                  else:
	   	     scores[SCORE_TYPE.FRCBAM] = minScore()
               elif validator.lower() == "ale":
                  scores[SCORE_TYPE.ALE] = runALE(inputSam, assembler, assembly, pairedMin, pairedMax, genomeSize)
               elif validator.lower() == "snp" or validator.lower() == "freebayes":
                  scores[SCORE_TYPE.SNP] = runSNP(inputSam, assembler, assembly, pairedMin, pairedMax, genomeSize)
               elif validator.lower() == "cgal":
                  scores[SCORE_TYPE.CGAL] = runCGAL(inputSam, assembler, assembly, pairedMin, pairedMax, genomeSize)
               elif validator.lower() == 'orf':
                  scores[SCORE_TYPE.ORF] = runORF(inputSam, assembler, assembly, pairedMin, pairedMax, genomeSize)
               elif validator.lower() == 'n50':
                  scores[SCORE_TYPE.N50] = runN50(inputSam, assembler, assembly, pairedMin, pairedMax, genomeSize)
               elif validator.lower() == 'reapr':
                  if pairedReads != "":
                     scores[SCORE_TYPE.REAPR] = runREAPR(pairedFiles, assembler, assembly, pairedMin, pairedMax, genomeSize)
                  else:
                     scores[SCORE_TYPE.REAPR] = minScore()
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

            # if a validator failed to run, skip it and dont record scores
            if scores[type] == None or scores[type] == minScore():
               continue

            # record for later if this is the tiebreaker
            if type == _tieBreakerType:
               tieBreakingScores[assembler] = float(scores[type])

            debugString = "%s %s:%s"%(debugString, SCORE_TYPE.reverse_mapping[type], scores[type])

            # record in db of all scores
            if type not in validationScores:
               validationScores[type] = dict()
            validationScores[type][assembler] = float(scores[type])

         if _settings.VERBOSE:
            print "%s"%(debugString)

         if needToRun: 
            if firstLine:
               lapfile.write("Assembly\t%s\n"%(header))
               firstLine = False
            asmName = getAssemblerName(assembler)
            lapfile.write("%s\t%s\n"%(asmName, scoreOutput))
            lapfile.flush()

   # output list of previously validated assemblers that could be resumed if needed
   if needToOutput:
      for assembler in validatedAsms:
         scoreOutput = ""
         asmScores = validatedAsms[assembler.lower()]
         for type in asmScores.keys():
            outputScore = asmScores[type]
            if asmScores[type] != None and scores[type] == minScore():
               outputScore = "None"
            elif asmScores[type] != None and (type == SCORE_TYPE.FRCBAM or type == SCORE_TYPE.SNP):
               outputScore = -1 * asmScores[type]

            if scoreOutput == "":
               scoreOutput = "%s"%(outputScore)
            else:
               scoreOutput = "%s\t%s"%(scoreOutput, outputScore)
         lapfile.write("%s\t%s\n"%(assembler.upper(), scoreOutput))

   # output failed assemblers
   for file in os.listdir("%s/Assemble/out/"%(_settings.rundir)): 
      if (file.endswith(".failed")):
          assembler = os.path.splitext(os.path.basename(file))[0]
          asmName = getAssemblerName(assembler)
          if asmName.lower() not in validatedAsms:
             lapfile.write("%s\t%s\n"%(asmName, failedOutput))
             lapfile.flush()

   # select best assembly
   global _scores
   if totalRun != 0 or bestAssembler == "": 
      if "%s"%(SCORE_TYPE.ALL) in _scores:
         _scores = SCORE_TYPE.reverse_mapping.keys()

      candidates = dict()
      numCandidatesPerType = int(float(len(input_file_names)) * TOP_PERCENTILE) + 1
      minNumToQualify = ((len(_scores) - len(_skipTypesForCandidates) - 1) / QUALIFY_RATIO)

      # first we select top candidates from each score then all the scores vote
      # get top percentile from all scores
      for type in _scores:
         type = int(type)
         counter = 0
         if type == "":
            continue
         if type == SCORE_TYPE.ALL:
            continue
         if type in _skipTypesForCandidates:
            continue
         if type not in validationScores:
            print "*** metAMOS: Error, type %s is not available, skipping"%(SCORE_TYPE.reverse_mapping[type])
            continue

         sorted_scores = sorted(validationScores[type].iteritems(), key=itemgetter(1), reverse=True)
          
         for tup in sorted_scores:
            if counter <= numCandidatesPerType:
               if tup[0] not in candidates:
                  candidates[tup[0]] = 0
               candidates[tup[0]] += 1
            else:
               break
            counter += 1

      # remove those candidates which did not qualify with at least the required scores (i.e. they have to be in top for several scores)
      bestScores = dict()
      bestAssemblers = dict()
      bestAssemblies = dict()

      for type in validationScores:
         for asm in validationScores[type]:
            if len(candidates) == 0 or (asm in candidates and candidates[asm] >= minNumToQualify):
               if type not in bestScores or float(validationScores[type][asm]) > bestScores[type]:
                  bestScores[type] = float(validationScores[type][asm])
                  bestAssemblers[type] = asm
                  bestAssemblies[type] = "%s/Assemble/out/%s.asm.contig"%(_settings.rundir, asm)

      assemblerVotes = dict()
      assemblyVotes = dict()
      currMaxVote = 0
      processedType = dict()
      ties = dict()

      # for each score, we will take the best assembly and vote for it, the one with the most votes after all scores wins
      for type in _scores:
         if type == "":
            continue

         type = int(type)
         if type == SCORE_TYPE.ALL:
            continue
         if type not in bestAssemblers.keys() or bestAssemblers[type] == None or bestAssemblers[type] == minScore():
            print "Warning: selected score %s was not available, skipping it"%(SCORE_TYPE.reverse_mapping[type].upper())
            continue
         if type in processedType.keys():
            continue
         processedType[type] = True

         # initialize scores to 0
         if bestAssemblers[type] not in assemblerVotes:
            assemblerVotes[bestAssemblers[type]] = 0
            assemblyVotes[bestAssemblies[type]] = 0

         # get the weight this score has
         weight = 1
         try:
            weight = SCORE_WEIGHTS[type]
         except KeyError:
            print "Uknown Type %s using weight %f"%(SCORE_TYPE.reverse_mapping[type], weight)
            weight = 1

         # finally vote for this assembler
         if _settings.VERBOSE:
            print "*** metAMOS: Score %s voted for assembler %s, weight %f"%(SCORE_TYPE.reverse_mapping[type], bestAssemblers[type], weight) 
         assemblerVotes[bestAssemblers[type]] += weight
         assemblyVotes[bestAssemblies[type]] += weight
     
      # now that everyone has voted, the winner is the asm with the most weighted votes
      for assembler in assemblerVotes:
         if float(assemblerVotes[assembler]) > currMaxVote:
            bestAssembler = assembler
            currMaxVote = float(assemblerVotes[assembler])

      # check for ties
      for assembler in assemblerVotes:
         if assemblerVotes[assembler] == currMaxVote:
            ties[assembler] = tieBreakingScores[assembler]
      if len(ties) > 1:
         if _settings.VERBOSE:
            print "*** metAMOS: Assemblers %s are tied after votes, tie breaking using %s."%(",".join(ties), SCORE_TYPE.reverse_mapping[_tieBreakerType])
         currMaxVote = None
         for asm in ties:            
            if currMaxVote == None or float(ties[asm]) > float(currMaxVote):
               bestAssembler = asm
               currMaxVote = float(ties[asm])
         if _settings.VERBOSE:
            print "*** metAMOS: Tie broken in favor of %s"%(bestAssembler)
            
      # sanity check, ensure we picked the same assembly file as assembler. This must always be true as we vote for both at the same time
      currMaxVote = 0
      for assembly in assemblyVotes:
         if assemblyVotes[assembly] > currMaxVote:
            bestAssembly = assembly
            currMaxVote = assemblyVotes[assembly]
         elif assemblyVotes[assembly] == currMaxVote and getAsmName(assembly) == bestAssembler:
            bestAssembly = assembly
            currMaxVote = assemblyVotes[assembly]
      if getAsmName(bestAssembly) != bestAssembler:
         print "Error: inconsistent assembly and assembler chosen %s %s"%(bestAssembly, bestAssembler)
         raise(JobSignalledBreak)

      if _settings.VERBOSE:
         print "*** metAMOS assembler %s selected."%(bestAssembler)   
      run_process(_settings, "unlink %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.PREFIX), "Validate") 

   # finally run quast
   if totalRun != 0:
      selectedReferences = open("%s/Validate/out/%s.ref.selected"%(_settings.rundir, _settings.PREFIX), 'w')
      if "quast" in _validators:
         # recruit a reference
         references = recruitGenomes(_settings,bestAssembly,"%s/refseq"%(_settings.DB_DIR), "%s/Validate/out/recruit"%(_settings.rundir), "Validate", 1)
         if len(references) > 0:
            for g in references:
               selectedReferences.write(os.path.splitext(os.path.basename(g))[0] + "\n")
            runQUAST(asmNames, asmFiles, pairedMin, pairedMax, genomeSize, references[0])
         else:
            print "Warning: could not recruit any references"
      selectedReferences.close()

   if not os.path.exists("%s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.PREFIX)):
      selectedAsm = open("%s/Validate/out/%s.asm.selected"%(_settings.rundir, _settings.PREFIX), 'w')
      selectedAsm.write("%s"%(bestAssembler))
      selectedAsm.close()

      # link the files for subsequent steps to the best assembler
      # this includes assemble/mapreads/validate results
      run_process(_settings, "unlink %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "ln %s %s/Assemble/out/%s.asm.contig"%(bestAssembly, _settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "unlink %s/Assemble/out/%s.linearize.scaffolds.final"%(_settings.rundir, _settings.PREFIX), "Validate")
      if os.path.exists("%s/Assemble/out/%s.linearize.scaffolds.final"%(_settings.rundir, bestAssembler)):
         run_process(_settings, "ln %s/Assemble/out/%s.linearize.scaffolds.final %s/Assemble/out/%s.linearize.scaffolds.final"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")
      else:
         run_process(_settings, "ln %s %s/Assemble/out/%s.linearize.scaffolds.final"%(bestAssembly, _settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "unlink %s/Assemble/out/%s.asm.tigr"%(_settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.asm.tigr %s/Assemble/out/%s.asm.tigr"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "unlink %s/Assemble/out/%s.contig.cnt"%(_settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.contig.cnt %s/Assemble/out/%s.contig.cnt"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "unlink %s/Assemble/out/%s.contig.cvg"%(_settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.contig.cvg %s/Assemble/out/%s.contig.cvg"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "unlink %s/Assemble/out/%s.afg"%(_settings.rundir, _settings.PREFIX), "Validate")
      if os.path.exists("%s/Assemble/out/%s.afg"%(_settings.rundir, bestAssembler)):
         run_process(_settings, "ln %s/Assemble/out/%s.afg %s/Assemble/out/%s.afg"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")

      # rather than linking orfs, we re-run orf finding without the --fast option
      setRunFast(False)
      run_process(_settings, "rm %s/Logs/findorfs.ok"%(_settings.rundir), "Validate")
      FindORFS("%s/Assemble/out/%s.contig.cvg"%(_settings.rundir, _settings.PREFIX), "%s/Assemble/out/%s.faa"%(_settings.rundir, _settings.PREFIX))
      setRunFast(True)
      run_process(_settings, "touch %s/Logs/findorfs.ok"%(_settings.rundir), "Validate")
      run_process(_settings, "unlink %s/Assemble/out/%s.faa"%(_settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "ln %s/FindORFS/out/%s.faa %s/Assemble/out/%s.faa"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "unlink %s/FindRepeats/in/%s.fna"%(_settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "ln %s/FindORFS/out/%s.fna %s/FindRepeats/in/%s.fna"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "unlink %s/FindRepeats/in/%s.faa"%(_settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "ln %s/FindORFS/out/%s.faa %s/FindRepeats/in/%s.faa"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "unlink %s/Annotate/in/%s.fna"%(_settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "ln %s/FindORFS/out/%s.fna %s/Annotate/in/%s.fna"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "unlink %s/Annotate/in/%s.faa"%(_settings.rundir, _settings.PREFIX), "Validate")
      run_process(_settings, "ln %s/FindORFS/out/%s.faa %s/Annotate/in/%s.faa"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Validate")

      for lib in _readlibs: 
         run_process(_settings, "unlink %s/Assemble/out/%s.lib%d.badmates"%(_settings.rundir, _settings.PREFIX, lib.id), "Validate")
         run_process(_settings, "ln %s/Assemble/out/%s.lib%d.badmates %s/Assemble/out/%s.lib%d.badmates"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
         run_process(_settings, "unlink %s/Assemble/out/%s.lib%d.hdr"%(_settings.rundir, _settings.PREFIX, lib.id), "Validate")
         run_process(_settings, "ln %s/Assemble/out/%s.lib%d.hdr %s/Assemble/out/%s.lib%d.hdr"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
         run_process(_settings, "unlink %s/Assemble/out/%s.lib%d.mappedmates"%(_settings.rundir, _settings.PREFIX,lib.id), "Validate")
         run_process(_settings, "ln %s/Assemble/out/%s.lib%d.mappedmates %s/Assemble/out/%s.lib%d.mappedmates"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
         run_process(_settings, "unlink %s/Assemble/out/%s.lib%d.mates_in_diff_contigs"%(_settings.rundir, _settings.PREFIX, lib.id), "Validate")
         run_process(_settings, "ln %s/Assemble/out/%s.lib%d.mates_in_diff_contigs %s/Assemble/out/%s.lib%d.mates_in_diff_contigs"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
         run_process(_settings, "unlink %s/Assemble/out/%s.lib%dcontig.reads"%(_settings.rundir, _settings.PREFIX, lib.id), "Validate")
         run_process(_settings, "ln %s/Assemble/out/%s.lib%dcontig.reads %s/Assemble/out/%s.lib%dcontig.reads"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
         run_process(_settings, "unlink %s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, lib.id), "Validate")
         run_process(_settings, "ln %s/Assemble/out/%s.lib%d.unaligned.fasta %s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, lib.id), "Validate")
         run_process(_settings, "unlink %s/Assemble/out/lib%d.unaligned.fastq"%(_settings.rundir, lib.id), "Validate")
         if os.path.exists("%s/Assemble/out/%s.lib%d.unaligned.fastq"%(_settings.rundir, bestAssembler, lib.id)):
            run_process(_settings, "ln %s/Assemble/out/%s.lib%d.unaligned.fastq %s/Assemble/out/lib%d.unaligned.fastq"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, lib.id), "Validate")

   lapfile.close()

