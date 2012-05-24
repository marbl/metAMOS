import os, sys, string, subprocess, distutils.util

print "<<Weclome to metAMOS install>>"

silentInstall=False
if (len(sys.argv) > 1):
  if sys.argv[1] == 'silent':
     silentInstall=True
     print "Running in silent mode"

ALLOW_FAST=True
OSTYPE="Linux"
OSVERSION="1"
MACHINETYPE="x86_64"

#identify machine type
p = subprocess.Popen("echo `uname`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(checkStdout, checkStderr) = p.communicate()
if checkStderr != "":
   print "Warning: Cannot determine OS, defaulting to %s"%(OSTYPE)
else:
    OSTYPE = checkStdout.strip()

p = subprocess.Popen("echo `uname -r`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(checkStdout, checkStderr) = p.communicate()
if checkStderr != "":
   print "Warning: Cannot determine OS version, defaulting to %s"%(OSVERSION)
else:
   OSVERSION = checkStdout.strip()

p = subprocess.Popen("echo `uname -m`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(checkStdout, checkStderr) = p.communicate()
if checkStderr != "":
   print "Warning: Cannot determine system type, defaulting to %s"%(MACHINETYPE)
else:
   MACHINETYPE = checkStdout.strip()

if OSTYPE == "Darwin":
   p = subprocess.Popen("echo `gcc --version`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if "Apple" not in checkStdout:
      ALLOW_FAST=False

#check for python version
if (sys.version_info[0] < 2) or (sys.version_info[0] == 2 and sys.version_info[1] < 6):
  print "Python version is %s. metAMOS requires at least 2.6"%(sys.version)
  sys.exit(1)

#check for DBs, etc

if not os.path.exists("./FastQC"):
    print "FastQC not found, optional for Preprocess, download now?"
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        archive = "fastqc_v0.10.0.zip"
        os.system("wget http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/%s" % archive)
        os.system("unzip %s" % archive)
        os.system("rm %s" % archive)
        os.system("chmod u+x FastQC/fastqc")

if not os.path.exists("./Utilities/models"):
    print "Genome models not found, optional for FCP/NB, download now?"
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        archive = "fcp_models.tar.gz"
        os.system("wget ftp://ftp.cbcb.umd.edu/pub/data/metamos/%s -O %s" %(archive, archive))
        #os.system("mv %s ./Utilities/models/." % archive)
        os.system("tar -C ./Utilities/ -xvf %s" % archive)
        os.system("rm %s"%archive)
        os.system("ln -s %s/Utilities/python/taxonomy.txt %s/Utilities/models/taxonomy.txt"%(sys.path[0], sys.path[0]))
        #os.system("chmod u+x Utlities/models")

if not os.path.exists("./Utilities/glimmer-mg"):
    print "Glimmer-MG not found, optional for FindORFS step. Caution, this will take approx. 24 hours to complete, including Phymm download & install. download & install now?"
    if silentInstall:
       dl = 'n'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        archive = "glimmer-mg-0.2.2.tar.gz"
        os.system("wget ftp://ftp.cbcb.umd.edu/pub/data/metamos/%s -O %s" %(archive, archive))
        #os.system("mv %s ./Utilities/models/." % archive)
        os.system("tar -C ./Utilities/ -xvf %s" % archive)
        os.system("rm %s"%archive)
        os.system("python ./Utilities/glimmer-mg/install_glimmer.py")
        #os.system("ln -s %s/Utilities/python/taxonomy.txt %s/Utilities/models/taxonomy.txt"%(sys.path[0], sys.path[0]))
        #os.system("chmod u+x Utlities/models")

if not os.path.exists("./Utilities/DB/refseq_protein.pal"):
    print "refseq protein DB not found, needed for Annotate step, download now?"
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        print "Download and install refseq protein DB.."
        os.system("perl ./Utilities/perl/update_blastdb.pl refseq_protein")
        os.system("mv refseq_protein.00.tar.gz ./Utilities/DB/.")
        os.system("mv refseq_protein.01.tar.gz ./Utilities/DB/.")
        os.system("mv refseq_protein.02.tar.gz ./Utilities/DB/.")
        os.system("mv refseq_protein.03.tar.gz ./Utilities/DB/.")
        
        os.system("tar -C ./Utilities/DB/ -xvf ./Utilities/DB/refseq_protein.00.tar.gz")
        os.system("tar -C ./Utilities/DB/ -xvf ./Utilities/DB/refseq_protein.01.tar.gz")
        os.system("tar -C ./Utilities/DB/ -xvf ./Utilities/DB/refseq_protein.02.tar.gz")
        os.system("tar -C ./Utilities/DB/ -xvf ./Utilities/DB/refseq_protein.03.tar.gz")
        print "    running fastacmd (might take a few min)..."
        os.system("fastacmd -d ./Utilities/DB/refseq_protein -p T -a T -D 1 -o ./Utilities/DB/allprots.faa")

if not os.path.exists("./Utilities/krona/taxonomy.tab"):
    print "ncbi taxonomy file not found, needed for Postprocess, download now?"
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        print "Download and install ncbi taxonomy.."
        os.system("./Utilities/krona/updateTaxonomy.sh")
        #os.system("rm *.dmp")

if not os.path.exists("./AMOS"):
    print "AMOS binaries not found, needed for all steps, download now?"
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        os.system("wget ftp://ftp.cbcb.umd.edu/pub/data/metamos/amos-%s-%s.binaries.tar.gz -O ./amos-binaries.tar.gz"%(OSTYPE, MACHINETYPE))
        os.system("tar -xvf amos-binaries.tar.gz")
        os.system("rm -rf amos-binaries.tar.gz")

if 0 or not os.path.exists("./phylosift"):
   print "PhyloSift binaries not found, optional for Annotate step, download now?"
   if silentInstall:
      dl = 'y'
   else:
      dl = raw_input("Enter Y/N: ")
   if dl == 'y' or dl == 'Y':
      os.system("wget ftp://ftp.cbcb.umd.edu/pub/data/metamos/phylosift-%s-%s-20120523.tar.bz2 -O ./phylosift.tar.bz2"%(OSTYPE, MACHINETYPE))
      os.system("tar -xvjf phylosift.tar.bz2")
      os.system("rm -rf phylosift.tar.bz2")

if not os.path.exists("./CA"):
   print "Celera Assembler binaries not found, optional for Assemble step, download now?"
   if silentInstall:
      dl = 'y'
   else:
      dl = raw_input("Enter Y/N: ")
   if dl == 'y' or dl == 'Y':
      if OSTYPE == 'Linux' and MACHINETYPE == "x86_64":
         os.system("wget http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-7.0/wgs-7.0-PacBio-Linux-amd64.tar.bz2")
         os.system("bunzip2 wgs-7.0-PacBio-Linux-amd64.tar.bz2")
         os.system("tar xvf wgs-7.0-PacBio-Linux-amd64.tar")
         os.system("rm -rf wgs-7.0-PacBio-Linux-amd64.tar")
      else:
         os.system("wget http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-7.0/wgs-7.0.tar.bz2 -O wgs-7.0.tar.bz2")
         os.system("bunzip2 wgs-7.0.tar.bz2")
         os.system("tar xvf wgs-7.0.tar")
         os.system("rm -rf wgs-7.0.tar")
         # patch CA to support PacBio sequences and non-apple compilers on OSX
         if not ALLOW_FAST:
            os.system("cd wgs-7.0/kmer/ && cp configure.sh configure.original")
            os.system("cd wgs-7.0/kmer/ && cat configure.original |sed s/\-fast//g > configure.sh")
         os.system("cd wgs-7.0/src/ && cp AS_global.h AS_global.original")
         os.system("cd wgs-7.0/src/ && cat AS_global.original | sed 's/AS_READ_MAX_NORMAL_LEN_BITS.*11/AS_READ_MAX_NORMAL_LEN_BITS      15/g' > AS_global.h")
         os.system("cd wgs-7.0/kmer && ./configure.sh && gmake install")
         os.system("cd wgs-7.0/src && gmake")
      os.system("mv wgs-7.0 CA")

# make sure we have setuptools available
sys.path.append(sys.path[0] + os.sep + "Utilities" + os.sep + "python")
from get_setuptools import use_setuptools
use_setuptools()

#os.system("
print "Run setup.py.."
os.system("python setup.py install_scripts --install-dir=`pwd` build_ext")
#print "Compile & optimize"
#distutils.util.byte_compile(['./runPipeline.py'],optimize=2,force=True)
#os.system("chmod a+wrx runPipeline.pyo")
os.system("mv runPipeline.py runPipeline")
os.system("mv initPipeline.py initPipeline")
