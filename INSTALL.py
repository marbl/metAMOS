import os, sys, string, subprocess, distutils.util, check_install, site

user_home = os.environ["HOME"]
print "<<Welcome to metAMOS install>>"


#add access to utils.py, for utils dir
INITIAL_SRC   = "%s%ssrc"%(sys.path[0], os.sep)
sys.path.append(INITIAL_SRC)
import utils
sys.path.append(utils.INITIAL_UTILS)

shellv = os.environ["SHELL"]
#add libs to pythonpath

#add site dir
site.addsitedir(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python")

if "PYTHONPATH" not in os.environ:
   os.environ["PYTHONPATH"] = ""
os.environ["PYTHONPATH"]+=utils.INITIAL_UTILS+os.sep+"python"+os.pathsep
os.environ["PYTHONPATH"] += utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.pathsep
os.environ["PYTHONPATH"] += utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"+os.pathsep
sys.path.append(utils.INITIAL_UTILS+os.sep+"python")
sys.path.append(utils.INITIAL_UTILS+os.sep+"python" + os.sep+"lib"+ os.sep+"python")

if 'bash' in shellv:
   os.system("export PYTHONPATH=%s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"))
   os.system("export PYTHONPATH=%s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"))
else:
   os.system("setenv PYTHONPATH %s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"))
   os.system("setenv PYTHONPATH %s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"))

if not os.path.exists("%s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"):
    os.system("mkdir %s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib")
if not os.path.exists("%s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"):
    os.system("mkdir %s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python")

silentInstall=False
if (len(sys.argv) > 1):
  if sys.argv[1] == 'silent':
     silentInstall=True
     print "Running in silent mode"

ALLOW_FAST=True
OSTYPE="Linux"
OSVERSION="1"
MACHINETYPE="x86_64"
kronaTools = "KronaTools-2.2"

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


if not os.path.exists("./Utilities/config/usage.ok"):
    print "MetAMOS would like to record anonymous usage statistics, is this ok ? "
    dl = 'n'
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        os.system("echo ok > ./Utilities/config/usage.ok")

#check for DBs, etc
if not os.path.exists("./Utilities/cpp/%s-%s/metaphylerClassify"%(OSTYPE, MACHINETYPE)) or not os.path.exists("./Utilities/perl/metaphyler/markers/markers.protein") or not os.path.exists("./Utilities/perl/metaphyler/markers/markers.dna"):
    print "Metaphyler (latest version) not found, optional for Annotate, download now?"
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        os.system("wget http://metaphyler.cbcb.umd.edu/MetaPhylerV1.25.tar.gz -O metaphyler.tar.gz")
        os.system("tar -C ./Utilities/perl/ -xvf metaphyler.tar.gz")
        os.system("mv ./Utilities/perl/MetaPhylerV1.25 ./Utilities/perl/metaphyler")
        os.system("perl ./Utilities/perl/metaphyler/installMetaphyler.pl")
        os.system("cp ./Utilities/perl/metaphyler/metaphylerClassify ./Utilities/cpp/%s-%s/metaphylerClassify"%(OSTYPE, MACHINETYPE))
        #os.system("cp ./Utilities/perl/metaphyler/metaphylerClassify ./Utilities/cpp/%s-%s/metaphylerClassify"%(OSTYPE, MACHINETYPE))
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

#not needed (for now)
#if not os.path.exists("./Utilities/cpp/%s-%s/samtools"%(OSTYPE, MACHINETYPE)):
#       os.system("wget http://sourceforge.net/projects/samtools/files/samtools/0.1.17/samtools-0.1.17.tar.bz2 -O samtools.tar.bz2")
#       os.system("tar -C ./Utilities/cpp/%s-%s/ -xvf samtools.tar.bz2"%(OSTYPE, MACHINETYPE))
#       os.system("cd ./Utilities/cpp/%s-%s"%(OSTYPE,MACHINETYPE))
#       os.system("make")
#       os.system("rm samtools.tar.bz2")
 
if not os.path.exists("./Utilities/DB/uniprot_sprot.fasta"):
    print "Uniprot/Swissprot DB not found, optional for Functional Annotation, download now?"
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        archive = "uniprot.tar.gz"
        os.system("wget ftp://ftp.cbcb.umd.edu/pub/data/metamos/%s -O %s" %(archive, archive))
        os.system("tar -C ./Utilities/DB/ -xvf %s" % archive)
        os.system("rm %s"%archive)


if not os.path.exists("./Utilities/models") or not os.path.exists("./Utilities/DB/blast_data"):
    print "Genome models not found, optional for FCP/NB, download now?"
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        archive = "fcp_models.tar.gz"
        os.system("wget ftp://ftp.cbcb.umd.edu/pub/data/metamos/%s -O %s" %(archive, archive))
        os.system("rm -rf ./Utilities/DB/blast_data")
        os.system("rm -rf ./Utilities/models")
        os.system("tar -C ./Utilities/ -xvf %s" % archive)
        os.system("rm %s"%archive)

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

if not os.path.exists("./AMOS") or 0:
    print "AMOS binaries not found, needed for all steps, download now?"
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        os.system("wget ftp://ftp.cbcb.umd.edu/pub/data/metamos/amos-3.2-BETA-%s-%s.binaries.tar.gz -O ./amos-binaries.tar.gz"%(OSTYPE, MACHINETYPE))
        os.system("tar -xvf amos-binaries.tar.gz")
        os.system("rm -rf amos-binaries.tar.gz")

if 1:
   # or not os.path.exists("./Utilities/python/psutil"):
   fail = 0
   try:
       import psutil
   except ImportError:
       print "psutil not found, required for memory usage estimation, download now?"
       fail = 1
   if not fail or silentInstall:
       dl = 'y'
   else:
       dl = raw_input("Enter Y/N: ")
   if fail and (dl == 'y' or dl == "Y"):

       os.system("wget http://psutil.googlecode.com/files/psutil-0.6.1.tar.gz -O ./psutil.tar.gz")
       os.system("tar -C ./Utilities/python -xvf psutil.tar.gz")
       os.system("mv ./Utilities/python/psutil-0.6.1 ./Utilities/python/psutil")
       os.chdir("./Utilities/python/psutil")
       os.system("python setup.py install --home=%spython"%(utils.INITIAL_UTILS+os.sep))
       os.chdir("%s"%(sys.path[0]))
       os.system("rm -rf psutil.tar.gz")
if 1:
   #not os.path.exists("./Utilities/python/cython"):
   fail = 0
   try:
       import cython
   except ImportError:
       print "cython modules not found, necessary for c-compiling python code, download now?"
       fail = 1
   if not fail or silentInstall:
       dl = 'y'
   else:
       dl = raw_input("Enter Y/N: ")
   if fail and (dl == 'y' or dl == "Y"):
       os.system("wget https://github.com/cython/cython/archive/master.zip -O ./cython.zip")
       os.system("unzip ./cython.zip")
       os.system("mv ./cython-master ./Utilities/python/cython")
       os.chdir("./Utilities/python/cython")
       os.system("python setup.py install --home=%spython"%(utils.INITIAL_UTILS+os.sep))
       os.chdir(sys.path[0])
       os.system("rm -rf cython.zip")
       #os.system("tar -C ./Utilities/python -xvf cython.tar.gz")
       #os.system("mv ./Utilities/python/pysam-0.6 ./Utilities/python/pysam")

if 1:
   # or not os.path.exists("./Utilities/python/pysam"):
   fail = 0
   try:
       import pysam
   except ImportError:
       print "pysam python modules not found, necessary for bowtie2 alignments, download now?"
       fail = 1

   if not fail or silentInstall:
       dl = 'y'
   else:
       dl = raw_input("Enter Y/N: ")
   if fail and (dl == 'y' or dl == "Y"):
       os.system("wget http://pysam.googlecode.com/files/pysam-0.6.tar.gz -O ./pysam.tar.gz")
       os.system("tar -C ./Utilities/python -xvf pysam.tar.gz")
       os.system("mv ./Utilities/python/pysam-0.6 ./Utilities/python/pysam")
       #for root install
       #os.system("sudo python ./Utilities/python/pysam/setup.py install")
       os.chdir("./Utilities/python/pysam")
       os.system("python setup.py install --home=%spython"%(utils.INITIAL_UTILS+os.sep))
       os.chdir(sys.path[0])
       os.system("rm -rf pysam.tar.gz")
       #os.system("ln -s %s/Utilities/python/taxonomy.txt %s/Utilities/models/taxonomy.txt"%(sys.path[0], sys.path[0]))
if 0 or not os.path.exists("./phylosift"):
   print "PhyloSift binaries not found, optional for Annotate step, download now?"
   if silentInstall:
      dl = 'y'
   else:
      dl = raw_input("Enter Y/N: ")
   if dl == 'y' or dl == 'Y':
      #os.system("wget ftp://ftp.cbcb.umd.edu/pub/data/metamos/phylosift-Linux-x86_6-20120523.tar.bz2 -O ./phylosift.tar.bz2"%(OSTYPE, MACHINETYPE))
      #phylosift OSX binaries included inside Linux X86_64 tarball..
      os.system("wget ftp://ftp.cbcb.umd.edu/pub/data/metamos/phylosift-Linux-x86_64-20120523.tar.bz2 -O ./phylosift.tar.bz2")
      os.system("tar -xvjf phylosift.tar.bz2")
      os.system("rm -rf phylosift.tar.bz2")

if not os.path.exists("./CA") or 0:
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

if not os.path.exists("KronaTools") or 0:
    print "KronaTools not found, needed for Postprocess, download now?"
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        # TODO: KronaTools should be on the FTP site for robustness to URL changes
        os.system("wget 'ftp://ftp.cbcb.umd.edu/pub/data/metamos/" + kronaTools + ".tar' -O %s.tar"%(kronaTools))
        os.system("tar -xvf %s.tar"%(kronaTools))
        os.system("rm -rf %s.tar"%(kronaTools))
        os.system("mv %s KronaTools"%(kronaTools))
        os.system("cd KronaTools && ./install.pl --prefix=.")

if not os.path.exists("KronaTools/taxonomy/taxonomy.tab") or 0:
    print "KronaTools taxonomy data not found, needed for Postprocess, download now (will take around 20 minutes)?"
    if silentInstall:
       dl = 'y'
    else:
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        os.system("cd KronaTools && ./updateTaxonomy.sh")

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

#print sys.path[0]
validate_install = 1
if validate_install:
    rt = check_install.validate_dir(sys.path[0].strip(),'required_file_list.txt')
    if rt == -1:
        print "MetAMOS not properly installed, please reinstall or contact development team for assistance"
        sys.exit(1)
    
