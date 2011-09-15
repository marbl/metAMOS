import os, sys, string, distutils.util

print "<<Weclome to metAMOS install>>"
#check for DBs, etc
if not os.path.exists("./Utilities/DB/refseq_protein.pal"):
    print "refseq protein DB not found, needed for Annotate step, download now?"
    dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        print "Download and install refseq protein DB.."
        os.system("perl ./Utilities/perl/update_blastdb.pl refseq_protein")
        os.system("mv refseq_protein.00.tar.gz ./Utilities/DB/.")
        os.system("mv refseq_protein.01.tar.gz ./Utilities/DB/.")
        os.system("mv refseq_protein.02.tar.gz ./Utilities/DB/.")
        os.system("mv refseq_protein.03.tar.gz ./Utilities/DB/.")
        
        os.system("tar -C ./Utilities/DB/ -xvf ./Utilities/DB/refseq_protein.00.tar")
        os.system("tar -C ./Utilities/DB/ -xvf ./Utilities/DB/refseq_protein.01.tar")
        os.system("tar -C ./Utilities/DB/ -xvf ./Utilities/DB/refseq_protein.02.tar")
        os.system("tar -C ./Utilities/DB/ -xvf ./Utilities/DB/refseq_protein.03.tar")

        os.system("fastacmd -d ./Utilities/DB/refseq_protein -p T -a T -D 1 -o ./Utilities/DB/allprots.faa")

if not os.path.exists("./Utilities/krona/taxonomy.tab")
    print "ncbi taxonomy file not found, needed for Postprocess, download now?"
    dl = raw_input("Enter Y/N: ")
    print "Download and install ncbi taxonomy.."
    os.system("./Utilities/krona/updateTaxonomy.sh")
    os.system("rm *.dmp")
    
#os.system("
print "Run setup.py.."
os.system("python setup.py install_scripts --install-dir=`pwd`")
#print "Compile & optimize"
#distutils.util.byte_compile(['./runPipeline.py'],optimize=2,force=True)
#os.system("chmod a+wrx runPipeline.pyo")
#os.system("mv runPipeline.pyo runPipeline")
