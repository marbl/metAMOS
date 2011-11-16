../createProject -f -m carsonella_pe_filt.fna -d test12  -i 150:450
cat test12/pipeline.ini
rm -rf test12
../createProject -d test12 -q -1 test12.fq -1 test1.fq -2 test2.fq -i 150:450
cat test12/pipeline.ini
rm -rf test12
../createProject -d test12 -q -1 test12.fq -1 test1.fq -2 test2.fq -i 150:450
cat test12/pipeline.ini
rm -rf test12
../createProject -d test12 -i 150:450 -q -1 test1.fq,test12.fq -2 test2.fq 
cat test12/pipeline.ini
rm -rf test12
