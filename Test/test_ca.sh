#bzip2 -d test4.sff.bz2
../initPipeline -l flx -s -m test4.sff -d test_ca -i 3000:4000 
../runPipeline -c amphora -a ca -p 4 -d test_ca -k 55 -t -f Preprocess
