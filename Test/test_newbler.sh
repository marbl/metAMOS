#bzip2 -d test4.sff.bz2
../initPipeline -l flx -s -m test4.sff -d test_newbler -i 3000:4000 
../runPipeline -c amphora -a newbler -p 4 -d test_newbler -k 55 -t -f Preprocess
