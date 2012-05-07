#bzip2 -d test4.sff.bz2
../initPipeline -l flx -q -m lib1.seq -d test_soap -i 3000:4000 
../runPipeline -c phylosift -p 8 -a soap -d test_soap -k 55 -t -f Preprocess,Abundance
