#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_pe_filt.fna -d test_velvet -i 3000:4000 
../runPipeline -c phylosift -p 8 -a velvet -d test_velvet -k 55 -t -f Preprocess
