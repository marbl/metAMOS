#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_pe_filt.fna  -d test_newbler -i 3000:4000 
../runPipeline -c phylosift -p 8 -a newbler -d test_newbler -k 55 -t -f Preprocess
