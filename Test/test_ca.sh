#bzip2 -d test4.sff.bz2
../initPipeline -f -m carsonella_pe_filt.fna -d test_ca -i 3000:4000 
../runPipeline -c kraken -a ca -p 4 -d test_ca -k 55 -f Preprocess
