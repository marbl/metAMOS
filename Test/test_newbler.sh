#bzip2 -d test4.sff.bz2
../initPipeline -f -l flx -m carsonella_pe_filt.fna  -d test_newbler -i 3000:4000 
../runPipeline -c fcp -p 8 -a newbler -d test_newbler -k 55 -f Preprocess
