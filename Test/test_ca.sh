#bzip2 -d test4.sff.bz2
../initPipeline -l flx -q -m carsonella_pe_filt.fq -d test_ca -i 3000:4000 
../runPipeline -c fcp -a ca -p 4 -d test_ca -k 55 -t -f Preprocess
