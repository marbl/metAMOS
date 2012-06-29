#bzip2 -d test4.sff.bz2
../initPipeline -l flx -q -m carsonella_pe_filt.fq -d test_soap -i 3000:4000 
../runPipeline -c fcp -p 8 -d test_soap -k 55 -t -f Preprocess,Abundance
