#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_pe_filt.fna -d test_fcp -i 3000:4000 
../runPipeline -v -c fcp -p 8 -a soap -d test_fcp -k 55 -t -f Preprocess -z phylum
