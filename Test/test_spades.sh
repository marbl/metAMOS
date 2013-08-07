#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_pe_filt.fna -d test_spades -i 3000:4000 
../runPipeline -r -v -c fcp -p 8 -a spades -d test_spades -k 55 -t -f FunctionalAnnotation,Postprocess -z phylum
