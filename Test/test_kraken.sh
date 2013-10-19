#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_pe_filt.fna -d test_kraken -i 3000:4000 
../runPipeline -u -r -v -c kraken -p 8 -a soap,velvet -d test_kraken -k 55 -t -n FunctionalAnnotation -f Postprocess -z phylum
