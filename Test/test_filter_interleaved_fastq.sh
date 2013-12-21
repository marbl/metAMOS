#bzip2 -d test4.sff.bz2
../initPipeline -q -m carsonella_pe_filt.fq -d test_filter_ifq -i 3000:4000 
../runPipeline -r -v -c kraken -p 8 -a soap -d test_filter_ifq -k 55 -t metamos -f FunctionalAnnotation,Postprocess -n Classify,Propagate -z phylum
