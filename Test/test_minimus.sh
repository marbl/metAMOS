../initPipeline -f -m carsonella_pe_filt.fna -d test_amos -i 500:3500
../runPipeline -a amos -c phylosift -g fraggenescan -p 15 -d test_amos -k 55 
