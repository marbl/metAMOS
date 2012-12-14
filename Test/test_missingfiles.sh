../initPipeline -l flx -f -m carsonella_pe_filt.fna.fake -d test_missing -i 3000:4000 
../initPipeline -l flx -f -1 carsonella_pe_filt.1.fna.fake -2 carsonella_pe_filt.2.fna.fake -i 3000:4000
../initPipeline -l flx -f -1 carsonella_pe_filt.fna.fake -d test_missing
