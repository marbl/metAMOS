../initPipeline -q -1 test1.fq -2 test2.fq  -d test_phylosift  -i 150:450
../runPipeline -a soap -c phylosift -p 4 -d test_phylosift -k 45 -t -f Postprocess
