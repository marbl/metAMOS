#bzip2 -d test4.sff.bz2
../initPipeline -d test_isolate -W isolate_test
../runPipeline -p 8 -k 51 -d test_isolate
