[CONFIG]
input CONTIGS
name Kraken
output [PREFIX].hits
location cpp/[MACHINE]/kraken/bin
threads --threads
commands kraken [THREADS] -db [DB]/kraken [INPUT] \
         | kraken-filter --db [DB]/kraken --threshold 0.05 > [OUTPUT]
[kraken]
-preload
