#!/bin/bash

# Based on pcawg14 "htseq_docker_run.sh" script, v0.9.4p4 & docker repo: dyndna/pcawg14_htseq:0.9.2
# README at http://dyndna.github.io/pcawg14_htseq

# cd to workdir.
# it is safe to keep same workdir across multiple containers given storage PATH can handle significant I/O overhead.
cd /data/vol2/pcawg14_htseq

echo -e "\n#### Batchrun on VM ####\n"

./scripts/batchrun/batchrun_docker.sh 1 100 ./scripts/batchrun/docker_batchrun_all.txt > ./logs/batch_1_100.log 2>&1 &

echo "Spawned batch_1_100 on $(date)"

# Wait until previous instance of docker completes downloading bam from cghub and sort it using samtools. Apprx. 1 hour.
# You may edit wait time to optimize pseudo parallelization!
sleep 4800

./scripts/batchrun/batchrun_docker.sh  101 200 ./scripts/batchrun/docker_batchrun_all.txt > ./logs/batch_101_200.log 2>&1 &

echo "Spawned batch_101_200 on $(date)"
sleep 4800

./scripts/batchrun/batchrun_docker.sh 201 300 ./scripts/batchrun/docker_batchrun_all.txt > ./logs/batch_201_300.log 2>&1 &

echo "Spawned batch_201_300 on $(date)"

echo -e "\n#### batch scripts spawned on VM ####\n"

## On your host system, execute this script to run htseq-count on first 300 samples:
# Use GNU Screen or nohup command to keep script alive after exiting terminal session. 
# cd <workdir>
# ./sample_batchscript.sh | tee -a logs/sample_batchscripts_1_300.log

#END
