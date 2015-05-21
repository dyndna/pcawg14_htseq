### PCAWG-14 HTSeq RNAseq analysis pipeline

[![Version 0.9.4p4](https://img.shields.io/badge/version-0.9.4p4-brightgreen.svg)](https://github.com/dyndna/pcawg14_htseq/releases/tag/0.9.4p4) [![Status Stable](https://img.shields.io/badge/status-stable-brightgreen.svg)](https://github.com/dyndna/pcawg14_htseq/releases/tag/0.9.4p4) [![Dockerfile](https://img.shields.io/badge/Dockerfile-v0.9.2-brightgreen.svg "Build Instructions")](https://registry.hub.docker.com/u/dyndna/docker-for-pcawg14-htseq/dockerfile/) [![Docker v0.9.2](https://img.shields.io/badge/docker-dyndna/pcawg14_htseq:0.9.2-brightgreen.svg "docker pull dyndna/pcawg14_htseq:0.9.2")](https://registry.hub.docker.com/u/dyndna/pcawg14_htseq)

[//]: # (https://img.shields.io/badge/status-offline-red.svg)

#### Acknowledgement:

*   Morten Muhlig Nielsen, PhD & Prof. Jakob Skou Pederson, PhD
*   [Laurent Jourdren](https://github.com/jourdren), Genomic Paris Centre
*   [André Kahles](https://github.com/akahles), PhD  
*   [Nuno Fonseca](https://github.com/nunofonseca), PhD   

***

### <span class="octicon octicon-tools"></span> Required configurations host unix system:

#### <span class="octicon octicon-repo-pull"></span> Install and config docker

To avoid running out of space on root drive, it is recommended to edit `/etc/sysconfig/docker` in CentOS or `/etc/default/docker` in Ubuntu, and set docker container PATH to large storage drive. It is also suggested to use local DNS instead of default Google DNS in docker config to avoid issues with downloading remote data. Sample docker configs are given within git repository (below) under `./info/docker_host_configs` directory.

```
docker pull dyndna/pcawg14_htseq:0.9.2
docker images
```

If downloaded docker image has name other than *dyndna/pcawg14_htseq:0.9.2*, rename or create alias image as follows:

    docker tag <IMAGE ID> dyndna/pcawg14_htseq:0.9.2

#### Set *MYWORKDIR* variable

This is a PATH within docker container and not necessarily should be present in host OS. You can keep in host `~/.bash_profile`

    export MYWORKDIR="/scratch"

#### Set up base directory

Go to host directory under which all data needs to be stored. Please read WARNING note below. Ideally, base directory should be freshly created directory under a non-system storage or scratch drive.

    cd /data/vol1

##### <span class="octicon octicon-alert"></span> WARNING <span class="octicon octicon-stop"></span>

Current docker version 1.6.0 continues to have inherent security weakness by not separating root privileges between host system and docker container(s). This will results in allowing docker container to gain root privileges for host system. You may force ownership to some dummy user, *foo* by passing `-e USER=$USER -e USERID=$UID` variables in `docker run` argument at `./scripts/batchrun/docker_batchrun_all.txt`, and thus disallow docker container to gain root privileges. If you do so or in any other case, **PLEASE DO NOT MOUNT** root block device or `\` or $HOME `/~` as base directory, and instead use separate drive or block device other than disk where host system is mounted. However, keep in mind that binding non-root UID will *recursively and irreversibly change owner permission* of all files and directories under base directory in host system. It goes without saying that **this can easily break your system, including lock you out of system if you mount system-level directories or entire user home directory as docker base directory.**

docker works best if used with care! More on security issues and how to minimize vulnerabilities: https://docs.docker.com/articles/security

Dockerfile for image, *dyndna/pcawg14_htseq:0.9.2* can be seen at https://registry.hub.docker.com/u/dyndna/docker-for-pcawg14-htseq/dockerfile/ If you build *dyndna/pcawg14_htseq:0.9.2* from source, image may have different name. If so, rename it as follows:

    docker tag <IMAGE ID> dyndna/pcawg14_htseq:0.9.2

#### <span class="octicon octicon-file-directory"></span> Set up work dir

Your work dir will be mounted as `/scratch` in docker container. Configure work directory as follows without changing any dir/file names:

```
pwd # your base directory, e.g., /data/vol1
git clone https://github.com/dyndna/pcawg14_htseq.git
cd pcawg14_htseq
ls
```

Directory structure under work dir: `<basedir>/pcawg14_htseq/`

>info  LICENSE  README.md  sample_batchscript.sh  scripts  set_env

##### Required edits in workdir

*   Move or symlink your cghub download credential file at `<workdir>/info/cgkey`

*   Use default extended annotation GTF, *rnaseq.gc19_extNc.gtf* or move/symlink your GTF file at `<workdir>/info/htseq.gtf`

```
cd <workdir>/info
tar xvzf htseq.gtf.tar.gz
mv rnaseq.gc19_extNc.gtf htseq.gtf
md5sum htseq.gtf
```

If you are using default extended annotation GTF derived from https://www.synapse.org/#!Synapse:syn3606092, MD5 checksum for `htseq.gtf` must match *48245eeff794b9e686466a80e66905c9*

bed to gtf conversion for syn3606092 file, *rnaseq.gc19_extNc.bed* was done as follows:

Ref.: http://onetipperday.blogspot.com/2012/08/convert-bed-to-gtf.html

    #wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToGenePred
    #wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
    <workdir>/scripts/bedToGenePred rnaseq.gc19_extNc.bed rnaseq.gc19_extNc.btgp.genepred
    <workdir>/scripts/genePredToGtf file rnaseq.gc19_extNc.btgp.genepred rnaseq.gc19_extNc.gptg.gtf

    cat rnaseq.gc19_extNc.gptg.gtf | awk 'BEGIN{FS="::|transcript_id"}{printf($1"::"$2"::"$3"\"; transcript_id"$5"::"$6"::"$7"::"$8"::"$9"::"$10"::"$11"\n")}' | sed -e "s/^chrM/MT/g;s/^chr//g"  > rnaseq.gc19_extNc.gtf

*   Obtain updated CGHub bam file summary tsv file from https://browser.cghub.ucsc.edu

#### docker command:

##### <span class="octicon octicon-playback-play"></span> Single sample run:

Using bam file related attributes from downloaded cghub summary file, run docker command as follows:

>docker run -d -v ***basedir***/pcawg14_htseq:$MYWORKDIR -e MYWORKDIR=$MYWORKDIR -w=$MYWORKDIR dyndna/pcawg14_htseq:0.9.2 /bin/bash -c "source ${MYWORKDIR}/set_env/bash_profile && source ${MYWORKDIR}/scripts/htseq_docker_run.sh ***"analysis_id" "filename" "checksum"*** | tee -a ${MYWORKDIR}/logs/htseq_docker_brc.log"

Example:

    docker run -d -v /data/vol1/pcawg14_htseq:$MYWORKDIR -e MYWORKDIR=$MYWORKDIR -w=$MYWORKDIR dyndna/pcawg14_htseq:0.9.2 /bin/bash -c "source ${MYWORKDIR}/set_env/bash_profile && source ${MYWORKDIR}/scripts/htseq_docker_run.sh "b227b026-ef3b-4194-b833-d6386e906587" "PCAWG.057da4ba-421e-4f39-afa8-c7de2ca665e2.TopHat2.v1.bam" "38c067f8289e9c0689fed2c54e9b569e" | tee -a ${MYWORKDIR}/logs/htseq_docker_brc.log"

##### <span class="octicon octicon-playback-fast-forward"></span> Batch run:

First, make format for `docker run` command as per `./scripts/batchrun/docker_batchrun_all.txt` using R script, `./scripts/batchrun/pcawg14_htseq_docker_batchrun.Rmd`

Then, go to workdir and execute following command to run sample 1 to 100.

    cd <workdir>
    ./scripts/batchrun/batchrun_docker.sh 1 100 ./scripts/batchrun/docker_batchrun_all.txt > ./logs/batch_1_100.log 2>&1 &

You can do pseudo parallelization by executing above script after fixed wait time, e.g., after 1 hour or so for allowing preceding sample to be downloaded and name sorted. Here is an example batch script to run three samples in parallel and process upto 300 samples. Please check `./sample_batchscript.sh` to confirm/edit your workdir. You may use GNU Screen or `nohup` command to keep script alive after exiting terminal session. 

    ./sample_batchscript.sh | tee -a logs/brc1_batchscripts_11_300.log

#### logging:

Although not clean, this pipeline will keep master log in `./logs/htseq_docker_brc.log` and individual batch run summary under `./logs/batch_*.log`

To know how many samples have been processed successfully, run `ls ./processed/*.gto | wc -l`.

END

