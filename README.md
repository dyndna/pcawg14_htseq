## pcawg14_htseq

[![Version 0.9.4p1](https://img.shields.io/badge/version-0.9.4p1-brightgreen.svg)](https://github.com/dyndna/pcawg14_htseq) [![Status Stable](https://img.shields.io/badge/status-stable-brightgreen.svg)](https://github.com/dyndna/pcawg14_htseq) [![Docker v0.9.2](https://img.shields.io/badge/docker-dyndna/pcawg14_htseq:0.9.2-brightgreen.svg "docker pull dyndna/pcawg14_htseq:0.9.2")](https://registry.hub.docker.com/u/dyndna/pcawg14_htseq)

[//]: # (https://img.shields.io/badge/status-offline-red.svg)

PCAWG-14 HTSeq RNAseq analysis pipeline

### In host OS:

#### Install and config docker

To avoid running out of space on root drive, it is recommended to edit `/etc/sysconfig/docker` in CentOS or `/etc/default/docker` in Ubuntu, and set docker container PATH to large storage drive. It is also suggested to use local DNS instead of default Google DNS in docker config to avoid issues with downloading remote data. Sample docker configs are given within git repository (below) under `./docker` directory.

```
docker pull dyndna/pcawg14_htseq:0.9.2
docker images
```

#### Set *MYWORKDIR* variable

This is a PATH within docker container and not necessarily should be present in host OS. You can keep in host `~/.bash_profile`

    export MYWORKDIR="/scratch"

#### Set up base directory

Go to host directory under which all data needs to be stored.

    cd /data/vol1

#### Set up work dir

Your work dir will be mounted as `/scratch` in docker container. Configure work directory as follows without changing any dir/file names:

```
pwd # your base directory, e.g., /data/vol1
git clone https://github.com/dyndna/pcawg14_htseq.git
cd pcawg14_htseq
ls
```

Directory structure under work dir: `<basedir>/pcawg14_htseq/`

>info  LICENSE  README.md  scripts  set_env

##### Required edits in workdir

*   Move or symlink your cghub download credential file at `<workdir>/info/cgkey`

*   Use default GENCODE v19 GTF, *rnaseq.gc19_extNc.gtf* or move/symlink your GTF file at `<workdir>/info/htseq.gtf`

```
cd <workdir>/info
tar xvzf htseq.gtf.tar.gz
mv rnaseq.gc19_extNc.gtf htseq.gtf
md5sum htseq.gtf
```

If you are using extended annotation GTF derived from https://www.synapse.org/#!Synapse:syn3606092, MD5 checksum for `htseq.gtf` must match *48245eeff794b9e686466a80e66905c9*

bed to gtf conversion for syn3606092 file, *rnaseq.gc19_extNc.bed* was done as follows:

Ref.: http://onetipperday.blogspot.com/2012/08/convert-bed-to-gtf.html

    #wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToGenePred
    #wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
    <workdir>/scripts/bedToGenePred rnaseq.gc19_extNc.bed rnaseq.gc19_extNc.btgp.genepred
    <workdir>/scripts/genePredToGtf file rnaseq.gc19_extNc.btgp.genepred rnaseq.gc19_extNc.gptg.gtf

    cat rnaseq.gc19_extNc.gptg.gtf | awk 'BEGIN{FS="::|transcript_id"}{printf($1"::"$2"::"$3"\"; transcript_id"$5"::"$6"::"$7"::"$8"::"$9"::"$10"::"$11"\n")}' | sed -e "s/^chrM/MT/g;s/^chr//g"  > rnaseq.gc19_extNc.gtf

*   Obtain updated CGHub bam file summary tsv file from https://browser.cghub.ucsc.edu

#### docker command:

Using bam file related attributes from downloaded cghub summary file, run docker command as follows:

>docker run -d -v ***basedir***/pcawg14_htseq:$MYWORKDIR -e MYWORKDIR=$MYWORKDIR -w=$MYWORKDIR dyndna/pcawg14_htseq:0.9.2 /bin/bash -c "source ${MYWORKDIR}/set_env/bash_profile && source ${MYWORKDIR}/scripts/htseq_docker_run.sh ***"analysis_id" "filename" "checksum"*** | tee -a ${MYWORKDIR}/logs/htseq_docker_brc.log"

Example:

    docker run -d -v /data/vol1/pcawg14_htseq:$MYWORKDIR -e MYWORKDIR=$MYWORKDIR -w=$MYWORKDIR dyndna/pcawg14_htseq:0.9.2 /bin/bash -c "source ${MYWORKDIR}/set_env/bash_profile && source ${MYWORKDIR}/scripts/htseq_docker_run.sh "b227b026-ef3b-4194-b833-d6386e906587" "PCAWG.057da4ba-421e-4f39-afa8-c7de2ca665e2.TopHat2.v1.bam" "38c067f8289e9c0689fed2c54e9b569e" | tee -a ${MYWORKDIR}/logs/htseq_docker_brc.log"

END
