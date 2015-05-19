## pcawg14_htseq

[![Version 0.9.4p1](https://img.shields.io/badge/version-0.9.4p1-brightgreen.svg)](https://github.com/dyndna/pcawg14_htseq) [![Status Stable](https://img.shields.io/badge/status-stable-brightgreen.svg)](https://github.com/dyndna/pcawg14_htseq) [![Docker v0.9.2](https://img.shields.io/badge/docker-dyndna/pcawg14_htseq:0.9.2-brightgreen.svg "docker pull dyndna/pcawg14_htseq:0.9.2")](https://registry.hub.docker.com/u/dyndna/pcawg14_htseq)

[//]: # (https://img.shields.io/badge/status-offline-red.svg)

PCAWG-14 HTSeq RNAseq analysis pipeline

### In host OS:

Install and config docker. To avoid running out of space on root drive, it is recommended to edit `/etc/sysconfig/docker` in CentOS or `/etc/default/docker` in Ubuntu, and set docker container PATH to large storage drive. It is also suggested to use local DNS instead of default Google DNS in docker config to avoid issues with downloading remote data. Sample docker configs are given within git repository (below) under `./docker` directory.

```
docker pull dyndna/pcawg14_htseq:0.9.2
docker images
```

Set *MYWORKDIR* variable. This is a PATH within docker container and not necessarily should be present in host OS. You can keep in host `~/.bash_profile`

    export MYWORKDIR="/scratch"

Go to host directory under which all data needs to be stored. This will be your base directory and mounted as `/scratch` in docker container.

    cd /data/vol1

Configure base directory as follows without changing any dir/file names:

```
pwd # your base directory, e.g., /data/vol1
git clone https://github.com/dyndna/pcawg14_htseq.git
ls -a
```

Directory structure under basedir: `/data/vol1`:

>. .. info scripts set_env

#### Required:

*   Move or symlink your cghub download credential file at `<basedir>/info/cgkey`

#### Optional:

*   Use default GENCODE v19 GTF or move/symlink GTF file at `<basedir>/info/htseq.gtf`

```
cd <basedir>/info
tar xvzf htseq.gtf.tar.gz
mv rnaseq.gc19_extNc.gtf htseq.gtf
```

*   Obtain updated CGHub bam file summary tsv file from https://browser.cghub.ucsc.edu

#### docker command:

Using bam file related attributes from downloaded cghub summary file, run docker command as follows:

>docker run -d -v /data/vol1:$MYWORKDIR -e MYWORKDIR=$MYWORKDIR -w=$MYWORKDIR dyndna/pcawg14_htseq:0.9.2 /bin/bash -c "source ${MYWORKDIR}/set_env/bash_profile && source ${MYWORKDIR}/scripts/htseq_docker_run.sh ***"analysis_id" "filename" "checksum"*** | tee -a ${MYWORKDIR}/logs/htseq_docker_brc.log"

Example:

```
docker run -d -v /data/vol1:$MYWORKDIR -e MYWORKDIR=$MYWORKDIR -w=$MYWORKDIR dyndna/pcawg14_htseq:0.9.2 /bin/bash -c "source ${MYWORKDIR}/set_env/bash_profile && source ${MYWORKDIR}/scripts/htseq_docker_run.sh "b227b026-ef3b-4194-b833-d6386e906587" "PCAWG.057da4ba-421e-4f39-afa8-c7de2ca665e2.TopHat2.v1.bam" "38c067f8289e9c0689fed2c54e9b569e" | tee -a ${MYWORKDIR}/logs/htseq_docker_brc.log"
```

END
