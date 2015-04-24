## pcawg14_htseq

>PCAWG-14 HTSeq RNAseq analysis pipeline
>v 0.9.2 | Passing; Unstable

#### Guide to run `HTSeq-count` using PCAWG-14 extended set of annotations

*   Run following in your docker-enabled host system, preferably Unix x64 system with minimum 16 VCPUs and 32 GB RAM (more is better!):

```
whoami #non-root user
cd /scratch #or your preferred storage path.
git clone https://github.com/dyndna/pcawg14_htseq.git
cd pcawg14_htseq
```

*   Set `MYWORKDIR` environment variable to `/home/samin1/pcawg14`

```
export MYWORKDIR="/home/samin1/pcawg14"
echo $MYWORKDIR
```

>You may keep this variable in your bashrc to avoid sourcing each time you login to host system. As of now, it is not possible to edit `MYWORKDIR` variable unless you prefer to rebuild docker image.

```
docker run -d --name="909c7cfb-8f65-4ab3-9189-fdce5cb6dc59" --cidfile="909c7cfb-8f65-4ab3-9189-fdce5cb6dc59.dockerid" -v /scratch/data/samir/pcawg/v3_htseq:$MYWORKDIR -e MYWORKDIR=$MYWORKDIR -e USER=$USER -e USERID=$UID -w=$MYWORKDIR dyndna/pcawg14_htseq:0.9.2 /bin/bash -c "source ${MYWORKDIR}/set_env/bash_profile && source ${MYWORKDIR}/htseq_docker.sh "rnaseq.gc19_extNc.gtf" "909c7cfb-8f65-4ab3-9189-fdce5cb6dc59" "PCAWG.057da4ba-421e-4f39-afa8-c7de2ca665e2.STAR.v1.bam" | tee -a ${MYWORKDIR}/htseq_docker.log"

```
