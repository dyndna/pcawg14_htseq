#!/bin/bash

#htseq.docker script: exon-level expression for pcawg14 extended annotations.
#v 0.9.3p3 | 04-15 | github.com/dyndna
#@sbamin

MYGTF="$1"
MYBAM="./ingest/${2}/${3}"
MYOUTID=$(basename "$MYBAM" .bam)
#MYSORTEDBAM="./ingest/${2}/${MYOUTID}.sorted.bam"

if [ -z "$MYWORKDIR" ]; then
	echo -e "\n##### ERROR: #####\nMissing ENV variable MYWORKDIR in docker container, unable to initiate analysis.\nSet ENV variable MYWORKDIR to point to mount point inside docker container.\nExiting." >&2
	exit 1
else
	cd "$MYWORKDIR"

	if [ $(pwd) != "$MYWORKDIR" ]; then
		echo -e "\n##### ERROR: #####\nWORKDIR "$MYWORKDIR" is not matching to current workdir: $(pwd)\nCheck if docker run arguments were set as per README file.\nExiting." >&2
		exit 1
	else
		echo "WORKDIR set to "$MYWORKDIR" and matches $(pwd)"

		echo "######################################################"
		echo -e "\n## SANITY CHECK ##\n"
		echo -e "\n### START_DATESTAMP: $(date) ###\n"
		echo "#### My ENV ####\n"
		env
		echo "PATH: "$PATH""
		echo -e "\n### SOME STATS ###\n"
		echo -e "\n#### SOFTWARE VERSIONS ####\n"
		samtools --version
		htseq-count --help | grep version
		echo -e "\n#### BAM HEADER ####"
		echo "BAM file is at "$MYBAM" with prefix "$MYOUTID""
		samtools view -H "$MYBAM"
		#echo -e "\n#### SORTED BAM HEADER ####"
		#echo "SORTED BAM file is at "$MYSORTEDBAM""
		#samtools view -H "$MYSORTEDBAM"
		echo -e "\n#### GTF HEADER ####"
		echo "GTF file is at "$MYGTF""
		head "$MYGTF"

		echo -e "\n## RUNNING htseq-count ##\n"
		#for position-sorted bams: default
		mkdir -p ${MYWORKDIR}/out
		samtools view -F 4 "$MYBAM" | htseq-count -m intersection-nonempty --stranded=no --idattr exon_id -r pos - "$MYGTF" > ./out/"${MYOUTID}".htseq.exon.log 2>&1 > ./out/"${MYOUTID}".htseq.exon.cnt
		#for name-sorted bams
		#samtools view -F 4 "$MYSORTEDBAM" | htseq-count -m intersection-nonempty --stranded=no --idattr gene_id -r name - "$MYGTF" > "$MYOUTID".htseq.log 2>&1 > "$MYOUTID".htseq.cnt

		echo -e "\n### END_DATESTAMP: $(date) ###\n"
		echo "######################################################"
		#END
	fi
fi

#Execute from any location, preferably outside mapped volume path within host OS:

#export MYWORKDIR="/home/samin1/pcawg14"

#docker run -d --name="exon-6e858862-703f-452d-889c-5d8cd03793b9" --cidfile="exon-6e858862-703f-452d-889c-5d8cd03793b9.dockerid" -v /scratch/data/samir/pcawg/v3_htseq:$MYWORKDIR -e MYWORKDIR=$MYWORKDIR -e USER=$USER -e USERID=$UID -w=$MYWORKDIR pcawg14/htseq:0.9.2 /bin/bash -c "source ${MYWORKDIR}/set_env/bash_profile && source ${MYWORKDIR}/htseq_docker_transid.sh "rnaseq.gc19_extNc.gtf" "6e858862-703f-452d-889c-5d8cd03793b9" "PCAWG.4588e050-e15f-476f-9065-ba747bdd3736.STAR.v1.bam" | tee -a ${MYWORKDIR}/htseq_docker.exon.log"

#END
