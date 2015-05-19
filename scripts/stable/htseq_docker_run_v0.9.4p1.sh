#!/bin/bash

#htseq.docker script: gene-level expression for pcawg14 extended annotations.
#v 0.9.4p1 | 05-15 | github.com/dyndna
#@sbamin

if [ "$1" == "-h" ] || [ "$1" == "--h" ] || [ "$#" -lt 3 ]; then
	echo "Require three arguments from cghub data in following order:"
	echo -e "Usage: $(basename $0) <analysis_id> <filename> <checksum>\n"
	echo -e "Example: $(basename $0) \"b117855d-bc00-4d34-aedf-e5fdfcc3942a\" \"PCAWG.0d97fa5e-8330-469b-b287-80a4bd881e95.STAR.v1.bam\" \"2a5c25bac2c2aa3eed65ae7fc656eb73\"\n"
	exit 0
else

	MYGTF="./info/htseq.gtf"
	MYBAMID="$1"
	MYBAM="./ingest/${MYBAMID}/${2}"
	MYOUTID=$(basename "$MYBAM" .bam)
	MYSORTEDBAM="./ingest/${MYBAMID}/${MYOUTID}.sorted.bam"
	MD5CGHUB="$3"
	MYKEY="./info/cgkey"

	if [ -z "$MYWORKDIR" ]; then
		echo -e "\n##### ERROR: #####\nMissing ENV variable MYWORKDIR in docker container, unable to initiate analysis.\nSet ENV variable MYWORKDIR to point to mount point inside docker container.\n#### Exiting with ERROR ####\n" >&2
		exit 1
	else
		cd "$MYWORKDIR"

		if [ $(pwd) != "$MYWORKDIR" ]; then
			echo -e "\n##### ERROR: #####\nWORKDIR "$MYWORKDIR" is not matching to current workdir: $(pwd)\nCheck if docker run arguments were set as per README file.\n#### Exiting with ERROR ####\n" >&2
			exit 1
		else
			echo "WORKDIR set to "$MYWORKDIR" and matches $(pwd)"
			mkdir -p {ingest,logs,out}
			ls -alh
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
			gtdownload --version
			echo "######################################################"
			echo -e "\n#### Download BAM from CGHub at "$MYBAM" ####"
			gtdownload -c "$MYKEY" -d "$MYBAMID" -p ingest -v

			MD5BAM=$(md5sum "$MYBAM" | cut -d ' ' -f 1)
			echo -e "\n### BAM_DOWNLOADED_DATESTAMP: $(date) ###\n"

			if [ "$MD5BAM" != "$MD5CGHUB" ]; then
				echo -e "\n#### ERROR: MD5 checksum failed for downloaded bam file. ####\n"
				echo "CGHUB specified checksum is "$MD5CGHUB""
				echo -e "Downloaded bam checksum is "$MD5BAM"\n#### Exiting with ERROR ####\n" >&2
				exit 1
			else
				echo -e "CGHUB specified checksum is "$MD5CGHUB""
				echo -e "Downloaded bam checksum is "$MD5BAM""
				echo -e "\n#### BAM HEADER ####"
				echo "BAM file is at "$MYBAM" with prefix "$MYOUTID""
				samtools view -H "$MYBAM"
				echo "######################################################"
				echo -e "\n#### NAME SORTING BAM ####"
				samtools sort -n "$MYBAM" "./ingest/${MYBAMID}/${MYOUTID}.sorted"
				# samtools index is not useful for name sorted bams.
				echo -e "\n#### SORTED BAM HEADER ####"
				echo -e "\n### SORTING_BAM_DONE_DATESTAMP: $(date) ###\n"

				if [ ! -f "$MYSORTEDBAM" ]; then
					echo -e "\n#### ERROR: Name sorting bam failed, can not continue htseq-count. ####\n#### Exiting with ERROR ####\n" >&2
					exit 1
				else
					echo "SORTED BAM file is at "$MYSORTEDBAM""
					samtools view -H "$MYSORTEDBAM"
					echo -e "\n#### GTF HEADER ####"
					echo "GTF file is at "$MYGTF""
					head "$MYGTF"
					echo "######################################################"
					echo -e "\n## RUNNING htseq-count ##\n"
					
					OUTTIME=$(date +%d%b%y_%H%M%S%Z)

					#for position-sorted bams: default
					#samtools view -F 4 "$MYBAM" | htseq-count -m intersection-nonempty --stranded=no --idattr gene_id -r pos - "$MYGTF" > ./out/"${MYOUTID}".htseq.log 2>&1 > ./out/"${MYOUTID}".htseq.cnt
					
					#for name-sorted bams
					samtools view -F 4 "$MYSORTEDBAM" | htseq-count -m intersection-nonempty --stranded=no --idattr gene_id -r name - "$MYGTF" > ./out/${MYOUTID}_${OUTTIME}.htseq.log 2>&1 > ./out/${MYOUTID}_${OUTTIME}.htseq.cnt

					echo -e "\n### END_DATESTAMP: $(date) ###\n"
					echo "######################################################"
					#END
				fi
			fi
		fi
	fi
fi

# Within host OS, execute following command from any location:
	
	#base directory in host OS, e.g., /data/vol1 or any other location where you cloned pcawg14_htseq git repository.
	#export MYWORKDIR="/scratch" # PATH for docker containter, do not edit it.
	#docker run -d -v /data/vol1:$MYWORKDIR -e MYWORKDIR=$MYWORKDIR -w=$MYWORKDIR dyndna/pcawg14_htseq:0.9.2 /bin/bash -c "source ${MYWORKDIR}/set_env/bash_profile && source ${MYWORKDIR}/scripts/htseq_docker_run.sh "b227b026-ef3b-4194-b833-d6386e906587" "PCAWG.057da4ba-421e-4f39-afa8-c7de2ca665e2.TopHat2.v1.bam" "38c067f8289e9c0689fed2c54e9b569e" | tee -a ${MYWORKDIR}/logs/htseq_docker_brc.log"

#END
