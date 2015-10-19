#!/bin/bash

# v 1.8.2 | 10-15 | @sbamin
#set -e # see MYEXIT variable below.

MYKEY="/scratch/genomic_med/gm_dep/samin1/orion/master/mycgkey"
MYBAMINDEX=${1}

#while read line; do
#	echo -e "\n#### Downloading $line at $(date) ####\n"
#	eval $(echo $line)
#done < ${MYBAMS}

echo -e "\n#### Batch begin $(date) ####\n"

while read myindex; do
	echo -e "\n############### Processing ${myindex} #################\n"
	cd /scratch/genomic_med/gm_dep/samin1/orion/workdirs/pcawg14/ingest/
	echo "Work dir is $(pwd)"

	MYGTCMD=$(sed -n "${myindex}"p /scratch/genomic_med/gm_dep/samin1/orion/master/mumbai_cmds.txt)
	echo -e "\n#### Downloading ${myindex} at $(date) ####\n"
	echo "${MYGTCMD}" 

	eval $(echo ${MYGTCMD})
	MYEXIT=$?

	if [ "$MYEXIT" == 0 ]; then

		echo -e "\n############### Download Completed for ${myindex} at $(date) #################\n"
		MYMSUB=$(echo "msub -N pcawg14_${myindex} -vARG1=${myindex},ARG2=${myindex},ARG3=pcawg14 /scratch/genomic_med/gm_dep/samin1/orion/scripts/pcawg14/htseq_nautilus_batchrun.pbs")
		echo -e "\n#### Running HTSeq on index ${myindex} at $(date) ####\n"
		echo "$MYMSUB"
		eval $(echo ${MYMSUB})
		echo -e "\n########### Submitted HTSeq msub job for index ${myindex} at $(date) #############\n"
	else
		echo -e "\n############### ERROR downloading index: ${myindex} at $(date) #################\n"
		echo -e "Following command returned exit status: "$MYEXIT""
		echo "${MYGTCMD}"
		echo -e "\n############### SKIPPING index: ${myindex} ###############\n" 
	fi
done < ${MYBAMINDEX}

echo -e "\n#### Batch end $(date) ####\n"

#END
