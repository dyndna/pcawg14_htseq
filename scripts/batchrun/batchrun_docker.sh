#!/bin/bash

#htseq docker batchrun script:
#v 0.9.4p3 | 05-15 | github.com/dyndna
#@sbamin

LINESTART="$1"
LINEEND="$2"
BATCHFILE="$3"

if [ "$1" == "-h" ] || [ "$1" == "--h" ] || [ "$#" -lt 3 ]; then
	echo "Require three arguments in following order:"
	echo -e "Usage: $(basename $0) <first line to read> <last line to read> <path to batchfile>\n"
	echo -e "Example: $(basename $0) 201 205 ./scripts/batchrun/docker_batchrun_all.txt\n"
	echo -e "If using default ./scripts/batchrun/docker_batchrun_all.txt, make sure you have workdir or git repo at /data/vol2/pcawg14_htseq\n"

	exit 0
else
	n=1
	mkdir -p ./logs

	while read line; do 
	    if [[ ($n -ge "$LINESTART") && ($n -le "$LINEEND")  ]]; then 
	        
	        MYTAG=$(date +%N | sed -e 's/000$//' -e 's/^0//')

	        echo -e "\n#### Attempting to run following docker command ####\n"
	        echo $line 
	        
	        eval $(echo $line)
	        echo -e "\n#### Now, running docker with endtag "$MYTAG" ####\n"
	        
	        CIDFILE=$(ls ./logs/*_${MYTAG}.cidfile)
	        echo "cidfile is at ${CIDFILE}"
	        docker wait $(cat ${CIDFILE})

	        echo -e "\n#### Docker run for ${CIDFILE} ended with exit code "$?" ####\n"
	        sleep 5

	    elif [[ $n -gt "$LINEEND" ]]; then
	        break
	    fi 
	    n=$(( $n + 1 ))
	done < ${BATCHFILE}
fi

# while read line
# do
#   echo -e "\n#### Downloading GNOS analysis ID "$line" ####\n"
#   START=$(date +%M)
#   docker run -i --name="$line" --cidfile="$line.dockerid" -v /data/vol1:/scratch -w /scratch dyndna/pcawg14_htseq:0.9.2 /bin/bash -c "cd bams && gtdownload -c ../.myconfs/.gcsa_Tdb34Ab6dcG -d "$line" -v" | tee -a ../logs/set1_${line}.log
#   END=$(date +%M)
#   DIFF=$(( $END - $START ))
#   echo "Download completed in $DIFF minutes"
# done
