## R script

# pcawg14_htseq flowr pipeline
# https://gitlab.com/gmapps/pcawg14_htseq
# 10-15 | @sbamin
# Based on flowr package, https://github.com/sahilseth/flowr

#### This script should be run from git code repository, ./pcawg14_htseq/flowr/ ####
## How To ##

## DO NOT RUN FROM R; USE /bin/bash SHELL environment

#1. load required modules: R, samtools, genetorrent, htseq-0.6.1.p2, git

#2. edit path to master sample info under "USER EDITABLE CONTENTS" in this script.

#3. clone git repo and checkout "shark" branch to run on HPC Shark at UTMDACC

	#mkdir -p ${HOME}/pipelines && cd ${HOME}/pipelines
	#git clone https://github.com/dyndna/pcawg14_htseq.git
	#cd pcawg14_htseq
	#git checkout shark
	#cd flowr

	# analysis parameters can be modified in pcawg14_htseq.R and pcawg14_htseq.conf ; 
	# CPU, memory, walltime, etc. requirements will be sourced from pcawg14_htseq.csv ;

#4. Run flowr dry run:
	#Rscript run_htseq.R batch1 FALSE 1 10

#Require four mandatory arguments: 
	#First argument: name of running batch without any spaces or special characters except underscore sign _. 
	#Second argument: Supply TRUE if you want to run pipeline without first checking for integrity. 
	#Third and fourth arguments show start and end range from master table, exclude header. 
	#Example: Rscript run_htseq.R batch1 TRUE 25 50 will run flowr on 25th through 50th sample.


#5: If all is well, execute flowr on cluster:

	#5.1 remove dummy directories created during dry run. No harm except cluttered directories if you do not delete those as flowr will create unique directories at each run.
	## rm -rf ${HOME}/runs/batch1* 

	#5.2 Actual run on cluster:
	#Rscript run_htseq.R batch1 TRUE 1 10

	# Done! See logs at ${HOME}/pipelines/pcawg14_htseq/reports/batch1/  ##

########## USER EDITABLE CONTENTS ##########
## fqdir is path to master data.frame in R rds which contains data in character mode with at least following colnames: samplename, analysis_id, file, and url. samplename and analysis ids must be unique across all libraries (column), but can be identical for an individual library(row).
	fqdir = "../master/pilot3_ebi.rds"
	if(!file.exists(fqdir))
		stop(sprintf("master table was not found at %s", fqdir))
######## END USER EDITABLE CONTENTS ########

#### UNLESS YOU KNOW FLOWR, DO NOT EDIT LINES BELOW ####

(WD <- getwd())
if (!is.null(WD)) setwd(WD)

print(sprintf("Current work dir is %s", getwd()))

if(!file.exists("../.git/description"))
  stop("You are not running this pipeline from a git directory, ./pcawg14_htseq/flowr/. Please clone repository first, cd to ./pcawg14_htseq/flowr/ and then run pipeline again.")

args <- commandArgs(TRUE)

if(length(args) < 4) stop(sprintf("You've specified nothing or %s; Require four mandatory arguments: name of running batch without any spaces or special characters except underscore sign _. Second argument: Supply TRUE if you want to run pipeline without first checking for integrity. Third and fourth arguments show start and end range from master table, exclude header. Example: Rscript run_htseq.R batch1 TRUE 25 50 will run flowr on 25th through 50th sample.", args[1]))

myexecute = as.logical(as.character(args[2]))
if(is.na(myexecute)) (myexecute = as.logical(as.character("FALSE")))

print(sprintf("User supplied flowname is: %s", args[1]))
print(sprintf("To run flowr without integrity check: %s", args[1]))

STARTTIME = make.names(format(Sys.time(),"X%d_%b_%y_%H%M%S%Z"))

myflowname = as.character(args[1])
myflowname = sub("\\s+", sprintf("batch_%s", STARTTIME), myflowname)

if(is.na(myflowname)) (myflowname = sprinf("batch_%s", STARTTIME))
print(sprintf("flowname is: %s", myflowname))

## store flowr run logs in batch specific directory ##
myoutdir = sprintf("../reports/%s", myflowname)
dir.create(path = myoutdir, recursive = TRUE, mode = "0775")
if(!file.exists(myoutdir))
	stop(sprintf("No such directory at %s", myoutdir))
print(sprintf("Reports will be archived at: %s", myoutdir))
## Following will source analysis steps from pcawg14_htseq.R ; 
	# analysis parameters can be modified in pcawg14_htseq.R and pcawg14_htseq.conf ; 
	# CPU, memory, walltime, etc. requirements will be sourced from pcawg14_htseq.csv ;
	
	source('pcawg14_htseq.R')

	tmpout = readRDS(file=fqdir)
	print(sprintf("#### master info file has %s samples ####", length(unique(tmpout$samplename))))

	sample_start = as.numeric(args[3])
	if(is.na(sample_start)) (sample_start = 1)

	sample_end = as.numeric(args[4])
	if(is.na(sample_end)) (sample_end = nrow(tmpout))

	out = tmpout[sample_start:sample_end,]

	write.table(out, sprintf("%s/flowr_run_%s_to_%s_%s_%s", myoutdir, sample_start, sample_end, STARTTIME, basename(fqdir)), 
		sep = "\t", quote = FALSE, row.names = FALSE)
	print(sprintf("#### Sample info to be processed was written at %s/flowr_run_%s_to_%s_%s_%s ####", myoutdir, 
		sample_start, sample_end, STARTTIME, basename(fqdir)))

	print(sprintf("#### Flowr will process %s samples, starting at index (header exlcluded) %s and ending at index (header exlcluded) %s of master_info file: %s ####", length(unique(out$samplename)), sample_start, sample_end, fqdir))

	print("###################### DOUBLE CHECK if sample count is correct above, else press CTRL C - 30 seconds wait #######################")

Sys.sleep(30)
myflowobj = vector("list", length(unique(out$samplename)))
myflowdata = vector("list", length(unique(out$samplename)))

for (i in 1:length(unique(out$samplename))) {
	s = unique(out$samplename)[i]

	fobj = get_flow(s, out, myflowname, myoutdir)
	myflowobj[[i]] = fobj

	if(is.flow(fobj)) {
		myflowdata[[i]] = submit_flow(fobj, execute = myexecute, verbose = TRUE)
		names(myflowdata)[i] <- as.character(s)
		names(myflowobj)[i] <- as.character(s)		
	} else {
		myflowdata[[i]] = sprintf("Unable to create flow object: See previous warnings for possible issue. Most likely reason being incorrect master table format for file, %s", fqdir)
		warning(sprintf("Unable to create flow object: See previous warnings for possible issue. Most likely reason being incorrect master table format for file, %s", fqdir))
		names(myflowdata)[i] <- as.character(s)
		names(myflowobj)[i] <- as.character(s)	
	}
}

sessionInfo()
mysessionInfo = sessionInfo()

## move Rplots.pdf to reports dir.
file.rename("Rplots.pdf", file.path(myoutdir, sprintf("flowr_%s_plots_%s.pdf", myflowname, STARTTIME)))

## save detailed log of flowr.
saveRDS(myflowobj, file.path(myoutdir, sprintf("flowr_fobj_%s_log_%s.rds", myflowname, STARTTIME)))
saveRDS(myflowdata, file.path(myoutdir, sprintf("flowr_submit_flow_%s_log_%s.rds", myflowname, STARTTIME)))
save.image(file = file.path(myoutdir, sprintf("flowr_logs_%s_log_%s.RData", myflowname, STARTTIME)))

#END#
