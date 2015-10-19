## R script
# v 2.0.2 | 10-15 | @sbamin

# pcawg14 htseq pipeline

#### Sample bash command ####
#to run on first 12 bam files with 2 files in parallel.
#Rscript script.R 1 12 pcawg14 2>&1 | tee -a err.log
################################

args <- commandArgs(TRUE)
print(sprintf("Supplied arguments are: %s", args))
print(paste0("##########################################"))

library(uuid)

myproject <- as.character(args[3])

batch_uuid <- UUIDgenerate()
outFile <- sprintf("/scratch/genomic_med/gm_dep/samin1/orion/workdirs/%s/logs/%s_%s_%s_%s.Rout", myproject, as.character(args[1]), as.character(args[2]), make.names(format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")), batch_uuid)
if (!isTRUE(file.info(dirname(outFile))$isdir)) dir.create(dirname(outFile), recursive=TRUE)

sink(outFile, split = TRUE)

print(sprintf("Project name is %s and batch ID is %s", myproject, batch_uuid))
print(sprintf("Work dir was %s", getwd()))

myworkdir <- sprintf("/scratch/genomic_med/gm_dep/samin1/orion/workdirs/%s", myproject)
#if (!isTRUE(file.info(myworkdir)$isdir)) dir.create(myworkdir, recursive=TRUE)

setwd(myworkdir)
print(sprintf("Work dir is now set to %s", getwd()))
print(sprintf("logs at: %s",outFile))

Sys.getenv("MYWORKDIR")
Sys.setenv(MYWORKDIR = sprintf("%s", getwd()))
print(sprintf("MYWORKDIR env variable set to: %s", Sys.getenv("MYWORKDIR")))
print(paste0("##########################################"))

if(!file.exists(sprintf("/scratch/genomic_med/gm_dep/samin1/orion/master/%s_master.RData", myproject))) stop(sprintf("Could not locate /scratch/genomic_med/gm_dep/samin1/orion/master/%s_master.RData", myproject))
load(sprintf("/scratch/genomic_med/gm_dep/samin1/orion/master/%s_master.RData", myproject))
ls()
dim(myinfo)
#myinfo[1:4,c(2,17,14,16)]
myinfo[1:4,c(1:4)]

print(paste0("##########################################"))

batchrun <- function(pairindex){
	pairindex <- as.numeric(pairindex)

	print(paste0("########## SAMPLE ENV VARIABLES ############"))
	print(sprintf("Work dir is set to %s", getwd()))
	lapply(c("ingest", "out", "logs" ,"processed"), function(mydir) dir.create(path=mydir))
	ls()

	print(myinfo[pairindex,c(1:4)])

	if(!file.exists(sprintf("/scratch/genomic_med/gm_dep/samin1/orion/master/%s.gtf", myproject))) stop(sprintf("Could not locate /scratch/genomic_med/gm_dep/samin1/orion/master/%s.gtf", myproject))
	MYGTF <- sprintf("/scratch/genomic_med/gm_dep/samin1/orion/master/%s.gtf", myproject)
	Sys.setenv(MYWORKDIR = sprintf("%s", getwd()))	
	Sys.setenv(MYGTF = sprintf("/scratch/genomic_med/gm_dep/samin1/orion/master/%s.gtf", myproject))
	Sys.setenv(MYKEY = "/scratch/genomic_med/gm_dep/samin1/orion/master/mycgkey")
	print(Sys.getenv("MYWORKDIR"))
	print(Sys.getenv("MYGTF"))
	print(Sys.getenv("MYKEY"))

	library(tools)
	if(!file.exists(sprintf("ingest/%s/%s", myinfo[pairindex,2], myinfo[pairindex,3]))) stop(sprintf("Could not locate ingest/%s/%s", myinfo[pairindex,2], myinfo[pairindex,3]))
#	command2 <- sprintf("samtools sort -n ingest/%s/%s ingest/%s/%s.sorted", myinfo[pairindex,2], myinfo[pairindex,3], myinfo[pairindex,2], basename(file_path_sans_ext(as.character(myinfo[pairindex,3]))))
	
	command2 <- sprintf("samtools sort -@ 3 -m 2G -n ingest/%s/%s ingest/%s/%s.sorted", myinfo[pairindex,2], myinfo[pairindex,3], myinfo[pairindex,2], basename(file_path_sans_ext(as.character(myinfo[pairindex,3]))))

	print(paste0("########## command2 command string ############"))
	print(paste0(command2))
	print(sprintf("######### RUNNING COMMAND2 ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
	system(command2,wait=TRUE)
	print(sprintf("######### COMMAND2 ENDED ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))	

	command3 <- sprintf("samtools view -F 4 ingest/%s/%s.sorted.bam | htseq-count -m intersection-nonempty --stranded=no --idattr gene_id -r name - ${MYGTF} > out/%s_%s_%s.htseq.log 2>&1 > out/%s_%s_%s.htseq.cnt", myinfo[pairindex,2], basename(file_path_sans_ext(as.character(myinfo[pairindex,3]))), pairindex, basename(file_path_sans_ext(as.character(myinfo[pairindex,3]))), make.names(format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")), pairindex, basename(file_path_sans_ext(as.character(myinfo[pairindex,3]))), make.names(format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
	print(paste0("########## command3 command string ############"))
	print(paste0(command3))
	print(sprintf("######### RUNNING COMMAND3 ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
	system(command3,wait=TRUE)
	print(sprintf("######### COMMAND3 ENDED ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))

	command4 <- sprintf("mv ingest/%s.gto processed/ & mv ingest/%s processed/ &", myinfo[pairindex,2], myinfo[pairindex,2])
	print(paste0("########## command4 command string ############"))
	print(paste0(command4))
	print(sprintf("######### RUNNING COMMAND4 ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
	system(command4,wait=TRUE)
	print(sprintf("######### COMMAND4 ENDED ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))

	print(paste0("########## Total Samples Processed ############"))
	print(system(sprintf("ls -alh %s/out/*.cnt | awk '{if($5 != 0) print $0}' | grep -Eo \"PCAWG.*v1\" | sort | uniq | wc -l", getwd())))
	print(sprintf("############################ DONE %03d/%s ######################################", pairindex, nrow(myinfo)))
}

#gtagent_batchrun(args[1] args[2])
library(parallel)
mclapply(seq(as.numeric(args[1]),as.numeric(args[2]),1),function(k) batchrun(k),mc.preschedule=FALSE,mc.cores=1)
Sys.sleep(5)
print(paste0("JOB DONE: ",batch_uuid," At: ",format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")," | LOGS at: ",outFile))
sink()
#end
