## htseq on flowr

### 10-2015 | Samir Amin @sbamin

require(flowr)
if(FALSE)
	stop("Requirements not met; unable to load library: flowr")
require(params)
if(FALSE)
	stop("Requirements not met; unable to load library: params")
require(ngsflows)
if(FALSE)
	stop("Requirements not met; unable to load library: ngsflows")
require(tools)
if(FALSE)
	stop("Requirements not met; unable to load library: tools")

## create a flow object
get_flow <- function(s, out, myflowname, myoutdir){
	sampout = subset(out, samplename == s)
	fqs = sampout$file
	dataurl = sampout$url
	sample_name = s
	flowname = myflowname
	myreports = myoutdir
	
	## create a flowmat
	load_opts("pcawg14_htseq.conf")
	flowmat = pcawg14_htseq(fqs = fqs, sample_name = sample_name, dataurl = dataurl)
	## read the flowdef		
	deffile = file.path("pcawg14_htseq.csv")
	def = as.flowdef(deffile)
	## saving plot as a pdf file
	plot_flow(def, pdf = TRUE, pdffile = file.path(myreports, sprintf("pcawg14_htseq_%s.pdf", flowname)))
	
	## make a plot to check
	plot_flow(def)
	
	## create a flow object
	fobj = to_flow(flowmat, def, flowname = flowname)
}

#' @param fqs a character vector of path to bam file for a sample
#' @param sample_name name of this sample (a character vector, length one.)
pcawg14_htseq <-function(fqs, sample_name, dataurl,
													my_gtdownload = get_opts("gtdownload_exe"),
													my_cgkey = get_opts("cgkey_path"),
													my_samtools = get_opts("samtools_exe"),
													my_htseq = get_opts("htseq_exe"),
													my_gtf=get_opts("gtf_path"),
													my_sort_mem = get_opts("sort_mem"),
													my_sort_nt = get_opts("sort_nt"),
													mycode_path=get_opts("mycode_path")
													){
	
	## we also want to make sure that all options are recovered. for example get_opts("NOTHING"), returns NULL
	## but does not fail.
	
	## this function will check ALL the arguments of this function and make sure none of them are null.
	check_args()
	
	time_tag <- "\"$(date +%d%b%y_%H%M%S%Z)\""

	cmd_download <- sprintf("gtdownload -c %s -v -d %s", 
		my_cgkey, dataurl)

	cmd_sort <- sprintf("samtools sort -@ %s -m %sG -n %s/%s %s/%s.sorted", 
											my_sort_nt, my_sort_mem,
											sample_name, fqs, sample_name, basename(file_path_sans_ext(as.character(fqs)))) 

	cmd_count <- sprintf("samtools view -F 4 %s/%s.sorted.bam | htseq-count -m intersection-nonempty --stranded=no --idattr gene_id -r name - %s > %s/%s_%s.htseq.cnt", 
											sample_name, basename(file_path_sans_ext(as.character(fqs))),
											my_gtf,
											sample_name, basename(file_path_sans_ext(as.character(fqs))), time_tag)

	cmd.list = list(download = cmd_download, 
									sort = cmd_sort,
									count = cmd_count)
	
	my.flowmat <- to_flowmat(x = cmd.list, samplename = sample_name)
	return(my.flowmat)
}

#END
