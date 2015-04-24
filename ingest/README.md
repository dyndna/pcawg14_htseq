## PCAWG-3 BAM PATH:

Download bams from CGHub/Annai GNOS to `<path_to_git_repo>/ingest/` using `gtdownload` command. More at *PCAWG Researcher's Guide*, https://wiki.oicr.on.ca/x/BKVCAw

Example gtdownload command running against cgHub for a single analysis uuid:

```
cd <path_to_git_repo>/ingest/
gtdownload -c ~/cghub.key -v -d <analysis_id>
gtdownload -c ~/cghub.key -v -d 6e858862-703f-452d-889c-5d8cd03793b9
```

>Where <analysis_id> can be fetched from `cgquery` or CGHub browser URL: 
<https://browser.cghub.ucsc.edu/search/?study=(%22PCAWG+2.0%22)&library_strategy=(RNA-Seq)&refassem_short_name=(GRCh37)&state=(live)> 

>PS: Depending on update in PCAWG data freeze version, search criteria for CGHub query may differ in future.

Downloaded bams should be stored in following order:

```
<path_to_git_repo>/ingest/<analysis_id>/<sample.bam>
<path_to_git_repo>/ingest/<analysis_id>/<sample.bam.bai>
```

>e.g., `ls /scratch/pcawg14_htseq/ingest/6e858862-703f-452d-889c-5d8cd03793b9` should show following output:

>PCAWG.4588e050-e15f-476f-9065-ba747bdd3736.STAR.v1.bam
>PCAWG.4588e050-e15f-476f-9065-ba747bdd3736.STAR.v1.bam.bai







