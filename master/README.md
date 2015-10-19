

### Annotations:

PCAWG14 extended annotations are available at https://www.synapse.org/#!Synapse:syn4590784 with direct link: https://www.synapse.org/#!Synapse:syn4590793

Original annotation file, `rnaseq.gc19_extNc.bed` was compiled by Morten Nielsen (Muhlig) of  Prof. Jakob Pedersen Lab, Aarhus University, Denmark

```
rnaseq.gc19_extNc.bed md5: 6a9937bf2a53a375765cebf3980690d1 ;
rnaseq.gc19_extNc.gtf md5: 48245eeff794b9e686466a80e66905c9 ;
```

Rename `rnaseq.gc19_extNc.gtf ` to `pcawg14.gtf`


### GNOS key:

Replace `mycgkey` with your GNOS key to download protected RNA-seq data from GNOS servers. This sets `MYKEY` environment variable referring to path to `mycgkey`.

You will also need to update absolute path to cgkey in file, `mumbai_cmds.txt`. You may use environment variable `$MYKEY` instead of absolute path.

### Library manifest:

Two files, `pcawg14_master.RData` and `mumbai_cmds.txt` were made using RNA-seq library manifest downloaded using `gnosquery` command as follows:

```
gnosquery -s https://gtrepo-ebi.annailabs.com  -a -o manifest_ebi.xml "center_name=MSKCC&state=live&filename=*STAR*"
```

##### Pending manifests:

```
gnosquery -s https://gtrepo-bsc.annailabs.com  -a -o manifest_bsc.xml "center_name=MSKCC&state=live&filename=*STAR*"
gnosquery -s https://gtrepo-osdc-icgc.annailabs.com -a -o manifest_osdc.xml "center_name=MSKCC&state=live&filename=*STAR*"
```

*   `gnosquery` using Gene Torrent installed on local (non-HPC) Ubuntu 14.04 machine with version: GeneTorrent gtdownload release 4.0.0 (SCM REV: git ref: 87efe569e7):

*   `gtdownload` using GeneTorrent version on HPC Nautilus: GeneTorrent gtdownload release 3.8.7, git ref: 476d7bd449, build: 207, type: tar.gz, target: CentOS6.4.x86_64

*   GeneTorrent download link: https://annaisystems.zendesk.com/hc/en-us/articles/205449647-GeneTorrent-Downloads


