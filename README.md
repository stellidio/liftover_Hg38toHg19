# liftover_Hg38toHg19
Change specific chromosome positions from hg38 to hg19 assembly

Liftover is the process of conversion of the genomic ranges from one reference genome to another.

https://genome.ucsc.edu/cgi-bin/hgLiftOver

At this project the initial genomics positions are in GRCh38(hg38) and they are converted to GRCh37(hg19).

For this purpose Bioconductor package "rtracklayer" is used. http://bioconductor.org/packages/release/bioc/html/rtracklayer.html

# Prerequisites

Download a chain file, according to the conversion that is made. In these case "hg38ToHg19.over.chain" is used.
(Download: http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/)
