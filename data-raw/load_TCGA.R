## code to prepare TCGA dataset

TCGA.rsem <- download.file(url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz", destfile = "TCGA_RSEM.tab")
TCGA.rsem <- read.table( "TCGA_RSEM.tab",sep="\t",header=T,row=1)
colnames(TCGA.rsem) <- gsub("-",".",colnames(TCGA.rsem))
TCGA.rsem <- TCGA.rsem[,colnames(TCGA.rsem)%in%patients.keep]

usethis::use_data(TCGA.rsem, overwrite = TRUE)
