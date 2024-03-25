## code to prepare TCGA dataset

TCGA.rsem <- download.file(url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz", destfile = "TCGA_RSEM.tab")
TCGA.rsem <- read.table( "TCGA_RSEM.tab",sep="\t",header=T,row=1)
colnames(TCGA.rsem) <- gsub("-",".",colnames(TCGA.rsem))
TCGA.rsem <- TCGA.rsem[,colnames(TCGA.rsem)%in%patients.keep]
TCGA.rsem <- TCGA.rsem[row.names(TCGA.rsem)%in%c(toupper(c("esr1","ar","foxa1","tox3","spdef","gata3","myb","msx2","tfap2b","esrrg")),toupper(c("foxc1","bcl11a","elf5","klf5","vgll1","nfib","id4","sox10","en1")),c("S100A8","IL8", "DKK1","HSD17B2",  "PERP","SOX11","AKR1B10","RARRES1","KRT7","KYNU", "PSAT1","KMO","CLCA2","SERHL2","CLDN8","UGT2B28"),c("ESR1","CA12","BCL2","GFRA1","GREB1","FAM134B","IGF1R","NPY1R","ANXA9","SERPINA5","SCCPDH","IRS1","ABAT","SERPINA3","TFF1","AGR3","NAT1","GATA3","ERBB4","MTL5"),RA_genes,"AR","ERBB2"),]

usethis::use_data(TCGA.rsem, overwrite = TRUE)
