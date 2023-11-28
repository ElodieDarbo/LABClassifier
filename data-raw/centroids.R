## code to prepare `centroids` dataset goes here

TCGA.rsem <- download.file(url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz", destfile = "TCGA_RSEM.tab")
TCGA.rsem <- read.table(TCGA.rsem,sep="\t",header=T,row=1)
colnames(TCGA.rsem) <- gsub("-",".",colnames(TCGA.rsem))
TCGA.rsem <- TCGA.rsem[,colnames(TCGA.rsem)%in%patients.keep]

TCGA.noGrey <- lab.classifier(as.matrix(TCGA.rsem),prefix="TCGA_RSEM")

clean.annots <- TCGA.noGrey$LABclass[TCGA.noGrey$LABclass$pred!="unknown",]
row.names(clean.annots) <- clean.annots$ID
clean.annots$ss.split <- "sensory"
clean.annots$ss.split[clean.annots$pred=="Basal"] <- "secretory"
clean.annots$ss.la <- clean.annots$pred
clean.annots$ss.la[clean.annots$pred=="Basal"] <- NA

clean.data <- TCGA.rsem[,clean.annots$ID]

centroides <- make.LABclassifier(clean.data,clean.annots)

usethis::use_data(centroides, overwrite = TRUE)
