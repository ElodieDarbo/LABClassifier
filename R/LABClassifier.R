rotate.45 <- function(score1,score2){
  a<-c(0.7,-0.7)
  b<-c(0.7,0.7)
  R<-rbind(a,b)
  v<-rbind(score1,score2)
  lp.tf<-R %*% v
  return(lp.tf)
}

plot.splits <- function(filename,LP.dat,LA.dat){
  mean.hsc <- mean.msc <- mean.lsc <- mean.asc <- classif <- NULL
  axis.min <- min(LP.dat[,1:2])
  axis.max <- max(LP.dat[,1:2])
  g.sp <- ggplot(LP.dat,aes(x=mean.hsc,y=mean.msc)) + theme_bw() +
    geom_point(aes(color=classif),size=0.7) +
    geom_abline(slope = 1,intercept = 0, linetype="dashed",size=0.4) +
    scale_colour_manual(values=c("sensory"="black","secretory"="red")) +
    labs(x="Sensory score",y="Secretory score") +
    coord_cartesian(xlim = c(axis.min,axis.max),ylim=c(axis.min,axis.max)) +
    ggtitle("All tumours by LP-split TFs") +
    theme(legend.position = "top",aspect.ratio = 1,legend.title=element_blank())

  axis.min <- min(LA.dat[,3:4])
  axis.max <- max(LA.dat[,3:4])

  g.la <- ggplot(LA.dat,aes(x=mean.lsc,y=mean.asc)) + theme_bw() +
    geom_point(aes(color=classif,fill=classif),size=0.7,shape=21) +
    geom_abline(slope = 1,intercept = 0, linetype="dashed",size=0.4) +
    scale_colour_manual(values=c("Luminal"="darkblue","MA"="mediumorchid")) +
    scale_fill_manual(values=c("Luminal"="darkblue","MA"="pink")) +
    labs(x="Luminal score",y="Apocrine score") +
    coord_cartesian(xlim = c(axis.min,axis.max),ylim=c(axis.min,axis.max)) +
    ggtitle("Sensory tumours by ER metagene") +
    theme(legend.position = "top",aspect.ratio = 1,legend.title=element_blank())

  multiplot(g.sp,g.la,cols=2)
  export.plot(filename,width=8,height=5)
}




ssGSEA.classif <- function(num, filename,sensor.genes,secretor.genes,asc.genes,lsc.genes) {
  #Use cell identity TFs to find the lumapo/basal (luminal progenitor) split
  sensor <- num[row.names(num)%in%sensor.genes,]
  if (sum(row.names(num)%in%sensor.genes)<length(sensor.genes)){
    message("Sensor genes")
    message("WARNING: genes ",paste(sensor.genes[!sensor.genes%in%row.names(num)],collapse=", ")," are not in your input data")
  }
  secretor <- num[row.names(num)%in%secretor.genes,]
  if (sum(row.names(num)%in%secretor.genes)<length(secretor.genes)){
    message("Secretor genes")
    message("WARNING: genes ",paste(secretor.genes[!secretor.genes%in%row.names(num)],collapse=", ")," are not in your input data")
  }
  # Genes for LA split
  #la.list <- c("ESR1", "CA12", "BCL2", "GFRA1", "GREB1", "FAM134B", "IGF1R", "NPY1R", "ANXA9", "SERPINA5", "SCCPDH", "IRS1", "ABAT", "SERPINA3", "MTL5", "IL8", "LIMCH1", "DKK1", "HSD17B2", "PERP", "TFAP2B", "SOX11", "AKR1B10", "RARRES1", "KRT7", "KYNU", "PSAT1", "PAPSS2", "KMO", "CLCA2")
  if (sum(row.names(num)%in%c(asc.genes,lsc.genes))<length(c(asc.genes,lsc.genes))){
    message("Luminal/Apocrine split genes")
    message("WARNING: genes ",paste(c(asc.genes,lsc.genes)[!c(asc.genes,lsc.genes)%in%row.names(num)],collapse=", ")," are not in your input data")
  }

  v <- F
  if (ncol(num)==1){
    v <- T
  }

  gene.split.list <- list(mean.hsc=sensor.genes,mean.msc=secretor.genes,mean.asc=asc.genes,mean.lsc=lsc.genes,RA_activ=RA_genes)
  message("Computing scores ...")
  all.scores <- gsva(as.matrix(num),gset.idx.list = gene.split.list,method="ssgsea", ssgsea.norm = FALSE, verbose = FALSE)/1000


  all.scores <- as.data.frame(t(all.scores))

  LP.dat <- all.scores

  lp.tf <- rotate.45(LP.dat$mean.hsc,LP.dat$mean.msc)
  lpsplit<-lp.tf[1,]
  names(lpsplit) <- row.names(all.scores)

  ##############################

  if (v){
    is.basal <- lpsplit<0
    if (is.basal){
      LABclass <- data.frame(row.names=colnames(num),score=lpsplit,pred="Basal",AR_activity=all.scores$RA_activ)
    } else {
      mal <- rotate.45(all.scores$mean.lsc,all.scores$mean.asc)
      lasplit<-mal[1,]
      if (lasplit>=0){
        LABclass <- data.frame(row.names=colnames(num),score=lasplit,pred="Luminal",AR_activity=all.scores$RA_activ)
      } else {
        LABclass <- data.frame(row.names=colnames(num),score=lasplit,pred="MA",AR_activity=all.scores$RA_activ)
      }
    }
  }
  else {
    #sensory cell tumours
    hormone <- which(lpsplit>=0)

    LA.dat <- all.scores[hormone,]
    #Rotate backwards by 45 degree to put the variation on one axis
    mal <- rotate.45(LA.dat$mean.lsc,LA.dat$mean.asc)
    lasplit<-mal[1,]
    names(lasplit) <- row.names(LA.dat)

    LP.dat$classif <- ifelse(lpsplit>=0,"sensory","secretory")
    LA.dat$classif <- ifelse(lasplit>=0,"Luminal","MA")

    message("Plotting classification split steps")
    plot.splits(filename,LP.dat,LA.dat)
    clean <- lasplit
    clean.lp <- lpsplit

    LABlum<-clean[which(clean>=0)]
    LABapo<-clean[which(clean<0)]
    LABbas<-clean.lp[which(clean.lp<0)]

    LUM<-rep("Luminal",length(LABlum))
    APO<-rep("MA",length(LABapo))
    BAS<-rep("Basal",length(LABbas))
    LABclass <- data.frame(row.names=c(names(LABlum),names(LABapo),names(LABbas)),score=c(as.numeric(LABlum),as.numeric(LABapo),as.numeric(LABbas)),pred=c(LUM,APO,BAS),AR_activity=all.scores[c(names(LABlum),names(LABapo),names(LABbas)),"RA_activ"])

  }

  return(LABclass)

}


prepare.data <- function(data,id.type,log2T,raw.counts,PAM50){
  v <- F
  if (is.vector(data)){
    if (is.null(names(data))) stop("ERROR: You submitted a vector that is not named.")
    data <- data.frame(row.names=names(data),sample=data)
    v <- T
    if (PAM50){
      message("WARNING: PAM50 can't be applied on a simgle sample. PAM50 is set to FALSE")
      PAM50 <- F
    }
  }
  switch(id.type,
         ensEMBL={
           ensg <- row.names(data)
           ensg <- sub("[.][0-9]+$","",ensg)
           m1 <- match(ensg,gene.length$v86)
           m2 <- match(ensg,gene.length$v75)
           symbols1 <- gene.length$SYMBOL[m1]
           symbols2 <- gene.length$SYMBOL[m2]
           symbols1[is.na(symbols1)] <- symbols2[is.na(symbols1)]
           data <- data[!is.na(symbols1),,drop=F]
           symbols1 <- na.omit(symbols1)
           dup <- duplicated(symbols1)
           data <- data[!dup,,drop=F]
           row.names(data) <- symbols1[!dup]
           message("You have submitted ",length(ensg)," EnsEMBL genes.")
           message("We found ",nrow(data)," corresponding gene symbols.")
         },
         entrezID={
           id <- row.names(data)
           m1 <- match(id,gene.length$entrezid)
           m2 <- match(id,gene.length$entrez75)
           symbols1 <- gene.length$SYMBOL[m1]
           symbols2 <- gene.length$SYMBOL[m2]
           symbols1[is.na(symbols1)] <- symbols2[is.na(symbols1)]
           data <- data[!is.na(symbols1),,drop=F]
           symbols1 <- na.omit(symbols1)
           dup <- duplicated(symbols1)
           data <- data[!dup,,drop=F]
           row.names(data) <- symbols1[!dup]
           message("You have submitted ",length(id)," ENTREZ genes.")
           message("We found ",nrow(data)," corresponding gene symbols.")
         },
         SYMBOL={
           # Do nothing
         },
         stop(paste(id.type, "not recognized"))
  )
  if (log2T){
    data <- log2(data+1)
  }
  if (raw.counts){
    message("Transforming RNA-seq raw counts into log2(TPM+1)")
    m <- match(row.names(data),gene.length$SYMBOL)
    ldata <- data[!is.na(m),,drop=F]
    gene.length <- gene.length[na.omit(m),,drop=F]
    ldata <- (ldata / (gene.length$effective_length+1))*1000
    TPM <- apply(ldata,2,function(x){
      x <- (x/sum(x))*10^6
    })
    log2T <- F
    data <- log2(TPM+1)
  }
  return(data)
}


#' \code{PAMgenefu} uses \code{molecular.subtyping} from genefu package.
#' Please cite the corresponding paper (See reference).
#' @param data A data.frame containing gene (gene symbols as rows) expression
#' from the samples (sample IDs as columns) to be classified or a  vector
#' containing gene expression named with gene SYMBOL
#' @references
#'
#' genefu:  Deena M.A. Gendoo et al. (2021). genefu: Computation of Gene
#' Expression-Based Signatures in Breast Cancer. R package version 2.26.0.
#' http://www.pmgenomics.ca/bhklab/software/genefu
#' @importFrom genefu molecular.subtyping
#' @return \code{PAMgenefu} Returns a data.frame with 1 column
#' \item{PAM50}{PAM50 classification: LumA and LumB (Luminal A and B), Basal, Her2.}


PAMgenefu <- function(data) {
  annots <- gene.length[gene.length$entrezid%in%pam50.robust$centroids.map$EntrezGene.ID,c("SYMBOL","entrezid")]
  annots <- merge(annots, data, by.x = "SYMBOL",by.y=0)
  annots <- unique(annots)
  row.names(annots) <- annots$SYMBOL

  colnames(annots) <- sub("SYMBOL","Gene.Symbol",colnames(annots))
  dat <- annots[,-c(1:2),drop=F]
  dat <- t(dat)

  annots <- annots[,1:2]

  # Subtype with PAM50
  output <- molecular.subtyping(
    sbt.model="pam50",
    data=dat,
    annot=annots,
    do.mapping=FALSE)

  pam50 <- output$subtype
  pam50 <- data.frame(row.names=names(pam50),PAM50=as.character(pam50))

  return(pam50)

}


#' \code{expression.dotplot} helps in visualising the sample classification
#'
#' @param data A data.frame containing gene (rows) expression from the samples (columns)
#' to be classified
#' @param predictions A data.frame containing the predictions obtained with
#' \code{classify_splits} function;
#' @param g1 A gene symbol to be visualized on x-axis;
#' @param g2 A gene symbol to be visualized on y-axis;
#' @param PAM50 If TRUE, the colors correspond to PAM50 classification
#'
#' @import ggplot2
#' @return \code{classify_splits} Returns a ggplot object


expression.dotplot <- function(data, predictions,g1,g2,PAM50=F){
  prediction <- NULL
  colnames(data) <- gsub("-",".",colnames(data))
  classes <- c("Basal","Luminal","MA")
  predictions$pred <- factor(as.vector(predictions$pred),levels=c("Luminal","Basal","MA"))
  if (PAM50){
    classes <- c("Basal","LumA","LumB","Her2")
    predictions <- predictions[predictions$PAM50%in%classes,]
    predictions$PAM50 <- factor(as.vector(predictions$PAM50),levels=c("LumA","LumB","Basal","Her2"))
    data <- data[,predictions$Row.names]
  }
  data <- as.data.frame(t(data))
  data$patientID <- row.names(data)
  #data <- data[!grepl("^TCGA",data$patientID),]
  data <- data[,c("patientID",g1,g2)]
  colnames(data) <- c("patientID","g1","g2")
  data$prediction <- predictions$pred[match(data$patientID,row.names(predictions))]
  if (PAM50){
    data$prediction <- predictions$PAM50[match(data$patientID,row.names(predictions))]
  }
  data <- data[order(data$prediction),]
  color.LAB <- c(LumA="darkblue",Her2="deeppink",LumB="lightblue",Basal="red",unclassified="grey",Luminal="darkblue",MA="mediumorchid")
  fill.LAB <- c(Her2="pink",MA="pink",LumA="darkblue",LumB="lightblue",Basal="red",unclassified="grey",Luminal="darkblue")
  color.LAB<- color.LAB[names(color.LAB)%in%classes]
  fill.LAB<- fill.LAB[names(fill.LAB)%in%classes]
  g <- ggplot(data,aes(x=g1,y=g2)) + theme_bw() +
    labs(x=paste(g1,"expression"),y=paste(g2,"expression")) +
    scale_colour_manual(values=color.LAB) +
    scale_fill_manual(values=fill.LAB) +
    theme(legend.position="top",
          legend.title=element_blank(),
          aspect.ratio = 1
    ) +
    geom_point(aes(color=prediction,fill=prediction),shape=21)

  return(g)
}


#' \code{LABclassifier} computes splits according expression
#' of TFs defining the splits enrichment according to ssGSEA
#' \code{LABclassifier} helps in visualising the sample classification
#'
#' @param data A data.frame containing gene (gene symbols as rows) expression
#' from the samples (sample IDs as columns) to be classified or a  vector
#' containing gene expression named with gene SYMBOL
#' @param dir.path A character string defining the path to where the results
#' will be stored. Default: current directory.
#' @param prefix A character string: prefix of the outputs. Default: myClassif
#' @param raw.counts A logical. Set it to TRUE if you are using raw sequencing
#' data, then they will be transformed in TPM. Default: FALSE;
#' @param log2T A logical. If true, data are log2 transformed. Default: FALSE;
#' @param id.type A character string. Gene identifier type in the input matrix. Possible values are
#' ensEMBL, entrezID, SYMBOL (default)
#' @param PAM50 A logical. Compute PAM50 classification using \code{molecular.subtyping}
#' function from genefu package. If used, please cite the corresponding paper (See reference).
#' @param plot A logical. Set it to TRUE to display classification
#' diagnostic plots. Default: FALSE
#' @param sensor.genes,secretor.genes,asc.genes,lsc.genes Vectors: Gene list to compute splits.
#' If no gene list is given, pre-defined ones are used. Default: NULL
#' @references
#' genefu:  Deena M.A. Gendoo et al. (2021). genefu: Computation of Gene
#' Expression-Based Signatures in Breast Cancer. R package version 2.26.0.
#' http://www.pmgenomics.ca/bhklab/software/genefu
#' @import ggplot2
#' @importFrom utils write.table
#' @export
#' @import grDevices
#' @import pheatmap
#' @import viridis
#' @importFrom GSVA gsva
#' @return A data.frame containing LAB classification
#' @examples
#' data(TCGA.rsem)
#' LABclassifier(TCGA.rsem,plot=TRUE)


LABclassifier <- function(data,dir.path=".",prefix="myClassif",raw.counts=F,log2T=F,id.type="SYMBOL",PAM50=F,plot=T,sensor.genes=NULL,secretor.genes=NULL,asc.genes=NULL,lsc.genes=NULL,colorBlind=F){
  AR_activity <- pred <- NULL
  message("Creating an Output folder in working directory: ",dir.path)
  dir.create(file.path(dir.path,"Output"), recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(dir.path,"Output",prefix)
  message("Classifying according to LAB classifier")
  if (is.null(sensor.genes)){
    sensor.genes <- toupper(c("esr1","ar","foxa1","tox3","spdef","gata3","myb","msx2","tfap2b","esrrg"))
  }
  if (is.null(secretor.genes)){
    secretor.genes <- toupper(c("foxc1","bcl11a","elf5","klf5","vgll1","nfib","id4","sox10","en1"))
  }
  if (is.null(asc.genes)){
    asc.genes <- c("S100A8","IL8", "DKK1","HSD17B2",  "PERP","SOX11","AKR1B10","RARRES1","KRT7","KYNU",
                   "PSAT1","KMO","CLCA2","SERHL2","CLDN8","UGT2B28")
  }
  if (is.null(lsc.genes)){
    lsc.genes <- c("ESR1","CA12","BCL2","GFRA1","GREB1","FAM134B","IGF1R","NPY1R","ANXA9","SERPINA5",
                   "SCCPDH","IRS1","ABAT","SERPINA3","TFF1","NAT1","ERBB4","MTL5")
  }
  data <- prepare.data(data,id.type,log2T,raw.counts,PAM50)
  if (ncol(data)==1){
    PAM50 <- F
  }
  LABclass <- ssGSEA.classif(data,prefix,sensor.genes,secretor.genes,asc.genes,lsc.genes)

  if (PAM50){
    message("Computing PAM50 classification")
    pam50.pred <- PAMgenefu(data)
    LABclass <- merge(LABclass,pam50.pred,by=0)
    row.names(LABclass) <- LABclass$Row.names
  }
  pred.final <- LABclass

  if (plot) {
    w <- 18
    h <- 4
    if (ncol(data)>1){
      annot.col <- pred.final[,c("pred","AR_activity")]
      colnames(annot.col)[1] <- "LABClassif"
      annot.colors <- list(LABClassif=c(MA="pink",Basal="red",Luminal="darkblue"), AR_activity = inferno(10, begin = 0, end = 0.8))
      g3 <- expression.dotplot(data, pred.final,"ESR1","FOXA1")
      g5 <- expression.dotplot(data, pred.final,"AR","FOXA1")
      g6 <- expression.dotplot(data, pred.final,"AR","ERBB2")
      g4 <- ggplot(pred.final,aes(x=AR_activity)) + theme_bw() +
        geom_density(aes(color=pred)) + geom_rug(aes(color=pred),length = unit(0.1, "npc")) +
        scale_color_manual(values=c(MA="pink",Basal="red",Luminal="darkblue")) +
        theme(aspect.ratio = 3/4,legend.title=element_blank(),legend.position = "top")

      if (PAM50){
        annot.col$PAM50 <- pred.final$PAM50
        annot.colors$PAM50 <- c(Her2="pink",Basal="red",LumA="darkblue",LumB="lightblue",Normal="grey")
        gpam1 <- expression.dotplot(data, pred.final,"ESR1","FOXA1",PAM50=T)
        gpam3 <- expression.dotplot(data, pred.final,"AR","FOXA1",PAM50=T)
        gpam4 <- expression.dotplot(data, pred.final,"AR","ERBB2",PAM50=T)
        gpam2 <- ggplot(pred.final,aes(x=AR_activity)) + theme_bw() +
          geom_density(aes(color=PAM50)) + geom_rug(aes(color=PAM50),length = unit(0.1, "npc")) +
          scale_color_manual(values=c(Her2="pink",Basal="red",LumA="darkblue",LumB="lightblue")) +
          theme(aspect.ratio = 3/4,legend.position = "top",legend.title=element_blank())
        print(multiplot(g3,gpam1,g5,gpam3,g6,gpam4,g4,gpam2,cols=4))
        h <- 8
        export.plot(paste0(prefix,"_LAB_predictions"),width=w,height=h)
        confusion <- as.data.frame.matrix(table(pred.final$pred,pred.final$PAM50))
        pheatmap(round(confusion[,colnames(confusion)%in%c("Basal","Normal","LumA","LumB","Her2")]),
                 cluster_rows = F,
                 cluster_cols = F,
                 color=colorRampPalette(c("white",colours()[72]))(100),
                 display_numbers = T,
                 number_format = "%.0f"
        )
        export.plot(paste0(prefix,"_LAB_PAM50_confusion"),width=8,height=6)

      } else {
       multiplot(g3,g5,g6,g4,cols=4)
       export.plot(paste0(prefix,"_LAB_predictions"),width=w,height=h)
      }
      data <- data[row.names(data)%in%c(sensor.genes,secretor.genes,asc.genes,lsc.genes),]
      data <- data - apply(data,1,mean)
      thr <- max(abs(data))
      if (colorBlind){
        palette <- colorRampPalette(colours()[c(121,121,121,24,142,142,142)])(100)
      } else {
        palette <- colorRampPalette(c("green","green","green","black","red","red","red"))(100)
      }
      pheatmap(data,
               show_colnames = F,
               #show_rownames = F,
               border_color = NA,
               annotation_col = annot.col,
               annotation_colors = annot.colors,
               breaks=seq(-thr,thr,thr*2/100),
               color=palette,
               clustering_method = "ward.D",
               clustering_distance_cols = "correlation",
               clustering_distance_rows = "correlation"
               )
      export.plot(paste0(prefix,"_LAB_predictions_heatmap"),width=8,height=12)

    } else {
      message("A unique sample can't be plotted")
    }
  }

  pred.final$Row.names <- pred.final$score <-  NULL
  write.table(pred.final,paste0(prefix,"_LAB_predictions.tab"),sep="\t",quote=F)
  # RA_genes <- c("UGT2B28","SEC14L2","SERHL2","FKBP5","MYBPC1","AQP3","CLDN8","ZBTB16","IQGAP2","GGT1","ECHDC2","AZGP1","FASN","MCCC2","FMO5","SORD","SLC15A2","PIP","EAF2","CROT","ALCAM","RND1")
  return(pred.final)
}
