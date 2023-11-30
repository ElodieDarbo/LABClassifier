rotate.45 <- function(score1,score2){
  a<-c(0.7,-0.7)
  b<-c(0.7,0.7)
  R<-rbind(a,b)
    v<-rbind(score1,score2)
    lp.tf<-R %*% v
    return(lp.tf)
}

split.distrib <- function(lpsplit){
  mix.data <- mixgroup(lpsplit,breaks=100)
  mixpar <- mixparam(mu = c(-1,1), sigma = c(1), pi= c(0.4,0.6))
  fitdata <- mix(mix.data, mixpar=mixpar,dist="norm")
  m1<-fitdata$parameters[1,"mu"]
  m2<-fitdata$parameters[2,"mu"]
  cutoff.lp<-(m1+m2)/2
  scale.lp<-(m2-m1)/2
  clean.lp<-(lpsplit-cutoff.lp)/scale.lp

  #Express the 0.2 as a multiple of the sd
  lp.sd<-mean(fitdata$parameters$sigma)/scale.lp
  times.sd<-0.8/lp.sd

  return(list(fitdata=fitdata,clean=clean.lp,cutoff=cutoff.lp,times.sd=times.sd))
}

plot.splits <- function(filename,LP.dat,lp.tf,LP.split,LA.dat,mal,LA.split){
    name<-paste(filename,".pdf",sep="")
    pdf(file=name,width=11,height=5,pointsize=10,family="Helvetica")

    par(mfrow=c(2,5))

    plot(LP.dat, main="All tumours by LP-split TFs",xlab="Sensory score",ylab="Secretory score")
    mtext(side = 3, text = "A",adj=(-0.32),line=0.7, cex=1.75)

    plot(lp.tf[1,],lp.tf[2,], main="45 degree rotation of LP-split TFs",xlab="Luminal Progenitor Score",ylab="Luminal Progenitor noise")
    mtext(side = 3, text = "B",adj=(-0.32),line=0.7, cex=1.75)

    hist(lp.tf[1,],breaks=50, main="Luminal Progenitor Scores",xlab="Luminal Progenitor Score")
    mtext(side = 3, text = "C",adj=(-0.32),line=0.7, cex=1.75)

    plot(LP.split$fitdata,xlab="Luminal Progenitor Score", main="Distribution of LP Scores")
    mtext(side = 3, text = "D",adj=(-0.32),line=0.7, cex=1.75)
    abline(v=LP.split$cutoff)

    plot(density(LP.split$clean),main="Normalised LP Scores",xlab=NA)
    abline(v=c(-0.2,0.2),col="grey")
    abline(v=c(-1,1))
    mtext(side = 3, text = "E",adj=(-0.32),line=0.7, cex=1.75)
    mtext(side = 1, line= 3, text = "Luminal Progenitor score",cex=0.7)
    mtext(side = 1, line= 4, text = paste("grey zone = means +/-",round(LP.split$times.sd,2),"sd"),cex=0.7)

    plot(LA.dat, main="Sensory tumours by ER metagene",xlab="Luminal score",ylab="Apocrine score")
    mtext(side = 3, text = "F",adj=(-0.32),line=0.7, cex=1.75)

    plot(mal[1,],mal[2,], main="45 degree rotation of ER metagene",xlab="Luminal-Apocrine score",ylab="Luminal-Apocrine noise")
    mtext(side = 3, text = "G",adj=(-0.32),line=0.7, cex=1.75)

    hist(mal[1,],breaks=50, main="Luminal-Apocrine Scores",xlab="Luminal-Apocrine score")
    mtext(side = 3, text = "H",adj=(-0.32),line=0.7, cex=1.75)


    plot(LA.split$fitdata,xlab="Luminal-Apocrine Score", main="Distribution of LA Scores")
    abline(v=LA.split$cutoff)
    mtext(side = 3, text = "I",adj=(-0.32),line=0.7, cex=1.75)

    plot(density(LA.split$clean),main="Normalised LA Scores",xlab=NA)
    abline(v=c(-0.2,0.2),col="grey")
    abline(v=c(-1,1))
    mtext(side = 3, text = "J",adj=(-0.32),line=0.7, cex=1.75)
    mtext(side = 1, line= 3, text = "Luminal-Apocrine score",cex=0.7)
    mtext(side = 1, line= 4, text = paste("grey zone = means +/-",round(LA.split$times.sd,2),"sd"),cex=0.7)

    dev.off()
}


initial.classif <- function(num, filename) {
  #Use cell identity TFs to find the lumapo/basal (luminal progenitor) split
  sensor.genes <- toupper(c("esr1","ar","foxa1","tox3","spdef","gata3","myb","msx2","tfap2b","esrrg"))
  sensor <- num[row.names(num)%in%sensor.genes,]
  if (sum(row.names(num)%in%sensor.genes)<length(sensor.genes)){
    message("Sensor genes")
    message("WARNING: genes ",paste(sensor.genes[!sensor.genes%in%row.names(num)],collapse=", ")," are not in your input data")
  }
  secretor.genes <- toupper(c("foxc1","bcl11a","elf5","klf5","vgll1","nfib","id4","sox10","en1"))
  secretor <- num[row.names(num)%in%secretor.genes,]
  if (sum(row.names(num)%in%secretor.genes)<length(secretor.genes)){
    message("Secretor genes")
    message("WARNING: genes ",paste(secretor.genes[!secretor.genes%in%row.names(num)],collapse=", ")," are not in your input data")
  }
 # Genes for LA split
  la.list <- c("ESR1", "CA12", "BCL2", "GFRA1", "GREB1", "FAM134B", "IGF1R", "NPY1R", "ANXA9", "SERPINA5", "SCCPDH", "IRS1", "ABAT", "SERPINA3", "MTL5", "IL8", "LIMCH1", "DKK1", "HSD17B2", "PERP", "TFAP2B", "SOX11", "AKR1B10", "RARRES1", "KRT7", "KYNU", "PSAT1", "PAPSS2", "KMO", "CLCA2")
  if (sum(row.names(num)%in%la.list)<length(secretor.genes)){
    message("Luminal/Apocrine split genes")
    message("WARNING: genes ",paste(la.list[!la.list%in%row.names(num)],collapse=", ")," are not in your input data")
  }


  message("Compute sensory / secretory split")

  hsc<- (sensor - apply(sensor,1,mean))/apply(sensor,1,function(x){ IQR(x)/(qnorm(0.75) - qnorm(0.25)) })
  msc<- (secretor - apply(secretor,1,median))/apply(secretor,1,function(x){ IQR(x)/(qnorm(0.75) - qnorm(0.25)) }) #milk secreting cell

  mean.hsc<-apply(hsc,2,mean)
  mean.msc<-apply(msc,2,mean)

  LP.dat <- data.frame(mean.hsc,mean.msc)

  lp.tf <- rotate.45(mean.hsc,mean.msc)
  lpsplit<-lp.tf[1,]


  #Find the means of the distributions, put the cutoff half way between them
  #Convert to a more easily interpretable scale (set the means to -1 and 1)
  LP.split <- split.distrib(lpsplit)

  ##############################
  message("Compute MA / luminal split")

  #sensory cell tumours
  hormone <- which(LP.split$clean>=0)

  sct<-num[,hormone]

  lsc.genes <- la.list[1:15]
  lsc <- sct[row.names(sct)%in%lsc.genes,]
  lsc <- (lsc - apply(lsc,1,mean))/apply(lsc,1,function(x){ IQR(x)/(qnorm(0.75) - qnorm(0.25)) })

  asc.genes <- la.list[16:30]
  asc <- sct[row.names(sct)%in%asc.genes,]
  asc <-  (asc - apply(asc,1,mean))/apply(asc,1,function(x){ IQR(x)/(qnorm(0.75) - qnorm(0.25)) })

  mean.lsc<-apply(lsc,2,mean)
  mean.asc<-apply(asc,2,mean)


  LA.dat <- data.frame(mean.lsc,mean.asc)

  #Rotate backwards by 45 degree to put the variation on one axis
  mal <- rotate.45(mean.lsc,mean.asc)
  lasplit<-mal[1,]


  #Find the means of the distributions, put the cutoff half way between them
  LA.split <- split.distrib(lasplit)

  message("Plotting classification split steps")
  plot.splits(filename,LP.dat,lp.tf,LP.split,LA.dat,mal,LA.split)

  ##############################################################################

  #Make labels
  clean <- LA.split$clean
  clean.lp <- LP.split$clean


  LABlum<-clean[which(clean>=0)]
  LABapo<-clean[which(clean<0)]
  LABbas<-clean.lp[which(clean.lp<0)]

  LUM<-rep("Luminal",length(LABlum))
  APO<-rep("MA",length(LABapo))
  BAS<-rep("Basal",length(LABbas))

  LABclass <- data.frame(ID=c(names(LABlum),names(LABapo),names(LABbas)),score=c(as.numeric(LABlum),as.numeric(LABapo),as.numeric(LABbas)),pred=c(LUM,APO,BAS))


  LABclass$pred[LABclass$score >=(-0.2) & LABclass$score<=0.2] <- "unknown"

 save(LABclass,file=paste0(filename,".RData"))

 return(LABclass)

}


lab.classifier <- function(data,dir.path=".",prefix){
  setwd(dir.path)
  message("Creating an LABmodel folder in working directory: ",getwd())
  dir.create("LABmodel", recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path("LABmodel",prefix)
  message("Classifying according to LAB classifier")
  LABclass <- initial.classif(data,prefix)
  return(LABclass)
}


cit.dfAggregate  <-  function (data, partition, MARGIN = 2, fAggreg = mean.na){
  mean.na <- NULL
  cMARGIN <- setdiff(c(1, 2), MARGIN)
  n <- length(partition)
  N <- dim(data)[MARGIN]
  p <- dim(data)[cMARGIN]
  if (n != N)
    stop("Error - function cit.dfAggregate : size of partition doesn't correspond to data dimension")
  l <- split(1:N, partition)
  d <- data
  if (MARGIN == 2)
    d <- t(data)
  d <- matrix(sapply(l, function(i) if (length(i) == 1) {
    unlist(d[i, ])
  }
  else {
    apply(d[i, ], 2, fAggreg)
  }), ncol = p, byrow = TRUE)
  d <- as.data.frame(d)
  rownames(d) <- names(l)
  colnames(d) <- dimnames(data)[[cMARGIN]]
  if (MARGIN == 2)
    d <- as.data.frame(t(d))
  d
}

computeCentroids <- function(annots,DF,col_interest) {
  stopifnot(length(which(colnames(annots)==col_interest))==1)
  stopifnot(length(which(rownames(annots) %in% names(DF)))>0)
  cl <- rep(NA,ncol(DF))
  names(cl) <- names(DF)
  nom <- names(DF)
  for (i in nom) {
    cl[i]<- annots[i,col_interest]
  }
  dd <- t(scale(t(DF),scale=F))
  L <- list()
  L$centroids <- cit.dfAggregate(dd, cl , MARGIN = 2, fAggreg = mean)
  return(L)
}


make.model <- function(annot,df,spl){
  switch(spl,
         ss={
           col_of_interest <- "ss.split"
           genes <- toupper(c("esr1","ar","foxa1","tox3","spdef","gata3","myb","msx2","tfap2b","esrrg","foxc1","bcl11a","elf5","klf5","vgll1","nfib","id4","sox10","en1"))
         },
         la={
           col_of_interest <- "ss.la"
           genes <- c("ESR1", "CA12", "BCL2", "GFRA1", "GREB1", "FAM134B", "IGF1R", "NPY1R", "ANXA9", "SERPINA5", "SCCPDH", "IRS1", "ABAT", "SERPINA3", "MTL5", "IL8", "LIMCH1", "DKK1", "HSD17B2", "PERP", "TFAP2B", "SOX11", "AKR1B10", "RARRES1", "KRT7", "KYNU", "PSAT1", "PAPSS2", "KMO", "CLCA2")
         }
  )
  df <- df[row.names(df)%in%genes,]
  df.smote <- as.data.frame(t(df))
  df.smote$class <- factor(annot[colnames(df),col_of_interest])
  df.smote <- smote(df.smote, var="class",over_ratio = 0.5)
  smote.annot <- data.frame(pred=df.smote$class,tmp=df.smote$class,row.names=paste0("s",1:nrow(df.smote)))
  colnames(smote.annot)[1] <- col_of_interest
  df.smote$class <- NULL
  df.smote <- as.data.frame(t(df.smote))
  colnames(df.smote) <- row.names(smote.annot)
  centroids <- computeCentroids(annots=annot,DF=df,col_interest=col_of_interest)
  return(centroids)
}

#' \code{make.LABclassifier} compute centroids
#' @param data A data.frame containing expression values from classified samples;
#' @param annot A data.frame containing the Luminal/Basal/MA annotations
#' for samples
#' @export
#' @import grDevices
#' @import graphics
#' @import stats
#' @importFrom themis smote
#' @import mixdist
#' @return A object of class \code{LABclassifier} containing a list
#' \item{centroid.ss}{Centroids to classify according to the sensory/secretor split}
#' \item{centroid.la}{Centroids to classify according to the Luminal/Apocrine split}

make.LABclassifier <- function(data,annot){
  centroid.ss <- make.model(annot,data,spl="ss")
  centroid.la <- make.model(annot,data,spl="la")
  centroides <- list(centroid.ss=centroid.ss,centroid.la=centroid.la)
  class(centroides) <- "LABclassifier"
  return(centroides)
}
