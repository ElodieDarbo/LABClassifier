#' \code{classify_splits} classifies samples using centroids
#'
#' @param d2 A data.frame containing gene (rows) expression from the samples (columns)
#' to be classified;
#' @param L A list containing centroids produced by \code{computeCentroids}
#' function;
#' @param opt Determine which split to be computed: ss: sensory / secretory split;
#' la : Luminal / Apocrine split;
#'
#' @return \code{classify_splits} Returns a data.frame with columns
#' \item{scores}{Distance to centroids}
#' \code{"sensory","secretory"} if \code{opt} set to ss or
#' \code{"Luminal","MA"} if \code{opt} set to la
#' \item{predictions}{Centroid association}
#' \code{"pred.ss","pred.ss.strict"} if \code{opt} set to ss : sensory / Basal
#' \code{"pred.la","pred.la.strict"} if \code{opt} set to la : Luminal,
#'  MA: Molecular Apocrine
#' The strict version annotates the samples too far from any
#' centroid as unclassified


classify_splits <- function(d2,L,opt=NULL) {
  if (opt == "ss"){
    classes <- c("secretory","sensory")
  } else if (opt == "la"){
    classes <- c("Luminal","MA")
  }
  if (is.vector(d2)){
    N <- 1
    d2 <- na.omit(d2)
    d2 <- scale(d2)
    m <- match(row.names(L$centroids),row.names(d2))
    d2 <- cbind(d2=d2[na.omit(m)],L$centroids[!is.na(m),])
    nb.genes <- nrow(d2)
    n <- ncol(L$centroids)
    d2 <- t(d2)
    tdist <- as.matrix(cor(t(d2),method="spearman"))
    tdist <- as.data.frame(matrix(1-c(tdist[1,2],tdist[1,3]),ncol=2))
  } else {
    N <- ncol(d2)
    #d2 <- t(scale(t(d2),scale=F))
    d2=merge(d2,L$centroids,by="row.names")
    rownames(d2)=d2[,1]
    d2[,1]=NULL
    n <- ncol(L$centroids)
    d2 <- t(d2)
    tdist <- as.matrix(cor(t(d2),method="spearman"))
    tdist <- 1-as.data.frame(tdist[1:N,(N+1):(N+n)])
  }
  colnames(tdist) <- classes
  if (opt == "ss"){
    pred <- ifelse(tdist$sensory <= 1,"sensory","secretory")
    tdist$pred.ss <- pred
    tdist$pred.ss.strict <- pred
    tdist$pred.ss.strict[tdist$sensory>0.9 & tdist$sensory<1.1] <- "unclassified"
  } else if (opt == "la"){
    pred <- ifelse(tdist$Luminal <= 1,"Luminal","MA")
    tdist$pred.la <- pred
    tdist$pred.la.strict <- pred
    tdist$pred.la.strict[tdist$Luminal>0.9 & tdist$Luminal<1.1] <- "unclassified"
  }
  return(tdist)
}

#' \code{expression.dotplot} helps in visualising the sample classification
#'
#' @param data A data.frame containing gene (rows) expression from the samples (columns)
#' to be classified
#' @param predictions A data.frame containing the predictions obtained with
#' \code{classify_splits} function;
#' @param g1 A gene symbol to be visualized on x-axis;
#' @param g2 A gene symbol to be visualized on y-axis;
#' @param confidence If TRUE, the colors correspond to classification
#' confidence rather than predicted classes;
#'
#' @import ggplot2
#' @return \code{classify_splits} Returns a ggplot object


expression.dotplot <- function(data, predictions,g1,g2,confidence=F){
  prediction <- NULL
  colnames(data) <- gsub("-",".",colnames(data))
  classes <- c("Basal","Luminal","MA","unclassified")
  predictions$pred <- factor(as.vector(predictions$prediction.strict),levels=c("Luminal","Basal","MA","unclassified"))
  if (confidence){
    classes <- c("high","medium","low")
  }
  data <- as.data.frame(t(data))
  data$patientID <- row.names(data)
  #data <- data[!grepl("^TCGA",data$patientID),]
  data <- data[,c("patientID",g1,g2)]
  colnames(data) <- c("patientID","g1","g2")
  data$prediction <- predictions$pred[match(data$patientID,row.names(predictions))]
  data <- data[order(data$prediction),]
  color.LAB <- c(Basal="red",unclassified="grey",Luminal="darkblue",MA="pink",high="red",medium="orange",low="grey")
  color.LAB<- color.LAB[names(color.LAB)%in%classes]
  g <- ggplot(data,aes(x=g1,y=g2)) + theme_bw() +
    labs(x=paste(g1,"expression"),y=paste(g2,"expression")) +
    scale_colour_manual(values=color.LAB) +
    theme(legend.position="top",
          legend.title=element_blank(),
          aspect.ratio = 1
    )
  if (confidence){
    data$confidence <- predictions$confidence[match(data$patientID,row.names(predictions))]
    g <- g + geom_point(data=data,aes(color=confidence))
  }
  else {
    g <- g + geom_point(aes(color=prediction))
    g <- g + geom_point(data=data[data$prediction=="MA",],shape=1,color="mediumorchid")
  }

  return(g)
}

#' \code{LABclassifier} helps in visualising the sample classification
#'
#' @param data A data.frame containing gene (gene symbols as rows) expression
#' from the samples (sample IDs as columns) to be classified or a  vector
#' containing gene expression named with gene SYMBOL
#' @param LABClassif An object of class \code{LABclassifier} obtained with
#' \code{make.LABclassifier} function. If NULL, the precomputed classifier
#' is used (default);
#' @param raw.counts A logical. Set it to TRUE if you are using raw sequencing
#' data, then they will be transformed in TPM. Default: FALSE;
#' @param plot A logical. Set it to TRUE to display classification
#' diagnostic plots. Default: FALSE
#'
#' @import ggplot2
#' @return \code{LABclassifier} Returns a data.frame with columns
#' \item{sensory2secretory}{Distance to sensory centroid}
#' \item{Luminal2Apocrine}{Distance to Luminal centroid}
#' \item{prediction}{Centroid association: Luminal, Basal, MA: Molecular Apocrine}
#' \item{prediction.strict}{Centroid association with annotation of the samples
#' too far from any centroid as unclassified}
#' \item{confidence}{Summarizes the distance to centroids, it is just an indicator:
#' low confidence : 0.9 > scores < 1.1 (unclassified samples in strict predictions)
#' medium confidence : 0.75 > scores < 1.25
#' high confidence : 0.75 < scores > 1.15}
#' @export
#' @examples
#' data(TCGA.rsem)
#' LABclassifier(TCGA.rsem,plot=TRUE)


LABclassifier <- function(data,LABClassif=NULL,raw.counts=F,plot=F){
  if (raw.counts){
    message("Transforming RNA-seq raw counts into log2(TPM+1)")
    m <- match(gene.length$SYMBOL,row.names(data))
    ldata <- data[na.omit(m),]
    gene.length <- gene.length[!is.na(m),]
    ldata <- (ldata / (gene.length$effective_length+1))*1000
    TPM <- apply(ldata,2,function(x){
      x <- (x/sum(x))*10^6
    })
    data <- log2(TPM+1)
  }
  if (!is.null(LABClassif)){
    if (class(LABClassif)!="LABclassifier"){
      stop("Error - LAB classifier must be created using make.LABclassifier function.")
    }
    centroid.ss <- LABClassif$centroid.ss
    centroid.la <- LABClassif$centroid.la
  }
  sensory2secretory <- Luminal2Apocrine <- prediction.strict <- xend <- yend <- confidence <- NULL
  message("Applying sensory/secretory splitting ...")
  pred.ss <- classify_splits(data,centroid.ss,"ss")
  message("Applying Luminal/Apocrine splitting ...")
  pred.la <- classify_splits(data,centroid.la,"la")

  message("Computing final predictions ...")
  final.pred <- ifelse(pred.ss$pred.ss=="secretory","Basal",pred.la$pred.la)
  final.pred.strict <- final.pred
  final.pred.strict[pred.ss$pred.ss.strict=="unclassified" | (pred.la$pred.la.strict=="unclassified" & pred.ss$pred.ss.strict=="sensory")] <- "unclassified"

  pred.final <- data.frame(sensory2secretory=pred.ss$sensory,Luminal2Apocrine=pred.la$Luminal,prediction=final.pred,prediction.strict=final.pred.strict,row.names=colnames(data))
  pred.final$prediction.strict <- factor(as.vector(pred.final$prediction.strict),levels=c("Luminal","Basal","MA","unclassified"))

  message("Computing confidence ...")
  s2s <- cut(pred.final$sensory2secretory,breaks = c(0,0.75,0.9,1.1,1.25,2))
  levels(s2s) <- c("high","medium","low","medium","high")
  l2a <- cut(pred.final$Luminal2Apocrine,breaks = c(0,0.75,0.9,1.1,1.25,2))
  levels(l2a) <- c("high","medium","low","medium","high")
  lab <- s2s
  lab[pred.ss$pred.ss=="sensory"] <- l2a[pred.ss$pred.ss=="sensory"]
  lab[lab!="low" & s2s=="low"] <- "low"
  lab[lab=="high" & s2s=="medium"] <- "medium"
  pred.final$confidence <- lab
  if (plot){
    message("Start plotting ...")
    thrs.lines <- data.frame(sensory2secretory=rep(0,3),xend=rep(1,3),Luminal2Apocrine=c(0.9,1,1.1),yend=c(0.9,1,1.1))
    g <- ggplot(pred.final,aes(x=sensory2secretory,y=Luminal2Apocrine)) + theme_bw() +
      lims(x=c(0,2),y=c(0,2)) +
      geom_vline(xintercept = c(0.9,1,1.1),linetype=c("dashed","solid","dashed"),size=rep(0.5,3),color=rep("darkgrey",3)) +
      geom_segment(data=thrs.lines,aes(xend=xend,yend=yend),linetype=c("dashed","solid","dashed"),size=0.5,color="darkgrey") +
      labs(x="Sensory to Secretory score",y="Luminal to Apocrine score") + theme(aspect.ratio = 1,legend.position="top")
      g1 <- g + geom_point(aes(color=prediction.strict)) + scale_color_manual(values=c(Basal="red",Luminal="darkblue",MA="pink",unclassified="grey")) + geom_point(data=pred.final[pred.final$prediction.strict=="MA",],shape=1,color="purple")
      g9 <- g + geom_point(data=pred.final,aes(color=confidence)) + scale_color_manual(values=c(high="red",medium="orange",low="grey"))
      if (is.data.frame(data)){
        g3 <- expression.dotplot(data, pred.final,"ESR1","FOXA1")
        g10 <- expression.dotplot(data, pred.final,"ESR1","FOXA1",confidence=T)
        g5 <- expression.dotplot(data, pred.final,"ESR1","AR")
        g11 <- expression.dotplot(data, pred.final,"ESR1","AR",confidence=T)
        g7 <- expression.dotplot(data, pred.final,"ESR1","ERBB2")
        g12 <- expression.dotplot(data, pred.final,"ESR1","ERBB2",confidence=T)
        multiplot(g9,g1,g10,g3,g11,g5,g12,g7,cols=4)
      } else {
        multiplot(g9,g1,cols = 2)
      }
  }
  return(pred.final)
}




