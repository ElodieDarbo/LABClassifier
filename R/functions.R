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
    #tdist <- as.data.frame(matrix(1-c(tdist[1,2],tdist[1,3]),ncol=2))
    tdist <- as.data.frame(matrix(1-as.vector(tdist[1,2:(length(classes)+1)]),ncol=length(classes)))
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
    tdist$pred.ss.strict[tdist$sensory>=0.9 & tdist$sensory<=1.1] <- "unclassified"
  } else if (opt == "la"){
    pred <- ifelse(tdist$Luminal <= 1,"Luminal","MA")
    tdist$pred.la <- pred
    tdist$pred.la.strict <- pred
    tdist$pred.la.strict[tdist$Luminal>=0.9 & tdist$Luminal<=1.1] <- "unclassified"
  }
  return(tdist)
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
#' @param confidence If TRUE, the colors correspond to classification
#' confidence rather than predicted classes;
#' @param PAM50 If TRUE, the colors correspond to PAM50 classification
#'
#' @import ggplot2
#' @return \code{classify_splits} Returns a ggplot object


expression.dotplot <- function(data, predictions,g1,g2,confidence=F,PAM50=F){
  prediction <- NULL
  colnames(data) <- gsub("-",".",colnames(data))
  classes <- c("Basal","Luminal","MA","unclassified")
  predictions$pred <- factor(as.vector(predictions$prediction.strict),levels=c("Luminal","Basal","MA","unclassified"))
  if (confidence){
    classes <- c("high","medium","low")
  }
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
  #scale_color_manual(values=c()) + geom_point(data=pred.final[pred.final$PAM50=="Her2",],shape=1,color="deeppink")
  color.LAB <- c(LumA="darkblue",Her2="pink",LumB="lightblue",Basal="red",unclassified="grey",Luminal="darkblue",MA="pink",high="red",medium="orange",low="grey")
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
    g <- g + geom_point(data=data,aes(color=confidence)) +
      scale_color_gradientn(colours = c("yellow","orange","red","black"),limits=c(0,0.5),values = c(0,0.02,0.1,0.2,0.5,1))

  }
  else {
    g <- g + geom_point(aes(color=prediction))
    if (PAM50){
      g <- g + geom_point(data=data[data$prediction=="Her2",],shape=1,color="deeppink")
    } else {
      g <- g + geom_point(data=data[data$prediction=="MA",],shape=1,color="mediumorchid")
    }
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
#' @param log2T A logical. If true, data are log2 transformed. Default: FALSE;
#' @param id.type A character string. Gene identifier type in the input matrix. Possible values are
#' ensEMBL, entrezID, SYMBOL (default)
#' @param PAM50 A logical. Compute PAM50 classification using \code{molecular.subtyping}
#' function from genefu package. If used, please cite the corresponding paper (See reference).
#' @param plot A logical. Set it to TRUE to display classification
#' diagnostic plots. Default: FALSE
#' @param prefix A character string. If set, a directory Output is created in the
#' working directory (if not exists), and the prediction result is saved in as a
#' tabulated file. The filename is prefix_LAB_prefix.tab
#' @references
#'
#' genefu:  Deena M.A. Gendoo et al. (2021). genefu: Computation of Gene
#' Expression-Based Signatures in Breast Cancer. R package version 2.26.0.
#' http://www.pmgenomics.ca/bhklab/software/genefu
#' @import ggplot2
#' @importFrom utils write.table
#' @return \code{LABclassifier} Returns a data.frame with columns
#' \item{sensory2secretory}{Distance to sensory centroid}
#' \item{Luminal2Apocrine}{Distance to Luminal centroid}
#' \item{prediction}{Centroid association: Luminal, Basal, MA: Molecular Apocrine}
#' \item{prediction.strict}{Centroid association with annotation of the samples
#' too far from any centroid as unclassified}
#' \item{confidence}{Summarizes the distance to centroids, it is just an indicator:
#' \itemize{
#' \item{high confidence}{ : 0.75 < scores > 1.15}
#' \item{medium confidence}{ : 0.75 > scores < 1.25}
#' \item{low confidence{ : 0.9 > scores < 1.1 (unclassified samples in strict predictions)}
#' }}}
#' @export
#' @examples
#' data(TCGA.rsem)
#' LABclassifier(TCGA.rsem,plot=TRUE)


LABclassifier <- function(data,LABClassif=NULL,raw.counts=F,log2T=F,id.type="SYMBOL",PAM50=F,plot=F,prefix=NULL){

  v <- F
  if (is.vector(data)){
    if (is.null(names(data))) stop("ERROR: You submitted a vector that is not named.")
    data <- data.frame(row.names=names(data),sample=data)
    v <- T
    if (PAM50){
      message("WARNING: PAM50 can't be applied on a simgle sample.")
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
  if (!is.null(LABClassif)){
    if (class(LABClassif)!="LABclassifier"){
      stop("Error - LAB classifier must be created using make.LABclassifier function.")
    }
    centroid.ss <- LABClassif$centroid.ss
    centroid.la <- LABClassif$centroid.la
  }
  sensory2secretory <- Luminal2Apocrine <- prediction.strict <- xend <- yend <- confidence <- NULL

  slipts.genes <- list(ss=row.names(centroid.ss$centroids),la=row.names(centroid.la$centroids))
  slipts.genes.pc <- list(ss=sum(slipts.genes$ss%in%row.names(data))/length(slipts.genes$ss),la=sum(slipts.genes$la%in%row.names(data))/length(slipts.genes$la))
  all.symbols <- row.names(data)
  if (v){
    g <- row.names(data)
    data <- data$sample
    names(data) <- g
  }

  message("Applying sensory/secretory splitting ...")
  if (slipts.genes.pc$ss!=1){
    message("WARNING: ",paste(slipts.genes$ss[!slipts.genes$ss%in%all.symbols],collapse=", "), " (",round(1 - slipts.genes.pc$ss,2)*100,"%) ",
            "are absent from your data.")
  }
  pred.ss <- classify_splits(data,centroid.ss,"ss")
  message("Applying Luminal/Apocrine splitting ...")
  if (slipts.genes.pc$la!=1){
    message("WARNING: ",paste(slipts.genes$la[!slipts.genes$la%in%all.symbols],collapse=", "), " (",round(1 - slipts.genes.pc$la,2)*100,"%) ",
            "are absent from your data.")
  }

  pred.la <- classify_splits(data,centroid.la,"la")


  # s2s <- cut(pred.final$sensory2secretory,breaks = c(0,0.75,0.9,1.1,1.25,2))
  # levels(s2s) <- c("high","medium","low","medium","high")
  # l2a <- cut(pred.final$Luminal2Apocrine,breaks = c(0,0.75,0.9,1.1,1.25,2))
  # levels(l2a) <- c("high","medium","low","medium","high")
  # lab <- s2s
  # lab[pred.ss$pred.ss=="sensory"] <- l2a[pred.ss$pred.ss=="sensory"]
  # lab[lab!="low" & s2s=="low"] <- "low"
  # lab[lab=="high" & s2s=="medium"] <- "medium"
  # pred.final$confidence <- lab

  #pred.final$combined.score <- pred.final$Luminal2Apocrine + pred.final$sensory2secretory

  compute.pvalue <- function(x,y,class){
    corner <- c(0,0)
    switch(class,
           luminal={
             wa <- 0.8
             wb <- 0.7
           },
           basal={
             x <- 2 - x
             y <- 2 - y
             wa <- 1.2
             wb <- 0
           },
           apocrine={
             y <- 2 - y
             wa <- 0.6
             wb <- 0.7
           }
           )

    distances <- sqrt(wa*(x-corner[1])^2+wb*(y-corner[2])^2)

    compute_slope <- function(x1, y1, x2, y2) {
      slope <- (y2 - y1) / (x2 - x1)
      return(slope)
    }

    # Function to compute the angle at the crossing of two lines
    compute_angle <- function(slope1, slope2) {
      # Calculate the angle in radians
      angle_rad <- atan(abs((slope2 - slope1) / (1 + slope1 * slope2)))

      # Convert the angle to degrees
      angle_deg <- angle_rad * 180 / pi

      return(angle_deg)
    }

    distances <- apply(data.frame(x,y,distances),1,function(p){
      x <- p[1]
      y <- p[2]
      distances <- p[3]
      corner <- c(0,0)
      slope1 <- compute_slope(0,0,1,0)
      slope2 <- compute_slope(0,0,x,y)
      angle <- compute_angle(slope1,slope2)
      b <- 1
      if (angle<45){
        max.dist <- sqrt((b-corner[1])^2+(b*slope2-corner[2])^2)
      } else {
        max.dist <- sqrt((b/slope2-corner[1])^2+(b-corner[2])^2)
      }
      distances <- distances / max.dist

      return(distances)
    })
    #print(plot(density(distances)))
    round(pnorm(distances,mean = 0,sd = 0.31,lower.tail = F),7)

    #pgamma(distances,shape = 1,scale = 0.25,lower.tail = F)
  }

  pred.final <- data.frame(sensory2secretory=pred.ss$sensory,Luminal2Apocrine=pred.la$Luminal,row.names=colnames(data))

  message("Computing confidence ...")
  #sprintf("%.10f", data$Value)
  pred.final$test.bas <- compute.pvalue(pred.final$sensory2secretory,pred.final$Luminal2Apocrine,"basal")
  pred.final$test.lum <- compute.pvalue(pred.final$sensory2secretory,pred.final$Luminal2Apocrine,"luminal")
  pred.final$test.apo <- compute.pvalue(pred.final$sensory2secretory,pred.final$Luminal2Apocrine,"apocrine")

  pred.final$confidence <- apply(pred.final[,c("test.lum","test.apo","test.bas")],1,max)

  message("Computing final predictions ...")
  #final.pred <- ifelse(pred.ss$pred.ss=="secretory","Basal",pred.la$pred.la)
  #final.pred.strict <- final.pred
  #final.pred.strict[pred.ss$pred.ss.strict=="unclassified" | (pred.la$pred.la.strict=="unclassified" & pred.ss$pred.ss.strict=="sensory")] <- "unclassified"

  pred.final$prediction <- apply(pred.final[,c("test.bas","test.lum","test.apo")],1,function(x){
    c("Basal","Luminal","MA")[which.max(x)]
  })

  pred.final$prediction.strict <- ifelse(pred.final$test.apo>=0.01,"MA","unclassified")
  pred.final$prediction.strict[pred.final$test.lum>=0.01] <- "Luminal"
  pred.final$prediction.strict[pred.final$test.bas>=0.01] <- "Basal"

  #pred.final <- data.frame(sensory2secretory=pred.ss$sensory,Luminal2Apocrine=pred.la$Luminal,prediction=final.pred,prediction.strict=final.pred.strict,row.names=colnames(data))
  pred.final$prediction.strict <- factor(as.vector(pred.final$prediction.strict),levels=c("Luminal","Basal","MA","unclassified"))

  if (PAM50){
    message("Computing PAM50 classification")
    pam50.pred <- PAMgenefu(data)
    pred.final <- merge(pred.final,pam50.pred,by=0)
    row.names(pred.final) <- pred.final$Row.names
  }

  if (plot){
    message("Start plotting ...")
    thrs.lines <- data.frame(sensory2secretory=rep(0,3),xend=rep(1,3),Luminal2Apocrine=c(0.9,1,1.1),yend=c(0.9,1,1.1))
    g <- ggplot(pred.final,aes(x=sensory2secretory,y=Luminal2Apocrine)) + theme_bw() +
      #lims(x=c(0,2),y=c(0,2)) +
      geom_vline(xintercept = c(0.9,1,1.1),linetype=c("dashed","solid","dashed"),size=rep(0.5,3),color=rep("darkgrey",3)) +
      geom_segment(data=thrs.lines,aes(xend=xend,yend=yend),linetype=c("dashed","solid","dashed"),size=0.5,color="darkgrey") +
      labs(x="Distance to sensory centroid",y="Distance to Luminal centroid") + theme(aspect.ratio = 1,legend.position="top")
    g1 <- g + geom_point(aes(color=prediction.strict)) + scale_color_manual(values=c(Basal="red",Luminal="darkblue",MA="pink",unclassified="grey")) +
        geom_point(data=pred.final[pred.final$prediction.strict=="MA",],shape=1,color="purple") +
        coord_flip() + scale_y_reverse() + scale_x_reverse()
      #g9 <- g + geom_point(data=pred.final,aes(color=confidence)) + scale_color_manual(values=c(high="red",medium="orange",low="grey"))
    g9 <- g + geom_point(aes(color=confidence)) +
        scale_color_gradientn(colours = c("yellow","orange","red","black"),limits=c(0,0.5),values = c(0,0.02,0.1,0.2,0.5,1)) +
        coord_flip() + scale_y_reverse() + scale_x_reverse()
    if (PAM50){
        gpam <- g + geom_point(data=pred.final,aes(color=PAM50)) + scale_color_manual(values=c(Basal="red",LumA="darkblue",Her2="pink",LumB="lightblue")) + geom_point(data=pred.final[pred.final$PAM50=="Her2",],shape=1,color="deeppink")+
          coord_flip() + scale_y_reverse() + scale_x_reverse()
    }
    if (!is.vector(data)){
        g3 <- expression.dotplot(data, pred.final,"ESR1","FOXA1")
        g10 <- expression.dotplot(data, pred.final,"ESR1","FOXA1",confidence=T)
        g5 <- expression.dotplot(data, pred.final,"ESR1","AR")
        g11 <- expression.dotplot(data, pred.final,"ESR1","AR",confidence=T)
        g7 <- expression.dotplot(data, pred.final,"ESR1","ERBB2")
        g12 <- expression.dotplot(data, pred.final,"ESR1","ERBB2",confidence=T)
        if (PAM50){
          gpam1 <- expression.dotplot(data, pred.final,"ESR1","FOXA1",PAM50=T)
          gpam2 <- expression.dotplot(data, pred.final,"ESR1","AR",PAM50=T)
          gpam3 <- expression.dotplot(data, pred.final,"ESR1","ERBB2",PAM50=T)
          multiplot(g9,g1,gpam,g10,g3,gpam1,g11,g5,gpam2,g12,g7,gpam3,cols=4)
        } else {
          multiplot(g9,g1,g10,g3,g11,g5,g12,g7,cols=4)
        }
      } else {
        multiplot(g9,g1,cols = 2)
      }
  }
  pred.final$Row.names <- NULL
  if (!is.null(prefix)){
    message("If not existing, creating a folder 'Output' in the working directory: ",getwd())
    suppressWarnings(dir.create("./Output",recursive = T))
    pred.final$test.bas <- sprintf("%.7f", pred.final$test.bas)
    pred.final$test.lum <- sprintf("%.7f", pred.final$test.lum)
    pred.final$test.apo <- sprintf("%.7f", pred.final$test.apo)

    pred.final$confidence <- sprintf("%.7f", pred.final$confidence)

    write.table(pred.final,paste0("./Output/",prefix,"_LAB_predictions.tab"),sep="\t",quote=F)
    if (plot){
      w <- 16
      h <- 8
      if (PAM50){
        h <- 12
      }
      if (is.vector(data)){
        h <- 4
        w <- 8
      }
      export.plot(paste0("./Output/",prefix,"_LAB_predictions"),width=w,height=h)
    }
  }
  return(invisible(pred.final))
}




