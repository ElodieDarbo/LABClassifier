import.data.zenodo <- function(dataset="all"){
  dir.create("ext_data",recursive = T,showWarnings = FALSE)
  data <- c(EORTC="./ext_data/EORTC_x3p_matrix.txt",ICGC="./ext_data/ICGC_matrix.txt",
               METABRIC= "./ext_data/METABRIC_matrix.txt", TCGA="./ext_data/TCGA_matrix.txt")
  if (dataset[1]!="all" & sum(dataset%in%names(data))!=length(dataset)){
    dataset <- dataset[!dataset%in%names(data)]
    stop(paste("The dataset(s)",paste(dataset,collapse=" "),"are not part of the available data.\n dataset must contain: all, TCGA, METABRIC, EORTC, ICGC."))
  }
  l <- list()
  if (length(dataset)==1 & dataset == "all"){
    dataset <- names(data)
    message(paste(dataset, collapse=", ")," will be dowloaded ...")
    download_zenodo(doi="10.5281/zenodo.10935179",path="./ext_data",files=basename(data),timeout=1000)
  }
  else {
    dataset <- dataset[dataset%in%names(data)]
    data <- data[dataset]
    message(paste(dataset, collapse=", ")," will be dowloaded ...")
    download_zenodo(doi="10.5281/zenodo.10935179",path="./ext_data",files=basename(data),timeout=1000)
  }
  for (d in names(data)){
    message("Loading ",d, " ...")
    tmp <- read.table(data[d],sep="\t",head=T,row=1)
    print(head(tmp[,1:3]))
    l[[d]] <- tmp
  }
  return(l)
}


# Export current plot in different formats and sizes; from Jacques van Helden
export.plot <- function (file.prefix="PlotExport",
                         export.formats="pdf", # supported: postscript, jpg, png, bmp, pdf
                         width=11, # in inches
                         height=8, # in inches
                         horizontal=T,
                         ... ## Additional parameters are passed to the export method
) {

  ppi <- 72
  file.ext <- c(
    postscript = "ps",
    pdf = "pdf",
    ps = "ps",
    eps = "eps",
    jpeg="jpg",
    jpg="jpg",
    bmp="bmp",
    png="png",
    svg="svg",
    tiff="tiff")
  for (f in export.formats) {
    from.dev <- dev.cur();

    file.name <- paste(file.prefix,file.ext[f], sep=".")

    if ((f == "postscript") || (f == "ps")) {
      postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal, ...)
    } else if (f == "eps") {
      postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal,onefile=F, ...)
    } else if (f == "pdf") {
      pdf(file.name, paper="special",width=width,height=height, ...)
    } else if ((f == "jpg") || (f == "jpeg")) {
      jpeg(file.name,width=(width*ppi),height=(height*ppi),quality=100, ...)
    } else if (f == "png") {
      png(file.name,width=width*ppi,height=height*ppi, ...)
    } else if (f == "bmp") {
      bitmap(file.name,width=width*ppi,height=height*ppi, ...)
    } else if (f == "svg") {
      svg(file.name,width=width*ppi,height=height*ppi, ...)
    } else if (f == "tiff") {
      #tiff(filename = "Rplot%03d.tiff", width = 480, height = 480, units = "px", pointsize = 12, compression = c("none", "rle", "lzw", "jpeg", "zip"), bg = "white", res = NA,  ..., type = c("cairo", "Xlib", "quartz"), antialias)
      tiff(file.name,width=width*ppi,height=height*ppi, compression = 'none', ...)
    }
    else {
      print(paste("Error: format ", f, " is not supported", sep=""))
      return()
    }
    to.dev <- dev.cur()
    dev.set(which=from.dev)
    dev.copy(which=to.dev)
    dev.set(which=to.dev)
    dev.off()
    dev.set(which=from.dev) ## This is required because dev.off() returns to the first, not the last, device
  }
}



#' \code{multiplot} helps to aggregate multiple ggplot objects in one plot
#'
#' @param ... ggplot objects comma separated;
#' @param plotlist A list of ggplot objects;
#' @param cols Number of columns in layout
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored.
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' @author http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
#' @export
#' @import grDevices
#' @import graphics
#' @import stats
#' @import grid


multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


