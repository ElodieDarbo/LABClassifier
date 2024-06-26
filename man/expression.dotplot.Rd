% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LABClassifier.R
\name{expression.dotplot}
\alias{expression.dotplot}
\title{\code{expression.dotplot} helps in visualising the sample classification}
\usage{
expression.dotplot(data, predictions, g1, g2, PAM50 = F)
}
\arguments{
\item{data}{A data.frame containing gene (rows) expression from the samples (columns)
to be classified}

\item{predictions}{A data.frame containing the predictions obtained with
\code{classify_splits} function;}

\item{g1}{A gene symbol to be visualized on x-axis;}

\item{g2}{A gene symbol to be visualized on y-axis;}

\item{PAM50}{If TRUE, the colors correspond to PAM50 classification}
}
\value{
\code{classify_splits} Returns a ggplot object
}
\description{
\code{expression.dotplot} helps in visualising the sample classification
}
