% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_PCA.R
\name{PCAPlot}
\alias{PCAPlot}
\title{Principal Component Analysis (PCA)}
\usage{
PCAPlot(
  data,
  group,
  show.PCA.labels = FALSE,
  cols = NULL,
  prefix = "",
  scale = FALSE
)
}
\arguments{
\item{data}{a dataframe of expression/abundance counts}

\item{group}{a factor specifying group for each sample (e.g. could be
represented by a column from a metadata file)}

\item{show.PCA.labels}{a boolean value (TRUE or FALSE) specifying if elements
(e.g. samples) should be labelled in the PCA plot. Labeling is based on
column names of the input data. Default value is FALSE.}

\item{cols}{a vector of colors (one color for each group)}

\item{prefix}{a character defining a prefix of output file}

\item{scale}{a boolean, if TRUE then data are scaled to unit variance. Default is FALSE}
}
\value{
1) scree plot; 2) plot of contributions of variables to PC1; 3) plot
of contributions of variables to PC2; 4) 2D PCA plot projecting samples over
first 2 principal components
}
\description{
A function calculating PCA. PCA plot shows the relationships
between samples (based on squared euclidean distances) in the data set.
Each dot represents one sample. This function provides also a scree plot
describing percentages of explained variances by each PCA component and
plots with contributions of variables to PCA1 and PCA2.
}
\examples{
{
PCAPlot(campp2_brca_1_batchCorrected, as.factor(campp2_brca_1_meta$subtype), show.PCA.labels =FALSE, cols=NULL, prefix="test_PCA")
}
}
