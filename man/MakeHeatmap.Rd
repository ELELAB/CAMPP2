% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Heatmap.R
\name{MakeHeatmap}
\alias{MakeHeatmap}
\title{Make Heatmap}
\usage{
MakeHeatmap(
  data,
  groups,
  prefix,
  viridis.palette = "turbo",
  data.type = "expression/abundance"
)
}
\arguments{
\item{data}{a feature count matrix (ideally normalized and batch corrected)
from "seq", "array", "ms" or "other" technology (with feature IDs as row
names and sample IDs as columns). It's recommended to import feature counts
using function "import_counts".}

\item{groups}{a factor specifying group for each sample (e.g. could be
represented by a column from a metadata file). Group's information
will be used in clustering the samples.}

\item{prefix}{a character vector defining prefix of output file name.}

\item{viridis.palette}{a character vector specifying color gradient to use for
the heatmap. Default option for viridis color palettes is 'turbo'. For more
options, see viridis vignette.}

\item{data.type}{a character vector (used in the heatmap legend) describing
type of the data being visualized. Default is "expression/abundance".}
}
\value{
heatmap
}
\description{
A function for making heatmaps to showcase the difference in
expression/abundance of features across samples and sample groups and to
visualise relations between them.
}
\examples{
{
MakeHeatmap(data=campp2_brca_1_batchCorrected,
groups=campp2_brca_1_meta$subtype, prefix="test", viridis.palette = "turbo",
data.type = "Feature counts")
}
}
