% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot_distributions.R
\name{PlotDistributions}
\alias{PlotDistributions}
\title{Plot Distributions}
\usage{
PlotDistributions(data, fitted.data)
}
\arguments{
\item{data}{a data.frame of feature counts, the same which was
used as an input for "FitDistributions" function.}

\item{fitted.data}{a list of distribution descriptions obtained as an output
from the function "FitDistributions" where the same feature counts are used
as an input.}
}
\value{
a .pdf file including multiple plots:
1) Cullen and Frey graph
2) Histogram and theoretical densities
3) Empirical and theoretical CDFs
4) Q-Q plot
5) P-P plot
}
\description{
Plotting the feature counts' Distributions.
As results, multiple plots are obtained for each feature:
1) Cullen and Frey graph
2) Histogram and theoretical densities
3) Empirical and theoretical CDFs
4) Q-Q plot
5) P-P plot

For some of the figures, plotted distributions depends on the characteristics
of the input data:
In case of only integer numbers are present in the dataset, poison and normal
distributions are tested. In case of presence of float values,
Weibull, gamma, log normal and normal distributions are tested. In case of
negative float values, only normal distribution is calculated.

Resulting figures are saved into the .pdf file.
}
\examples{
{
PlotDistributions(campp2_brca_1_batchCorrected[1:10,], campp2_brca_1_distributionsFit)
}
}
