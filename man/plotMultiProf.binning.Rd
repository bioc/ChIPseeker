% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotTagMatrix.R
\name{plotMultiProf.binning}
\alias{plotMultiProf.binning}
\title{internal function}
\usage{
plotMultiProf.binning(
  tagMatrix,
  xlab = "Genomic Region (5'->3')",
  ylab = "Peak Count Frequency",
  conf,
  facet = "none",
  free_y = TRUE,
  upstream = NULL,
  downstream = NULL,
  label,
  ...
)
}
\arguments{
\item{tagMatrix}{tagMatrix}

\item{xlab}{xlab}

\item{ylab}{ylab}

\item{conf}{confidence interval}

\item{facet}{one of 'none', 'row' and 'column'}

\item{free_y}{if TRUE, y will be scaled by AvgProf}

\item{upstream}{the upstream extension}

\item{downstream}{the downstream extension}

\item{label}{the label of the center}

\item{...}{additional parameter}
}
\description{
internal function
}
