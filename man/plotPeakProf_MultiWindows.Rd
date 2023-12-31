% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotTagMatrix.R
\name{plotPeakProf_MultiWindows}
\alias{plotPeakProf_MultiWindows}
\title{plotPeakProf_MultiWindows}
\usage{
plotPeakProf_MultiWindows(
  peak,
  upstream,
  downstream,
  conf,
  by,
  type,
  windows_name = NULL,
  weightCol = NULL,
  TxDb = NULL,
  xlab = "Genomic Region (5'->3')",
  ylab = "Peak Count Frequency",
  facet = "row",
  free_y = TRUE,
  verbose = TRUE,
  nbin = NULL,
  ignore_strand = FALSE,
  ...
)
}
\arguments{
\item{peak}{peak file or GRanges object}

\item{upstream}{upstream position}

\item{downstream}{downstream position}

\item{conf}{confidence interval}

\item{by}{feature of interest}

\item{type}{one of "start_site", "end_site", "body"}

\item{windows_name}{the name for each window, which will also be showed in the picture as labels}

\item{weightCol}{column name of weight}

\item{TxDb}{TxDb object or self-made granges objects}

\item{xlab}{xlab}

\item{ylab}{ylab}

\item{facet}{one of 'none', 'row' and 'column'}

\item{free_y}{if TRUE, y will be scaled by AvgProf}

\item{verbose}{print message or not}

\item{nbin}{the amount of bines}

\item{ignore_strand}{ignore the strand information or not}

\item{...}{additional parameter}
}
\value{
ggplot object
}
\description{
plot the profile of peaks in two or more windows
}
\details{
This function comes from \url{https://github.com/YuLab-SMU/ChIPseeker/issues/189}
`
\code{plotPeakProf_MultiWindows()} is almost the same as \code{plotPeakProf2()}, having
the main difference of accepting two or more granges objects. Accepting more
granges objects can help compare the same peaks in different windows.

\code{TxDb} parameter can accept txdb object.
But many regions can not be obtained by txdb object. In this case,
Users can provide self-made granges served the same role 
as txdb object and pass to \code{TxDb} object.

\code{by} the features of interest. 

(1) if users use \code{txdb}, \code{by} can be one of 'gene', 'transcript', 'exon', 
'intron' , '3UTR' , '5UTR', 'UTR'. These features can be obtained by functions from txdb object.

(2) if users use self-made granges object, \code{by} can be everything. Because this \code{by}
will not pass to functions to get features, which is different from the case of using 
txdb object. This \code{by} is only used to made labels showed in picture.

\code{type} means the property of the region. one of the "start site",
"end site" and "body".

\code{upstream} and \code{downstream} parameter have different usages:

(1) if \code{type == 'body'}, \code{upstream} and \code{downstream} can use to extend 
the flank of body region.

(2) if \code{type == 'start_site'/'end_site'}, \code{upstream} and \code{downstream} refer to
the upstream and downstream of the start_site or the end_site.

\code{weightCol} refers to column in peak file. This column acts as a weight value. Details
see \url{https://github.com/YuLab-SMU/ChIPseeker/issues/15}

\code{nbin} refers to the number of bins. \code{getTagMatrix()} provide a binning method
to get the tag matrix.

There are two ways input a list of window.

(1) Users can input a list of self-made granges objects

(2) Users can input a list of \code{by} and only one \code{type}. In this way, 
\code{plotPeakProf_MultiWindows()} can made a list of window from txdb object based on \code{by} and \code{type}.

Warning: 

(1) All of these window should be the same type. It means users can only
compare a list of "start site"/"end site"/"body region" with the same upstream
and downstream.

(2) So it will be only one \code{type} and several \code{by}.

(3) Users can make window by txdb object or self-made granges object. Users can only
choose one of 'gene', 'transcript', 'exon', 'intron' , '3UTR' , '5UTR' or 'UTR' in the
way of using txdb object. User can input any \code{by} in the way of using 
self-made granges object.

(4) Users can mingle the \code{by} designed for the two ways. \code{plotPeakProf_MultiWindows} can
accpet the hybrid \code{by}. But the above rules should be followed.
}
