% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matcher.R
\name{MXmatcher}
\alias{MXmatcher}
\title{specific matcher for MXE}
\usage{
MXmatcher(
  i,
  below_thresh = 0,
  redExon,
  minOverlap = 0.05,
  gtf_transcripts,
  gtf_exons,
  protein_coding_transcripts,
  transcript_starts
)
}
\arguments{
\item{i}{index of redExon}

\item{below_thresh}{threshold to make skipped form match less than}

\item{redExon}{dataframe generated initially in getForeground}

\item{minOverlap}{minimum overlap to count as same exon}

\item{gtf_transcripts}{transcripts gtf dataframe}

\item{gtf_exons}{exons gtf dataframe}

\item{protein_coding_transcripts}{protein coding transcripts}

\item{transcript_starts}{start location of transcript}
}
\value{
matched transcript rownumber
}
\description{
specific matcher for MXE
}
\keyword{internal}
