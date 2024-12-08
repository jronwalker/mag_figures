% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculations.R
\name{comm_long}
\alias{comm_long}
\title{Generate table ready for ggplot2.}
\usage{
comm_long(abund_df, rank, min_cutoff = 0.01)
}
\arguments{
\item{abund_df}{Data frame object containing taxonomic information and proportional/relative abundance information.
Works best with the frac_abund output but any data frame object should work. Taxonomic column names should
be one or multiple of the following: "bin", "domain", "phylum", "class", "order", "family", "genus", "species".
Sample column names can be anything but the columns should be numeric.}

\item{rank}{Taxonomic rank at which you would like to make the graph. If none chosen the lowest will be used.}

\item{min_cutoff}{The cutoff to aggregate a group into the "Others". The cutoff is compared with the maximum
relative abundance of that group in your each of your samples. If that group does not surpass the cutoff in
any sample it is added to "Others". Default is 0.01 (i.e. 1\% in proportional data).}
}
\value{
The output will be a data frame containing the fractional abundance of each taxonomic rank chosen.
}
\description{
The functions takes a community data frame where columns are taxonomic information and abundances
and rows are unique taxa (MAGs, ASVs, OTUs) and returns a long formatted table for ggplot2.
}
