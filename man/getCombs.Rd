% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dStruct.R
\name{getCombs}
\alias{getCombs}
\title{Identifies subgroupings of replicates for assessing within-group and between-group variation.}
\usage{
getCombs(reps_A, reps_B, batches = F, between_combs = NULL,
  within_combs = NULL)
}
\arguments{
\item{reps_A}{Number of replicates of group A.}

\item{reps_B}{Number of replicates of group B.}

\item{batches}{Logical suggesting if replicates of group A and B were performed in batches and are labelled accordingly. If TRUE, a heterogeneous/homogeneous subset may not have multiple samples from the same batch.}

\item{between_combs}{Dataframe with each column containing groupings of replicates of groups A and B, which will be used to assess between-group variation.}

\item{within_combs}{Data.frame with each column containing groupings of replicates of groups A or B, which will be used to assess within-group variation.}
}
\value{
List of two dataframes, containing groupings for within-group and between-group variation.
}
\description{
Identifies subgroupings of replicates for assessing within-group and between-group variation.
}
