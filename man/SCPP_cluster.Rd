\name{SCPP_cluster}
\alias{SCPP_cluster}
\title{Spectral Clustering using Projection Pursuit}
\usage{
SCPP_cluster(X, K, v0, ndim, nMicro, betamax, betamin, smult, minsize, minprop, omega, type)
}
\description{
Generates a binary partitioning tree by recursively partitioning a dataset
using spectral clustering within an optimal subspace.
}
\arguments{
  \item{X}{a numeric matrix (num_data x num_dimensions); the dataset to be clustered.}
  \item{K}{the number of clusters to extract.}
  \item{v0}{(optional) initialisations for projection puprsuit. a function(X) of the data being split, which returns a matrix with ncol(X)*ndim rows. each column of the output of v0(X) is used as an initialisation for projection pursuit. That is, the i-th initialisation is via the projection matrix: matrix(v0[,i], ncol = ndim). the solution with the minimum spectral connectivity is used within the final model. initialisations are determined separately for each cluster being split. if omitted then a single initialisation is used; the first ndim principal components.}
  \item{ndim}{(optional) dimension of subspace.if omitted then ndim = 2}
  \item{nMicro}{(optional) number of microclusters. running time is quadratic in the number of microclusters. approximation of true objective is better the more microclusters used. if omitted then nMicro = 200}
  \item{betamax}{(optional) initial value of beta, which affects pairwise similarities by increasing the similarity of points outside beta standard deviations from the mean. This reduces the effect of outliers. if omitted then betamax = 5.}
  \item{betamin}{(optional) smallest value of beta considered. beta is reduced by 0.1 from initial value until a desired balance in cluster size is met. if omitted then betamin = 0.5.}
  \item{smult}{(optional) multiplicative factor applied to scaling parameter used in pairwise similarities. scaling parameter is determined for each cluster (C) being split as sqrt(mean(eigen(cov(C))$values[1:intr]))*smult*(4/3/n)^(1/(4+intr)), where intr is an estimate of the intrinsic dimensionality of $C$. if omitted then smult = 1.}
  \item{minsize}{(optional) the minimum cluster size allowable. if omitted then minsize = nrow(X)/K/5, i.e., 20 percent of the average cluster size.}
  \item{minprop}{(optional) the minimum cluster expressed as a proportion of the cluster being split. if omitted then minimum cluster size is detemined by minsize.}
  \item{omega}{(optional) parameter controlling the orthogonality of the columns of the projection vector. if omitted then omega = 1.}
  \item{type}{(optional) type of Laplacian ("standard" or "normalised") to use. if omitted then type = "normalised".}
}
\value{
  a named list containing
  \item{$cluster}{cluster assignment vector.}
  \item{$model}{matrix containing the would-be location of each node (depth and position at depth) within a complete tree of appropriate depth.}
  \item{$nodes}{unnamed list each element of which is a named list containing details of the binary partitions at each node in the model.}
  \item{$data}{the data matrix being clustered.}
  \item{$args}{named list of arguments passed to SCPP_cluster.}
}
\references{
  Hofmeyr, D., Pavlidis, N., Eckley, I. (2018) Minimum spectral connectivity projection pursuit. \emph{Statistics and Computing}, to appear.
}
\examples{
## load synthetic control chart dataset
data(chart)

## obtain clustering solution using SCPP
sol <- SCPP_cluster(chart$x, 6)

## plot cluster model
SCPP_plot(sol)

## evaluate performance using external cluster validity metrics
cluster_performance(sol$cluster, chart$c)
}
