\name{concProbKMeans}

\alias{concProbKMeans}

\title{Concordance probability based on k-means clustering in discrete setting}

\description{This function calculates the concordance probability based on the k-means clustering algorithm in a discrete setting.}

\usage{concProbKMeans(inputDT, nClusMax, letsTime = TRUE)}

\arguments{
\item{inputDT}{data.table object consisting of the columns 'observed' and 'predicted'. This is an obligatory argument, without default value.}

\item{nClusMax}{Numeric integer indicating the maximum number of clusters.}

\item{letsTime}{Logical vector of length 1 indicating whether or not the computing time of this function should be measured. Default value is TRUE.}
}

\value{
\code{concProbContKMeans} returns a list that contains the following items:
\item{concProbGlobal}{the concordance probability, obtained following the supplied values for the arguments}
\item{time}{if \code{letsTime} is TRUE, the computing time of this function}
}

\references{Van Oirbeek, R., Ponnet, J., Verdonck, T. (2021). Concordance Probability in a Big Data Setting.}

\author{Robin Van Oirbeek, Jolien Ponnet and Tim Verdonck}

\examples{
library(concProb)
library(data.table)

nSamples <- 10000
obsValues <- rbinom(nSamples, 1, 0.75)
inputDT <- as.data.table(data.frame(observed = obsValues, 
 predicted = runif(nSamples)))
 
concProbKMeans(inputDT, 30, letsTime = FALSE)$concProbGlobal
concProbKMeans(inputDT, 20)$time
concProbKMeans(inputDT, 10)
}
