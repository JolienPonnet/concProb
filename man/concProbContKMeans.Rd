\name{concProbContKMeans}

\alias{concProbContKMeans}

\title{Concordance probability based on k-means clustering in continuous setting}

\description{This function calculates the concordance probability based on the k-means clustering algorithm in a continuous setting.}

\usage{concProbContKMeans(inputDT, nu = 0, nClus, letsTime = TRUE)}

\arguments{
\item{inputDT}{data.table object consisting of the columns 'observed' and 'predicted'. This is an obligatory argument, without default value.}

\item{nu}{Numeric vector of length 1 referring to the value that the observations have to differ at least from one another in order to be taken into account by the algorithm. Default value is 0.}

\item{nClus}{Numeric integer indicating the number of clusters.}

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
library(MASS)
library(data.table)

nSamples <- 10000
bvn1 <- mvrnorm(nSamples, mu =  c(0,0), Sigma = matrix(c(1, 0.75, 0.75, 1), 2) )
inputDT <- as.data.table(data.frame(observed = bvn1[,1], predicted = bvn1[,2]))

concProbContKMeans(inputDT, nClus = 10, letsTime = FALSE)$concProbGlobal
concProbContKMeans(inputDT, nu = 0.4, nClus = 100)$time
concProbContKMeans(inputDT, nu = 0.2, nClus = 50)
}
