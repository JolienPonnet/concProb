\name{concProbGrid}

\alias{concProbGrid}

\title{Concordance probability based on marginal approximation in discrete setting}

\description{This function calculates the concordance probability based on a marginal approximation in a discrete setting.}

\usage{concProbGrid(inputDT, quantSplits, letsTime = TRUE)}

\arguments{
\item{inputDT}{data.table object consisting of the columns 'observed' and 'predicted'. This is an obligatory argument, without default value.}

\item{quantSplits}{Numeric vector representing the quantiles to define the boundary values of the grid.}

\item{letsTime}{Logical vector of length 1 indicating whether or not the computing time of this function should be measured. Default value is TRUE.}
}

\value{
\code{concProbContGrid} returns a list that contains the following items:
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
 
concProbGrid(inputDT, quantSplits = seq(0.1,0.9,0.1))
concProbGrid(inputDT, quantSplits = 0.5)$concProbGlobal
concProbGrid(inputDT, quantSplits = seq(0.1,0.9,0.1), letsTime = TRUE)$time
}
