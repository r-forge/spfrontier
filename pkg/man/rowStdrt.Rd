\name{rowStdrt}
\alias{rowStdrt}
\title{Spatial contiguity matrix standartisation}
\usage{
rowStdrt(W)
}
\arguments{
  \item{W}{a spatial contiguity matrix to be standatised}
}
\description{
\code{rowStdrt} standartizes spatial contiguity matrix by
rows
}
\details{
The function divides every element in an argument matrix by
the sum of elements in its row. Some spatial estimation
requires this standartisation (generally - for faster
calculations)
}
\examples{
# Completely filled 10x10 queen contiguity matrix
queenW <- genW(100, type="queen")
queenW

# Standartisation
stQueenW <- rowStdrt(queenW)
stQueenW
}

