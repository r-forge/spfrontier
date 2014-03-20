\name{genW}
\alias{genW}
\title{Standard spatial contiguity matrixes}
\usage{
genW(n, type = "rook")
}
\arguments{
  \item{n}{a number of objects with spatial interaction to
  be arranged.See 'Details' for objects arranging
  principle}

  \item{type}{an optional type of spatial interaction.
  Currently 'rook' and ''queen' values are supported, to
  produce Rook and Queen Contiguity matrix. See references
  for more info. By default set to rook.}
}
\description{
\code{genW} generates an spatial contiguity matrix (rook or
queen)
}
\details{
To generate spatial interaction between \code{n} objects
the function arranges them on a chess board. A number of
columns is calculated as a square root of \code{n}, rounded
to the top. The last row contains empty cells, if n is not
quadratic
}
\examples{
#Completely filled 10x10 rook contiguity matrix
rookW <- genW(100)
rookW

#Partly filled 10x10 rook contiguity matrix
rookW <- genW(90)
rookW

#' #Completely filled 10x10 queen contiguity matrix
rookW <- genW(90)
rookW
}
\references{
Anselin, L. (1988). Spatial Econometrics: Methods and
Models. Kluwer Academic Publishers, Dordrecht, The
Netherlands.
}
\keyword{frontier}
\keyword{spatial}
\keyword{stochastic}

