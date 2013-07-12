\name{gsymeigen}
\alias{gsymeigen}

\title{
Returns the eigenvalues/eigenvectors of a symetric matrix.
}
\description{
This function is like the function \code{\link{eigen}},
 but it only works on symetric matrices. No attemt is made to check if the matrix is symetric.
}
\usage{
gsymeigen(x)
}

\arguments{
  \item{x}{
 A square matrix usually stored as an object of class  \code{"gmatrix"}.
}
}
\details{
No attempt to check if the matrix is symetric is made.
}
\value{
Returns and object of class \code{"geigen"}.
}

\author{
Nathan Morris
}

\seealso{
\code{\link{geigen-class}}
}
\examples{
	A=gmatrix(grnorm(10*10),10,10)
	B = crossprod(A)
	geig = gsymeigen(B)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
