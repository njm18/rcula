\name{geigen-class}
\Rdversion{1.1}
\docType{class}
\alias{geigen-class}
\alias{values-methods}
\alias{values,geigen-method}
\alias{vectors-methods}
\alias{vectors,geigen-method}

\title{Class \code{"geigen"}}
\description{
Contains eigenvectors and eigenvalues of a matrix.
}
\section{Objects from the Class}{
Objects of this class are greated by a call to \code{\link{gsymeigen}}
}
\section{Slots}{
  \describe{
    \item{\code{vectors}:}{The eigenvectors stored as the columns of an object of class \code{"gmatrix"}. }
    \item{\code{values}:}{The eigenvalues stored as and object of class \code{"gvector"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{values}{\code{signature(a = "geigen")}: extracts the eigenvalues }
    \item{vectors}{\code{signature(a = "geigen")}: extracts the eigenvectors}
	 }
}

\author{
	Nathan Morris
}

\seealso{
\code{\link{gsymeigen}}, 
}

\keyword{classes}
