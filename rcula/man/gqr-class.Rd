\name{gqr-class}
\Rdversion{1.1}
\docType{class}
\alias{gqr-class}

\alias{qr.coef}
\alias{qr.coef-methods}
\alias{qr.coef,ANY-method}
\alias{qr.coef,gqr-method}
\alias{solve}
\alias{solve-methods}
\alias{solve,ANY-method}
\alias{solve,gqr-method}


\title{Class \code{"gqr"}}
\description{
Conatians the results of a QR decomposition from the \code{qr()} method. Note that the \code{solve} function uses QR decpompostion. Also, \code{qr.coef} acts on the \code{gqr} object.
}
\section{Objects from the Class}{
Objects can be created by calls of the \code{qr()} method.
}
\section{Slots}{
  \describe{
    \item{\code{qr}:}{Information about the Q and R matrices stored as and bject of class \code{"gmatrix"}.}
    \item{\code{qraux}:}{Information about the Q matrix stored as an object of class \code{"gvector"}. }
  }
}
\section{Methods}{
  \describe{
    \item{qr.coef}{Simular to the \code{\link{qr.coef}} function. }
    \item{solve}{Simular to the \code{\link{solve.qr}} function. }
	 }
}

\author{
Nathan Morris
}


\seealso{
\code{\link{qr}}
}

\keyword{classes}
