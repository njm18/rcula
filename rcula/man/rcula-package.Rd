\name{rcula-package}
\alias{rcula-package}
%\alias{rcula}
\docType{package}
\title{
A conection between are and CULA.
}
\description{
This package implements a few simple functions to make function from the CULA library for performing matrix decompostions on a GPU.
}
\details{
	\tabular{ll}{
		Package: \tab rcula\cr
		Type: \tab Package\cr
		Version: \tab 0.1\cr
		Date: \tab 2013-7-13\cr
		License: \tab GPL-2\cr
		LazyLoad: \tab yes\cr
		Depends: \tab methods, gmatrix\cr
	}
This package is fairly simple. It implements what is essentually base R functiontionality for the function \code{solve},
\code{qr} and \code{qr.coef}. Note that unlike R these functions do not check the rank of the matrix, and solve will not 
through and error for ill conditioned matrices. We use QR decomposition for the \code{solve} functionality. Only symetric eigen
decomposition is implemented here with the function \code{\link{gsymeigen}}. 
}
\author{
Nathan Morris
}
\references{
None
}

\keyword{NVIDIA}
\keyword{GPU}

\seealso{
 \code{\link{gqr-class}}, \code{\link{gsymeigen}}, \code{\link{geigen-class}}, \code{\link{culaTest}}
}
\examples{
#solve examples
A=gmatrix(grnorm(400),20,20)
b=grnorm(20)
Ainverse = solve(A)
z=solve(A,b)

#qr
B=qr(A)
Ainverse2=solve(B)

#symetric eigen
C=crossprod(A)
Ceig=gsymeigen(C)
vals = h(values(Ceig))
vecs = h(vectors(Ceig))

}
