\name{culaTest}
\alias{culaTest}
\title{
Test the rcula installation.
}
\description{
This function will error out if the GPU and CPU versions of the operations differ by some tolerance level.
}
\usage{
culaTest(n = 100, tol = 1e-10)
}

\arguments{
  \item{n}{
  Size of the matrix to simulate for testin.
}
  \item{tol}{
  Tolerance lever for differences between the CPU and GPU versions of various operations
}
}

\author{
Nathan Morris
}

\examples{
gtest()
}

