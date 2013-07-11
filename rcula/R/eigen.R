
setClass("geigen",
		representation(
				vectors = "gmatrix",
				values  = "gvector"),
		prototype = list(
				vectors = NULL,
				values  = NULL)
) 


symeigen = function (x) 
{
	x <- as.gmatrix(x,dup=FALSE)
	if(x@type>1L)
		type(x)=0L
	x =gnamestrip(x)
	n <- nrow(x)
	if (!n) 
		stop("0 x 0 matrix")
	if (n != ncol(x)) 
		stop("non-square matrix in 'eigen'")
	
	if (sum(!is.finite(x),retgpu=FALSE)>0) 
		stop("infinite or missing values in 'x'")
	#rcula_eigen_nonsymm(SEXP A_in, SEXP val_in)
	
	ret=new("geigen", vectors=gdup(x), values=gvector(n))
	z <- .Call("rcula_eigen_symm", ret@vectors, ret@values)

	ord <- order(h(ret@values), decreasing = TRUE)
	ret@values =ret@values[ord]
	ret@vectors=ret@vectors[,ord, drop=FALSE]
	return(ret)
}

#eigen = function (x, symmetric, only.values = FALSE) 
#{
#	x <- as.gmatrix(x,dup=FALSE)
#	if(x@type>1L)
#		type(x)=0L
#	x =gnamestrip(x)
#	n <- nrow(x)
#	if (!n) 
#		stop("0 x 0 matrix")
#	if (n != ncol(x)) 
#		stop("non-square matrix in 'eigen'")
#	
#	if (sum(is.finite(x),retgpu=FALSE)>0) 
#		stop("infinite or missing values in 'x'")
#	if (missing(symmetric)) 
#		symmetric <- isSymmetric.matrix(x)
#	
#	dbl.n <- double(n)
#	if (symmetric) {
#		
#		z <- .Fortran("rs", n, n, x, values = dbl.n, !only.values, 
#				vectors = x, dbl.n, dbl.n, ierr = integer(1L), 
#				PACKAGE = "base")
#		if (z$ierr) 
#			stop(gettextf("'rs' returned code %d in 'eigen'", 
#							z$ierr), domain = NA)
#		
#		ord <- sort.list(z$values, decreasing = TRUE)
#	}
#	else {
#		
#		z <- .Fortran("rg", n, n, x, values = dbl.n, ivalues = dbl.n, 
#				!only.values, vectors = x, integer(n), dbl.n, 
#				ierr = integer(1L), PACKAGE = "base")
#		if (z$ierr) 
#			stop(gettextf("'rg' returned code %d in 'eigen'", 
#							z$ierr), domain = NA)
#		ind <- z$ivalues > 0L
#		if (any(ind)) {
#			ind <- seq.int(n)[ind]
#			z$values <- complex(real = z$values, imaginary = z$ivalues)
#			if (!only.values) {
#				z$vectors[, ind] <- complex(real = z$vectors[, ind], imaginary = z$vectors[, ind + 1])
#				z$vectors[, ind + 1] <- Conj(z$vectors[, ind])
#			}
#		}
#		
#		ord <- sort.list(Mod(z$values), decreasing = TRUE)
#	}
#	list(
#			values  = z$values[ord],
#			vectors = if (!only.values) z$vectors[, ord, drop = FALSE])
#}