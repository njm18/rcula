# TODO: Add comment
# 
# Author: nmorris
###############################################################################


solve.gmatrix = function (a, b, useQr=FALSE, skipCheck=FALSE, tol=1E-7 )
{
	tol=as.numeric(tol)[1]
	if(a@type>1L)
		type(a)=0L
	if(!missing(b))
		if(b@type>1L)
			type(b)=type(a)
	
	if(useQr) {
		a <- qr(a, tol = tol)
		nc <- ncol(a$qr)
		if (a$rank != nc) 
			stop("singular matrix 'a' in 'solve'")
		if (missing(b)) {
			if (nc != nrow(a$qr)) 
				stop("only square matrices can be inverted")
			b <- diag(1, nc)
			colnames(b) <- rownames(a$qr)
		}
		return(qr.coef(a, b))
	} else {
		if (missing(b)) {
			if ( ncol(a) != nrow(a)) 
				stop("only square matrices can be inverted")
			dupb <- gident(ncol(a))
		} else
			dupb=gdup(b)
		#static SEXP modLa_dgesv(SEXP A_in, SEXP B_in, SEXP ipiv_in)
		dupa=gdup(a)
		browser()
		hperm=.Call("rcula_dgesv",dupa, dupb)
		if(!skipCheck) {
			ha=h(a)
			hdupa=h(dupa)
			.Call("check_inverse_condition",ha, hdupa,hperm, tol)
		}
		return(dupb)
	}

}


solve.qr = function (a, b, ...) 
{
	if (!is.qr(a)) 
		stop("this is the \"qr\" method for the generic function solve()")
	nc <- ncol(a$qr)
	nr <- nrow(a$qr)
	if (a$rank != min(nc, nr)) 
		if (a$rank != nc) 
			stop("singular matrix 'a' in 'solve'")
	if (missing(b)) {
		if (nc != nr) 
			stop("only square matrices can be inverted")
		b <- diag(1, nc)
	}
	res <- qr.coef(a, b)
	res[is.na(res)] <- 0
	res
}