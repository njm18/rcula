# TODO: Add comment
# 
# Author: nmorris
###############################################################################

setMethod("solve", signature(a = "gmatrix", b = "ANY"),
		function (a, b)
		{
			if(a@type>1L)
				type(a)=0L
			if(!missing(b))
				if(b@type>1L)
					type(b)=type(a)
			#browser()
			a <- qr(a)
			nc <- ncol(a@qr)
			#if (a$rank != nc) 
			#	stop("singular matrix 'a' in 'solve'")
			if (missing(b)) {
				if (nc != nrow(a@qr)) 
					stop("only square matrices can be inverted")
				b <- gident(nc, type=a@qr@type)
				colnames(b) <- rownames(a@qr)
			}
			return(qr.coef(a, b))
		}
)


setMethod("solve", signature(a = "gqr", b = "ANY"),
		function (a, b, ...) 
		{
			nc <- ncol(a@qr)
			nr <- nrow(a@qr)
#			if (a$rank != min(nc, nr)) 
#				if (a$rank != nc) 
#					stop("singular matrix 'a' in 'solve'")
			if (missing(b)) {
				if (nc != nr) 
					stop("only square matrices can be inverted")
				b <- gident(nc, type=a@qr@type)
				colnames(b) <- rownames(a@qr)
			}
			res <- qr.coef(a, b)
			#res[is.na(res)] <- 0
			res
		}
)

#slv = function (a, b, useQr=FALSE, skipCheck=FALSE)
#{
#	if(a@type>1L)
#		type(a)=0L
#	if(!missing(b))
#		if(b@type>1L)
#			type(b)=type(a)
#	browser()
#	if(useQr) {
#		a <- qr(a)
#		nc <- ncol(a$qr)
#		if (a$rank != nc) 
#			stop("singular matrix 'a' in 'solve'")
#		if (missing(b)) {
#			if (nc != nrow(a$qr)) 
#				stop("only square matrices can be inverted")
#			b <- diag(1, nc)
#			colnames(b) <- rownames(a$qr)
#		}
#		return(qr.coef(a, b))
#	} else {
#		if (missing(b)) {
#			if ( ncol(a) != nrow(a)) 
#				stop("only square matrices can be inverted")
#			dupb <- gident(ncol(a))
#		} else
#			dupb=gdup(b)
#		#static SEXP modLa_dgesv(SEXP A_in, SEXP B_in, SEXP ipiv_in)
#		dupa=gdup(a)
#		browser()
#		hperm=.Call("rcula_dgesv",dupa, dupb)
#		if(!skipCheck) {
#			ha=h(a)
#			hdupa=h(dupa)
#			.Call("check_inverse_condition",ha, hdupa,hperm, tol)
#		}
#		return(dupb)
#	}
#	
#}
