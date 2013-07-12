# TODO: Add comment
# 
# Author: nmorris
###############################################################################


culaTest = function(n=100, tol= 1e-10) {
	
	myCheck = function(d, msg, tol) {
		if(max(abs(d))>tol)
			stop(paste("Testing dectected problems with the", msg))
	}
	A=gmatrix(grnorm(n*n),n,n)
	b=grnorm(n)
	
	#test solve (solve calls qr and qr.coef - so we skip testing those)
	myCheck(h(solve(A))-solve(h(A)), "solve(a) function.", tol)
	myCheck(h(solve(A,b))-solve(h(A),h(b)), "solve(a,b) function.", tol)
	
	#test eigen
	B = crossprod(A)
	geig = gsymeigen(B)
	heig = eigen(B)
	
	myCheck(h(values(geig))-heig$val, "gsymeig() function.", tol)
	vecg=h(vectors(geig))
	vech=heig$vec
	d=vecg-vech
	d=colSums(d*d)
	vecg = t( ifelse(d>.1, -1, 1) * t(vecg)) #account for the fact that eigenvec is only spcified upt to +/-1 factor
	myCheck(vecg-vech, "gsymeig() function.", tol)
	return("No Problems Detected")
}