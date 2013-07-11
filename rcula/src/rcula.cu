
#include "rcula.h"
#include <R_ext/Lapack.h>


struct gmat get_gmat_struct(SEXP A_in) {
	struct gmat A;
	struct gptr *gpu_ptr = (struct gptr*) R_ExternalPtrAddr(GET_SLOT(A_in, install("ptr")));
	A.d_vec = gpu_ptr->d_vec;
	A.nrow = INTEGER(GET_SLOT(A_in, install("nrow")))[0];
	A.ncol = INTEGER(GET_SLOT(A_in, install("ncol")))[0];
	A.type = INTEGER(GET_SLOT(A_in, install("type")))[0];
	A.device = INTEGER(GET_SLOT(A_in, install("device")))[0];
	return A;
}

struct gvec get_gvec_struct(SEXP A_in) {
	struct gvec A;
	struct gptr *gpu_ptr = (struct gptr*) R_ExternalPtrAddr(GET_SLOT(A_in, install("ptr")));
	A.d_vec = gpu_ptr->d_vec;
	A.length = INTEGER(GET_SLOT(A_in, install("length")))[0];
	A.type = INTEGER(GET_SLOT(A_in, install("type")))[0];
	A.device = INTEGER(GET_SLOT(A_in, install("device")))[0];
	return A;
}

/*
struct matrix create_gvector(SEXP a) {
	struct matrix A;
	struct gpuvec *gpu_ptr = (struct gpuvec*) R_ExternalPtrAddr(GET_SLOT(A_in, install("ptr")));
	A.d_vec = gpu_ptr->d_vec;
	A.rows = INTEGER(GET_SLOT(A_in, install("nrow")))[0];
	A.cols = INTEGER(GET_SLOT(A_in, install("ncol")))[0];
	return A;
}



SEXP creat_gmatrix(struct matrix A) {
	SEXP ret_final;
	PROTECT(ret_final = NEW_OBJECT(MAKE_CLASS("gmatrix")));
	SET_SLOT(ret_final, install("nrow"), PROTECT(asSEXPint(A.rows)));
	SET_SLOT(ret_final, install("ncol"), PROTECT(asSEXPint(A.cols)));
	SET_SLOT(ret_final, install("ptr"), PROTECT(gpu_register(A.d_vec)));
	SET_SLOT(ret_final, install("type"), PROTECT(asSEXPint(A.type)));
	SET_SLOT(ret_final, install("device"), PROTECT(asSEXPint(A.device)));
	UNPROTECT(8L);
	return A;
}*/



void call_error(culaStatus s) {
	int info;
	char buf[256];
	info = culaGetErrorInfo();
	culaGetErrorInfoString(s, info, buf, sizeof(buf));
	error("Culbas returned error: %s", buf);
}

SEXP asSEXPint(int myint) { //wraps an integer as a sexp
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = myint;
	UNPROTECT(1);
	return ans;
}

SEXP asSEXPreal(double myint) { //wraps an integer as a sexp
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, 1));
	REAL(ans)[0] = myint;
	UNPROTECT(1);
	return ans;
}

SEXP rcula_initialize() {
	culaStatus s;
	s=culaInitialize();
	CHECK_CULA(s);
    return asSEXPint(1L);
}

SEXP rcula_shutdown() {
	culaShutdown();
    return asSEXPint(1L);
}


SEXP rcula_qr(SEXP A_in, SEXP qraux_in)
{
    int m, n;
    culaStatus s;


	struct gmat A = get_gmat_struct(A_in);
	struct gptr *qraux = (struct gptr*) R_ExternalPtrAddr(qraux_in);

	if (A.type >1L)
		error("'a' must be a of type 'double' or 'single'");
    //PROTECT(A = duplicate(Ain));
    //Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    m = A.nrow;//Adims[0];
    n = A.ncol;//Adims[1];

    if(A.type==0L)
    	s = culaDeviceDgeqrfp(m, n, (double *) A.d_vec, m,(double *) qraux->d_vec);
    else
    	s = culaDeviceSgeqrfp(m, n, (float *) A.d_vec, m,(float *)qraux->d_vec);
    CHECK_CULA(s);

    return asSEXPint(1L);
}



SEXP rcula_modqr_coef(SEXP qr_in, SEXP qraux_in, SEXP B_in)
{
    int n, nrhs, k;
    culaStatus s;
	struct gmat qr = get_gmat_struct(qr_in);
	struct gmat B = get_gmat_struct(B_in);
	struct gvec qraux = get_gvec_struct(qraux_in);



    k =qraux.length;
//    if (!(isMatrix(Bin) && isReal(Bin)))
//	error(_("'b' must be a numeric matrix"));

//    PROTECT(B = duplicate(Bin));
//    Qdims = INTEGER(coerceVector(getAttrib(qr, R_DimSymbol), INTSXP));
    n = qr.nrow;
  //  Bdims = INTEGER(coerceVector(getAttrib(B, R_DimSymbol), INTSXP));
    if(B.nrow != n)
    	error("right-hand side should have %d not %d rows", n,B.nrow );
    nrhs = B.ncol;
    if(qr.type==0L)
    	s=culaDeviceDormqr('L', 'T', n, nrhs, k,
        		(double *) qr.d_vec, n, (double *)qraux.d_vec, (double *) B.d_vec, n);
    else
    	s=culaDeviceSormqr('L', 'T', n, nrhs, k,
        		(float *) qr.d_vec, n, (float *)qraux.d_vec, (float *) B.d_vec, n);
    CHECK_CULA(s);

    if(qr.type==0L)
    	s=culaDeviceDtrtrs('U', 'N', 'N', k, nrhs,
    			(double *) qr.d_vec, n, (double *) B.d_vec, n);
    else
    	s=s=culaDeviceStrtrs('U', 'N', 'N', k, nrhs,
    			(float *) qr.d_vec, n, (float *) B.d_vec, n);
    CHECK_CULA(s);

    return asSEXPint(1L);
}


SEXP rcula_dgesv(SEXP A_in, SEXP B_in)
{
    int n;
    culaStatus s;

	struct gmat A = get_gmat_struct(A_in);
	struct gmat B = get_gmat_struct(B_in);


    if (A.type>1L)
    	error(("'a' must be of type 'single' or 'double.'"));
    if (B.type>1L)
    	error(("'b' must be of type 'single' or 'double.'"));

    n = A.nrow;
    if(n == 0) error(("'a' is 0-diml"));
    int p = B.ncol;
    if(p == 0) error(("no right-hand side in 'b'"));
    if(A.ncol != n)
    	error(("'a' (%d x %d) must be square"), n, A.ncol);
    if(B.nrow != n)
    	error(("'b' (%d x %d) must be compatible with 'a' (%d x %d)"), B.nrow, p, n, n);
    SEXP ret;
    PROTECT(ret = allocVector(INTSXP, n));
    int *ipiv = INTEGER(ret);
    for(int tmp=0;tmp<n;tmp++)
    	ipiv[tmp]=1;
    if(A.type==0L)
    	s = culaDgesv(n, p, (double *)A.d_vec, n, ipiv, (double *)B.d_vec, n);
    else
    	s = culaSgesv(n, p, (float *) A.d_vec, n, ipiv, (float *)B.d_vec, n);

    CHECK_CULA(s);
    UNPROTECT(1);
    return ret;
}


SEXP check_inverse_condition(SEXP Ain, SEXP Avalsin, SEXP permin, SEXP tolin) {
	double tol = REAL(tolin)[0];
	double *A = REAL(Ain);
	double *Avals= REAL(Avalsin);
	int *ipiv = INTEGER(permin);
	double rcond;
	int info;

	int *Adims = INTEGER(coerceVector(getAttrib(Ain, R_DimSymbol), INTSXP));
    int n = Adims[0];
	double anorm = F77_CALL(dlange)("1", &n, &n, A, &n, (double*) NULL);
	double *work = (double *) R_alloc(4*n, sizeof(double));
	F77_CALL(dgecon)("1", &n, Avals, &n, &anorm, &rcond, work, ipiv, &info);
	if (rcond < tol)
		error(("system is computationally singular: reciprocal condition number = %g"),
				rcond);
	return(asSEXPreal(rcond));

}

SEXP rcula_eigen_symm(SEXP A_in, SEXP val_in)
{
	culaStatus s;
	struct gmat A = get_gmat_struct(A_in);
	struct gvec val = get_gvec_struct(val_in);

	if(A.nrow != val.length)
		error("dim mismatch");

	if(A.type==0L)
		s=culaDeviceDsyev('V', 'U', A.ncol, (double *) A.d_vec , A.nrow, (double *) val.d_vec );
	else
		s=culaDeviceSsyev('V', 'U', A.ncol, (float *) A.d_vec , A.nrow, (float *) val.d_vec );


	CHECK_CULA(s);

	return asSEXPint(1L);
}

/*
SEXP rcula_eigen_nonsymm(SEXP A_in, SEXP val_in)
{
    culaStatus s;
	struct gmat A = get_gmat_struct(A_in);
	struct gvec xaux = get_gvec_struct(val_in);
//	int symm=INTEGER(symm_in)[0];

	if(A.nrow != val.length)
		error("dim mismatch");


	if(A.type==0L)
		s=culaDeviceDsgeev(�V�,�U�, A.ncol,
				(double *) A.d_vec , A.nrow, (double *) val.d_vec );
	else
		s=culaDeviceSsgeev(�N�,�V�, A.ncol,
				(float *) A.d_vec , A.nrow, (float *) val.d_vec );

	CHECK_CULA(s);
	return asSEXPint(1L);
}
//culaDeviceSsyev(�V�,�U�, cols, ptr(stores matrix/eigvec out), rows, prt(vals))
*/
