

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <math.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Arith.h>
#include <R_ext/Boolean.h>
//#include <cutil_inline.h>
//include "cublas.h"
#include <stdio.h>
#include <cublas_v2.h>
#include <curand_kernel.h>
#include <cula.h>



#define CHECK_CULA(s) \
	if( s != culaNoError )\
	call_error(s);\


#define IDX2(i ,j ,ld) (((j)*(ld))+(i))


#define GET_BLOCKS_PER_GRID(n)  \
	int blocksPerGrid = (n + (threads_per_block[currentDevice]) - 1) / (threads_per_block[currentDevice]); \
	int operations_per_thread = 1;  \
	if(blocksPerGrid>MAX_BLOCKS) {  \
		blocksPerGrid = MAX_BLOCKS;  \
		int total_threads = blocksPerGrid*(threads_per_block[currentDevice]); \
		operations_per_thread = (n + total_threads -1) / total_threads; \
	}


#define DECERROR0 cudaError_t  cudaStat
#define DECERROR1 cudaError_t  cudaStat, status1
#define DECERROR2 cudaError_t  cudaStat, status1, status2
#define DECERROR3 cudaError_t  cudaStat, status1, status2, status3

//Macros for malloc
#define CUDA_MALLOC(MPTR,MN)  \
		cudaStat = cudaMalloc( (void **)&(MPTR),MN) ;\
		if (cudaStat != cudaSuccess ){\
			R_gc();\
			cudaStat = cudaMalloc( (void **)&(MPTR),MN) ;\
			if (cudaStat != cudaSuccess ){\
				error("CUDA memory allocation error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
			}\
		}

#define CUDA_MALLOC_CLEAN_1(MPTR,MN,MCLEANPTR)  \
		cudaStat = cudaMalloc( (void **)&(MPTR),MN) ;\
		if (cudaStat != cudaSuccess ) {\
			R_gc();\
			cudaStat = cudaMalloc( (void **)&(MPTR),MN) ;\
			if (cudaStat != cudaSuccess ){\
				status1=cudaFree(MCLEANPTR);\
				if (status1 != cudaSuccess) {\
					error("CUDA memory allocation and free error (potential memory leak) in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
				}\
				error("CUDA memory allocation error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
			}\
		}

#define CUDA_MALLOC_CLEAN_2(MPTR,MN,MCLEANPTR1,MCLEANPTR2)  \
		cudaStat = cudaMalloc( (void **)&(MPTR),MN ) ;\
		if (cudaStat != cudaSuccess ) {\
			R_gc();\
			cudaStat = cudaMalloc( (void **)&(MPTR),MN) ;\
			if (cudaStat != cudaSuccess ){\
				status1=cudaFree(MCLEANPTR1);\
				status2=cudaFree(MCLEANPTR2);\
				if (status1 != cudaSuccess || status2 != cudaSuccess) {\
					error("CUDA memory allocation error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
				}\
				error("CUDA memory allocation error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
			}\
		}

//Macros for checking errors
#define CUDA_ERROR \
		if (cudaStat != cudaSuccess ) {\
		 error("Error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
		}

//Macros for checking the kernal and cleaning up if there are errors
#define CUDA_CHECK_KERNAL  \
		cudaStat = cudaDeviceSynchronize(); \
		if (cudaStat != cudaSuccess ) {\
			error("Kernal error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
		}

#define CUDA_CHECK_KERNAL_CLEAN_1(MCLEANPTR1)  \
		cudaStat = cudaDeviceSynchronize(); \
		if (cudaStat != cudaSuccess ) {\
			status1=cudaFree(MCLEANPTR1);\
			if (status1 != cudaSuccess ) {\
				error("Kernal error and memory free errors (potential memory leak) in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
			}\
		 error("Kernal error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
		}
#define CUDA_CHECK_KERNAL_CLEAN_2(MCLEANPTR1,MCLEANPTR2)  \
		cudaStat = cudaDeviceSynchronize(); \
		if (cudaStat != cudaSuccess ) {\
			status1=cudaFree(MCLEANPTR1);\
			status2=cudaFree(MCLEANPTR2);\
			if (status1 != cudaSuccess || status2 != cudaSuccess) {\
				error("Kernal error and memory free errors (potential memory leak) in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
			}\
		 error("Kernal error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
		}

#define CUDA_CHECK_KERNAL_CLEAN_3(MCLEANPTR1,MCLEANPTR2,MCLEANPTR3)  \
		cudaStat = cudaDeviceSynchronize(); \
		if (cudaStat != cudaSuccess ) {\
			status1=cudaFree(MCLEANPTR1);\
			status2=cudaFree(MCLEANPTR2);\
			status3=cudaFree(MCLEANPTR3);\
			if (status1 != cudaSuccess || status2 != cudaSuccess || status3 != cudaSuccess ) {\
				error("Kernal error and memory free errors (potential memory leak) in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
			}\
		 error("Kernal error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
		}


//macros for cleaning up
#define CUDA_CLEAN_1(MCLEANPTR1)  \
		status1=cudaFree(MCLEANPTR1);\
		if (status1 != cudaSuccess) {\
			error("Memory free errors (potential memory leak) in '%s.' (%s)\n", __func__, cudaGetErrorString(status1));\
		}

#define CUDA_CLEAN_2(MCLEANPTR1,MCLEANPTR2)  \
		status1=cudaFree(MCLEANPTR1);\
		status2=cudaFree(MCLEANPTR2);\
		if (status1 != cudaSuccess || status2 != cudaSuccess ) {\
			if (status1 != cudaSuccess && status2 == cudaSuccess )\
				error("Memory free errors (potential memory leak) in '%s.' (%s)\n", __func__, cudaGetErrorString(status1));\
			if (status1 == cudaSuccess && status2 != cudaSuccess )\
				error("Memory free errors (potential memory leak) in '%s.' (%s)\n", __func__, cudaGetErrorString(status2));\
			if (status1 != cudaSuccess && status2 != cudaSuccess )\
				error("Memory free errors (potential memory leak) in '%s.' (%s) (%s)\n", __func__, cudaGetErrorString(status1), cudaGetErrorString(status2));\
		}


//memory copy macros (this assumes that DST is on the device and needs to be cleaned up when errors arise)
#define CUDA_MEMCPY_CLEAN(DST, SRC, COUNT,KIND)  \
    cudaStat=cudaMemcpy(DST, SRC, COUNT, KIND) ;\
    if (cudaStat != cudaSuccess) {\
		status1=cudaFree(DST);\
		if (status1 != cudaSuccess  ) {\
			error("Memory copy and memory free errors (potential memory leak) in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
		}\
		error("Memory copy error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
    }

#define CUDA_MEMCPY_CLEAN_1(DST, SRC, COUNT,KIND, MCLEANPTR)  \
    cudaStat=cudaMemcpy(DST, SRC, COUNT, KIND) ;\
    if (cudaStat != cudaSuccess) {\
		status1=cudaFree(DST);\
		status2=cudaFree(MCLEANPTR);\
		if (status1 != cudaSuccess || status2 != cudaSuccess ) {\
			error("Memory copy and memory free errors (potential memory leak) in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
		}\
		error("Memory copy error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
    }

#define CUDA_MEMCPY_CLEAN_2(DST, SRC, COUNT,KIND, MCLEANPTR1, MCLEANPTR2)  \
    cudaStat=cudaMemcpy(DST, SRC, COUNT, KIND) ;\
    if (cudaStat != cudaSuccess) {\
		status1=cudaFree(DST);\
		status2=cudaFree(MCLEANPTR1);\
		status3=cudaFree(MCLEANPTR2);\
		if (status1 != cudaSuccess || status2 != cudaSuccess || status3 != cudaSuccess ) {\
			error("Memory copy and memory free errors (potential memory leak) in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
		}\
		error("Memory copy error in '%s.' (%s)\n", __func__, cudaGetErrorString(cudaStat));\
    }

#define 	PROCESS_TYPE_NO_SIZE\
	int type = INTEGER(in_type)[0];\
	if(type>3)\
		error("Incorrect type passed to '%s.'", __func__);

#define 	PROCESS_TYPE\
	int type = INTEGER(in_type)[0];\
	int mysizeof;\
	if(type==0)\
		mysizeof=sizeof(double);\
	else if(type==1)\
		mysizeof=sizeof(float);\
	else if(type==2 || type==3)\
		mysizeof=sizeof(int);\
	else\
		error("Incorrect type passed to '%s.'", __func__);

#define 	PROCESS_TYPE_SF\
	int type = INTEGER(in_type)[0];\
	int mysizeof;\
	if(type==0)\
		mysizeof=sizeof(double);\
	else if(type==1)\
		mysizeof=sizeof(float);\
	else \
		error("Incorrect type passed to '%s.' Type must be 'double' or 'float.'", __func__);

#define PTR_DBL(A) \
		(double *) A->d_vec

#define PTR_FLOAT(A) \
		(float *) A->d_vec

#define PTR_INT(A) \
		(int *) A->d_vec

#define CALL_KERNAL\
		if(type==0)\
			KERNAL(PTR_DBL, double)\
		else if(type==1)\
			KERNAL(PTR_FLOAT, float)\
		else \
			KERNAL(PTR_INT, int)\

#define CALL_KERNAL_SF\
		if(type==0)\
			KERNAL(PTR_DBL, double)\
		else if(type==1)\
			KERNAL(PTR_FLOAT, float)\


#define MAX_DEVICE 20
#define MAX_BLOCKS 65000

#ifdef DEFINEGLOBALSHERE
#define GLOBAL
#else
#define GLOBAL extern
#endif

GLOBAL __device__ int CUDA_R_Na_int;
GLOBAL __device__ double CUDA_R_Na_double;
GLOBAL __device__ float CUDA_R_Na_float;
GLOBAL cublasHandle_t handle[MAX_DEVICE];
GLOBAL int total_states[MAX_DEVICE] ;
GLOBAL curandState* dev_states[MAX_DEVICE];
GLOBAL int threads_per_block[MAX_DEVICE] ;
GLOBAL int dev_state_set[MAX_DEVICE] ;
GLOBAL int dev_cublas_set[MAX_DEVICE];
GLOBAL int currentDevice ;



union ieee_float {
  unsigned int myint;
  float myfloat;
};
union ieee_double {
  unsigned long mylong;
  double mydouble;
};


#define RNAREAL 0x7ff00000000007A2
#define MYNAFLOAT 0x7F8000FF

//make NA
template <typename T>
__forceinline__ __device__ void MAKE_NA(T *ret) {
	((ieee_double*) ret)->mylong=RNAREAL;
}
template <>
__forceinline__ __device__ void MAKE_NA(float *ret) {
	((ieee_float*) ret)->myint=MYNAFLOAT;
}

template <>
__forceinline__ __device__ void MAKE_NA(int *ret) {
	ret[0]=INT_MIN;
}


//is NA
template <typename T>
__forceinline__ __device__ int IS_NA(T *val) {
	return(((ieee_double*) val)->mylong==RNAREAL);
}
template <>
__forceinline__ __device__ int IS_NA(float *val) {
	return(((ieee_float*) val)->myint==MYNAFLOAT);
}

template <>
__forceinline__ __device__ int IS_NA(int *val) {
	return(INT_MIN==val[0]);
}

//is NAN
template <typename T>
__forceinline__ __device__ int IS_NAN(T val) {
	return(isnan(val));
}
template <>
__forceinline__ __device__ int IS_NAN(float val) {
	return(isnan(val));
}

template <>
__forceinline__ __device__ int IS_NAN(int val) {
	return(INT_MIN==val);
}

//return nan
template <typename T>
__forceinline__ __device__ T RET_NAN(void) {
	return(NAN);
}
template <>
__forceinline__ __device__ float RET_NAN(void) {
	return(NAN);
}

template <>
__forceinline__ __device__ int RET_NAN(void) {
	return(INT_MIN);
}


extern "C" {

struct gptr {
	 void *d_vec;
};


struct gvec {
	 void *d_vec;
	 int length;
	 int type;
	 int device;
};

struct gmat {
   void *d_vec;
   int nrow;
   int ncol;
   int type;
   int device;
};

SEXP rcula_qr(SEXP A_in, SEXP qraux_in);
SEXP rcula_modqr_coef(SEXP qr_in, SEXP qraux_in, SEXP B_in);
SEXP rcula_initialize();
SEXP rcula_shutdown();
SEXP rcula_eigen_symm(SEXP A_in, SEXP val_in);
SEXP rcula_dgesv(SEXP A_in, SEXP B_in);
SEXP check_inverse_condition(SEXP Ain, SEXP Avalsin, SEXP permin, SEXP tolin) ;
}


/*
__global__ void VecAdd(const double* A, const double* B, double* C, int N);
__global__ void myFKernal(const double* bigY, const double* bigX,
		const double* EmpPriorR_x_Sqrt2SigmaInv, double* myF, double x, int N);
__global__ void logKernal(double* y, int N);
__global__ void vec_logspace_add_addKernal(double *a, double *b, double add, int N);
__global__ void vec_logspace_addKernal(double *a, double *b, int N);
*/
