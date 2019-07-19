/*
 * blas_tmpl.h
 *
 *  Created on: 31.01.2018
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef BLAS_TMPL_H_
#define BLAS_TMPL_H_

#define BLAS_NEEDBUNDERSCORE

#ifdef BLAS_NEEDBUNDERSCORE
#define BLASFUNC(FUNC) FUNC##_
#else
#define BLASFUNC(FUNC) FUNC
#endif

typedef int BLAS_INTEGER;

extern "C"
{
// LDLT solve
void BLASFUNC(dsysv)( const char* uplo, const BLAS_INTEGER* n, const BLAS_INTEGER* nrhs, double* a, const BLAS_INTEGER* lda,
                BLAS_INTEGER* ipiv, double* b, const BLAS_INTEGER* ldb, double* work, const BLAS_INTEGER* lwork, BLAS_INTEGER* info );
void BLASFUNC(ssysv)( const char* uplo, const BLAS_INTEGER* n, const BLAS_INTEGER* nrhs, float* a, const BLAS_INTEGER* lda,
                BLAS_INTEGER* ipiv, float* b, const BLAS_INTEGER* ldb, float* work, const BLAS_INTEGER* lwork, BLAS_INTEGER* info );
// General linear solve
void BLASFUNC(dgesv)(const BLAS_INTEGER *N, const BLAS_INTEGER *NRHS, double *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, double *B, const BLAS_INTEGER *LDB, BLAS_INTEGER *INFO);
void BLASFUNC(sgesv)(const BLAS_INTEGER *N, const BLAS_INTEGER *NRHS, float *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, float *B, const BLAS_INTEGER *LDB, BLAS_INTEGER *INFO);
double BLASFUNC(dasum)(const BLAS_INTEGER *N, const double *DX, const BLAS_INTEGER *INCX);
double BLASFUNC(sasum)(const BLAS_INTEGER *N, const float *DX, const BLAS_INTEGER *INCX);
// scales a vector
void BLASFUNC(dscal)(const BLAS_INTEGER *N, const double *DA, double *DX, const BLAS_INTEGER *INCX);
void BLASFUNC(sscal)(const BLAS_INTEGER *N, const float *DA, float *DX, const BLAS_INTEGER *INCX);
// copy vector
BLAS_INTEGER BLASFUNC(dcopy)(const BLAS_INTEGER *nn, const double *B, const BLAS_INTEGER *incb, double *C, const BLAS_INTEGER *incC);
BLAS_INTEGER BLASFUNC(scopy)(const BLAS_INTEGER *nn, const float *B, const BLAS_INTEGER *incb, float *C, const BLAS_INTEGER *incC);
// copy matrix
void BLASFUNC(dlacpy)(const char *uplo, const BLAS_INTEGER *mm, const BLAS_INTEGER *nn, const double *A, const BLAS_INTEGER *LDA, double *B, const BLAS_INTEGER *LDB);
void BLASFUNC(slacpy)(const char *uplo, const BLAS_INTEGER *mm, const BLAS_INTEGER *nn, const float *A, const BLAS_INTEGER *LDA, float *B, const BLAS_INTEGER *LDB);
// scalar vector multiplication
double BLASFUNC(ddot)(const BLAS_INTEGER *n, const double *B, const BLAS_INTEGER *incb, const double *C, const BLAS_INTEGER *incC);
double BLASFUNC(sdot)(const BLAS_INTEGER *n, const float *B, const BLAS_INTEGER *incb, const float *C, const BLAS_INTEGER *incC);

void BLASFUNC(dger) (const BLAS_INTEGER *M, const BLAS_INTEGER *N, const double *ALPHA, const double *X, const BLAS_INTEGER *INCX, const double *Y, const BLAS_INTEGER *INCY, double *A, const BLAS_INTEGER *LDA);
void BLASFUNC(sger) (const BLAS_INTEGER *M, const BLAS_INTEGER *N, const float *ALPHA, const float *X, const BLAS_INTEGER *INCX, const float *Y, const BLAS_INTEGER *INCY, float *A, const BLAS_INTEGER *LDA);
// Y = a*X+Y
void BLASFUNC(daxpy)(const BLAS_INTEGER *N, const double *SA, const double *SX, const BLAS_INTEGER *INCX, double *SY, const BLAS_INTEGER *INCY);
void BLASFUNC(saxpy)(const BLAS_INTEGER *N, const float *SA, const float *SX, const BLAS_INTEGER *INCX, float *SY, const BLAS_INTEGER *INCY);
// mul matrices
BLAS_INTEGER BLASFUNC(dgemm)(const char *transA, const char *transB, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const BLAS_INTEGER *common, const double *alpha, const double *a, const BLAS_INTEGER *lda,
              const double *b, const BLAS_INTEGER *ldb, const double *beta, double *c, const BLAS_INTEGER *ldc);
BLAS_INTEGER BLASFUNC(sgemm)(const char *transA, const char *transB, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const BLAS_INTEGER *common, const float *alpha, const float *a, const BLAS_INTEGER *lda,
              const float *b, const BLAS_INTEGER *ldb, const float *beta, float *c, const BLAS_INTEGER *ldc);
BLAS_INTEGER BLASFUNC(dgemv)(const char *TRANS, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const double *ALPHA, const double *A, const BLAS_INTEGER *LDA,
		const double *X, const BLAS_INTEGER *INCX, const double *BETA, double *Y, const BLAS_INTEGER *INCY);
BLAS_INTEGER BLASFUNC(sgemv)(const char *TRANS, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const float *ALPHA, const float *A, const BLAS_INTEGER *LDA,
		const float *X, const BLAS_INTEGER *INCX, const float *BETA, float *Y, const BLAS_INTEGER *INCY);
BLAS_INTEGER BLASFUNC(dsymm)(const char *SIDE, const char *UPLO, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const double *ALPHA, double *A, const BLAS_INTEGER *LDA, double *B, const BLAS_INTEGER *LDB,
		const double *BETA, double *C, const BLAS_INTEGER *LDC);
BLAS_INTEGER BLASFUNC(ssymm)(const char *SIDE, const char *UPLO, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const float *ALPHA, float *A, const BLAS_INTEGER *LDA, float *B, const BLAS_INTEGER *LDB,
		const float *BETA, float *C, const BLAS_INTEGER *LDC);
BLAS_INTEGER BLASFUNC(dsymv)(const char *UPLO, const BLAS_INTEGER *N, const double *ALPHA, double *A, const BLAS_INTEGER *LDA, double *X, const BLAS_INTEGER *INCX,
		const double *BETA, double *Y, const BLAS_INTEGER *INCY);
BLAS_INTEGER BLASFUNC(ssymv)(const char *UPLO, const BLAS_INTEGER *N, const float *ALPHA, float *A, const BLAS_INTEGER *LDA, float *X, const BLAS_INTEGER *INCX,
		const float *BETA, float *Y, const BLAS_INTEGER *INCY);
BLAS_INTEGER BLASFUNC(dsbmv)(const char *UPLO, const BLAS_INTEGER *N, const BLAS_INTEGER *K, const double *ALPHA, const double *A, const BLAS_INTEGER *LDA, const double *X, const BLAS_INTEGER *INCX,
		const double *BETA, double *Y, const BLAS_INTEGER *INCY);
BLAS_INTEGER BLASFUNC(ssbmv)(const char *UPLO, const BLAS_INTEGER *N, const BLAS_INTEGER *K, const float *ALPHA, const float *A, const BLAS_INTEGER *LDA, const float *X, const BLAS_INTEGER *INCX,
		const float *BETA, float *Y, const BLAS_INTEGER *INCY);
// Сholesky decomposition
BLAS_INTEGER BLASFUNC(dpotrf)(
	    const char *uplo, const BLAS_INTEGER *n, double *a, const BLAS_INTEGER *lda, BLAS_INTEGER *info);
BLAS_INTEGER BLASFUNC(spotrf)(
	    const char *uplo, const BLAS_INTEGER *n, float *a, const BLAS_INTEGER *lda, BLAS_INTEGER *info);
// Inverse based on Cholesky decomposition
BLAS_INTEGER BLASFUNC(dpotri)(const char *uplo, const BLAS_INTEGER *n, double *a, const BLAS_INTEGER *lda, BLAS_INTEGER *info);
BLAS_INTEGER BLASFUNC(spotri)(const char *uplo, const BLAS_INTEGER *n, float *a, const BLAS_INTEGER *lda, BLAS_INTEGER *info);
// Solve. Input - Cholesky
BLAS_INTEGER BLASFUNC(dpotrs)(const char *uplo, const BLAS_INTEGER *n, const BLAS_INTEGER *nrhs, double *a, const BLAS_INTEGER *lda,
	    double *b, const BLAS_INTEGER *ldb, BLAS_INTEGER *info);
BLAS_INTEGER BLASFUNC(spotrs)(const char *uplo, const BLAS_INTEGER *n, const BLAS_INTEGER *nrhs, float *a, const BLAS_INTEGER *lda,
	    float *b, const BLAS_INTEGER *ldb, BLAS_INTEGER *info);
// Gets eigenvalues and eigenvectors for symmentric matrix A
void BLASFUNC(dsyev)(const char *JOBZ, const char *UPLO, BLAS_INTEGER *N, double *A, BLAS_INTEGER *LDA, double *W, double *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO);
void BLASFUNC(ssyev)(const char *JOBZ, const char *UPLO, BLAS_INTEGER *N, float *A, BLAS_INTEGER *LDA, float *W, float *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO);
// Solve symm
BLAS_INTEGER BLASFUNC(dposv)(const char *uplo, const BLAS_INTEGER *n, const BLAS_INTEGER *nrhs, double *a, const BLAS_INTEGER *lda, double *b, const BLAS_INTEGER *ldb, BLAS_INTEGER *info);
BLAS_INTEGER BLASFUNC(sposv)(const char *uplo, const BLAS_INTEGER *n, const BLAS_INTEGER *nrhs, float *a, const BLAS_INTEGER *lda, float *b, const BLAS_INTEGER *ldb, BLAS_INTEGER *info);

void BLASFUNC(dgesvd)(const char *JOBU, const char *JOBVT, const BLAS_INTEGER *M, const BLAS_INTEGER *N,
		double* A, const BLAS_INTEGER *LDA, double *S, double *U, const BLAS_INTEGER *LDU, double *VT, const BLAS_INTEGER *LDVT,
		const double *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO);
void BLASFUNC(sgesvd)(const char *JOBU, const char *JOBVT, const BLAS_INTEGER *M, const BLAS_INTEGER *N,
		float* A, const BLAS_INTEGER *LDA, float *S, float *U, const BLAS_INTEGER *LDU, float *VT, const BLAS_INTEGER *LDVT,
		const float *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO);
double BLASFUNC(dlange)(const char *NORM, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const double *A, const BLAS_INTEGER *LDA, double *WORK);
double BLASFUNC(slange)(const char *NORM, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const float *A, const BLAS_INTEGER *LDA, float *WORK);

void BLASFUNC(dlasrt)(const char *ID, const BLAS_INTEGER *N, double *D, BLAS_INTEGER *INFO);
void BLASFUNC(slasrt)(const char *ID, const BLAS_INTEGER *N, float *D, BLAS_INTEGER *INFO);

void BLASFUNC(dgetrf)(const BLAS_INTEGER *M, const BLAS_INTEGER *N, double *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, BLAS_INTEGER *INFO);
void BLASFUNC(sgetrf)(const BLAS_INTEGER *M, const BLAS_INTEGER *N, float *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, BLAS_INTEGER *INFO);
void BLASFUNC(dgetrs)(const char *Trans, const BLAS_INTEGER *N, const BLAS_INTEGER *NRHS, const double *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, double *B, const BLAS_INTEGER *LDB, BLAS_INTEGER *INFO);
void BLASFUNC(sgetrs)(const char *Trans, const BLAS_INTEGER *N, const BLAS_INTEGER *NRHS, const float *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, float *B, const BLAS_INTEGER *LDB, BLAS_INTEGER *INFO);
void BLASFUNC(dsyev)(const char *JOBZ, const char *UPLO, BLAS_INTEGER *N, double *A, BLAS_INTEGER *LDA, double *W, double *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO);
void BLASFUNC(ssyev)(const char *JOBZ, const char *UPLO, BLAS_INTEGER *N, float *A, BLAS_INTEGER *LDA, float *W, float *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO);
// converting float-double
void BLASFUNC(slag2d)(const BLAS_INTEGER *M, const BLAS_INTEGER *N, const float *SA, const BLAS_INTEGER *LDSA, double *A, const BLAS_INTEGER *LDA, const BLAS_INTEGER *INFO);
void BLASFUNC(dlag2s)(const BLAS_INTEGER *M, const BLAS_INTEGER *N, const double *A, const BLAS_INTEGER *LDA, float *SA, const BLAS_INTEGER *LDSA, const BLAS_INTEGER *INFO);

}

namespace iTTL
{

template <typename T>
struct type_constants
{
	static const T m_one;
	static const T zero;
	static const T one;
};

template <typename T>
const T type_constants<T>::m_one=-1;
template <typename T>
const T type_constants<T>::zero=0;
template <typename T>
const T type_constants<T>::one=1;

template <typename T>
T dot(const BLAS_INTEGER *n, const T *B, const BLAS_INTEGER *incb, const T *C, const BLAS_INTEGER *incC);

template <>
inline double dot<double>(const BLAS_INTEGER *n, const double *B, const BLAS_INTEGER *incB, const double *C, const BLAS_INTEGER *incC)
{
	return BLASFUNC(ddot)(n,B,incB,C,incC);
}
template <>
inline float dot<float>(const BLAS_INTEGER *n, const float *B, const BLAS_INTEGER *incB, const float *C, const BLAS_INTEGER *incC)
{
	return BLASFUNC(sdot)(n,B,incB,C,incC);
}

template <typename T>
void ger(const BLAS_INTEGER *M, const BLAS_INTEGER *N, const T *ALPHA, const T *X, const BLAS_INTEGER *INCX, const T *Y, const BLAS_INTEGER *INCY, T *A, const BLAS_INTEGER *LDA);

template <>
inline void ger<double>(const BLAS_INTEGER *M, const BLAS_INTEGER *N, const double *ALPHA, const double *X, const BLAS_INTEGER *INCX, const double *Y, const BLAS_INTEGER *INCY, double *A, const BLAS_INTEGER *LDA)
{
	BLASFUNC(dger)(M, N, ALPHA, X, INCX, Y, INCY, A, LDA);
}
template <>
inline void ger<float>(const BLAS_INTEGER *M, const BLAS_INTEGER *N, const float *ALPHA, const float *X, const BLAS_INTEGER *INCX, const float *Y, const BLAS_INTEGER *INCY, float *A, const BLAS_INTEGER *LDA)
{
	BLASFUNC(sger)(M, N, ALPHA, X, INCX, Y, INCY, A, LDA);
}

template <typename T>
BLAS_INTEGER copy(const BLAS_INTEGER *nn, const T *B, const BLAS_INTEGER *incb, T *C, const BLAS_INTEGER *incC);

template <>
inline BLAS_INTEGER copy<double>(const BLAS_INTEGER *nn, const double *B, const BLAS_INTEGER *incb, double *C, const BLAS_INTEGER *incC)
{
	return BLASFUNC(dcopy)(nn, B, incb, C, incC);
}
template <>
inline BLAS_INTEGER copy<float>(const BLAS_INTEGER *nn, const float *B, const BLAS_INTEGER *incb, float *C, const BLAS_INTEGER *incC)
{
	return BLASFUNC(scopy)(nn, B, incb, C, incC);
}

template <typename T>
void lacpy(const char *uplo, const BLAS_INTEGER *mm, const BLAS_INTEGER *nn, const T *A, const BLAS_INTEGER *LDA, T *B, const BLAS_INTEGER *LDB);

template <>
inline void lacpy<double>(const char *uplo, const BLAS_INTEGER *mm, const BLAS_INTEGER *nn, const double *A, const BLAS_INTEGER *LDA, double *B, const BLAS_INTEGER *LDB)
{
	BLASFUNC(dlacpy)(uplo, mm, nn, A, LDA, B, LDB);
}
template <>
inline void lacpy<float>(const char *uplo, const BLAS_INTEGER *mm, const BLAS_INTEGER *nn, const float *A, const BLAS_INTEGER *LDA, float *B, const BLAS_INTEGER *LDB)
{
	BLASFUNC(slacpy)(uplo, mm, nn, A, LDA, B, LDB);
}


template <typename T>
void axpy(const BLAS_INTEGER *N, const T *SA, const T *SX, const BLAS_INTEGER *INCX, T *SY, const BLAS_INTEGER *INCY);

template <>
inline void axpy<double>(const BLAS_INTEGER *N, const double *SA, const double *SX, const BLAS_INTEGER *INCX, double *SY, const BLAS_INTEGER *INCY)
{
	BLASFUNC(daxpy)(N, SA, SX, INCX, SY, INCY);
}
template <>
inline void axpy<float>(const BLAS_INTEGER *N, const float *SA, const float *SX, const BLAS_INTEGER *INCX, float *SY, const BLAS_INTEGER *INCY)
{
	BLASFUNC(saxpy)(N, SA, SX, INCX, SY, INCY);
}

template <typename T>
BLAS_INTEGER gemm(const char *transA, const char *transB, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const BLAS_INTEGER *common, const T *alpha, const T *a, const BLAS_INTEGER *lda,
        const T *b, const BLAS_INTEGER *ldb, const T *beta, T *c, const BLAS_INTEGER *ldc);


template <>
inline BLAS_INTEGER gemm<double>(const char *transA, const char *transB, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const BLAS_INTEGER *common, const double *alpha, const double *a, const BLAS_INTEGER *lda,
        const double *b, const BLAS_INTEGER *ldb, const double *beta, double *c, const BLAS_INTEGER *ldc)
{
	return BLASFUNC(dgemm)(transA, transB, M, N, common, alpha, a, lda, b, ldb, beta, c, ldc);
}
template <>
inline BLAS_INTEGER gemm<float>(const char *transA, const char *transB, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const BLAS_INTEGER *common, const float *alpha, const float *a, const BLAS_INTEGER *lda,
        const float *b, const BLAS_INTEGER *ldb, const float *beta, float *c, const BLAS_INTEGER *ldc)
{
	return BLASFUNC(sgemm)(transA, transB, M, N, common, alpha, a, lda, b, ldb, beta, c, ldc);
}

template <typename T>
BLAS_INTEGER gemv(const char *TRANS, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const T *ALPHA, const T *A, const BLAS_INTEGER *LDA,
		const T *X, const BLAS_INTEGER *INCX, const T *BETA, T *Y, const BLAS_INTEGER *INCY);


template <>
inline BLAS_INTEGER gemv<double>(const char *TRANS, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const double *ALPHA, const double *A, const BLAS_INTEGER *LDA,
		const double *X, const BLAS_INTEGER *INCX, const double *BETA, double *Y, const BLAS_INTEGER *INCY)
{
	return BLASFUNC(dgemv)(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY);
}
template <>
inline BLAS_INTEGER gemv<float>(const char *TRANS, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const float *ALPHA, const float *A, const BLAS_INTEGER *LDA,
		const float *X, const BLAS_INTEGER *INCX, const float *BETA, float *Y, const BLAS_INTEGER *INCY)
{
	return BLASFUNC(sgemv)(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY);
}


template <typename T>
BLAS_INTEGER symm(const char *SIDE, const char *UPLO, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const T *ALPHA, T *A, const BLAS_INTEGER *LDA, T *B, const BLAS_INTEGER *LDB,
		const T *BETA, T *C, const BLAS_INTEGER *LDC);


template <>
inline BLAS_INTEGER symm<double>(const char *SIDE, const char *UPLO, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const double *ALPHA, double *A, const BLAS_INTEGER *LDA, double *B, const BLAS_INTEGER *LDB,
		const double *BETA, double *C, const BLAS_INTEGER *LDC)
{
	return BLASFUNC(dsymm)(SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}
template <>
inline BLAS_INTEGER symm<float>(const char *SIDE, const char *UPLO, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const float *ALPHA, float *A, const BLAS_INTEGER *LDA, float *B, const BLAS_INTEGER *LDB,
		const float *BETA, float *C, const BLAS_INTEGER *LDC)
{
	return BLASFUNC(ssymm)(SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}

template <typename T>
BLAS_INTEGER symv(const char *UPLO, const BLAS_INTEGER *N, const T *ALPHA, T *A, const BLAS_INTEGER *LDA, T *X, const BLAS_INTEGER *INCX,
		const T *BETA, T *Y, const BLAS_INTEGER *INCY);


template <>
inline BLAS_INTEGER symv<double>(const char *UPLO, const BLAS_INTEGER *N, const double *ALPHA, double *A, const BLAS_INTEGER *LDA, double *X, const BLAS_INTEGER *INCX,
		const double *BETA, double *Y, const BLAS_INTEGER *INCY)
{
	return BLASFUNC(dsymv)(UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY);
}
template <>
inline BLAS_INTEGER symv<float>(const char *UPLO, const BLAS_INTEGER *N, const float *ALPHA, float *A, const BLAS_INTEGER *LDA, float *X, const BLAS_INTEGER *INCX,
		const float *BETA, float *Y, const BLAS_INTEGER *INCY)
{
	return BLASFUNC(ssymv)(UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY);
}

template <typename T>
BLAS_INTEGER sbmv(const char *UPLO, const BLAS_INTEGER *N, const BLAS_INTEGER *K, const T *ALPHA, const T *A, const BLAS_INTEGER *LDA, const T *X, const BLAS_INTEGER *INCX,
		const T *BETA, T *Y, const BLAS_INTEGER *INCY);


template <>
inline BLAS_INTEGER sbmv<double>(const char *UPLO, const BLAS_INTEGER *N, const BLAS_INTEGER *K, const double *ALPHA, const double *A, const BLAS_INTEGER *LDA, const double *X, const BLAS_INTEGER *INCX,
		const double *BETA, double *Y, const BLAS_INTEGER *INCY)
{
	return BLASFUNC(dsbmv)(UPLO, N, K, ALPHA, A, LDA, X, INCX, BETA, Y, INCY);
}
template <>
inline BLAS_INTEGER sbmv<float>(const char *UPLO, const BLAS_INTEGER *N, const BLAS_INTEGER *K, const float *ALPHA, const float *A, const BLAS_INTEGER *LDA, const float *X, const BLAS_INTEGER *INCX,
		const float *BETA, float *Y, const BLAS_INTEGER *INCY)
{
	return BLASFUNC(ssbmv)(UPLO, N, K, ALPHA, A, LDA, X, INCX, BETA, Y, INCY);
}


template <typename T>
void scal(const BLAS_INTEGER *N, const T *DA, T *DX, const BLAS_INTEGER *INCX);

template<>
inline void scal<double>(const BLAS_INTEGER *N, const double *DA, double *DX, const BLAS_INTEGER *INCX)
{
	BLASFUNC(dscal)(N, DA, DX, INCX);
}
template<>
inline void scal<float>(const BLAS_INTEGER *N, const float *DA, float *DX, const BLAS_INTEGER *INCX)
{
	BLASFUNC(sscal)(N, DA, DX, INCX);
}

template <typename T>
void potrf(const char *uplo, const BLAS_INTEGER *n, T *a, const BLAS_INTEGER *lda, BLAS_INTEGER *info);

template<>
inline void potrf<double>(const char *uplo, const BLAS_INTEGER *n, double *a, const BLAS_INTEGER *lda, BLAS_INTEGER *info)
{
	BLASFUNC(dpotrf)(uplo, n, a, lda, info);
}
template<>
inline void potrf<float>(const char *uplo, const BLAS_INTEGER *n, float *a, const BLAS_INTEGER *lda, BLAS_INTEGER *info)
{
	BLASFUNC(spotrf)(uplo, n, a, lda, info);
}

template <typename T>
void potri(const char *uplo, const BLAS_INTEGER *n, T *a, const BLAS_INTEGER *lda, BLAS_INTEGER *info);

template<>
inline void potri<double>(const char *uplo, const BLAS_INTEGER *n, double *a, const BLAS_INTEGER *lda, BLAS_INTEGER *info)
{
	BLASFUNC(dpotri)(uplo, n, a, lda, info);
}
template<>
inline void potri<float>(const char *uplo, const BLAS_INTEGER *n, float *a, const BLAS_INTEGER *lda, BLAS_INTEGER *info)
{
	BLASFUNC(spotri)(uplo, n, a, lda, info);
}

template <typename T>
void potrs(const char *uplo, const BLAS_INTEGER *n, const BLAS_INTEGER *nrhs, T *a, const BLAS_INTEGER *lda,
	    T *b, const BLAS_INTEGER *ldb, BLAS_INTEGER *info);

template<>
inline void potrs<double>(const char *uplo, const BLAS_INTEGER *n, const BLAS_INTEGER *nrhs, double *a, const BLAS_INTEGER *lda,
	    double *b, const BLAS_INTEGER *ldb, BLAS_INTEGER *info)
{
	BLASFUNC(dpotrs)(uplo, n, nrhs, a, lda, b, ldb, info);
}
template<>
inline void potrs<float>(const char *uplo, const BLAS_INTEGER *n, const BLAS_INTEGER *nrhs, float *a, const BLAS_INTEGER *lda,
	    float *b, const BLAS_INTEGER *ldb, BLAS_INTEGER *info)
{
	BLASFUNC(spotrs)(uplo, n, nrhs, a, lda, b, ldb, info);
}

template <typename T>
void syev(const char *JOBZ, const char *UPLO, BLAS_INTEGER *N, T *A, BLAS_INTEGER *LDA, T *W, T *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO);

template<>
inline void syev<double>(const char *JOBZ, const char *UPLO, BLAS_INTEGER *N, double *A, BLAS_INTEGER *LDA, double *W, double *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO)
{
	BLASFUNC(dsyev)(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO);
}
template<>
inline void syev<float>(const char *JOBZ, const char *UPLO, BLAS_INTEGER *N, float *A, BLAS_INTEGER *LDA, float *W, float *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO)
{
	BLASFUNC(ssyev)(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO);
}


template <typename T>
T asum(const BLAS_INTEGER *N, const T *DX, const BLAS_INTEGER *INCX);

template <>
inline double asum<double>(const BLAS_INTEGER *N, const double *DX, const BLAS_INTEGER *INCX)
{
	return BLASFUNC(dasum)(N, DX, INCX);
}
template <>
inline float asum<float>(const BLAS_INTEGER *N, const float *DX, const BLAS_INTEGER *INCX)
{
	return BLASFUNC(sasum)(N, DX, INCX);
}

template <typename T>
void gesvd(const char *JOBU, const char *JOBVT, const BLAS_INTEGER *M, const BLAS_INTEGER *N,
		T* A, const BLAS_INTEGER *LDA, T *S, T *U, const BLAS_INTEGER *LDU, T *VT, const BLAS_INTEGER *LDVT,
		const T *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO);

template <>
inline void gesvd<double>(const char *JOBU, const char *JOBVT, const BLAS_INTEGER *M, const BLAS_INTEGER *N,
		double* A, const BLAS_INTEGER *LDA, double *S, double *U, const BLAS_INTEGER *LDU, double *VT, const BLAS_INTEGER *LDVT,
		const double *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO)
{
	BLASFUNC(dgesvd)(JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO);
}
template <>
inline void gesvd<float>(const char *JOBU, const char *JOBVT, const BLAS_INTEGER *M, const BLAS_INTEGER *N,
		float* A, const BLAS_INTEGER *LDA, float *S, float *U, const BLAS_INTEGER *LDU, float *VT, const BLAS_INTEGER *LDVT,
		const float *WORK, BLAS_INTEGER *LWORK, BLAS_INTEGER *INFO)
{
	BLASFUNC(sgesvd)(JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO);
}

template <typename T>
T lange(const char *NORM, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const T *A, const BLAS_INTEGER *LDA, T *WORK);

template<>
inline double lange<double>(const char *NORM, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const double *A, const BLAS_INTEGER *LDA, double *WORK)
{
	return BLASFUNC(dlange)(NORM, M, N, A, LDA, WORK);
}
template<>
inline float lange<float>(const char *NORM, const BLAS_INTEGER *M, const BLAS_INTEGER *N, const float *A, const BLAS_INTEGER *LDA, float *WORK)
{
	return BLASFUNC(slange)(NORM, M, N, A, LDA, WORK);
}

template <typename T>
void gesv(const BLAS_INTEGER *N, const BLAS_INTEGER *NRHS, T *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, T *B, const BLAS_INTEGER *LDB, BLAS_INTEGER *INFO);

template<>
inline void gesv<double>(const BLAS_INTEGER *N, const BLAS_INTEGER *NRHS, double *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, double *B, const BLAS_INTEGER *LDB, BLAS_INTEGER *INFO)
{
	BLASFUNC(dgesv)(N, NRHS, A, LDA, IPIV, B, LDB, INFO);
}
template<>
inline void gesv<float>(const BLAS_INTEGER *N, const BLAS_INTEGER *NRHS, float *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, float *B, const BLAS_INTEGER *LDB, BLAS_INTEGER *INFO)
{
	BLASFUNC(sgesv)(N, NRHS, A, LDA, IPIV, B, LDB, INFO);
}

template <typename T>
void lasrt(const char *ID, const BLAS_INTEGER *N, T *D, BLAS_INTEGER *INFO);

template<>
inline void lasrt<double>(const char *ID, const BLAS_INTEGER *N, double *D, BLAS_INTEGER *INFO)
{
	BLASFUNC(dlasrt)(ID, N, D, INFO);
}
template<>
inline void lasrt<float>(const char *ID, const BLAS_INTEGER *N, float *D, BLAS_INTEGER *INFO)
{
	BLASFUNC(slasrt)(ID, N, D, INFO);
}

template <typename T>
void getrf(const BLAS_INTEGER *M, const BLAS_INTEGER *N, T *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, BLAS_INTEGER *INFO);

template <>
inline void getrf<double>(const BLAS_INTEGER *M, const BLAS_INTEGER *N, double *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, BLAS_INTEGER *INFO)
{
	BLASFUNC(dgetrf)(M, N, A, LDA, IPIV, INFO);
}
template <>
inline void getrf<float>(const BLAS_INTEGER *M, const BLAS_INTEGER *N, float *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, BLAS_INTEGER *INFO)
{
	BLASFUNC(sgetrf)(M, N, A, LDA, IPIV, INFO);
}

template <typename T>
void getrs(const char *Trans, const BLAS_INTEGER *N, const BLAS_INTEGER *NRHS, const T *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, T *B, const BLAS_INTEGER *LDB, BLAS_INTEGER *INFO);

template <>
inline void getrs<double>(const char *Trans, const BLAS_INTEGER *N, const BLAS_INTEGER *NRHS, const double *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, double *B, const BLAS_INTEGER *LDB, BLAS_INTEGER *INFO)
{
	BLASFUNC(dgetrs)(Trans, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
}
template <>
inline void getrs<float>(const char *Trans, const BLAS_INTEGER *N, const BLAS_INTEGER *NRHS, const float *A, const BLAS_INTEGER *LDA, BLAS_INTEGER *IPIV, float *B, const BLAS_INTEGER *LDB, BLAS_INTEGER *INFO)
{
	BLASFUNC(sgetrs)(Trans, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
}

template <typename T0, typename T1>
void lag2(const BLAS_INTEGER *M, const BLAS_INTEGER *N, const T0 *A0, const BLAS_INTEGER *LDA0, T1 *A1, const BLAS_INTEGER *LDA1, const BLAS_INTEGER *INFO);

template <>
inline void lag2<float, double>(const BLAS_INTEGER *M, const BLAS_INTEGER *N, const float *A0, const BLAS_INTEGER *LDA0, double *A1, const BLAS_INTEGER *LDA1, const BLAS_INTEGER *INFO)
{
	BLASFUNC(slag2d)(M, N, A0, LDA0, A1, LDA1, INFO);
}
template <>
inline void lag2<double, float>(const BLAS_INTEGER *M, const BLAS_INTEGER *N, const double *A0, const BLAS_INTEGER *LDA0, float *A1, const BLAS_INTEGER *LDA1, const BLAS_INTEGER *INFO)
{
	BLASFUNC(dlag2s)(M, N, A0, LDA0, A1, LDA1, INFO);
}


} // namespace iTTL



#endif /* BLAS_TMPL_H_ */
