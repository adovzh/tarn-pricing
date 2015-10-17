/** \file  InterfaceCLAPACK.cpp
    \brief C++ source file implementing interface routines between Blitz++ classes and CLAPACK.
           Copyright 2005, 2011 by Erik Schlögl

		   Redistribution and use in source and binary forms, with or without modification, are permitted provided
		   that the following conditions are met:
		
           -# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

           -# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

           -# The name of the copyright holder may not be used to endorse or promote products derived from this software without specific prior written permission.
		
           THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
		   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
		   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
		   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
		   OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
		   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
		   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
		   POSSIBILITY OF SUCH DAMAGE.
           */

#include <stdexcept>
#include "InterfaceCLAPACK.hpp"

extern "C" {
  /* DSPEV computes all the eigenvalues and, optionally, eigenvectors of a 
     real symmetric Array A in packed storage. Linked in from CLAPACK.lib */
  extern int dspev_(char *jobz, char *uplo, long int *n, double *
     ap, double *w, double *z, long int *ldz, double *work, 
     long int *info);

  /* DPOTRF computes the Cholesky factorization of a real symmetric   
    positive definite matrix A.   

    The factorization has the form   
       A = U**T * U,  if UPLO = 'U', or   
       A = L  * L**T,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular. 
	Linked in from CLAPACK.lib */
  extern int dpotrf_(char *uplo, long int *n, double *a, long int *
	lda, long int *info);

  /* DPOTRI computes the inverse of a real symmetric positive definite   
    matrix A using the Cholesky factorization A = U**T*U or A = L*L**T   
    computed by DPOTRF. Linked in from CLAPACK.lib */
  extern int dpotri_(char *uplo, long int *n, double *a, long int *
	lda, long int *info);
     
  /* DGESV computes the solution to a real system of linear equations   
       A * X = B,   
     where A is an N-by-N Array and X and B are N-by-NRHS matrices.
     Linked in from CLAPACK.lib */
  extern int dgesv_(long int *n, long int *nrhs, double *a, long int 
	*lda, long int *ipiv, double *b, long int *ldb, long int *info);

  /* DGESVD computes the singular value decomposition (SVD) of a real   
    M-by-N matrix A, optionally computing the left and/or right singular   
    vectors. The SVD is written   

         A = U * SIGMA * transpose(V)   

    where SIGMA is an M-by-N matrix which is zero except for its   
    min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and   
    V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA   
    are the singular values of A; they are real and non-negative, and   
    are returned in descending order.  The first min(m,n) columns of   
    U and V are the left and right singular vectors of A.
	Linked in from CLAPACK.lib  */
  extern int dgesvd_(char *jobu, char *jobvt, long int *m, long int *n, 
	double *a, long int *lda, double *s, double *u, long int *ldu, double *vt, long int *ldvt, double *work, long int *lwork, 
	long int *info);

  /* DGTTRS solves one of the systems of equations   
       A*X = B  or  A'*X = B,   
    with a tridiagonal matrix A using the LU factorization computed   
    by DGTTRF. */
  extern int dgttrs_(char *trans, long int *n, long int *nrhs, 
	double *dl, double *d__, double *du, double *du2, 
	long int *ipiv, double *b, long int *ldb, long int *info);

  /* DGTTRF computes an LU factorization of a real tridiagonal matrix A   
    using elimination with partial pivoting and row interchanges.   

    The factorization has the form   
       A = L * U   
    where L is a product of permutation and unit lower bidiagonal   
    matrices and U is upper triangular with nonzeros in only the main   
    diagonal and first two superdiagonals. */
  extern int dgttrf_(long int *n, double *dl, double *d__, 
	double *du, double *du2, long int *ipiv, long int *info);
}

inline double my_abs(double x) { return (x>0.0) ? x : -x; }

/// Determinant of a real symmetric positive definite matrix.
double quantfin::interfaceCLAPACK::PositiveSymmetricMatrixDeterminant(const Array<double,2>& A)
{
  int i;
  Array<double,2> L(A.extent(blitz::firstDim),A.extent(blitz::secondDim));
  quantfin::interfaceCLAPACK::Cholesky(A,L,'L');
  double det = 1.0;
  for (i=0;i<A.extent(blitz::firstDim);i++) det *= L(i,i);
  return det*det;
}

/// Calculate Cholesky factorization of a real symmetric positive definite matrix
void quantfin::interfaceCLAPACK::Cholesky(const Array<double,2>& A,Array<double,2>& triangular,char LorU)
{
  int i,j;
  long int n = A.rows();
  if (n!=A.columns()) throw(std::logic_error("Array must be square"));
  double* ap  = new double[n*n];
  double* pos = ap;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) *pos++ = A(j,i); }
  long int info = 0;
  dpotrf_(&LorU,&n,ap,&n,&info);
  triangular = 0.0;
  pos = ap;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      if ((LorU=='L')&&(j>=i)) triangular(j,i) = *pos;
      if ((LorU=='U')&&(j<=i)) triangular(j,i) = *pos;
	  pos++; }}
  delete[] ap;
  if (info) throw(std::logic_error("Cholesky factorization failed"));
}

/// Determinant of a real symmetric positive definite matrix.
void quantfin::interfaceCLAPACK::PositiveSymmetricMatrixInverse(const Array<double,2>& A,Array<double,2>& inverseA)
{
  int i,j;
  long int n = A.rows();
  if (n!=A.columns()) throw(std::logic_error("Array must be square"));
  double* ap  = new double[n*n];
  double* pos = ap;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) *pos++ = A(j,i); }
  long int info = 0;
  char LorU = 'L';
  dpotrf_(&LorU,&n,ap,&n,&info);
  if (info) {
    delete[] ap;
	throw(std::logic_error("Cholesky factorization failed")); }
  dpotri_(&LorU,&n,ap,&n,&info);
  pos = ap;
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  if (j>=i) inverseA(j,i) = inverseA(i,j) = *pos; 
	  pos++; }}
  delete[] ap;
  if (info) throw(std::logic_error("Matrix inversion failed"));
}

 void quantfin::interfaceCLAPACK::MoorePenroseInverse(const Array<double,2>& A,Array<double,2>& inverseA,double eps)
{
  int m = A.rows();
  int n = A.columns();
  Array<double,2> U(m,n),V(n,n);
  Array<double,1> sigma(n);
  SingularValueDecomposition(A,U,sigma,V);
  int rank = 0;
  while ((rank<sigma.extent(blitz::firstDim))&&(std::abs(sigma(rank))>eps)) {
	sigma(rank) = 1.0/sigma(rank);
	rank++; }
  while (rank<sigma.extent(blitz::firstDim)) {
	sigma(rank) = 0.0;
	rank++; }
  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::thirdIndex k;
  V = V(i,j) * sigma(j);
  inverseA = blitz::sum(V(i,k)*U(j,k),k);
}

void quantfin::interfaceCLAPACK::SingularValueDecomposition(const Array<double,2>& A,Array<double,2>& U,Array<double,1>& sigma,Array<double,2>& V)
{
  int i,j;
  long int m = A.rows();
  long int n = A.columns();
  long int lwork = 5 * std::max(m,n);
  double* ap  = new double[n*m];
  double* s   = new double[std::min(n,m)];
  double* u   = new double[m*m];
  double* vt  = new double[n*n];
  double* w   = new double[lwork];
  double* pos = ap;
  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) *pos++ = A(j,i); }
  long int info = 0;
  char jobu  = 'S';
  char jobvt = 'A';
  dgesvd_(&jobu,&jobvt,&m,&n,ap,&m,s,u,&m,vt,&n,w,&lwork,&info);
  sigma = 0.0;
  for (i=0;i<std::min(n,m);i++) sigma(i) = s[i];
  pos = u;
  for (i=0;i<n;i++) {
	for (j=0;j<m;j++) U(j,i) = *pos++; }
  pos = vt;
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) V(i,j) = *pos++; }
  delete[] ap;
  delete[] s;
  delete[] u;
  delete[] vt;
  delete[] w;
  if (info) throw(std::logic_error("Singular value decomposition failed"));
}

/// Solve system of linear equations A X = B using CLAPACK routines.
void quantfin::interfaceCLAPACK::SolveLinear(const Array<double,2>& A,Array<double,2>& X,const Array<double,2>& B)
{
  int i,j;
  long int n = A.rows();
  if (n!=A.columns()) throw(std::logic_error("Array must be square"));
  long int nrhs = B.columns();
  double* ap  = new double[n*n];
  double* pos = ap;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) *pos++ = A(j,i); }
  double* bp  = new double[nrhs*n];
  pos = bp;
  for (i=0;i<nrhs;i++) {
    for (j=0;j<n;j++) *pos++ = B(j,i); }
  long int* ipiv = new long int[n];
  long int info = 0;
  dgesv_(&n,&nrhs,ap,&n,ipiv,bp,&n,&info);
  if (!info) {
    pos = bp;
    for (i=0;i<nrhs;i++) {
      for (j=0;j<n;j++) X(j,i) = *pos++; }}
  delete[] ap;
  delete[] bp;
  delete[] ipiv;
  if (info) throw(std::logic_error("Linear equation solve failed"));
}

/// Solve symmetric eigenvalue problem using CLAPACK routines.
void quantfin::interfaceCLAPACK::SymmetricEigenvalueProblem(
                                   const Array<double,2>& A,      ///< Symmetric Array to be decomposed.
                                   Array<double,1>& eigval,       ///< Array (vector) containing all non-zero eigenvalues.
                                   Array<double,2>& eigvec,       ///< Array of eigenvectors (Array, each column is an eigenvector)
                                   double eps                     ///< Threshold for comparison to zero, default 1e-12
                                   )
{
  int i,j;
  long int n = A.rows();
  if (n!=A.columns()) throw(std::logic_error("Array must be square"));
  double* ap  = new double[(n*(n+1))/2];
  double* pos = ap;
  for (i=0;i<n;i++) {
    for (j=0;j<=i;j++) *pos++ = A(j,i); }
  double* w = new double[n];
  double* z = new double[n*n];
  double* work = new double[3*n];
  long int ldz  = n;
  long int info = 0;
  char jobz = 'V';
  char uplo = 'U';
  dspev_(&jobz,&uplo,&n,ap,w,z,&ldz,work,&info);
  if (!info) {
    int k = n;
    for (i=0;i<n;i++) {
      if (my_abs(w[i])<=eps) k--; }
    Array<double,1> val(k);
    Array<double,2> vec(n,k);
    int l = 0;
    pos = z;
    for (i=0;i<n;i++) {
      if (my_abs(w[i])>eps) {
        val(l) = w[i];
        for (j=0;j<n;j++) vec(j,l) = *pos++;
        l++; }
      else pos += n; }
    eigval.resize(k);
    eigvec.resize(n,k);
    eigval = val;
    eigvec = vec; }
  delete[] ap;
  delete[] w;
  delete[] z;
  delete[] work;
  if (info) throw(std::logic_error("Eigenvalue decomposition failed"));
}

/// Solve system of linear equations A X = B using CLAPACK routines, where A is a triagonal matrix.
void quantfin::interfaceCLAPACK::SolveTridiagonal(const Array<double,2>& A,Array<double,2>& X,const Array<double,2>& B)
{
  int i,j;
  long int n = A.rows();
  if (n!=A.columns()) throw(std::logic_error("Array must be square"));
  long int nrhs = B.columns();
  // subdiagonal
  double* dl  = new double[n-1];
  double* pos = dl;
  for (i=0;i<n-1;i++) *pos++ = A(i+1,i); 
  // diagonal
  double* d   = new double[n];
  pos = d;
  for (i=0;i<n;i++) *pos++ = A(i,i); 
  // superdiagonal
  double* du  = new double[n-1];
  pos = du;
  for (i=0;i<n-1;i++) *pos++ = A(i,i+1); 
  double* bp  = new double[nrhs*n];
  // right hand side
  pos = bp;
  for (i=0;i<nrhs;i++) {
    for (j=0;j<n;j++) *pos++ = B(j,i); }
  long int* ipiv = new long int[n];
  long int info = 0;
  // LU factorization
  double* du2 = new double[n-1];
  dgttrf_(&n,dl,d,du,du2,ipiv,&info);
  // Solve tridiagonal system
  static char trans = 'N';
  dgttrs_(&trans,&n,&nrhs,dl,d,du,du2,ipiv,bp,&n,&info);
  // prepare return value
  if (!info) {
    pos = bp;
    for (i=0;i<nrhs;i++) {
      for (j=0;j<n;j++) X(j,i) = *pos++; }}
  // clean up
  delete[] dl;
  delete[] d;
  delete[] du;
  delete[] du2;
  delete[] ipiv;
  if (info) throw(std::logic_error("Linear equation solve failed"));
}

/// Solve system of linear equations A X = B using CLAPACK routines, where A is a triagonal matrix.
void quantfin::interfaceCLAPACK::SolveTridiagonalSparse(const Array<double,2>& A,Array<double,2>& X,const Array<double,2>& B)
{
  int i,j;
  long int n = A.rows();
  if (3!=A.columns()) throw(std::logic_error("Tridiagonal matrix must be represented as N x 3 matrix"));
  long int nrhs = B.columns();
  // subdiagonal
  double* dl  = new double[n-1];
  double* pos = dl;
  for (i=0;i<n-1;i++) *pos++ = A(i+1,0); 
  // diagonal
  double* d   = new double[n];
  pos = d;
  for (i=0;i<n;i++) *pos++ = A(i,1); 
  // superdiagonal
  double* du  = new double[n-1];
  pos = du;
  for (i=0;i<n-1;i++) *pos++ = A(i,2); 
  double* bp  = new double[nrhs*n];
  // right hand side
  pos = bp;
  for (i=0;i<nrhs;i++) {
    for (j=0;j<n;j++) *pos++ = B(j,i); }
  long int* ipiv = new long int[n];
  long int info = 0;
  // LU factorization
  double* du2 = new double[n-1];
  dgttrf_(&n,dl,d,du,du2,ipiv,&info);
  // Solve tridiagonal system
  static char trans = 'N';
  dgttrs_(&trans,&n,&nrhs,dl,d,du,du2,ipiv,bp,&n,&info);
  // prepare return value
  if (!info) {
    pos = bp;
    for (i=0;i<nrhs;i++) {
      for (j=0;j<n;j++) X(j,i) = *pos++; }}
  // clean up
  delete[] dl;
  delete[] d;
  delete[] du;
  delete[] du2;
  delete[] ipiv;
  if (info) throw(std::logic_error("Linear equation solve failed"));
}

