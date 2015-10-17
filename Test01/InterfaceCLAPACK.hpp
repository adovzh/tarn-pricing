/** \file  InterfaceCLAPACK.hpp
    \brief Header file declaring interface routines between Blitz++ classes and CLAPACK.
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

#ifndef INTERFACECLAPACK_HPP
#define INTERFACECLAPACK_HPP
#include "MSWarnings.hpp"
#include <blitz/array.h>

namespace quantfin { namespace interfaceCLAPACK {
          
using blitz::Array;

/// Calculate Cholesky factorization of a real symmetric positive definite matrix.
void Cholesky(const Array<double,2>& A,Array<double,2>& triangular,char LorU);

void SingularValueDecomposition(const Array<double,2>& A,Array<double,2>& U,Array<double,1>& sigma,Array<double,2>& V);

/// Solve symmetric eigenvalue problem.
void SymmetricEigenvalueProblem(const Array<double,2>& A,Array<double,1>& eigval,Array<double,2>& eigvec,double eps = 1e-12);

/// Solve system of linear equations A X = B using CLAPACK routines.
void SolveLinear(const Array<double,2>& A,Array<double,2>& X,const Array<double,2>& B);

/// Solve system of linear equations A X = B using CLAPACK routines, where A is a tridiagonal matrix.
void SolveTridiagonal(const Array<double,2>& A,Array<double,2>& X,const Array<double,2>& B);
void SolveTridiagonalSparse(const Array<double,2>& A,Array<double,2>& X,const Array<double,2>& B);

/// Determinant of a real symmetric positive definite matrix.
double PositiveSymmetricMatrixDeterminant(const Array<double,2>& A);

/// Inverse of a real symmetric positive definite matrix.
void PositiveSymmetricMatrixInverse(const Array<double,2>& A,Array<double,2>& inverseA);

void MoorePenroseInverse(const Array<double,2>& A,Array<double,2>& inverseA,double eps = 1e-6);
}}

#endif
