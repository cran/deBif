/***
   NAME
     globals
   DESCRIPTION
     Header file with global definitions, such as error codes, return
     values and function prototypes.

    Copyright (C) 2018, Andre M. de Roos, University of Amsterdam

    This file is part of the deBif software package.

    This is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software. If not, see <http://www.gnu.org/licenses/>.

    Last modification: AMdR - Apr 08, 2022
***/
#ifndef GLOBALS
#define GLOBALS

#if defined(__APPLE__) && !defined(__clang__)
// Work-arounds for using GCC with Accelerate framework on Mac OS Yosemite
// First address a bug in <os/base.h>. It should protect against __has_extension() being undefined. Provide a phony definition of __has_extension()
#ifndef __has_extension
#define __has_extension(x)        0
#endif
// Second, GCC doesn't support blocks. The following prevents using vImage features (simplified interoperability with Core Graphics and Core Video) that can not be used by GCC anyway
#define vImage_Utilities_h
#define vImage_CVUtilities_h
#endif

#include "ctype.h"
#include "float.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdarg.h"
#include "stdint.h"
#include "string.h"
#include <sys/stat.h>

#if !defined(TRUE)
#if defined(true)
#define TRUE                      true
#else
#define TRUE                      1
#endif
#endif

/*
 *====================================================================================================================================
 *  Various other macro definitions
 *====================================================================================================================================
 */

#ifdef _MSC_VER
#include <float.h>
#define issane(a)                 ((_fpclass(a) == _FPCLASS_NN) || (_fpclass(a) == _FPCLASS_NZ) || (_fpclass(a) == _FPCLASS_PZ) || (_fpclass(a) == _FPCLASS_PN))
#else
#define issane(a)                 ((fpclassify(a) == FP_ZERO) || (fpclassify(a) == FP_NORMAL))
#endif

#define max(a, b)                 (((a) > (b)) ? (a) : (b))
#define min(a, b)                 (((a) < (b)) ? (a) : (b))
#define sign(a)                   (((a) < 0.0) ? (-1.0) : (1.0))

#if (defined(CURVE))
#undef EXTERN
#define EXTERN                    extern
#else
#undef EXTERN
#define EXTERN
#endif

#define R_USE_C99_IN_CXX
// See https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Fortran-character-strings
#define USE_FC_LEN_T
#include "R.h"
#include <Rdefines.h>
#include <Rinternals.h>

#include "R_ext/Print.h"
#include <R_ext/Utils.h>


/*
 *====================================================================================================================================
 *  Macros used in curve.c
 *====================================================================================================================================
 */

#define FAILURE                   0
#define SUCCES                    1

#define SINGULARITY               100
#define NORM_OVERFLOW             101
#define NO_CONVERGENCE            102
#define ILLEGAL_INPUT             103
#define FAILED_EVALUATION         104

#define FORWARD                   0                                                 // Computational methods for derivatives
#define CENTRAL                   1

#define UNDEFINED                 1000
#define BP                        1001
#define EQ                        1002
#define HP                        1003
#define LP                        1004
#define LC                        1005

#define MAX_STR_LEN               1024


/*
 *====================================================================================================================================
 *  Global variable definitions
 *====================================================================================================================================
 */

EXTERN double                     Jacobian_Step;


/*
 *====================================================================================================================================
 *  Function prototypeing: curve.c
 *====================================================================================================================================
 */

int                               ErrorMsg(const char *msg);
int                               checkInterrupt();
int                               FindPoint(const int pntdim, const int freeparsdim, double *guess, double *tanvec,
                                            double rhstol, double vartol, const int max_iter, int *niter,
                                            int (*fnc)(double *, double *),
                                            int (jacfun)(const int, double *, const int, double *, int (*fnc)(double *, double *), int));
int                               Jacobian(const int pntdim, double *pnt, const int fncdim, double *jac,
                                           int (*fnc)(double *, double *), int method);
int                               TangentVec(const int pntdim, double *sol, double *JacExport, double *tanvec,
                                             int (*fnc)(double *, double *),
                                             int (jacfun)(const int, double *, const int, double *, int (*fnc)(double *, double *), int),
                                             double *det);
int                               Determinant(const int N, double *M, double *det, double *cond);
int                               SolveLinearSystem(const int N, double *A, double *B);


/*
 *====================================================================================================================================
 *  Single/double precision and BLAS/LAPACK function mapping
 *====================================================================================================================================
 */

#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
// See https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Fortran-character-strings
#ifndef FCONE
#define FCONE
#endif

// Lapack index type

#define LAPACK_SIZE_T             int

/*
 * Routines and operations for which Lapack functions are used:
 *
 * Solving a system of Linear equations:
 *        FindPoint()   - Compute adjustment in Newton iteration
 *        TangentVec()  - Compute tangent vector of solution curve using the full Jacobian
 *
 * Determine matrix determinant and/or condition:
 *        TangentVec()  - Determine matrix determinant for BP detection
 */

#define dgetrf                    F77_CALL(dgetrf)        // Used in: Determinant()
#define dgecon                    F77_CALL(dgecon)        // Used in: Determinant()
#define dgesvx                    F77_CALL(dgesvx)        // Used in: SolveLinearSystem()


/*==================================================================================================================================*/

#endif /* GLOBALS */
