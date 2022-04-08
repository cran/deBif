/***
  NAME
    curve
  DESCRIPTION
    Routines that are used for locating points on steady state bifurcation curves
    of non-linear ODE systems and curves determining dynamic regimes of such systems
    in two-parameter domains. More generally, the routines allow for computing fixed
    points of a system of non-linear equations.

    Copyright (C) 2021, Andre M. de Roos, University of Amsterdam

    This file is part of the deBif software package.

    deBif is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    deBif is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with deBif. If not, see <http://www.gnu.org/licenses/>.

    Last modification: AMdR - Mar 08, 2021
***/
#ifndef CURVE
#define CURVE
#endif
#include "globals.h"


/*
 *====================================================================================================================================
 *      Some numerical settings
 *====================================================================================================================================
 */

#ifndef RHSMAX
#define RHSMAX                    0.999
#endif

#ifndef FASTNUMERICS
#define FASTNUMERICS              1
#endif

#ifndef JACOBIAN_UPDATES
#define JACOBIAN_UPDATES          4
#endif

#ifndef JACOBIAN_MIN_STEP
#define JACOBIAN_MIN_STEP         1.0E-5
#endif

static int                        FastNumerics = FASTNUMERICS;
static char                       errstr[MAX_STR_LEN];

/*==================================================================================================================================*/

static double anorm(int rows, int cols, double *a)
{
  double       tmp, maxval = 0.0;

  for (int ii = 0; ii < rows; ii++)
    {
      tmp = 0.0;
      for (int jj = 0; jj < cols; jj++) tmp += *(a + ii * cols + jj);
      maxval = max(maxval, tmp);
    }

  return maxval;
}


static inline double dot(int dim, double *x, double *y)
{
  double                          ss = 0;

  for (int ii = 0; ii < dim; ii++) ss += x[ii] * y[ii];

   return ss;
}


static inline double nrm2(int dim, double *x)
{
  double                          ss = 0;

  for (int ii = 0; ii < dim; ii++) ss += x[ii] * x[ii];

   return sqrt(ss);
}


/*==================================================================================================================================*/

int ErrorMsg(const char *msg)

{
#if (defined(CMDLINEDEBUG))
  fprintf(stderr, "\n%s\n", msg);
#else

  REprintf("%s\n", msg);
  warning(msg);
  R_FlushConsole();
  R_ProcessEvents();

#endif

  return FAILURE;
}


/*==================================================================================================================================*/

int FindPoint(const int pntdim, const int freeparsdim, double *guess, double *tanvec,
              double rhstol, double vartol, const int max_iter, int *niter,
              int (*fnc)(double *, double *),
              int (jacfun)(const int, double *, const int, double *, int (*fnc)(double *, double *), int))

/*
 * FindPoint -  Routine locates a point on a curve determined by a
 *              system of non-linear, algebraic equations.
 *              The iteration adjusts the vector-elements following a simple
 *              Newton-Chord method with Broyden update (see Kuznetsov pg. 418).
 *              Pseudo-arclength continuation is used to continue past curve folds.
 *
 * Arguments -  pntdim      : The dimension of the solution point on the curve.
 *                            The dimension of the system of equations is
 *                            assumed to be exactly 1 less.
 *              freeparsdim : Number of free parameters
 *              guess       : Pointer to an array containing the initial point
 *                            to start the iteration from. The first element of
 *                            the vector is assumed to be non-adjustable parameter.
 *              tanvec      : Tangent vector along the solution curve.
 *              rhstol      : Tolerance for right-hand side of defining equations
 *              vartol      : Tolerance for variable change
 *              max_iter    : Maximum stepnumber allowed in iteration.
 *              fnc         : Pointer to function specifying the system of
 *                            equations. The function must have a (double)
 *                            pointer as first argument, containing the point
 *                            in which to evaluate the system and a (double)
 *                            pointer as second argument, containing the
 *                            results after evaluation of the equations.
 */

{
  register int  iter, i, j;
  int           rhsdim;
  int           pntdim2 = pntdim*pntdim;
  int           retcode = NO_CONVERGENCE;
  int           dyconverged = TRUE, rhsconverged = TRUE;
  double        ynorm, dynorm, rhsnorm;
  double        *dBaseMem, *y, *tv, *dy, *rhs;
  double        *Jac, *JacCopy;

  y = dBaseMem = calloc((4 * pntdim + 2 * pntdim2), sizeof(double));
  if (!dBaseMem) return ErrorMsg("Memory allocation error in FindPoint()");

  tv      = y + pntdim;
  dy      = tv + pntdim;
  rhs     = dy + pntdim;
  Jac     = rhs + pntdim;
  JacCopy = Jac + pntdim2;

  // If tangent vector is given, we are doing pseudo-arc length continuation
  rhsdim = pntdim - (tanvec != NULL);

  memcpy(y, guess, pntdim * sizeof(double));
  memset((void *)dy, 0, pntdim*sizeof(double));

  *niter = max_iter;
  // The iteration loop
  for (iter = 0; iter < max_iter; iter++)
    {
      // Compute norm of Y and of dY
      ynorm  = nrm2(pntdim, y);
      dynorm = nrm2(pntdim, dy);
      if (!issane(ynorm) || !issane(dynorm))
        {
          ErrorMsg("Point norm overflow in FindPoint");
          retcode = NORM_OVERFLOW;
          break;
        }

      // Compute RHS and its norm
      memset((void *)rhs, 0, pntdim*sizeof(double));
      if ((*fnc)(y, rhs) == FAILURE)
        {
          ErrorMsg("Right-hand side computation failed");
          retcode = FAILED_EVALUATION;
          break;
        }
      rhsnorm = nrm2(rhsdim, rhs);
      // Return if diverged
      if ((!issane(rhsnorm)) || (rhsnorm/(1.0 + rhsnorm) > RHSMAX))
        {
          ErrorMsg("RHS norm overflow in FindPoint");
          retcode = NORM_OVERFLOW;
          break;
        }

      dyconverged = ((dynorm/pntdim) < vartol);
      rhsconverged = ((rhsnorm/pntdim) < rhstol);
      // The dimension of rhs is pntdim, for example in case of BP localisation or PGR calculations
      if (dyconverged && rhsconverged)
        {
          memcpy(guess, y, pntdim * sizeof(double));

          retcode = SUCCES;
          *niter = iter;
          break;
        }

      // Compute Jacobian every JACOBIAN_UPDATES steps, otherwise the Jacobian is updated
      // via a Broyden update (see below)
      if (!(iter % JACOBIAN_UPDATES))
        {
          memset((void *)Jac, 0, (pntdim*pntdim)*sizeof(double));
          /*
           * Notice that the Jacobian is stored as
           *
           *                    |dF1/dy1 ... dFn/dy1|
           *                    |dF1/dy2 ... dFn/dy2|
           *               J =  |   .          .   |
           *                    |   .          .   |
           *                    |dF1/dyn ... dFn/dyn|
           *
           * The matrix is hence stored in column-wise (fortran) style.
           * From a C perspective this means that all coefficients pertaining to yi are to be found
           * in ROW i (as opposed to column i).
           *
           * Solving J.dy = -F(y) with dy = (dy1 ... dyn) and F(y) = (F1(y) ... Fn(y)) requires the variable
           * trans[1] to be defined as {"N"} (see programs/various/testlapack.c for details).
          */
          (*jacfun)(pntdim, y, rhsdim, Jac, fnc, FORWARD);                          // Compute J = F_x(X^k)
        }
      else                                                                          // Broyden update of Jacobian
        {
          dynorm = dot(pntdim, dy, dy);
          // See 10.7 on pg. 419 in Kuznetsov. Notice though that Jac is the transposed jacobian
          for (i = 0; i < pntdim; i++)
            for (j = 0; j < rhsdim; j++) Jac[j + i*rhsdim] += rhs[j]*dy[i]/dynorm;
        }

      // Extend the Jacobian matrix to include an additional row for the tangent
      // vector if we are doing pseudo-archlength continuation
      // If tangent is not given rhsdim == pntdim and we follow simple Newton
      memset((void *)JacCopy, 0, (pntdim*pntdim)*sizeof(double));
      for (i = 0; i < pntdim; i++)
        memcpy(JacCopy + i*pntdim, Jac + i*rhsdim, rhsdim*sizeof(double));          // Extract dF/dx
      for (int jj = 0; jj < rhsdim; jj++) dy[jj] = -rhs[jj];

      if (tanvec != NULL)
        {
          // When tangent is present, find new point via pseudo-arclength continuation
          for (i = 0; i < pntdim; i++)
            {
              *(JacCopy + i*pntdim + rhsdim) = tanvec[i];
              tv[i] = y[i] - guess[i];
            }
          dy[rhsdim] = dot(pntdim, tv, tanvec);
        }

      // Solve the linear system
      retcode = SolveLinearSystem(pntdim, JacCopy, dy);
      if (retcode != SUCCES)
        {
          sprintf(errstr, "Failed to solve linear system in FindPoint() on iteration %d", iter);
          ErrorMsg(errstr);
          break;
        }

      // Adjust point
      for (int jj = 0; jj < pntdim; jj++) y[jj] += dy[jj];
      retcode = NO_CONVERGENCE;
    }

  free(dBaseMem);
  dBaseMem = NULL;

  return retcode;
}


/*==================================================================================================================================*/

int TangentVec(const int pntdim, double *sol, double *JacExport, double *tanvec, int (*fnc)(double *, double *),
               int (jacfun)(const int, double *, const int, double *, int (*fnc)(double *, double *), int),
               double *det)

/*
 * TangentVec - routine determines the direction of the curve defined by the
 *              system of equations
 *
 *                    F(y) = 0
 *
 *              The point y is considered to have a dimension of exactly 1
 *              larger than the number of equations (i.e. the dimension of
 *              F(y)).
 *
 * Arguments -  pntdim  : The dimension of the solution point on the curve.
 *              y       : Pointer to an array containing the fixed point
 *              tanvec  : Pointer to return tangent vector
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 */

{
  register int  j;
  int           rhsdim = pntdim - 1, pntdim2 = pntdim*pntdim, retcode;
  double        norm;
  double        *dBaseMem, *y, *Jac, *JacCopy;

  y = dBaseMem = calloc((pntdim + 2*pntdim2), sizeof(double));
  if (!dBaseMem) return ErrorMsg("Memory allocation error in TangentVec()");

  Jac     = y + pntdim;
  JacCopy = Jac + pntdim2;

  // Initialize
  memcpy(y, sol, pntdim * sizeof(double));
  norm = nrm2(pntdim, y);
  if (!issane(norm))
    {
      ErrorMsg("Norm overflow in curvedir");
      free(dBaseMem);
      dBaseMem = NULL;
      return NORM_OVERFLOW;
    }

  // Determine the Jacobian of the extended system (variable plus parameter
  // dependence).
  (*jacfun)(pntdim, y, rhsdim, JacCopy, fnc, CENTRAL);
  if (JacExport) memcpy(JacExport, JacCopy, pntdim*rhsdim*sizeof(double));

  // Append the current tangent vector as the last row to the jacobian to
  // preserve direction. See the matcont manual at
  // http://www.matcont.ugent.be/manual.pdf, page 10 & 11
  // Notice, however, it is here added as the last COLUMN because of the
  // Fortran column-wise storage!

  for (j = 0; j < pntdim; j++)
    {
      memcpy(Jac + j*pntdim, JacCopy + j*rhsdim, rhsdim * sizeof(double));
      *(Jac + j*pntdim + rhsdim) = tanvec[j];
    }

  memset((void *)JacCopy, 0, (pntdim*pntdim)*sizeof(double));
  memcpy(JacCopy, Jac, pntdim2 * sizeof(double));
  memset((void *)tanvec, 0, pntdim*sizeof(double));
  tanvec[rhsdim] = 1.0;

  // Solve the linear system
  retcode = SolveLinearSystem(pntdim, JacCopy, tanvec);
  if (retcode != SUCCES)
    {
      ErrorMsg("Failed to solve for tangent vector in TangentVec()");
      memset((void *)tanvec, 0, pntdim*sizeof(double));
      tanvec[0] = 1.0;
      free(dBaseMem);
      dBaseMem = NULL;
      return retcode;
    }

  if (det)
    {
      // Replace the last row of the (saved) Jacobian with the newly computed tangent vector
      // to compute the determinant for BP detection
      for (j = 0; j < pntdim; j++)
        {
          memcpy(JacCopy + j*pntdim, Jac + j*pntdim, rhsdim * sizeof(double));
          *(JacCopy + j*pntdim + rhsdim) = tanvec[j];
        }
      Determinant(pntdim, JacCopy, det, NULL);
    }
  norm = nrm2(pntdim, tanvec);    // Normalize and store
  for (int ii = 0; ii < pntdim; ii++) tanvec[ii] /= norm;

  free(dBaseMem);
  dBaseMem = NULL;

  return SUCCES;
}


/*==================================================================================================================================*/

int CentralDerivative(int fncdim, int (*fnc)(double *, double *), double *farg, double *frhs, double *x, double h0, double *feq,
                      double *result, int fast)
{
  // Compute the derivative using the 5-point rule (x-h, x-h/2, x, x+h/2, x+h). Note that the central point is not used.
  // Compute the error using the difference between the 5-point and the 3-point rule (x-h,x,x+h). Again the central point is not used.
  int    i;
  double old, h, r[2], trunc[2], round[2], error[2];
  double fm1, fp1, fmh, fph, r3, r5, e3, e5, dy;

  old =*x;
  h   = h0;
  for (i = 0; i < 2; i++)
    {
      *x = old + h;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;

      // Trick to reduce precision errors if non-optimized computation is used. See Num. Recipes 9.7, pg. 388.
      if (fast)
        {
          memcpy(result, feq, fncdim * sizeof(double));
          h =*x - old;
        }
      else
        fp1 =*feq;

      *x = old - h;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;

      if (fast)
        {
          for (int jj = 0; jj < fncdim; jj++)
            {
              result[jj] -= feq[jj];
              result[jj] /= 2.0 * h;
            }
          *x = old;
          return SUCCES;
        }
      else
        fm1 =*feq;

      *x = old - h/2;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;
      fmh =*feq;

      *x = old + h/2;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;
      fph =*feq;

      r3 = 0.5*(fp1 - fm1);
      r5 = (4.0/3.0)*(fph - fmh) - (1.0/3.0)*r3;

      e3 = (fabs(fp1) + fabs(fm1))*DBL_EPSILON;
      e5 = 2.0*(fabs(fph) + fabs(fmh))*DBL_EPSILON + e3;

      // The next term is due to finite precision in x+h = O (eps * x)
      dy = max(fabs(r3/h), fabs(r5/h))*(fabs(*x)/h)*DBL_EPSILON;

      // The truncation error in the r5 approximation itself is O(h^4). However, for safety, we estimate the error from r5-r3,
      // which is O(h^2).  By scaling h we will minimise this estimated error, not the actual truncation error in r5.
      r[i]     = r5/h;
      trunc[i] = fabs((r5 - r3)/h);                                                 // Estimated truncation error O(h^2)
      round[i] = fabs(e5/h) + dy;                                                   // Rounding error (cancellations)
      error[i] = round[i] + trunc[i];

      // Compute an optimised stepsize to minimize the total error, using the scaling of the truncation error (O(h^2)) and
      // rounding error (O(1/h)).
      if ((i == 0) && ((round[0] < trunc[0]) && (round[0] > 0 && trunc[0] > 0)))
        h = (h0)*pow(round[0]/(2.0*trunc[0]), 1.0/3.0);
      else
        break;
    }

  // Check that the new error is smaller, and that the new derivative is consistent with the error bounds of the original estimate.
  if ((i == 1) && (error[1] < error[0]) && (fabs(r[1] - r[0]) < 4.0*error[0]))
    *result = r[1];
  else
    *result = r[0];

  *x = old;
  return SUCCES;
}


/*==================================================================================================================================*/

static int ForwardDerivative(int fncdim, int (*fnc)(double *, double *), double *farg, double *frhs, double *x, double h0,
                             double *feq, double *result, int fast)
{
  int    i;
  double old, h, r[2], trunc[2], round[2], error[2];
  double f1, f2, f3, f4, r2, r4, e4, dy;

  old =*x;
  h   = h0;
  for (i = 0; i < 2; i++)
    {
      *x = old + h;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;

      // Trick to reduce precision errors if non-optimized computation is used. See Num. Recipes 9.7, pg. 388.
      if (fast)
        {
          h =*x - old;
          memcpy(result, feq, fncdim * sizeof(double));
          *x = old;
          if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;

          for (int jj = 0; jj < fncdim; jj++)
            {
              result[jj] -= feq[jj];
              result[jj] /= h;
            }
          return SUCCES;
        }
      else
        f4 =*feq;

      *x = old + (3.0/4.0)*h;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;
      f3 =*feq;

      *x = old + h/2.0;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;
      f2 =*feq;

      *x = old + h/4.0;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;
      f1 =*feq;

      r2 = 2.0*(f4 - f2);
      r4 = (22.0/3.0)*(f4 - f3) - (62.0/3.0)*(f3 - f2) + (52.0/3.0)*(f2 - f1);

      // Estimate the rounding error for r4
      e4 = 2*20.67*(fabs(f4) + fabs(f3) + fabs(f2) + fabs(f1))*DBL_EPSILON;

      // The next term is due to finite precision in x+h = O (eps * x)
      dy = max(fabs(r2/h), fabs(r4/h))*fabs(*x/h)*DBL_EPSILON;

      // The truncation error in the r4 approximation itself is O(h^3). However, for safety, we estimate the error from r4-r2,
      // which is O(h).  By scaling h we will minimise this estimated error, not the actual truncation error in r4.
      r[i]     = r4/h;
      trunc[i] = fabs((r4 - r2)/h);                                                 // Estimated truncation error O(h)
      round[i] = fabs(e4/h) + dy;
      error[i] = round[i] + trunc[i];

      // Compute an optimised stepsize to minimize the total error, using the scaling of the estimated truncation error
      // (O(h)) and rounding error (O(1/h)).
      if ((i == 0) && ((round[0] < trunc[0]) && (round[0] > 0 && trunc[0] > 0)))
        h = (h0)*pow(round[0]/(trunc[0]), 1.0/2.0);
      else
        break;
    }
  *x = old;

  // Check that the new error is smaller, and that the new derivative is consistent with the error bounds of the original estimate.
  if ((i == 1) && (error[1] < error[0]) && (fabs(r[1] - r[0]) < 4.0*error[0]))
    *result = r[1];
  else
    *result = r[0];

  return SUCCES;
}


/*==================================================================================================================================*/

int Jacobian(const int pntdim, double *pnt, const int fncdim, double *jac, int (*fnc)(double *, double *), int method)
/*
 * Routine determines the Jacobian of the n-dimensional function F(y) w.r.t. the m-dimensional
 * variable y at the current point given by 'pnt'. The routine hence returns in 'jac' the
 * following matrix of partial derivatives:
 *
 *                    |dF1/dy1 ... dFn/dy1|
 *                    |   .          .    |
 *               Df = |   .          .    |
 *                    |   .          .    |
 *                    |dF1/dym ... dFn/dym|
 *
 * Notice that all coefficients pertaining to yi are to be found in ROW i (as opposed to column i).
 * The matrix is hence stored in column-wise (fortran) style.
 */


{
  register int                    j, k;
  double                          *dBaseMem, *y, *rhs, ydif;

  y = dBaseMem = calloc(2*pntdim, sizeof(double));
  if (!dBaseMem) return ErrorMsg("Memory allocation error in Jacobian()");

  rhs = y + pntdim;

  // Initialize
  memcpy(y, pnt, pntdim * sizeof(double));
  memset((void *)jac, 0, (pntdim*fncdim)*sizeof(double));

  for (j = 0; j < pntdim; j++)
    {
      ydif = fabs(Jacobian_Step * y[j]);
      ydif = max(ydif, JACOBIAN_MIN_STEP);

      if (FastNumerics == 1)
        {
          if (method == FORWARD)
            {
              if (ForwardDerivative(fncdim, fnc, y, rhs, y + j, ydif, rhs, jac + j*fncdim, 1) == FAILURE)
                {
                  ErrorMsg("Right-hand side computation failed");
                  free(dBaseMem);
                  dBaseMem = NULL;
                  return FAILURE;
                }
            }
          else
            {
              if (CentralDerivative(fncdim, fnc, y, rhs, y + j, ydif, rhs, jac + j*fncdim, 1) == FAILURE)
                {
                  ErrorMsg("Right-hand side computation failed");
                  free(dBaseMem);
                  dBaseMem = NULL;
                  return FAILURE;
                }
            }
        }
      else
        {
          for (k = 0; k < fncdim; k++)
            {
              if (method == FORWARD)
                {
                  if (ForwardDerivative(fncdim, fnc, y, rhs, y + j, ydif, rhs + k, jac + j*fncdim + k, 0) == FAILURE)
                    {
                      ErrorMsg("Right-hand side computation failed");
                      free(dBaseMem);
                      dBaseMem = NULL;
                      return FAILURE;
                    }
                }
              else
                {
                  if (CentralDerivative(fncdim, fnc, y, rhs, y + j, ydif, rhs + k, jac + j*fncdim + k, 0) == FAILURE)
                    {
                      ErrorMsg("Right-hand side computation failed");
                      free(dBaseMem);
                      dBaseMem = NULL;
                      return FAILURE;
                    }
                }
            }
        }
    }

  free(dBaseMem);
  dBaseMem = NULL;
  return SUCCES;
}


/*==================================================================================================================================*/

int Determinant(const int N, double *M, double *det, double *cond)

{
  char                            whichnorm;
  int                             j;
  int                             retval = SUCCES;
  double                          *A, *work, norm;
  LAPACK_SIZE_T                   nc = N, lwork = 4*N, *ipiv, *iwork, liwork = N, info;
  double                          *dBaseMem;
  LAPACK_SIZE_T                   *iBaseMem;

  // Allocate temporarily minimally allowed size for workspace arrays
  A    = dBaseMem = calloc((N*N + lwork), sizeof(double));
  if (dBaseMem == NULL) return ErrorMsg("Memory allocation error in Determinant()");

  work = A + N*N;

  ipiv  = iBaseMem = calloc(N + liwork, sizeof(LAPACK_SIZE_T));
  if (iBaseMem == NULL)
    {
      free(dBaseMem);
      dBaseMem = NULL;
      return ErrorMsg("Memory allocation error in Determinant()");
    }

  iwork = ipiv + N;

  // Copy the matrix
  memcpy(A, M, (N * N) * sizeof(double));

  dgetrf(&nc, &nc, A, &nc, ipiv, &info);
  if (info < 0)
    {
      sprintf(errstr, "Illegal value for parameter %d in dgetrf()", abs((int)info));
      ErrorMsg(errstr);
      return ILLEGAL_INPUT;
    }

  if (det)
    {
      *det = 1.0;
      if (!info)
        {
          for (j = 0; j < N; j++)
            {
              if (ipiv[j] != (j + 1))
                *det *= -A[j*N + j];
              else
                *det *= A[j*N + j];
            }
        }
    }

  if (info > 0) return SINGULARITY;

  if (cond)
    {
      norm      = anorm(N, N, M);
      whichnorm = '1';
      dgecon(&whichnorm, &nc, A, &nc, &norm, cond, work, iwork, &info FCONE);
      if (info < 0)
        {
          sprintf(errstr, "Illegal value for parameter %d in dgecon()", abs((int)info));
          ErrorMsg(errstr);
          return ILLEGAL_INPUT;
        }
    }

  free(dBaseMem);
  free(iBaseMem);
  dBaseMem = NULL;
  iBaseMem = NULL;

  return retval;
}


/*==================================================================================================================================*/

int SolveLinearSystem(const int N, double *A, double *B)

{
/*
 * This function solves the linear equation system A*x = B, where A is a NxN
 * matrix and B a N-dimensional vector.
 *
 * The matrix A has to be in Fortran column-wise format
 *
 * The solution is returned in *B, whereas the content of A is not changed.
 */

  char                            fact, trans, equed;
  double                          *Ac, *Af, *r, *c, *Bc, *x, *work, rcond, ferr = 0, berr;
  int                             retval = SUCCES;
  LAPACK_SIZE_T                   nc = N, nrhs = 1, lwork = 4*N, *ipiv, *iwork, info;
  double                          *dBaseMem;
  LAPACK_SIZE_T                   *iBaseMem;

  // Allocate temporarily minimally allowed size for workspace arrays
  Ac   = dBaseMem = calloc((2*N*N + 4*N + lwork), sizeof(double));
  if (dBaseMem == NULL) return ErrorMsg("Memory allocation error in SolveLinearSystem()");

  Ac   = dBaseMem;
  Af   = Ac + N*N;
  r    = Af + N*N;
  c    = r + N;
  Bc   = c + N;
  x    = Bc + N;
  work = x + N;

  ipiv  = iBaseMem = calloc(2*N, sizeof(LAPACK_SIZE_T));
  if (iBaseMem == NULL)
    {
      free(dBaseMem);
      dBaseMem = NULL;
      return ErrorMsg("Memory allocation error in SolveLinearSystem()");
    }

  iwork = ipiv + N;

  fact  = 'E';
  trans = 'N';

  // Fill the matrix and the right-hand side vector
  memcpy(Ac, A, (N * N) * sizeof(double));
  memcpy(Bc, B, N * sizeof(double));

  dgesvx(&fact, &trans, &nc, &nrhs, Ac, &nc, Af, &nc, ipiv, &equed, r, c, Bc, &nc, x, &nc, &rcond, &ferr, &berr, work, iwork, &info FCONE FCONE FCONE);

  // Check for singularity of the matrix
  if (info < 0)
    {
      sprintf(errstr, "Illegal value for parameter %d in dgesvx()", abs((int)info));
      ErrorMsg(errstr);
      retval = ILLEGAL_INPUT;
    }
  else if (info > 0)
    {
      ErrorMsg("(Nearly) Singular matrix in SolveLinearSystem(), solving the linear system A*x = B:\n");
      retval = SINGULARITY;
    }
  else
    {
      memcpy(B, x, N * sizeof(double));
    }

  free(dBaseMem);
  free(iBaseMem);
  dBaseMem = NULL;
  iBaseMem = NULL;

  return retval;
}


/*==================================================================================================================================*/
