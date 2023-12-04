/*
  NAME
    deBif

  PURPOSE
    Interface to the routines that are used for locating points on steady state
    bifurcation curves of non-linear ODE systems and curves determining dynamic
    regimes of such systems in two-parameter domains. More generally, the routines
    allow for computing fixed points of a system of non-linear equations.

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

    Last modification: AMdR - Feb 21, 2023
*/

#include "globals.h"

#define FEMTO                     1.0E-15

/*
 *====================================================================================================================================
 *  Definition of global variables and parameters
 *====================================================================================================================================
 */

static int                        pntDim, sysOdeDim = 0, freeparsdim;
static int                        CurveType;

static SEXP                       R_UserFunc;
static SEXP                       R_PntVals;
static SEXP                       R_FixedPars;

static long int                   N_Protected;

static int                        glorder, ninterval, CDfinemeshdim;
static double                     *CDupoldp = NULL, *CDwi = NULL, *CDwt = NULL, *CDwpvec = NULL, *CDwp = NULL;

static int                        dGlobalMemDim = 0;
static double                     *dGlobalMem = NULL;
static double                     *state0, *rhs, *ups, *ficdmat, *ficd, *xp, *t, *blockjac, *partjac, *ic;

static int                        state0Dim, rhsDim, upsRows, upsCols, ficdmatRows, ficdmatCols, ficdDim, xpRows, xpCols,
                                  tRows, tCols, blockjacRows, blockjacCols, partjacRows, partjacCols, icDim;

/*
 *====================================================================================================================================
 *  Implementation of problem specification routines
 *====================================================================================================================================
 */

int EQsystem(double *argument, double *result)

{
  SEXP                            Result;                                           // Protected
  SEXP                            R_pnt;                                            // Unprotected

  //******************************************************************************
  // Map current estimate of solution to R variables

  memcpy(REAL(R_PntVals), argument, (freeparsdim + sysOdeDim) * sizeof(double));

  R_pnt = LCONS(R_UserFunc, LCONS(R_NilValue, LCONS(R_PntVals, LCONS(R_FixedPars, R_NilValue))));
  N_Protected++;

  Result = PROTECT(R_forceAndCall(R_pnt, 3, R_GlobalEnv));

  if (!sysOdeDim) sysOdeDim = length(coerceVector(VECTOR_ELT(Result, 0), REALSXP));

  memcpy(result, REAL(coerceVector(VECTOR_ELT(Result, 0), REALSXP)), sysOdeDim * sizeof(double));
  UNPROTECT(1);
  N_Protected--;

  return SUCCES;
}


/*==================================================================================================================================*/

static inline double dot(int dim, double *x, double *y)
{
  double                          ss = 0;

  for (int ii = 0; ii < dim; ii++) ss += x[ii] * y[ii];

   return ss;
}

static inline void bialt2AI(int matdim, double *A, double *result)
{
  double                          *dblpnt = result;

  for (int p = 1; p < matdim; p++)
    for (int q = 0; q < p; q++)
      for (int r = 1; r < matdim; r++)
        for (int s = 0; s < r; s++)
          {
            if (r == q) *dblpnt = -A[p*matdim + s];
            else if ((r != p) && (s == q)) *dblpnt = A[p*matdim + r];
            else if ((r == p) && (s == q)) *dblpnt = A[p*matdim + p] + A[q*matdim + q];
            else if ((r == p) && (s != q)) *dblpnt = A[q*matdim + s];
            else if (s == p) *dblpnt = -A[q*matdim + r];
            else *dblpnt = 0.0;
            dblpnt++;
          }
  return;
}

// Computes the matrix product A x B  (A : (rowsA * cArB) matrix; B : (cArB * colsB) matrix)
static inline void matXmat(int rowsA, int cArB, int colsB, double *A, double *B, double *AxB)
{
  memset(AxB, 0, (rowsA * colsB) * sizeof(double));
  for (int ii = 0; ii < rowsA; ii++)
    for (int jj = 0; jj < colsB; jj++)
      for (int kk = 0; kk < cArB; kk++)
        AxB[ii * colsB + jj] += A[ii * cArB + kk] * B[kk * colsB + jj];

   return;
}

// Multiplies the vector V with the constant c
static inline void conXvec(int dim, double c, double *V)
{
  for (int ii = 0; ii < dim; ii++) V[ii] *= c;

   return;
}

/*==================================================================================================================================*/

int BPcondition(const int pntdim, double *y, int (*fnc)(double *, double *), int method, double *res)

/*
 * BPcondition - Routine computes the additional conditions determining the location of
 *               a branching point, see the Matcont documentation (Branch point locator,
 *               page 36, eq. 41)
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'. Notice that the
 *                        point consists of:
 *                          p:  the first bifurcation parameter
 *                          q:  the second bifurcation parameter
 *                          x:  the solution point
 *                          b:  an additional value
 *                          v:  the eigenvector (same dimension as x)
 *              y       : Pointer to an array containing the values of the unknowns
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              method  : Method to use for differential computation: FORWARD or CENTRAL
 *              res     : Pointer to the result vector to adjust and/or compute. If NULL
 *                        the routine is used to initialize the additional variables for
 *                        the BP continuation
 */
{
  int                             resindx, retcode  = SUCCES;
  double                          *dBaseMem, *Jac;
  double                          eval, *evec;

  // Some of the following vectors are over-sized, but that's OK
  Jac = dBaseMem = calloc((pntdim * pntdim), sizeof(double));
  if (!dBaseMem) return ErrorMsg("Memory allocation error in BPcondition()");

  // Determine the Jacobian of the extended system (variable plus parameter dependence).
  Jacobian(pntdim, y, sysOdeDim, Jac, fnc, method);

  /*
   * The resulting Jacobian equals the following matrix of partial derivatives:
   *
   *           |dF1/dp1 ... dFn/dp1|
   *           |dF1/dp2 ... dFn/dp2|
   *           |dF1/dy1 ... dFn/dy1|
   *      Df = |   .           .   |
   *           |   .           .   |
   *           |   .           .   |
   *           |dF1/dyn ... dFn/dyn|
   *
   * In which n = pntdim-2. Notice that all coefficients pertaining to yi are to be found
   * in ROW i + 2 (as opposed to column i). The Jacobian matrix is hence stored in column-wise
   * (fortran-style) form.
   */

  /*
   *  The extended system to solve for is:
   *
   *     F(x, p) + b*v   = 0
   *     (F_x(x, p))^T v = 0
   *     v^T F_p(x, p)   = 0
   *     v^T v - 1       = 0
   *
   *     with initial conditions b = 0 and v the eigenvector of the matrix (F_x(x, p))^T pertaining to
   *     the eigenvalue with the smallest norm. The unknowns are:
   *
   *     p:  the first bifurcation parameter
   *     q:  the second bifurcation parameter
   *     x:  the solution point
   *     b:  an additional value
   *     v:  the eigenvector (same dimension as x)
   *
   */
  // Adjust the base equations
  eval = y[freeparsdim + sysOdeDim];
  evec = y + freeparsdim + sysOdeDim + 1;

  //  F(x, p) + b*v   = 0
  for (resindx = 0; resindx < sysOdeDim; resindx++) res[resindx] += eval*evec[resindx];

  // (F_x(x, p))^T v = 0
  for (int ii = 0; ii < sysOdeDim; ii++, resindx++) res[resindx] = dot(sysOdeDim, evec, Jac + (freeparsdim + ii) * sysOdeDim);

  // v^T F_p(x, p)   = 0 : Take the second parameter in the list of arguments
  res[resindx++] = dot(sysOdeDim, evec, Jac + sysOdeDim);

  // v^T v - 1       = 0
  res[resindx++] = dot(sysOdeDim, evec, evec) - 1;

  free(dBaseMem);
  dBaseMem = NULL;

  return retcode;
}


/*==================================================================================================================================*/

int HPcondition(const int pntdim, double *y, int (*fnc)(double *, double *), int method, double *res)

/*
 * HPcondition - Routine computes the additional conditions determining the location of
 *               a Hopf bifurcation point, see equations (10.40) and (10.41) on page 485-486
 *               in Kuznetsov (1996), Elements of applied bifurcation analysis.
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'. Notice that this
 *                        equals 2 (for the bifurcation parameters) plus the dimension of
 *                        the vector of state variables
 *              y       : Pointer to an array containing as first one or two element
 *                        the value of the free parameter(s) and as subsequent elements
 *                        the values of the state variables y.
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              method  : Method to use for differential computation: FORWARD or CENTRAL
 *              res     : Pointer to the result
 */
{
  const int                       bialtdim = sysOdeDim*(sysOdeDim - 1)/2;
  int                             retcode = SUCCES;
  double                          *dBaseMem, *Jac, *Jx, *JI;

  // Some of the following vectors are over-sized, but that's OK
  Jac = dBaseMem = calloc(pntdim*pntdim + sysOdeDim*sysOdeDim + bialtdim*bialtdim, sizeof(double));
  if (!dBaseMem) return ErrorMsg("Memory allocation error in HPcondition()");

  Jx    = Jac + pntdim*pntdim;
  JI    = Jx  + sysOdeDim*sysOdeDim;

  // Determine the Jacobian of the extended system (variable plus parameter dependence).
  Jacobian(pntdim, y, sysOdeDim, Jac, fnc, method);

  /*
   * The resulting Jacobian equals the following matrix of partial derivatives:
   *
   *           |dF1/dp1 ... dFn/dp1|
   *           |dF1/dp2 ... dFn/dp2|
   *           |dF1/dy1 ... dFn/dy1|
   *      Df = |   .           .   |
   *           |   .           .   |
   *           |   .           .   |
   *           |dF1/dyn ... dFn/dyn|
   *
   * In which n = pntdim-2. Notice that all coefficients pertaining to yi are to be found
   * in ROW i + 2 (as opposed to column i). The Jacobian matrix is hence stored in column-wise
   * (fortran-style) form. Below it is transposed into a row-wise (C-style) form.
   */

  // Extract the restricted Jacobian in transposed (row-wise, C-style) form
  for (int ii = 0; ii < sysOdeDim; ii++)
    for (int jj = 0; jj < sysOdeDim; jj++) Jx[ii*sysOdeDim + jj] = Jac[(freeparsdim + jj)*sysOdeDim + ii];

  // Now construct the bialternate matrix product of 2J and I
  bialt2AI(sysOdeDim, Jx, JI);

  retcode = Determinant(bialtdim, JI, res + sysOdeDim, NULL);
  if ((retcode != SUCCES) && (retcode != SINGULARITY))
    {
      ErrorMsg("Failed to compute determinant of bialternate matrix product in HPcondition()");
      free(dBaseMem);
      dBaseMem = NULL;
      return retcode;
    }

  free(dBaseMem);
  dBaseMem = NULL;

  return retcode;
}


/*==================================================================================================================================*/

int LPcondition(const int pntdim, double *y, int (*fnc)(double *, double *), const int method, double *res)

/*
 * LPcondition -  Routine computes the factor determining the location of a limit point, i.e.
 *                the parameter component of the tangent vector. This component always has
 *                index 0 in the vector of the solution point.
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'. Notice that this
 *                        equals 2 (for the bifurcation parameters) plus the dimension of
 *                        the vector of state variables
 *              y       : Pointer to an array containing as first one or two element
 *                        the value of the free parameter(s) and as subsequent elements
 *                        the values of the state variables y.
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              method  : Method to use for differential computation: FORWARD or CENTRAL
 *
 *   Continuation of a limitpoint is carried out using the defining system
 *   (10.97) on page 515 of Kuznetsov (1998):
 *
 *    F(x, p)                 = 0
 *    F_x(x, p) q             = 0
 *    (F_x(x, p))^T p - eps*p = 0
 *    q^T q - 1               = 0
 *    p^T p - 1               = 0
 *
 *    with initial conditions b = 0 and v the eigenvector of the matrix
 *    (F_x(x, p))^T pertaining to the eigenvalue with the smallest norm.
 *    The unknowns are:
 *
 *    p:   the bifurcation parameter
 *    x:   the solution point
 *    eps: an additional value
 *    q:   the eigenvector of F_x pertaining to the 0 eigenvalue
 *    p:   the eigenvector of (F_x)^T pertaining to the 0 eigenvalue
 *
 *   The advantage of this defining system is that detection of Bogdanov-Takens
 *   points and cusp point are straightforward
 */
{
  register int                    resindx;
  double                          *dBaseMem, *Jac, *Jx, *JxT, eps, *pv, *qv;

  Jac = dBaseMem = calloc((2 * pntdim * sysOdeDim), sizeof(double));
  if (!dBaseMem) return ErrorMsg("Memory allocation error in LPcondition()");

  Jx = Jac + pntdim * sysOdeDim;

  // Determine the Jacobian of the extended system (variable plus parameter dependence).
  Jacobian(pntdim, y, sysOdeDim, Jac, fnc, method);

  /*
   * The Jacobian equals the following (n+2)x(n) matrix of partial derivatives:
   *
   *           |dF1/dp1 ... dFn/dp1|
   *           |dF1/dp2 ... dFn/dp2|
   *           |dF1/dy1 ... dFn/dy1|
   *      Df = |   .           .   |
   *           |   .           .   |
   *           |   .           .   |
   *           |dF1/dyn ... dFn/dyn|
   *
   * In which n = pntdim-2 (i.e. equal to the number of state variables).
   * Notice that all coefficients pertaining to yi are to be found in ROW i + 2
   * (as opposed to column i). The matrix is hence stored in column-wise (fortran) style.
   */

  // Extract the restricted Jacobian in transposed (row-wise, C-style) form
  for (int ii = 0; ii < sysOdeDim; ii++)
    for (int jj = 0; jj < sysOdeDim; jj++) Jx[ii * sysOdeDim + jj] = Jac[(freeparsdim + jj)*sysOdeDim + ii];

  // The transposed and restricted Jacobian of the system is found at row 2..(sysOdeDim+2)
  JxT = Jac + freeparsdim * sysOdeDim;

  // Extract the additional variables from the state
  eps  = y[freeparsdim + sysOdeDim];
  qv   = y + freeparsdim +     sysOdeDim + 1;
  pv   = y + freeparsdim + 2 * sysOdeDim + 1;

  resindx = sysOdeDim;
  // Add the additional values
  // A qv = 0
  for (int ii = 0; ii < sysOdeDim; ii++, resindx++) res[resindx] = dot(sysOdeDim, qv, Jx + ii * sysOdeDim);

  // A^T pv - eps*pv = 0
  for (int ii = 0; ii < sysOdeDim; ii++, resindx++) res[resindx] = dot(sysOdeDim, pv, JxT + ii * sysOdeDim) - eps * pv[ii];

  // <qv, qv) - 1 =
  res[resindx++] = dot(sysOdeDim, qv, qv) - 1;

  // <pv, pv) - 1 =
  res[resindx++] = dot(sysOdeDim, pv, pv) - 1;

  free(dBaseMem);
  dBaseMem = NULL;

  return SUCCES;
}


/*==================================================================================================================================*/

static void ExtSystemLCblockjac(double *xp, double *state0, double lcperiod, double *partjacout, int method)

// ExtSystemLCblockjac <- function(xp, state, parms, curveData, nopts = NULL) {

{
  int                             CDwpRows  = (sysOdeDim * (glorder + 1));
  int                             CDwpCols  = (sysOdeDim * (glorder));
  double                          *sysjac;

  sysjac = blockjac + sysOdeDim;

  /* The first term on the right-hand side of the blockjac assignment below amounts to a direct copy
   * of CDwp into rows 1..(1 + CDwpRows) of partjac
   *
   *   wploc = curveData$wp/dt;
   *   blockjac[range1,1+(1:(jaccol-2))] <- wploc[range1,] - state["LCperiod"]*fastkron(glorder, sysOdeDim, t(CDwt[,j]), sysjac)
   */
  double dt = 1.0 / ninterval;

  for (int ii = 0; ii < (CDwpRows * CDwpCols); ii++) partjacout[partjacCols + ii] = CDwp[ii] / dt;

  // Evaluate function value on each collocation point
  for (int ii = 0; ii < glorder; ii++)
    {
      memcpy(state0 + 1, xp + ii * sysOdeDim, sysOdeDim * sizeof(double));
      memset(blockjac, 0, (blockjacRows * blockjacCols) * sizeof(double));
      Jacobian((freeparsdim + sysOdeDim), state0, sysOdeDim, blockjac, EQsystem, method);

      // Derivatives w.r.t. the bifurcation parameter in first row
      memcpy(partjacout + (ii * sysOdeDim), blockjac, sysOdeDim * sizeof(double));

      // blockjac[range1,1+(1:(jaccol-2))] <- wploc[range1,] - state["LCperiod"]*fastkron(glorder, sysOdeDim, t(CDwt[,j]), sysjac)
      for (int jj = 0; jj < (glorder + 1); jj++)
        for(int kk = 0; kk < sysOdeDim; kk++)
          for(int ll = 0; ll < sysOdeDim; ll++)
            {
              partjacout[(1 + jj * sysOdeDim + kk) * partjacCols + (ii * sysOdeDim) + ll]  -= lcperiod * CDwt[ii * (glorder + 1) + jj] * sysjac[kk * sysOdeDim + ll];
            }

      // Derivatives w.r.t. the bifurcation parameter in last row
      EQsystem(state0, partjacout + ((partjacRows - 1) * partjacCols) + (ii * sysOdeDim));
    }

  // Multiply the parameter derivatives with -state["LCperiod"]
  for (int jj = 0; jj < partjacCols; jj++) partjacout[jj] *= -lcperiod;

  // Multiply the last row with -1
  for (int jj = 0, kk = ((partjacRows - 1) * partjacCols); jj < partjacCols; jj++) partjacout[kk + jj] *= -1;

  return;
}


/*==================================================================================================================================*/

static int ExtSystemLCjac(const int pntdim, double *y, const int fncdim, double *fulljac, int (*fnc)(double *, double *), int method)

{
  double                          dt = 1.0 / ((double)ninterval);

  memset(dGlobalMem, 0, dGlobalMemDim * sizeof(double));
  memcpy(state0, y, (freeparsdim + sysOdeDim) * sizeof(double));                    // Copy bifurcation parameter and state at t=0

  // In LCfunc.R the limit cycle period is added, but that seems unnecessary
  // state0[(freeparsdim + sysOdeDim)] = y[pntdim -1];
  double LCperiod = y[pntdim -1];

  // Extract a matrix with state variable values.
  // In R:  ups: sysOdeDim rows; finemesh  columns (fortran column-wise storage)
  // In C:  ups: finemesh  rows; sysOdeDim columns (row-wise storage)
  memcpy(ups, y + freeparsdim, (upsRows * upsCols) * sizeof(double));

  /* The jacobian has pntdim rows and fncdim columns and has the following layout:
   *
   *           |dF1/dp1   ... dFn/dp1  |
   *           |dF1/dy1   ... dFn/dy1  |
   *           |   .             .     |
   *      Df = |   .             .     |
   *           |   .             .     |
   *           |dF1/dyn-1 ... dFn/dyn-1|
   *           |dF1/dyn   ... dFn/dyn  |
   *
   * Here p1 is the bifurcation parameter and y1..yn-1 are the nodal points
   * at which the LC is approximated (dimension: finemeshdim * sysOdeDim)
   * and yn is the LC period.
   *
   * F1 .. Fn are the values returned by LCcondition, which include:
   *
   *  - (ninterval * glorder * sysOdeDim) values determined by the
   *    right-hand side of the model at the nodal points at which the LC
   *    is approximated
   *  - sysOdeDim values determined by the circular boundary condition
   *  - 1 value (the last) corresponding to the integral condition
   *
   */

  memset(fulljac, 0, pntdim * fncdim * sizeof(double));

  // Evaluate the block jacobians
  for (int ii = 0; ii < ninterval; ii++)
    {
      // Value of polynomial on each collocation point
      //     xp <- ups[, range1] %*% curveData$wt
      // Notice that CDwt (curveData$wt) should be stored as a matrix with glorder
      // rows and (glorder + 1) columns
      matXmat(xpRows, (glorder + 1), xpCols, CDwt, ups + ii * glorder * sysOdeDim, xp);

      /* partjac dimensions:
       *
       * Columns:   partjacCols = (sysOdeDim * glorder) (= blockrow in R)
       * Rows   :   partjacRows = (sysOdeDim * (glorder + 1) + 2) rows (= blockcol in R)
       *
       * Notice that just like for fulljac, in partjac all the derivative values
       * w.r.t. a particular variable should be in the same ROW.
       */

      // partjac <- ExtSystemLCblockjac(xp, state0, parms, curveData, nopts)
      memset(partjac, 0, (partjacRows * partjacCols) * sizeof(double));
     ExtSystemLCblockjac(xp, state0, LCperiod, partjac, method);

      // The derivative w.r.t. the bifurcation parameter
      memcpy(fulljac + ii * partjacCols, partjac, partjacCols * sizeof(double));

      //  fulljac[rowrange, 1 + colrange] <- partjac[(1:blockrow), 1 + (1:blockcol)]
      for (int jj = 0; jj < sysOdeDim * (glorder+1); jj++)
        memcpy(fulljac + (1 + ii * (sysOdeDim * glorder) + jj) * fncdim + ii * partjacCols, partjac + (1 + jj) * partjacCols, partjacCols * sizeof(double));

      // Last row of jacobian
      memcpy(fulljac + (pntdim - 1) * fncdim + ii * partjacCols, partjac + (partjacRows - 1) * partjacCols, partjacCols * sizeof(double));

      // Derivative of the integral constraint
      // p <- dt*(CDupoldp[,range1]*curveData$pwi)                                  //pwi is sysOdeDim stacked rows of wi
      // ic[1 + colrange] <- ic[1 + colrange] + p[1:blockcol]
      for (int jj =0; jj < (glorder + 1); jj++)
        for (int kk = 0; kk < sysOdeDim; kk++)
          ic[1 + (ii * glorder + jj) * sysOdeDim + kk] += dt * (CDupoldp[(ii * glorder + jj) * sysOdeDim + kk] * CDwi[jj]);
    }

  // Derivative of the boundary condition into the one-to-last column of fulljac
  for (int ii = 0; ii < sysOdeDim; ii++)
    {
      fulljac[(ii + 1) * fncdim + (fncdim - (sysOdeDim + 1) + ii)] = 1;                            // Derivative w.r.t. ups[kk] in row 2 and further
      fulljac[(pntDim - (sysOdeDim + 1) + ii) * fncdim + (fncdim - (sysOdeDim + 1) + ii)] = -1;   // Derivative w.r.t. ups[(CDfinemeshdim - 1) * sysOdeDim + kk] in end rows (but not the last)
    }

  // Derivative of the integral constraint into the last column of fulljac
  for (int ii = 0; ii < pntdim; ii++) fulljac[ii * fncdim + (fncdim - 1)] = ic[ii];

  return SUCCES;
}


/*==================================================================================================================================*/

static int LCcondition(double *argument, double *result)

// ExtSystemLC <- function(t, state, parms, curveData, nopts = NULL) {

{

  memset(dGlobalMem, 0, dGlobalMemDim * sizeof(double));
  memcpy(state0, argument, (freeparsdim + sysOdeDim) * sizeof(double));                    // Copy bifurcation parameter and state at t=0

  // In LCfunc.R the limit cycle period is added, but that seems unnecessary
  // state0[(freeparsdim + sysOdeDim)] = argument[pntDim -1];
  double LCperiod = argument[pntDim -1];

  // Extract a matrix with state variable values.
  // In R:  ups: sysOdeDim rows; finemesh  columns (fortran column-wise storage)
  // In C:  ups: finemesh  rows; sysOdeDim columns (row-wise storage)
  memcpy(ups, argument + freeparsdim, (upsRows * upsCols) * sizeof(double));

  // Compute the values for the integral condition: element-wise multiplication
  // of ups and upoldp and summing over the columns
  memset(ficd, 0, CDfinemeshdim * sizeof(double));
  for (int ii = 0; ii < upsRows; ii++)
    for (int jj = 0; jj < upsCols; jj++)
      ficd[ii] += ups[ii * upsCols + jj] * CDupoldp[ii * upsCols + jj];

  // Evaluate the derivatives at all the nodal points of the mesh
  double dt = 1.0 / ((double)ninterval);
  for (int ii = 0; ii < ninterval; ii++)
    {
      // Value of polynomial on each collocation point
      //     xp <- ups[, range1] %*% curveData$wt
      // Notice that CDwt (curveData$wt) should be stored as a matrix with glorder
      // rows and (glorder + 1) columns
      matXmat(xpRows, (glorder + 1), xpCols, CDwt, ups + ii * glorder * sysOdeDim, xp);

      // Derivative of polynomial on each collocation point
      // t  <- (ups[, range1] %*% CDwpvec)/dt
      // Notice that CDwpvec (curveData$wpvec) should be stored as a matrix with glorder
      // rows and (glorder + 1) columns
      matXmat(tRows, (glorder + 1), tCols, CDwpvec, ups + ii * glorder * sysOdeDim, t);
      conXvec((tRows  * tCols), ninterval, t);

      // Evaluate function value on each collocation point
      for (int jj = 0; jj < glorder; jj++)
        {
          memcpy(state0 + 1, xp + jj * sysOdeDim, sysOdeDim * sizeof(double));
          EQsystem(state0, rhs);

          for (int kk = 0; kk < sysOdeDim; kk++)
            result[(ii * glorder + jj) * sysOdeDim + kk] = t[jj * tCols + kk] - LCperiod * rhs[kk];
        }

      // Put the appropriate values of the integral condition in the matrix
      memcpy(ficdmat + ii * ficdmatCols, ficd + ii * glorder, ficdmatCols * sizeof(double));
    }

  // Ciruclar boundary conditions
  for (int kk = 0; kk < sysOdeDim; kk++)
    result[(CDfinemeshdim - 1) * sysOdeDim + kk] = ups[kk] - ups[(CDfinemeshdim - 1) * sysOdeDim + kk];

  //  Integral constraint
  result[CDfinemeshdim * sysOdeDim] = 0.0;
  for (int ii = 0; ii < ninterval; ii++)
    result[CDfinemeshdim * sysOdeDim] += dt * dot(ficdmatCols, CDwi, ficdmat + ii * ficdmatCols);

  return SUCCES;
}


/*==================================================================================================================================*/

int AllEquations(double *argument, double *result)

{
  int                             retval = SUCCES;

  // Compute the basic system of equations
  EQsystem(argument, result);

  //==================================================================================================================================
  // Add the final value in case of BP, HP or LP continuation

  if (CurveType == BP)
    retval = BPcondition(pntDim, argument, EQsystem, CENTRAL, result);
  else if (CurveType == HP)
    retval = HPcondition(pntDim, argument, EQsystem, CENTRAL, result);
  else if (CurveType == LP)
    retval = LPcondition(pntDim, argument, EQsystem, CENTRAL, result);

  return retval;
}


/*==================================================================================================================================*/

SEXP deBif(SEXP curveType, SEXP userFunc, SEXP initVals, SEXP fixedParVals, SEXP tanVec,
           SEXP rhsTol, SEXP varTol, SEXP jacStep, SEXP maxIter, SEXP glOrder, SEXP nInterval, SEXP cData)

{
  int                             MaxIter, retcode = FAILURE, listel = 0, nIter = 10;
  double                          RhsTol, VarTol;
  double                          *dBaseMem, *point, *tanvec, *JacExport, *dblPnt;
  SEXP                            outputList = R_NilValue, nms, outputListEl[4], R_VarNames;
  SEXP                            R_VarNamesLC, R_cDataNames = R_NilValue;

  N_Protected = 0L;
  glorder = ninterval = CDfinemeshdim = 0;

  //============================== Process the curve type argument ===================================================================

  if (!strcmp(CHAR(STRING_ELT(curveType, 0)), "BP"))      CurveType = BP;
  else if (!strcmp(CHAR(STRING_ELT(curveType, 0)), "EQ")) CurveType = EQ;
  else if (!strcmp(CHAR(STRING_ELT(curveType, 0)), "HP")) CurveType = HP;
  else if (!strcmp(CHAR(STRING_ELT(curveType, 0)), "LP")) CurveType = LP;
  else if (!strcmp(CHAR(STRING_ELT(curveType, 0)), "LC")) CurveType = LC;

  if ((CurveType == EQ) || (CurveType == LC)) freeparsdim = 1;
  else freeparsdim = 2;

  //============================== Process the user funcion argument =================================================================

  R_UserFunc = PROTECT(duplicate(userFunc));
  N_Protected++;

  //============================== Process the numerical option argument =============================================================

  RhsTol        = asReal(rhsTol);
  VarTol        = asReal(varTol);
  Jacobian_Step = asReal(jacStep);
  MaxIter       = asInteger(maxIter);
  glorder       = asInteger(glOrder);
  ninterval     = asInteger(nInterval);

  if (CurveType == LC)
    {
      char    optname[MAX_STR_LEN];
      int     ncols = length(cData);

      // Get the dimensions from curveData
      R_cDataNames = PROTECT(getAttrib(cData, R_NamesSymbol));
      N_Protected++;

      for (int ii = 0, nset = 0; ii < ncols; ii++)
        {
          strcpy(optname, CHAR(STRING_ELT(R_cDataNames, ii)));
          if (strcmp(optname, "finemeshdim") == 0)
            {
              CDfinemeshdim = asInteger(VECTOR_ELT(cData, ii));
              nset++;
            }
          else if (strcmp(optname, "statedim") == 0)
            {
              sysOdeDim = asInteger(VECTOR_ELT(cData, ii));
              nset++;
            }
          if (nset == 2) break;
        }
    }

  //============================== Process the initial point argument ================================================================

  pntDim = length(initVals);

  point = dBaseMem = calloc((2 * pntDim + (pntDim * pntDim)), sizeof(double));      // point, tanvec, JacExport
  if (!dBaseMem)
    {
      ErrorMsg("Memory allocation error in deBif()");
      UNPROTECT(N_Protected);

      return outputList;
    }

  tanvec    = point     + pntDim;
  JacExport = tanvec    + pntDim;

  memcpy(point, REAL(initVals), pntDim * sizeof(double));

  R_VarNames = PROTECT(getAttrib(initVals, R_NamesSymbol));
  N_Protected++;

  if (CurveType == LC)
    {
      char    optname[MAX_STR_LEN];
      int     ncols = length(cData);

      // Allocate the global memory
      state0Dim    = freeparsdim + sysOdeDim + 1;
      rhsDim       = sysOdeDim;
      upsRows      = CDfinemeshdim;
      upsCols      = sysOdeDim;
      ficdmatRows  = ninterval;
      ficdmatCols  = glorder+1;
      ficdDim      = CDfinemeshdim;
      xpRows       = glorder;
      xpCols       = sysOdeDim;
      tRows        = glorder;
      tCols        = sysOdeDim;
      blockjacRows = freeparsdim + sysOdeDim;
      blockjacCols = sysOdeDim;
      partjacRows  = sysOdeDim * (glorder+1) + freeparsdim + 1;
      partjacCols  = sysOdeDim * glorder;
      icDim        = pntDim;

      dGlobalMemDim  = 0;
      dGlobalMemDim += state0Dim;
      dGlobalMemDim += rhsDim;
      dGlobalMemDim += upsRows      * upsCols;
      dGlobalMemDim += ficdmatRows  * ficdmatCols;
      dGlobalMemDim += ficdDim;
      dGlobalMemDim += xpRows       * xpCols;
      dGlobalMemDim += tRows        * tCols;
      dGlobalMemDim += blockjacRows * blockjacCols;
      dGlobalMemDim += partjacRows  * partjacCols;
      dGlobalMemDim += icDim;

      int dGlobalMemDim2 = 0;
      dGlobalMemDim2 += (pntDim - 2);                                                // CDupoldp (exclude par & period)
      dGlobalMemDim2 += (glorder + 1);                                               // CDwi
      dGlobalMemDim2 += (glorder)                   * (glorder + 1);                 // CDwt
      dGlobalMemDim2 += (glorder)                   * (glorder + 1);                 // CDwpvec
      dGlobalMemDim2 += (sysOdeDim * (glorder + 1)) * (sysOdeDim * glorder);         // CDwp

      state0 = dGlobalMem = calloc(dGlobalMemDim + dGlobalMemDim2, sizeof(double));
      if (!dGlobalMem)
        {
          ErrorMsg("Memory allocation error in deBif()");
          UNPROTECT(N_Protected);

          free(dBaseMem);
          dBaseMem = NULL;
          dGlobalMem = NULL;
          return outputList;
        }

      rhs       = state0    + state0Dim;
      ups       = rhs       + rhsDim;
      ficdmat   = ups       + upsRows      * upsCols;
      ficd      = ficdmat   + ficdmatRows  * ficdmatCols;
      xp        = ficd      + ficdDim;
      t         = xp        + xpRows       * xpCols;
      blockjac  = t         + tRows        * tCols;
      partjac   = blockjac  + blockjacRows * blockjacCols;
      ic        = partjac   + partjacRows  * partjacCols;

      CDupoldp  = ic        + icDim;
      CDwi      = CDupoldp  + (pntDim - 2);
      CDwt      = CDwi      + (glorder + 1);
      CDwpvec   = CDwt      + ((glorder) * (glorder + 1));
      CDwp      = CDwpvec   + ((glorder) * (glorder + 1));

      for (int ii = 0, nset = 0; ii < ncols; ii++)
        {
          strcpy(optname, CHAR(STRING_ELT(R_cDataNames, ii)));
          if (strcmp(optname, "upoldp") == 0)
            {
              memcpy(CDupoldp, REAL(VECTOR_ELT(cData, ii)), (pntDim - 2) * sizeof(double));
              nset++;
            }
          else if (strcmp(optname, "wi") == 0)
            {
              memcpy(CDwi, REAL(VECTOR_ELT(cData, ii)), (glorder + 1) * sizeof(double));
              nset++;
            }
          else if (strcmp(optname, "wt") == 0)
            {
              memcpy(CDwt, REAL(VECTOR_ELT(cData, ii)), ((glorder) * (glorder + 1)) * sizeof(double));
              nset++;
            }
          else if (strcmp(optname, "wpvec") == 0)
            {
              memcpy(CDwpvec, REAL(VECTOR_ELT(cData, ii)), ((glorder) * (glorder + 1)) * sizeof(double));
              nset++;
            }
          else if (strcmp(optname, "wp") == 0)
            {
              int matsize = ((sysOdeDim * (glorder + 1)) * (sysOdeDim * glorder));
              memcpy(CDwp, REAL(VECTOR_ELT(cData, ii)), matsize * sizeof(double));
              nset++;
            }
          if (nset == 5) break;
        }

      R_PntVals = PROTECT(allocVector(REALSXP, (freeparsdim + sysOdeDim + 1)));
      N_Protected++;

      // create names for single state vector
      R_VarNamesLC = PROTECT(allocVector(STRSXP, (freeparsdim + sysOdeDim + 1)));
      N_Protected++;

      for (int ii = 0; ii < (freeparsdim + sysOdeDim); ii++)
        SET_STRING_ELT(R_VarNamesLC, ii, STRING_ELT(R_VarNames, ii));

      SET_STRING_ELT(R_VarNamesLC, (freeparsdim + sysOdeDim), mkChar("LCperiod"));

      setAttrib(R_PntVals, R_NamesSymbol, R_VarNamesLC);

      UNPROTECT(1);
      N_Protected--;
    }
  else
    {
      R_PntVals = PROTECT(duplicate(initVals));
      N_Protected++;
    }

  //================================ Process the parameters argument =================================================================

  R_FixedPars = PROTECT(duplicate(fixedParVals));
  N_Protected++;

  //================================ Process the tangent argument ====================================================================

  memcpy(tanvec, REAL(tanVec), pntDim * sizeof(double));

  //============================= Compute the solution point =========================================================================

  retcode = FAILURE;
  if (CurveType == LC)
    {
      retcode = FindPoint(pntDim, freeparsdim, point, tanvec, RhsTol, VarTol, MaxIter, &nIter, LCcondition, ExtSystemLCjac);
    }
  else
    {
      retcode = FindPoint(pntDim, freeparsdim, point, tanvec, RhsTol, VarTol, MaxIter, &nIter, AllEquations, Jacobian);
    }

  //============================= Return the solution point =========================================================================

  if (retcode == SUCCES)
    {
      PROTECT(outputListEl[listel] = allocVector(REALSXP, pntDim));
      N_Protected++;

      dblPnt = REAL(outputListEl[listel]);
      memcpy(dblPnt, point, pntDim * sizeof(double));

      setAttrib(outputListEl[listel], R_NamesSymbol, R_VarNames);
      listel++;

      PROTECT(outputListEl[listel] = allocVector(INTSXP, 1));
      N_Protected++;

      // return the number of iterations
      INTEGER(outputListEl[listel])[0] = nIter;
      listel++;

      if (CurveType == LC)
        retcode = TangentVec(pntDim, point, JacExport, tanvec, LCcondition, ExtSystemLCjac, NULL);
      else
        retcode = TangentVec(pntDim, point, JacExport, tanvec, AllEquations, Jacobian, NULL);

      if (retcode == SUCCES)
        {
          PROTECT(outputListEl[listel] = allocVector(REALSXP, pntDim));
          N_Protected++;

          dblPnt = REAL(outputListEl[listel]);
          memcpy(dblPnt, tanvec, pntDim * sizeof(double));

          setAttrib(outputListEl[listel], R_NamesSymbol, R_VarNames);
          listel++;

          PROTECT(outputListEl[listel] = allocMatrix(REALSXP, pntDim, pntDim));
          N_Protected++;

          dblPnt = REAL(outputListEl[listel]);
          // memcpy(dblPnt, JacExport, pntDim * sysOdeDim * sizeof(double));
          for (int ii = 0; ii < pntDim; ii++)
            {
              memcpy((dblPnt + ii * pntDim), (JacExport + ii * (pntDim - 1)), (pntDim - 1) * sizeof(double));
              dblPnt[ii * pntDim + pntDim - 1] = tanvec[ii];
            }
          listel++;
        }

      outputList = PROTECT(allocVector(VECSXP, listel));
      N_Protected++;

      // create names
      nms = PROTECT(allocVector(STRSXP, listel));
      N_Protected++;

      switch (listel)
        {
          case 4:
            SET_STRING_ELT(nms, 3, mkChar("Jacobian"));
          case 3:
            SET_STRING_ELT(nms, 2, mkChar("tanvec"));
          case 2:
            SET_STRING_ELT(nms, 1, mkChar("niter"));
          default:
            SET_STRING_ELT(nms, 0, mkChar("y"));
        }

      for (int ii = 0; ii < listel; ii++) SET_VECTOR_ELT(outputList, ii, outputListEl[ii]);
      setAttrib(outputList, R_NamesSymbol, nms);
    }

  UNPROTECT(N_Protected);

  if (dGlobalMem) free(dGlobalMem);
  free(dBaseMem);
  dBaseMem = NULL;
  dGlobalMem = NULL;

  return outputList;
}
