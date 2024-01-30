/*CPD_JACOBIAN_INCOMPLETE CPD Jacobian of incomplete tensor.
 *  See cpd_jacobian_incomplete.m for documentation.
 *
 * Author(s):  Nico Vervliet       (Nico.Vervliet@kuleuven.be)
 *             Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
 *
 * Version History:
 * - 2017/04/04   NV      Complex data and order(U) > order(T) added.
 * - 2017/03/06   NV      Initial version.
 */

/* INCLUDES -------------------------------------------------------------------*/
#include "mex.h"
#include "matrix.h"

/* Computational routines -----------------------------------------------------*/

void fill_jacobian_real (const mxArray *jacobian,  /* Pointer to J */
                         const mxArray *factors,   /* Pointer to transposed
                                                      factors matrices U */
                         const mxArray *subs,      /* Pointer to indices of
                                                      known entries */
                         const size_t NZ,          /* Length of U */
                         const size_t NT,          /* Order of T (<= NZ !) */
                         const size_t R,           /* Rank = rows of U{n} */
                         const size_t nke)         /* Number of known entries */    
{

  /* Auxiliary variables */ 
  size_t m, n, r, i, ind;
  double tmp[NZ][R];

  /* Pointers to data */
  double *jacobianptr = mxGetPr(jacobian);
  double *factorptr[NZ];
  int    *subptr[NT];


  for (n = 0; n < NT; n++)
    subptr[n] = (int *) mxGetData(mxGetCell(subs,n));
  
  for (n = 0; n < NZ; n++) 
    factorptr[n] = mxGetPr(mxGetCell(factors,n));
  
  /* Loop over all known entries */
  for (i = 0; i < nke; i++) {
    /* Retrieve relevant factor rows from memory */
    for (n = 0; n < NT; n++) {
      for (r = 0; r < R; r++) {
        ind = subptr[n][i];
        tmp[n][r] = factorptr[n][(ind-1)*R + r];
      }
    }
    for (n = NT; n < NZ; n++) {
      for (r = 0; r < R; r++) {
        tmp[n][r] = factorptr[n][r];
      }
    }

    /* Create Jacobian entries using element-wise products */
    for (n = 0; n < NZ; n++) {
      for (r = 0; r < R; r++) {
        *jacobianptr = 1;
        for (m = 0; m < NZ; m++) {
          if (m == n)
            continue;
          *jacobianptr *= tmp[m][r];
        } 
        jacobianptr++;
      } /* end r loop */
    } /* end n loop */
  } /* end i loop */
}

void fill_jacobian_complex (const mxArray *jacobian,  /* Pointer to J */
                            const mxArray *factors,   /* Pointer to transposed
                                                         factors matrices U */
                            const mxArray *subs,      /* Pointer to indices of
                                                         known entries */
                            const size_t NZ,          /* Length of U */
                            const size_t NT,          /* Order of T (<= NZ !) */
                            const size_t R,           /* Rank = rows of U{n} */
                            const size_t nke)         /* Number of known entries */  
{

  /* Auxiliary variables */ 
  size_t m, n, r, i, ind;
  double tmpR[NZ][R];
  double tmpI[NZ][R];
  double a, b;

  /* Pointers to data */
  double *jacobianptrR = mxGetPr(jacobian);
  double *factorptrR[NZ];
  double *jacobianptrI = mxGetPi(jacobian);
  double *factorptrI[NZ];
  int    *subptr[NT];

  for (n = 0; n < NT; n++)
    subptr[n] = (int *) mxGetData(mxGetCell(subs,n));
  
  for (n = 0; n < NZ; n++) { 
    factorptrR[n] = mxGetPr(mxGetCell(factors,n));
    factorptrI[n] = mxGetPi(mxGetCell(factors,n));
  }

  /* Loop over all known entries */
  for (i = 0; i < nke; i++) {
    /* Retrieve relevant factor rows from memory */
    for (n = 0; n < NT; n++) {
      for (r = 0; r < R; r++) {
        ind = subptr[n][i];
        tmpR[n][r] = factorptrR[n][(ind-1)*R + r];
        tmpI[n][r] = factorptrI[n][(ind-1)*R + r];
      }      
    }
    for (n = NT; n < NZ; n++) {
      for (r = 0; r < R; r++) {
        tmpR[n][r] = factorptrR[n][r];
        tmpI[n][r] = factorptrI[n][r];
      }
    }

    /* Create Jacobian entries using element-wise products */
    for (n = 0; n < NZ; n++) {
      for (r = 0; r < R; r++) {
        *jacobianptrR = 1;
        *jacobianptrI = 0;
        for (m = 0; m < NZ; m++) {
          if (m == n)
            continue;
          a = *jacobianptrR;
          b = *jacobianptrI;
          *jacobianptrR = a * tmpR[m][r] - b*tmpI[m][r];
          *jacobianptrI = b * tmpR[m][r] + a*tmpI[m][r];
        } 
        jacobianptrR++;
        jacobianptrI++;
      } /* end r loop */
    } /* end n loop */
  } /* end i loop */
}

/*CPD_JACOBIAN_INCOMPLETE CPD Jacobian of incomplete tensor.
 *  Gateway routine for matlab. 
 *
 *  See cpd_jacobian_incomplete.m for documentation.
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (nrhs < 3)
    mexErrMsgIdAndTxt("cpd_jacobian_incomplete:input", "The number of inputs should be 3.\n");
  if (!mxIsDouble(prhs[0]))
    mexErrMsgIdAndTxt("cpd_jacobian_incomplete:input", "Jacobian should be a matrix.\n");
  if (!mxIsCell(prhs[1]))
    mexErrMsgIdAndTxt("cpd_jacobian_incomplete:input", "Factors should be a cell.\n");
  if (!mxIsCell(prhs[2]))
    mexErrMsgIdAndTxt("cpd_jacobian_incomplete:input", "Subs should be a cell.\n");
  
  size_t NZ  = mxGetNumberOfElements(prhs[1]);
  size_t NT  = mxGetNumberOfElements(prhs[2]);
  size_t nke = mxGetNumberOfElements(mxGetCell(prhs[2],0));
  size_t R   = mxGetM(mxGetCell(prhs[1],0));
  size_t n;

  if (NT > NZ)
    mexErrMsgIdAndTxt("cpd_jacobian_incomplete:input", "Length(subs) should smaller than or equal to length(factors).\n");

  size_t iscomplex = mxIsComplex(mxGetCell(prhs[1],0));
  for (n = 1; n < NZ; n++) {
    /* Test factors */
    if (mxGetM(mxGetCell(prhs[1],n)) != R)
      mexErrMsgIdAndTxt("cpd_jacobian_incomplete:input", "The number of rows of the (transposed!) factors should be R for all n.\n");
    if (mxIsComplex(mxGetCell(prhs[1],0)) != iscomplex)
      mexErrMsgIdAndTxt("cpd_jacobian_incomplete:input", "All factor matrices should be either real or complex.\n");
  }

  for (n = 1; n < NT; n++) {
    if (mxGetNumberOfElements(mxGetCell(prhs[2],n)) != nke)
      mexErrMsgIdAndTxt("cpd_jacobian_incomplete:input", "The length(subs(n)) should be equal for all n.\n");
    if (!mxIsInt32(mxGetCell(prhs[2],n)))
      mexErrMsgIdAndTxt("cpd_jacobian_incomplete:input", "subs{n} should be int32 for all n.\n");
  }

  /* Test number of entries in jacobian */
  if (mxGetJc(prhs[0])[mxGetN(prhs[0])] != nke*R*NZ)
    mexErrMsgIdAndTxt("cpd_jacobian_incomplete:input", "The number of entries in Jacobian should match the length of subs{n}*R*N.\n");

  if (mxIsComplex(prhs[0]) != iscomplex)
    mexErrMsgIdAndTxt("cpd_jacobian_incomplete:input", "If U{n} are real (complex), J should be real (complex) as well.\n");

  /* Test output: should be none */
  if (nlhs > 0)
    mexErrMsgIdAndTxt("cpd_jacobian_incomplete:output", "No output should be requested.\n");

  /* Call algorithm */
  if (mxIsComplex(mxGetCell(prhs[1],0)))
    fill_jacobian_complex(prhs[0], prhs[1], prhs[2], NZ, NT, R, nke);
  else 
    fill_jacobian_real(prhs[0], prhs[1], prhs[2], NZ, NT, R, nke);
}




