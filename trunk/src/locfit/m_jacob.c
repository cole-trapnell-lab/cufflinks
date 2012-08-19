/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include <stdlib.h>
#include "math.h"
#include "stdio.h"
#include "mutil.h"

#define DEF_METH JAC_EIGD

int jac_reqd(int p) { return(2*p*(p+1)); }

double *jac_alloc(J,p,wk)
jacobian *J;
int p;
double *wk;
{ if (wk==NULL)
    wk = (double *)calloc(2*p*(p+1),sizeof(double));
  J->Z = wk; wk += p*p;
  J->Q = wk; wk += p*p;
  J->wk= wk; wk += p;
  J->dg= wk; wk += p;
  return(wk);
}

void jacob_dec(J, meth)
jacobian *J;
int meth;
{ int i, j, p;

  if (J->st != JAC_RAW) return;

  J->sm = J->st = meth;
  switch(meth)
  { case JAC_EIG:
      eig_dec(J->Z,J->Q,J->p);
      return;
    case JAC_EIGD:
      p = J->p;
      for (i=0; i<p; i++)
        J->dg[i] = (J->Z[i*(p+1)]<=0) ? 0.0 : 1/sqrt(J->Z[i*(p+1)]);
      for (i=0; i<p; i++)
        for (j=0; j<p; j++)
          J->Z[i*p+j] *= J->dg[i]*J->dg[j];
      eig_dec(J->Z,J->Q,J->p);
      J->st = JAC_EIGD;
      return;
    case JAC_CHOL:
      chol_dec(J->Z,J->p);
      return;
    default: printf("jacob_dec: unknown method %d",meth);
  }
}

int jacob_solve(J,v) /* (X^T W X)^{-1} v */
jacobian *J;
double *v;
{ int i, rank;

  if (J->st == JAC_RAW) jacob_dec(J,DEF_METH);

  switch(J->st)
  { case JAC_EIG:
      return(eig_solve(J,v));
    case JAC_EIGD:
      for (i=0; i<J->p; i++) v[i] *= J->dg[i];
      rank = eig_solve(J,v);
      for (i=0; i<J->p; i++) v[i] *= J->dg[i];
      return(rank);
    case JAC_CHOL:
      return(chol_solve(J->Z,v,J->p));
  }
  printf("jacob_solve: unknown method %d",J->st);
  return(0);
}

int jacob_hsolve(J,v) /*  J^{-1/2} v */
jacobian *J;
double *v;
{ int i;

  if (J->st == JAC_RAW) jacob_dec(J,DEF_METH);

  switch(J->st)
  { case JAC_EIG:
      return(eig_hsolve(J,v));
    case JAC_EIGD: /* eigenvalues on corr matrix */
      for (i=0; i<J->p; i++) v[i] *= J->dg[i];
      return(eig_hsolve(J,v));
    case JAC_CHOL:
      return(chol_hsolve(J->Z,v,J->p));
  }
  printf("jacob_hsolve: unknown method %d",J->st);
  return(0);
}

double jacob_qf(J,v)  /* vT J^{-1} v */
jacobian *J;
double *v;
{ int i;

  if (J->st == JAC_RAW) jacob_dec(J,DEF_METH);

  switch (J->st)
  { case JAC_EIG:
      return(eig_qf(J,v));
    case JAC_EIGD:
      for (i=0; i<J->p; i++) v[i] *= J->dg[i];
      return(eig_qf(J,v));
    case JAC_CHOL:
      return(chol_qf(J->Z,v,J->p));
    default:
      printf("jacob_qf: invalid method\n");
      return(0.0);
  }
}
