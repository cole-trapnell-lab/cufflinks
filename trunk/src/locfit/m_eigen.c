/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include <math.h>
#include "mutil.h"
#define E_MAXIT 20
#define E_TOL 1.0e-8
#define SQR(x) ((x)*(x))

double e_tol(D,p)
double *D;
int p;
{ double mx;
  int i;
  if (E_TOL <= 0.0) return(0.0);
  mx = D[0];
  for (i=1; i<p; i++) if (D[i*(p+1)]>mx) mx = D[i*(p+1)];
  return(E_TOL*mx);
}

void eig_dec(X,P,d)
double *X, *P;
int d;
{ int i, j, k, iter, ms;
  double c, s, r, u, v;

  for (i=0; i<d; i++)
    for (j=0; j<d; j++) P[i*d+j] = (i==j);

  for (iter=0; iter<E_MAXIT; iter++)
  { ms = 0;
    for (i=0; i<d; i++)
      for (j=i+1; j<d; j++)
        if (SQR(X[i*d+j]) > 1.0e-15*fabs(X[i*d+i]*X[j*d+j]))
        { c = (X[j*d+j]-X[i*d+i])/2;
          s = -X[i*d+j];
          r = sqrt(c*c+s*s);
          c /= r;
          s = sqrt((1-c)/2)*(2*(s>0)-1);
          c = sqrt((1+c)/2);
          for (k=0; k<d; k++)
          { u = X[i*d+k]; v = X[j*d+k];
            X[i*d+k] = u*c+v*s;
            X[j*d+k] = v*c-u*s;
          }
          for (k=0; k<d; k++)
          { u = X[k*d+i]; v = X[k*d+j];
            X[k*d+i] = u*c+v*s;
            X[k*d+j] = v*c-u*s;
          }
          X[i*d+j] = X[j*d+i] = 0.0;
          for (k=0; k<d; k++)
          { u = P[k*d+i]; v = P[k*d+j];
            P[k*d+i] = u*c+v*s;
            P[k*d+j] = v*c-u*s;
          }
          ms = 1;
        }
    if (ms==0) return;
  }
  //printf("eig_dec not converged\n");
}

int eig_solve(J,x)
jacobian *J;
double *x;
{ int d, i, j, rank;
  double  *D, *P, *Q, *w;
  double tol;

  D = J->Z;
  P = Q = J->Q;
  d = J->p;
  w = J->wk;

  tol = e_tol(D,d);

  rank = 0;
  for (i=0; i<d; i++)
  { w[i] = 0.0;
    for (j=0; j<d; j++) w[i] += P[j*d+i]*x[j];
  }
  for (i=0; i<d; i++)
    if (D[i*d+i]>tol)
    { w[i] /= D[i*(d+1)];
      rank++;
    }
  for (i=0; i<d; i++)
  { x[i] = 0.0;
    for (j=0; j<d; j++) x[i] += Q[i*d+j]*w[j];
  }
  return(rank);
}

int eig_hsolve(J,v)
jacobian *J;
double *v;
{ int i, j, p, rank;
  double *D, *Q, *w;
  double tol;

    rank = 0;
  D = J->Z;
  Q = J->Q;
  p = J->p;
  w = J->wk;

  tol = e_tol(D,p);

  for (i=0; i<p; i++)
  { w[i] = 0.0;
    for (j=0; j<p; j++) w[i] += Q[j*p+i]*v[j];
  }
  for (i=0; i<p; i++)
    if (D[i*p+i]>tol)
    { w[i] /= sqrt(D[i*(p+1)]);
      rank++;
    }
  return(rank);
}

double eig_qf(J,v)
jacobian *J;
double *v;
{ int i, j, p;
  double sum, tol;

  p = J->p;
  sum = 0.0;
  tol = e_tol(J->Z,p);

  for (i=0; i<p; i++)
    if (J->Z[i*p+i]>tol)
    { J->wk[i] = 0.0;
      for (j=0; j<p; j++) J->wk[i] += J->Q[j*p+i]*v[j];
      sum += J->wk[i]*J->wk[i]/J->Z[i*p+i];
    }
  return(sum);
}
