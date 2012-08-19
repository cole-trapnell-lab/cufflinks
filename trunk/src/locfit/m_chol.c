/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include <math.h>
#include "mutil.h"

void chol_dec(A,n)
double *A;
int n;
{ int i, j, k;
  for (j=0; j<n; j++)
  { k = n*j+j;
    for (i=0; i<j; i++) A[k] -= A[n*i+j]*A[n*i+j];
    if (A[k]<=0)
    { for (i=j; i<n; i++) A[n*j+i] = 0.0; }
    else
    { A[k] = sqrt(A[k]);
      for (i=j+1; i<n; i++)
      { for (k=0; k<j; k++)
          A[n*j+i] -= A[n*k+i]*A[n*k+j];
        A[n*j+i] /= A[n*j+j];
      }
    }
  }
  for (j=0; j<n; j++)
    for (i=j+1; i<n; i++) A[n*i+j] = 0.0;
}

int chol_solve(A,v,p)
double *A, *v;
int p;
{ int i, j;

  for (i=0; i<p; i++)
  { for (j=0; j<i; j++) v[i] -= A[j*p+i]*v[j];
    v[i] /= A[i*p+i];
  }
  for (i=p-1; i>=0; i--)
  { for (j=i+1; j<p; j++) v[i] -= A[i*p+j]*v[j];
    v[i] /= A[i*p+i];
  }
  return(p);
}

int chol_hsolve(A,v,p)
double *A, *v;
int p;
{ int i, j;

  for (i=0; i<p; i++)
  { for (j=0; j<i; j++) v[i] -= A[j*p+i]*v[j];
    v[i] /= A[i*p+i];
  }
  return(p);
}

double chol_qf(A,v,p)
double *A, *v;
int p;
{ int i, j;
  double sum;
 
  sum = 0.0;
  for (i=0; i<p; i++)
  { for (j=0; j<i; j++) v[i] -= A[j*p+i]*v[j];
    v[i] /= A[i*p+i];
    sum += v[i]*v[i];
  }
  return(sum);
}
