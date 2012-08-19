/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

void svd(x,p,q,d,mxit)  /* svd of square matrix */
double *x, *p, *q;
INT d, mxit;
{ INT i, j, k, iter, ms, zer;
  double r, u, v, cp, cm, sp, sm, c1, c2, s1, s2, mx;
  for (i=0; i<d; i++)
    for (j=0; j<d; j++) p[i*d+j] = q[i*d+j] = (i==j);
  for (iter=0; iter<mxit; iter++)
  { ms = 0;
    for (i=0; i<d; i++)
      for (j=i+1; j<d; j++)
      { mx = MAX(fabs(x[i*d+j]),fabs(x[j*d+i]));
        zer = 1;
        if (mx*mx>1.0e-15*fabs(x[i*d+i]*x[j*d+j]))
        { if (fabs(x[i*(d+1)])<fabs(x[j*(d+1)]))
          { for (k=0; k<d; k++)
            { u = x[i*d+k]; x[i*d+k] = x[j*d+k]; x[j*d+k] = u;
              u = p[k*d+i]; p[k*d+i] = p[k*d+j]; p[k*d+j] = u;
            }
            for (k=0; k<d; k++)
            { u = x[k*d+i]; x[k*d+i] = x[k*d+j]; x[k*d+j] = u;
              u = q[k*d+i]; q[k*d+i] = q[k*d+j]; q[k*d+j] = u;
            }
          }
          cp = x[i*(d+1)]+x[j*(d+1)];
          sp = x[j*d+i]-x[i*d+j];
          r = sqrt(cp*cp+sp*sp);
          if (r>0) { cp /= r; sp /= r; }
              else { cp = 1.0; zer = 0;}
          cm = x[i*(d+1)]-x[j*(d+1)];
          sm = x[i*d+j]+x[j*d+i];
          r = sqrt(cm*cm+sm*sm);
          if (r>0) { cm /= r; sm /= r; }
              else { cm = 1.0; zer = 0;}
          c1 = cm+cp;
          s1 = sm+sp;
          r = sqrt(c1*c1+s1*s1);
          if (r>0) { c1 /= r; s1 /= r; }
              else { c1 = 1.0; zer = 0;}
          if (fabs(s1)>ms) ms = fabs(s1);
          c2 = cm+cp;
          s2 = sp-sm;
          r = sqrt(c2*c2+s2*s2);
          if (r>0) { c2 /= r; s2 /= r; }
              else { c2 = 1.0; zer = 0;}
          for (k=0; k<d; k++)
          { u = x[i*d+k]; v = x[j*d+k];
            x[i*d+k] = c1*u+s1*v;
            x[j*d+k] = c1*v-s1*u;
            u = p[k*d+i]; v = p[k*d+j];
            p[k*d+i] = c1*u+s1*v;
            p[k*d+j] = c1*v-s1*u;
          }
          for (k=0; k<d; k++)
          { u = x[k*d+i]; v = x[k*d+j];
            x[k*d+i] = c2*u-s2*v;
            x[k*d+j] = s2*u+c2*v;
            u = q[k*d+i]; v = q[k*d+j];
            q[k*d+i] = c2*u-s2*v;
            q[k*d+j] = s2*u+c2*v;
          }
          if (zer) x[i*d+j] = x[j*d+i] = 0.0;
          ms = 1;
        }
      }
    if (ms==0) iter=mxit+10;
  }
  if (iter==mxit) WARN(("svd not converged"));
  for (i=0; i<d; i++)
    if (x[i*d+i]<0)
    { x[i*d+i] = -x[i*d+i];
      for (j=0; j<d; j++) p[j*d+i] = -p[j*d+i];
    }
}

INT svdsolve(x,w,P,D,Q,d,tol) /* original X = PDQ^T; comp. QD^{-1}P^T x */
double *x, *w, *P, *D, *Q, tol;
INT d;
{ INT i, j, rank;
  double mx;
  if (tol>0)
  { mx = D[0];
    for (i=1; i<d; i++) if (D[i*(d+1)]>mx) mx = D[i*(d+1)];
    tol *= mx;
  }
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

void hsvdsolve(x,w,P,D,Q,d,tol) /* original X = PDQ^T; comp. D^{-1/2}P^T x */
double *x, *w, *P, *D, *Q, tol;
INT d;
{ INT i, j;
  double mx;
  if (tol>0)
  { mx = D[0];
    for (i=1; i<d; i++) if (D[i*(d+1)]>mx) mx = D[i*(d+1)];
    tol *= mx;
  }
  for (i=0; i<d; i++)
  { w[i] = 0.0;
    for (j=0; j<d; j++) w[i] += P[j*d+i]*x[j];
  }
  for (i=0; i<d; i++) if (D[i*d+i]>tol) w[i] /= sqrt(D[i*(d+1)]);
  for (i=0; i<d; i++) x[i] = w[i];
}

void QRupd(X,u,p,w,y)
double *X, *u, *w, y;
INT p;
{ INT i, j;
  double s, c, r, t;
  for (i=0; i<p; i++)
  { c = X[i*p+i]; s = u[i];
    if (s!=0)
    { r = sqrt(c*c+s*s);
      c = c/r; s = s/r;
      for (j=i; j<p; j++)
      { t = c*X[i*p+j]+s*u[j];
        u[j] = c*u[j]-s*X[i*p+j];
        X[i*p+j] = t;
      }
      if (w!=NULL)
      { t = c*w[i]+s*y;
        y = c*y-s*w[i];
        w[i] = t;
      }
    }
  }
}

void QR1(X,n,p,w)   /* QR of X (n*p); take w for the ride, if not NULL */
double *X, *w;
INT n, p;
{ INT i, j, k, mi;
  double c, s, mx, nx, t;
  for (j=0; j<p; j++)
  { mi = j; mx = fabs(X[(p+1)*j]); nx = mx*mx;
    for (i=j+1; i<n; i++)
    { nx += X[p*i+j]*X[p*i+j];
      if (fabs(X[p*i+j])>mx) { mi = i; mx = fabs(X[p*i+j]); }
    }
    for (i=0; i<p; i++)
    { t = X[p*j+i]; X[p*j+i] = X[p*mi+i]; X[p*mi+i] = t; }
    if (w != NULL)
    { t = w[j]; w[j] = w[mi]; w[mi] = t; }
    if (X[(p+1)*j]>0)
    { for (i=j; i<p; i++) X[p*j+i] = -X[p*j+i];
      if (w != NULL) w[j] = -w[j];
    }
    nx = sqrt(nx);
    c = nx*(nx-X[(p+1)*j]);
    if (c!=0)
    { for (i=j+1; i<p; i++)
      { s = 0;
        for (k=j; k<n; k++)
          s += X[p*k+i]*X[p*k+j];
        s = (s-nx*X[p*j+i])/c;
        for (k=j; k<n; k++)
          X[p*k+i] -= s*X[p*k+j];
        X[p*j+i] += s*nx;
      }
      if (w != NULL)
      { s = 0;
        for (k=j; k<n; k++)
          s += w[k]*X[p*k+j];
        s = (s-nx*w[j])/c;
        for (k=j; k<n; k++)
          w[k] -= s*X[p*k+j];
        w[j] += s*nx;
      }
      X[p*j+j] = nx;
   }
} }

void bacK(R,f,p)   /* R^{-1} f */
double *R, *f;
INT p;
{ INT i, j;
  for (i=p-1; i>=0; i--)
  { for (j=i+1; j<p; j++)
      f[i] -= R[i*p+j]*f[j];
    f[i] /= R[i*(p+1)];
  }
}

void bacT(R,f,p,i0,i1)   /* R^{-1} (R^T)^{-1} f; p columns; use i0-(i1-1) */
double *R, *f;
INT p;
int i0;
int i1;
{ INT i, j;
  for (i=i0; i<i1; i++)
  { for (j=i0; j<i; j++)
      f[i-i0] -= R[j*p+i]*f[j-i0];
    f[i-i0] /= R[i*(p+1)];
  }
  for (i=i1-1; i>=i0; i--)
  { for (j=i+1; j<i1; j++)
      f[i-i0] -= R[i*p+j]*f[j-i0];
    f[i-i0] /= R[i*(p+1)];
  }
}

void solve(A,b,d) /* this is crude! A organized by column. */
double *A, *b;
INT d;
{ INT i, j, k;
  double piv;
  for (i=0; i<d; i++)
  { piv = A[(d+1)*i];
    for (j=i; j<d; j++) A[j*d+i] /= piv;
    b[i] /= piv;
    for (j=0; j<d; j++) if (j != i)
    { piv = A[i*d+j];
      A[i*d+j] = 0;
      for (k=i+1; k<d; k++)
        A[k*d+j] -= piv*A[k*d+i];
      b[j] -= piv*b[i];
    }
  }
}

void grsc(X,d)
double *X;
INT d;
{ INT i, j, k;
  double nx;
  for (i=0; i<d; i++)
  { nx = 0;
    for (j=0; j<d; j++) nx += X[i*d+j]*X[i*d+j];
    nx = sqrt(nx);
    for (j=0; j<d; j++) X[i*d+j] /= nx;
    for (j=i+1; j<d; j++)
    { nx = 0;
      for (k=0; k<d; k++) nx += X[i*d+k]*X[j*d+k];
      for (k=0; k<d; k++) X[j*d+k] -= nx*X[i*d+k];
    }
  }
}

void setzero(v,p)
double *v;
INT p;
{ int i;
  for (i=0; i<p; i++) v[i] = 0.0;
}

void unitvec(x,k,p)
double *x;
INT k, p;
{ setzero(x,p);
  x[k] = 1.0;
}

double innerprod(v1,v2,p)
double *v1, *v2;
INT p;
{ int i;
  double s;
  s = 0;
  for (i=0; i<p; i++) s += v1[i]*v2[i];
  return(s);
}

void addouter(re,v1,v2,p,c)
double *re, *v1, *v2, c;
INT p;
{ int i, j;
  for (i=0; i<p; i++)
    for (j=0; j<p; j++)
      re[i*p+j] += c*v1[i]*v2[j];
}

void multmatscal(A,z,n)
double *A, z;
INT n;
{ INT i;
  for (i=0; i<n; i++) A[i] *= z;
}

/*
 *  transpose() transposes an m*n matrix in place.
 *  At input, the matrix has n rows, m columns and
 *    x[0..n-1] is the is the first column.
 *  At output, the matrix has m rows, n columns and
 *    x[0..m-1] is the first column.
 */
void transpose(x,m,n)
double *x;
INT m, n;
{ INT t0, t, ti, tj;
  double z;
  for (t0=1; t0<m*n-2; t0++)
    { ti = t0%m; tj = t0/m;
      do
      { t = ti*n+tj;
        ti= t%m;
        tj= t/m;
      } while (t<t0);
      z = x[t];
      x[t] = x[t0];
      x[t0] = z;
    }
}

/* trace of an n*n square matrix. */
double m_trace(x,n)
double *x;
int n;
{ int i;
  double sum;
  sum = 0;
  for (i=0; i<n; i++)
    sum += x[i*(n+1)];
  return(sum);
}
