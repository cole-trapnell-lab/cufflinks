/*
 *   Copyright (c) 1996-2001 Jiayang Sun, Catherine Loader.
 *   This file is used by the simultaneous confidence band
 *   additions to Locfit.
 *
 */

#include "local.h"

extern INT cvi;

static double *fd, *ft, *lij, *d1a;
static INT par;

void assignk0(z,d,n) /* z should be n*(2*d*d+2*d+2); */
double *z;
INT d, n;
{ d1a= z; z += d*d*n;
  ft = z; z += n*(d*(d+1)+1);
  fd = z; z += n*(d+1);
}

void christ(d,nn,nl)  /* lij[i][j] = res proj. of Tij to (T1...Td) */
double nl;
INT d, nn;
{ INT i, j, k, l;
  double p4, *ll, v[1+MXDIM];
  for (i=0; i<d; i++)
    for (j=i; j<d; j++)
    { ll = &lij[(i*d+j)*nn];
      for (k=0; k<=d; k++) v[k] = innerprod(&ft[k*nn],ll,nn);
      bacT(fd,v,d+1,0,d+1);
      for (k=0; k<nn; k++)
        for (l=0; l<=d; l++)
          ll[k] -= ft[l*nn+k]*v[l];
      p4 = 0;
      for (k=0; k<=i+1; k++)
        p4 += fd[k*(d+1)+i+1]*fd[k*(d+1)+j+1];
      p4 = (fd[i+1]*fd[j+1]-p4)/(nl*nl);
      for (k=0; k<nn; k++)
        ll[k] = lij[(j*d+i)*nn+k] = ll[k] + p4*ft[k];
    }
}

void d1(n,d)   /* d1[i][j] = e_i^T (A^T A)^{-1} B_j^T */
INT n, d;
{ INT a, b, i, j;
  double *dd, v[MXDIM];
  for (i=0; i<d; i++)
  { for (j=0; j<d; j++) v[j] = 0;
    v[i] = 1;
    bacT(fd,v,d+1,1,d+1);
    for (j=0; j<d; j++)
    { dd = &d1a[(i*d+j)*n];
      for (a=0; a<n; a++)
      { dd[a] = 0;
        for (b=0; b<d; b++)
          dd[a] += v[b] * lij[(j*d+b)*n+a]; 
} } } }

void k2x(lf,des,kap)
lfit *lf;
design *des;
double *kap;
{ double det, s;
  INT i, j, k, l, d, m;
  d = lf->mi[MDIM];
  m = wdiag(lf,des,ft,1+(d>1),2,0);
  lij = &ft[(d+1)*m];
  for (i=0; i<m; i++)
    for (j=0; j<=d; j++)
      fd[i*(d+1)+j] = ft[j*m+i];
  QR1(fd,m,d+1,NULL);
  s = 0;
  if (d>1)
  { christ(d,m,fd[0]);
    d1(m,d);
    for (j=0; j<d; j++)
      for (k=0; k<j; k++)
        for (l=0; l<m; l++)
          s += d1a[(j*d+k)*m+l]*d1a[(k*d+j)*m+l]
             - d1a[(j*d+j)*m+l]*d1a[(k*d+k)*m+l];
  }
  det = 1;
  for (j=1; j<=d; j++)
    det *= fd[j*(d+2)]/fd[0];
  kap[0] = det;
  kap[2] = s*det*fd[0]*fd[0];
}

void l1x(lf,des,lap,re)
lfit *lf;
design *des;
double *lap;
INT re;
{ double det, t, sumcj, nu, *u, v[MXDIM];
  INT i, j, j1, d, m;
  d = lf->mi[MDIM]; u = des->res;
  m = wdiag(lf,des,ft,2,2,0);
  lij = &ft[(d+1)*m];
  for (i=0; i<m; i++)
  { t = ft[(re+1)*m+i];
    ft[(re+1)*m+i] = ft[d*m+i];
    ft[d*m+i] = t;
    for (j=0; j<d; j++) /* don't copy last column */
      fd[i*d+j] = ft[j*m+i];
    u[i] = ft[d*m+i];
  }
  QR1(fd,m,d,&ft[d*m]);
  bacK(fd,&ft[d*m],d);
  nu = 0;
  for (i=0; i<m; i++)
  { for (j=0; j<d; j++)
      u[i] -= ft[j*m+i]*ft[d*m+j];
    nu += u[i]*u[i];
  }    /* now u is outward vector, nu = ||u|| */
  sumcj = 0;
  for (i=0; i<d; i++) /* copy l(d-1,i) to l(re,i) */
    for (j=0; j<m; j++)
      lij[(re*d+i)*m+j] = lij[((d-1)*d+i)*m+j];
  for (j=0; j<d; j++)
  { if (j != re)
    { j1 = (j==(d-1)) ? re : j;
      for (i=0; i<d-1; i++)
        v[i] = innerprod(&lij[(i*d+j)*m],u,m);
      bacT(fd,v,d,1,d);
      sumcj += -v[j1];
    }
  }                                   /* stage 3,4 now complete */
  det = 1;
  for (j=1; j<d; j++)
    det *= fd[j*(d+1)]/fd[0];
  lap[0] = det;
  lap[1] = sumcj*det*fd[0]/sqrt(nu);
}

void m0x(lf,des,m0,re,rg)
lfit *lf;
design *des;
double *m0;
INT re, rg;
{ double det, t;
  INT d, m, i, j;
  d = lf->mi[MDIM];
  m = wdiag(lf,des,ft,1,2,0);
  for (i=0; i<m; i++)
  { t=ft[(rg+1)*m+i]; ft[(rg+1)*m+i]=ft[d*m+i]; ft[d*m+i]=t;
    t=ft[(re+1)*m+i]; ft[(re+1)*m+i]=ft[(d-1)*m+i]; ft[(d-1)*m+i]=t;
    for (j=0; j<=d; j++)
      fd[i*(d+1)+j] = ft[j*m+i];
  }
  det = 1;
  QR1(fd,m,d+1,NULL);
  for (j=1; j<d-1; j++)
    det *= fd[j*(d+2)]/fd[0];
  m0[0] = det*atan2(fd[d*(d+2)],-par*fd[d*(d+1)-1]);
}

INT constants(des,lf,kap)
design *des;
lfit *lf;
double *kap;
{ double h, k0[3], k1[3], l0[2], l1[2], m0[1], m1[1];
  double z[MXDIM], delt[MXDIM], mk, ml, mm;
  INT d, i, j, nnn, wt, index[MXDIM], *mi, pe, re, rg;
  cvi = -1; /* avoid cross valid */
  mi = lf->mi;
  d = mi[MDIM];
  if (lf_error) return(0);
  if ((lf->mi[MKER] != WPARM) && (lf->dp[DALP]>0))
    WARN(("constants are approximate for varying h"));
  mi[MP] = calcp(mi,mi[MDEG]);
  deschk(des,mi[MN],mi[MP]);
  preproc(des,lf,mi[MKER]!=WPARM);
  nnn = (ident==1) ? lf->mi[MP] : lf->mi[MN];
  lf->L = checkvarlen(lf->L,2*nnn*(d*d+d+1),"_hatmat",VDOUBLE);
  assignk0(vdptr(lf->L),d,nnn);
  mi[MDC] = 1;
  des->xev = z;

  mk = 1.0;
  for (i=0; i<d; i++)
  { index[i] = 0;
    z[i] = lf->fl[i];
    delt[i] = (lf->fl[i+d]-z[i])/(3*mi[MMINT]);
    mk *= delt[i];
  }
  i = 0;
  
  k0[0] = k0[1] = k0[2] = 0.0;
  l0[0] = l0[1] = 0.0;
  m0[0] = 0.0;

#ifdef CVERSION
  if (mi[MIT]==IMONT)
  { for (i=0; i<mi[MMINT]; i++)
    { for (j=0; j<d; j++) z[j] = lf->fl[j]+(lf->fl[j+d]-lf->fl[j])*runif();
      if ((mi[MKER]!=WPARM) | (!hasparcomp(lf)))
      { h = nbhd(lf,des,(INT)(mi[MN]*lf->dp[DALP]),lf->dp[DFXH],0);
        locfit(lf,des,h,1);
      }
      k2x(lf,des,k1);
      k0[0] += k1[0];
    }
    for (j=0; j<d; j++) k0[0] *= lf->fl[j+d]-lf->fl[j];
    kap[0] = k0[0]/mi[MMINT];
    return(1);
  }
#endif

  while(1)
  {
    wt = 1;
    for (i=0; i<d; i++)
      wt *= (4-2*(index[i]%2==0)-(index[i]==0)-(index[i]==mi[MMINT]));
    if ((mi[MKER]!=WPARM) | (!hasparcomp(lf)))
    { h = nbhd(lf,des,(INT)(mi[MN]*lf->dp[DALP]),lf->dp[DFXH],0);
      locfit(lf,des,h,1);
    }
    k2x(lf,des,k1);
    k0[0] += wt*mk*k1[0];
    k0[2] += wt*mk*k1[2];

    for (re=0; re<d; re++) if ((index[re]==0) | (index[re]==mi[MMINT]))
    { l1x(lf,des,l1,re);
      ml = 1;
      for (i=0; i<d; i++) if (i!=re) ml *= delt[i];
      pe = 1-2*(index[re]==0);
      l0[0] += wt*ml*l1[0];
      l0[1] += wt*ml*pe*l1[1];

      for (rg=re+1; rg<d; rg++) if ((index[rg]==0) | (index[rg]==mi[MMINT]))
      { par = pe*(1-2*(index[rg]==0));
        m0x(lf,des,m1,re,rg);
        mm = 1;
        for (i=0; i<d; i++) if ((i!=re) & (i!=rg)) mm *= delt[i];
        m0[0] += wt*mm*m1[0];
      }
    }

    /* compute next grid point */
    for (i=0; i<d; i++)
    { index[i]++;
      z[i] = lf->fl[i]+3*delt[i]*index[i];
      if (index[i]>mi[MMINT])
      { index[i] = 0;
        z[i] = lf->fl[i];
        if (i==d-1) /* done */
        { kap[0] = k0[0];
          kap[1] = l0[0]/2;
          if (d==1) return(2);
          k0[2] = -k0[2] - d*(d-1)*k0[0]/2;
          if (mi[MDEB]>0)
          { printf("constants:\n");
            printf("  k0: %8.5f  k2: %8.5f\n",k0[0],k0[2]);
            printf("  l0: %8.5f  l1: %8.5f\n",l0[0],l1[1]);
            printf("  m0: %8.5f\n",m0[0]);
            printf("  check: %8.5f\n",(k0[0]+k0[2]+l0[1]+m0[0])/(2*PI));
          }
          kap[2] = (k0[2]+l0[1]+m0[0])/(2*PI);
          return(3);
        }
      }
      else i = d;
    }

  }
}

double tailp(c,k0,m,d,nu)
double c, *k0, nu;
INT m, d;
{ INT i;
  double p;
  p = 0;
  if (nu==0)
  { for (i=0; i<m; i++) if (k0[i]>0)
      p += k0[i]*exp(LGAMMA((d+1-i)/2.0)-(d+1-i)*LOGPI/2)
          *(1-pchisq(c*c,(double) d+1-i));
  }
  else
  { for (i=0; i<m; i++) if (k0[i]>0)
      p += k0[i]*exp(LGAMMA((d+1-i)/2.0)-(d+1-i)*LOGPI/2)
          *(1-pf(c*c/(d+1-i),(double) (d+1-i), nu));
  }
  return(p);
}

double taild(c,k0,m,d,nu)
double c, *k0, nu;
INT m, d;
{ double p;
  INT i;
  p = 0;
  if (nu==0)
  { for (i=0; i<m; i++) if (k0[i]>0)
      p += k0[i]*exp(LGAMMA((d+1-i)/2.0)-(d+1-i)*LOGPI/2)
          *2*c*dchisq(c*c,(double) (d+1-i));
  }
  else
  { for (i=0; i<m; i++) if (k0[i]>0)
      p += k0[i]*exp(LGAMMA((d+1-i)/2.0)-(d+1-i)*LOGPI/2)
          *2*c*df(c*c/(d+1-i),(double) (d+1-i), nu)/(d+1-i);
  }
  return(-p);
}

double critval(k0,m,d,al,it,s,nu)
double *k0, al, nu;
INT m, d, it, s;
{ double c, cn, c0, c1, tp, td;
  INT j;
  if (m<0) ERROR(("critval: no terms?"));
  if (m>d+1) m = d+1;
  if ((al<=0) | (al>=1)) ERROR(("critval: invalid alpha %8.5f",al));
  if (lf_error) return(0.0);
  if (al>0.5) WARN(("critval: A mighty large tail probability al=%8.5f",al));
  if (s==1) al = 2*al;
  if (m==0) { d = 0; k0[0] = 1; m = 1; }
  c = 2.0; c0 = 0.0; c1 = 0.0;
  for (j=0; j<it; j++)
  { tp = tailp(c,k0,m,d,nu)-al;
    td = taild(c,k0,m,d,nu);
    if (tp>0) c0 = c;
    if (tp<0) c1 = c;
    cn = c - tp/td;
    if (cn<c0) cn = (c+c0)/2;
    if ((c1>0.0) && (cn>c1)) cn = (c+c1)/2;
    c = cn;
    if (fabs(tp/al)<1.0e-10) return(c);
  }
  return(c);
}

#ifdef SVERSION
void scritval(k0,d,cov,m,rdf,z)
double *k0, *z, *cov, *rdf;
INT *d, *m;
{ lf_error = 0;
  *z = critval(k0,*m,*d,1-*cov,10,2,*rdf);
}
#endif
