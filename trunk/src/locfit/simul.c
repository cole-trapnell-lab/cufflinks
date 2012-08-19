/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

static double pen, sig2;

void goldensec(f,des,tr,eps,xm,ym,meth)
double (*f)(), eps, *xm, *ym;
INT meth;
design *des;
lfit *tr;
{ double x[4], y[4], xx[11], yy[11];
  INT i, im;
  xx[0] = tr->dp[DFXH];
  if (xx[0]<=0)
  { ERROR(("regband: initialize h>0"));
    return;
  }
  for (i=0; i<=10; i++)
  { if (i>0) xx[i] = (1+GOLDEN)*xx[i-1];
    yy[i] = f(xx[i],des,tr,meth);
    if ((i==0) || (yy[i]<yy[im])) im = i;
  }
  if (im==0) im = 1;
  if (im==10)im = 9;
  x[0] = xx[im-1]; y[0] = yy[im-1];
  x[1] = xx[im];   y[1] = yy[im];
  x[3] = xx[im+1]; y[3] = yy[im+1];
  x[2] = GOLDEN*x[3]+(1-GOLDEN)*x[0];
  y[2] = f(x[2],des,tr,meth);
  while (x[3]-x[0]>eps)
  { if (y[1]<y[2])
    { x[3] = x[2]; y[3] = y[2];
      x[2] = x[1]; y[2] = y[1];
      x[1] = GOLDEN*x[0]+(1-GOLDEN)*x[3];
      y[1] = f(x[1],des,tr,meth);
    }
    else
    { x[0] = x[1]; y[0] = y[1];
      x[1] = x[2]; y[1] = y[2];
      x[2] = GOLDEN*x[3]+(1-GOLDEN)*x[0];
      y[2] = f(x[2],des,tr,meth);
    }
  }
  im = 0;
  for (i=1; i<4; i++) if (y[i]<y[im]) im = i;
  *xm = x[im]; *ym = y[im];
}

double dnk(x,k)
double x;
INT k;
{ double f;
  switch(k)
  { case 0: f = 1; break;
    case 1: f = -x; break;
    case 2: f = x*x-1; break;
    case 3: f = x*(x*x-3); break;
    case 4: f = 3-x*x*(6-x*x); break;
    case 5: f = -x*(15-x*x*(10-x*x)); break;
    case 6: f = -15+x*x*(45-x*x*(15-x*x)); break;
    default: ERROR(("dnk: k=%d too large",k)); return(0.0);
  }
  return(f*exp(-x*x/2)/S2PI);
}

double locai(h,des,tr)
double h;
design *des;
lfit *tr;
{ double cp;
  tr->dp[DALP] = h;
  startlf(des,tr,procv,0);
  ressumm(tr,des);
  cp = -2*tr->dp[DLK]+pen*tr->dp[DT0];
  return(cp);
}

double loccp(h,des,tr,m) /* m=1: cp    m=2: gcv */
double h;
design *des;
lfit *tr;
int m;
{ double cp;
  INT dg;
  tr->dp[DALP] = 0;
  tr->dp[DFXH] = h;
  dg = tr->mi[MDEG]; tr->mi[MDEG] = tr->mi[MDEG0];
  startlf(des,tr,procv,0);
  ressumm(tr,des);
  if (m==1)
    cp = -2*tr->dp[DLK]/sig2 - tr->mi[MN] + 2*tr->dp[DT0];
  else cp = -2*tr->mi[MN]*tr->dp[DLK]/((tr->mi[MN]-tr->dp[DT0])*(tr->mi[MN]-tr->dp[DT0]));
  printf("h %8.5f  deg %2d  rss %8.5f  trl %8.5f  cp: %8.5f\n",h,tr->mi[MDEG],-2*tr->dp[DLK],tr->dp[DT0],cp);
  tr->mi[MDEG0] = tr->mi[MDEG]; tr->mi[MDEG] = dg;
  return(cp);
}

double cp(des,tr,meth)
design *des;
lfit *tr;
INT meth;
{ double hm, ym;
  goldensec(loccp,des,tr,0.001,&hm,&ym,meth);
  return(hm);
}

double gkk(des,tr)
design *des;
lfit *tr;
{ double h, h5, nf, th;
  INT i, j, n, dg0, dg1;
  tr->mi[MEV]=EDATA;
  tr->dp[DALP] = 0;
  n = tr->mi[MN];
  dg0 = tr->mi[MDEG0];     /* target degree */
  dg1 = dg0+1+(dg0%2==0);  /* pilot degree */
  nf = exp(log(1.0*n)/10); /* bandwidth inflation factor */
  h = tr->dp[DFXH];        /* start bandwidth */
  for (i=0; i<=10; i++)
  { tr->mi[MDEG] = dg1;
    tr->dp[DFXH] = h*nf;
    startlf(des,tr,procv,0);
    th = 0;
    for (j=10; j<n-10; j++)
      th += tr->coef[dg1*n+j]*tr->coef[dg1*n+j];
th *= n/(n-20.0);
    h5 = sig2*Wikk(tr->mi[MKER],dg0)/th;
    h = exp(log(h5)/(2*dg1+1));
/* printf("pilot %8.5f  sel %8.5f\n",tr->dp[DFXH],h); */
  }
  return(h);
}

double rsw(des,tr,kk)
design *des;
lfit *tr;
INT *kk;
{ INT i, j, k, nmax, nvm, n, mk, ev, dg0, dg1;
  double rss[6], cp[6], th22, dx, d2, hh;
  nmax = 5;
  ev = tr->mi[MEV];  tr->mi[MEV] = EGRID;
  mk = tr->mi[MKER]; tr->mi[MKER]= WRECT;
  dg0 = tr->mi[MDEG0];
  dg1 = 1 + dg0 + (dg0%2==0);
  tr->mi[MDEG]= 4;
  for (k=nmax; k>0; k--)
  { tr->mg[0] = k;
    tr->fl[0] = 1.0/(2*k); tr->fl[1] = 1-1.0/(2*k);
    tr->dp[DALP] = 0; tr->dp[DFXH] = 1.0/(2*k);
    startlf(des,tr,procv,0);
    nvm = tr->nvm;
    rss[k] = 0;
    for (i=0; i<k; i++) rss[k] += -2*tr->lik[i];
  }
  n = tr->mi[MN]; k = 1;
  for (i=1; i<=nmax; i++)
  { /* cp[i] = (n-5*nmax)*rss[i]/rss[nmax]-(n-10*i); */
    cp[i] = rss[i]/sig2-(n-10*i);
    if (cp[i]<cp[k]) k = i;
  }
  *kk = k;
  tr->mg[0] = k;
  tr->fl[0] = 1.0/(2*k); tr->fl[1] = 1-1.0/(2*k);
  tr->dp[DALP] = 0; tr->dp[DFXH] = 1.0/(2*k);
  startlf(des,tr,procv,0);
  tr->mi[MKER] = mk; tr->mi[MEV] = ev;
  nvm = tr->nvm;
  th22 = 0;
  for (i=10; i<n-10; i++)
  { j = floor(k*datum(tr,0,i));
    if (j>=k) j = k-1;
    dx = datum(tr,0,i)-evptx(tr,0,j);
    if (dg1==2)
      d2 = tr->coef[2*nvm+j]+dx*tr->coef[3*nvm+j]+dx*dx*tr->coef[4*nvm+j]/2;
    else d2 = tr->coef[4*nvm+j];
    th22 += d2*d2;
  }
  hh = Wikk(mk,dg0)*sig2/th22*(n-20.0)/n;
  return(exp(log(hh)/(2*dg1+1)));
}

void rband(des,tr,hhat,meth,nmeth,kk)
design *des;
lfit *tr;
double *hhat;
INT *meth, *nmeth, *kk;
{ INT i, deg;
  double h0;

  /* first, estimate sigma^2 */
  deg = tr->mi[MDEG]; tr->mi[MDEG] = 2;
  h0 = tr->dp[DFXH];  tr->dp[DFXH] = 0.05;
printf("alp: %8.5f  h: %8.5f  deg %2d  ev %2d\n",tr->dp[DALP],tr->dp[DFXH],tr->mi[MDEG],tr->mi[MEV]);
  startlf(des,tr,procv,0);
  ressumm(tr,des);
  tr->mi[MDEG] = deg; tr->dp[DFXH] = h0;
  sig2 = tr->dp[DRV]; 
  printf("sd est: %8.5f\n",sqrt(tr->dp[DRV]));

  for (i=0; i<*nmeth; i++)
  { switch(meth[i])
    { case 1: hhat[i] = cp(des,tr,1);
              break;
      case 2: hhat[i] = cp(des,tr,2);
              break;
      case 3: hhat[i] = gkk(des,tr);
              break;
      case 4: hhat[i] = rsw(des,tr,kk);
              break;
      default: hhat[i] = 0;
    }
    tr->dp[DFXH] = h0;
    tr->mi[MDEG] = deg;
  }
}
