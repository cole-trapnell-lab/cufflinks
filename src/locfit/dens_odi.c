/*
 *   Copyright (c) 1996-200 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *  Routines for one-dimensional numerical integration
 *  in density estimation. The entry point is
 *
 *  onedint(cf,mi,l0,l1,resp)
 *
 *  which evaluates int W(u)u^j exp( P(u) ), j=0..2*deg.
 *  P(u) = cf[0] + cf[1]u + cf[2]u^2/2 + ... + cf[deg]u^deg/deg!
 *  l0 and l1 are the integration limits.
 *  The results are returned through the vector resp.
 *
 */

#include "local.h"

static int debug;

INT exbctay(b,c,n,z) /* n-term taylor series of e^(bx+cx^2) */
double b, c, *z;
INT n;
{ double ec[20];
  INT i, j;
  z[0] = 1;
  for (i=1; i<=n; i++) z[i] = z[i-1]*b/i;
  if (c==0.0) return(n);
  if (n>=40)
  { WARN(("exbctay limit to n<40"));
    n = 39;
  }
  ec[0] = 1;
  for (i=1; 2*i<=n; i++) ec[i] = ec[i-1]*c/i;
  for (i=n; i>1; i--)
    for (j=1; 2*j<=i; j++)
      z[i] += ec[j]*z[i-2*j];
  return(n);
}

double explinjtay(l0,l1,j,cf)
/* int_l0^l1 x^j e^(a+bx+cx^2); exbctay aroud l1 */
double l0, l1, *cf;
INT j;
{ double tc[40], f, s;
  INT k, n;
  if ((l0!=0.0) | (l1!=1.0)) WARN(("explinjtay: invalid l0, l1"));
  n = exbctay(cf[1]+2*cf[2]*l1,cf[2],20,tc);
  s = tc[0]/(j+1);
  f = 1/(j+1);
  for (k=1; k<=n; k++)
  { f *= -k/(j+k+1.0);
    s += tc[k]*f;
  }
  return(f);
}

void explint1(l0,l1,cf,I,p) /* int x^j exp(a+bx); j=0..p-1 */
double l0, l1, *cf, *I;
INT p;
{ double y0, y1, f;
  INT j, k, k1;
  y0 = lf_exp(cf[0]+l0*cf[1]);
  y1 = lf_exp(cf[0]+l1*cf[1]);
  if (p<2*fabs(cf[1])) k = p; else k = (INT)fabs(cf[1]);

  if (k>0)
  { I[0] = (y1-y0)/cf[1];
    for (j=1; j<k; j++) /* forward steps for small j */
    { y1 *= l1; y0 *= l0;
      I[j] = (y1-y0-j*I[j-1])/cf[1];
    }
    if (k==p) return;
    y1 *= l1; y0 *= l0;
  }

  f = 1; k1 = k;
  while ((k<50) && (f>1.0e-8)) /* initially Ik = diff(x^{k+1}e^{a+bx}) */
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
    if (k>=p) f *= fabs(cf[1])/(k+1);
    k++;
  }
  if (k==50) WARN(("explint1: want k>50"));
  I[k] = 0.0;
  for (j=k-1; j>=k1; j--) /* now do back step recursion */
    I[j] = (I[j]-cf[1]*I[j+1])/(j+1);
}

void explintyl(l0,l1,cf,I,p) /* small c, use taylor series and explint1 */
double l0, l1, *cf, *I;
INT p;
{ INT i;
  double c;
  explint1(l0,l1,cf,I,p+8);
  c = cf[2];
  for (i=0; i<p; i++)
    I[i] = (((I[i+8]*c/4+I[i+6])*c/3+I[i+4])*c/2+I[i+2])*c+I[i];
}

void solvetrid(X,y,m)
double *X, *y;
INT m;
{ INT i;
  double s;
  for (i=1; i<m; i++)
  { s = X[3*i]/X[3*i-2];
    X[3*i] = 0; X[3*i+1] -= s*X[3*i-1];
    y[i] -= s*y[i-1];
  }
  for (i=m-2; i>=0; i--)
  { s = X[3*i+2]/X[3*i+4];
    X[3*i+2] = 0;
    y[i] -= s*y[i+1];
  }
  for (i=0; i<m; i++) y[i] /= X[3*i+1];
}

void initi0i1(I,cf,y0,y1,l0,l1)
double *I, *cf, y0, y1, l0, l1;
{ double a0, a1, c, d, bi;
  d = -cf[1]/(2*cf[2]); c = sqrt(2*fabs(cf[2]));
  a0 = c*(l0-d); a1 = c*(l1-d);
  if (cf[2]<0)
  { bi = lf_exp(cf[0]+cf[1]*d+cf[2]*d*d)/c;
    if (a0>0)
    { if (a0>6) I[0] = (y0*ptail(-a0)-y1*ptail(-a1))/c;
      else I[0] = S2PI*(pnorm(-a0,0.0,1.0)-pnorm(-a1,0.0,1.0))*bi;
    }
    else
    { if (a1< -6) I[0] = (y1*ptail(a1)-y0*ptail(a0))/c;
      else I[0] = S2PI*(pnorm(a1,0.0,1.0)-pnorm(a0,0.0,1.0))*bi;
    }
  }
  else
    I[0] = (y1*daws(a1)-y0*daws(a0))/c;
  I[1] = (y1-y0)/(2*cf[2])+d*I[0];
}

void explinsid(l0,l1,cf,I,p) /* large b; don't use fwd recursion */
double l0, l1, *cf, *I;
INT p;
{ INT k, k0, k1, k2;
  double y0, y1, Z[150];
if (debug) printf("side: %8.5f %8.5f %8.5f    limt %8.5f %8.5f  p %2d\n",cf[0],cf[1],cf[2],l0,l1,p);
 
  k0 = 2;
  k1 = (INT)(fabs(cf[1])+fabs(2*cf[2]));
  if (k1<2) k1 = 2;
  if (k1>p+20) k1 = p+20;
  k2 = p+20;

  if (debug) printf("k0 %2d  k1 %2d  k2 %2d  p %2d\n",k0,k1,k2,p);

  y0 = lf_exp(cf[0]+l0*(cf[1]+l0*cf[2]));
  y1 = lf_exp(cf[0]+l1*(cf[1]+l1*cf[2]));
  initi0i1(I,cf,y0,y1,l0,l1);
if (debug) printf("i0 %8.5f  i1 %8.5f\n",I[0],I[1]);

  y1 *= l1; y0 *= l0; /* should be x^(k1)*exp(..) */
  if (k0<k1) /* center steps; initially x^k*exp(...) */
    for (k=k0; k<k1; k++)
    { y1 *= l1; y0 *= l0;
      I[k] = y1-y0;
      Z[3*k] = k; Z[3*k+1] = cf[1]; Z[3*k+2] = 2*cf[2];
    }
   
  y1 *= l1; y0 *= l0; /* should be x^(k1)*exp(..) */
if (debug) printf("k1 %2d  y0 %8.5f  y1 %8.5f\n",k1,y0,y1);
  for (k=k1; k<k2; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }
  I[k2] = I[k2+1] = 0.0;
  for (k=k2-1; k>=k1; k--)
    I[k] = (I[k]-cf[1]*I[k+1]-2*cf[2]*I[k+2])/(k+1);

  if (k0<k1)
  { I[k0] -= k0*I[k0-1];
    I[k1-1] -= 2*cf[2]*I[k1];
    Z[3*k0] = Z[3*k1-1] = 0;
    solvetrid(&Z[3*k0],&I[k0],k1-k0);
  }
if (debug)
{ printf("explinsid:\n");
  for (k=0; k<p; k++) printf("  %8.5f\n",I[k]);
}
}

void explinbkr(l0,l1,cf,I,p) /* small b,c; use back recursion */
double l0, l1, *cf, *I;
INT p;
{ INT k, km;
  double y0, y1;
  y0 = lf_exp(cf[0]+l0*(cf[1]+cf[2]*l0));
  y1 = lf_exp(cf[0]+l1*(cf[1]+cf[2]*l1));
  km = p+10;
  for (k=0; k<=km; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }
  I[km+1] = I[km+2] = 0;
  for (k=km; k>=0; k--)
    I[k] = (I[k]-cf[1]*I[k+1]-2*cf[2]*I[k+2])/(k+1);
}

void explinfbk0(l0,l1,cf,I,p) /* fwd and bac recur; b=0; c<0 */
double l0, l1, *cf, *I;
INT p;
{ double y0, y1, f1, f2, f, ml2;
  INT k, ks;

  y0 = lf_exp(cf[0]+l0*l0*cf[2]);
  y1 = lf_exp(cf[0]+l1*l1*cf[2]);
  initi0i1(I,cf,y0,y1,l0,l1);

  ml2 = MAX(l0*l0,l1*l1);
  ks = 1+(INT)(2*fabs(cf[2])*ml2);
  if (ks<2) ks = 2;
  if (ks>p-3) ks = p;

  /* forward recursion for k < ks */
  for (k=2; k<ks; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = (y1-y0-(k-1)*I[k-2])/(2*cf[2]);
  }
  if (ks==p) return;

  y1 *= l1*l1; y0 *= l0*l0;
  for (k=ks; k<p; k++) /* set I[k] = x^{k+1}e^(a+cx^2) | {l0,l1} */
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }

  /* initialize I[p-2] and I[p-1] */
  f1 = 1.0/p; f2 = 1.0/(p-1);
  I[p-1] *= f1; I[p-2] *= f2;
  k = p; f = 1.0;
  while (f>1.0e-8)
  { y1 *= l1; y0 *= l0;
    if ((k-p)%2==0) /* add to I[p-2] */
    { f2 *= -2*cf[2]/(k+1);
      I[p-2] += (y1-y0)*f2;
    }
    else /* add to I[p-1] */
    { f1 *= -2*cf[2]/(k+1);
      I[p-1] += (y1-y0)*f1;
      f *= 2*fabs(cf[2])*ml2/(k+1);
    }
    k++;
  }
  
  /* use back recursion for I[ks..(p-3)] */
  for (k=p-3; k>=ks; k--)
    I[k] = (I[k]-2*cf[2]*I[k+2])/(k+1);
}

void explinfbk(l0,l1,cf,I,p) /* fwd and bac recur; b not too large */
double l0, l1, *cf, *I;
INT p;
{ double y0, y1;
  INT k, ks, km;

  y0 = lf_exp(cf[0]+l0*(cf[1]+l0*cf[2]));
  y1 = lf_exp(cf[0]+l1*(cf[1]+l1*cf[2]));
  initi0i1(I,cf,y0,y1,l0,l1);

  ks = (INT)(3*fabs(cf[2]));
  if (ks<3) ks = 3;
  if (ks>0.75*p) ks = p; /* stretch the forward recurs as far as poss. */
  /* forward recursion for k < ks */
  for (k=2; k<ks; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = (y1-y0-cf[1]*I[k-1]-(k-1)*I[k-2])/(2*cf[2]);
  }
  if (ks==p) return;

  km = p+15;
  y1 *= l1*l1; y0 *= l0*l0;
  for (k=ks; k<=km; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }
  I[km+1] = I[km+2] = 0.0;
  for (k=km; k>=ks; k--)
    I[k] = (I[k]-cf[1]*I[k+1]-2*cf[2]*I[k+2])/(k+1);
}

void recent(I,resp,wt,p,s,x)
double *I, *resp, *wt, x;
INT p, s;
{ INT i, j;

  /* first, use W taylor series I -> resp */
  for (i=0; i<=p; i++)
  { resp[i] = 0.0;
    for (j=0; j<s; j++) resp[i] += wt[j]*I[i+j];
  }

  /* now, recenter x -> 0 */
  if (x==0) return;
  for (j=0; j<=p; j++) for (i=p; i>j; i--) resp[i] += x*resp[i-1];
}

void recurint(l0,l2,cf,resp,p,ker)
double l0, l2, *cf, *resp;
INT p, ker;
{ INT i, s;
  double l1, d0, d1, d2, dl, z0, z1, z2, wt[20], ncf[3], I[50], r1[5], r2[5];
if (debug) printf("\nrecurint: %8.5f %8.5f %8.5f   %8.5f %8.5f\n",cf[0],cf[1],cf[2],l0,l2);

  if (cf[2]==0) /* go straight to explint1 */
  { s = wtaylor(wt,0.0,ker);
if (debug) printf("case 1\n");
    explint1(l0,l2,cf,I,p+s);
    recent(I,resp,wt,p,s,0.0);
    return;
  }

  dl = l2-l0;
  d0 = cf[1]+2*l0*cf[2];
  d2 = cf[1]+2*l2*cf[2];
  z0 = cf[0]+l0*(cf[1]+l0*cf[2]);
  z2 = cf[0]+l2*(cf[1]+l2*cf[2]);

  if ((fabs(cf[1]*dl)<1) && (fabs(cf[2]*dl*dl)<1))
  { ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
if (debug) printf("case 2\n");
    s = wtaylor(wt,l0,ker);
    explinbkr(0.0,dl,ncf,I,p+s);
    recent(I,resp,wt,p,s,l0);
    return;
  }

  if (fabs(cf[2]*dl*dl)<0.001) /* small c, use explint1+tay.ser */
  { ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
if (debug) printf("case small c\n");
    s = wtaylor(wt,l0,ker);
    explintyl(0.0,l2-l0,ncf,I,p+s);
    recent(I,resp,wt,p,s,l0);
    return;
  }

  if (d0*d2<=0) /* max/min in [l0,l2] */
  { l1 = -cf[1]/(2*cf[2]);
    z1 = cf[0]+l1*(cf[1]+l1*cf[2]);
    d1 = 0.0;
    if (cf[2]<0) /* peak, integrate around l1 */
    { s = wtaylor(wt,l1,ker);
      ncf[0] = z1; ncf[1] = 0.0; ncf[2] = cf[2];
if (debug) printf("case peak  p %2d  s %2d\n",p,s);
      explinfbk0(l0-l1,l2-l1,ncf,I,p+s);
      recent(I,resp,wt,p,s,l1);
      return;
    }
  }

  if ((d0-2*cf[2]*dl)*(d2+2*cf[2]*dl)<0) /* max/min is close to [l0,l2] */
  { l1 = -cf[1]/(2*cf[2]);
    z1 = cf[0]+l1*(cf[1]+l1*cf[2]);
    if (l1<l0) { l1 = l0; z1 = z0; }
    if (l1>l2) { l1 = l2; z1 = z2; }

    if ((z1>=z0) & (z1>=z2)) /* peak; integrate around l1 */
    { s = wtaylor(wt,l1,ker);
if (debug) printf("case 4\n");
      d1 = cf[1]+2*l1*cf[2];
      ncf[0] = z1; ncf[1] = d1; ncf[2] = cf[2];
      explinfbk(l0-l1,l2-l1,ncf,I,p+s);
      recent(I,resp,wt,p,s,l1);
      return;
    }

    /* trough; integrate [l0,l1] and [l1,l2] */
    for (i=0; i<=p; i++) r1[i] = r2[i] = 0.0;
    if (l0<l1)
    { s = wtaylor(wt,l0,ker);
if (debug) printf("case 5\n");
      ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
      explinfbk(0.0,l1-l0,ncf,I,p+s);
      recent(I,r1,wt,p,s,l0);
    }
    if (l1<l2)
    { s = wtaylor(wt,l2,ker);
if (debug) printf("case 6\n");
      ncf[0] = z2; ncf[1] = d2; ncf[2] = cf[2];
      explinfbk(l1-l2,0.0,ncf,I,p+s);
      recent(I,r2,wt,p,s,l2);
    }
    for (i=0; i<=p; i++) resp[i] = r1[i]+r2[i];
    return;
  }

  /* Now, quadratic is monotone on [l0,l2]; big b; moderate c */
  if (z2>z0+3) /* steep increase, expand around l2 */
  { s = wtaylor(wt,l2,ker);
if (debug) printf("case 7\n");


    ncf[0] = z2; ncf[1] = d2; ncf[2] = cf[2];
    explinsid(l0-l2,0.0,ncf,I,p+s);
    recent(I,resp,wt,p,s,l2);
if (debug) printf("7 resp: %8.5f %8.5f %8.5f %8.5f\n",resp[0],resp[1],resp[2],resp[3]);
    return;
  }

  /* bias towards expansion around l0, because it's often 0 */
if (debug) printf("case 8\n");
  s = wtaylor(wt,l0,ker);
  ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
  explinsid(0.0,l2-l0,ncf,I,p+s);
  recent(I,resp,wt,p,s,l0);
  return;
}

INT onedexpl(cf,mi,resp)
double *cf, *resp;
INT *mi;
{ INT i;
  double f0, fr, fl;
  if (mi[MDEG]>=2) ERROR(("onedexpl only valid for deg=0,1"));
  if (fabs(cf[1])>=EFACT) return(LF_BADP);

  f0 = exp(cf[0]); fl = fr = 1.0;
  for (i=0; i<=2*mi[MDEG]; i++)
  { f0 *= i+1;
    fl /=-(EFACT+cf[1]);
    fr /=  EFACT-cf[1];
    resp[i] = f0*(fr-fl);
  }
  return(LF_OK);
}

INT onedgaus(cf,mi,resp)
double *cf, *resp;
INT *mi;
{ INT i;
  double f0, mu, s2;
  if (mi[MDEG]>=3)
  { ERROR(("onedgaus only valid for deg=0,1,2"));
    return(LF_ERR);
  }
  if (2*cf[2]>=GFACT*GFACT) return(LF_BADP);

  s2 = 1/(GFACT*GFACT-2*cf[2]);
  mu = cf[1]*s2;
  resp[0] = 1.0;
  if (mi[MDEG]>=1)
  { resp[1] = mu;
    resp[2] = s2+mu*mu;
    if (mi[MDEG]==2)
    { resp[3] = mu*(3*s2+mu*mu);
      resp[4] = 3*s2*s2 + mu*mu*(6*s2+mu*mu);
    }
  }
  f0 = S2PI * exp(cf[0]+mu*mu/(2*s2))*sqrt(s2);
  for (i=0; i<=2*mi[MDEG]; i++) resp[i] *= f0;
  return(LF_OK);
}

INT onedint(cf,mi,l0,l1,resp) /* int W(u)u^j exp(..), j=0..2*deg */
double *cf, l0, l1, *resp;
INT *mi;
{ double u, uj, y, ncf[4], rr[5];
  INT deg, i, j;
    memset(rr, 0, sizeof(rr));
if (debug) printf("onedint: %f %f %f   %f %f\n",cf[0],cf[1],cf[2],l0,l1);
  deg = mi[MDEG];

  if (deg<=2)
  { for (i=0; i<3; i++) ncf[i] = (i>deg) ? 0.0 : cf[i];
    ncf[2] /= 2;

    if (mi[MKER]==WEXPL) return(onedexpl(ncf,mi,resp));
    if (mi[MKER]==WGAUS) return(onedgaus(ncf,mi,resp));

    if (l1>0)
      recurint(MAX(l0,0.0),l1,ncf,resp,2*deg,mi[MKER]);
    else for (i=0; i<=2*deg; i++) resp[i] = 0;

    if (l0<0)
    { ncf[1] = -ncf[1];
      l0 = -l0; l1 = -l1;
      recurint(MAX(l1,0.0),l0,ncf,rr,2*deg,mi[MKER]);
    }
    else for (i=0; i<=2*deg; i++) rr[i] = 0.0;

    for (i=0; i<=2*deg; i++)
      resp[i] += (i%2==0) ? rr[i] : -rr[i];

    return(LF_OK);
  }

  /* For degree >= 3, we use Simpson's rule. */
  for (j=0; j<=2*deg; j++) resp[j] = 0.0;
  for (i=0; i<=mi[MMINT]; i++)
  { u = l0+(l1-l0)*i/mi[MMINT];
    y = cf[0]; uj = 1;
    for (j=1; j<=deg; j++)
    { uj *= u;
      y += cf[j]*uj/fact[j];
    }
    y = (4-2*(i%2==0)-(i==0)-(i==mi[MMINT])) *
          W(fabs(u),mi[MKER])*exp(MIN(y,300.0));
    for (j=0; j<=2*deg; j++)
    { resp[j] += y;
      y *= u;
    }
  }
  for (j=0; j<=2*deg; j++) resp[j] = resp[j]*(l1-l0)/(3*mi[MMINT]);
  return(LF_OK);
}

