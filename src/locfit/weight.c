/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *  Defines the weight functions and related quantities used
 *  in LOCFIT.
 */

#include "local.h"


/* The weight functions themselves.  Used everywhere. */
double W(u,ker)
double u;
INT ker;
{ u = fabs(u);
  switch(ker)
  { case WRECT: return((u>1) ? 0.0 : 1.0);
    case WEPAN: return((u>1) ? 0.0 : 1-u*u);
    case WBISQ: if (u>1) return(0.0);
                u = 1-u*u; return(u*u);
    case WTCUB: if (u>1) return(0.0);
                u = 1-u*u*u; return(u*u*u);
    case WTRWT: if (u>1) return(0.0);
                u = 1-u*u; return(u*u*u);
    case WQUQU: if (u>1) return(0.0);
                u = 1-u*u; return(u*u*u*u);
    case WTRIA: if (u>1) return(0.0);
                return(1-u);
    case W6CUB: if (u>1) return(0.0);
                u = 1-u*u*u; u = u*u*u; return(u*u);
    case WGAUS: return(exp(-SQR(GFACT*u)/2.0));
    case WEXPL: return(exp(-EFACT*u));
    case WMACL: return(1/((u+1.0e-100)*(u+1.0e-100)));
    case WMINM: ERROR(("WMINM in W"));
                return(0.0);
    case WPARM: return(1.0);
  }
  return(0.0);
}

INT iscompact(ker)
INT ker;
{ if ((ker==WEXPL) | (ker==WGAUS) | (ker==WMACL) | (ker==WPARM)) return(0);
  return(1);
}

double weightprod(lf,u,h)
lfit *lf;
double *u, h;
{ INT i, ker;
  double sc, w;
  w = 1.0;
  ker = lf->mi[MKER];
  for (i=0; i<lf->mi[MDIM]; i++)
  { sc = lf->sca[i];
    switch(lf->sty[i])
    { case STLEFT:
        if (u[i]>0) return(0.0);
        w *= W(-u[i]/(h*sc),ker);
        break;
      case STRIGH:
        if (u[i]<0) return(0.0);
        w *= W(u[i]/(h*sc),ker);
        break;
      case STANGL:
        w *= W(2*fabs(sin(u[i]/(2*sc)))/h,ker);
        break;
      case STCPAR:
        break;
      default:
        w *= W(fabs(u[i])/(h*sc),ker);
    }
    if (w==0.0) return(w);
  }
  return(w);
}

double weightsph(lf,u,h,hasdi,di)
lfit *lf;
double *u, h, di;
INT hasdi;
{ INT i;

  if (!hasdi) di = rho(u,lf->sca,lf->mi[MDIM],lf->mi[MKT],lf->sty);

  for (i=0; i<lf->mi[MDIM]; i++)
  { if ((lf->sty[i]==STLEFT) && (u[i]>0.0)) return(0.0);
    if ((lf->sty[i]==STRIGH) && (u[i]<0.0)) return(0.0);
  }
  if (h==0) return((di==0.0) ? 1.0 : 0.0);

  return(W(di/h,lf->mi[MKER]));
}

double weight(lf,x,t,h,hasdi,di)
lfit *lf;
double *x, *t, h, di;
INT hasdi;
{ double u[MXDIM];
  INT i;
  for (i=0; i<lf->mi[MDIM]; i++) u[i] = (t==NULL) ? x[i] : x[i]-t[i];
  switch(lf->mi[MKT])
  { case KPROD: return(weightprod(lf,u,h));
    case KSPH:  return(weightsph(lf,u,h,hasdi,di));
  }
  ERROR(("weight: unknown kernel type %d",lf->mi[MKT]));
  return(1.0);
}

double sgn(x)
double x;
{ if (x>0) return(1.0);
  if (x<0) return(-1.0);
  return(0.0);
}

double WdW(u,ker) /* W'(u)/W(u) */
double u;
INT ker;
{ double eps=1.0e-10;
  if (ker==WGAUS) return(-GFACT*GFACT*u);
  if (ker==WPARM) return(0.0);
  if (fabs(u)>=1) return(0.0);
  switch(ker)
  { case WRECT: return(0.0);
    case WTRIA: return(-sgn(u)/(1-fabs(u)+eps));
    case WEPAN: return(-2*u/(1-u*u+eps));
    case WBISQ: return(-4*u/(1-u*u+eps));
    case WTRWT: return(-6*u/(1-u*u+eps));
    case WTCUB: return(-9*sgn(u)*u*u/(1-u*u*fabs(u)+eps));
    case WEXPL: return((u>0) ? -EFACT : EFACT);
  }
  ERROR(("WdW: invalid kernel"));
  return(0.0);
}

/* deriv. weights .. spherical, product etc
   u, sc, sty needed only in relevant direction
   Acutally, returns (d/dx W(||x||/h) ) / W(.)
*/
double weightd(u,sc,d,ker,kt,h,sty,di)
double u, sc, h, di;
INT d, ker, kt, sty;
{ if (sty==STANGL)
  { if (kt==KPROD)
      return(-WdW(2*sin(u/(2*sc)),ker)*cos(u/(2*sc))/(h*sc));
    if (di==0.0) return(0.0);
    return(-WdW(di/h,ker)*sin(u/sc)/(h*sc*di));
  }
  if (sty==STCPAR) return(0.0);
  if (kt==KPROD)
    return(-WdW(u/(h*sc),ker)/(h*sc));
  if (di==0.0) return(0.0);
  return(-WdW(di/h,ker)*u/(h*di*sc*sc));
}

double weightdd(u,sc,d,ker,kt,h,sty,di,i0,i1)
double *u, *sc, h, di;
INT d, ker, kt, *sty, i0, i1;
{ double w;
  w = 1;
  if (kt==KPROD)
  {
    w = WdW(u[i0]/(h*sc[i0]),ker)*WdW(u[i1]/(h*sc[i1]),ker)/(h*h*sc[i0]*sc[i1]);
  }
  return(0.0);
}

/* Derivatives W'(u)/u.
   Used in simult. conf. band computations,
   and kernel density bandwidth selectors. */
double Wd(u,ker)
double u;
INT ker;
{ double v;
  if (ker==WGAUS) return(-SQR(GFACT)*exp(-SQR(GFACT*u)/2));
  if (ker==WPARM) return(0.0);
  if (fabs(u)>1) return(0.0);
  switch(ker)
  { case WEPAN: return(-2.0);
    case WBISQ: return(-4*(1-u*u));
    case WTCUB: v = 1-u*u*u;
                return(-9*v*v*u);
    case WTRWT: v = 1-u*u;
                return(-6*v*v);
    default: ERROR(("Invalid kernel %d in Wd",ker));
  }
  return(0.0);
}

/* Second derivatives W''(u)-W'(u)/u.
   used in simult. conf. band computations in >1 dimension. */
double Wdd(u,ker)
double u;
INT ker;
{ double v;
  if (ker==WGAUS) return(SQR(u*GFACT*GFACT)*exp(-SQR(u*GFACT)/2));
  if (ker==WPARM) return(0.0);
  if (u>1) return(0.0);
  switch(ker)
  { case WBISQ: return(12*u*u);
    case WTCUB: v = 1-u*u*u;
                return(-9*u*v*v+54*u*u*u*u*v);
    case WTRWT: return(24*u*u*(1-u*u));
    default: ERROR(("Invalid kernel %d in Wdd",ker));
  }
  return(0.0);
}

/* int u1^j1..ud^jd W(u) du.
   Used for local log-linear density estimation.
   Assume all j_i are even.
   Also in some bandwidth selection.
*/
double wint(d,j,nj,ker)
INT d, *j, nj, ker;
{ double I, z;
  int k, dj;
  dj = d;
    I = 0.0;
  for (k=0; k<nj; k++) dj += j[k];
  switch(ker) /* int_0^1 u^(dj-1) W(u)du  */
  { case WRECT: I = 1.0/dj; break;
    case WEPAN: I = 2.0/(dj*(dj+2)); break;
    case WBISQ: I = 8.0/(dj*(dj+2)*(dj+4)); break;
    case WTCUB: I = 162.0/(dj*(dj+3)*(dj+6)*(dj+9)); break;
    case WTRWT: I = 48.0/(dj*(dj+2)*(dj+4)*(dj+6)); break;
    case WTRIA: I = 1.0/(dj*(dj+1)); break;
    case WQUQU: I = 384.0/(dj*(dj+2)*(dj+4)*(dj+6)*(dj+8)); break;
    case W6CUB: I = 524880.0/(dj*(dj+3)*(dj+6)*(dj+9)*(dj+12)*(dj+15)*(dj+18)); break;
    case WGAUS: switch(d)
                { case 1: I = S2PI/GFACT; break;
                  case 2: I = 2*PI/(GFACT*GFACT); break;
                  default: I = exp(d*log(S2PI/GFACT)); /* for nj=0 */
                }
                for (k=0; k<nj; k++) /* deliberate drop */
                  switch(j[k])
                  { case 4: I *= 3.0/(GFACT*GFACT);
                    case 2: I /= GFACT*GFACT;
                  }
                return(I);
    case WEXPL: I = factorial(dj-1)/ipower(EFACT,dj); break;
    default: ERROR(("Unknown kernel %d in exacint",ker));
  }
  if ((d==1) && (nj==0)) return(2*I); /* common case quick */
  z = (d-nj)*LOGPI/2-LGAMMA(dj/2.0);
  for (k=0; k<nj; k++) z += LGAMMA((j[k]+1)/2.0);
  return(2*I*exp(z));
}

/* taylor series expansion of weight function around x.
   0 and 1 are common arguments, so are worth programming
   as special cases.
   Used in density estimation.
*/
INT wtaylor(f,x,ker)
double *f, x;
INT ker;
{ double v;
  switch(ker)
  { case WRECT:
      f[0] = 1.0;
      return(1);
    case WEPAN:
      f[0] = 1-x*x; f[1] = -2*x; f[2] = -1;
      return(3);
    case WBISQ:
      v = 1-x*x;
      f[0] = v*v;   f[1] = -4*x*v; f[2] = 4-6*v;
      f[3] = 4*x;   f[4] = 1;
      return(5);
    case WTCUB:
      if (x==1.0)
      { f[0] = f[1] = f[2] = 0; f[3] = -27; f[4] = -81; f[5] = -108;
        f[6] = -81; f[7] = -36; f[8] = -9; f[9] = -1; return(10); }
      if (x==0.0)
      { f[1] = f[2] = f[4] = f[5] = f[7] = f[8] = 0;
        f[0] = 1; f[3] = -3; f[6] = 3; f[9] = -1; return(10); }
      v = 1-x*x*x;
      f[0] = v*v*v; f[1] = -9*v*v*x*x; f[2] = x*v*(27-36*v);
      f[3] = -27+v*(108-84*v);         f[4] = -3*x*x*(27-42*v);
      f[5] = x*(-108+126*v);           f[6] = -81+84*v;
      f[7] = -36*x*x; f[8] = -9*x;     f[9] = -1;
      return(10);
    case WTRWT:
      v = 1-x*x;
      f[0] = v*v*v; f[1] = -6*x*v*v; f[2] = v*(12-15*v);
      f[3] = x*(20*v-8); f[4] = 15*v-12; f[5] = -6; f[6] = -1;
      return(7);
    case WTRIA:
      f[0] = 1-x; f[1] = -1;
      return(2);
    case WQUQU:
      v = 1-x*x;
      f[0] = v*v*v*v; f[1] = -8*x*v*v*v; f[2] = v*v*(24-28*v);
      f[3] = v*x*(56*v-32); f[4] = (70*v-80)*v+16; f[5] = x*(32-56*v);
      f[6] = 24-28*v; f[7] = 8*x; f[8] = 1;
      return(9);
    case W6CUB:
      v = 1-x*x*x;
      f[0] = v*v*v*v*v*v;
      f[1] = -18*x*x*v*v*v*v*v;
      f[2] = x*v*v*v*v*(135-153*v);
      f[3] = v*v*v*(-540+v*(1350-816*v));
      f[4] = x*x*v*v*(1215-v*(4050-v*3060));
      f[5] = x*v*(-1458+v*(9234+v*(-16254+v*8568)));
      f[6] = 729-v*(10206-v*(35154-v*(44226-v*18564)));
      f[7] = x*x*(4374-v*(30132-v*(56862-v*31824)));
      f[8] = x*(12393-v*(61479-v*(92664-v*43758)));
      f[9] = 21870-v*(89100-v*(115830-v*48620));
      f[10]= x*x*(26730-v*(69498-v*43758));
      f[11]= x*(23814-v*(55458-v*31824));
      f[12]= 15849-v*(34398-v*18564);
      f[13]= x*x*(7938-8568*v);
      f[14]= x*(2970-3060*v);
      f[15]= 810-816*v;
      f[16]= 153*x*x;
      f[17]= 18*x;
      f[18]= 1;
      return(19);
  }
  ERROR(("Invalid kernel %d in wtaylor",ker));
  return(0);
}

/* convolution int W(x)W(x+v)dx.
   used in kde bandwidth selection.
*/
double Wconv(v,ker)
double v;
INT ker;
{ double v2;
  switch(ker)
  { case WGAUS: return(SQRPI/GFACT*exp(-SQR(GFACT*v)/4));
    case WRECT:
      v = fabs(v);
      if (v>2) return(0.0);
      return(2-v);
    case WEPAN:
      v = fabs(v);
      if (v>2) return(0.0);
      return((2-v)*(16+v*(8-v*(16-v*(2+v))))/30);
    case WBISQ:
      v = fabs(v);
      if (v>2) return(0.0);
      v2 = 2-v;
      return(v2*v2*v2*v2*v2*(16+v*(40+v*(36+v*(10+v))))/630);
  }
  ERROR(("Wconv not implemented for kernel %d",ker));
  return(0.0);
}

/* derivative of Wconv.
   1/v d/dv int W(x)W(x+v)dx
   used in kde bandwidth selection.
*/
double Wconv1(v,ker)
double v;
INT ker;
{ double v2;
  v = fabs(v);
  switch(ker)
  { case WGAUS: return(-0.5*SQRPI*GFACT*exp(-SQR(GFACT*v)/4));
    case WRECT:
      if (v>2) return(0.0);
      return(1.0);
    case WEPAN:
      if (v>2) return(0.0);
      return((-16+v*(12-v*v))/6);
    case WBISQ:
      if (v>2) return(0.0);
      v2 = 2-v;
      return(-v2*v2*v2*v2*(32+v*(64+v*(24+v*3)))/210);
  }
  ERROR(("Wconv1 not implemented for kernel %d",ker));
  return(0.0);
}

/* 4th derivative of Wconv.
   used in kde bandwidth selection (BCV, SJPI, GKK)
*/
double Wconv4(v,ker)
double v;
INT ker;
{ double gv;
  switch(ker)
  { case WGAUS:
      gv = GFACT*v;
      return(exp(-SQR(gv)/4)*GFACT*GFACT*GFACT*(12-gv*gv*(12-gv*gv))*SQRPI/16);
  }
  ERROR(("Wconv4 not implemented for kernel %d",ker));
  return(0.0);
}

/* 5th derivative of Wconv.
   used in kde bandwidth selection (BCV method only)
*/
double Wconv5(v,ker) /* (d/dv)^5 int W(x)W(x+v)dx */
double v;
INT ker;
{ double gv;
  switch(ker)
  { case WGAUS:
      gv = GFACT*v;
      return(-exp(-SQR(gv)/4)*GFACT*GFACT*GFACT*GFACT*gv*(60-gv*gv*(20-gv*gv))*SQRPI/32);
  }
  ERROR(("Wconv5 not implemented for kernel %d",ker));
  return(0.0);
}

/* 6th derivative of Wconv.
   used in kde bandwidth selection (SJPI)
*/
double Wconv6(v,ker)
double v;
INT ker;
{ double gv, z;
  switch(ker)
  { case WGAUS:
      gv = GFACT*v;
      gv = gv*gv;
      z = exp(-gv/4)*(-120+gv*(180-gv*(30-gv)))*0.02769459142;
      gv = GFACT*GFACT;
      return(z*gv*gv*GFACT);
  }
  ERROR(("Wconv6 not implemented for kernel %d",ker));
  return(0.0);
}

/* int W(v)^2 dv / (int v^2 W(v) dv)^2
   used in some bandwidth selectors
*/
double Wikk(ker,deg)
INT ker, deg;
{ switch(deg)
  { case 0:
    case 1: /* int W(v)^2 dv / (int v^2 W(v) dv)^2 */
      switch(ker)
      { case WRECT: return(4.5);
        case WEPAN: return(15.0);
        case WBISQ: return(35.0);
        case WGAUS: return(0.2820947918*GFACT*GFACT*GFACT*GFACT*GFACT);
        case WTCUB: return(34.15211105);
        case WTRWT: return(66.08391608);
      }
    case 2:
    case 3: /* 4!^2/8*int(W1^2)/int(v^4W1)^2
               W1=W*(n4-v^2n2)/(n0n4-n2n2) */
      switch(ker)
      { case WRECT: return(11025.0);
        case WEPAN: return(39690.0);
        case WBISQ: return(110346.9231);
        case WGAUS: return(14527.43412);
        case WTCUB: return(126500.5904);
        case WTRWT: return(254371.7647);
      }
  }
  ERROR(("Wikk not implemented for kernel %d",ker));
  return(0.0);
}
