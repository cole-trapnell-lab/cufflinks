/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   The function dens_integrate(lf,des,z) is used to integrate a density
 *   estimate (z=1) or the density squared (z=2). This is used to renormalize
 *   the estimate (function dens_renorm) or in the computation of LSCV
 *   (function dnes_lscv). The implementation is presently for d=1.
 *
 *   The computation orders the fit points selected by locfit, and
 *   integrates analytically over each interval. For the log-link,
 *   the interpolant used is peicewise quadratic (with one knot in
 *   the middle of each interval); this differs from the cubic interpolant
 *   used elsewhere in Locfit.
 *
 *   TODO: allow for xlim. What can be done simply in >=2 dimensions?
 *         fix df computation (in lscv) for link=IDENT.
 */

#include "local.h"

/*
 * Finds the order of observations in the array x, and
 * stores in integer array ind.
 * At input, lset l=0 and r=length(x)-1.
 * At output, x[ind[0]] <= x[ind[1]] <= ...
 */
void lforder(ind,x,l,r)
INT *ind, l, r;
double *x;
{ double piv;
  INT i, i0, i1;
  piv = (x[ind[l]]+x[ind[r]])/2;
  i0 = l; i1 = r;
  while (i0<=i1)
  { while ((i0<=i1) && (x[ind[i0]]<=piv)) i0++;
    while ((i0<=i1) && (x[ind[i1]]>piv))  i1--;
    if (i0<i1)
    { ISWAP(ind[i0],ind[i1]);
      i0++; i1--;
    }
  }
  /* now, x[ind[l..i1]] <= piv < x[ind[i0..r]].
     put the ties in the middle */
  while ((i1>=l) && (x[ind[i1]]==piv)) i1--;
  for (i=l; i<=i1; i++)
    if (x[ind[i]]==piv)
    { ISWAP(ind[i],ind[i1]);
      while (x[ind[i1]]==piv) i1--;
    }

  if (l<i1) lforder(ind,x,l,i1);
  if (i0<r) lforder(ind,x,i0,r);
}

/*
 *  estdiv integrates the density between fit points x0 and x1.
 *  f0, f1 are function values, d0, d1 are derivatives.
 */
double estdiv(x0,x1,f0,f1,d0,d1,link)
double x0, x1, f0, f1, d0, d1;
INT link;
{ double cf[4], I[2], dlt, e0, e1;

  if (x0==x1) return(0.0);

  if (link==LIDENT)
  {
/* cf are integrals of hermite polynomials.
 * Then adjust for x1-x0.
 */
    cf[0] = cf[1] = 0.5;
    cf[2] = 1.0/12.0; cf[3] = -cf[2];
    return( (cf[0]*f0+cf[1]*f1)*(x1-x0)
          + (cf[2]*d0+cf[3]*d1)*(x1-x0)*(x1-x0) );
  }

/*
 * this is for LLOG
 */

  dlt = (x1-x0)/2;
  cf[0] = f0;
  cf[1] = d0;
  cf[2] = ( 2*(f1-f0) - dlt*(d1+3*d0) )/(4*dlt*dlt);
  recurint(0.0,dlt,cf,I,0,WRECT);
  e0 = I[0];

  cf[0] = f1;
  cf[1] = -d1;
  cf[2] = ( 2*(f0-f1) + dlt*(d0+3*d1) )/( 4*dlt*dlt );
  recurint(0.0,dlt,cf,I,0,WRECT);
  e1 = I[0];

  return(e0+e1);
}

/*
 *   Evaluate the integral of the density estimate to the power z.
 *   This would be severely messed up, if I ever implement parcomp
 *   for densities.
 */
double dens_integrate(lf,des,z)
lfit *lf;
design *des;
INT z;
{ INT has_deriv, i, i0, i1, link, nv, *ind;
  double *xev, *fit, *deriv, sum, term;
  double d0, d1, f0, f1;

  if (lf->mi[MDIM]>=2)
  { WARN(("dens_integrate requires d=1"));
    return(0.0);
  }

  link = lf->mi[MLINK];
  has_deriv = (lf->mi[MDEG] > 0); /* not right? */
  fit = lf->coef;
  if (has_deriv)
    deriv = &lf->coef[lf->nvm];
  xev = vdptr(lf->xxev);

  /*
   * order the vertices
   */
  nv = lf->nv;
  if (lf->mi[MN]<nv) return(0.0);
  ind = des->ind;
  for (i=0; i<nv; i++) ind[i] = i;
  lforder(ind,xev,0,nv-1);
  sum = 0.0;

  /*
   * Estimate the contribution of the boundaries.
   * should really check flim here.
   */
  i0 = ind[0]; i1 = ind[1];
  f1 = fit[i0];
  d1 = (has_deriv) ? deriv[i0] :
         (fit[i1]-fit[i0])/(xev[i1]-xev[i0]);
  if (d1 <= 0) WARN(("dens_integrate - ouch!"));
  if (z==2)
  { if (link==LLOG)
    { f1 *= 2; d1 *= 2; }
    else
    { d1 = 2*d1*f1; f1 = f1*f1; }
  }
  term = (link==LIDENT) ? f1*f1/(2*d1) : exp(f1)/d1;
  sum += term;

  i0 = ind[nv-2]; i1 = ind[nv-1];
  f0 = fit[i1];
  d0 = (has_deriv) ? deriv[i1] :
         (fit[i1]-fit[i0])/(xev[i1]-xev[i0]);
  if (d0 >= 0) WARN(("dens_integrate - ouch!"));
  if (z==2)
  { if (link==LLOG)
    { f0 *= 2; d0 *= 2; }
    else
    { d0 = 2*d0*f0; f0 = f0*f0; }
  }
  term = (link==LIDENT) ? -f0*f0/(2*d0) : exp(f0)/d0;
  sum += term;
  
  for (i=1; i<nv; i++)
  { i0 = ind[i-1]; i1 = ind[i];
    f0 = fit[i0]; f1 = fit[i1];
    d0 = (has_deriv) ? deriv[i0] :
              (f1-f0)/(xev[i1]-xev[i0]);
    d1 = (has_deriv) ? deriv[i1] : d0;
    if (z==2)
    { if (link==LLOG)
      { f0 *= 2; f1 *= 2; d0 *= 2; d1 *= 2; }
      else
      { d0 *= 2*f0; d1 *= 2*f1; f0 = f0*f0; f1 = f1*f1; }
    }
    term = estdiv(xev[i0],xev[i1],f0,f1,d0,d1,link);
    sum += term;
  }

  return(sum);
}

void dens_renorm(lf,des)
lfit *lf;
design *des;
{ INT i;
  double sum;
  sum = dens_integrate(lf,des,1);
  if (sum==0.0) return;
  sum = log(sum);
  for (i=0; i<lf->nv; i++) lf->coef[i] -= sum;
}

void dens_lscv(des,lf)
lfit *lf;
design *des;
{ double df, fh, fh_cv, infl, z0, z1, x[MXDIM];
  int i, n, j, ev;
  z1 = df = 0.0;
  ev = lf->mi[MEV];
  n = lf->mi[MN];
  if ((ev==EDATA) | (ev==ECROS)) ev = EFITP;

  z0 = dens_integrate(lf,des,2);

  for (i=0; i<n; i++)
  { for (j=0; j<lf->mi[MDIM]; j++) x[j] = datum(lf,j,i);
    fh = base(lf,i)+dointpoint(lf,des,x,PCOEF,ev,i);
    if (lf->mi[MLINK]==LLOG) fh = exp(fh);
    infl = dointpoint(lf,des,x,PT0,ev,i);
    infl = infl * infl;
    if (infl>1) infl = 1;
    fh_cv = (lf->mi[MLINK] == LIDENT) ?
       (n*fh - infl) / (n-1.0) : fh*(1-infl)*n/(n-1.0);
    z1 += fh_cv;
    df += infl;
  }

  vdptr(lf->L)[0] = z0-2*z1/n;
  vdptr(lf->L)[1] = df;
}
