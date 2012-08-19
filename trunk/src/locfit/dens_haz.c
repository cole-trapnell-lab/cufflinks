/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *   Integration for hazard rate estimation. The functions in this
 *   file are used to evaluate
 *      sum int_0^{Ti} W_i(t,x) A()A()' exp( P() ) dt
 *   for hazard rate models.
 *
 *   These routines assume the weight function is supported on [-1,1].
 *   hasint_sph multiplies by exp(base(lf,i)), which allows estimating
 *   the baseline in a proportional hazards model, when the covariate
 *   effect base(lf,i) is known.
 *
 *   TODO:
 *     hazint_sph, should be able to reduce mint in some cases with
 *       small integration range. onedint could be used for beta-family
 *       (RECT,EPAN,BISQ,TRWT) kernels.
 *     hazint_prod, restrict terms from the sum based on x values.
 *       I should count obs >= max, and only do that integration once.
 */

#include "local.h"

static double ilim[2*MXDIM], *ff, tmax;

/*
 *  hrao returns 0 if integration region is empty.
 *               1 otherwise.
 */
INT haz_sph_int(lf,dfx,cf,h,r1)
lfit *lf;
double *dfx, *cf, h, *r1;
{ double s, t0, t1, wt, th;
  INT dim, j, p, *mi;
  mi = lf->mi;
  s = 0; p = mi[MP];
  dim = mi[MDIM];
  for (j=1; j<dim; j++) s += SQR(dfx[j]/(h*lf->sca[j]));
  if (s>1) return(0);

  setzero(r1,p*p);
  t1 = sqrt(1-s)*h*lf->sca[0];
  t0 = -t1;
  if (t0<ilim[0])   t0 = ilim[0];
  if (t1>ilim[dim]) t1 = ilim[dim];
  if (t1>dfx[0]) t1 = dfx[0];
  if (t1<t0) return(0);

/*  Numerical integration by Simpson's rule.
 */
  for (j=0; j<=mi[MMINT]; j++)
  { dfx[0] = t0+(t1-t0)*j/mi[MMINT];
    wt = weight(lf,dfx,NULL,h,0,0.0);
    fitfun(lf,dfx,NULL,ff,NULL,0);
    th = innerprod(cf,ff,p);
    if (mi[MLINK]==LLOG) th = exp(th);
    wt *= 2+2*(j&1)-(j==0)-(j==mi[MMINT]);
    addouter(r1,ff,ff,p,wt*th);
  }
  multmatscal(r1,(t1-t0)/(3*mi[MMINT]),p*p);

  return(1);
}

INT hazint_sph(t,resp,r1,lf,cf,h)
lfit *lf;
double *t, *resp, *r1, *cf, h;
{ INT i, j, p, st;
  double dfx[MXDIM], eb, sb;
  p = lf->mi[MP];
  setzero(resp,p*p);
  sb = 0.0;

  for (i=0; i<=lf->mi[MN]; i++)
  {
    if (i==lf->mi[MN])
    { dfx[0] = tmax-t[0];
      for (j=1; j<lf->mi[MDIM]; j++) dfx[j] = 0.0;
      eb = exp(sb/lf->mi[MN]);
    }
    else
    { eb = exp(base(lf,i)); sb += base(lf,i);
      for (j=0; j<lf->mi[MDIM]; j++) dfx[j] = datum(lf,j,i)-t[j];
    }

    st = haz_sph_int(lf,dfx,cf,h,r1);
    if (st)
      for (j=0; j<p*p; j++) resp[j] += eb*r1[j];
  }
  return(LF_OK);
}

INT hazint_prod(t,resp,x,lf,cf,h)
lfit *lf;
double *t, *resp, *x, *cf, h;
{ INT d, p, deg, i, j, k, st;
  double dfx[MXDIM], t_prev,
         hj, hs, ncf[MXDEG], ef, il1;
  double prod_wk[MXDIM][2*MXDEG+1], eb, sb;

  p = lf->mi[MP]; d = lf->mi[MDIM];
  deg = lf->mi[MDEG];
  setzero(resp,p*p);
  hj = hs = h*lf->sca[0];
    memset(dfx, 0.0, sizeof(dfx));
  ncf[0] = cf[0];
  for (i=1; i<=deg; i++)
  { ncf[i] = hj*cf[(i-1)*d+1]; hj *= hs;
  }

/*   for i=0..n....
 *     First we compute prod_wk[j], j=0..d.
 *     For j=0, this is int_0^T_i (u-t)^k W((u-t)/h) exp(b0*(u-t)) du
 *     For remaining j,   (x(i,j)-x(j))^k Wj exp(bj*(x..-x.))
 *
 *     Second, we add to the integration (exp(a) incl. in integral)
 *     with the right factorial denominators.
 */
  t_prev = ilim[0]; sb = 0.0;
  for (i=0; i<=lf->mi[MN]; i++)
  { if (i==lf->mi[MN])
    { dfx[0] = tmax-t[0];
      for (j=1; j<d; j++) dfx[j] = 0.0;
      eb = exp(sb/lf->mi[MN]);
    }
    else
    { eb = exp(base(lf,i)); sb += base(lf,i);
      for (j=0; j<d; j++) dfx[j] = datum(lf,j,i)-t[j];
    }

    if (dfx[0]>ilim[0]) /* else it doesn't contribute */
    {
/* time integral */
      il1 = (dfx[0]>ilim[d]) ? ilim[d] : dfx[0];
      if (il1 != t_prev) /* don't repeat! */
      { st = onedint(ncf,lf->mi,ilim[0]/hs,il1/hs,prod_wk[0]);
        if (st>0) return(st);
        hj = eb;
        for (j=0; j<=2*deg; j++)
        { hj *= hs;
          prod_wk[0][j] *= hj;
        }
        t_prev = il1;
      }

/* covariate terms */
      for (j=1; j<d; j++)
      {
        ef = 0.0;
        for (k=deg; k>0; k--) ef = (ef+dfx[j])*cf[1+(k-1)*d+j];
        ef = exp(ef);
        prod_wk[j][0] = ef * W(dfx[j]/(h*lf->sca[j]),lf->mi[MKER]);
        for (k=1; k<=2*deg; k++)
          prod_wk[j][k] = prod_wk[j][k-1] * dfx[j];
      }

/*  add to the integration.  */
      prodint_resp(resp,prod_wk,d,deg,p);
    } /* if dfx0 > ilim0 */
  } /* n loop */

/* symmetrize */
  for (k=0; k<p; k++)
    for (j=k; j<p; j++)
      resp[j*p+k] = resp[k*p+j];
  return(LF_OK);
}

INT hazint(t,resp,resp1,lf,cf,h)
lfit *lf;
double *t, *resp, *resp1, *cf, h;
{ if (lf->mi[MDIM]==1) return(hazint_prod(t,resp,resp1,lf,cf,h));
  if (lf->mi[MKT]==KPROD) return(hazint_prod(t,resp,resp1,lf,cf,h));

  return(hazint_sph(t,resp,resp1,lf,cf,h));
}

void haz_init(lf,des,il)
lfit *lf;
design *des;
double *il;
{ int i;
  tmax = datum(lf,0,0);
  for (i=1; i<lf->mi[MN]; i++) tmax = MAX(tmax,datum(lf,0,i));
  ff = des->xtwx.wk;
  for (i=0; i<2*lf->mi[MDIM]; i++) ilim[i] = il[i];
}
