/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

/*
  Functions implementing the adaptive bandwidth selection.
  Will make the final call to nbhd() to set smoothing weights
  for selected bandwidth, But will **not** make the
  final call to locfit().
*/

#include "local.h"

static double hmin;

double acri(lk,t0,t2,pen)
double lk, t0, t2, pen;
{ double y;
/* return(-2*lk/(t0*exp(pen*log(1-t2/t0)))); */
  /* return((-2*lk+pen*t2)/t0); */
  y = (MAX(-2*lk,t0-t2)+pen*t2)/t0;
  return(y);
}

double mmse(lf,des)
lfit *lf;
design *des;
{ int i, ii, j, p, p1;
  double sv, sb, *l, dp;
  l = des->wd;
  wdiag(lf,des,l,(INT)0,(INT)1,(INT)0);
  sv = sb = 0;
  p = lf->mi[MP];
  for (i=0; i<des->n; i++)
  { sv += l[i]*l[i];
    ii = des->ind[i];
    dp = des->di[ii];
    for (j=0; j<lf->mi[MDEG]; j++) dp *= des->di[ii];
    sb += fabs(l[i])*dp;
  }
  p1 = factorial((int)lf->mi[MDEG]+1);
  return(sv+sb*sb*lf->dp[DADP]*lf->dp[DADP]/(p1*p1));
}

static double mcp, clo, cup;

/*
  Initial bandwidth will be (by default)
  k-nearest neighbors for k small, just lage enough to
  get defined estimate (unless user provided nonzero DALP
  or DFXH components)
*/

INT ainitband(des,lf)
design *des;
lfit   *lf;
{ INT lf_status = LF_OK, p, z, cri, noit, redo;
  double h, ho, t[6];
  p = des->p;
  cri = lf->mi[MACRI];
  noit = !((cri==AOK) | (cri==ANONE));
  z = (INT)(lf->mi[MN]*lf->dp[DALP]);
  if ((noit) && (z<p+2)) z = p+2;
  redo = 0; ho = -1;
  do
  { h = nbhd(lf,des,z,lf->dp[DFXH],redo);
    if (z<des->n) z = des->n;
    if (h>ho) lf_status = locfit(lf,des,h,noit);
    if (cri==ANONE) return(lf_status);
    z++;
    redo = 1;
  } while ((z<=lf->mi[MN]) && ((h==0)||(lf_status!=LF_OK)));
  hmin = h;

  switch(lf->mi[MACRI])
  { case ACP:
      local_df(lf,des,t);
      mcp = acri(des->llk,t[0],t[2],lf->dp[DADP]);
      return(lf_status);
    case AKAT:
      local_df(lf,des,t);
      clo = des->cf[0]-lf->dp[DADP]*t[5];
      cup = des->cf[0]+lf->dp[DADP]*t[5];
      return(lf_status);
    case AMDI:
      mcp = mmse(lf,des);
      return(lf_status);
    case AOK: return(lf_status);
  }
  ERROR(("aband1: unknown criterion"));
  return(LF_ERR);
}

/*
  aband2 increases the initial bandwidth until lack of fit results,
  or the fit is close to a global fit. Increase h by 1+0.3/d at
  each iteration.
*/

double aband2(des,lf,h0)
design *des;
lfit   *lf;
double h0;
{ double t[6], h, h1, nu1, cp, ncp, tlo, tup;
  INT d, inc, n, p, done;
  d = lf->mi[MDIM]; n = lf->mi[MN]; p = lf->mi[MP];
  h1 = h = h0;
  done = 0; nu1 = 0.0;
  inc = 0; ncp = 0.0;
  while ((!done) & (nu1<(n-p)*0.95))
  { h = nbhd(lf,des,0,(1+0.3/d)*h,1);
    if (locfit(lf,des,h,1)>0) WARN(("aband2: failed fit"));
    local_df(lf,des,t);
    nu1 = t[0]-t[2]; /* tr(A) */
    switch(lf->mi[MACRI])
    { case AKAT:
        tlo = des->cf[0]-lf->dp[DADP]*t[5];
        tup = des->cf[0]+lf->dp[DADP]*t[5];
/* printf("h %8.5f  tlo %8.5f  tup %8.5f\n",h,tlo,tup); */
        done = ((tlo>cup) | (tup<clo));
        if (!done)
        { clo = MAX(clo,tlo);
          cup = MIN(cup,tup);
          h1 = h;
        }
        break;
      case ACP:
        cp = acri(des->llk,t[0],t[2],lf->dp[DADP]);
/* printf("h %8.5f  lk %8.5f  t0 %8.5f  t2 %8.5f  cp %8.5f\n",h,des->llk,t[0],t[2],cp); */
        if (cp<mcp) { mcp = cp; h1 = h; }
        if (cp>=ncp) inc++; else inc = 0;
        ncp = cp;
        done = (inc>=10) | ((inc>=3) & ((t[0]-t[2])>=10) & (cp>1.5*mcp));
        break;
      case AMDI:
        cp = mmse(lf,des);
        if (cp<mcp) { mcp = cp; h1 = h; }
        if (cp>ncp) inc++; else inc = 0;
        ncp = cp;
        done = (inc>=3);
        break;
    }
  }
  return(h1);
}

/*
  aband3 does a finer search around best h so far. Try
  h*(1-0.2/d), h/(1-0.1/d), h*(1+0.1/d), h*(1+0.2/d)
*/
double aband3(des,lf,h0)
design *des;
lfit   *lf;
double h0;
{ double t[6], h, h1, cp, tlo, tup;
  INT i, i0, d, n;
  d = lf->mi[MDIM]; n = lf->mi[MN];

  h1 = h0;
  i0 = (lf->mi[MACRI]==AKAT) ? 1 : -2;
  if (h0==hmin) i0 = 1;
  for (i=i0; i<=2; i++)
  { if (i==0) i++;
    h = h0*(1+0.1*i/lf->mi[MDIM]);
    h = nbhd(lf,des,0,h,1);
    if (locfit(lf,des,h,1)>0) WARN(("aband3: failed fit"));
    local_df(lf,des,t);
    switch (lf->mi[MACRI])
    { case AKAT:
        tlo = des->cf[0]-lf->dp[DADP]*t[5];
        tup = des->cf[0]+lf->dp[DADP]*t[5];
        if ((tlo>cup) | (tup<clo)) /* done */
          i = 2;
        else
        { h1 = h;
          clo = MAX(clo,tlo);
          cup = MIN(cup,tup);
        }
        break;
      case ACP:
        cp = acri(des->llk,t[0],t[2],lf->dp[DADP]);
        if (cp<mcp) { mcp = cp; h1 = h; }
        else
        { if (i>0) i = 2; }
        break;
      case AMDI:
        cp = mmse(lf,des);
        if (cp<mcp) { mcp = cp; h1 = h; }
        else
        { if (i>0) i = 2; }
    }
  }
  return(h1);
}
