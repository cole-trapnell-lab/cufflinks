/*
 *   Copyright (c) 1996-2001 Jiayang Sun, Catherine Loader.
 *   This file is used by the simultaneous confidence band
 *   additions to Locfit.
 *
 */

#include "local.h"
extern int cvi;
static double scb_crit, *x, c[10], kap[5], kaq[5], max_p2;
static int type, side;

double covar_par(lf,des,x1,x2)
lfit *lf;
design *des;
double x1, x2;
{ double *v1, *v2, *wk;
  paramcomp *pc;
  int i, j, p, ispar;

  v1 = des->f1; v2 = des->ss; wk = des->oc;
  ispar = (lf->mi[MKER]==WPARM) && (hasparcomp(lf));
  p = lf->mi[MP];

/*  for parametric models, the covariance is
 *  A(x1)^T (X^T W V X)^{-1} A(x2)
 *  which we can find easily from the parametric component.
 */
  if (ispar)
  { pc = &lf->pc;
    fitfun(lf,&x1,pc->xbar,v1,NULL,0);
    fitfun(lf,&x2,pc->xbar,v2,NULL,0);
    jacob_hsolve(&lf->pc.xtwx,v1);
    jacob_hsolve(&lf->pc.xtwx,v2);
  }

/*  for non-parametric models, we must use the cholseky decomposition
 *  of M2 = X^T W^2 V X. Courtesy of comp_vari, we already have
 *  des->P = M2^{1/2} M1^{-1}.
 */
  if (!ispar)
  { fitfun(lf,&x1,des->xev,wk,NULL,0);
    for (i=0; i<p; i++)
    { v1[i] = 0;
      for (j=0; j<p; j++) v1[i] += des->P[i*p+j]*wk[j];
    }
    fitfun(lf,&x2,des->xev,wk,NULL,0);
    for (i=0; i<p; i++)
    { v2[i] = 0;
      for (j=0; j<p; j++) v2[i] += des->P[i*p+j]*wk[j];
    }
  }

  return(innerprod(v1,v2,p));
}

void cumulant(lf,des,sd)
lfit *lf;
design *des;
double sd;
{ double b2i, b3i, b3j, b4i;
  double ss, si, sj, uii, uij, ujj, k1;
  INT ii, i, j, jj, *mi;
  for (i=1; i<10; i++) c[i] = 0.0;
  k1 = 0;
  mi = lf->mi;

  /* ss = sd*sd; */
  ss = covar_par(lf,des,des->xev[0],des->xev[0]);

/*
 * this isn't valid for nonparametric models. At a minimum,
 * the sums would have to include weights. Still have to work
 * out the right way.
 */
  for (i=0; i<mi[MN]; i++)
  { ii = des->ind[i];
    b2i = b2(des->th[i],mi[MTG],prwt(lf,ii));
    b3i = b3(des->th[i],mi[MTG],prwt(lf,ii));
    b4i = b4(des->th[i],mi[MTG],prwt(lf,ii));
    si = covar_par(lf,des,des->xev[0],datum(lf,0,ii));
    uii= covar_par(lf,des,datum(lf,0,ii),datum(lf,0,ii));
    if (lf_error) return;

    c[2] += b4i*si*si*uii;
    c[6] += b4i*si*si*si*si;
    c[7] += b3i*si*uii;
    c[8] += b3i*si*si*si;
    /* c[9] += b2i*si*si*si*si;
       c[9] += b2i*b2i*si*si*si*si; */
    k1 += b3i*si*(si*si/ss-uii);

    /* i=j components */
    c[1] += b3i*b3i*si*si*uii*uii;
    c[3] += b3i*b3i*si*si*si*si*uii;
    c[4] += b3i*b3i*si*si*uii*uii;

    for (j=i+1; j<mi[MN]; j++)
    { jj = des->ind[j];
      b3j = b3(des->th[j],mi[MTG],prwt(lf,jj));
      sj = covar_par(lf,des,des->xev[0],datum(lf,0,jj));
      uij= covar_par(lf,des,datum(lf,0,ii),datum(lf,0,jj));
      ujj= covar_par(lf,des,datum(lf,0,jj),datum(lf,0,jj));

      c[1] += 2*b3i*b3j*si*sj*uij*uij;
      c[3] += 2*b3i*b3j*si*si*sj*sj*uij;
      c[4] += b3i*b3j*uij*(si*si*ujj+sj*sj*uii);
      if (lf_error) return;
    }
  }
  c[5] = c[1];
  c[7] = c[7]*c[8];
  c[8] = c[8]*c[8];

  c[1] /= ss; c[2] /= ss; c[3] /= ss*ss; c[4] /= ss;
  c[5] /= ss; c[6] /= ss*ss; c[7] /= ss*ss;
  c[8] /= ss*ss*ss; c[9] /= ss*ss;

/* constants used in p(x,z) computation */
  kap[1] = k1/(2*sqrt(ss));
  kap[2] = 1 + 0.5*(c[1]-c[2]+c[4]-c[7]) - 3*c[3] + c[6] + 1.75*c[8];
  kap[4] = -9*c[3] + 3*c[6] + 6*c[8] + 3*c[9];

/* constants used in q(x,u) computation */
  kaq[2] = c[3] - 1.5*c[8] - c[5] - c[4] + 0.5*c[7] + c[6] - c[2];
  kaq[4] = -3*c[3] - 6*c[4] - 6*c[5] + 3*c[6] + 3*c[7] - 3*c[8] + 3*c[9];
}

/* q2(u) := u+q2(x,u) in paper */
double q2(u)
double u;
{ return(u-u*(36.0*kaq[2] + 3*kaq[4]*(u*u-3) + c[8]*((u*u-10)*u*u+15))/72.0);
}

/*  p2(u) := p2(x,u) in paper */
double p2(u)
double u;
{ return( -u*( 36*(kap[2]-1+kap[1]*kap[1])
     + 3*(kap[4]+4*kap[1]*sqrt(kap[3]))*(u*u-3)
     + c[8]*((u*u-10)*u*u+15) ) / 72 );
}

void procvscb2(des,lf,v)
design *des;
lfit *lf;
INT v;
{ double thhat, sd, *lo, *hi, u;
  int err, tmp;
  x = des->xev = evpt(lf,v);
  tmp = lf->mi[MPC];
  if ((lf->mi[MKER]==WPARM) && (hasparcomp(lf)))
  { lf->coef[v] = thhat = addparcomp(lf,des->xev,PCOEF);
    lf->nlx[v] = sd = addparcomp(lf,des->xev,PNLX);
  }
  else
  { lf->mi[MPC] = 0;
    procv(des,lf,v);
    thhat = lf->coef[v];
    sd = lf->nlx[v];
  }
  if (type >= 2)
  { if (lf->mi[MKER] != WPARM)
      WARN(("nonparametric fit; correction is invalid"));
    cumulant(lf,des,sd);
  }
  lf->mi[MPC] = tmp;
  lo = vdptr(lf->L);
  hi = &lo[lf->nvm];
  switch(type)
  { case 0:
    case 1: /* basic scr */
      lo[v] = thhat - scb_crit * sd;
      hi[v] = thhat + scb_crit * sd;
      return;
    case 2: /* centered scr */
      lo[v] = thhat - kap[1]*sd - scb_crit*sd*sqrt(kap[2]);
      hi[v] = thhat - kap[1]*sd + scb_crit*sd*sqrt(kap[2]);
      return;
    case 3: /* corrected 2 */
      u = solve_secant(q2,scb_crit,0.0,2*scb_crit,0.000001,BDF_NONE,&err);
      lo[v] = thhat - u*sd;
      hi[v] = thhat + u*sd;
      return;
    case 4: /* corrected 2' */
      u = fabs(p2(scb_crit));
      max_p2 = MAX(max_p2,u);
      lo[v] = thhat;
      hi[v] = thhat;
      return;
  }
  ERROR(("procvscb2: invalid type"));
}

void scb(des,lf)
design *des;
lfit *lf;
{ double kap[10], *lo, *hi;
  INT i, *mi, nterms;
  mi = lf->mi;
  mi[MP] = calcp(mi,mi[MDEG]);
  type = mi[MGETH] - 70;
  deschk(des,mi[MN],mi[MP]);
  des->pref = 0;
  cvi = -1; /* inhibit cross validation */
  mi[MLINK] = defaultlink(mi[MLINK],mi[MTG]);

  if (type==0)
  { kap[0] = 1;
    scb_crit = critval(kap,1,0,0.05,10,2,0.0);
  }
  else
  { compparcomp(des,lf,0);
    nterms = constants(des,lf,kap);
    scb_crit = critval(kap,nterms,mi[MDIM],0.05,10,2,0.0);
  }

  max_p2 = 0.0;
  startlf(des,lf,procvscb2,0);
  if (type==4)
  { lo = vdptr(lf->L);
    hi = &lo[lf->nvm];
    for (i=0; i<lf->nv; i++)
    {
      lo[i] -= (scb_crit-max_p2)*lf->nlx[i];
      hi[i] += (scb_crit-max_p2)*lf->nlx[i];
    }
  }
}

#ifdef CVERSION
extern lfit lf;
extern design des;
extern vari *aru;

lfit *lf_sim;
design *des_sim;

double scbsim_fun(x)
double x;
{ double y;
  evptx(lf_sim,0,0) = x;
  procv(des_sim,lf_sim,0);

  if (type>=2)
  { if (lf_sim->mi[MKER] != WPARM)
      WARN(("nonparametric fit; correction is invalid"));
    cumulant(lf_sim,des_sim,lf_sim->nlx[0]);
  }

  y = lf_link(dareval(aru,0,&x),lf_sim->mi[MLINK]);
  y = (lf_sim->coef[0] - y) / lf_sim->nlx[0];

  switch(type)
  {
    case 2:
      y = (y-kap[1]) / sqrt(kap[2]);
      break;
    case 3:
      y = (y-kap[1])/sqrt(kap[2]);
      y = (y>0) ? y+q2(y) : y - q2(y);
      break;
  }

  switch(side)
  { case -1: return(-y);
    case  1: return(y);
    default: return(fabs(y));
  }
}

static double max;

void do_scbsim(des,lf)
design *des;
lfit *lf;
{ double y;
  int err;

  lf_sim = lf;
  des_sim = des;

  trchck(lf,1,1,lf->mi[MDIM],lf->mi[MP],1);
  y = max_quad(scbsim_fun,lf->fl[0],lf->fl[1],10,0.00001,&err,'y');
  max = y;
}

void scbsim(lf,des)
lfit *lf;
design *des;
{ double kap[5];
  int nterms;

  lf->mi[MEV] = 100;
  startlf(des,lf,scbsim_fun,1);

  nterms = constants(des,lf,kap);
  printf("xmx: %10.6f  max: %10.6f  k0 %10.6f %10.6f  pr %10.6f\n",0.0,max,kap[0],kap[1],tailp(max,kap,nterms,lf->mi[MDIM],0.0));
}

void cscbsim(v)
vari *v;
{ int i;
  side = 0; type = 1;
  fitoptions(&lf,v,0);

  i = getarg(v,"mean",1);
  if (i==0)
  { WARN(("cscbsim: no mean function; setting = 0"));
    aru = arbuild("0",0,0,NULL,0,1);
  }
  else
  { aru = arbuild(argval(v,i),0,strlen(argval(v,i))-1,NULL,0,1);
    setvarname(aru,"_aru");
  }

  i = getarg(v,"corr",1);
  if (i>0) type = getlogic(v,i);
  if (lf_error) return;

  i = getarg(v,"side",1);
  if (i>0) sscanf(argval(v,i),"%d",&side);
  if (lf_error) return;

  scbsim(&lf,&des);
}
#endif
