/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
  Routines for computing weight diagrams.
  wdiag(lf,des,lx,deg,ty,exp)
  Must locfit() first, unless ker==WPARM and has par. comp.
 
  also, vertex processing function procvhatm().

  cwdiag() entry for CLOCFIT.
 */

#include "local.h"

static double *wd;
extern double robscale;
extern void unitvec();

void nnresproj(lf,des,u,m,p,mi)
lfit *lf;
design *des;
double *u;
INT m, p, *mi;
{ INT i, j;
  double link[LLEN];
  setzero(des->f1,p);
  for (j=0; j<m; j++)
  { stdlinks(link,lf,des->ind[j],des->th[j],robscale);
    for (i=0; i<p; i++) des->f1[i] += link[ZDDLL]*d_xij(des,j,i)*u[j];
  }
  jacob_solve(&des->xtwx,des->f1);
  for (i=0; i<m; i++)
    u[i] -= innerprod(des->f1,d_xi(des,i),p)*des->w[i];
}

void wdexpand(l,n,ind,m)
double *l;
INT n, m, *ind;
{ INT i, j, t;
  double z;
  for (j=m; j<n; j++) { l[j] = 0.0; ind[j] = -1; }
  j = m-1;
  while (j>=0)
  { if (ind[j]==j) j--;
    else
    { i = ind[j];
      z = l[j]; l[j] = l[i]; l[i] = z;
      t = ind[j]; ind[j] = ind[i]; ind[i] = t;
      if (ind[j]==-1) j--;
    }
  }

/*  for (i=n-1; i>=0; i--)
  { l[i] = ((j>=0) && (ind[j]==i)) ? l[j--] : 0.0; } */
}

INT wdiagp(lf,des,lx,deg,ty,exp)
lfit *lf;
design *des;
double *lx;
INT deg, ty, exp;
{ INT i, j, p, *mi, *deriv, nd;
  double *l1;
  mi = lf->mi; p = des->p;
  deriv = lf->deriv; nd = lf->nd;
  fitfun(lf,des->xev,lf->pc.xbar,des->f1,deriv,nd);
  if (exp)
  { jacob_solve(&lf->pc.xtwx,des->f1);
    for (i=0; i<mi[MN]; i++)
      lx[i] = innerprod(des->f1,d_xi(des,i),p);
    return(mi[MN]);
  }
  jacob_hsolve(&lf->pc.xtwx,des->f1);
  for (i=0; i<p; i++) lx[i] = des->f1[i];

  if (deg>=1)
    for (i=0; i<mi[MDIM]; i++)
    { deriv[nd] = i;
      l1 = &lx[(i+1)*p];
      fitfun(lf,des->xev,lf->pc.xbar,l1,deriv,nd+1);
      jacob_hsolve(&lf->pc.xtwx,l1);
    }

  if (deg>=2)
    for (i=0; i<mi[MDIM]; i++)
    { deriv[nd] = i;
      for (j=0; j<mi[MDIM]; j++)
      { deriv[nd+1] = j;
        l1 = &lx[(i*mi[MDIM]+j+mi[MDIM]+1)*p];
        fitfun(lf,des->xev,lf->pc.xbar,l1,deriv,nd+2);
        jacob_hsolve(&lf->pc.xtwx,l1);
    } }
  return(p);
}

INT wdiag(lf,des,lx,deg,ty,exp)
lfit *lf;
design *des;
double *lx;
INT deg, ty, exp;
/* deg=0: l(x) only.
   deg=1: l(x), l'(x) (approx/exact ? mi[MDC] );
   deg=2: l(x), l'(x), l''(x);
   ty = 1: e1 (X^T WVX)^{-1} X^T W        -- hat matrix
   ty = 2: e1 (X^T WVX)^{-1} X^T WV^{1/2} -- scb's
*/
{ double w, *X, *lxd, *lxdd, wdd, wdw, *ulx, link[LLEN], h;
  double dfx[MXDIM], hs[MXDIM];
  INT i, ii, j, k, l, m, d, p, *mi, *deriv, nd;
  if ((lf->mi[MKER]==WPARM) && (hasparcomp(lf)))
    return(wdiagp(lf,des,lx,deg,ty,exp));
  mi = lf->mi; h = des->h;
  deriv = lf->deriv; nd = lf->nd;
  wd = des->wd;
  d = mi[MDIM]; p = des->p; X = d_x(des);
  ulx = des->res;
  m = des->n;
  for (i=0; i<d; i++) hs[i] = h*lf->sca[i];
  if (deg>0)
  { lxd = &lx[m];
    setzero(lxd,m*d);
    if (deg>1)
    { lxdd = &lxd[d*m];
      setzero(lxdd,m*d*d);
  } }
  if (nd>0) fitfun(lf,des->xev,des->xev,des->f1,deriv,nd); /* c(0) */
    else unitvec(des->f1,0,p);
  jacob_solve(&des->xtwx,des->f1);   /* c(0) (X^TWX)^{-1} */
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    lx[i] = innerprod(des->f1,&X[i*p],p); /* c(0)(XTWX)^{-1}X^T */
    if ((deg>0) && (mi[MDC]))
    { wd[i] = Wd(des->di[ii]/h,mi[MKER]);
      for (j=0; j<d; j++)
      { dfx[j] = datum(lf,j,ii)-des->xev[j];
        lxd[j*m+i] = lx[i]*des->w[i]*weightd(dfx[j],lf->sca[j],
          d,mi[MKER],mi[MKT],h,lf->sty[j],des->di[ii]);
             /* c(0) (XTWX)^{-1}XTW' */
      }
      if (deg>1)
      { wdd = Wdd(des->di[ii]/h,mi[MKER]);
        for (j=0; j<d; j++)
          for (k=0; k<d; k++)
          { w = (des->di[ii]==0) ? 0 : h/des->di[ii];
            w = wdd * (des->xev[k]-datum(lf,k,ii)) * (des->xev[j]-datum(lf,j,ii))
                  * w*w / (hs[k]*hs[k]*hs[j]*hs[j]);
            if (j==k) w += wd[i]/(hs[j]*hs[j]);
            lxdd[(j*d+k)*m+i] = lx[i]*w;
              /* c(0)(XTWX)^{-1}XTW'' */
          }
      }
    }
    lx[i] *= des->w[i];
  }
  if ((deg==2) && (mi[MDC]))
  { for (i=0; i<d; i++)
    { deriv[nd] = i;
      fitfun(lf,des->xev,des->xev,des->f1,deriv,nd+1);
      for (k=0; k<m; k++)
      { stdlinks(link,lf,des->ind[k],des->th[k],robscale);
        for (j=0; j<p; j++)
          des->f1[j] -= link[ZDDLL]*lxd[i*m+k]*X[k*p+j];
        /* c'(x)-c(x)(XTWX)^{-1}XTW'X */
      }
      jacob_solve(&des->xtwx,des->f1); /* (...)(XTWX)^{-1} */
      for (j=0; j<m; j++)
        ulx[j] = innerprod(des->f1,&X[j*p],p); /* (...)XT */
      for (j=0; j<d; j++)
        for (k=0; k<m; k++)
        { ii = des->ind[k];
          dfx[j] = datum(lf,j,ii)-des->xev[j];
          wdw = des->w[k]*weightd(dfx[j],lf->sca[j],d,mi[MKER],mi[MKT],h,
            lf->sty[j],des->di[ii]);
          lxdd[(i*d+j)*m+k] += ulx[k]*wdw;
          lxdd[(j*d+i)*m+k] += ulx[k]*wdw;
        }
        /* + 2(c'-c(XTWX)^{-1}XTW'X)(XTWX)^{-1}XTW' */
    }
    for (j=0; j<d*d; j++) nnresproj(lf,des,&lxdd[j*m],m,p,mi);
        /* * (I-X(XTWX)^{-1} XTW */
  }
  if (deg>0)
  { if (mi[MDC]) for (j=0; j<d; j++) nnresproj(lf,des,&lxd[j*m],m,p,mi);
             /* c(0)(XTWX)^{-1}XTW'(I-X(XTWX)^{-1}XTW) */
    for (i=0; i<d; i++)
    { deriv[nd]=i;
      fitfun(lf,des->xev,des->xev,des->f1,deriv,nd+1);
      jacob_solve(&des->xtwx,des->f1);
      for (k=0; k<m; k++)
        for (l=0; l<p; l++)
          lxd[i*m+k] += des->f1[l]*X[k*p+l]*des->w[k];
            /* add c'(0)(XTWX)^{-1}XTW */
    }
  }
  if (deg==2)
  { for (i=0; i<d; i++)
    { deriv[nd]=i;
      for (j=0; j<d; j++)
      { deriv[nd+1]=j;
        fitfun(lf,des->xev,des->xev,des->f1,deriv,nd+2);
        jacob_solve(&des->xtwx,des->f1);
        for (k=0; k<m; k++)
          for (l=0; l<p; l++)
            lxdd[(i*d+j)*m+k] += des->f1[l]*X[k*p+l]*des->w[k];
        /* + c''(x)(XTWX)^{-1}XTW */
      }
    }
  }
  k = 1+d*(deg>0)+d*d*(deg==2);

  if (exp) wdexpand(lx,mi[MN],des->ind,m);
 
  if (ty==1) return(m);
  for (i=0; i<m; i++)
  { stdlinks(link,lf,des->ind[i],des->th[i],robscale);
    link[ZDDLL] = sqrt(fabs(link[ZDDLL]));
    for (j=0; j<k; j++) lx[j*m+i] *= link[ZDDLL];
  }
  return(m);
}

INT procvhatm(des,lf,v)
design *des;
lfit *lf;
INT v;
{ INT k = 0, n;
  double *l;
  n = (ident==0) ? lf->mi[MN] : des->p;
  if ((lf->mi[MKER]!=WPARM) | (!hasparcomp(lf))) k = procvraw(des,lf,v);
  l = (double *)viptr(lf->L,v*n);
  wdiag(lf,des,l,(INT)0,(INT)1,1);
  return(k);
}

#ifdef CVERSION
extern lfit lf;
extern design des;

void cwdiag(v)
vari *v;
{ INT i;
  vari *ve, *vr;
  i = getarg(v,"ev",0);
  if (i==0) ERROR(("wdiag: no ev point"));
  fitoptions(&lf,v,0);
  if (lf_error) return;
  ve = varith(argval(v,i),"wdev",STPLOTVAR);
  vr = createvar("wdres",STHIDDEN,lf.mi[MN],VDOUBLE);
  lf.xxev = ve;
  lf.nv = lf.nvm = 1;
  lf.L = vr;
  lf.mi[MEV] = EPRES;
  startlf(&des,&lf,procvhatm,1);
  deletevar(ve);
  saveresult(vr,argarg(v,0),STREGULAR);
}

#endif
