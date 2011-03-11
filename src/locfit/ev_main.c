/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

void kdtre_guessnv(int* nvm, int* ncm, int* vc, double* dp, int* mi);

INT cvi=-1;

/*
 * the C version has its own createvar, file src-c/vari.c.
 * For other versions, use the simple function here.
 */
#ifndef CVERSION
vari *createvar(name,type,n,mode)
varname name;
INT type, n, mode;
{ vari *v;
printf("creating variable\n");
  v = (vari *)malloc(sizeof(vari));
  vdptr(v) = (double *)calloc( vlength(v)=n ,
              (mode==VDOUBLE) ? sizeof(double) : sizeof(INT));
  return(v);
}
#endif

/*
 *  guessnv is used to guess the number of vertices etc in the
 *  evaluation structure.
 */
void guessnv(nvm,ncm,vc,dp,mi)
double *dp;
int *nvm, *ncm, *vc, *mi;
{ switch(mi[MEV])
  { case ETREE:
      atree_guessnv(nvm,ncm,vc,dp,mi);
      return;
    case EPHULL:
      *nvm = *ncm = mi[MK]*mi[MDIM];
      *vc = mi[MDIM]+1;
      return;
    case EDATA:
    case ECROS:
      *nvm = mi[MN];
      *ncm = *vc = 0;
      return;
    case EKDTR:
    case EKDCE:
      kdtre_guessnv(nvm,ncm,vc,dp,mi);
      return;
    case EXBAR:
    case ENONE:
      *nvm = 1;
      *ncm = *vc = 0;
      return;
    case EPRES:
      /* preset -- leave everything unchanged. */
      return;
    default:
      ERROR(("guessnv: I don't know this evaluation structure."));
  }
  return;
}

/*
 * trchck checks the working space on the lfit structure 
 * has space for nvm vertices and ncm cells.
 */
int lfit_reqd(d,nvm,ncm)
INT d, nvm, ncm;
{ return(nvm*(3*d+8)+ncm);
}
int lfit_reqi(nvm,ncm,vc)
INT nvm, ncm, vc;
{ return(ncm*vc+3*MAX(ncm,nvm));
}

void trchck(tr,nvm,ncm,d,p,vc)
lfit *tr;
INT nvm, ncm, d, p, vc;
{ INT rw, *k;
  double *z;

  tr->xxev = checkvarlen(tr->xxev,d*nvm,"_lfxev",VDOUBLE);

  rw = nvm*(3*d+8)+ncm;
  tr->tw = checkvarlen(tr->tw,lfit_reqd(d,nvm,ncm),"_lfwork",VDOUBLE);
  z = (double *)vdptr(tr->tw);
  tr->coef= z; z += nvm*(d+1);
  tr->nlx = z; z += nvm*(d+1);
  tr->t0  = z; z += nvm*(d+1);
  tr->lik = z; z += 3*nvm;
  tr->h   = z; z += nvm;
  tr->deg = z; z += nvm;
  tr->sv  = z; z += ncm;
  if (z != (double *)vdptr(tr->tw)+rw)
    WARN(("trchck: double assign problem"));

  rw = lfit_reqi(nvm,ncm,vc);
  tr->iw = checkvarlen(tr->iw,rw,"_lfiwork",VINT);
  k = (INT *)vdptr(tr->iw);
  tr->ce = k; k += vc*ncm;
  tr->s  = k; k += MAX(ncm,nvm);
  tr->lo = k; k += MAX(ncm,nvm);
  tr->hi = k; k += MAX(ncm,nvm);
  if (k != (INT *)vdptr(tr->iw)+rw)
    WARN(("trchck: int assign problem"));

  tr->nvm = nvm; tr->ncm = ncm; tr->mi[MDIM] = d; tr->mi[MP] = p; tr->vc = vc;
}

#ifdef CVERSION
void reassign(lf)
lfit *lf;
{ INT i, nvm, ncm, vc, d, k, p, *iw;
  double *tw, *ntw;
  setvarname(lf->tw,"__lfwork"); /* prevent overwrite */
  setvarname(lf->iw,"__lfiwork");
  nvm = lf->nvm; ncm = lf->ncm; vc = lf->vc;
  tw = (double *)vdptr(lf->tw);
  iw = (INT *)vdptr(lf->iw);
  d = lf->mi[MDIM];
  p = lf->mi[MP];
  trchck(lf,2*nvm,ncm,d,p,vc);
  ntw = vdptr(lf->tw);
/*
  xev is stored in blocks of d. other matrices by blocks on nvm
*/
  k = nvm*d;
  memcpy(vdptr(lf->xxev),tw,k*sizeof(double));
  tw += k; ntw += 2*k;
  for (i=0; i<2*p+2*d+6; i++)
  { memcpy(ntw,tw,nvm*sizeof(double));
    tw += nvm; ntw += 2*nvm;
  }
  k = ncm;       memcpy(lf->sv,tw,k*sizeof(double));   tw += k;

  k = vc*ncm;       memcpy(lf->ce,iw,k*sizeof(INT)); iw += k;
  k = MAX(ncm,nvm); memcpy(lf->s,iw,k*sizeof(INT));  iw += k;
  k = MAX(ncm,nvm); memcpy(lf->lo,iw,k*sizeof(INT)); iw += k;
  k = MAX(ncm,nvm); memcpy(lf->hi,iw,k*sizeof(INT)); iw += k;
  deletename("__lfwork");
  deletename("__lfiwork");
}
#endif

void dataf(des,lf)
design *des;
lfit *lf;
{ INT d, i, j, ncm, nv, vc;

  d = lf->mi[MDIM];
  guessnv(&nv,&ncm,&vc,lf->dp,lf->mi);
  trchck(lf,nv,0,d,des->p,0);

  for (i=0; i<nv; i++)
    for (j=0; j<d; j++) evptx(lf,i,j) = datum(lf,j,i);
  for (i=0; i<nv; i++)
  { des->vfun(des,lf,i);
    lf->s[i] = 0;
  }
  lf->nv = lf->nvm = nv; lf->nce = 0;
}

void xbarf(des,lf)
design *des;
lfit *lf;
{ int i, d, nvm, ncm, vc;
  d = lf->mi[MDIM];
  guessnv(&nvm,&ncm,&vc,lf->dp,lf->mi);
  trchck(lf,1,0,d,des->p,0);
  for (i=0; i<d; i++) evptx(lf,0,i) = lf->pc.xbar[i];
  des->vfun(des,lf,0);
  lf->s[0] = 0;
  lf->nv = 1; lf->nce = 0;
}

#ifndef GR
void preset(des,lf)
design *des;
lfit *lf;
{ INT i, nv;
  double *tmp;
  nv = lf->nvm;
  tmp = vdptr(lf->xxev);
  trchck(lf,nv,0,lf->mi[MDIM],des->p,0);
  lf->xxev->dpr = tmp;
  for (i=0; i<nv; i++)
  { 
    des->vfun(des,lf,i);
    lf->s[i] = 0;
  }
  lf->nv = nv; lf->nce = 0;
}
#endif

void crossf(des,lf)
design *des;
lfit *lf;
{ INT d, i, j, n, nv, ncm, vc;

  n = lf->mi[MN]; d = lf->mi[MDIM];
  guessnv(&nv,&ncm,&vc,lf->dp,lf->mi);
  trchck(lf,n,0,d,des->p,0);

  for (i=0; i<n; i++)
    for (j=0; j<d; j++) evptx(lf,i,j) = datum(lf,j,i);
  for (cvi=0; cvi<n; cvi++)
  { lf->s[cvi] = 0;
    des->vfun(des,lf,cvi);
  }
  cvi = -1;
  lf->nv = n; lf->nce = 0; lf->mi[MN] = n;
}

void gridf(des,tr)
design *des;
lfit *tr;
{ INT d, i, j, nv, u0, u1, z;
  nv = 1; d = tr->mi[MDIM];
  for (i=0; i<d; i++)
  { if (tr->mg[i]==0)
      tr->mg[i] = 2+(INT)((tr->fl[i+d]-tr->fl[i])/(tr->sca[i]*tr->dp[DCUT]));
    nv *= tr->mg[i];
  }
  trchck(tr,nv,0,d,des->p,1<<d);
  for (i=0; i<nv; i++)
  { z = i;
    for (j=0; j<d; j++)
    { u0 = z%tr->mg[j];
      u1 = tr->mg[j]-1-u0;
      evptx(tr,i,j) = (tr->mg[j]==1) ? tr->fl[j] :
                      (u1*tr->fl[j]+u0*tr->fl[j+d])/(tr->mg[j]-1);
      z = z/tr->mg[j];
    }
    tr->s[i] = 0;
    des->vfun(des,tr,i);
  }
  tr->nv = nv; tr->nce = 0;
}

/*
  add a new vertex at the midpoint of (x[i0],x[i1]).
  return the vertex number.
*/
INT newsplit(des,lf,i0,i1,pv)
design *des;
lfit *lf;
INT i0, i1, pv;
{ INT i, nv;

  /* first, check to see if the new point already exists */
  if (i0>i1) ISWAP(i0,i1);
  nv = lf->nv;
  for (i=i1+1; i<nv; i++)
    if ((lf->lo[i]==i0) && (lf->hi[i]==i1)) return(i);
  
  /* the point is new. Now check we have space for the new point. */
  if (nv==lf->nvm)
  {
#ifdef CVERSION
    reassign(lf);
#else
    ERROR(("newsplit: out of vertex space"));
    return(-1);
#endif
  }

  /* compute the new point, and evaluate the fit */
  lf->lo[nv] = i0;
  lf->hi[nv] = i1;
  for (i=0; i<lf->mi[MDIM]; i++)
    evptx(lf,nv,i) = (evptx(lf,i0,i)+evptx(lf,i1,i))/2;
  if (pv) /* pseudo vertex */
  { lf->h[nv] = (lf->h[i0]+lf->h[i1])/2;
    lf->s[nv] = 1; /* pseudo-vertex */
  }
  else /* real vertex */
  { des->vfun(des,lf,nv);
    lf->s[nv] = 0;
  }
  lf->nv++;

  return(nv);
}
