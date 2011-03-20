/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *   This file contains functions for constructing and
 *   interpolating the adaptive tree structure. This is
 *   the default evaluation structure used by Locfit.
 */

#include "local.h"

/*
  Guess the number of fitting points.
  Needs improving!
*/
void atree_guessnv(nvm,ncm,vc,dp,mi)
double *dp;
INT *nvm, *ncm, *vc, *mi;
{ double a0, cu, ifl;
  int i, nv, nc;

  *ncm = 1<<30; *nvm = 1<<30;
  *vc = 1 << mi[MDIM];

  if (dp[DALP]>0)
  { a0 = (dp[DALP] > 1) ? 1 : 1/dp[DALP];
    if (dp[DCUT]<0.01)
    { WARN(("guessnv: cut too small."));
      dp[DCUT] = 0.01;
    }
    cu = 1;
    for (i=0; i<mi[MDIM]; i++) cu *= MIN(1.0,dp[DCUT]);
    nv = (int)((5*a0/cu+1)**vc);  /* this allows 10*a0/cu splits */
    nc = (int)(10*a0/cu+1);      /* and 10*a0/cu cells */
    if (nv<*nvm) *nvm = nv;
    if (nc<*ncm) *ncm = nc;
  }

  if (*nvm == 1<<30) /* by default, allow 100 splits */
  { *nvm = 102**vc;
    *ncm = 201;
  }

  /* inflation based on mi[MK] */
  ifl = mi[MK]/100.0;
  *nvm = (int)(ifl**nvm);
  *ncm = (int)(ifl**ncm);
  
}

/*
  Determine whether a cell in the tree needs splitting.
  If so, return the split variable (0..d-1).
  Otherwise, return -1.
*/
INT atree_split(lf,ce,le,ll,ur)
lfit *lf;
INT *ce;
double *le, *ll, *ur;
{ INT d, vc, i, is;
  double h, hmin, score[MXDIM];
    memset(score, 0, sizeof(score));
  d = lf->mi[MDIM]; vc = 1<<d;

  hmin = 0.0;
  for (i=0; i<vc; i++)
  { h = lf->h[ce[i]];
    if ((h>0) && ((hmin==0)|(h<hmin))) hmin = h;
  }

  is = 0;
  for (i=0; i<d; i++)
  { le[i] = (ur[i]-ll[i])/lf->sca[i];
    if ((lf->sty[i]==STCPAR) || (hmin==0))
      score[i] = 2*(ur[i]-ll[i])/(lf->fl[i+d]-lf->fl[i]);
    else
      score[i] = le[i]/hmin;
    if (score[i]>score[is]) is = i;
  }
  if (lf->dp[DCUT]<score[is]) return(is);
  return(-1);
}

/*
  recursively grow the tree structure, begining with the parent cell.
*/
void atree_grow(des,lf,ce,ct,term,ll,ur)
design *des;
lfit *lf;
INT *ce, *ct, *term;
double *ll, *ur;
{ INT i, i0, i1, d, vc, ns, tk, nce[1024], pv;
  double le[MXDIM], z;
  d = lf->mi[MDIM]; vc = 1<<d;

  /* does this cell need splitting?
     If not, wrap up and return.
  */
  ns = atree_split(lf,ce,le,ll,ur);
  if (ns==-1)
  { if (ct != NULL) /* reconstructing terminal cells */
    { for (i=0; i<vc; i++) term[*ct*vc+i] = ce[i];
      (*ct)++;
    }
    return;
  }

  /* split the cell at the midpoint on side ns */
  tk = 1<<ns;
  for (i=0; i<vc; i++)
  { if ((i&tk)==0) nce[i] = ce[i];
    else
    { i0 = ce[i];
      i1 = ce[i-tk];
      pv = (lf->sty[i]!=STCPAR) &&
           (le[ns] < (lf->dp[DCUT]*MIN(lf->h[i0],lf->h[i1])));
      nce[i] = newsplit(des,lf,i0,i1,pv);
      if (lf_error) return;
    }
  }
  z = ur[ns]; ur[ns] = (z+ll[ns])/2;
  atree_grow(des,lf,nce,ct,term,ll,ur);
  if (lf_error) return;
  ur[ns] = z;
  for (i=0; i<vc; i++)
    nce[i] = ((i&tk)== 0) ? nce[i+tk] : ce[i];
  z = ll[ns]; ll[ns] = (z+ur[ns])/2;
  atree_grow(des,lf,nce,ct,term,ll,ur);
  ll[ns] = z;
}

void atree_start(des,lf)
design *des;
lfit *lf;
{ INT i, j, vc, d, ncm, nvm, k;
  double ll[MXDIM], ur[MXDIM];

  d = lf->mi[MDIM];
  atree_guessnv(&nvm,&ncm,&vc,lf->dp,lf->mi);
  trchck(lf,nvm,ncm,d,des->p,vc);

  /* Set the lower left, upper right limits. */
  for (j=0; j<d; j++)
  { ll[j] = lf->fl[j];
    ur[j] = lf->fl[j+d];
  }

  /* Set the initial cell; fit at the vertices. */
  for (i=0; i<vc; i++)
  { j = i;
    for (k=0; k<d; ++k)
    { evptx(lf,i,k) = (j%2) ? ur[k] : ll[k];
      j >>= 1;
    }
    lf->ce[i] = i;
    des->vfun(des,lf,i);
    if (lf_error) return;
    lf->s[i] = 0;
  }
  lf->nv = vc;

  /* build the tree */
  atree_grow(des,lf,lf->ce,NULL,NULL,ll,ur);
  lf->nce = 1;
}

double atree_int(tr,x,what)
lfit *tr;
double *x;
INT what;
{ double vv[64][64], *ll, *ur, h, xx[MXDIM];
  INT d, i, lo, tk, ns, nv, nc, vc, ce[64];
  d = tr->mi[MDIM];
  vc = 1<<tr->mi[MDIM];
  for (i=0; i<vc; i++)
  { setzero(vv[i],vc);
    nc = exvval(tr,vv[i],i,d,what,1);
    ce[i] = tr->ce[i];
  }
  ns = 0;
  while(ns!=-1)
  { ll = evpt(tr,ce[0]); ur = evpt(tr,ce[vc-1]);
    ns = atree_split(tr,ce,xx,ll,ur);
    if (ns!=-1)
    { tk = 1<<ns;
      h = ur[ns]-ll[ns];
      lo = (2*(x[ns]-ll[ns])) < h;
      for (i=0; i<vc; i++) if ((tk&i)==0)
      { nv = newsplit((design *)NULL,tr,ce[i],ce[i+tk],0);
        if (lf_error) return(0.0);
        if (lo)
        { ce[i+tk] = nv;
          if (tr->s[nv]) exvvalpv(vv[i+tk],vv[i],vv[i+tk],d,ns,h,nc);
                    else exvval(tr,vv[i+tk],nv,d,what,1);
        }
        else
        { ce[i] = nv;
          if (tr->s[nv]) exvvalpv(vv[i],vv[i],vv[i+tk],d,ns,h,nc);
                    else exvval(tr,vv[i],nv,d,what,1);
      } }
  } }
  ll = evpt(tr,ce[0]); ur = evpt(tr,ce[vc-1]);
  return(rectcell_interp(x,vdptr(tr->xxev),vv,ll,ur,d,nc));
}
