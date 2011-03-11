/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *   Routines for building and interpolating the kd tree.
 *   Initially, this started from the loess code.
 *
 *   Todo: EKDCE isn't working.
 */

#include "local.h"

void newcell();
static INT nterm;

void kdtre_guessnv(nvm,ncm,vc,dp,mi)
double *dp;
INT *nvm, *ncm, *vc, *mi;
{ int k;
  if (mi[MEV] == EKDTR)
  { nterm = (INT) (dp[DCUT]/4 * mi[MN] * MIN(dp[DALP],1.0) );
    k = 2*mi[MN]/nterm;
    *vc = 1<<mi[MDIM];
    *ncm = 2*k+1;
    *nvm = (k+2)**vc/2;
    return;
  }
  if (mi[MEV] == EKDCE)
  { nterm = (INT) (mi[MN] * dp[DALP]);
    *vc = 1;
    *nvm = 1+(INT)(2*mi[MN]/nterm);
    *ncm = 2**nvm+1;
    return;
  }
  *nvm = *ncm = *vc = 0;
  return;
}

/*
  Split x[pi[l..r]] into two equal sized sets.

  Let m=(l+r)/2.
  At return,
    x[pi[l..m]] < x[pi[m+1..r]].
    Assuming no ties:
      If l+r is odd, the sets have the same size.
      If l+r is even, the low set is larger by 1.
    If there are ties, all ties go in the low set.
*/      
int ksmall(l, r, m, x, pi)
INT l, r, m, *pi;
double *x;
{
  int il, ir, jl, jr;
  double t;


  while(l<r)
  { t = x[pi[m]];

    /*
      permute the observations so that
        x[pi[l..il]] < t <= x[pi[ir..r]].
    */
    ir = l; il = r;
    while (ir<il)
    { while ((ir<=r) && (x[pi[ir]] < t)) ir++;
      while ((il>=l) && (x[pi[il]]>= t)) il--;
      if (ir<il) ISWAP(pi[ir],pi[il]);
    }

    /*
      move  = t to the middle:
        x[pi[l..il]]  < t
        x[pi[jl..jr]] = t
        x[pi[ir..r]] > t
    */
    jl = ir; jr = r;
    while (ir<jr)
    { while ((ir<=r)  && (x[pi[ir]]== t)) ir++;
      while ((jr>=jl) && (x[pi[jr]] > t)) jr--;
      if (ir<jr) ISWAP(pi[ir],pi[jr]);
    }

    /*
      we're done if m is in the middle, jl <= m <= jr.
    */
    if ((jl<=m) & (jr>=m)) return(jr);

    /*
      update l or r.
    */
    if (m>=ir) l = ir;
    if (m<=il) r = il;
  }
  if (l==r) return(l);
  ERROR(("ksmall failure"));
  return(0);
}

INT terminal(tr,p,pi,fc,d,m,split_val)
lfit *tr;
INT p, *pi, d, fc, *m;
double *split_val;
{ INT i, k, lo, hi, split_var;
  double max, min, score, max_score, t;

  /*
    if there are fewer than fc points in the cell, this cell
    is terminal.
  */
  lo = tr->lo[p]; hi = tr->hi[p];
  if (hi-lo < fc) return(-1);

  /* determine the split variable */
  max_score = 0.0; split_var = 0;
  for (k=0; k<d; k++)
  { max = min = datum(tr, k, pi[lo]);
    for (i=lo+1; i<=hi; i++)
    { t = datum(tr,k,pi[i]);
      if (t<min) min = t;
      if (t>max) max = t;
    }
    score = (max-min) / tr->sca[k];
    if (score > max_score)
    { max_score = score;
      split_var = k;
    }
  }
  if (max_score==0) /* all points in the cell are equal */
    return(-1);

  *m = ksmall(lo,hi,(lo+hi)/2, dvari(tr,split_var), pi);
  *split_val = datum(tr, split_var, pi[*m]);

  if (*m==hi) /* all observations go lo */
    return(-1);
  return(split_var);
}

void kdtre_start(des,tr)
design *des;
lfit *tr;
{ INT i, j, vc, d, nc, nv, ncm, nvm, k, m, n, p, *pi;
  double sv;
  d = tr->mi[MDIM]; n = tr->mi[MN]; pi = des->ind;
  kdtre_guessnv(&nvm,&ncm,&vc,tr->dp,tr->mi);
  trchck(tr,nvm,ncm,d,des->p,vc);

  nv = 0;
  if (tr->mi[MEV] != EKDCE)
  { for (i=0; i<vc; i++)
    { j = i;
      for (k=0; k<d; ++k)
      { evptx(tr,i,k) = tr->fl[d*(j%2)+k];
        j >>= 1;
      }
    }
    nv = vc;
    for (j=0; j<vc; j++) tr->ce[j] = j;
  }

  for (i=0; i<n; i++) pi[i] = i;
  p = 0; nc = 1;
  tr->lo[p] = 0; tr->hi[p] = n-1;
  tr->s[p] = -1;
  while (p<nc)
  { k = terminal(tr,p,pi,nterm,d,&m,&sv);
    if (k>=0)
    {
      if ((ncm<nc+2) | (2*nvm<2*nv+vc))
      { WARN(("Insufficient space for full tree"));
        tr->nce = nc; tr->nv = nv;
        return;
      }

      /* new lo cell has obsn's tr->lo[p]..m */
      tr->lo[nc] = tr->lo[p];
      tr->hi[nc] = m;
      tr->s[nc] = -1;

      /* new hi cell has obsn's m+1..tr->hi[p] */
      tr->lo[nc+1] = m+1;
      tr->hi[nc+1] = tr->hi[p];
      tr->s[nc+1] = -1;

      /* cell p is split on variable k, value sv */
      tr->s[p] = k;
      tr->sv[p] = sv;
      tr->lo[p] = nc; tr->hi[p] = nc+1;

      nc=nc+2; i = nv;

      /* now compute the new vertices. */
      if (tr->mi[MEV] != EKDCE)
        newcell(&nv,vc,vdptr(tr->xxev), d, k, sv,
             &tr->ce[p*vc], &tr->ce[(nc-2)*vc], &tr->ce[(nc-1)*vc]);

    }
    else if (tr->mi[MEV]==EKDCE) /* new vertex at cell center */
    { sv = 0;
      for (i=0; i<d; i++) evptx(tr,nv,i) = 0;
      for (j=tr->lo[p]; j<=tr->hi[p]; j++)
      { sv += prwt(tr,pi[j]);
        for (i=0; i<d; i++)
          evptx(tr,nv,i) += datum(tr,i,pi[j])*prwt(tr,pi[j]);
      }
      for (i=0; i<d; i++) evptx(tr,nv,i) /= sv;
      tr->mi[MN] = tr->hi[p]-tr->lo[p]+1;
      des->ind = &pi[tr->lo[p]];
      des->vfun(des,tr,nv);
      tr->mi[MN] = n; des->ind = pi;
      nv++;
    }
    p++;
  }

  /* We've built the tree. Now do the fitting. */
  if (tr->mi[MEV]==EKDTR)
    for (i=0; i<nv; i++) des->vfun(des,tr,i);

  tr->nce = nc; tr->nv = nv;
  return;
}

void newcell(nv,vc,xev, d, k, split_val, cpar, clef, crig)
double *xev, split_val;
INT *nv, vc, d, k, *cpar, *clef, *crig;
{ INT i, ii, j, j2, tk, match;
  tk = 1<<k;
  for (i=0; i<vc; i++)
  { if ((i&tk) == 0)
    { for (j=0; j<d; j++) xev[*nv*d+j] = xev[d*cpar[i]+j];
      xev[*nv*d+k] = split_val;
      match = 0; j = vc; /* no matches in first vc points */
      while ((j<*nv) && (!match))
      { j2 = 0;
        while ((j2<d) && (xev[*nv*d+j2] == xev[j*d+j2])) j2++;
        match = (j2==d);
        if (!match) j++;
      }
      ii = i+tk;
      clef[i] = cpar[i];
      clef[ii]= crig[i] = j;
      crig[ii]= cpar[ii];
      if (!match) (*nv)++;
    }
  }
  return;
}

extern void hermite2();

double blend(lf,s,x,ll,ur,j,nt,t,what)
lfit *lf;
double s, *x, *ll, *ur;
INT j, nt, *t, what;
{ INT k, k1, m, nc, j0, j1, *ce;
  double v0, v1, xibar, g0[3], g1[3], gg[4], gp[4], phi[4];
  ce = lf->ce;
  for (k=0; k<4; k++)  /* North South East West */
  { k1 = (k>1);
    v0 = ll[k1]; v1 = ur[k1];
    j0 = ce[j+2*(k==0)+(k==2)];
    j1 = ce[j+3-2*(k==1)-(k==3)];
    xibar = (k%2==0) ? ur[k<2] : ll[k<2];
    m = nt;
    while ((m>=0) && ((lf->s[t[m]] != (k<=1)) | (lf->sv[t[m]] != xibar))) m--;
    if (m >= 0)
    { m = (k%2==1) ? lf->lo[t[m]] : lf->hi[t[m]];
      while (lf->s[m] != -1)
        m = (x[lf->s[m]] < lf->sv[m]) ? lf->lo[m] : lf->hi[m];
      if (v0 < evptx(lf,ce[4*m+2*(k==1)+(k==3)],k1))
      { j0 = ce[4*m+2*(k==1)+(k==3)];
        v0 = evptx(lf,j0,k1);
      }
      if (evptx(lf,ce[4*m+3-2*(k==0)-(k==2)],k1) < v1)
      { j1 = ce[4*m+3-2*(k==0)-(k==2)];
        v1 = evptx(lf,j1,k1);
      }
    }
    nc = exvval(lf,g0,j0,2,what,0);
    nc = exvval(lf,g1,j1,2,what,0);
    if (nc==1)
      gg[k] = linear_interp((x[(k>1)]-v0),v1-v0,g0[0],g1[0]);
    else
    { hermite2(x[(k>1)]-v0,v1-v0,phi);
      gg[k] = phi[0]*g0[0]+phi[1]*g1[0]+(phi[2]*g0[1+k1]+phi[3]*g1[1+k1])*(v1-v0);
      gp[k] = phi[0]*g0[2-k1] + phi[1]*g1[2-k1];
    }
  }
  s = -s;
  if (nc==1)
    for (k=0; k<2; k++)
      s += linear_interp(x[k]-ll[k],ur[k]-ll[k],gg[3-2*k],gg[2-2*k]);
    else
    for (k=0; k<2; k++) /* EW NS */
    { hermite2(x[k]-ll[k],ur[k]-ll[k],phi);
      s += phi[0]*gg[3-2*k] + phi[1]*gg[2-2*k]
          +(phi[2]*gp[3-2*k] + phi[3]*gp[2-2*k]) * (ur[k]-ll[k]);
    }
  return(s);
}

double kdtre_int(lf,x,what)
lfit *lf;
double *x;
INT what;
{ INT d, vc, k, t[20], nt, nc, *ce, j;
  double *ll, *ur, ff, vv[64][64];
  d = lf->mi[MDIM];
  vc = lf->vc;
  if (d > 6) ERROR(("d too large in kdint"));

  /* descend the tree to find the terminal cell */
  nt = 0; t[nt] = 0; k = 0;
  while (lf->s[k] != -1)
  { nt++;
    if (nt>=20) { ERROR(("Too many levels in kdint")); return(NOSLN); }
    k = t[nt] = (x[lf->s[k]] < lf->sv[k]) ? lf->lo[k] : lf->hi[k];
  }

  ce = &lf->ce[k*vc];
  ll = evpt(lf,ce[0]);
  ur = evpt(lf,ce[vc-1]);
  nc = 0;
  for (j=0; j<vc; j++) nc = exvval(lf,vv[j],ce[j],d,what,0);
  ff = rectcell_interp(x,NULL,vv,ll,ur,d,nc);

  if (d==2) ff = blend(lf,ff,x,ll,ur,k*vc,nt,t,what);
  return(ff);
}
