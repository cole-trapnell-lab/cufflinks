/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *  Functions for determining bandwidth; smoothing neighborhood
 *  and smoothing weights.
 */

#include "local.h"

extern vari *vb;

double rho(x,sc,d,kt,sty) /* ||x|| for appropriate distance metric */
double *x, *sc;
INT d, kt, *sty;
{ double rhoi[MXDIM], s;
  INT i;
  for (i=0; i<d; i++)
  { if (sty!=NULL)
    { switch(sty[i])
      { case STANGL:  rhoi[i] = 2*sin(x[i]/(2*sc[i])); break;
        case STCPAR: rhoi[i] = 0; break;
        default: rhoi[i] = x[i]/sc[i];
    } }
    else rhoi[i] = x[i]/sc[i];
  }

  if (d==1) return(fabs(rhoi[0]));

  s = 0;
  if (kt==KPROD)
  { for (i=0; i<d; i++)
    { rhoi[i] = fabs(rhoi[i]);
      if (rhoi[i]>s) s = rhoi[i];
    }
    return(s);
  }

  if (kt==KSPH)
  { for (i=0; i<d; i++)
      s += rhoi[i]*rhoi[i];
    return(sqrt(s));
  }

  ERROR(("rho: invalid kt"));
  return(0.0);
}

double kordstat(x,k,n,ind)
double *x;
INT k, n, *ind;
{ INT i, i0, i1, l, r;
  double piv;
  if (k<1) return(0.0);
  i0 = 0; i1 = n-1;
  while (1)
  { piv = x[ind[(i0+i1)/2]];
    l = i0; r = i1;
    while (l<=r)
    { while ((l<=i1) && (x[ind[l]]<=piv)) l++;
      while ((r>=i0) && (x[ind[r]]>piv)) r--;
      if (l<=r) ISWAP(ind[l],ind[r]);
    } /* now, x[ind[i0..r]] <= piv < x[ind[l..i1]] */
    if (r<k-1) i0 = l;  /* go right */
    else /* put pivots in middle */
    { for (i=i0; i<=r; )
        if (x[ind[i]]==piv) { ISWAP(ind[i],ind[r]); r--; }
        else i++;
      if (r<k-1) return(piv);
      i1 = r;
    }
  }
}

/* check if i'th data point is in limits */
INT inlim(lf,xlim,i,d)
lfit *lf;
double *xlim;
INT i, d;
{ INT j, k;
  k = 1;
  for (j=0; j<d; j++)
  { if (xlim[j]<xlim[j+d])
      k &= ((datum(lf,j,i)>=xlim[j]) & (datum(lf,j,i)<=xlim[j+d]));
  }
  return(k);
}

double compbandwid(di,ind,x,n,d,nn,fxh)
double *di, *x, fxh;
INT n, d, *ind, nn;
{ INT i;
  double nnh;

#ifdef CVERSION
  if (nn<0)
    return(dareval(vb,0,x));
#endif

  if (nn<=0) return(fxh);

  if (nn<n)
    nnh = kordstat(di,nn,n,ind);
  else
  { nnh = 0;
    for (i=0; i<n; i++) nnh = MAX(nnh,di[i]);
    nnh = nnh*exp(log(1.0*nn/n)/d);
  }
  return(MAX(fxh,nnh));
}

/*
  fast version of nbhd for ordered 1-d data
*/
double nbhd1(lf,des,k,fxh)
lfit *lf;
design *des;
INT k;
double fxh;
{ double x, h, *xd;
  INT i, l, r, m, n, z;
  n = lf->mi[MN];
  x = des->xev[0];
  xd = dvari(lf,0);

  /* find closest data point to x */
  if (x<=xd[0]) z = 0;
  else
  if (x>=xd[n-1]) z = n-1;
  else
  { l = 0; r = n-1;
    while (r-l>1)
    { z = (r+l)/2;
      if (xd[z]>x) r = z;
              else l = z;
    }
    /* now, xd[0..l] <= x < x[r..n-1] */
    if ((x-xd[l])>(xd[r]-x)) z = r; else z = l;
  }
  /* closest point to x is xd[z] */

  if (k>0) /* set h to nearest neighbor bandwidth */
  { l = r = z;
    if (l==0) r = k-1;
    if (r==n-1) l = n-k;
    while (r-l<k-1)
    { if ((x-xd[l-1])<(xd[r+1]-x)) l--; else r++;
      if (l==0) r = k-1;
      if (r==n-1) l = n-k;
    }
    h = x-xd[l];
    if (h<xd[r]-x) h = xd[r]-x;
  } else h = 0;

  if (h<fxh) h = fxh;

  m = 0;
  if (xd[z]>x) z--; /* so xd[z]<=x */
  /* look left */
  for (i=z; i>=0; i--) if (inlim(lf,lf->xl,i,1))
  { des->di[i] = x-xd[i];
    des->w[m] = weight(lf,&xd[i],&x,h,1,des->di[i]);
    if (des->w[m]>0)
    { des->ind[m] = i;
      m++; 
    } else i = 0;
  }
  /* look right */
  for (i=z+1; i<n; i++) if (inlim(lf,lf->xl,i,1))
  { des->di[i] = xd[i]-x;
    des->w[m] = weight(lf,&xd[i],&x,h,1,des->di[i]);
    if (des->w[m]>0)
    { des->ind[m] = i;
      m++; 
    } else i = n;
  }

  des->n = m;
  return(h);
}

double nbhd(lf,des,nn,fxh,redo)
lfit *lf;
design *des;
double fxh;
INT redo, nn;
{ INT d, i, j, m, n, *mi;
  double h, u[MXDIM];

  mi = lf->mi;

  if (mi[MKT]==KCE) return(0.0);
  d = mi[MDIM]; n = mi[MN];

  /* ordered 1-dim; use fast searches */
  if ((nn<=n) & (lf->ord) & (mi[MKER]!=WMINM) & (lf->sty[0]!=STANGL))
    return(nbhd1(lf,des,nn,fxh));

  if (!redo)
  { for (i=0; i<n; i++)
    { for (j=0; j<d; j++) u[j] = datum(lf,j,i)-des->xev[j];
      des->di[i] = rho(u,lf->sca,d,mi[MKT],lf->sty);
      des->ind[i] = i;
    }
  }
  else
    for (i=0; i<n; i++) des->ind[i] = i;

  if (mi[MKER]==WMINM) return(minmax(lf,des));

  h = compbandwid(des->di,des->ind,des->xev,n,mi[MDIM],nn,fxh);

  m = 0;
  for (i=0; i<n; i++) if (inlim(lf,lf->xl,i,d))
  { for (j=0; j<d; j++) u[j] = datum(lf,j,i);
    des->w[m] = weight(lf,u,des->xev,h,1,des->di[i]);
    if (des->w[m]>0)
    { des->ind[m] = i;
      m++;
    }
  }
  des->n = m;
  return(h);
}
