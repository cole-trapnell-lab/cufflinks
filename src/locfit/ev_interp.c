/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

double linear_interp(h,d,f0,f1)
double h, d, f0, f1;
{ if (d==0) return(f0);
  return( ( (d-h)*f0 + h*f1 ) / d );
}

void hermite2(x,z,phi)
double x, z, *phi;
{ double h;
  if (z==0)
  { phi[0] = 1.0; phi[1] = phi[2] = phi[3] = 0.0;
    return;
  }
  h = x/z;
  if (h<0)
  { phi[0] = 1; phi[1] = 0;
    phi[2] = h; phi[3] = 0;
    return;
  }
  if (h>1)
  { phi[0] = 0; phi[1] = 1;
    phi[2] = 0; phi[3] = h-1;
    return;
  }
  phi[1] = h*h*(3-2*h);
  phi[0] = 1-phi[1];
  phi[2] = h*(1-h)*(1-h);
  phi[3] = h*h*(h - 1);
}

double cubic_interp(h,f0,f1,d0,d1)
double h, f0, f1, d0, d1;
{ double phi[4];
  hermite2(h,1.0,phi);
  return(phi[0]*f0+phi[1]*f1+phi[2]*d0+phi[3]*d1);
}

double cubintd(h,f0,f1,d0,d1)
double h, f0, f1, d0, d1;
{ double phi[4];
  phi[1] = 6*h*(1-h);
  phi[0] = -phi[1];
  phi[2] = (1-h)*(1-3*h);
  phi[3] = h*(3*h-2);
  return(phi[0]*f0+phi[1]*f1+phi[2]*d0+phi[3]*d1);
}

/*
  interpolate over a rectangular cell.
    x = interpolation point. 
    xev = evaluation points (do I need this?)
    vv = array of vertex values.
    ll = lower left corner.
    ur = upper right corner.
    d = dimension.
    nc = no of coefficients.
*/
double rectcell_interp(x,xev,vv,ll,ur,d,nc)
double *x, *xev, vv[64][64], *ll, *ur;
INT d, nc;
{ double phi[4];
  INT i, j, k, tk;

  tk = 1<<d;
  for (i=0; i<tk; i++) if (vv[i][0]==NOSLN) return(NOSLN);

  /* no derivatives - use multilinear interpolation */
  if (nc==1)
  { for (i=d-1; i>=0; i--)
    { tk = 1<<i;
      for (j=0; j<tk; j++)
        vv[j][0] = linear_interp(x[i]-ll[i],ur[i]-ll[i],vv[j][0],vv[j+tk][0]);
    }
    return(vv[0][0]);
  }

  /* with derivatives -- use cubic */
  if (nc==d+1)
  { for (i=d-1; i>=0; i--)
    { hermite2(x[i]-ll[i],ur[i]-ll[i],phi);
      tk = 1<<i;
      phi[2] *= ur[i]-ll[i];
      phi[3] *= ur[i]-ll[i];
      for (j=0; j<tk; j++)
      { vv[j][0] = phi[0]*vv[j][0] + phi[1]*vv[j+tk][0]
                 + phi[2]*vv[j][i+1] + phi[3]*vv[j+tk][i+1];
        for (k=1; k<=i; k++)
          vv[j][k] = phi[0]*vv[j][k] + phi[1]*vv[j+tk][k];
      }
    }
    return(vv[0][0]); 
  }

  /* with all coefs -- use multicubic */
  for (i=d-1; i>=0; i--)
  { hermite2(x[i]-ll[i],ur[i]-ll[i],phi);
    tk = 1<<i;
    phi[2] *= ur[i]-ll[i];
    phi[3] *= ur[i]-ll[i];
    for (j=0; j<tk; j++)
      for (k=0; k<tk; k++)
        vv[j][k] = phi[0]*vv[j][k] + phi[1]*vv[j+tk][k]
                 + phi[2]*vv[j][k+tk] + phi[3]*vv[j+tk][k+tk];
  }
  return(vv[0][0]);
}

INT exvval(lf,vv,nv,d,what,z)
lfit *lf;
double *vv;
INT nv, d, what, z;
{ INT i, k;
  double *values;
  k = (z) ? 1<<d : d+1;
  for (i=1; i<k; i++) vv[i] = 0.0;
  switch(what)
  { case PCOEF:
      values = lf->coef;
      break;
    case PVARI:
    case PNLX:
      values = lf->nlx;
      break;
    case PT0:
      values = lf->t0;
      break;
    case PBAND:
      vv[0] = lf->h[nv];
      return(1);
    case PDEGR:
      vv[0] = lf->deg[nv];
      return(1);
    case PLIK:
      vv[0] = lf->lik[nv];
      return(1);
    case PRDF:
      vv[0] = lf->lik[2*lf->nvm+nv];
      return(1);
    default:
      ERROR(("Invalid what in exvval"));
      return(0);
  }
  vv[0] = values[nv];
  if ((lf->mi[MDEG]==0) && (lf->mi[MDC]==0))return(1);
  if (z)
  { for (i=0; i<d; i++) vv[1<<i] = values[(i+1)*lf->nvm+nv];
    return(1<<d);
  }
  else
  { for (i=1;i<=d; i++) vv[i] = values[i*lf->nvm+nv];
    return(d+1);
  }
}

void exvvalpv(vv,vl,vr,d,k,dl,nc)
double *vv, *vl, *vr, dl;
INT d, k, nc;
{ INT i, tk, td;
  double f0, f1;
  if (nc==1)
  { vv[0] = (vl[0]+vr[0])/2;
    return;
  }
  tk = 1<<k;
  td = 1<<d;
  for (i=0; i<td; i++) if ((i&tk)==0)
  { f0 = (vl[i]+vr[i])/2 + dl*(vl[i+tk]-vr[i+tk])/8;
    f1 = 1.5*(vr[i]-vl[i])/dl - (vl[i+tk]+vr[i+tk])/4;
    vv[i] = f0;
    vv[i+tk] = f1;
  }
} 

double gridint(tr,x,what)
lfit *tr;
double *x;
INT what;
{ INT d, i, j, jj, v[MXDIM], vc, z0, sk, nce[1024], nc;
  double *ll, *ur, vv[64][64];
  d = tr->mi[MDIM];
  ll = evpt(tr,0); ur = evpt(tr,tr->nv-1);
  z0 = 0; vc = 1<<d;
  for (j=d-1; j>=0; j--)
  { v[j] = (INT)((tr->mg[j]-1)*(x[j]-ll[j])/(ur[j]-ll[j]));
    if (v[j]<0) v[j]=0;
    if (v[j]>=tr->mg[j]-1) v[j] = tr->mg[j]-2;
    z0 = z0*tr->mg[j]+v[j];
  }
  nce[0] = z0; nce[1] = z0+1; sk = jj = 1; 
    memset(nce, 0, sizeof(nce));
  for (i=1; i<d; i++)
  { sk *= tr->mg[i-1];
    jj<<=1;
    for (j=0; j<jj; j++)
      nce[j+jj] = nce[j]+sk;
  }
    nc = 0;
  for (i=0; i<vc; i++)
    nc = exvval(tr,vv[i],nce[i],d,what,1);
  ll = evpt(tr,nce[0]);
  ur = evpt(tr,nce[vc-1]);
  return(rectcell_interp(x,vdptr(tr->xxev),vv,ll,ur,d,nc));
}

double fitpint(lf,x,what,i)
lfit *lf;
double *x;
INT what, i;
{ double vv[1+MXDIM];
  exvval(lf,vv,i,lf->mi[MDIM],what,0);
  return(vv[0]);
}

double dointpointpf(lf,des,x,what)
lfit *lf;
design *des;
double *x;
INT what;
{ locfit(lf,des,0.0,0);
  if (what==PCOEF) return(des->cf[0]);
  if ((what==PNLX)|(what==PT0)) return(sqrt(comp_infl(lf,des)));
  ERROR(("dointpointpf: invalid what"));
  return(0.0);
}

double xbarint(lf,x,what)
lfit *lf;
double *x;
INT what;
{ INT i, nc;
  double vv[1+MXDIM], f;
  nc = exvval(lf,vv,0,lf->mi[MDIM],what,0);
  f = vv[0];
  if (nc>1)
    for (i=0; i<lf->mi[MDIM]; i++)
      f += vv[i+1]*(x[i]-evptx(lf,0,i));
  return(f);
}

double dointpoint(lf,des,x,what,ev,j)
lfit *lf;
design *des;
double *x;
INT what, ev, j;
{ double xf, f;
  INT i;
  for (i=0; i<lf->mi[MDIM]; i++) if (lf->sty[i]==STANGL)
  { xf = floor(x[i]/(2*PI*lf->sca[i]));
    x[i] -= xf*2*PI*lf->sca[i];
  }
    f = 0;
  if (ident==1) return(dointpointpf(lf,des,x,what));
  switch(ev)
  { case EGRID: f = gridint(lf,x,what); break;
    case EKDTR: f = kdtre_int(lf,x,what); break;
    case ETREE: f = atree_int(lf,x,what); break;
    case EPHULL: f = triang_int(lf,x,what); break;
    case EFITP: f = fitpint(lf,x,what,j); break;
    case EXBAR: f = xbarint(lf,x,what); break;
    case ENONE: f = 0; break;
    default: ERROR(("dointpoint: cannot interpolate this structure"));
  }
  if (((what==PT0)|(what==PNLX)) && (f<0)) f = 0.0;
  f += addparcomp(lf,x,what);
  return(f);
}
