/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

extern int lf_status;
static double u[MXDIM], ilim[2*MXDIM], *ff;
static lfit   *den_lf;
static design *den_des;
INT fact[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800};

INT multint(), prodint(), gausint(), mlinint();

void prresp(coef,resp,p)
double *coef, *resp;
INT p;
{ INT i, j;
  printf("Coefficients:\n");
  for (i=0; i<p; i++) printf("%8.5f ",coef[i]);
  printf("\n");
  printf("Response matrix:\n");
  for (i=0; i<p; i++)
  { for (j=0; j<p; j++) printf("%9.6f, ",resp[i+j*p]);
    printf("\n");
  }
}

INT multint(t,resp1,resp2,lf,cf,h)
lfit *lf;
double *t, *resp1, *resp2, *cf, h;
{ INT d, p, i, j, k, m, m1, w, z, z1;
  double th, wt, dj[MXDIM];
  d = lf->mi[MDIM]; p = lf->mi[MP];
  setzero(resp1,p*p);
  m = 1; m1 = lf->mi[MMINT]+1;
  for (i=0; i<d; i++)
  { m *= m1;
    dj[i] = (ilim[i+d]-ilim[i])/lf->mi[MMINT];
  }
  for (i=0; i<m; i++)
  { z = i; w = 1;
    for (j=d-1; j>=0; j--)
    { z1 = z%m1;
      u[j] = t[j]+ilim[j]+dj[j]*z1;
      w *= (4-2*(z1%2==0)-(z1==0)-(z1==lf->mi[MMINT]));
      z /= m1;
    }
    wt = w*weight(lf,u,t,h,0,0.0);
    if (wt>0)
    { fitfun(lf,u,t,ff,NULL,0);
      th = innerprod(ff,cf,p);
      switch(lf->mi[MLINK])
      { case LLOG:
          addouter(resp1,ff,ff,p,wt*lf_exp(th));
          break;
        case LIDENT:
          addouter(resp1,ff,ff,p,wt);
          break;
        default:
          ERROR(("multint: Invalid link"));
          return(LF_LNK);
      }
    }
  }
  wt = 1;
  for (j=0; j<d; j++) wt *= dj[j]/3;
  for (j=0; j<p; j++)
    for (k=j; k<p; k++)
      resp1[p*k+j] = resp1[p*j+k] = resp1[p*j+k]*wt;
  return(LF_OK);
}

INT mlinint(t,resp1,resp2,lf,cf,h)
lfit *lf;
double *t, *resp1, *resp2, *cf, h;
{
  double hd, nb, wt, wu, g[4], w0, w1, v, *sca;
  INT d, p, i, j, jmax, k, l, z, jj[2];
  d = lf->mi[MDIM]; p = lf->mi[MP]; sca = lf->sca;
  hd = 1;
  for (i=0; i<d; i++) hd *= h*sca[i];

  if (lf->mi[MLINK]==LIDENT)
  { setzero(resp1,p*p);
    resp1[0] = wint(d,NULL,0,lf->mi[MKER])*hd;
    if (lf->mi[MDEG]==0) return(LF_OK);
    jj[0] = 2; w0 = wint(d,jj,1,lf->mi[MKER])*hd*h*h;
    for (i=0; i<d; i++) resp1[(i+1)*p+i+1] = w0*sca[i]*sca[i];
    if (lf->mi[MDEG]==1) return(LF_OK);
    for (i=0; i<d; i++)
    { j = p-(d-i)*(d-i+1)/2;
      resp1[j] = resp1[p*j] = w0*sca[i]*sca[i]/2;
    }
    if (d>1) { jj[1] = 2; w0 = wint(d,jj,2,lf->mi[MKER])*hd*h*h*h*h; }
    jj[0] = 4; w1 = wint(d,jj,1,lf->mi[MKER])*hd*h*h*h*h/4;
    z = d+1;
    for (i=0; i<d; i++)
    { k = p-(d-i)*(d-i+1)/2;
      for (j=i; j<d; j++)
      { l = p-(d-j)*(d-j+1)/2;
        if (i==j) resp1[z*p+z] = w1*SQR(sca[i])*SQR(sca[i]);
        else
        { resp1[z*p+z] = w0*SQR(sca[i])*SQR(sca[j]);
          resp1[k*p+l] = resp1[k+p*l] = w0/4*SQR(sca[i])*SQR(sca[j]);
        }
        z++;
    } }
    return(LF_OK);
  }
  switch(lf->mi[MDEG])
  { case 0:
      resp1[0] = lf_exp(cf[0])*wint(d,NULL,0,lf->mi[MKER])*hd;
      return(LF_OK);
    case 1:
      nb = 0.0;
      for (i=1; i<=d; i++)
      { v = h*cf[i]*sca[i-1];
        nb += v*v;
      }
      if (lf->mi[MKER]==WGAUS)
      { w0 = 1/(GFACT*GFACT);
        g[0] = lf_exp(cf[0]+w0*nb/2+d*log(S2PI/2.5));
        g[1] = g[3] = g[0]*w0;
        g[2] = g[0]*w0*w0;
      }
      else
      { wt = wu = lf_exp(cf[0]);
        w0 = wint(d,NULL,0,lf->mi[MKER]); g[0] = wt*w0;
        g[1] = g[2] = g[3] = 0.0;
        j = 0; jmax = (d+2)*lf->mi[MMINT];
        while ((j<jmax) && (wt*w0/g[0]>1.0e-8))
        { j++;
          jj[0] = 2*j; w0 = wint(d,jj,1,lf->mi[MKER]);
          if (d==1) g[3] += wt * w0;
          else
          { jj[0] = 2; jj[1] = 2*j-2; w1 = wint(d,jj,2,lf->mi[MKER]);
            g[3] += wt*w1;
            g[2] += wu*(w0-w1);
          }
          wt /= (2*j-1.0); g[1] += wt*w0;
          wt *= nb/(2*j); g[0] += wt*w0;
          wu /= (2*j-1.0)*(2*j);
          if (j>1) wu *= nb;
        }
        if (j==jmax) WARN(("mlinint: series not converged"));
      }
      g[0] *= hd; g[1] *= hd;
      g[2] *= hd; g[3] *= hd;
      resp1[0] = g[0];
      for (i=1; i<=d; i++)
      { resp1[i] = resp1[(d+1)*i] = cf[i]*SQR(h*sca[i-1])*g[1];
        for (j=1; j<=d; j++)
        { resp1[(d+1)*i+j] = (i==j) ? g[3]*SQR(h*sca[i-1]) : 0;
          resp1[(d+1)*i+j] += g[2]*SQR(h*h*sca[i-1]*sca[j-1])*cf[i]*cf[j];
        }
      }
      return(LF_OK);
  }
  ERROR(("mlinint: deg=0,1 only"));
  return(LF_ERR);
}

void prodint_resp(resp,prod_wk,dim,deg,p)
double *resp, prod_wk[MXDIM][2*MXDEG+1];
int dim, deg, p;
{ double prod;
  int i, j, k, j1, k1;

  prod = 1.0;
  for (i=0; i<dim; i++) prod *= prod_wk[i][0];
  resp[0] += prod;
  if (deg==0) return;

  for (j1=1; j1<=deg; j1++)
  { for (j=0; j<dim; j++)
    { prod = 1.0;
      for (i=0; i<dim; i++) prod *= prod_wk[i][j1*(j==i)];
      prod /= fact[j1];
      resp[1 + (j1-1)*dim +j] += prod;
    }
  }

  for (k1=1; k1<=deg; k1++)
    for (j1=k1; j1<=deg; j1++)
    { for (k=0; k<dim; k++)
        for (j=0; j<dim; j++)
        { prod = 1.0;
          for (i=0; i<dim; i++) prod *= prod_wk[i][k1*(k==i) + j1*(j==i)];
          prod /= fact[k1]*fact[j1];
          resp[ (1+(k1-1)*dim+k)*p + 1+(j1-1)*dim+j] += prod;
        }
    }
}

INT prodint(t,resp,resp2,lf,coef,h)
lfit *lf;
double *t, *resp, *resp2, *coef, h;
{ INT dim, deg, p, i, j, k, st;
  double cf[MXDEG+1], hj, hs, prod_wk[MXDIM][2*MXDEG+1];

    st = 0;
  dim = lf->mi[MDIM];
  deg = lf->mi[MDEG];
  p = lf->mi[MP];
  for (i=0; i<p*p; i++) resp[i] = 0.0;
  cf[0] = coef[0];

/*  compute the one dimensional terms
 */
  for (i=0; i<dim; i++)
  { hj = 1; hs = h*lf->sca[i];
    for (j=0; j<deg; j++)
    { hj *= hs;
      cf[j+1] = hj*coef[ j*dim+i+1 ];
    }
    st = onedint(cf,lf->mi,ilim[i]/hs,ilim[i+dim]/hs,prod_wk[i]);
    if (st==LF_BADP) return(st);
    hj = 1;
    for (j=0; j<=2*deg; j++)
    { hj *= hs;
      prod_wk[i][j] *= hj;
    }
    cf[0] = 0.0; /* so we only include it once, when d>=2 */
  }

/*  transfer to the resp array
 */
  prodint_resp(resp,prod_wk,dim,deg,p);

/* Symmetrize.
*/
  for (k=0; k<p; k++)
    for (j=k; j<p; j++)
      resp[j*p+k] = resp[k*p+j];

  return(st);
}

INT gausint(t,resp,C,cf,h,mi,sca)
double *t, *resp, *C, *cf, h, *sca;
INT *mi;
{ double nb, det, z, *P;
  INT d, p, i, j, k, l, m1, m2, f;
  d = mi[MDIM]; p = mi[MP];
  m1 = d+1; nb = 0;
  P = &C[d*d];
  resp[0] = 1;
  for (i=0; i<d; i++)
  { C[i*d+i] = SQR(GFACT/(h*sca[i]))-cf[m1++];
    for (j=i+1; j<d; j++) C[i*d+j] = C[j*d+i] = -cf[m1++];
  }
  eig_dec(C,P,d);
  det = 1;
  for (i=1; i<=d; i++)
  { det *= C[(i-1)*(d+1)];
    if (det <= 0) return(LF_BADP);
    resp[i] = cf[i];
    for (j=1; j<=d; j++) resp[j+i*p] = 0;
    resp[i+i*p] = 1;
    svdsolve(&resp[i*p+1],u,P,C,P,d,0.0);
  }
  svdsolve(&resp[1],u,P,C,P,d,0.0);
  det = sqrt(det);
  for (i=1; i<=d; i++)
  { nb += cf[i]*resp[i];
    resp[i*p] = resp[i];
    for (j=1; j<=d; j++)
      resp[i+p*j] += resp[i]*resp[j];
  }
  m1 = d;
  for (i=1; i<=d; i++)
    for (j=i; j<=d; j++)
    { m1++; f = 1+(i==j);
      resp[m1] = resp[m1*p] = resp[i*p+j]/f;
      m2 = d;
      for (k=1; k<=d; k++)
      { resp[m1+k*p] = resp[k+m1*p] =
        ( resp[i]*resp[j*p+k] + resp[j]*resp[i*p+k]
        + resp[k]*resp[i*p+j] - 2*resp[i]*resp[j]*resp[k] )/f;
        for (l=k; l<=d; l++)
        { m2++; f = (1+(i==j))*(1+(k==l));
          resp[m1+m2*p] = resp[m2+m1*p] = ( resp[i+j*p]*resp[k+l*p]
            + resp[i+k*p]*resp[j+l*p] + resp[i+l*p]*resp[j+k*p]
            - 2*resp[i]*resp[j]*resp[k]*resp[l] )/f;
    } } }
  z = lf_exp(d*0.918938533+cf[0]+nb/2)/det;
  multmatscal(resp,z,p*p);
  return(LF_OK);
}

int likeden(coef, lk0, f1, A)
double *coef, *lk0, *f1, *A;
{ double lk, r;
  int i, j, p, rstat;

    lk = 0;
  lf_status = LF_OK;
  p = den_des->p;
  if ((den_lf->mi[MLINK]==LIDENT) && (coef[0] != 0.0)) return(NR_BREAK);
  lf_status = (den_des->itype)(den_des->xev,A,den_des->xtwx.Q,den_lf,coef,den_des->h);
  if (lf_error) lf_status = LF_ERR;
  if (lf_status==LF_BADP)
  { *lk0 = -1.0e300;
    return(NR_REDUCE);
  }
  if (lf_status!=LF_OK) return(NR_BREAK);
  if (den_lf->mi[MDEB]>2) prresp(coef,A,(INT)p);

  den_des->xtwx.p = p;
  rstat = NR_OK;
  switch(den_lf->mi[MLINK])
  { case LLOG:
      r = den_des->ss[0]/A[0];
      coef[0] += log(r);
      multmatscal(A,r,(INT)p*p);
      A[0] = den_des->ss[0];
      lk = -A[0];
      if (fabs(coef[0]) > 700)
      { lf_status = LF_OOB;
        rstat = NR_REDUCE;
      }
      for (i=0; i<p; i++)
      { lk += coef[i]*den_des->ss[i];
        f1[i] = den_des->ss[i]-A[i];
      }
      break;
    case LIDENT:
      lk = 0.0;
      for (i=0; i<p; i++)
      { f1[i] = den_des->ss[i];
        for (j=0; j<p; j++)
          den_des->res[i] -= A[i*p+j]*coef[j];
      }
      break;
  }
  *lk0 = den_des->llk = lk;

  return(rstat);
}

INT inre(x,bound,d)
double *x, *bound;
INT d;
{ INT i, z;
  z = 1;
  for (i=0; i<d; i++)
    if (bound[i]<bound[i+d])
      z &= (x[i]>=bound[i]) & (x[i]<=bound[i+d]);
  return(z);
}

INT setintlimits(lf, x, h, ang, lset)
lfit *lf;
INT *ang, *lset;
double *x, h;
{ INT d, i;
  d = lf->mi[MDIM];
  *ang = *lset = 0;
  for (i=0; i<d; i++)
  { if (lf->sty[i]==STANGL)
    { ilim[i+d] = ((h<2) ? 2*asin(h/2) : PI)*lf->sca[i];
      ilim[i] = -ilim[i+d];
      *ang = 1;
    }
    else
    { ilim[i+d] = h*lf->sca[i];
      ilim[i] = -ilim[i+d];

      if (lf->sty[i]==STLEFT) { ilim[i+d] = 0; *lset = 1; }
      if (lf->sty[i]==STRIGH) { ilim[i] = 0;   *lset = 1; }

      if (lf->xl[i]<lf->xl[i+d]) /* user limits for this variable */
      { if (lf->xl[i]-x[i]> ilim[i])
        { ilim[i] = lf->xl[i]-x[i]; *lset=1; }
        if (lf->xl[i+d]-x[i]< ilim[i+d])
        { ilim[i+d] = lf->xl[i+d]-x[i]; *lset=1; }
      }
    }
    if (ilim[i]==ilim[i+d]) return(LF_DEMP); /* empty integration */
  }
  return(LF_OK);
}

INT selectintmeth(mi,lset,ang)
INT *mi, lset, ang;
{
  if (mi[MIT]==IDEFA) /* select the default method */
  { if (mi[MTG]==THAZ)
    { if (ang) return(IDEFA);
      return( IHAZD );
    }

    if (mi[MUBAS]) return(IMULT);

    if (ang) return(IMULT);

    if (iscompact(mi[MKER]))
    { if (mi[MKT]==KPROD) return(IPROD);
      if (lset)
        return( (mi[MDIM]==1) ? IPROD : IMULT );
      if (mi[MDEG]<=1) return(IMLIN);
      if (mi[MDIM]==1) return(IPROD);
      return(IMULT);
    }

    if (mi[MKER]==WGAUS)
    { if (lset) WARN(("Integration for Gaussian weights ignores limits"));
      if ((mi[MDIM]==1)|(mi[MKT]==KPROD)) return(IPROD);
      return(IMLIN);
    }

    return(IDEFA);
  }

  /* user provided an integration method, check it is valid */

  if (mi[MTG]==THAZ)
  { if (ang) return(INVLD);
    if (!iscompact(mi[MKER])) return(INVLD);
    return( ((mi[MKT]==KPROD) | (mi[MKT]==KSPH)) ? IHAZD : INVLD );
  }

  if ((ang) && (mi[MIT] != IMULT)) return(INVLD);

  switch(mi[MIT])
  { case IMULT: return( iscompact(mi[MKER]) ? IMULT : INVLD );
    case IPROD: return( ((mi[MDIM]==1) | (mi[MKT]==KPROD)) ? IPROD : INVLD );
    case IMLIN: return( ((mi[MKT]==KSPH) && (!lset) &&
      (mi[MDEG]<=1)) ? IMLIN : INVLD );
  }

  return(INVLD);
}

INT densinit(lf,des,h,cf,m)
lfit *lf;
design *des;
double h, *cf;
INT m;
{ INT deg, p, i, ii, j, nnz, rnz, lset, ang, status;
  double w;

  den_lf  = lf;
  den_des = des;

  p = des->p; deg = lf->mi[MDEG];
  ff = des->xtwx.wk;
  cf[0] = NOSLN;
  for (i=1; i<p; i++) cf[i] = 0.0;

  if (!inre(des->xev,lf->xl,lf->mi[MDIM])) return(LF_XOOR);

  status = setintlimits(lf,des->xev,h,&ang,&lset);
  if (status != LF_OK) return(status);

  switch(selectintmeth(lf->mi,lset,ang))
  { case IMULT: des->itype = multint; break;
    case IPROD: des->itype = prodint; break;
    case IMLIN: des->itype = mlinint; break;
    case IHAZD: des->itype = hazint; break;
    case INVLD: ERROR(("Invalid integration method %d",lf->mi[MIT]));
                break;
    case IDEFA: ERROR(("No integration type available for this model"));
                break;
    default: ERROR(("densinit: unknown integral type"));
  }

  switch(deg)
  { case 0: rnz = 1; break;
    case 1: rnz = 1; break;
    case 2: rnz = lf->mi[MDIM]+1; break;
    case 3: rnz = lf->mi[MDIM]+2; break;
    default: ERROR(("densinit: invalid degree %d",deg));
  }
  if (lf_error) return(LF_ERR);

  setzero(des->ss,p);
  nnz = 0;
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    if (!cens(lf,ii))
    { w = des->w[i]*prwt(lf,ii);
      for (j=0; j<p; j++) des->ss[j] += d_xij(des,i,j)*w;
      if (des->w[i]>0.00001) nnz++;
  } }

  if (lf->mi[MTG]==THAZ) haz_init(lf,des,ilim);

  if (lf->mi[MDEB]>2)
  { printf("    LHS: ");
    for (i=0; i<p; i++) printf(" %8.5f",des->ss[i]);
    printf("\n");
  }

  switch(lf->mi[MLINK])
  { case LIDENT:
      cf[0] = 0.0;
      return(LF_OK);
    case LLOG:
      if (nnz<rnz) { cf[0] = -1000; return(LF_DNOP); }
      cf[0] = 0.0;
      return(LF_OK);
    default:
      ERROR(("unknown link in densinit"));
      return(LF_ERR);
  }
}
