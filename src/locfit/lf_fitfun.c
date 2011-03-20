/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   Evaluate the locfit fitting functions.
 *     calcp(mi,deg)
 *       calculates the number of fitting functions.
 *     makecfn(des,lf)
 *       makes the coef.number vector.
 *     fitfun(lf,x,t,f,der,nd)
 *       lf is the local fit structure.
 *       x is the data point.
 *       t is the fitting point.
 *       f is a vector to return the results.
 *       deriv is a vector of derivative variables.
 *       nd is the number of derivatives.
 *     designmatrix(lf, des)
 *       is a wrapper for fitfun to build the design matrix.
 *
 */

#include "local.h"

INT calcp(mi,deg)
INT *mi, deg;
{ INT i, k;
    
    if (mi[MUBAS]) return(mi[MP]);
    
    switch (mi[MKT])
    { case KSPH:
        case KCE:
            k = 1;
            for (i=1; i<=deg; i++) k = k*(mi[MDIM]+i)/i;
            return(k);
        case KPROD: return(mi[MDIM]*deg+1);
        case KLM: return(mi[MDIM]);
    }
    ERROR(("calcp: invalid kt %d",mi[MKT]));
    return(0);
}

INT coefnumber(deriv,nd,kt,d,deg)
INT *deriv, nd, kt, d, deg;
{ INT d0, d1, t;
    
    if (d==1)
    { if (nd<=deg) return(nd);
        return(-1);
    }
    
    if (nd==0) return(0);
    if (deg==0) return(-1);
    if (nd==1) return(1+deriv[0]);
    if (deg==1) return(-1);
    if (kt==KPROD) return(-1);
    
    if (nd==2)
    { d0 = deriv[0]; d1 = deriv[1];
        if (d0<d1) { t = d0; d0 = d1; d1 = t; }
        return((d+1)*(d0+1)-d0*(d0+3)/2+d1);
    }
    if (deg==2) return(-1);
    
    ERROR(("coefnumber not programmed for nd>=3"));
    return(-1);
}

void makecfn(des,lf)
design *des;
lfit *lf;
{ int i;
    INT *mi, nd;
    
    nd = lf->nd;
    mi = lf->mi;
    
    des->cfn[0] = coefnumber(lf->deriv,nd,mi[MKT],mi[MDIM],mi[MDEG]);
    des->ncoef = 1;
    if (nd >= mi[MDEG]) return;
    if (mi[MDIM]>1)
    { if (nd>=2) return;
        if ((nd>=1) && (mi[MKT]==KPROD)) return;
    }
    
    for (i=0; i<mi[MDIM]; i++)
    { lf->deriv[nd] = i;
        des->cfn[i+1] = coefnumber(lf->deriv,nd+1,mi[MKT],mi[MDIM],mi[MDEG]);
    }
    des->ncoef = 1+mi[MDIM];
}

void fitfunangl(dx,ff,sca,cd,deg)
double dx, *ff, sca;
INT deg, cd;
{
    if (deg>=3) WARN(("Can't handle angular model with deg>=3"));
    
    switch(cd)
    { case 0:
            ff[0] = 1;
            ff[1] = sin(dx/sca)*sca;
            ff[2] = (1-cos(dx/sca))*sca*sca;
            return;
        case 1:
            ff[0] = 0;
            ff[1] = cos(dx/sca);
            ff[2] = sin(dx/sca)*sca;
            return;
        case 2:
            ff[0] = 0;
            ff[1] = -sin(dx/sca)/sca;
            ff[2] = cos(dx/sca);
            return;
        default: WARN(("Can't handle angular model with >2 derivs"));
    }
}

void fitfun(lf,x,t,f,deriv,nd)
lfit *lf;
double *x, *t, *f;
INT *deriv, nd;
{ 
    INT d, deg, m, i, j, k, ct_deriv[MXDIM];
    double ff[MXDIM][1+MXDEG], dx[MXDIM];
    
#ifdef SVERSION
    if (lf->mi[MUBAS])
    { if (nd>0) WARN(("User basis does not take derivatives"));
        basis(x,t,f,lf->mi[MDIM],lf->mi[MP]);
        return;
    }
#endif
    
    d = lf->mi[MDIM];
    deg = lf->mi[MDEG];
    m = 0;
    
    if (lf->mi[MKT]==KLM)
    { for (i=0; i<d; i++) f[m++] = x[i];
        return;
    }
    memset(ct_deriv, 0, sizeof(ct_deriv));
    for (i = 0; i < MXDIM; ++i)
    {
        for (j = 0; j < 1 + MXDEG; ++j)
        {
            ff[i][j] = 0;
        }
    }
    
    f[m++] = (nd==0);
    if (deg==0) return;
    
    for (i=0; i<d; i++)
    { ct_deriv[i] = 0;
        dx[i] = (t==NULL) ? x[i] : x[i]-t[i];
    }
    for (i=0; i<nd; i++) ct_deriv[deriv[i]]++;
    
    for (i=0; i<d; i++)
    { switch(lf->sty[i])
        {
            case STANGL:
                fitfunangl(dx[i],ff[i],lf->sca[i],ct_deriv[i],lf->mi[MDEG]);
                break;
            default:
                for (j=0; j<ct_deriv[i]; j++) ff[i][j] = 0.0;
                ff[i][ct_deriv[i]] = 1.0;
                for (j=ct_deriv[i]+1; j<=deg; j++)
                    ff[i][j] = ff[i][j-1]*dx[i]/(j-ct_deriv[i]);
        }
    }
    
    /*
     *  Product kernels. Note that if ct_deriv[i] != nd, that implies
     *  there is differentiation wrt another variable, and all components
     *  involving x[i] are 0.
     */
    if ((d==1) || (lf->mi[MKT]==KPROD))
    { for (j=1; j<=deg; j++)
        for (i=0; i<d; i++)
            f[m++] = (ct_deriv[i]==nd) ? ff[i][j] : 0.0;
        return;
    }
    
    /*
     *  Spherical kernels with the full polynomial basis.
     *  Presently implemented up to deg=3.
     */
    for (i=0; i<d; i++)
        f[m++] = (ct_deriv[i]==nd) ? ff[i][1] : 0.0;
    if (deg==1) return;
    
    for (i=0; i<d; i++)
    {
        /* xi^2/2 terms. */
        f[m++] = (ct_deriv[i]==nd) ? ff[i][2] : 0.0;
        
        /* xi xj terms */
        for (j=i+1; j<d; j++)
            f[m++] = (ct_deriv[i]+ct_deriv[j]==nd) ? ff[i][1]*ff[j][1] : 0.0;
    }
    if (deg==2) return;
    
    for (i=0; i<d; i++)
    { 
        /* xi^3/6 terms */
        f[m++] = (ct_deriv[i]==nd) ? ff[i][3] : 0.0;
        
        /* xi^2/2 xk terms */
        for (k=i+1; k<d; k++)
            f[m++] = (ct_deriv[i]+ct_deriv[k]==nd) ? ff[i][2]*ff[k][1] : 0.0;
        
        /* xi xj xk terms */
        for (j=i+1; j<d; j++)
        { f[m++] = (ct_deriv[i]+ct_deriv[j]==nd) ? ff[i][1]*ff[j][2] : 0.0;
            for (k=j+1; k<d; k++)
                f[m++] = (ct_deriv[i]+ct_deriv[j]+ct_deriv[k]==nd) ?
                ff[i][1]*ff[j][1]*ff[k][1] : 0.0;
        }
    }
    if (deg==3) return;
    
    ERROR(("fitfun: can't handle deg=%d for spherical kernels",deg));
}

/*
 *  Build the design matrix. Assumes des->ind contains the indices of
 *  the required data points; des->n the number of points; des->xev
 *  the fitting point.
 */
void designmatrix(lf,des)
lfit *lf;
design *des;
{ int i, ii, j, p;
    double *X, u[MXDIM];
    
    X = d_x(des);
    p = des->p;
    
    if (lf->mi[MUBAS])
    {
#ifdef SVERSION
        vbasis(lf->x,des->xev,lf->mi[MN],lf->mi[MDIM],des->ind,des->n,p,X);
#else
        ERROR(("user basis in S version only\n"));
#endif
        return;
    }
    
    for (i=0; i<des->n; i++)
    { ii = des->ind[i];
        for (j=0; j<lf->mi[MDIM]; j++) u[j] = datum(lf,j,ii);
        fitfun(lf,u,des->xev,&X[i*p],NULL,(INT)0);
    }
}
