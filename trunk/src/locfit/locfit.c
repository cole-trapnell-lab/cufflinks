/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

static double s0, s1, tol;
static lfit   *lf_lf;
static design *lf_des;
int lf_status;
int ident=0;
int (*like)();
extern double robscale;

int likereg(coef, lk0, f1, Z)
double *coef, *lk0, *f1, *Z;
{ INT i, ii, j, p, *mi;
    double lk, ww, link[LLEN], *X;
    lf_status = LF_OK;
    lk = 0.0; p = lf_des->p;
    mi = lf_lf->mi;
    setzero(Z,p*p);
    setzero(f1,p);
    for (i=0; i<lf_des->n; i++)
    { ii = lf_des->ind[i];
        X = d_xi(lf_des,i);
        lf_des->th[i] = base(lf_lf,ii)+innerprod(coef,X,p);
        lf_status = stdlinks(link,lf_lf,ii,lf_des->th[i],robscale);
        if (lf_status == LF_BADP)
        { *lk0 = -1.0e300;
            return(NR_REDUCE);
        }
        if (lf_error) lf_status = LF_ERR;
        if (lf_status != LF_OK) return(NR_BREAK);
        
        ww = lf_des->w[i];
        lk += ww*link[ZLIK];
        for (j=0; j<p; j++)
            f1[j] += X[j]*ww*link[ZDLL];
        addouter(Z, X, X, p, ww*link[ZDDLL]);
    }
    if (mi[MDEB]>2) prresp(coef,Z,p);
    if (mi[MDEB]>1) printf("  likelihood: %8.5f\n",lk);
    *lk0 = lf_des->llk = lk;
    
    switch (lf_lf->mi[MTG]&63) /* parameter checks */
    { case TGAUS: /* prevent iterations! */
            if ((mi[MLINK]==LIDENT)&((mi[MTG]&128)==0)) return(NR_BREAK);
            break;
        case TPOIS:
        case TGEOM:
        case TWEIB:
        case TGAMM:
            if ((mi[MLINK]==LLOG) && (fabs(coef[0])>700))
            { lf_status = LF_OOB;
                return(NR_REDUCE);
            }
            if (lk > -1.0e-5*s0)
            { lf_status = LF_PF;
                return(NR_REDUCE);
            }
            break;
        case TRBIN:
        case TLOGT:
            if (lk > -1.0e-5*s0)
            { lf_status = LF_PF;
                return(NR_REDUCE);
            }
            if (fabs(coef[0])>700)
            { lf_status = LF_OOB;
                return(NR_REDUCE);
            }
            break;
    }
    return(NR_OK);
}

INT robustinit(lf,des)
lfit *lf;
design *des;
{ int i;
    for (i=0; i<des->n; i++)
        des->res[i] = resp(lf,des->ind[i])-base(lf,des->ind[i]);
    des->cf[0] = median(des->res,des->n);
    for (i=1; i<des->p; i++) des->cf[i] = 0.0;
    tol = 1.0e-6;
    return(LF_OK);
}

INT circinit(lf,des)
lfit *lf;
design *des;
{ int i, ii;
    double s0, s1;
    s0 = s1 = 0.0;
    for (i=0; i<des->n; i++)
    { ii = des->ind[i];
        s0 += des->w[i]*prwt(lf,ii)*sin(resp(lf,ii)-base(lf,ii));
        s1 += des->w[i]*prwt(lf,ii)*cos(resp(lf,ii)-base(lf,ii));
    }
    des->cf[0] = atan2(s0,s1);
    for (i=1; i<des->p; i++) des->cf[i] = 0.0;
    tol = 1.0e-6;
    return(LF_OK);
}

INT reginit(lf,des)
lfit *lf;
design *des;
{ int i, ii;
    double sb, link[LLEN];
    s0 = s1 = sb = 0;
    for (i=0; i<des->n; i++)
    { ii = des->ind[i];
        links(base(lf,ii),resp(lf,ii),lf->mi[MTG],LINIT,link,cens(lf,ii),prwt(lf,ii),1.0);
        s1 += des->w[i]*link[ZDLL];
        s0 += des->w[i]*prwt(lf,ii);
        sb += des->w[i]*prwt(lf,ii)*base(lf,ii);
    }
    if (s0==0) return(LF_NOPT); /* no observations with W>0 */
    setzero(des->cf,des->p);
    tol = 1.0e-6*s0;
    switch(lf->mi[MLINK])
    { case LIDENT:
            des->cf[0] = (s1-sb)/s0;
            return(LF_OK);
        case LLOG:
            if (s1<=0.0)
            { des->cf[0] = -1000;
                return(LF_INFA);
            }
            des->cf[0] = log(s1/s0) - sb/s0;
            return(LF_OK);
        case LLOGIT:
            if (s1<=0.0)
            { des->cf[0] = -1000;
                return(LF_INFA);
            }
            if (s1>=s0)
            { des->cf[0] = +1000;
                return(LF_INFA);
            }
            des->cf[0] = logit(s1/s0)-sb/s0;
            return(LF_OK);
        case LINVER:
            if (s1<=0.0)
            { des->cf[0] = 1000;
                return(LF_INFA);
            }
            des->cf[0] = s0/s1-sb/s0;
            return(LF_OK);
        case LSQRT:
            des->cf[0] = sqrt(s1/s0)-sb/s0;
            return(LF_OK);
        case LASIN:
            des->cf[0] = asin(sqrt(s1/s0))-sb/s0;
            return(LF_OK);
        default:
            ERROR(("reginit: invalid link %d",lf->mi[MLINK]));
            return(LF_ERR);
    }
}

int lfinit(lf,des)
lfit *lf;
design *des;
{ 
    //double u[MXDIM];
    INT *mi;
    
    mi = lf->mi;
    des->xtwx.sm = (mi[MDEG0]<mi[MDEG]) ? JAC_CHOL : JAC_EIGD;
    
    designmatrix(lf,des);
    
    like = likereg;
    switch(mi[MTG]&63)
    { case TDEN:
        case TRAT:
        case THAZ:
            like = likeden;
            tol = (mi[MLINK]==LLOG) ? 1.0e-6 : 0.0;
            return(densinit(lf,des,des->h,des->cf,des->n));
        case TCAUC:
        case TROBT:
            return(robustinit(lf,des));
        case TCIRC:
            return(circinit(lf,des));
        default:
            return(reginit(lf,des));
    }
}

void lfiter(lf,des)
lfit *lf;
design *des;
{ int err;
    max_nr(like, des->cf, des->oc, des->res, des->f1,
           &des->xtwx, (int)des->p, (int)lf->mi[MMXIT], tol, &err);
    switch(err)
    { case NR_OK: return;
        case NR_NCON:
            WARN(("max_nr not converged"));
            return;
        case NR_NDIV:
            WARN(("max_nr reduction problem"));
            return;
    }
    WARN(("max_nr return status %d",err));
}

int use_robust_scale(int tg)
{ if ((tg&64)==0) return(0); /* not quasi - no scale */
    if (((tg&128)==0) & (((tg&63)!=TROBT) & ((tg&63)!=TCAUC))) return(0);
    return(1);
}

int locfit(lf,des,h,noit)
lfit *lf;
design *des;
double h;
int noit;
{ int i, p;
    
    if (lf->mi[MDEB]>0)
    { printf("locfit: ");
        for (i=0; i<lf->mi[MDIM]; i++) printf(" %10.6f",des->xev[i]);
        printf("  h = %8.5f\n",h);
    }
    
    lf_lf  = lf;
    lf_des = des;
    des->h = h;
    p = des->p;
    
    lf_status = lfinit(lf,des);
    if (lf_status != LF_OK) return(lf_status);
    
    if (use_robust_scale((int)lf->mi[MTG]))
        lf_robust(lf,des);
    else
    { robscale = 1.0;
        lfiter(lf,des);
    }
    
    if (lf_status == LF_OOB)
        for (i=1; i<p; i++) des->cf[i] = 0.0;
    
    if ((lf->mi[MTG]&63)==TDEN) /* convert from rate to density */
    { switch(lf->mi[MLINK])
        { case LLOG:
                des->cf[0] -= log(lf->dp[DSWT]);
                break;
            case LIDENT:
                multmatscal(des->cf,1.0/lf->dp[DSWT],des->p);
                break;
            default: ERROR(("Density adjustment; invalid link"));
        }
    }
    
    return(lf_status);
}
