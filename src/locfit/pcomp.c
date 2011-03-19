/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 * 
 *
 * functions for computing and subtracting, adding the
 * parametric component
 */

#include "local.h"

INT noparcomp(lf)
lfit *lf;
{ 
    INT tg;
    if (lf->mi[MDEG0]<lf->mi[MDEG]) return(1);
    if (lf->mi[MUBAS]) return(1);
    tg = lf->mi[MTG] & 63;
    if (tg<=THAZ) return(1);
    if (tg==TROBT) return(1);
    if (tg==TCAUC) return(1);
    return(0);
}

INT hasparcomp(lf)
lfit *lf;
{ 
    return(lf->mi[MPC]);
}

int pc_reqd(d,p)
INT d, p;
{ 
    return(d + 2*p + jac_reqd(p));
}

void pcchk(pc,d,p,lc)
paramcomp *pc;
INT d, p, lc;
{ 
    //INT k;
    double *z;
    pc->wk = checkvarlen(pc->wk,pc_reqd(d,p),"_pcwork",VDOUBLE);
    z = vdptr(pc->wk);
    //k = 0;
    
    pc->xbar = z; z += d;
    pc->coef = z; z += p;
    pc->f    = z; z += p;
    
    z = jac_alloc(&pc->xtwx,p,z);
    pc->xtwx.p = p;
}

void compparcomp(des,lf,nopc)
design *des;
lfit *lf;
INT nopc;
{ 
    INT i, j, k;
    double wt, sw;
    paramcomp *pc;
    pc = &lf->pc;
    pcchk(pc,lf->mi[MDIM],lf->mi[MP],1);
    
    if (pc == NULL)
    {
        fprintf(stderr, "Error: locfit cannot allocate pc working memory\n");
        pcchk(pc,lf->mi[MDIM],lf->mi[MP],1);
        return;
    }
    
    for (i=0; i<lf->mi[MDIM]; i++) pc->xbar[i] = 0.0;
    sw = 0.0;
    for (i=0; i<lf->mi[MN]; i++)
    { 
        wt = prwt(lf,i);
        sw += wt;
        for (j=0; j<lf->mi[MDIM]; j++)
            pc->xbar[j] += datum(lf,j,i)*wt;
        des->ind[i] = i;
        des->w[i] = 1.0;
    }
    for (i=0; i<lf->mi[MDIM]; i++) pc->xbar[i] /= sw;
    if ((nopc) || noparcomp(lf))
    { lf->mi[MPC] = 0;
        return;
    }
    lf->mi[MPC] = 1;
    des->xev = pc->xbar;
    k = locfit(lf,des,0.0,0);

    if (lf_error) return;
    switch(k)
    { case LF_NOPT:
            ERROR(("compparcomp: no points in dataset?"));
            return;
        case LF_INFA:
            ERROR(("compparcomp: infinite parameters in param. component"));
            return;
        case LF_NCON:
            ERROR(("compparcom: not converged"));
            return;
        case LF_OOB:
            ERROR(("compparcomp: parameters out of bounds"));
            return;
        case LF_PF:
            WARN(("compparcomp: perfect fit"));
        case LF_OK:
            for (i=0; i<lf->mi[MP]; i++)
            { pc->coef[i] = des->cf[i];
                pc->xtwx.dg[i] = des->xtwx.dg[i];
                pc->xtwx.wk[i] = des->xtwx.wk[i];
            }
            for (i=0; i<lf->mi[MP]*lf->mi[MP]; i++)
            { pc->xtwx.Z[i] = des->xtwx.Z[i];
                pc->xtwx.Q[i] = des->xtwx.Q[i];
            }
            pc->xtwx.sm = des->xtwx.sm;
            pc->xtwx.st = des->xtwx.st;
            return;
        default:
            ERROR(("compparcomp: locfit unknown return status %d",k));
            return;
    }
}

void subparcomp(des,lf,coef)
design *des;
lfit *lf;
double *coef;
{ INT i, *deriv, nd;
    
    if (!hasparcomp(lf)) return;
    
    deriv = lf->deriv;
    nd = lf->nd;
    fitfun(lf,des->xev,lf->pc.xbar,des->f1,deriv,nd);
    coef[0] -= innerprod(lf->pc.coef,des->f1,lf->mi[MP]);
    if (des->ncoef == 1) return;
    
    for (i=0; i<lf->mi[MDIM]; i++)
    { deriv[nd] = i;
        fitfun(lf,des->xev,lf->pc.xbar,des->f1,deriv,nd+1);
        coef[i+1] -= innerprod(lf->pc.coef,des->f1,lf->mi[MP]);
    }
}

void subparcomp2(des,lf,vr,il)
design *des;
lfit *lf;
double *vr, *il;
{ double t0, t1;
    INT i, *deriv, nd, *mi;
    
    if (!hasparcomp(lf)) return;
    
    mi = lf->mi;
    deriv = lf->deriv;
    nd = lf->nd;
    fitfun(lf,des->xev,lf->pc.xbar,des->f1,deriv,nd);
    for (i=0; i<mi[MP]; i++) lf->pc.f[i] = des->f1[i];
    jacob_solve(&lf->pc.xtwx,des->f1);
    t0 = sqrt(innerprod(lf->pc.f,des->f1,mi[MP]));
    vr[0] -= t0;
    il[0] -= t0;
    if ((t0==0) | (des->ncoef==1)) return;
    
    for (i=0; i<mi[MDIM]; i++)
    { deriv[nd] = i;
        fitfun(lf,des->xev,lf->pc.xbar,lf->pc.f,deriv,nd+1);
        t1 = innerprod(lf->pc.f,des->f1,mi[MP])/t0;
        vr[i+1] -= t1;
        il[i+1] -= t1;
    }
}

double addparcomp(lf,x,c)
lfit *lf;
double *x;
int c;
{ double y;
    if (!hasparcomp(lf)) return(0.0);
    fitfun(lf,x,lf->pc.xbar,lf->pc.f,lf->deriv,lf->nd);
    if (c==PCOEF) return(innerprod(lf->pc.coef,lf->pc.f,lf->mi[MP]));
    if ((c==PNLX)|(c==PT0)|(c==PVARI))
    { y = sqrt(jacob_qf(&lf->pc.xtwx,lf->pc.f));
        return(y);
    }
    return(0.0);
}
