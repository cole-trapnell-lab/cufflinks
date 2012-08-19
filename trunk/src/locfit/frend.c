/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

extern INT cvi;
extern double robscale;

double resp(lf,i)
lfit *lf;
INT i;
{ if (lf->y==NULL) return(0.0);
    return(lf->y[i]);
}

double prwt(lf,i)
lfit *lf;
INT i;
{ if (i==cvi) return(0.0);
    if (lf->w==NULL) return(1.0);
    return(lf->w[i]);
}

double base(lf,i)
lfit *lf;
INT i;
{ if (lf->base==NULL) return(0.0);
    return(lf->base[i]);
}

double cens(lf,i)
lfit *lf;
INT i;
{ if (lf->c==NULL) return(0.0);
    return(lf->c[i]);
}

double vocri(lk,t0,t2,pen)
double lk, t0, t2, pen;
{ if (pen==0) return(-2*t0*lk/((t0-t2)*(t0-t2)));
    return((-2*lk+pen*t2)/t0);
}

INT procvraw(des,lf,v)
design *des;
lfit *lf;
INT v;
{ INT lf_status;
    int i;
    double h, coef[1+MXDIM];
    des->xev = evpt(lf,v);
    
    lf_status = ainitband(des,lf);
    
    if (!lf_error) switch(lf->mi[MACRI])
    { case AKAT:
        case ACP:
        case AMDI:
            h = aband2(des,lf,des->h);
            h = aband3(des,lf,h);
            h = nbhd(lf,des,0,h,1);
            lf_status = locfit(lf,des,h,0);
            break;
        case ANONE:
        case AOK:
            break;
    }
    
    lf->h[v] = des->h;
    for (i=0; i<des->ncoef; i++) coef[i] = des->cf[cfn(des,i)];
    
    if (!lf_error)
    { if (lf->mi[MDC]) dercor(des,lf,coef);
        subparcomp(des,lf,coef);
        for (i=0; i<des->ncoef; i++) lf->coef[i*lf->nvm+v] =  coef[i];
    }
    
    lf->deg[v] = lf->mi[MDEG];
    
    return(lf_status);
}

/*
 * Set default values for the likelihood e.t.c. This
 * is called in cases where the optimization for the fit
 * has failed.
 */

void set_default_like(lf,nvm,v,d)
lfit *lf;
INT nvm, v;
int d;
{ INT i;
    lf->lik[v] = lf->lik[nvm+v] = 0;
    lf->lik[2*nvm+v] = 0; /* should use sum of weights here? */
    for (i=0; i<=d; i++)
        lf->t0[i*nvm+v] = lf->nlx[i*nvm+v] = 0.0;
}

INT procv(des,lf,v)
design *des;
lfit *lf;
INT v;
{ INT d, p, nvm, i, k;
    double trc[6], t0[1+MXDIM], vari[1+MXDIM];
    memset(vari, 0, sizeof(vari));
    k = procvraw(des,lf,v);
    if (lf_error) return(k);
    
    d = lf->mi[MDIM]; p = lf->mi[MP];
    nvm = lf->nvm;
    
    switch(k)
    { case LF_OK: break;
        case LF_NCON:
            WARN(("procv: locfit did not converge"));
            break;
        case LF_OOB:
            WARN(("procv: parameters out of bounds"));
            break;
        case LF_PF:
            if (lf->mi[MDEB]>1) WARN(("procv: perfect fit"));
            set_default_like(lf,nvm,v,d);
            return(k);
        case LF_NOPT:
            WARN(("procv: no points with non-zero weight"));
            set_default_like(lf,nvm,v,d);
            return(k);
        case LF_INFA:
            if (lf->mi[MDEB]>1) WARN(("procv: initial value problem"));
            set_default_like(lf,nvm,v,d);
            return(k);
        case LF_DEMP:
            WARN(("procv: density estimate, empty integration region"));
            set_default_like(lf,nvm,v,d);
            return(k);
        case LF_XOOR:
            WARN(("procv: fit point outside xlim region"));
            set_default_like(lf,nvm,v,d);
            return(k);
        case LF_DNOP:
            if (lf->mi[MDEB]>1)
                WARN(("density estimation -- insufficient points in smoothing window"));
            set_default_like(lf,nvm,v,d);
            return(k);
        case LF_FPROB:
            WARN(("procv: f problem; likelihood failure"));
            set_default_like(lf,nvm,v,d);
            return(k);
        default:
            WARN(("procv: unknown return code %d",k));
            set_default_like(lf,nvm,v,d);
            return(k);
    }
    
    comp_vari(lf,des,trc,t0);
    lf->lik[v] = des->llk;
    lf->lik[nvm+v] = trc[2];
    lf->lik[2*nvm+v] = trc[0]-trc[2];
    lf->nlx[v] = sqrt(des->V[0]);
    
    for (i=0; i<des->ncoef; i++)
        vari[i] = des->V[p*cfn(des,0) + cfn(des,i)];
    vari[0] = sqrt(vari[0]);
    if (vari[0]>0) for (i=1; i<des->ncoef; i++) vari[i] /= vari[0];
    t0[0] = sqrt(t0[0]);
    if (t0[0]>0) for (i=1; i<des->ncoef; i++) t0[i] /= t0[0];
    
    subparcomp2(des,lf,vari,t0);
    for (i=0; i<des->ncoef; i++)
    { lf->nlx[i*lf->nvm+v] = vari[i];
        lf->t0[i*lf->nvm+v]  = t0[i];
    }
    
    return(k);
}

double intvo(des,lf,c0,c1,a,p,t0,t20,t21)
design *des;
lfit *lf;
double *c0, *c1, a, t0, t20, t21;
INT p;
{ double th, lk, link[LLEN];
    INT i;
    lk = 0;
    for (i=0; i<des->n; i++)
    { th = (1-a)*innerprod(c0,&des->X[i*p],p) + a*innerprod(c1,&des->X[i*p],p);
        stdlinks(link,lf,des->ind[i],th,robscale);
        lk += des->w[i]*link[ZLIK];
    }
    des->llk = lk;
    return(vocri(des->llk,t0,(1-a)*t20+a*t21,lf->dp[DADP]));
}

INT procvvord(des,lf,v)
design *des;
lfit *lf;
INT v;
{ 
    static const int x_1 = 4;
    static const int y_1 = 10;
    double tr[6], gcv, g0, ap, coef[x_1][y_1], t2[4], th, md = 0.0;
    INT i, j, k = 0, d1, *mi, i0, p1, ip;
    mi = lf->mi;
    des->xev = evpt(lf,v);
    
    for (i = 0; i < x_1; ++i)
    {
        for (j = 0; j < y_1; ++j)
        {
            coef[i][j] = 0;
        }
    }
    
    lf->h[v] = nbhd(lf,des,(INT)(mi[MN]*lf->dp[DALP]),lf->dp[DFXH],0);
    if (lf->h[v]<=0) WARN(("zero bandwidth in procvvord"));
    
    ap = lf->dp[DADP];
    if ((ap==0) & ((mi[MTG]&63)!=TGAUS)) ap = 2.0;
    d1 = mi[MDEG]; p1 = mi[MP];
    for (i=0; i<p1; i++) coef[0][i] = coef[1][i] = coef[2][i] = coef[3][i] = 0.0;
    i0 = 0; g0 = 0;
    ip = 1;
    for (i=mi[MDEG0]; i<=d1; i++)
    { mi[MDEG] = i; des->p = mi[MP] = calcp(mi,i);
        k = locfit(lf,des,lf->h[v],0);
        
        local_df(lf,des,tr);
        gcv = vocri(des->llk,tr[0],tr[2],ap);
        if ((i==mi[MDEG0]) || (gcv<g0)) { i0 = i; g0 = gcv; md = i; }
        
        for (j=0; j<des->p; j++) coef[i][j] = des->cf[j];
        t2[i] = tr[2];
        
#ifdef RESEARCH
        printf("variable order\n");
        if ((ip) && (i>mi[MDEG0]))
        { for (j=1; j<10; j++)
        { gcv = intvo(des,lf,coef[i-1],coef[i],j/10.0,des->p,tr[0],t2[i-1],t2[i]);
            if (gcv<g0) { g0 = gcv; md = i-1+j/10.0; }
        }
        }
#endif
    }
    
    if (i0<d1) /* recompute the best fit */
    { mi[MDEG] = i0; des->p = mi[MP] = calcp(mi,i0);
        k = locfit(lf,des,lf->h[v],0);
        for (i=mi[MP]; i<p1; i++) des->cf[i] = 0.0;
        i0 = (INT)md; if (i0==d1) i0--;
        th = md-i0;
        for (i=0; i<p1; i++) des->cf[i] = (1-th)*coef[i0][i]+th*coef[i0+1][i];
        mi[MDEG] = d1; mi[MP] = p1;
    }
    
    for (i=0; i<p1; i++) lf->coef[i*lf->nvm+v] = des->cf[i];
    lf->deg[v] = md;
    return(k);
}

/* special version of ressumm to estimate sigma^2, with derivative estimation */
void ressummd(lf,des)
lfit *lf;
design *des;
{ INT i;
    double s0, s1;
    s0 = s1 = 0.0;
    if ((lf->mi[MTG]&64)==0)
    { lf->dp[DRV] = 1.0;
        return;
    }
    for (i=0; i<lf->nv; i++)
    { s0 += lf->lik[2*lf->nvm+i];
        s1 += lf->lik[i];
    }
    if (s0==0.0)
        lf->dp[DRV] = 0.0;
    else
        lf->dp[DRV] = -2*s1/s0;
}

void ressumm(lf,des)
lfit *lf;
design *des;
{ INT i, j, ev, tg, orth;
    double *dp, *oy, pw, r1, r2, rdf, t0, t1, u[MXDIM], link[LLEN];
    dp = lf->dp;
    dp[DLK] = dp[DT0] = dp[DT1] = 0;
    if ((lf->mi[MEV]==EKDCE) | (lf->mi[MEV]==EPRES))
    { dp[DRV] = 1.0;
        return;
    }
    if (lf->nd>0)
    { ressummd(lf,des);
        return;
    }
    r1 = r2 = 0.0;
    ev = lf->mi[MEV];
    if ((ev==EDATA) | (ev==ECROS)) ev = EFITP;
    orth = (lf->mi[MGETH]==4) | (lf->mi[MGETH]==5);
    for (i=0; i<lf->mi[MN]; i++)
    { for (j=0; j<lf->mi[MDIM]; j++) u[j] = datum(lf,j,i);
        des->th[i] = base(lf,i)+dointpoint(lf,des,u,PCOEF,ev,i);
        des->wd[i] = resp(lf,i) - des->th[i];
        des->w[i] = 1.0;
        des->ind[i] = i;
    }
    
    tg = lf->mi[MTG];
    lf->dp[DRSC] = 1.0;
    if ((tg==TROBT+64) | (tg==TCAUC+64)) /* global robust scale */
    { oy = lf->y; lf->y = des->wd;
        des->xev = lf->pc.xbar;
        locfit(lf,des,0.0,1);
        lf->y = oy;
        lf->dp[DRSC] = robscale;
    }
    
    if (orth) /* orthog. residuals */
    { int od, op;
        des->n = lf->mi[MN];
        od = lf->mi[MDEG]; op = lf->mi[MP];
        lf->mi[MDEG] = 1;
        lf->mi[MP] = des->p = 1+lf->mi[MDIM];
        oy = lf->y; lf->y = des->wd;
        des->xev = lf->pc.xbar;
        locfit(lf,des,0.0,1);
        for (i=0; i<lf->mi[MN]; i++) oy[i] = resp(lf,i) - des->th[i];
        lf->y = oy;
        lf->mi[MDEG] = od; lf->mi[MP] = op;
    }
    
    for (i=0; i<lf->mi[MN]; i++)
    { for (j=0; j<lf->mi[MDIM]; j++) u[j] = datum(lf,j,i);
        t0 = dointpoint(lf,des,u,PT0,ev,i);
        t1 = dointpoint(lf,des,u,PNLX,ev,i);
        stdlinks(link,lf,i,des->th[i],lf->dp[DRSC]);
        t1 = t1*t1*link[ZDDLL];
        t0 = t0*t0*link[ZDDLL];
        if (t1>1) t1 = 1;
        if (t0>1) t0 = 1; /* no observation gives >1 deg.free */
        dp[DLK] += link[ZLIK];
        dp[DT0] += t0;
        dp[DT1] += t1;
        pw = prwt(lf,i);
        if (pw>0)
        { r1 += link[ZDLL]*link[ZDLL]/pw;
            r2 += link[ZDDLL]/pw;
        }
        if (orth) des->di[i]  = t1;
    }
    
    if (orth) return;
    
    dp[DRV] = 1.0;
    if ((lf->mi[MTG]&64)==64) /* quasi family */
    { rdf = lf->mi[MN]-2*dp[DT0]+dp[DT1];
        if (rdf<1.0)
        { WARN(("Estimated rdf < 1.0; not estimating variance"));
        }
        else
            dp[DRV] = r1/r2 * lf->mi[MN] / rdf;
    }
    
    /* try to ensure consistency for family="circ"! */
    if (((lf->mi[MTG]&63)==TCIRC) & (lf->mi[MDIM]==1))
    { INT *ind, nv;
        double dlt, th0, th1;
        ind = des->ind;
        nv = lf->nv;
        for (i=0; i<nv; i++) ind[i] = i;
        lforder(ind,vdptr(lf->xxev),0,nv-1);
        for (i=1; i<nv; i++)
        { dlt = evptx(lf,ind[i],0)-evptx(lf,ind[i-1],0);
            th0 = lf->coef[ind[i]]-dlt*lf->coef[ind[i]+nv]-lf->coef[ind[i-1]];
            th1 = lf->coef[ind[i]]-dlt*lf->coef[ind[i-1]+nv]-lf->coef[ind[i-1]];
            if ((th0>PI)&(th1>PI))
            { for (j=0; j<i; j++)
                lf->coef[ind[j]] += 2*PI;
                i--;
            }
            if ((th0<(-PI))&(th1<(-PI)))
            { for (j=0; j<i; j++)
                lf->coef[ind[j]] -= 2*PI;
                i--;
            }
        }
    }
}

double rss(lf,des,df)
lfit *lf;
design *des;
double *df;
{ double ss;
    INT i;
    ss = 0;
    if (ident==1)
    { for (i=0; i<lf->mi[MN]; i++)
        ss += SQR(resp(lf,i)-lf->coef[i]);
        *df = lf->mi[MN]-lf->mi[MP];
        return(ss);
    }
    ressumm(lf,des);
    *df = lf->mi[MN] - 2*lf->dp[DT0] + lf->dp[DT1];
    return(-2*lf->dp[DLK]);
}
