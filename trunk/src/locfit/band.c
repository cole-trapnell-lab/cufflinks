/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

extern lfit lf;
extern design des;
extern void fitoptions();

static double hmin, gmin, sig2, pen, vr, tb;

#define BGCV 1
#define BCP  2
#define BIND 3

INT procvbind(des,lf,v)
design *des;
lfit *lf;
INT v;
{ double s0, s1, bi;
    int i, ii, k;
    k = procvraw(des,lf,v);
    wdiag(lf,des,des->wd,0,1,0);
    s0 = s1 = 0.0;
    for (i=0; i<des->n; i++)
    { ii = des->ind[i];
        s0+= prwt(lf,ii)*des->wd[i]*des->wd[i];
        bi = prwt(lf,ii)*fabs(des->wd[i]*ipower(des->di[ii],lf->mi[MDEG]+1));
        s1+= bi*bi;
    }
    vr += s0;
    tb += s1;
    return(k);
}

double bcri(h,c,cri)
double h;
INT c, cri;
{ double num, den;
    INT (*pv)();
    lf.dp[c] = h;
    if ((cri&63)==BIND)
    { pv = procvbind;
        vr = tb = 0.0;
    }
    else pv = procv;
    if (cri<64) startlf(&des,&lf,pv,0);
    switch(cri&63)
    { case BGCV:
            ressumm(&lf,&des);
            num = -2*lf.mi[MN]*lf.dp[DLK];
            den = lf.mi[MN]-lf.dp[DT0];
            return(num/(den*den));
        case BCP:
            ressumm(&lf,&des);
            return(-2*lf.dp[DLK]/sig2-lf.mi[MN]+pen*lf.dp[DT0]);
        case BIND:
            return(vr+pen*pen*tb);
    } 
    ERROR(("bcri: unknown criterion"));
    return(0.0);
}

void bsel2(h0,g0,ifact,c,cri)
double h0, g0, ifact;
INT c, cri;
{ INT done, inc;
    double h1, g1;
    h1 = h0; g1 = g0;
    done = inc = 0;
    while (!done)
    { h1 *= 1+ifact;
        g0 = g1;
        g1 = bcri(h1,c,cri);
        if (g1<gmin) { hmin = h1; gmin = g1; }
        if (g1>g0) inc++; else inc = 0;
        switch(cri)
        { case BIND:
                done = (inc>=4) & (vr<lf.nv);
                break;
            default:
                done = (inc>=4);
        }
    }
}

void bsel3(h0,g0,ifact,c,cri)
double h0, g0, ifact;
INT c, cri;
{ double h1, g1;
    INT i;
    hmin = h0; gmin = g0;
    for (i=-1; i<=1; i++) if (i!=0)
    { h1 = h0*(1+i*ifact);
        g1 = bcri(h1,c,cri);
        if (g1<gmin) { hmin = h1; gmin = g1; }
    }
    return;
}

void bselect(c,cri,pn)
INT c, cri;
double pn;
{ double h0, g0, ifact;
    INT i;
    pen = pn;
    if (cri==BIND) pen /= factorial((int)lf.mi[MDEG]+1);
    hmin = h0 = lf.dp[c];
    if (h0==0) ERROR(("bselect: initial bandwidth is 0"));
    if (lf_error) return;
    sig2 = 1.0;
    
    gmin = g0 = bcri(h0,c,cri);
    if (cri==BCP)
    { sig2 = lf.dp[DRV];
        g0 = gmin = bcri(h0,c,cri+64);
    }
    
    ifact = 0.3;
    bsel2(h0,g0,ifact,c,cri);
    
    for (i=0; i<5; i++)
    { ifact = ifact/2;
        bsel3(hmin,gmin,ifact,c,cri);
    }
    lf.dp[c] = hmin;
    startlf(&des,&lf,procv,0);
    ressumm(&lf,&des);
}

double compsda(x,h,n)
double *x, h;
INT n;
/* n/(n-1) * int( fhat''(x)^2 dx ); bandwidth h */
{ INT i, j;
    double ik, sd, z;
    ik = wint(1,NULL,0,WGAUS);
    sd = 0;
    
    for (i=0; i<n; i++)
        for (j=i; j<n; j++)
        { z = (x[i]-x[j])/h;
            sd += (2-(i==j))*Wconv4(z,WGAUS)/(ik*ik);
        }
    sd = sd/(n*(n-1)*h*h*h*h*h);
    return(sd);
}

double widthsj(x,lambda,n)
double *x, lambda;
INT n;
{ double ik, a, b, td, sa, z, c, c1, c2, c3;
    INT i, j;
    a = GFACT*0.920*lambda*exp(-log((double)n)/7)/SQRT2;
    b = GFACT*0.912*lambda*exp(-log((double)n)/9)/SQRT2;
    ik = wint(1,NULL,0,WGAUS);
    
    td = 0;
    for (i=0; i<n; i++)
        for (j=i; j<n; j++)
        { z = (x[i]-x[j])/b;
            td += (2-(i==j))*Wconv6(z,WGAUS)/(ik*ik);
        }
    
    td = -td/(n*(n-1));
    j = 2.0;
    c1 = Wconv4(0.0,WGAUS);
    c2 = wint(1,&j,1,WGAUS);
    c3 = Wconv(0.0,WGAUS);  /* (2*c1/(c2*c3))^(1/7)=1.357 */
    sa = compsda(x,a,n);
    c = b*exp(log(c1*ik/(c2*c3*GFACT*GFACT*GFACT*GFACT)*sa/td)/7)*SQRT2;
    return(c);
}

void kdecri(x,h,res,c,k,ker,n)
double *x, h, *res, c;
INT k, ker, n;
{ INT i, j;
    double degfree, dfd, pen, s, r0, r1, d0, d1, ik, wij;
    
    if (h<=0) WARN(("kdecri, h = %6.4f",h));
    
    res[0] = res[1] = 0.0;
    ik = wint(1,NULL,0,ker);
    switch(k)
    { case 1: /* aic */
            pen = 2.0;
            for (i=0; i<n; i++)
            { r0 = d0 = 0.0;
                for (j=0; j<n; j++)
                { s = (x[i]-x[j])/h;
                    r0 += W(s,ker);
                    d0 += s*s*Wd(s,ker);
                }
                d0 = -(d0+r0)/(n*h*h*ik);  /* d0 = d/dh fhat(xi) */
                r0 /= n*h*ik;              /* r0 = fhat(xi) */
                res[0] += -2*log(r0)+pen*W(0.0,ker)/(n*h*ik*r0);
                res[1] += -2*d0/r0-pen*W(0.0,ker)/(n*h*ik*r0)*(d0/r0+1.0/h);
            }
            return;
        case 2: /* ocv */
            for (i=0; i<n; i++)
            { r0 = 0.0; d0 = 0.0;
                for (j=0; j<n; j++) if (i!=j)
                { s = (x[i]-x[j])/h;
                    r0 += W(s,ker);
                    d0 += s*s*Wd(s,ker);
                }
                d0 = -(d0+r0)/((n-1)*h*h);
                r0 = r0/((n-1)*h);
                res[0] -= log(r0);
                res[1] -= d0/r0;
            }
            return;
        case 3: /* lscv */
            r0 = r1 = d0 = d1 = degfree = 0.0;
            for (i=0; i<n; i++)
            { dfd = 0.0;
                for (j=0; j<n; j++)
                { s = (x[i]-x[j])/h;
                    wij = W(s,ker);
                    dfd += wij;
                    /* 
                     *  r0 = \int fhat * fhat = sum_{i,j} W*W( (Xi-Xj)/h ) / n^2 h.
                     *  d0 is it's derivative wrt h.
                     *
                     *  r1 = 1/n sum( f_{-i}(X_i) ).
                     *  d1 is  it's derivative wrt h.
                     *
                     *  degfree = sum_i ( W_0 / sum_j W( (Xi-Xj)/h ) ) is fitted d.f.
                     */
                    r0 += Wconv(s,ker);
                    d0 += -s*s*Wconv1(s,ker);
                    if (i != j)
                    { r1 += wij;
                        d1 += -s*s*Wd(s,ker);
                    }
                }
                degfree += 1.0/dfd;
            }
            d1 -= r1;
            d0 -= r0;
            res[0] = r0/(n*n*h*ik*ik)   - 2*r1/(n*(n-1)*h*ik);
            res[1] = d0/(n*n*h*h*ik*ik) - 2*d1/(n*(n-1)*h*h*ik);
            res[2] = degfree;
            return;
        case 4: /* bcv */
            r0 = d0 = 0.0;
            for (i=0; i<n; i++)
                for (j=i+1; j<n; j++)
                { s = (x[i]-x[j])/h;
                    r0 += 2*Wconv4(s,ker);
                    d0 += 2*s*Wconv5(s,ker);
                }
            d0 = (-d0-r0)/(n*n*h*h*ik*ik);
            r0 = r0/(n*n*h*ik*ik);
            j = 2.0;
            d1 = wint(1,&j,1,ker);
            r1 = Wconv(0.0,ker);
            res[0] = (d1*d1*r0/4+r1/(n*h))/(ik*ik);
            res[1] = (d1*d1*d0/4-r1/(n*h*h))/(ik*ik);
            return;
        case 5: /* sjpi */
            s = c*exp(5*log(h)/7)/SQRT2;
            d0 = compsda(x,s,n);
            res[0] = d0; /* this is S(alpha) in SJ */
            res[1] = exp(log(Wikk(WGAUS,0)/(d0*n))/5)-h;
            return;
        case 6: /* gas-k-k */
            s = exp(log(1.0*n)/10)*h;
            d0 = compsda(x,s,n);
            res[0] = d0;
            res[1] = exp(log(Wikk(WGAUS,0)/(d0*n))/5)-h;
            return;
    }
    ERROR(("kdecri: what???"));
    return;
}

double esolve(x,j,h0,h1,k,c,ker,n)
double *x, h0, h1, c;
INT j, k, ker, n;
{ 
    double h[7], d[7], r[7], res[4], min, minh, fact;
    INT i, nc;
    memset(h, 0, sizeof(h));
    memset(d, 0, sizeof(d));
    memset(r, 0, sizeof(r));
    
    min = 1.0e30; minh = 0.0;
    fact = 1.00001;
    h[6] = h0; kdecri(x,h[6],res,c,j,ker,n);
    r[6] = res[0]; d[6] = res[1];
    if (lf_error) return(0.0);
    nc = 0;
    for (i=0; i<k; i++)
    { h[5] = h[6]; r[5] = r[6]; d[5] = d[6];
        h[6] = h0*exp((i+1)*log(h1/h0)/k);
        kdecri(x,h[6],res,c,j,ker,n);
        r[6] = res[0]; d[6] = res[1];
        if (lf_error) return(0.0);
        if (d[5]*d[6]<0)
        { h[2] = h[0] = h[5]; d[2] = d[0] = d[5]; r[2] = r[0] = r[5];
            h[3] = h[1] = h[6]; d[3] = d[1] = d[6]; r[3] = r[1] = r[6];
            while ((h[3]>fact*h[2])|(h[2]>fact*h[3]))
            { h[4] = h[3]-d[3]*(h[3]-h[2])/(d[3]-d[2]);
                if ((h[4]<h[0]) | (h[4]>h[1])) h[4] = (h[0]+h[1])/2;
                kdecri(x,h[4],res,c,j,ker,n);
                r[4] = res[0]; d[4] = res[1];
                if (lf_error) return(0.0);
                h[2] = h[3]; h[3] = h[4];
                d[2] = d[3]; d[3] = d[4];
                r[2] = r[3]; r[3] = r[4];
                if (d[4]*d[0]>0) { h[0] = h[4]; d[0] = d[4]; r[0] = r[4]; }
                else { h[1] = h[4]; d[1] = d[4]; r[1] = r[4]; }
            }
            if (j>=4) return(h[4]); /* first min for BCV etc */
            if (r[4]<=min) { min = r[4]; minh = h[4]; }
            nc++;
        }
    }
    if (nc==0) minh = (r[5]<r[6]) ? h0 : h1;
    return(minh);
}

void kdeselect(band,x,ind,h0,h1,meth,nm,ker,n)
double h0, h1, *band, *x;
INT *ind, *meth, nm, ker, n;
{ double scale, c;
    INT i, k;
    k = n/4;
    for (i=0; i<n; i++) ind[i] = i;
    scale = kordstat(x,n+1-k,n,ind) - kordstat(x,k,n,ind);
    c = widthsj(x,scale,n);
    for (k=0; k<nm; k++)
        band[k] = esolve(x,meth[k],h0,h1,10,c,ker,n);
}

#ifdef CVERSION
void band(v)
vari *v;
{ INT i, c, cri;
    double pen;
    char *z;
    fitoptions(&lf,v,0);
    
    c = DALP;
    i = getarg(v,"comp",1);
    if (i>0)
    { z = argval(v,i);
        if (z[0]=='h') c = DFXH;
    }
    
    cri = BGCV;
    i = getarg(v,"bcri",1);
    if (i>0)
    { z = argval(v,i);
        if (z[0]=='c') cri = BCP;
        if (z[0]=='i') cri = BIND;
    }
    
    pen = 2.0;
    i = getarg(v,"pen",1);
    if (i>0)
        pen = darith(argval(v,i));
    
    bselect(c,cri,pen);
}
#endif

#ifdef SVERSION
void slscv(x,n,h,z)
double *x, *h, *z;
int *n;
{ INT i;
    double res[4];
    kdecri(x,*h,res,0.0,3,WGAUS,*n);
    z[0] = res[0];
    z[1] = res[2];
}
#endif
