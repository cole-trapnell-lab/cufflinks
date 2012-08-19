/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   Compute minimax weights for local regression.
 */

#include "local.h"

int mmsm_ct;

static int debug=0;
#define CONVTOL 1.0e-8
#define SINGTOL 1.0e-10
#define NR_SINGULAR 100

static lfit *mm_lf;
static design *mm_des;
static double mm_gam;

double ipower(x,n) /* use for n not too large!! */
double x;
int n;
{ if (n==0) return(1.0);
    if (n<0) return(1/ipower(x,-n));
    return(x*ipower(x,n-1));
}

double setmmwt(des,lf,a,gam)
design *des;
lfit *lf;
double *a, gam;
{ double ip, w0, w1, sw, wt;
    INT i, p;
    sw = 0.0;
    p = lf->mi[MP];
    for (i=0; i<lf->mi[MN]; i++)
    { ip = innerprod(a,d_xi(des,i),p);
        wt = prwt(lf,i);
        w0 = ip - gam*des->wd[i];
        w1 = ip + gam*des->wd[i];
        des->w[i] = 0.0;
        if (w0>0) { des->w[i] = w0; sw += wt*w0*w0; }
        if (w1<0) { des->w[i] = w1; sw += wt*w1*w1; }
    }
    return(sw/2-a[0]);
}

/* compute sum_{w!=0} AA^T; e1-sum wA  */
int mmsums(coef,f,z,J)
double *coef, *f, *z;
jacobian *J;
{ int i, j, p, sing;
    double *A;
    
    mmsm_ct++;
    A = J->Z;
    *f = setmmwt(mm_des,mm_lf,coef,mm_gam);
    
    p = mm_lf->mi[MP];
    setzero(A,p*p);
    setzero(z,p);
    z[0] = 1.0;
    
    for (i=0; i<mm_lf->mi[MN]; i++)
        if (mm_des->w[i]!=0.0)
        { addouter(A,d_xi(mm_des,i),d_xi(mm_des,i),p,prwt(mm_lf,i));
            for (j=0; j<p; j++) z[j] -= prwt(mm_lf,i)*mm_des->w[i]*mm_des->X[i*p+j];
        }
    
    J->st = JAC_RAW;
    jacob_dec(J,JAC_EIGD);
    
    sing = 0;
    for (i=0; i<p; i++) sing |= (J->Z[i*p+i]<SINGTOL);
    if ((debug) & (sing)) printf("SINGULAR!!!!\n");
    
    /* printf("%8.5f %8.5f  %8.5f   f %8.5f  z %8.5f %8.5f\n",
     coef[0],coef[1],mm_gam,*f,z[0],z[1]); */
    return((sing) ? NR_SINGULAR : NR_OK);
}

double updatesd(des,lf,z,p,a,a0,sw0,gam)
design *des;
lfit *lf;
int p;
double *z, *a, *a0, sw0, gam;
{ double f, sw, c0, c1, tmp[10];
    INT i, j, sd;
    
    if (debug) printf("updatesd\n");
    for (i=0; i<p; i++) if (des->xtwx.Z[i*p+i]<SINGTOL) sd = i;
    if (des->xtwx.dg[sd]>0)
        for (i=0; i<p; i++) tmp[i] = des->xtwx.Q[p*i+sd]*des->xtwx.dg[i];
    else
    { for (i=0; i<p; i++) tmp[i] = 0.0;
        tmp[sd] = 1.0;
    }
    
    mmsums(a0,&sw,z,&des->xtwx);
    
    c0 = c1 = 0.0;
    for (i=0; i<p; i++)
    { c0 += tmp[i]*z[i];
        for (j=0; j<p; j++)
            c1 += tmp[i]*des->xtwx.Z[i*p+j]*tmp[j];
    }
    if (debug) printf("sdir: c0 %8.5f  c1 %8.5f  z %8.5f %8.5f  tmp %8.5f %8.5f\n",c0,c1,z[0],z[1],tmp[0],tmp[1]);
    if (c0<0) for (i=0; i<p; i++) tmp[i] = -tmp[i];
    
    f = 1.0;
    for (i=0; i<p; i++) a[i] = a0[i]+tmp[i];
    sw = setmmwt(des,lf,a,gam);
    
    if (sw<sw0) /* double till we drop */
    { while(1)
    { f *= 2;
        sw0 = sw;
        for (i=0; i<p; i++) a[i] = a0[i]+f*tmp[i];
        sw = setmmwt(des,lf,a,gam);
        if (sw>sw0-CONVTOL) /* go back one step */
        { f /= 2;
            for (i=0; i<p; i++) a[i] = a0[i]+f*tmp[i];
            sw0 = setmmwt(des,lf,a,gam);
            return(sw0);
        }
    }
    }
    
    /* halve till success */
    while (1)
    { f *= 0.5;
        for (i=0; i<p; i++) a[i] = a0[i]+f*tmp[i];
        sw = setmmwt(des,lf,a,gam);
        if (sw<sw0+CONVTOL) return(sw);
    }
}

int mm_initial(des,lf,z,p,coef)
design *des;
lfit *lf;
int p;
double *z, *coef;
{ int st;
    double f;
    
    setzero(coef,p);
    coef[0] = 1;
    while (1)
    {
        st = mmsums(coef,&f,z,&des->xtwx);
        if (st==NR_OK) return(0);
        coef[0] *= 2;
        if (coef[0]>1e8) return(1);
    } 
}

int mmax(coef, old_coef, f1, delta, J, p, maxit, tol, err)
double *coef, *old_coef, *f1, *delta, tol;
int p, maxit, *err;
jacobian *J;
{ double f, old_f, lambda;
    int i, j, fr, sing;
    
    *err = NR_OK;
    J->p = p;
    J->st = JAC_RAW;
    fr = mmsums(coef,&f,f1,J);
    
    for (j=0; j<maxit; j++)
    { memcpy(old_coef,coef,p*sizeof(double));
        old_f = f;
        
        /* compute delta = Newton-Raphson increment */
        /* jacob_dec(J,JAC_EIGD); */
        
        sing = (fr==NR_SINGULAR);
        
        if (fr == NR_SINGULAR)
        { J->st = JAC_RAW;
            if (j==0) printf("init singular\n");
            f = updatesd(mm_des,mm_lf,delta,p,coef,old_coef,f,mm_gam);
            fr = mmsums(coef,&f,f1,J);
        }
        else
        { 
            jacob_solve(J,f1);
            memcpy(delta,f1,p*sizeof(double));
            /* printf("delta %8.5f %8.5f\n",f1[0],f1[1]); */
            lambda = 1.0;
            do
            {
                for (i=0; i<p; i++) coef[i] = old_coef[i] + lambda*delta[i];
                J->st = JAC_RAW;
                fr = mmsums(coef,&f,f1,J);
                
                lambda = lambda/2.0;
                /* if (fr==NR_SINGULAR) printf("singular\n"); */
            } while (((lambda>0.000000001) & (f > old_f+0.001)) /* | (fr==NR_SINGULAR) */ );
            
            if (f>old_f+0.001) { printf("lambda prob\n"); *err = NR_NDIV; return(f); }
            
        }
        if (f==0.0)
        { if (sing) printf("final singular - conv\n");
            return(f);
        }
        
        if (debug)
        { for (i=0; i<p; i++) printf("%8.5f ",coef[i]);
            printf(" f %8.5f\n",f);
        }
        
        if ((j>0) & (fabs(f-old_f)<tol)) return(f);
    }
    if (sing) printf("final singular\n");
    WARN(("findab not converged"));
    *err = NR_NCON;
    return(f);
}

double findab(double gam)
{ double *coef, sl;
    int i, p, nr_stat;
    
    mm_gam = gam;
    p = mm_lf->mi[MP];
    
    /* starting values for nr iteration */
    coef = mm_des->cf;
    for (i=0; i<p; i++) coef[i] = 0.0;
    if (mm_initial(mm_des, mm_lf, mm_des->f1, p, coef))
    { WARN(("findab: initial value divergence"));
        return(0.0);
    }
    else
        mmax(coef, mm_des->oc, mm_des->res, mm_des->f1,
             &mm_des->xtwx, p, mm_lf->mi[MMXIT], CONVTOL, &nr_stat);
    
    if (nr_stat != NR_OK) return(0.0);
    
    sl = 0.0;
    for (i=0; i<mm_lf->mi[MN]; i++) sl += fabs(mm_des->w[i])*mm_des->wd[i];
    
    return(sl-gam);
}

double weightmm(coef,di,ff,mi,gam)
double *coef, di, *ff, gam;
INT *mi;
{ double y1, y2, ip;
    ip = innerprod(ff,coef,mi[MP]);
    y1 = ip-gam*di; if (y1>0) return(y1/ip);
    y2 = ip+gam*di; if (y2<0) return(y2/ip);
    return(0.0);
}

double minmax(lf,des)
lfit *lf;
design *des;
{ double h, u[MXDIM], gam;
    int i, j, m, d1, p1, err_flag;
    
    mmsm_ct = 0;
    d1 = lf->mi[MDEG]+1;
    p1 = factorial(d1);
    for (i=0; i<lf->mi[MN]; i++)
    { for (j=0; j<lf->mi[MDIM]; j++) u[j] = datum(lf,j,i);
        des->wd[i] = lf->dp[DALP]/p1*ipower(des->di[i],d1);
        des->ind[i] = i;
        fitfun(lf,u,des->xev,d_xi(des,i),NULL,(INT)0);
    }
    /* designmatrix(lf,des); */
    
    /* find gamma (i.e. solve eqn 13.17 from book), using the secant method.
     * As a side effect, this finds the other minimax coefficients.
     * First, set some global pointers a, mm_lf, mm_des.
     * Note that 13.17 is rewritten as
     *   g2 = sum |l_i(x)| (||xi-x||^(p+1) M/(s*(p+1)!))
     * where g2 = gamma * s * (p+1)! / M. The gam variable below is g2.
     * The smoothing parameter is lf->dp[DALP] == M/s.
     */
    mm_lf = lf;
    mm_des = des;
    gam = solve_secant(findab, 0.0, 0.0,1.0, 0.0000001, BDF_EXPRIGHT, &err_flag);
    
    /*
     * Set the smoothing weights, in preparation for the actual fit.
     */
    h = 0.0; m = 0;
    for (i=0; i<lf->mi[MN]; i++)
    { des->w[m] = weightmm(des->cf, des->wd[i],&des->X[i*lf->mi[MP]],lf->mi,gam);
        if (des->w[m]>0)
        { if (des->di[i]>h) h = des->di[i];
            des->ind[m] = i;
            m++;
        }
    }
    des->n = m;
    return(h);
}
