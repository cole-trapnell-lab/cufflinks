/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

/*
 preplot():  interpolates the fit to a new set of points.
 lf  -- the fit structure.
 des -- design structure; really only needed for parfit.
 x   -- the points to predict at.
 f   -- vector to return the predictions.
 se  -- vector to return std errors (NULL if not req'd)
 band-- char for conf band type. ('n'=none, 'g'=global etc.)
 n   -- no of predictions (or vector of margin lengths for grid)
 where -- where to predict:
 1 = points in the array x.
 2 = grid defined by margins in x.
 3 = data points from lf (ignore x).
 4 = fit points from lf (ignore x).
 what -- what to predict.
 (PCOEF etc; see lfcons.h file)
 
 cpreplot(): C version front end.
 setpppoints(): (C version) used to set preplot points.
 */

static char cb;
double *sef, *fit, sigmahat;

void predptall(lf,des,x,what,ev,i)
lfit *lf;
design *des;
double *x;
INT what, ev, i;
{ double lik, rdf;
    fit[i] = dointpoint(lf,des,x,what,ev,i);
    if (cb=='n') return;
    sef[i] = dointpoint(lf,des,x,PNLX,ev,i);
    if (cb=='g')
    { sef[i] *= sigmahat;
        return;
    }
    if (cb=='l')
    { lik = dointpoint(lf,des,x,PLIK,ev,i);
        rdf = dointpoint(lf,des,x,PRDF,ev,i);
        sef[i] *= sqrt(-2*lik/rdf);
        return;
    }
    if (cb=='p')
    { sef[i] = sigmahat*sqrt(1+sef[i]*sef[i]);
        return;
    }
}

void prepvector(lf,des,x,n,what) /* interpolate a vector */
lfit *lf;
design *des;
double **x;
INT n, what;
{ INT i, j;
    double xx[MXDIM];
    for (i=0; i<n; i++)
    { for (j=0; j<lf->mi[MDIM]; j++) xx[j] = x[j][i];
        predptall(lf,des,xx,what,lf->mi[MEV],i);
        if (lf_error) return;
    }
}

void prepfitp(lf,des,what)
lfit *lf;
design *des;
INT what;
{ 
    INT i;
    //d = lf->mi[MDIM];
    for (i=0; i<lf->nv; i++)
    { predptall(lf,des,evpt(lf,i),what,EFITP,i);
        if (lf_error) return;
    }
}

void prepgrid(lf,des,x,mg,n,what) /* interpolate a grid given margins */
design *des;
lfit *lf;
double **x;
INT *mg, n, what;
{ INT i, ii, j, d;
    double xv[MXDIM];
    d = lf->mi[MDIM];
    for (i=0; i<n; i++)
    { ii = i;
        for (j=0; j<d; j++)
        { xv[j] = x[j][ii%mg[j]];
            ii /= mg[j];
        }
        predptall(lf,des,xv,what,lf->mi[MEV],i);
        if (lf_error) return;
    }
}

void preplot(lf,des,x,f,se,band,mg,where,what)
lfit *lf;
design *des;
double **x, *f, *se;
INT *mg, where, what;
char band;
{ INT d = 0, i, n;
    double *xx[MXDIM];
    d = lf->mi[MDIM];
    fit = f;
    sef = se;
    cb = band;
    if (cb!='n') sigmahat = sqrt(lf->dp[DRV]);
    
    switch(where)
    { case 1: /* vector */
            n = mg[0];
            prepvector(lf,des,x,n,what);
            break;
        case 2: /* grid */
            n = 1;
            for (i=0; i<d; i++) n *= mg[i];
            prepgrid(lf,des,x,mg,n,what);
            break;
        case 3: /* data */
            n = lf->mi[MN];
            if ((lf->mi[MEV]==EDATA) | (lf->mi[MEV]==ECROS))
                prepfitp(lf,des,what);
            else
            { for (i=0; i<d; i++) xx[i] = dvari(lf,i);
                prepvector(lf,des,xx,n,what);
            }
            break;
        case 4: /* fit points */
            n = lf->nv;
            prepfitp(lf,des,what);
            break;
        default:
            ERROR(("unknown where in preplot"));
            return;
    }
    
    if ((what==PT0)|(what==PVARI))
        for (i=0; i<n; i++) f[i] = f[i]*f[i];
}

#ifdef CVERSION
extern lfit lf;
extern design des;

void cpreplot(pp,vc,band)
pplot *pp;
vari *vc;
char band;
{ double *data[MXDIM];
    INT j, mg[MXDIM];
    for (j=0; j<pp->d; j++)
    { data[j] = vdptr(pp->data[j]);
        mg[j] = pp->data[j]->n;
    }
    j = getarg(vc,"what",0);
    pp->wh = (j>0) ? ppwhat(argval(vc,j)) : PCOEF;
    
    preplot(&lf,&des,data,vdptr(pp->fit),vdptr(pp->se),band,mg,pp->gr,pp->wh);
}

INT setpppoints(pp,where,mg,xl)
pplot *pp;
char *where;
INT *mg;
double *xl;
{ INT d, i, j, n, m;
    varname vn;
    d = pp->d = lf.mi[MDIM];
    if (strcmp(where,"fitp")==0)
    { n = lf.nv;
        for (j=0; j<d; j++)
        { sprintf(vn,"_pred%d",j);
            pp->data[j] = createvar(vn,STPLOTVAR,n,VDOUBLE);
            if (lf_error) return(0);
            for (i=0; i<n; i++)
                vassn(pp->data[j],i,evptx(&lf,i,j));
        }
        pp->gr = 4;
        return(n);
    }
    if (strcmp(where,"data")==0)
    { recondat(1,&n);
        for (j=0; j<d; j++)
        { sprintf(vn,"_pred%d",j);
            pp->data[j] = createvar(vn,STPLOTVAR,n,VDOUBLE);
            if (lf_error) return(0);
            for (i=0; i<n; i++)
                vassn(pp->data[j],i,datum(&lf,j,i));
        }
        pp->gr = 3;
        return(n);
    }
    if (strcmp(where,"grid")==0)
    { n = 1;
        for (j=0; j<d; j++)
        { sprintf(vn,"_pred%d",j);
            m = (mg==NULL) ? 40 : mg[j];
            pp->data[j] = createvar(vn,STPLOTVAR,m,VDOUBLE);
            if (lf_error) return(0);
            if (m==1)
                vassn(pp->data[j],0,(xl[d+j]+xl[j])/2);
            else
                for (i=0; i<m; i++)
                    vassn(pp->data[j],i,xl[j]+i*(xl[d+j]-xl[j])/(m-1));
            n *= m;
            pp->gr = 2;
        }
        return(n);
    }
    ERROR(("setpppoints: invalid where=%s",where));
    return(0);
}

#endif
