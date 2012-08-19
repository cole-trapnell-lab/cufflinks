/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 startlf(des,lf,vfun,nopc) -- starting point for locfit.
 des and lf are pointers to the design and fit structures.
 vfun is the vertex processing function.
 nopc=1 inhibits computation of parametric component.
 fitdefault(lf,n,d) -- fit default parameters.
 lf is pointer to fit; n and d are sample size and dimension.
 deschk()  -- assignment function for the design structure.
 preproc() -- fit preprocessing (limits, scales, paramcomp etc.)
 bbox()    -- compute bounding box.
 
 fitoptions()
 clocfit()  -- start point for CLocfit - interpret cmd line etc.
 */

#include "local.h"

extern INT cvi;

void fitdefault(lf,n,d)
lfit *lf;
INT n, d;
{ INT i, *mi;
    double *dp;
    
    dp = lf->dp;
    mi = lf->mi;
    
    mi[MTG] = TNUL;
    mi[MTG] = (lf->y==NULL) ? TDEN : (64+TGAUS);
    mi[MLINK] = LDEFAU;
    mi[MACRI] = ANONE;
    mi[MDEG] = mi[MDEG0] = 2;
    mi[MEV] = (ident==1) ? EDATA : ETREE;
    mi[MKT] = KSPH; mi[MKER] = WTCUB;
    mi[MIT] = IDEFA; mi[MDC] = mi[MREN] = 0;
    mi[MK] = 100; mi[MMINT] = 20;
    mi[MMXIT] = 20;
    mi[MN] = n; mi[MDIM] = d;
    mi[MDEB] = 0;
    mi[MUBAS] = 0;
    
    
    dp[DALP] = 0.7; dp[DFXH] = dp[DADP] = 0.0;
    dp[DCUT] = 0.8;
    
    if (d<=0)
        ERROR(("must set MDIM before calling fitdefault"));
    for (i=0; i<d; i++)
    { lf->sca[i] = 1.0;
        lf->xl[i] = lf->xl[i+d] = 0.0;
        lf->fl[i] = lf->fl[i+d] = 0.0;
    }
}

int des_reqd(n,p)
INT n, p;
{ 
    return (n*(p+5)+2*p*p+4*p + jac_reqd(p));
}
int des_reqi(INT n) { return(n); }

void deschk(des,n,p)
design *des;
INT n, p;
{ 
    double *z;
    des->dw = checkvarlen(des->dw,des_reqd(n,p),"_deswork",VDOUBLE);
    z = vdptr(des->dw);
    des->X = z; z += n*p;
    setzero(des->X, n*p);
    
    des->w = z; z += n;
    setzero(des->w, n);
    
    des->res=z; z += n;
    setzero(des->res, n);
    
    des->di =z; z += n;
    setzero(des->di, n);
    
    des->th =z; z += n;
    setzero(des->th, n);
    
    des->wd =z; z += n;
    setzero(des->wd, n);
    
    des->V  =z; z += p*p;
    setzero(des->V, p*p);
    
    des->P  =z; z += p*p;
    setzero(des->P, p*p);
    
    des->f1 =z; z += p;
    setzero(des->f1, p);
    
    des->ss =z; z += p;
    setzero(des->ss, p);
    
    des->oc =z; z += p;
    setzero(des->oc, p);
    
    des->cf =z; z += p;
    setzero(des->cf, p);
    
    z = jac_alloc(&des->xtwx,p,z);
    
    des->index = checkvarlen(des->index,des_reqi(n),"_desidx",VINT);
    des->ind = (INT *)vdptr(des->index);
    des->n = n; 
    des->p = p;
    des->xtwx.p = p;
}

void bbox(lf,bx)
lfit *lf;
double *bx;
{ INT i, j, d, n;
    double z, mx, mn;
    d = lf->mi[MDIM]; n = lf->mi[MN];
    for (i=0; i<d; i++)
        if (bx[i]==bx[i+d])
        { if (lf->sty[i]==STANGL)
        { bx[i] = 0.0; bx[i+d] = 2*PI*lf->sca[i];
        }
        else
        { mx = mn = datum(lf,i,0);
            for (j=1; j<n; j++)
            { mx = MAX(mx,datum(lf,i,j));
                mn = MIN(mn,datum(lf,i,j));
            }
            if (lf->xl[i]<lf->xl[i+d]) /* user set xlim; maybe use them. */
            { z = mx-mn;
                if (mn-0.2*z < lf->xl[i]) mn = lf->xl[i];
                if (mx+0.2*z > lf->xl[i+d]) mx = lf->xl[i+d];
            }
            bx[i] = mn;
            bx[i+d] = mx;
        }
        }
}

void preproc(des,lf,nopc)
design *des;
lfit *lf;
INT nopc;
{ INT d, i, j, n;
    double xb;
    d = lf->mi[MDIM]; n = lf->mi[MN];
    lf->mi[MLINK] = defaultlink(lf->mi[MLINK],lf->mi[MTG]);
    if (!validlinks(lf->mi[MLINK],lf->mi[MTG]))
    { ERROR(("Invalid family/link combination"));
        return;
    }
    compparcomp(des,lf,nopc);
    if (lf->w==NULL)
        lf->dp[DSWT] = lf->mi[MN];
    else
    { lf->dp[DSWT] = 0;
        for (i=0; i<lf->mi[MN]; i++) lf->dp[DSWT] += prwt(lf,i);
    }
    for (i=0; i<d; i++)
        if (lf->sca[i]<=0) /* set automatic scales */
        { if (lf->sty[i]==STANGL) lf->sca[i] = 1.0;
        else
        { xb = lf->sca[i] = 0.0;
            for (j=0; j<n; j++) xb += datum(lf,i,j);
            xb /= n;
            for (j=0; j<n; j++) lf->sca[i] += SQR(datum(lf,i,j)-xb);
            lf->sca[i] = sqrt(lf->sca[i]/(n-1));
        }
        }
    bbox(lf,lf->fl);
}

#ifdef CVERSION
extern void do_scbsim();
#endif

void startlf(des,lf,vfun,nopc)
design *des;
lfit *lf;
INT (*vfun)(), nopc;
{ 
    INT i, *mi;
    des->vfun = vfun;
    mi = lf->mi;
    mi[MP] = calcp(mi,mi[MDEG]);
    des->pref = 0;
    cvi = -1; /* inhibit cross validation */
    deschk(des,mi[MN],mi[MP]);
    if (mi[MDEB]>0) printf("preprocess\n");
        preproc(des,lf,nopc);
        if (mi[MDEB]>0) printf("preprocess ok\n");
            if (lf_error) return;
    lf->ord = 0;
    makecfn(des,lf);
    if ((mi[MDIM]==1) && (lf->sty[0]!=STANGL))
    { i = 1;
        while ((i<mi[MN]) && (datum(lf,0,i)>=datum(lf,0,i-1))) i++;
        lf->ord = (i==mi[MN]);
    }
    
    if (mi[MDEB]>0) printf("call eval structure\n");
        switch(mi[MEV])
    { case EPHULL: triang_start(des,lf); break;
        case EDATA:  dataf(des,lf); break;
        case ECROS:  crossf(des,lf); break;
        case EGRID:  gridf(des,lf); break;
        case ETREE:  atree_start(des,lf); break;
        case EKDCE:  mi[MKT] = KCE;
        case EKDTR:  kdtre_start(des,lf); break;
        case EPRES:  preset(des,lf); break;
        case EXBAR:  xbarf(des,lf); break;
        case ENONE:  lf->nv = lf->nce = 0;
            return;
#ifdef CVERSION
        case 100: do_scbsim(des,lf); break;
#endif
        default: ERROR(("startlf: Invalid evaluation structure"));
    }
    
    /* renormalize for family=density */
    if ((mi[MREN]) && (mi[MTG]==TDEN)) dens_renorm(lf,des);
        }

#ifdef CVERSION
extern lfit lf;
extern design des;
extern plots pl[];
int curwin;
vari *vb;

INT nofit()
{ if (lf.mi==NULL) return(1);
    return(lf.mi[MEV]==ENULL);
}

void endfit()
{ INT i;
    for (i=0; i<MAXWIN; i++)
        if (pl[i].track != NULL)
        { curwin = i;
            cmdint(pl[i].track);
        }
}

INT drl(key,dv,mi)
char *key;
INT *dv, *mi;
{ INT i, nd;
    nd = readilist(dv,key,0,mi[MDEG],0);
    for (i=0; i<nd; i++)
    { if ((dv[i]<1) | (dv[i]>mi[MDIM]))
        ERROR(("drl: Invalid derivatives %s",key));
        dv[i]--;
    }
    return(nd);
}

void fitoptions(lf,vc,re)
lfit *lf;
vari *vc;
INT re;
{ 
    INT d = 0, n, i, i0, i1, *mi;
    char kc, *key;
    vari *v;
    
    re &= (!nofit());
    i0 = getarg(vc,"formula",1);
    if ((!re) && (i0==0)) { ERROR(("no formula")); return; }
    i1 = getarg(vc,"data",1);
    if (i1>0) doreaddata(argval(vc,i1),(INT)0);
    if (re)
        recondat(0,&lf->mi[MN]);
    else
    { lf->base = lf->y = lf->c = lf->w = NULL;
        lf->nd = 0;
        strcpy(lf->yname,"_NuLl");
        strcpy(lf->wname,"_NuLl");
        strcpy(lf->bname,"_NuLl");
        strcpy(lf->cname,"_NuLl");
    }
    if (i0>0) /* interpret formula */
    { key = argval(vc,i0);
        n = -1;
        i0 = i1 = 0; d = 0;
        while ((i0<strlen(key)) && (key[i0]!='~')) i0++;
        if (key[i0] != '~') { ERROR(("invalid formula %s",key)); return; }
        if (i0>0)
        { key[i0] = '\0';
            lf->y = vdptr(findvar(key,1,&n));
            strcpy(lf->yname,key);
            key[i0] = '~';
        }
        i1 = i0 = i0+1;
        while (i1<strlen(key))
        { while ((i1<strlen(key)) && (key[i1]!='+')) i1++;
            kc = key[i1]; key[i1] = '\0';
            lf->sty[d] = KPROD;
            if (stm(&key[i0],"left(",5))
            { lf->sty[d] = STLEFT;
                i0 = i0+5; key[i1-1] = '\0';
            }
            else if (stm(&key[i0],"right(",6))
            { lf->sty[d] = STRIGH;
                i0 = i0+6; key[i1-1] = '\0';
            }
            else if (stm(&key[i0],"ang(",4))
            { lf->sty[d] = STANGL;
                i0 = i0+4; key[i1-1] = '\0';
            }
            else if (stm(&key[i0],"cpar(",5))
            { lf->sty[d] = STCPAR;
                i0 = i0+5; key[i1-1] = '\0';
            }
            dvari(lf,d) = vdptr(findvar(&key[i0],1,&n));
            strcpy(lf->xname[d],&key[i0]);
            if (lf->sty[d]!=KPROD) key[i1-1] = ')';
            d++; key[i1] = kc;
            i0 = i1 = i1+1;
        }
        fitdefault(lf,n,d);
    }
    mi = lf->mi;
    
    i = getarg(vc,"weights",1);
    if (i>0)
    { lf->w = vdptr(findvar(argval(vc,i),1,&mi[MN]));
        strcpy(lf->wname,argval(vc,i));
    }
    i = getarg(vc,"cens",1);
    if (i>0)
    { lf->c = vdptr(findvar(argval(vc,i),1,&mi[MN]));
        strcpy(lf->cname,argval(vc,i));
    }
    i = getarg(vc,"base",1);
    if (i>0)
    { lf->base = vdptr(findvar(argval(vc,i),1,&mi[MN]));
        strcpy(lf->bname,argval(vc,i));
    }
    
    i = getarg(vc,"scale",1);
    if (i>0)
    { if (argvalis(vc,i,"T"))
        for (i=0; i<d; i++) lf->sca[i] = 0;
    else if (argvalis(vc,i,"F"))
        for (i=0; i<d; i++) lf->sca[i] = 1;
    else
        arvect(argval(vc,i),lf->sca,d,0);
    }
    
    i = getarg(vc,"vb",0);
    if (i>0)
    { lf->dp[DALP] = -1;
        vb = arbuild(argval(vc,i),0,strlen(argval(vc,i))-1,NULL,0,1);
        setvarname(vb,"_varband");
    }
    else
    { i = getarg(vc,"alpha",1);
        if (i>0) arvect(argval(vc,i),&lf->dp[DALP],3,1);
    }
    
    i = getarg(vc,"deg",1);
    if (i>0)
    { i =  readilist(&mi[MDEG0],argval(vc,i),1,2,0);
        if (i==1) mi[MDEG] = mi[MDEG0];
    }
    
    i = getarg(vc,"family",1);if (i>0) setstrval(mi,MTG,argval(vc,i));
    i = getarg(vc,"link",1);  if (i>0) setstrval(mi,MLINK,argval(vc,i));
    i = getarg(vc,"ev",1); 
    if (i>0)
    { v = findvar(argval(vc,i),0,NULL);
        if (v!=NULL)
        { mi[MEV] = EPRES;
            lf->xxev= v;
            lf->nvm = v->n;
        }
        else
            setstrval(mi,MEV,argval(vc,i));
    }
    i = getarg(vc,"acri",1);  if (i>0) setstrval(mi,MACRI,argval(vc,i));
    
    i = getarg(vc,"mg",1);
    if (i>0) readilist(lf->mg,argval(vc,i),1,MXDIM,1);
    
    i = getarg(vc,"kt",1);   if (i>0) setstrval(mi,MKT, argval(vc,i));
    i = getarg(vc,"kern",1); if (i>0) setstrval(mi,MKER,argval(vc,i));
    i = getarg(vc,"itype",1);if (i>0) setstrval(mi,MIT, argval(vc,i));
    
    i = getarg(vc,"cut",1);
    if (i>0) lf->dp[DCUT] = darith(argval(vc,i));
    
    i = getarg(vc,"flim",1);
    if (i>0) arvect(argval(vc,i),lf->fl,2*d,2);
    
    i = getarg(vc,"xlim",1);
    if (i>0) arvect(argval(vc,i),lf->xl,2*d,2);
    
    i = getarg(vc,"deriv",0);
    if (i>0) lf->nd = drl(argval(vc,i),lf->deriv,lf->mi);
    i = getarg(vc,"dc",1); if (i>0) mi[MDC] = getlogic(vc,i);
    i = getarg(vc,"maxk",1); if (i>0) readilist(&mi[MK],argval(vc,i),1,1,0);
    i = getarg(vc,"mint",1); if (i>0) readilist(&mi[MMINT],argval(vc,i),1,1,0);
    i = getarg(vc,"maxit",1); if (i>0) readilist(&mi[MMXIT],argval(vc,i),1,1,0);
    i = getarg(vc,"renorm",1);if (i>0) mi[MREN] = getlogic(vc,i);
    i = getarg(vc,"debug",1); if (i>0) readilist(&mi[MDEB],argval(vc,i),1,1,0);
}

void clocfit(v,re)
INT re;
vari *v;
{
    lf.ord = 0;
    lf.kap[0] = lf.kap[1] = lf.kap[2] = 0.0; lf.nk = 0;
    fitoptions(&lf,v,re);
    if (lf_error)
    { if (lf.mi!=NULL) lf.mi[MEV] = ENULL;
        return;
    }
    
    
    lf.nv = 0;
    if (lf.mi[MDEG0]==lf.mi[MDEG])
    { startlf(&des,&lf,procv,0);
        if (!lf_error) ressumm(&lf,&des);
    }
    else
        startlf(&des,&lf,procvvord,0);
    if (lf_error)
    { if (!re) lf.mi[MEV] = ENULL;
        return;
    }
    
    //printf("Evaluation structure %d, %d points.\n",lf.mi[MEV],lf.nv);
    if (argarg(v,0) != NULL) dosavefit(&lf,argarg(v,0),"wb",(INT)0);
    endfit();
}

#endif
