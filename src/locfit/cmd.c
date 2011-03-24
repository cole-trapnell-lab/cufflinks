/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

#include <unistd.h>
#ifdef DOS
#include <dos.h>
#endif

#include "local.h"

#ifdef CVERSION

#define MAXK 20

FILE *ofile;

device devps, devwin;
design des;
lfit lf;
vari *aru;

extern plots pl[];
pplot pp;
struct lfcol mycol[MAXCOLOR];
char *lfhome;
extern char filename[100];

INT lf_error, lfcm[10];

vari *curstr;
void cmdint();
void del_lines();

/*
 INDEX  data input and output functions
 savefit:   f/end for user savefit.
 readdata:  f/end for readdata
 savedata:  f/end for savedata
 recondat:  reconnect data to fit
 */

void savefit(v,mode)
vari *v;
char *mode;
{ INT j, fp;
    char *filename;
    filename = getargval(v,"file",1);
    if (filename==NULL)
    { ERROR(("savefit: no filename"));
        return;
    }
    j = getarg(v,"fp",1);
    fp = (j>0) ? getlogic(v,j) : 0;
    dosavefit(&lf,filename,mode,fp);
    if (mode[0]=='r') endfit();
}

void readdata(v)
vari *v;
{ INT i, fp;
    i = getarg(v,"data",1);
    if (i==0) i = getarg(v,"file",1);
    if (i==0) { ERROR(("readdata: no file name")); return; }
    fp = getarg(v,"fp",1);
    fp = (fp>0) ? getlogic(v,fp) : 0;
    doreaddata(argval(v,i),fp);
}

void savedata(v)
vari *v;
{ INT fp;
    if (argarg(v,0)==NULL) { ERROR(("savedata: no file name")); return; }
    fp = getarg(v,"fp",0);
    fp = (fp>0) ? getlogic(v,fp) : 0;
    dosavedata(v,fp);
}

void recondat(xonly,n)
INT xonly, *n;
{ INT i;
    *n = -1;
    for (i=0; i<lf.mi[MDIM]; i++) dvari(&lf,i) = vdptr(findvar(lf.xname[i],1,n));
    if (lf_error | xonly) return;
    lf.y = vdptr(findvar(lf.yname,1,n));
    lf.c = vdptr(findvar(lf.cname,1,n));
    lf.w = vdptr(findvar(lf.wname,1,n));
    lf.base=vdptr(findvar(lf.bname,1,n));
}

/*
 INDEX  Call fitting functions e.t.c.
 ckap():       compute SCB constants.
 crband():     regression bandwidths
 ckdeb():      kde bandwidths
 */

void ckap(v)
vari *v;
{ INT i, nd;
    nd = 0;
    if (v->n==1) /* compute for existing fit */
    { if (nofit()) { ERROR(("ckap: no fit, no arguments")); }
    else recondat(0,&lf.mi[MN]);
    }
    else      /* new fit specification */
        fitoptions(&lf,v,0);
    if (lf_error) return;
    lf.nk = constants(&des,&lf,lf.kap);
    if (lf_error) { lf.nk=0; return; }
    printf("kappa0:");
    for (i=0; i<lf.nk; i++) printf(" %8.5f",lf.kap[i]);
    printf("\n");
}

void crband(v)
vari *v;
{ double h[4];
    INT i, kk, meth[4], nm;
    meth[0] = 1; meth[1] = 2; meth[2] = 3; meth[3] = 4;
    nm = 4;
    fitoptions(&lf,v,0);
    lf.mi[MDEG0] = lf.mi[MDEG]; lf.mi[MDEG] = 4;
    rband(&des,&lf,h,meth,&nm,&kk);
    for (i=0; i<nm; i++)
        printf("%8.5f ",h[i]);
    printf("\n");
    return;
}

void ckdeb(v)
vari *v;
{ INT i, mm[6], nm, n;
    double *x, band[6], h0, h1;
    char meth[6][5];
    strcpy(meth[0],"AIC");  strcpy(meth[1],"LCV");
    strcpy(meth[2],"LSCV"); strcpy(meth[3],"BCV");
    strcpy(meth[4],"SJPI"); strcpy(meth[5],"GKK");
    n = -1;
    i = getarg(v,"x",1);
    x = vdptr(findvar(argval(v,i),1,&n));
    if (lf_error) return;
    h0 = 0.02; h1 = 1.0;
    i = getarg(v,"h0",1); if (i>0) h0 = darith(argval(v,i));
    i = getarg(v,"h1",1); if (i>0) h1 = darith(argval(v,i));
    
    deschk(des,n,1);
    mm[0]=1; mm[1]=2; mm[2]=3; mm[3]=4; mm[4]=5; mm[5]=6; nm=6;
    kdeselect(band,x,des.ind,h0,h1,mm,nm,WGAUS,n);
    for (i=0; i<nm; i++)
        printf("%s: %8.6f ",meth[mm[i]-1],band[i]);
    printf("\n");
}

/*
 INDEX post-fitting functions
 docrit():    compute c for scb's.
 crit():      f/end to docrit();
 backtr():    back transform theta in likelihood models.
 predict():   interpolate the fit.
 printdata(): print the current dataset.
 printfit():  print the current fit.
 summdata():  summarize current dataset.
 summfit():   summarize current fit.
 */

double docrit(v)
vari *v;
{ double df, al;
    INT i;
    df = 0; al = 0.05;
    i = getarg(v,"df",1); if (i>0) sscanf(argval(v,i),"%lf",&df);
    i = getarg(v,"al",1); if (i>0) sscanf(argval(v,i),"%lf",&al);
    return(critval(lf.kap,lf.nk,lf.mi[MDIM],al,10,2,df));
}

void crit(v)
vari *v;
{ vari *vr;
    vr = createvar("crit",STHIDDEN,1,VDOUBLE);
    if (lf_error) return;
    vassn(vr,0,docrit(v));
    saveresult(vr,argarg(v,0),STREGULAR);
}

double backtr(th,mi,nd)
double th;
INT *mi, nd;
{ if (nd>0) return(th);
    return(invlink(th,mi[MLINK]));
}

void predict(vc)
vari *vc;
{ 
    double *data[MXDIM];
    varname vn;
    INT i, k, j, gr, n, z, mg[MXDIM];
    memset(mg, 0, sizeof(mg));
    dosavefit(&lf,getargval(vc,"fit",0),"rb",(INT)0);
    if (nofit()) ERROR(("predict: no fit to interpolate\n"));
    if (lf_error)  return;
    
    gr=0;
    i = getarg(vc,"grid",0);
    if (i>0) gr = getlogic(vc,i);
    
    i = getarg(vc,"where",0);
    if (i>0) n = setpppoints(&pp,argval(vc,i),NULL,lf.fl);
    else
    { 
        for (j=0; j<lf.mi[MDIM]; j++)
        { i = getarg(vc,lf.xname[j],1);
            if (i==0)
            { 
                ERROR(("predict: missing variables"));
                return;
            }
            if (gr) n = 0;
            sprintf(vn,"_pred%d",j);
            pp.data[j] = varith(argval(vc,i),vn,STPLOTVAR);
            if (lf_error) return;
            if (gr) mg[j] = pp.data[j]->n;
        }
        n = pp.data[0]->n;
        pp.gr = 1+gr;
    }
    
    for (j=0; j<lf.mi[MDIM]; j++) data[j] = vdptr(pp.data[j]);
    pp.d = lf.mi[MDIM];
    
    switch(pp.gr)
    { case 1:
            n = pp.data[0]->n;
            break;
        case 2:
            n = 1;
            for (i=0; i<lf.mi[MDIM]; i++)
            {
                mg[i] = pp.data[i]->n;
                n *= mg[i];
            }
            break;
        case 3:
            n = lf.mi[MN];
            break;
        case 4:
            n = lf.nv;
            break;
        default:
            ERROR(("cpreplot where problem"));
    }
    
    if (argarg(vc,0)==NULL)
        pp.fit = createvar("predict",STHIDDEN,n,VDOUBLE);
    else
        pp.fit = createvar(argarg(vc,0),STREGULAR,n,VDOUBLE);
    if (lf_error) return;
    pp.se = NULL;
    cpreplot(&pp,vc,'n');
    if (lf_error) return;
    for (j=0; j<n; j++)
        if (vitem(pp.fit,j)!=NOSLN)
            vassn(pp.fit,j,backtr(vitem(pp.fit,j),lf.mi,lf.nd));
    if (argarg(vc,0)!=NULL) return;
    for (j=0; j<n; j++)
    { for (i=0; i<lf.mi[MDIM]; i++)
    { z = j;
        if (pp.gr==2)
        { for (k=0; k<i; k++) z /= mg[k];
            z = z%mg[i];
        }
        //printf("%10.6f ",data[i][z]);
    }
        if (vitem(pp.fit,j)==NOSLN) printf("   Not computed\n");
        else printf("   %10.6f\n",vitem(pp.fit,j));
    }
    deletevar(pp.fit);
}

void printfit(v)
vari *v;
{ INT d, i = 0, j, k, cs, ck, nk, wh[MAXK];
    double rs, alp, c = 0, fh;
    cs = ck = nk = 0; rs = 1.0;
    
    dosavefit(&lf,getargval(v,"fit",i),"rb",(INT)0);
    for (i=1; i<v->n; i++) if (!argused(v,i))
    { if (argvalis(v,i,"x"))     { setused(v,i); wh[nk++]=1; }
        if (argvalis(v,i,"fhat"))  { setused(v,i); wh[nk++]=2; }
        if (argvalis(v,i,"coef"))  { setused(v,i); wh[nk++]=2; }
        if (argvalis(v,i,"nlx"))   { setused(v,i); wh[nk++]=3; }
        if (argvalis(v,i,"infl"))  { setused(v,i); wh[nk++]=4; }
        if (argvalis(v,i,"se"))    { setused(v,i); wh[nk++]=5; cs=1; }
        if (argvalis(v,i,"cband")) { setused(v,i); wh[nk++]=7; cs=ck=1; }
        if (argvalis(v,i,"h"))     { setused(v,i); wh[nk++]=8; }
        if (argvalis(v,i,"deg"))   { setused(v,i); wh[nk++]=9; }
    }
    if (nk==0) /* default: x and fhat */
    { wh[nk++] = 1; wh[nk++] = 2;
    }
    d = lf.mi[MDIM];
    alp = 0.95;
    
    if (cs) rs = sqrt(lf.dp[DRV]);
    if (ck)
    { c = critval(lf.kap,lf.nk,lf.mi[MDIM],1-alp,10,2,0.0);
        printf("using c = %8.5f\n",c);
    }
    
    for (i=0; i<lf.nv; i++) if (!lf.s[i])
    { fh = lf.coef[i]+addparcomp(&lf,evpt(&lf,i),PCOEF);
        for (j=0; j<nk; j++) switch(wh[j])
        { case 1:
                for (k=0; k<d; k++)
                    printf("%8.5f ",evptx(&lf,i,k));
                break;
            case 2: printf(" %12.6f ",backtr(fh,lf.mi,0)); break;
            case 3: printf(" %12.6f ",lf.nlx[i]); break;
            case 4: printf(" %12.6f ",lf.t0[i]); break;
            case 5: printf(" %12.6f ",rs*lf.nlx[i]); break;
            case 7: printf(" (%12.6f,%12.6f) ",fh-c*rs*lf.nlx[i],fh+c*rs*lf.nlx[i]);
                break;
            case 8: printf(" %12.6f ",lf.h[i]); break;
            case 9: printf(" %6.4f ",lf.deg[i]); break;
            default: ERROR(("prfit: what??"));
        }
        printf("\n");
    }
}

vari *knotsvar(name,n)
varname *name;
INT n;
{ vari *v;
    v = createvar("=knotv",STHIDDEN,n,VDOUBLE);
    if (lf_error) return(NULL);
    if (name!=NULL) v = saveresult(v,name,STREGULAR);
    return(v);
}

void knots(v)
vari *v;
{ INT i, j, k, n;
    vari *vr;
    if (nofit()) { ERROR(("knots: no fit")); return; }
    n = lf.nv; /* should delete pseudo vertices */
    for (k=0; k<v->n; k++)
    { vr = NULL;
        for (j=0; j<lf.mi[MDIM]; j++)
            if (argvalis(v,k,lf.xname[j]))
            { vr = knotsvar(argarg(v,k),n);
                for (i=0; i<n; i++) vassn(vr,i,evptx(&lf,i,j));
                setused(v,k);
            }
        if (argvalis(v,k,"fit")|argvalis(v,k,"coef"))
        { vr = knotsvar(argarg(v,k),n);
            for (i=0; i<n; i++) vassn(vr,i,backtr(lf.coef[i],lf.mi,0));
            setused(v,k);
        }
        if (argvalis(v,k,"h")|argvalis(v,k,"band"))
        { vr = knotsvar(argarg(v,k),n);
            for (i=0; i<lf.nv; i++) vassn(vr,i,lf.h[i]);
            setused(v,k);
        }
        if (argvalis(v,k,"deg"))
        { vr = knotsvar(argarg(v,k),n);
            for (i=0; i<lf.nv; i++) vassn(vr,i,lf.deg[i]);
            setused(v,k);
        }
        ((carg *)viptr(v,k))->result = vr;
    }
}

void summfit(v)
vari *v;
{ int i;
    dosavefit(&lf,getargval(v,"fit",1),"rb",0);
    printf("Response variable: %s\n",lf.yname);
    printf("Predictor variables: ");
    for (i=0; i<lf.mi[MDIM]; i++) printf("%s ",lf.xname[i]);
    printf("\nDegree of fit: %d\n",lf.mi[MDEG]);
    printf("Smoothing parameters: NN %f  fix %f  pen %f\n",
           lf.dp[DALP],lf.dp[DFXH],lf.dp[DADP]);
    printf("Fitting Family: ");
    switch(lf.mi[MTG]&63)
    { case TDEN: printf("Density Estimation\n"); break;
        case TRAT: printf("Poisson Process Rate Estimation\n"); break;
        case THAZ: printf("Hazard Rate Estimation\n"); break;
        case TGAUS:printf("Local Regression\n"); break;
        case TLOGT:printf("Binomial\n"); break;
        case TPOIS:printf("Poisson\n"); break;
        case TGAMM:printf("Exponential/Gamma\n"); break;
        case TGEOM:printf("Geometric/Negative Binomial\n"); break;
        case TCIRC:printf("Circular - Von Mises\n"); break;
    }
    printf("Fitted Degrees of Freedom: %8.5f\n",lf.dp[DT0]);
    printf("Number of fit points: %d\n",lf.nv);
    printf("Evaluation structure: ");
    switch(lf.mi[MEV])
    { case ENULL: printf("None\n"); break;
        case ETREE: printf("Rectangular tree\n"); break;
        case EPHULL:printf("Triangulation\n"); break;
        case EDATA: printf("Data\n"); break;
        case EGRID: printf("Grid\n"); break;
        case EKDTR: printf("K-d Tree\n"); break;
        case EKDCE: printf("K-d Tree (centers)\n"); break;
        case ECROS: printf("Data, Cross-Validation\n"); break;
        case EPRES: printf("User-provided\n"); break;
        default:    printf("Unknown\n");
    }
}

void AC(name,r,g,b,p)
char *name;
INT r, g, b, p;
{ devwin.AddColor(name,r,g,b,p);
    devps.AddColor(name,r,g,b,p);
}

INT getcolidx(cname, def)
char *cname;
int def;
{ int i;
    if (cname==NULL) return(def);
    for (i=0; i<8; i++)
        if (strcmp(cname,mycol[i].name)==0) return(i);
    WARN(("color %s not found",cname));
    return(def);
}

void greyscale(v)
vari *v;
{ INT i, j0, j1;
    j0 = getcolidx(getargval(v,"lo",1),0);
    j1 = getcolidx(getargval(v,"hi",1),1);
    for (i=0; i<=10; i++)
        AC("",((10-i)*mycol[j0].r+i*mycol[j1].r)/11,
           ((10-i)*mycol[j0].g+i*mycol[j1].g)/11,
           ((10-i)*mycol[j0].b+i*mycol[j1].b)/11,8+i);
}

void setcolor(v)
vari *v;
{
    return NULL;
//   int i;
//    lfcm[CBAK] = getcolidx(getargval(v,"back",0),lfcm[CBAK]);
//    
//    i = getarg(v,"fore",1);
//    if (i>0)
//    { lfcm[CAXI] = getcolidx(argval(v,i));
//        for (i=CTEX; i<CPA2; i++) lfcm[i] = lfcm[CAXI];
//    }
//    
//    lfcm[CAXI] = getcolidx(getargval(v,"axis",0),lfcm[CAXI]);
//    lfcm[CTEX] = getcolidx(getargval(v,"text",0),lfcm[CTEX]);
//    lfcm[CLIN] = getcolidx(getargval(v,"lines",0),lfcm[CLIN]);
//    lfcm[CPOI] = getcolidx(getargval(v,"points",0),lfcm[CPOI]);
//    lfcm[CCON] = getcolidx(getargval(v,"cont",0),lfcm[CCON]);
//    lfcm[CCLA] = getcolidx(getargval(v,"clab",0),lfcm[CCLA]);
//    lfcm[CSEG] = getcolidx(getargval(v,"cseg",0),lfcm[CSEG]);
//    lfcm[CPA1] = getcolidx(getargval(v,"patch1",0),lfcm[CPA1]);
//    lfcm[CPA2] = getcolidx(getargval(v,"patch2",0),lfcm[CPA2]);
//    if (lfcm[CAXI]==lfcm[0]) WARN(("axis color = background color"));
//    if (lfcm[CTEX]==lfcm[0]) WARN(("text color = background color"));
//    if (lfcm[CLIN]==lfcm[0]) WARN(("lines color = background color"));
//    if (lfcm[CPOI]==lfcm[0]) WARN(("points color = background color"));
//    if (lfcm[CCON]==lfcm[0]) WARN(("cont color = background color"));
//    if (lfcm[CCLA]==lfcm[0]) WARN(("clab color = background color"));
//    if (lfcm[CSEG]==lfcm[0]) WARN(("cseg color = background color"));
//    if (lfcm[CPA1]==lfcm[0]) WARN(("patch1 color = background color"));
//    if (lfcm[CPA2]==lfcm[0]) WARN(("patch2 color = background color"));
//    if (lfcm[CPA1]==lfcm[CPA2]) WARN(("patch1 color = patch2 color"));
}

void table(v)
vari *v;
{ INT i = 0, j = 0, ix, iy, m, mx, my, n, nx[15], ny[15], count[100];
    double xl[2], yl[2], xs[15], ys[15], *x, *y;
    i = getarg(v,"x",1);
    if (i==0)
    { ERROR(("table: no x variable"));
        return;
    }
    n = -1;
    x = vdptr(findvar(argval(v,i),1,&n));
    xl[0] = xl[1] = x[0];
    for (i=1; i<n; i++)
    { if (x[i]<xl[0]) xl[0] = x[i];
        if (x[i]>xl[1]) xl[1] = x[i];
    }
    i = getarg(v,"m",0);
    if (i>0) sscanf(argval(v,i),"%d",&m); else m = 5;
    mx = pretty(xl,m,xs);
    if (lf_error) return;
    
    i = getarg(v,"y",1);
    if (i>0)
    { y = vdptr(findvar(argval(v,i),1,&n));
        yl[0] = yl[1] = y[0];
        for (i=1; i<n; i++)
        { if (y[i]<yl[0]) yl[0] = y[i];
            if (y[i]>yl[1]) yl[1] = y[i];
        }
        my = pretty(yl,m,ys);
    }
    else { y = NULL; my = 0; }
    if (lf_error) return;
    
    for (i=0; i<15; i++) nx[i] = ny[i] = 0;
    for (i=0; i<=(mx+1)*(my+1); i++) count[i] = 0;
    for (i=0; i<n; i++)
    { if (x[i]<xs[0]) ix = 0;
        if (x[i]>=xs[mx-1]) ix = mx;
        if ((x[i]>=xs[0]) & (x[i]<xs[mx-1]))
            for (j=1; j<mx; j++)
                if ((x[i]>=xs[j-1]) & (x[i]<xs[j])) ix = j;
        if (my>0)
        { if (y[i]<ys[0]) iy = 0;
            if (y[i]>=ys[my-1]) iy = my;
            if ((y[i]>=ys[0]) & (y[i]<ys[my-1]))
                for (j=1; j<my; j++)
                    if ((y[i]>=ys[j-1]) & (y[i]<ys[j])) iy = j;
        } else iy = 0;
        nx[ix] = ny[iy] = 1;
        count[ix*(my+1)+iy]++;
    }
    if (my>0) printf("          ");
    for (i=0; i<=mx; i++) if (nx[i]>0)
        printf("  %4g-",(i==0) ? xl[0] : xs[i-1]);
    printf("\n");
    if (my>0) printf("          ");
    for (i=0; i<=mx; i++) if (nx[i]>0)
        printf("  %4g ",(i==mx) ? xl[1] : xs[i]);
    printf("\n\n");
    for (j=0; j<=my; j++) if (ny[j]>0)
    { if (my>0)
        printf("%4g-%4g ",(j==0) ? yl[0] : ys[j-1],
               (j==my) ? yl[1] : ys[j]);
        for (i=0; i<=mx; i++)
            if (nx[i]>0) printf("%6d ",count[i*(my+1)+j]);
        printf("\n");
    }
}

/*
 INDEX control functions:
 setout(): set output file.
 cmdint(): send off the command...
 locfit_dispatch(): called by the main program.
 */

void setout(v)
vari *v;
{ INT i, i0;
    char md;
    i0 = getarg(v,"file",1);
    if (i0==0)
    { if (ofile!=NULL) fclose(ofile);
        ofile = NULL;
        printf("Output set to stdout\n");
        return;
    }
    
    md = 'w';
    i = getarg(v,"mode",1);
    if ((i>0) && (argval(v,i)[0]=='a')) md = 'a';
    
    setfilename(argval(v,i0),"",&md,0);
    if (ofile != NULL) fclose(ofile);
    ofile = fopen(filename,&md);
    if (ofile == NULL)
        ERROR(("setout: can't open %s for writing",filename));
    else
        printf("Output set to file %s\n",filename);
}

void dosleep(v)
vari *v;
{ INT i;
    i = getarg(v,"time",1);
    if (i==0) return;
    sscanf(argval(v,i),"%d",&i);
    (void)sleep(i);
}

void setdef(v)
vari *v;
{ INT i, n;
    carg *ca;
    vari *vd;
    
    if (argarg(v,0)==NULL)
    { ERROR(("Unnamed Defintion"));
        return;
    }
    n = vlength(v)-1;
    vd = createvar(argarg(v,0),STSYSTEM,n,VARGL);
    if (lf_error) return;
    
    for (i=0; i<n; i++)
    { ca = (carg *)viptr(vd,i);
        ca->arg = argarg(v,i+1);
        ca->val = argval(v,i+1);
        setused(v,i+1);
    }
    sprintf(curstr->name,"=%s",argarg(v,0));
}

extern void cscbsim();

void dcmdint(v)
vari *v;
{ INT i;
    if (v==NULL)
    { ERROR(("dcmdint received NULL"));
        return;
    }
    if (argvalis(v,0,"band"))  { band(v); return; }
    if (argvalis(v,0,"crit"))  { crit(v); return; }
    if (argvalis(v,0,"def"))   { setdef(v); return; }
    if (argvalis(v,0,"endfor")) { dec_forvar(); return; }
    if (argvalis(v,0,"for"))    { inc_forvar(); return; }
    if (argvalis(v,0,"example")){example(v); return; }
    if (argvalis(v,0,"help"))   {example(v); return; }
    if (argvalis(v,0,"?"))      {example(v); return; }
    if (argvalis(v,0,"exit")) exit(0);
    if (argvalis(v,0,"quit")) exit(0);
    if (argvalis(v,0,"q()"))  exit(0);
    if (argvalis(v,0,"fitted")){ cfitted(v,RMEAN); return; }
    if (argvalis(v,0,"greyscale"))   { greyscale(v); return; }
    if (argvalis(v,0,"kappa")) { ckap(v); return; }
    if (argvalis(v,0,"kdeb"))  { ckdeb(v); return; }
    if (argvalis(v,0,"knots")) { knots(v); return; }
    if (argvalis(v,0,"locfit"))   { clocfit(v,0); return; }
    if (argvalis(v,0,"relocfit")) { clocfit(v,1); return; }
    if (argvalis(v,0,"plot"))     { printf("use plotfit or plotdata\n"); return; }
    if (argvalis(v,0,"plotdata")) { plotdata(v); return; }
    if (argvalis(v,0,"plotfit"))  { plotfit(v); return; }
    if (argvalis(v,0,"replot"))   { plotopt(v,1); return; }
    if (argvalis(v,0,"predict"))  { predict(v); return; }
    if (argvalis(v,0,"prfit"))    { printfit(v); return; }
    if (argvalis(v,0,"rband"))    { crband(v); return; }
    if (argvalis(v,0,"readdata")) { readdata(v); return; }
    if (argvalis(v,0,"readfile")) { readfile(v); return; }
    if (argvalis(v,0,"readfit"))  { savefit(v,"rb"); return; }
    if (argvalis(v,0,"residuals")){ cfitted(v,RDEV); return; }
    if (argvalis(v,0,"run"))      return;
    if (argvalis(v,0,"savedata")) { savedata(v); return; }
    if (argvalis(v,0,"savefit"))  { savefit(v,"wb"); return; }
    if (argvalis(v,0,"scbmax"))   { cscbsim(v); return; }
    if (argvalis(v,0,"scbsim"))   { cscbsim(v); return; }
    if (argvalis(v,0,"seed"))     { rseed(argval(v,1)); setused(v,1); return; }
    if (argvalis(v,0,"setcolor")) { setcolor(v); return; }
    if (argvalis(v,0,"setout"))   { setout(v); return; }
    if (argvalis(v,0,"outf"))     { setout(v); return; }
    if (argvalis(v,0,"setplot"))  { setplot(v); return; }
    if (argvalis(v,0,"sleep"))    { dosleep(v); return; }
    if (argvalis(v,0,"summfit"))  { summfit(v); return; }
    if (argvalis(v,0,"table"))    { table(v); return; }
    if (argvalis(v,0,"track"))    { plottrack(v); return; }
    if (argvalis(v,0,"wdiag"))    { cwdiag(v); return; }
    for (i=0; i<vlength(v); i++)
    { ((carg *)viptr(v,i))->result = varith(argval(v,i),argarg(v,i),STREGULAR);
        setused(v,i);
        if (lf_error) return;
    }
}

void cmdint(v)
vari *v;
{ vari *vv, *vr;
    INT i, j, mn, nr;
    if (v==NULL) return;
    
    for (i=0; i<vlength(v); i++)
    { setunused(v,i);
        /* ((carg *)viptr(v,i))->used = 0; */
        ((carg *)viptr(v,i))->result = NULL;
    }
    
    setused(v,0);
    if (vlength(v)==1)
    { j = 0;
        vv = findvar(argval(v,0),0,&j);
        if ((vv!=NULL) && ((vv->mode==VARGL) & (!argvalis(v,0,"=cline"))))
        { 
            cmdint(vv);
            return;
        }
    }
    
    /* dcmdint processes command */
    dcmdint(v);
    
    /* print the results of unassigned expression.
     * First, determine mn = maximum number of rows in the
     * output. Note that vr->stat==STHIDDEN determines whether
     * the result was unassigned.
     */
    mn = 0; nr = 0;
    for (i=0; i<vlength(v); i++)
    { vr = ((carg *)viptr(v,i))->result;
        if ((vr != NULL) && (vr->stat==STHIDDEN))
            switch(vr->mode)
        { case VCHAR: if (mn<1) mn = 1;
                break;
            case VINT:
            case VDOUBLE: if (mn<vr->n) mn = vr->n;
                break;
        }
    }
    
    /* now, print the unassigned variables.
     
     for (i=0; i<mn; i++)
     { for (j=0; j<vlength(v); j++)
     { vr = ((carg *)viptr(v,j))->result;
     if ((vr != NULL) && (vr->stat==STHIDDEN))
     switch(vr->mode)
     { case VDOUBLE: printf("%8.5f  ",vitem(vr,i)); break;
     case VCHAR:   printf("%s  ",vdptr(vr)); break;
     case VINT:    printf("%4d  ", vitem(vr,i)); break;
     }
     }
     printf("\n");
     }
     */
    
    for (i=0; i<vlength(v); i++)
        deleteifhidden(((carg *)viptr(v,i))->result);
}

INT locfit_dispatch(char *z)

{ vari *v;
    
    makecmd(z);
    while (1)
    { lf_error = 0;
        v = getcmd();
        if (v==NULL)
        { del_lines();
            return(0);
        }
        cmdint(v);
    }
}

void setuplf()
{ INT i;
    char command[100];
    vari *v;
    
    lfhome = getenv("LFHOME");
    initdb();
    
    ofile = NULL;
    lf.tw = lf.xxev = lf.L = lf.iw = des.dw = lf.pc.wk = NULL;
    des.index = NULL;
    lf.mg = calloc(MXDIM,sizeof(INT));
    
    v = createvar("mi",STSYSPEC,LENM,VINT); v->dpr = (double *)lf.mi;
    v = createvar("dp",STSYSPEC,LEND,VDOUBLE); v->dpr = lf.dp;
    v = createvar("alpha",STSYSPEC,1,VDOUBLE); v->dpr = &lf.dp[DALP];
    v = createvar("h",    STSYSPEC,1,VDOUBLE); v->dpr = &lf.dp[DFXH];
    v = createvar("pen",  STSYSPEC,1,VDOUBLE); v->dpr = &lf.dp[DADP];
    v = createvar("infl", STSYSPEC,1,VDOUBLE); v->dpr = &lf.dp[DT0];
    v = createvar("vari", STSYSPEC,1,VDOUBLE); v->dpr = &lf.dp[DT1];
    v = createvar("like", STSYSPEC,1,VDOUBLE); v->dpr = &lf.dp[DLK];
    v = createvar("resv", STSYSPEC,1,VDOUBLE); v->dpr = &lf.dp[DRV];
    
    for (i=0; i<MAXWIN; i++)
    { pl[i].xyzs = NULL;
        pl[i].id = i;
        pl[i].ty = PLNONE;
        pl[i].track = NULL;
    }
    //SetWinDev(&devwin);
    //SetPSDev(&devps);
    //  AC("white",255,255,255,0);
    //  AC("black",  0,  0,  0,1);
    //  AC(  "red",255,  0,  0,2);
    //  AC("green",  0,255,  0,3);
    //  AC( "blue",  0,  0,255,4);
    //  AC("magenta",255,0,255,5);
    //  AC("yellow",255,255, 0,6);
    //  AC( "cyan",  0,255,255,7);
    lfcm[0] = 0;
    for (i=CAXI; i<=CPA1; i++) lfcm[i] = 1;
    lfcm[CPA2] = 2;
    rseed("LocalFit");
    if (setfilename("LFInit","cmd","r",0))
    { sprintf(command,"run %s",filename);
        locfit_dispatch(command);
    }
}

#endif
