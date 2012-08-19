/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

/*
  Functions for computing residuals and fitted values from
  the locfit object.

  fitted(lf,des,fit,what,cv,ty) computes fitted values from the
    fit structure in lf. 
  resid(y,c,w,th,mi,ty) converts fitted values to residuals
  cfitted(v,ty) is CVERSION front end, interpreting command
    line arguments and computing th.
  vfitted(ty) for use by arithmetic interpreter.
*/

#include "local.h"

double resid(y,w,th,mi,ty,res)
INT *mi, ty;
double y, w, th, *res;
{ double raw;
  INT tg;

  tg = mi[MTG] & 63;
  if ((tg==TGAUS) | (tg==TROBT) | (tg==TCAUC))
    raw = y-res[ZMEAN];
  else
    raw = y-w*res[ZMEAN];
  switch(ty)
  { case RDEV:
      if (res[ZDLL]>0) return(sqrt(-2*res[ZLIK]));
            else return(-sqrt(-2*res[ZLIK]));
    case RPEAR:
      if (res[ZDDLL]<=0)
      { if (res[ZDLL]==0) return(0);
        return(NOSLN);
      }
      return(res[ZDLL]/sqrt(res[ZDDLL]));
    case RRAW:  return(raw);
    case RLDOT: return(res[ZDLL]);
    case RDEV2: return(-2*res[ZLIK]);
    case RLDDT: return(res[ZDDLL]);
    case RFIT:  return(th);
    case RMEAN: return(res[ZMEAN]);
    default: ERROR(("resid: unknown residual type %d",ty));
  }
  return(0.0);
}

double studentize(res,inl,var,ty,link)
double res, inl, var, *link;
int ty;
{ double den;
  inl *= link[ZDDLL];
  var = var*var*link[ZDDLL];
  if (inl>1) inl = 1;
  if (var>inl) var = inl;
  den = 1-2*inl+var;
  if (den<0) return(0.0);
  switch(ty)
  { case RDEV:
    case RPEAR:
    case RRAW:
    case RLDOT:
      return(res/sqrt(den));
    case RDEV2:
      return(res/den);
    default: return(res);
  }
}

void fitted(lf,des,fit,what,cv,st,ty)
lfit *lf;
design *des;
double *fit;
INT what, cv, st, ty;
{ INT i, j, d, n, ev;
  double xx[MXDIM], th, inl, var, link[LLEN];
    inl = 0.0;
    var = 0.0;
  n = lf->mi[MN];
  d = lf->mi[MDIM];
  ev = lf->mi[MEV];
  cv &= (ev!=ECROS);
  if ((lf->mi[MEV]==EDATA)|(lf->mi[MEV]==ECROS)) ev = EFITP;
  for (i=0; i<n; i++)
  { for (j=0; j<d; j++) xx[j] = datum(lf,j,i);
    th = dointpoint(lf,des,xx,what,ev,i);
    if ((what==PT0)|(what==PVARI)) th = th*th;
    if (what==PCOEF)
    { th += base(lf,i);
      stdlinks(link,lf,i,th,lf->dp[DRSC]);
      if ((cv)|(st))
      { inl = dointpoint(lf,des,xx,PT0,ev,i);
        inl = inl*inl;
        if (cv)
        { th -= inl*link[ZDLL];
          stdlinks(link,lf,i,th,lf->dp[DRSC]);
        }
        if (st) var = dointpoint(lf,des,xx,PNLX,ev,i);
      }
      fit[i] = resid(resp(lf,i),prwt(lf,i),th,lf->mi,ty,link);
      if (st) fit[i] = studentize(fit[i],inl,var,ty,link);
    } else fit[i] = th;
    if (lf_error) return;
  }
}

#ifdef CVERSION
extern lfit lf;
extern design des;

vari *vfitted(type)
INT type;
{ vari *v;
  INT n;
  n = lf.mi[MN];
  v = createvar("vfitted",STHIDDEN,n,VDOUBLE);
  recondat(1,&n);
  if (lf_error) return(NULL);

  fitted(&lf,&des,vdptr(v),PCOEF,0,0,type);
  return(v);
}

void cfitted(v,ty)
vari *v;
INT ty;
{ double *f;
  vari *vr;
  INT i, n, cv, st, wh;

  i = getarg(v,"type",1);
  if (i>0) ty = restyp(argval(v,i));

  i = getarg(v,"cv",1); cv = (i>0) ? getlogic(v,i) : 0;
  i = getarg(v,"studentize",1); st = (i>0) ? getlogic(v,i) : 0;

  wh = PCOEF;
  i = getarg(v,"what",1);
  if (i>0) wh = ppwhat(argval(v,i));

  recondat(ty==5,&n);
  if (lf_error) return;

  vr = createvar("fitted",STHIDDEN,n,VDOUBLE);
  f = vdptr(vr);
  fitted(&lf,&des,f,wh,cv,st,ty);

  saveresult(vr,argarg(v,0),STREGULAR);
}
#endif
