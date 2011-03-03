/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

#define HUBERC 2.0

extern double rs, log();

INT defaultlink(link,family)
INT link, family;
{ if (link==LDEFAU)
    switch(family&63)
    { case TDEN:
      case TRAT:
      case THAZ:
      case TGAMM:
      case TGEOM:
      case TPROB:
      case TPOIS: return(LLOG);
      case TCIRC:
      case TGAUS:
      case TCAUC:
      case TROBT: return(LIDENT);
      case TRBIN:
      case TLOGT: return(LLOGIT);
    }
  if (link==LCANON)
    switch(family&63)
    { case TDEN:
      case TRAT:
      case THAZ:
      case TPROB:
      case TPOIS: return(LLOG);
      case TGEOM:
        WARN(("Canonical link unavaialable for geometric family; using inverse"));
      case TGAMM: return(LINVER);
      case TCIRC:
      case TGAUS:
      case TCAUC:
      case TROBT: return(LIDENT);
      case TRBIN:
      case TLOGT: return(LLOGIT);
    }
  return(link);
}

INT validlinks(link,family)
INT link, family;
{ switch(family&63)
  { case TDEN:
    case TRAT:
    case THAZ:
      return((link==LLOG) | (link==LIDENT));
    case TGAUS:
      return((link==LIDENT) | (link==LLOG) | (link==LLOGIT));
    case TROBT:
    case TCAUC:
    case TCIRC:
      return(link==LIDENT);
    case TLOGT:
      return((link==LLOGIT) | (link==LIDENT) | (link==LASIN));
    case TRBIN:
      return(link==LLOGIT);
    case TGAMM:
      return((link==LLOG) | (link==LINVER) | (link==LIDENT));
    case TGEOM:
      return((link==LLOG) | (link==LIDENT));
    case TPOIS:
    case TPROB:
      return((link==LLOG) | (link==LSQRT) | (link==LIDENT));
  }
  ERROR(("Unknown family %d in validlinks",family));
  return(0);
}

INT famdens(mean,th,link,res,cens,w)
double mean, th, *res, w;
INT link, cens;
{ if (cens)
    res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
  else
  { res[ZLIK] = w*th;
    res[ZDLL] = res[ZDDLL] = w;
  }
  return(LF_OK);
}

INT famgaus(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ double z, pz, dp;
  if (link==LINIT)
  { res[ZDLL] = w*y;
    return(LF_OK);
  }
  z = y-mean;
  if (cens)
  { if (link!=LIDENT)
    { ERROR(("Link invalid for censored Gaussian family"));
      return(LF_LNK);
    }
    pz = pnorm(-z,0.0,1.0);
    dp = ((z>6) ? ptail(-z) : exp(-z*z/2)/pz)/2.5066283;
    res[ZLIK] = w*log(pz);
    res[ZDLL] = w*dp;
    res[ZDDLL]= w*dp*(dp-z);
    return(LF_OK);
  }
  res[ZLIK] = -w*z*z/2; 
  switch(link)
  { case LIDENT:
      res[ZDLL] = w*z;
      res[ZDDLL]= w;
      break;
    case LLOG:
      res[ZDLL] = w*z*mean;
      res[ZDDLL]= w*mean*mean;
      break;
    case LLOGIT:
      res[ZDLL] = w*z*mean*(1-mean);
      res[ZDDLL]= w*mean*mean*(1-mean)*(1-mean);
      break;
    default:
      ERROR(("Invalid link for Gaussian family"));
      return(LF_LNK);
  }
  return(LF_OK);
}

INT famrobu(y,mean,th,link,res,cens,w,rs)
double y, mean, th, *res, w, rs;
INT link, cens;
{ double z, sw;
  if (link==LINIT)
  { res[ZDLL] = w*y;
    return(LF_OK);
  }
  sw = (w==1.0) ? 1.0 : sqrt(w); /* don't want unnecess. sqrt! */
  z = sw*(y-mean)/rs;
  res[ZLIK] = (fabs(z)<HUBERC) ? -z*z/2 : HUBERC*(HUBERC/2.0-fabs(z));
  if (z< -HUBERC)
  { res[ZDLL] = -sw*HUBERC/rs;
    res[ZDDLL]= 0.0;
    return(LF_OK);
  }
  if (z> HUBERC)
  { res[ZDLL] = sw*HUBERC/rs;
    res[ZDDLL]= 0.0;
    return(LF_OK);
  }
  res[ZDLL] =  sw*z/rs;
  res[ZDDLL] = w/(rs*rs);
  return(LF_OK);
}

INT famcauc(y,p,th,link,res,cens,w,rs)
double y, p, th, *res, w, rs;
INT link, cens;
{ double z;
  if (link!=LIDENT)
  { ERROR(("Invalid link in famcauc"));
    return(LF_LNK);
  }
  z = w*(y-th)/rs;
  res[ZLIK] = -log(1+z*z);
  res[ZDLL] = 2*w*z/(rs*(1+z*z));
  res[ZDDLL] = 2*w*w*(1-z*z)/(rs*rs*(1+z*z)*(1+z*z));
  return(LF_OK);
}

INT famrbin(y,p,th,link,res,cens,w)
double y, p, th, *res, w;
INT link, cens;
{ double s2y;
  if (link==LINIT)
  { res[ZDLL] = y;
    return(LF_OK);
  }
  if ((y<0) | (y>w)) /* goon observation; delete it */
  { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
    return(LF_OK);
  }
  res[ZLIK] = (th<0) ? th*y-w*log(1+exp(th)) : th*(y-w)-w*log(1+exp(-th));
  if (y>0) res[ZLIK] -= y*log(y/w);
  if (y<w) res[ZLIK] -= (w-y)*log(1-y/w);
  res[ZDLL] = (y-w*p);
  res[ZDDLL]= w*p*(1-p);
  if (-res[ZLIK]>HUBERC*HUBERC/2.0)
  { s2y = sqrt(-2*res[ZLIK]);
    res[ZLIK] = HUBERC*(HUBERC/2.0-s2y);
    res[ZDLL] *= HUBERC/s2y;
    res[ZDDLL] = HUBERC/s2y*(res[ZDDLL]-1/(s2y*s2y)*w*p*(1-p));
  }
  return(LF_OK);
}

INT fambino(y,p,th,link,res,cens,w)
double y, p, th, *res, w;
INT link, cens;
{ double wp;
  if (link==LINIT)
  { if (y<0) y = 0;
    if (y>w) y = w;
    res[ZDLL] = y;
    return(LF_OK);
  }
  wp = w*p;
  if (link==LIDENT)
  { if ((p<=0) && (y>0)) return(LF_BADP);
    if ((p>=1) && (y<w)) return(LF_BADP);
    res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
    if (y>0)
    { res[ZLIK] += y*log(wp/y);
      res[ZDLL] += y/p;
      res[ZDDLL]+= y/(p*p);
    }
    if (y<w)
    { res[ZLIK] += (w-y)*log((w-wp)/(w-y));
      res[ZDLL] -= (w-y)/(1-p);
      res[ZDDLL]+= (w-y)/SQR(1-p);
    }
    return(LF_OK);
  }
  if (link==LLOGIT)
  { if ((y<0) | (y>w)) /* goon observation; delete it */
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
      return(LF_OK);
    }
    res[ZLIK] = (th<0) ? th*y-w*log(1+exp(th)) : th*(y-w)-w*log(1+exp(-th));
    if (y>0) res[ZLIK] -= y*log(y/w);
    if (y<w) res[ZLIK] -= (w-y)*log(1-y/w);
    res[ZDLL] = (y-wp);
    res[ZDDLL]= wp*(1-p);
    return(LF_OK);
  }
  if (link==LASIN)
  { if ((p<=0) && (y>0)) return(LF_BADP);
    if ((p>=1) && (y<w)) return(LF_BADP);
    if ((th<0) | (th>PI/2)) return(LF_BADP);
    res[ZDLL] = res[ZDDLL] = res[ZLIK] = 0;
    if (y>0)
    { res[ZDLL] += 2*y*sqrt((1-p)/p);
      res[ZLIK] += y*log(wp/y);
    }
    if (y<w)
    { res[ZDLL] -= 2*(w-y)*sqrt(p/(1-p));
      res[ZLIK] += (w-y)*log((w-wp)/(w-y));
    }
    res[ZDDLL] = 4*w;
    return(LF_OK);
  }
  ERROR(("link %d invalid for binomial family",link));
  return(LF_LNK);
}

INT fampois(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ double wmu, pt, dp, dq;
  if (link==LINIT)
  { res[ZDLL] = MAX(y,0.0);
    return(LF_OK);
  }
  wmu = w*mean;
  if (cens)
  { if (y<=0)
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
      return(LF_OK);
    }
    pt = igamma(wmu,y);
    dp = exp((y-1)*log(wmu)-wmu-LGAMMA(y))/pt;
    dq = dp*((y-1)/wmu-1);
    res[ZLIK] = log(pt);
    if (link==LLOG)
    { res[ZDLL] = dp*wmu;
      res[ZDDLL]= -(dq-dp*dp)*wmu*wmu-dp*wmu;
      return(LF_OK);
    }
    if (link==LIDENT)
    { res[ZDLL] = dp*w;
      res[ZDDLL]= -(dq-dp*dp)*w*w;
      return(LF_OK);
    }
    if (link==LSQRT)
    { res[ZDLL] = dp*2*w*th;
      res[ZDDLL]= -(dq-dp*dp)*(4*w*w*mean)-2*dp*w;
      return(LF_OK);
  } }
  if (link==LLOG)
  { if (y<0) /* goon observation - delete it */
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0;
      return(LF_OK);
    }
    res[ZLIK] = res[ZDLL] = y-wmu;
    if (y>0) res[ZLIK] += y*(th-log(y/w));
    res[ZDDLL] = wmu;
    return(LF_OK);
  }
  if (link==LIDENT)
  { if ((mean<=0) && (y>0)) return(LF_BADP);
    res[ZLIK] = y-wmu;
    res[ZDLL] = -w;
    res[ZDDLL] = 0;
    if (y>0)
    { res[ZLIK] += y*log(wmu/y);
      res[ZDLL] += y/mean;
      res[ZDDLL]= y/(mean*mean);
    }
    return(LF_OK);
  }
  if (link==LSQRT)
  { if ((mean<=0) && (y>0)) return(LF_BADP);
    res[ZLIK] = y-wmu;
    res[ZDLL] = -2*w*th;
    res[ZDDLL]= 2*w;
    if (y>0)
    { res[ZLIK] += y*log(wmu/y);
      res[ZDLL] += 2*y/th;
      res[ZDDLL]+= 2*y/mean;
    }
    return(LF_OK);
  }
  ERROR(("link %d invalid for Poisson family",link));
  return(LF_LNK);
}

INT famgamm(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ double pt, dg;
  if (link==LINIT)
  { res[ZDLL] = MAX(y,0.0);
    return(LF_OK);
  }
  if ((mean<=0) & (y>0)) return(LF_BADP);
  if (cens)
  { if (y<=0)
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
      return(LF_OK);
    }
    if (link==LLOG)
    { pt = 1-igamma(y/mean,w);
      dg = exp((w-1)*log(y/mean)-y/mean-LGAMMA(w));
      res[ZLIK] = log(pt);
      res[ZDLL] = y*dg/(mean*pt);
      res[ZDDLL]= dg*(w*y/mean-y*y/(mean*mean))/pt+SQR(res[ZDLL]);
      return(LF_OK);
    }
    if (link==LINVER)
    { pt = 1-igamma(th*y,w);
      dg = exp((w-1)*log(th*y)-th*y-LGAMMA(w));
      res[ZLIK] = log(pt);
      res[ZDLL] = -y*dg/pt;
      res[ZDDLL]= dg*y*((w-1)*mean-y)/pt+SQR(res[ZDLL]);
      return(LF_OK);
    }
  }
  else
  { if (y<0) WARN(("Negative Gamma observation"));
    if (link==LLOG)
    { res[ZLIK] = -y/mean+w*(1-th);
      if (y>0) res[ZLIK] += w*log(y/w);
      res[ZDLL] = y/mean-w;
      res[ZDDLL]= y/mean;
      return(LF_OK);
    }
    if (link==LINVER)
    { res[ZLIK] = -y/mean+w-w*log(mean);
      if (y>0) res[ZLIK] += w*log(y/w);
      res[ZDLL] = -y+w*mean;
      res[ZDDLL]= w*mean*mean;
      return(LF_OK);
    }
    if (link==LIDENT)
    { res[ZLIK] = -y/mean+w-w*log(mean);
      if (y>0) res[ZLIK] += w*log(y/w);
      res[ZDLL] = (y-mean)/(mean*mean);
      res[ZDDLL]= w/(mean*mean);
      return(LF_OK);
    }
  }
  ERROR(("link %d invalid for Gamma family",link));
  return(LF_LNK);
}

INT famgeom(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ double p, pt, dp, dq;
  if (link==LINIT)
  { res[ZDLL] = MAX(y,0.0);
    return(LF_OK);
  }
  p = 1/(1+mean);
  if (cens) /* censored observation */
  { if (y<=0)
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0;
      return(LF_OK);
    }
    pt = 1-ibeta(p,w,y);
    dp = -exp(LGAMMA(w+y)-LGAMMA(w)-LGAMMA(y)+(y-1)*th+(w+y-2)*log(p))/pt;
    dq = ((w-1)/p-(y-1)/(1-p))*dp;
    res[ZLIK] = log(pt);
    res[ZDLL] = -dp*p*(1-p);
    res[ZDDLL]= (dq-dp*dp)*p*p*(1-p)*(1-p)+dp*(1-2*p)*p*(1-p);
    res[ZDDLL]= -res[ZDDLL];
    return(LF_OK);
  }
  else
  { res[ZLIK] = (y+w)*log((y/w+1)/(mean+1));
    if (y>0) res[ZLIK] += y*log(w*mean/y);
    if (link==LLOG)
    { res[ZDLL] = (y-w*mean)*p;
      res[ZDDLL]= (y+w)*p*(1-p);
      return(LF_OK);
    }
    if (link==LIDENT)
    { res[ZDLL] = (y-w*mean)/(mean*(1+mean));
      res[ZDDLL]= w/(mean*(1+mean));
      return(LF_OK);
    }
  }
  ERROR(("link %d invalid for geometric family",link));
  return(LF_LNK);
}

INT famweib(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ double yy;
  yy = pow(y,w);
  if (link==LINIT)
  { res[ZDLL] = MAX(yy,0.0);
    return(LF_OK);
  }
  if (cens)
  { res[ZLIK] = -yy/mean;
    res[ZDLL] = res[ZDDLL] = yy/mean;
    return(LF_OK);
  }
  res[ZLIK] = 1-yy/mean-th;
  if (yy>0) res[ZLIK] += log(w*yy);
  res[ZDLL] = -1+yy/mean;
  res[ZDDLL]= yy/mean;
  return(LF_OK);
}

INT famcirc(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
INT link, cens;
{ if (link==LINIT)
  { res[ZDLL] = w*sin(y);
    res[ZLIK] = w*cos(y);
    return(LF_OK);
  }
  res[ZDLL] = w*sin(y-mean);
  res[ZDDLL]= w*cos(y-mean);
  res[ZLIK] = res[ZDDLL]-w;
  return(LF_OK);
}

void robustify(res,rs)
double *res, rs;
{ double sc, z;
  sc = rs*HUBERC;
  if (res[ZLIK] > -sc*sc/2) return;
  z = sqrt(-2*res[ZLIK]);
  res[ZDDLL]= -sc*res[ZDLL]*res[ZDLL]/(z*z*z)+sc*res[ZDDLL]/z;
  res[ZDLL]*= sc/z;
  res[ZLIK] = sc*sc/2-sc*z;
}

double lf_link(y,lin)
double y;
INT lin;
{ switch(lin)
  { case LIDENT: return(y);
    case LLOG:   return(log(y));
    case LLOGIT: return(logit(y));
    case LINVER: return(1/y);
    case LSQRT:  return(sqrt(fabs(y)));
    case LASIN:  return(asin(sqrt(y)));
  }
  ERROR(("link: unknown link %d",lin));
  return(0.0);
}

double invlink(th,lin)
double th;
INT lin;
{ switch(lin)
  { case LIDENT: return(th);
    case LLOG:   return(lf_exp(th));
    case LLOGIT: return(expit(th));
    case LINVER: return(1/th);
    case LSQRT:  return(th*fabs(th));
    case LASIN:  return(sin(th)*sin(th));
    case LINIT:  return(0.0);
  }
  ERROR(("invlink: unknown link %d",lin));
  return(0.0);
}

INT links(th,y,fam,lin,res,cd,w,rs) /* the link and various related functions */
double th, y, *res, w, cd, rs;
INT fam, lin;
{ double mean;
  INT c, link, st;
  c = (INT)cd; link = (INT)lin;

  mean = res[ZMEAN] = invlink(th,lin);
  if (lf_error) return(LF_LNK);

  switch(fam&63)
  { case THAZ:
    case TDEN:
    case TRAT: return(famdens(mean,th,link,res,c,w));
    case TGAUS: st = famgaus(y,mean,th,link,res,c,w);
                break;
    case TLOGT: st = fambino(y,mean,th,link,res,c,w);
                break;
    case TRBIN: return(famrbin(y,mean,th,link,res,c,w));
    case TPROB:
    case TPOIS: st = fampois(y,mean,th,link,res,c,w);
                break;
    case TGAMM: st = famgamm(y,mean,th,link,res,c,w);
                break;
    case TGEOM: st = famgeom(y,mean,th,link,res,c,w);
                break;
    case TWEIB: return(famweib(y,mean,th,link,res,c,w));
    case TCIRC: st = famcirc(y,mean,th,link,res,c,w);
                break;
    case TROBT: return(famrobu(y,mean,th,link,res,c,w,rs));
    case TCAUC: return(famcauc(y,mean,th,link,res,c,w,rs));
    default:
      ERROR(("links: invalid family %d",fam));
      return(LF_FAM);
  }
  if (st!=LF_OK) return(st);
  if (link==LINIT) return(st);
  if ((fam&128)==128) robustify(res,rs);
  return(st);
}

/*
  stdlinks is a version of links when family, link, response e.t.c
  all come from the standard places.
*/
INT stdlinks(res,lf,i,th,rs)
lfit *lf;
double th, rs, *res;
INT i;
{ return(links(th,resp(lf,i),lf->mi[MTG],lf->mi[MLINK],res,cens(lf,i),prwt(lf,i),rs));
}

/*
 *  functions used in variance, skewness, kurtosis calculations
 *  in scb corrections.
 */

double b2(th,tg,w)
double th, w;
INT tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(w);
    case TPOIS: return(w*lf_exp(th));
    case TLOGT:
      y = expit(th);
      return(w*y*(1-y));
  }
  ERROR(("b2: invalid family %d",tg));
  return(0.0);
}

double b3(th,tg,w)
double th, w;
INT tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(0.0);
    case TPOIS: return(w*lf_exp(th));
    case TLOGT:
      y = expit(th);
      return(w*y*(1-y)*(1-2*y));
  }
  ERROR(("b3: invalid family %d",tg));
  return(0.0);
}

double b4(th,tg,w)
double th, w;
INT tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(0.0);
    case TPOIS: return(w*lf_exp(th));
    case TLOGT:
      y = expit(th); y = y*(1-y);
      return(w*y*(1-6*y));
  }
  ERROR(("b4: invalid family %d",tg));
  return(0.0);
}
