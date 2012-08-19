/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *  setstrval() is a  function for converting string arguments to Locfit's
 *    numeric values.  A typical call will be setstrval(lf.mi,MKER,"gauss").
 *
 *  components that can be set in this manner are
 *    MKER  (weight function)
 *    MKT   (kernel type -- spherical or product)
 *    MTG   (local likelihood family)
 *    MLINK (link function)
 *    MIT   (integration type for density estimation)
 *    MEV   (evaluation structure)
 *    MACRI (adaptive criterion)
 *
 *  INT ppwhat(str) interprets the preplot what argument.
 *  INT restyp(str) interprets the residual type argument.
 *
 */

#include "local.h"

static char *famil[17] =
  { "density", "ate",   "hazard",    "gaussian", "binomial",
    "poisson", "gamma", "geometric", "circular", "obust", "huber",
    "weibull", "cauchy","probab",    "logistic", "nbinomial", "vonmises" };
static int   fvals[17] = 
  { TDEN,  TRAT,  THAZ,  TGAUS, TLOGT,
    TPOIS, TGAMM, TGEOM, TCIRC, TROBT, TROBT,
    TWEIB, TCAUC, TPROB, TLOGT, TGEOM, TCIRC };

INT lffamily(z)
char *z;
{ INT quasi, robu, f;
  quasi = robu = 0;
  while ((z[0]=='q') | (z[0]=='r'))
  { quasi |= (z[0]=='q');
    robu  |= (z[0]=='r');
    z++;
  }
  f = pmatch(z,famil,fvals,16,-1);
  if ((z[0]=='o') | (z[0]=='a')) robu = 0;
  if (f==-1)
  { WARN(("unknown family %s",z));
    f = TGAUS;
  }
  if (quasi) f += 64;
  if (robu)  f += 128;
  return(f);
}

void getlffam(z,x)
char **z;
INT *x;
{ *x = lffamily(z[0]);
}

static char *wfuns[13] = {
  "rectangular", "epanechnikov", "bisquare",    "tricube",
  "triweight",   "gaussian",     "triangular",  "ququ",
  "6cub",        "minimax",      "exponential", "maclean", "parametric" };
static int wvals[13] = { WRECT, WEPAN, WBISQ, WTCUB,
  WTRWT, WGAUS, WTRIA, WQUQU, W6CUB, WMINM, WEXPL, WMACL, WPARM };

static char *ktype[3] = { "spherical", "product", "center" };
static int   kvals[3] = { KSPH, KPROD, KCE };

static char *ltype[8] = { "default", "canonical", "identity", "log",
                          "logi",    "inverse",   "sqrt",     "arcsin" };
static int   lvals[8] = { LDEFAU, LCANON, LIDENT, LLOG,
                          LLOGIT, LINVER, LSQRT,  LASIN };

static char *etype[9] = { "tree",     "phull", "data", "grid", "kdtree",
                          "kdcenter", "cross", "xbar", "none" };
static int   evals[9] = { ETREE, EPHULL, EDATA, EGRID, EKDTR,
                          EKDCE, ECROS,  EXBAR, ENONE };

static char *itype[6] = { "default", "multi", "product", "mlinear",
                          "hazard",  "monte" };
static int   ivals[6] = { IDEFA, IMULT, IPROD, IMLIN, IHAZD, IMONT };

static char *atype[5] = { "none", "cp", "ici", "mindex", "ok" };
static int   avals[5] = { ANONE, ACP, AKAT, AMDI, AOK };

void setstrval(mi,v,z)
INT *mi, v;
char *z;
{ 
  switch(v)
  { case MKER:
      mi[v] = pmatch(z, wfuns, wvals, 13, WTCUB);
      return;

    case MKT:
      mi[v] = pmatch(z, ktype, kvals, 3, KSPH);
      return;

    case MTG:
      mi[v] = lffamily(z);
      return;

    case MLINK:
      mi[v] = pmatch(z, ltype, lvals, 8, LDEFAU);
      return;

    case MIT:
      mi[v] = pmatch(z, itype, ivals, 6, IDEFA);
      return;

    case MEV:
      mi[v] = pmatch(z, etype, evals, 9, ETREE);
      return;

    case MACRI:
      mi[v] = pmatch(z, atype, avals, 5, ANONE);
      return;
  }

  WARN(("setstrval: invalid value %d",v));
  return;
}

static char *rtype[8] = { "deviance", "d2",    "pearson", "raw",
                          "ldot",     "lddot", "fit",     "mean" };
static int   rvals[8] = { RDEV, RDEV2, RPEAR, RRAW, RLDOT, RLDDT, RFIT, RMEAN};

static char *whtyp[8] = { "coef", "nlx", "infl", "band",
                          "degr", "like", "rdf", "vari" };
static int   whval[8] = { PCOEF, PNLX, PT0, PBAND, PDEGR, PLIK, PRDF, PVARI };

INT restyp(z)
char *z;
{ int val;
  
  val = pmatch(z, rtype, rvals, 8, -1);
  if (val==-1) ERROR(("Unknown type = %s",z));
  return((INT)val);
}

INT ppwhat(z)
char *z;
{ int val;
  
  val = pmatch(z, whtyp, whval, 8, -1);
  if (val==-1) ERROR(("Unknown what = %s",z));
  return((INT)val);
}
