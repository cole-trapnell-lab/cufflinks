typedef struct {
  int n;
  double *dpr;
} vari;    

typedef struct {
  double *Z, *Q, *dg, *f2;
  int p, sm; } xtwxstruc;

typedef struct {
  vari *wk;
  double *coef, *xbar, *f;
  xtwxstruc xtwx; } paramcomp;

typedef struct {
  vari *dw, *index;
  double *xev, *X, *w, *di, *res, *th, *wd, h, xb[15];
  double *V, *P, *f1, *ss, *oc, *cf, llk;
  xtwxstruc xtwx;
  int *ind, n, p, pref, (*itype)();
  int (*vfun)(); } design;

typedef struct {
  vari *tw, *L, *iw, *xxev;
  double *x[15], *y, *w, *base, *c, *xl;
  double *coef, *nlx, *t0, *lik, *h, *deg;
  double *sv, *fl, *sca, *dp, kap[3];
  int *ce, *s, *lo, *hi, sty[15];
  int *mg, nvm, ncm, vc;
  int nl, nv, nnl, nce, nk, nn, *mi, ord, deriv[9], nd;
  paramcomp pc;
  varname yname, xname[15], wname, bname, cname; } lfit;

extern void mlbcall(
             double *x, double *y,
             double *xx, double *ff, int n);
