/*
 *   Copyright (c) 1998-2000 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *   Structures, typedefs etc used in Locfit
 */

typedef char varname[64];

/*
 *  Define the vari type for locfit variables and related macros.
 *  For the C version, an enhanced vari type is needed;
 *  for other versions a simple structure suffices.
 */
#ifdef CVERSION

typedef struct {
  varname name;
  INT n, bytes, mode, stat;
  double *dpr; } vari;
#define checkvarlen(v,n,name,mode) (createvar(name,STSYSTEM,n,mode))
#define vmode(v) ((v)->mode)

#else

typedef struct {
  INT n;
  double *dpr;
} vari;
#define viptr(v,i) (&(v)->dpr[i])
#define checkvarlen(v,len,name,mode) \
   ((((v)!=NULL) && (vlength(v) >= (len))) ? (v) : createvar((name),0,(len),(mode)))
#endif

#define vlength(v) ((v)->n)

typedef struct {
  char *arg, *val;
  vari *result;
  INT used; } carg;

typedef struct {
  void (*AddColor)(), (*SetColor)(), (*ClearScreen)(), (*TextDim)(), (*DoText)();
  void (*DrawPoint)(), (*DrawLine)(), (*DrawPatch)(), (*wrapup)();
  INT (*makewin)(), ticklength, defth, deftw;
} device;

typedef struct {
  vari *wk;
  double *coef, *xbar, *f;
  jacobian xtwx; } paramcomp;

typedef struct {
  vari *tw, *L, *iw, *xxev;
  double *x[MXDIM], *y, *w, *base, *c;
  double *coef, *nlx, *t0, *lik, *h, *deg;
  double *sv, dp[LEND], kap[3];
  double sca[MXDIM], fl[2*MXDIM], xl[2*MXDIM];
  INT *ce, *s, *lo, *hi, sty[MXDIM];
  INT *mg, nvm, ncm, vc;
  INT nl, nv, nnl, nce, nk, nn, mi[LENM], ord, deriv[MXDEG+2], nd;
  paramcomp pc;
  varname yname, xname[MXDIM], wname, bname, cname; } lfit;

#define datum(lf,i,j) (lf)->x[i][j]
#define dvari(lf,i)   (lf)->x[i]
#define evpt(lf,i) (&(lf)->xxev->dpr[(i)*(lf)->mi[MDIM]])
#define evptx(lf,i,k) ((lf)->xxev->dpr[(i)*(lf)->mi[MDIM]+(k)])

typedef struct {
  vari *data[MXDIM], *fit, *se;
  INT d, wh, gr;
} pplot;

typedef struct {
  char cmd;
  double x, *v, (*f)();
  INT m, nx[3];
  vari *vv; } arstruct;

typedef struct {
  vari *x, *y, *z;
  char type;
  INT id, t, n, nx, ny, pch; } plxyz;

typedef struct {
  double theta, phi, xl[2], yl[2], zl[2], sl[10];
  INT id, ty, nsl;
  char main[50], xlab[50], ylab[50], zlab[50];
  vari *track, *xyzs; } plots;

#define PLNONE 0
#define PLDATA 1
#define PLFIT  2
#define PLTRK  4

struct lfcol {
  char name[10];
  INT n, r, g, b;
};
