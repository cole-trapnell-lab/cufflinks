/*
 *   Copyright (c) 1998-2000 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   The design structure used in Locfit, and associated macro definitions.
 */

typedef struct {
  vari *dw, *index;
  double *xev;       /* fit point length p               */
  double *X;         /* design matrix, length n*p        */
  double *w, *di, *res, *th, *wd, h, xb[MXDIM];
  double *V, *P, *f1, *ss, *oc, *cf, llk;
  jacobian xtwx;     /* to store X'WVX and decomposition */
  int cfn[1+MXDIM], ncoef;
  int *fix;          /* indicator vector, showing fixed variables. */
  INT *ind, n, p, pref, (*itype)();
  INT (*vfun)();     /* pointer to the vertex processing function. */
} design;

#define cfn(des,i) (des->cfn[i])
#define d_x(des) ((des)->X)
#define d_xi(des,i) (&(des)->X[i*((des)->p)])
#define d_xij(des,i,j) ((des)->X[i*((des)->p)+j])
#define is_fixed(des,i) ((des)->fix[i]==1)

extern int des_reqd(), des_reqi();
