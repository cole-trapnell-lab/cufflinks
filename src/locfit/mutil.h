/*
 *   Copyright (c) 1998-2000 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   Headers for math utility functions.
 */

#ifndef I_MUT_H
#define I_MUT_H

typedef struct {
  double *Z;   /* jacobian matrix, length p*p          */
  double *Q;   /* eigenvalue matrix, length p*p        */
  double *wk;  /* work vector in eig_solve, length p   */
  double *dg;  /* diag vector in eigd, length p        */
  int p;       /* dimension */
  int st;      /* status    */
  int sm;      /* requested decomposition */
} jacobian;

/* m_jacob.c */
extern int jac_reqd();
extern double *jac_alloc();
extern void jacob_dec(),   chol_dec(),   eig_dec();
extern int  jacob_solve(), chol_solve(), eig_solve();
extern int  jacob_hsolve(),chol_hsolve(),eig_hsolve();
extern double jacob_qf(),  chol_qf(),    eig_qf();

/* m_max.c */
extern double max_grid(), max_golden(), max_quad(), max_nr();

/* solve.c */
extern double solve_secant(), solve_nr(), solve_fp();

#define BDF_NONE  0
#define BDF_EXPLEFT  1
#define BDF_EXPRIGHT 2

/* return codes for functions optimized by max_nr */
#define NR_OK 0
#define NR_INVALID 1
#define NR_BREAK   2
#define NR_REDUCE  3
#define NR_NCON  10
#define NR_NDIV  11


/* jacobian status definitions */
#define JAC_RAW 0
#define JAC_CHOL 1
#define JAC_EIG  2
#define JAC_EIGD 3


#endif  /* define I_MUT_H */
