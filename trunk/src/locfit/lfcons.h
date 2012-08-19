/*
 *   Copyright (c) 1998 Lucent Technologies.
 *   See README file for details.
 */

/*
  Numeric values for constants used in locfit
*/

/*
  MXDIM and MXDEG are maximum dimension and local polynomial
  degree for Locfit. Note that some parts of the code may be
  more restrictive.
*/
#define MXDIM 15
#define MXDEG 7

/*
  floating point constants
*/
#ifndef PI
#define PI    3.141592653589793238462643
#endif
#define S2PI  2.506628274631000502415765
#define SQRT2 1.4142135623730950488
#define SQRPI 1.77245385090552
#define LOGPI 1.144729885849400174143427
#define GOLDEN 0.61803398874989484820
#define HL2PI 0.91893853320467267 /* log(2pi)/2 */
#define SQRPI 1.77245385090552    /* sqrt(pi)   */

/*
  Criteria for adaptive local fitting  mi[MACRI]
  1: localized CP;  2: ICI (katkovnik);  3: curvature model index
  4: Increase bandwidth until locfit returns LF_OK
*/
#define ANONE 0
#define ACP  1
#define AKAT 2
#define AMDI 3
#define AOK  4

/*
  vector of double precision parameters.
  0, 1, 2 are the three components of the smoothing parameter.
  3 cut parameter for adaptive evaluation structures.
  4-8 are likelihood, degrees of freedom and residual variance,
  computed as part of the fit.
  Stored as the lf.dp vector.
*/
#define DALP 0
#define DFXH 1
#define DADP 2
#define DCUT 3
#define DLK  4
#define DT0  5
#define DT1  6
#define DRV  7
#define DSWT 8
#define DRSC 9
#define LEND 10

/*
  Evaluation structures mi[MEV]
  EFITP special for `interpolation' at fit points
*/
#define ENULL  0
#define ETREE  1
#define EPHULL 2
#define EDATA  3
#define EGRID  4
#define EKDTR  5
#define EKDCE  6
#define ECROS  7
#define EPRES  8
#define EXBAR  9
#define ENONE  10
#define EFITP  50

/*
  integer parameters: sample size; dimension; number of local parameters etc.
  stored as the lf.mi vector.
*/
#define MN     0
#define MP     1
#define MDEG0  2
#define MDEG   3
#define MDIM   4
#define MACRI  5
#define MKER   6
#define MKT    7
#define MIT    8
#define MMINT  9
#define MMXIT 10
#define MREN  11
#define MEV   12
#define MTG   13
#define MLINK 14
#define MDC   15
#define MK    16
#define MDEB  17
#define MGETH 18
#define MPC   19
#define MUBAS 20
#define LENM  21

/*
  Link functions mi[MLINK].
  Mostly as in table 4.1 of the book.
  LDEFAU and LCANON are used to select default and canonical
  links respectively. LINIT shouldn't be selected by user...
*/
#define LINIT  0
#define LDEFAU 1
#define LCANON 2
#define LIDENT 3
#define LLOG   4
#define LLOGIT 5
#define LINVER 6
#define LSQRT  7
#define LASIN  8

/*
  components of vector returned by the links() function
  in family.c. ZLIK the likelihood; ZMEAN = estimated mean;
  ZDLL = derivative of log-likelihood; ZDDLL = - second derivative
*/
#define LLEN  4
#define ZLIK  0
#define ZMEAN 1
#define ZDLL  2
#define ZDDLL 3

/*
  weight functions mi[MKER].
  see Table 3.1 or the function W() in weights.c for definitions.
*/
#define WRECT 1
#define WEPAN 2
#define WBISQ 3
#define WTCUB 4
#define WTRWT 5
#define WGAUS 6
#define WTRIA 7
#define WQUQU 8
#define W6CUB 9
#define WMINM 10
#define WEXPL 11
#define WMACL 12
#define WPARM 13

/*
  type of multivariate weight function mi[MKT]
  KSPH (spherical)  KPROD (product)
  others shouldn't be used at present.
*/
#define KSPH   1
#define KPROD  2
#define KCE    3
#define KLM    4

#define STANGL 4
#define STLEFT 5
#define STRIGH 6
#define STCPAR 7

/*
  Local likelihood family mi[MTG]
  for quasi-likelihood, add 64.
*/
#define TNUL 0
#define TDEN 1
#define TRAT 2
#define THAZ 3
#define TGAUS 4
#define TLOGT 5
#define TPOIS 6
#define TGAMM 7
#define TGEOM 8
#define TCIRC 9
#define TROBT 10
#define TRBIN 11
#define TWEIB 12
#define TCAUC 13
#define TPROB 14

/*
  Integration type mi[MIT] for integration in
  density estimation.
*/
#define INVLD 0
#define IDEFA 1
#define IMULT 2
#define IPROD 3
#define IMLIN 4
#define IHAZD 5
#define IMONT 7

/*
  For prediction functions, what to predict?
  PCOEF -- coefficients        PT0   -- influence function
  PNLX  -- ||l(x)||            PBAND -- bandwidth h(x)
  PDEGR -- local poly. degree  PLIK  -- max. local likelihood
  PRDF  -- local res. d.f.     PVARI -- ||l(x)||^2
*/
#define PCOEF 1
#define PT0   2
#define PNLX  3
#define PBAND 4
#define PDEGR 5
#define PLIK  6
#define PRDF  7
#define PVARI 8

/*
  Residual Types
*/
#define RDEV  1
#define RPEAR 2
#define RRAW  3
#define RLDOT 4
#define RDEV2 5
#define RLDDT 6
#define RFIT  7
#define RMEAN 8

/*
  components of the colour vector
*/
#define CBAK 0
#define CAXI 1
#define CTEX 2
#define CLIN 3
#define CPOI 4
#define CCON 5
#define CCLA 6
#define CSEG 7
#define CPA1 8
#define CPA2 9

/*
  variable types: double, INT, char, argument list
*/
#define VDOUBLE 0
#define VINT    1
#define VCHAR   2
#define VARGL   3
#define VPREP   4
#define VARC    5
#define VVARI   6
#define VXYZ    7

/*
  variable status
*/
#define STEMPTY   0
#define STREGULAR 1
#define STHIDDEN  3
#define STPLOTVAR 4
#define STSYSTEM  5
#define STSYSPEC  6
#define STREADFI  7

/*
  return status for the locfit() function
*/
#define LF_OK   0
#define LF_OOB  2   /* out of bounds, or large unstable parameter */
#define LF_PF   3   /* perfect fit; interpolation; deviance=0 */
#define LF_NCON 4   /* not converged */
#define LF_NOPT 6   /* no or insufficient points with non-zero wt */
#define LF_INFA 7   /* initial failure e.g. log(0) */
#define LF_DEMP 10  /* density -- empty integration region */
#define LF_XOOR 11  /* density -- fit point outside xlim region */
#define LF_DNOP 12  /* density version of 6 */
#define LF_FPROB 80
#define LF_BADP 81  /* bad parameters e.g. neg prob for binomial */
#define LF_LNK  82  /* invalid link */
#define LF_FAM  83  /* invalid family */
#define LF_ERR  99  /* error */
