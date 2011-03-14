/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 * 
 *
 *  Most of the changes formerly needed here are handled through
 *  the Makefiles and #ifdef's.
 *
 *
 */

#define CVERSION 1

#ifndef I_LF_H
#define I_LF_H

/*RVERSION*/

/*
 *   DIRSEP: '/' for unix; '\\' for DOS
 */
#ifdef DOS
#define DIRSEP '\\'
#else
#define DIRSEP '/'
#endif


/*
   Some older math libraries have no lgamma() function, and gamma(arg)
   actually returns log(gamma(arg)). If so, you need to change
   LGAMMA macro below.

   If all else fails, you can also use lflgamma().

   Use the definitions for erf, erfc and daws only if your
   math libraries don't include these functions.
 */
#ifdef DOS
#define LGAMMA(arg) lflgamma(arg)
#define erf(x) lferf(x)
#define erfc(x) lferfc(x)
#else
#define LGAMMA(arg) lgamma(arg)
#endif
#define daws(x) lfdaws(x)


/*
   Does your system support popen() and pclose()
   for pipes? For most flavours of unix, yes.
   For other OS's, usually not, and you'll need to
   uncomment the following line.

   (This is only relevant if making the C version).
*/
/* #define NOPIPES */


/*
   (the #ifdef's below should now take care of this).
   the INT type is used for all integers provided in .C() calls from S.
   For the S version, this should be long int.
   For the R version, should be int.
   For the C version, either is adequate.
   Usually this only makes a difference on 64 bit systems.
*/

//#ifndef SWINVERSION

//#ifdef RVERSION
//typedef int INT;
//#else
//#ifdef SVERSION
//typedef long int INT;
//#else
typedef int INT;
//#endif /* SVERSION */
//#endif /* RVERSION */

//#endif /* SWINVERSION */

/******** NOTHING BELOW HERE NEEDS CHANGING **********/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef RVERSION
#undef LGAMMA
#define LGAMMA(arg) Rf_lgammafn(arg)
extern double Rf_lgammafn();
#define SVERSION
#endif

#ifdef SWINVERSION
#define SVERSION
#include "newredef.h"
#endif

#include "mutil.h"
#include "lfcons.h"
#include "lfstruc.h"
#include "design.h"
#include "lffuns.h"

#ifdef CVERSION
//#undef printf
//#define printf lfprintf
//extern int lfprintf(const char *format, ...);
//extern int printf(const char *format, ...);
#endif

#ifdef SVERSION
#define printf printf
#endif

#ifdef INTERFACE
#define printf printf
#endif

#define ERROR(args) printf("Error: "), printf args , printf("\n"), lf_error=1
#define WARN(args)  printf("Warning: "),printf args, printf("\n")

#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define SGN(x) (((x)>0) ? 1 : -1)
#define SQR(x) ((x)*(x))
#define NOSLN 0.1278433
#define GFACT 2.5
#define EFACT 3.0

#define MAXCOLOR 20
#define MAXWIN 5

#ifdef SWINVERSION
#define ISWAP(a,b) { int zz; zz = a; a = b; b = zz; }
#else
#define ISWAP(a,b) { INT zz; zz = a; a = b; b = zz; }
extern INT lf_error;
#endif

extern INT lf_error;

double lf_exp(double x);

#endif /* I_LF_H */
