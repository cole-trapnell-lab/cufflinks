/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

#define LOG_2           0.6931471805599453094172321214581765680755
#define IBETA_LARGE     1.0e30
#define IBETA_SMALL     1.0e-30
#define IGAMMA_LARGE    1.0e30
#define DOUBLE_EP     2.2204460492503131E-16

double dchisq(x, df)
double x, df;
{ return(exp(log(x/2)*(df/2-1) - x/2 - LGAMMA(df/2) - LOG_2));
}

double df(x, df1, df2)
double x, df1, df2;
{ double p;
  p = exp(LGAMMA((df1+df2)/2) + df1/2*log(df1/df2) + (df1/2-1)*log(x)
               - LGAMMA(df1/2) - LGAMMA(df2/2) - (df1+df2)/2*log(1+x*df1/df2));
  return(p);
}

double ibeta(x, a, b)
double x, a, b;
{ int flipped = 0, i, k, count;
  double I = 0, temp, pn[6], ak, bk, next, prev, factor, val;
  if (x <= 0) return(0);
  if (x >= 1) return(1);
/* use ibeta(x,a,b) = 1-ibeta(1-x,b,z) */
  if ((a+b+1)*x > (a+1))
  { flipped = 1;
    temp = a;
    a = b;
    b = temp;
    x = 1 - x;
  }
  pn[0] = 0.0;
  pn[2] = pn[3] = pn[1] = 1.0;
  count = 1;
  val = x/(1.0-x);
  bk = 1.0;
  next = 1.0;
  do
  { count++;
    k = count/2;
    prev = next;
    if (count%2 == 0)
      ak = -((a+k-1.0)*(b-k)*val)/((a+2.0*k-2.0)*(a+2.0*k-1.0));
    else
      ak = ((a+b+k-1.0)*k*val)/((a+2.0*k)*(a+2.0*k-1.0));
    pn[4] = bk*pn[2] + ak*pn[0];
    pn[5] = bk*pn[3] + ak*pn[1];
    next = pn[4] / pn[5];
    for (i=0; i<=3; i++)
      pn[i] = pn[i+2];
    if (fabs(pn[4]) >= IBETA_LARGE)
      for (i=0; i<=3; i++)
        pn[i] /= IBETA_LARGE;
    if (fabs(pn[4]) <= IBETA_SMALL)
      for (i=0; i<=3; i++)
        pn[i] /= IBETA_SMALL;
  } while (fabs(next-prev) > DOUBLE_EP*prev);
  factor = a*log(x) + (b-1)*log(1-x);
  factor -= LGAMMA(a+1) + LGAMMA(b) - LGAMMA(a+b);
  I = exp(factor) * next;
  return(flipped ? 1-I : I);
}

/*
 * Incomplete gamma function.
 * int_0^x u^{df-1} e^{-u} du / Gamma(df).
 */
double igamma(x, df)
double x, df;
{ double factor, term, gintegral, pn[6], rn, ak, bk;
  int i, count, k;
  if (x <= 0.0) return(0.0);

  if (df < 1.0)
    return( exp(df*log(x)-x-LGAMMA(df+1.0)) + igamma(x,df+1.0) );

/* 
 * this is unstable for large df
 */
  factor = exp(df*log(x) - x - LGAMMA(df));

  if (x > 1.0 && x >= df)
  {
    pn[0] = 0.0;
    pn[2] = pn[1] = 1.0;
    pn[3] = x;
    count = 1;
    rn = 1.0 / x;
    do
    { count++;
      k = count / 2;
      gintegral = rn;
      if (count%2 == 0)
      { bk = 1.0;
        ak = (double)k - df;
      } else
      { bk = x;
        ak = (double)k;
      }
      pn[4] = bk*pn[2] + ak*pn[0];
      pn[5] = bk*pn[3] + ak*pn[1];
      rn = pn[4] / pn[5];
      for (i=0; i<4; i++)
        pn[i] = pn[i+2];
      if (pn[4] > IGAMMA_LARGE)
        for (i=0; i<4; i++)
          pn[i] /= IGAMMA_LARGE;
    } while (fabs(gintegral-rn) > DOUBLE_EP*rn);
    gintegral = 1.0 - factor*rn;
  }
  else
  { /*   For x<df, use the series
     *   dpois(df,x)*( 1 + x/(df+1) + x^2/((df+1)(df+2)) + ... )
     *   This could be slow if df large and x/df is close to 1.
     */
    gintegral = term = 1.0;
    rn = df;
    do
    { rn += 1.0;
      term *= x/rn;
      gintegral += term;
    } while (term > DOUBLE_EP*gintegral);
    gintegral *= factor/df;
  }
  return(gintegral);
}

double pf(q, df1, df2)
double q, df1, df2;
{ return(ibeta(q*df1/(df2+q*df1), df1/2, df2/2));
}

double pchisq(q, df)
double q, df;
{ return(igamma(q/2, df/2));
}

#ifdef RVERSION
extern double Rf_pnorm5();
double pnorm(x,mu,s)
double x, mu, s;
{ return(Rf_pnorm5(x, mu, s, 1L, 0L));
}
#else
double pnorm(x,mu,s)
double x, mu, s;
{ if(x == mu)
    return(0.5);
  x = (x-mu)/s;
  if(x > 0) return((1 + erf(x/SQRT2))/2);
  return(erfc(-x/SQRT2)/2);
}
#endif
