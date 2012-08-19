/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 * 
 *
 *  ar_setfunction sets a function pointer on the arithmetic structure,
 *    according to the first len compnents of the string z.
 *    Also sets the rt->cmd field, and returns the required number of
 *    arguments.
 */

#include "local.h"

double vrexp()  { return(rexp(1.0)); }
double vrnorm() { return(rnorm(0.0,1.0)); }
double vpnorm(double x){ return(pnorm(x,0.0,1.0)); }
double vdnorm(double x){ return(exp(-x*x/2)/S2PI); }
double dummyf() { return(0.0); }
double frac(double x) { return(x-floor(x)); }

double vmin(v)
vari *v;
{ int i;
  double z, x;
  z = vitem(v,0);
  for (i=1; i<vlength(v); i++)
  { x = vitem(v,i);
    if (x<z) z = x;
  }
  return(z);
}

double vmax(v)
vari *v;
{ int i;
  double z, x;
  z = vitem(v,0);
  for (i=1; i<vlength(v); i++)
  { x = vitem(v,i);
    if (x>z) z = x;
  }
  return(z);
}

double vsum(v)
vari *v;
{ int i;
  double z;
  z = 0.0;
  for (i=0; i<vlength(v); i++) z += vitem(v,i);
  return(z);
}
double vmean(vari *v) { return(vsum(v)/vlength(v)); }

int ar_setfunction(rt,z,len)
arstruct *rt;
int z;
int len;
{ int rargs;
  rt->f = NULL;
  rt->cmd = 'f';
  rargs = 1;
  if (len==3)
  { if (stm(z,"sin",3)) rt->f = sin;
    if (stm(z,"cos",3)) rt->f = cos;
    if (stm(z,"tan",3)) rt->f = tan;
    if (stm(z,"exp",3)) rt->f = exp;
    if (stm(z,"log",3)) rt->f = log;
    if (stm(z,"abs",3)) rt->f = fabs;
    if (stm(z,"seq",3)) { rt->f = dummyf; rt->cmd = 'Q'; rargs=3; }
    if (stm(z,"min",3)) { rt->f = vmin; rt->cmd = 'S'; rargs=1; }
    if (stm(z,"max",3)) { rt->f = vmax; rt->cmd = 'S'; rargs=1; }
    if (stm(z,"sum",3)) { rt->f = vsum; rt->cmd = 'S'; rargs=1; }
    if (stm(z,"rep",3)) { rt->f = dummyf; rt->cmd = 'R'; rargs=2; }
  }
  if (len==4)
  { if (stm(z,"frac",4)) rt->f = frac;
    if (stm(z,"sqrt",4)) rt->f = sqrt;
    if (stm(z,"rexp",4)) { rt->f = vrexp; rt->cmd = 'G'; }
    if (stm(z,"mean",4)) { rt->f = vmean; rt->cmd = 'S'; rargs=1; }
  }
  if (len==5)
  { if (stm(z,"floor",5)) rt->f = floor;
    if (stm(z,"pnorm",5)) rt->f = vpnorm;
    if (stm(z,"dnorm",5)) rt->f = vdnorm;
    if (stm(z,"logit",5)) rt->f = logit;
    if (stm(z,"expit",5)) rt->f = expit;
    if (stm(z,"runif",5)) { rt->f = runif; rt->cmd='G'; }
    if (stm(z,"rnorm",5)) { rt->f = vrnorm;rt->cmd='G'; }
    if (stm(z,"rpois",5)) { rt->f = rpois; rt->cmd='H'; rargs=2; }
  }
  if (len==6)
  { if (stm(z,"sample",6)) { rt->f = dummyf;rt->cmd='M'; rargs=2; }
    if (stm(z,"fitted",6)) { rt->f = dummyf; rt->cmd='Z'; rargs=0; }
  }
  if (len==9)
  { if (stm(z,"residuals",9))
    { rt->f= dummyf; rt->cmd='Y'; rargs=0; }
  }
  if (rt->f==NULL)
  { rt->cmd = 'e';
    ERROR(("unknown function"));
  }
  return(rargs);
}
