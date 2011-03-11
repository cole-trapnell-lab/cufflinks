/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 *
 *   Functions for interpreting and manipulating command line
 *   arguments.
 */

#include "local.h"

char *argval(vari *v,int i)
{ if (i<0) return(NULL);
  return(((carg *)viptr(v,i))->val);
}

int getarg(v,s,un) /* un=1: unnamed permitted un=2: next unused */
vari *v;
int un;
char *s;
{ int i;
  if (un==2)
  { for (i=1; i<vlength(v); i++)
    { if (!argused(v,i))
      { setused(v,i);
        return(i);
      }
    }
    return(0);
  }
  for (i=1; i<vlength(v); i++)
  { if ((!argused(v,i)) && (argarg(v,i)!=NULL))
    { if (strcmp(argarg(v,i),s)==0)
      { setused(v,i);
        return(i);
      }
    }
  }
  if (!un) return(0);
  for (i=1; i<vlength(v); i++)
  { if ((!argused(v,i)) && (argarg(v,i)==NULL))
    { setused(v,i);
      return(i);
    }
  }
  return(0);
}

char *getargval(vari *v, char *s, int un)
{ int i;
  i = getarg(v,s,un);
  if (i==0) return(NULL);
  return(argval(v,i));
}

int readilist(ivec,z,n0,n1,pad)
char *z;
int *ivec, n0, n1, pad;
{ int i, n, nd;
  n = 1;
  for (i=0; i<strlen(z); i++)
    if (z[i]==',') { z[i]=' '; n++; }
  if (n>n1)
  { WARN(("too many items in ilist"));
    n = n1;
  }
  for (i=0; i<n; i++)
  { nd = sscanf(z,"%d",&ivec[i]);
    //if (nd!=1) WARN(("problem scaning ilist %s",&ivec[i]));
    if (i<n-1) while (*z!=' ') z++;
  }
  if (pad)
  { for (i=n; i<n1; i++) ivec[i] = ivec[0]; 
    n = n1;
  }
  if (n<n0) WARN(("too few items in ilist"));
  return(n);
}

int getlogic(v,i)
vari *v;
int i;
{ char *z;
  if (argarg(v,i)==NULL) return(1);
  z = argval(v,i);
  if ((z[0]=='T') | (z[0]=='t') | (z[0]=='1')) return(1);
  if ((z[0]=='F') | (z[0]=='f') | (z[0]=='0')) return(0);
  ERROR(("getlogic: invalid logical argument %s",z));
  return(0);
}
