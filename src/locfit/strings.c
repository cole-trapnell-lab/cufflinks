/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *  miscellaneous string handling functions.
 *  used mostly in arith (C version) and lfstr.c
 *
 *  stm(u,v,k)        do the first k components of u, v match?
 *  ct_match(z1,z2)   counts number of matching components.
 *  pmatch(z,strings,vals,n,def)
 *                    finds the best match for z among the strings;
 *                    returns corresponding component of vals.
 *                    n = no of strings; def = value if no match.
 *  matchrt(z,i,i2,op,cl)
 *  matchlf(z,i1,i,op,cl)
 *                    Parenthesis matching. If op='(' and cl=')', matchrt
 *                    searches z, starting at z[i]='(' and ending and z[i2],
 *                    for the closing ')', taking care of nesting.
 *                    matchlf does the reverse.
 *
 *  checkltor(z,i1,i2,c)
 *                    Checks the string z, left to right from z[i1] to z[i2]
 *                    but skipping parenthesized () and [] parts, for the
 *                    occurence of any character from c. If a match is found,
 *                    the index is returned. If no match is found, return -1.
 *
 *  checkrtol(z,i1,i2,c) Same as checkltor(), but searches right to left.
 *  strip(z)          replaces underscores in z by spaces.
 */

#include "local.h"

/* do the first k components of u, v match? */
int stm(char *u, char *v, int k) { return((strncmp(u,v,k)==0)); }

int ct_match(z1, z2)
char *z1, *z2;
{ int ct = 0;
  while (z1[ct]==z2[ct])
  { if (z1[ct]=='\0') return(ct+1);
    ct++;
  }
  return(ct);
}

int pmatch(z, strings, vals, n, def)
char *z, **strings;
int *vals, n, def;
{ int i, ct, best, best_ct;
  best = -1;
  best_ct = 0;

  for (i=0; i<n; i++)
  { ct = ct_match(z,strings[i]);
    if (ct==strlen(z)+1) return(vals[i]);
    if (ct>best_ct) { best = i; best_ct = ct; }
  }
  if (best==-1) return(def);
  return(vals[best]);
}

int matchrt(z,i,i2,op,cl)
char *z, op, cl;
int i, i2;
{ int k;
  if (z[i] != op)
  { ERROR(("matchrt: wrong start character"));
    return(i);
  }
  k = 0;
  while (1)
  { if (z[i]==op) k++;
    if (z[i]==cl) k--;
    if (k==0) return(i);
    i++;
    if (i>i2)
    { ERROR(("matchrt: unbalanced %c%c: %s",op,cl,z));
      return(i);
    }
  }
}

int matchlf(z,i1,i,op,cl)
char *z, op, cl;
int i, i1;
{ int k;
  if (z[i] != cl)
  { ERROR(("matchlf: wrong end character"));
    return(i);
  }
  k = 0;
  while (1)
  { if (z[i]==op) k--;
    if (z[i]==cl) k++;
    if (k==0) return(i);
    i--;
    if (i<i1)
    { ERROR(("matchlf: unbalanced %c%c: %s",op,cl,z));
      return(i);
    }
  }
}

int checkltor(z,i1,i2,c)
char *z, *c;
int i1, i2;
{ int i;
  i = i1;
  while (i<=i2)
  {
    if (strchr(c,z[i]) != NULL) return(i);
    if (z[i]=='(') i = matchrt(z,i,i2,'(',')');
    if (z[i]=='[') i = matchrt(z,i,i2,'[',']');
    i++;
    if (lf_error) return(-1);
  }
  return(-1);
}

int checkrtol(z,i1,i2,c)
char *z, *c;
int i1, i2;
{ int i;
  i = i2;
  while (i >= i1)
  { if (strchr(c,z[i]) != NULL) return(i);
    if (z[i]==')') i = matchlf(z,i1,i,'(',')');
    if (z[i]==']') i = matchlf(z,i1,i,'[',']');
    i--;
    if (lf_error) return(-1);
  }
  return(-1);
}

void strip(z)
char *z;
{ do { if (*z=='_') *z=' '; } while (*(++z)!='\0');
}
