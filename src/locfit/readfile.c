/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *   Function to read and create Locfit variables from an ASCII file.
 *
 *   Command syntax:
 *     locfit> readfile filename v1 v2 v3
 *   filename is the file to be read (no default extension or path is added).
 *   v1 v2 v3  etc are the names of the variables to create.
 *
 *   File format: The file should be a plain ASCII file organized
 *     in matrix format, with one variable per column and one observation
 *     per row. Fields are separated by spaces.
 *     
 */

#include "local.h"

extern char filename[];
static FILE *aaa;

void readfile(vc)
vari *vc;
{ int i, j, k, n, nv;
  char wc[50], *fn;
  double *dpr;
  vari *v;

  i = getarg(vc,"file",1);
  if (i==0)
  { ERROR(("readfile: no file"));
    return;
  }

  fn = argval(vc,i);
  setfilename(fn,"","r",0);
  if (lf_error) return;

  i = getarg(vc,"arith",0); /* now automatic - leave for backward compat. */

  aaa = fopen(filename,"r");
  v = createvar("readfile",STREADFI,0,VDOUBLE);

  n = 0;
  do
  { k = fscanf(aaa,"%s",wc);
    if (k==1)
    { vassn(v,n,darith(wc));
      n++;
    }
  } while (k==1);
  fclose(aaa);
  dpr = vdptr(v);
  deletevar(v);

  nv = 0;
  for (i=1; i<vlength(vc); i++)
    if (!argused(vc,i)) nv++;
  if (nv==0) { ERROR(("readfile: no variables specified")); return; }
  if (n%nv != 0)
    WARN(("number of items not multiple of number of variables"));

  n /= nv;
  transpose(dpr,n,nv);
  nv = 0;
  for (i=1; i<vlength(vc); i++)
    if (!argused(vc,i))
    { v = createvar(argval(vc,i),STREGULAR,n,VDOUBLE);
      if (v==NULL) return;
      for (j=0; j<n; j++) vassn(v,j,dpr[nv*n+j]);
      nv++;
      if (lf_error) return;
    }
  if (argarg(vc,0)!=NULL)
    dosavedata(vc,0);
  else
    for (i=1; i<vlength(vc); i++) setused(vc,i);
}
