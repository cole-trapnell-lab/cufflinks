/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 *
 *  Functions for handling locfit variables in the C version.
 */

#include "local.h"

#define MAXV 1000
#define LF_WORK 1024

static char *db;
static INT lfwptr, lf_work;
vari root;

void initdb() /* initialize locfit's work space */
{ char *z;
  z = getenv("LFWORK");
  if (z==NULL) lf_work = LF_WORK;
    else sscanf(z,"%d",&lf_work);
  lf_work <<= 10;
  db = (char *)malloc(lf_work);
  root.stat = STSYSTEM;
  root.mode = VVARI;
  root.dpr = (double *)db;
  lfwptr = root.bytes = MAXV*sizeof(vari);
  root.n = 0;
}

vari *growvar(vold,n)
vari *vold;
INT n;
{ vari *vnew;
  INT reqd_bytes;

  if (vold==NULL)
  { ERROR(("growvar: NULL old"));
    return(NULL);
  }

  reqd_bytes = vbytes(n,vmode(vold));
  if (reqd_bytes <= vold->bytes) return(vold);

  vnew = createvar("_grow",vold->stat,n,vmode(vold));
  memcpy(vdptr(vnew),vdptr(vold),vbytes(vlength(vold),vmode(vold)));
  setvarname(vnew,vold->name);
  vlength(vnew) = vlength(vold);
  deletevar(vold);
  return(vnew);
}

INT vbytes(n,mode)
INT n, mode;
{ switch(mode)
  { case VDOUBLE: return(n*sizeof(double));
    case VINT:    return(n*sizeof(INT));
    case VCHAR:   return(n);
    case VARGL:   return(n*sizeof(carg));
    case VPREP:   return(sizeof(pplot));
    case VARC:    return(n*sizeof(arstruct));
    case VVARI:   return(n*sizeof(vari));
    case VXYZ:    return(n*sizeof(plxyz));
  }
  ERROR(("unknown mode %d in vbytes",mode));
  return(0);
}

void *viptr(v,i) /* return pointer to ith data item, take account of mode */
vari *v;
INT i;
{ switch(vmode(v))
  { case VDOUBLE: return(&v->dpr[i]);
    case VCHAR: return(&((char *)v->dpr)[i]);
    case VARGL: return(&((carg *)v->dpr)[i]);
    case VARC:  return(&((arstruct *)v->dpr)[i]);
    case VVARI: return(&((vari *)v->dpr)[i]);
    case VXYZ:  return(&((plxyz *)v->dpr)[i]);
  }
  ERROR(("Unknown mode %d in viptr",vmode(v)));
  return(NULL);
}

/* vdptr with NULL check */
double *vdptr(v)
vari *v;
{ if (v==NULL) return(NULL);
  return(v->dpr);
}

/* return the i'th data item. Cyclic. */
double vitem(v,i)
vari *v;
INT i;
{ int index;
  if ((v==NULL) || (vlength(v)==0)) return(0.0);
  index = i % vlength(v);
  switch(v->mode)
  { case VDOUBLE: return( vdptr(v)[index] );
    case VINT:
    { INT *z;
      z = (INT *)vdptr(v);
      return(z[index]);
    }
    case VCHAR:
    { char *z;
      z = (char *)vdptr(v);
      return(z[index]);
    }
  }
  ERROR(("Invalid mode in vitem()"));
  return(0.0);
}

void vassn(v,i,x)
vari *v;
INT i;
double x;
{ vdptr(v)[i] = x;
}

void setvarname(v,name)
vari *v;
varname name;
{ if (strcmp(v->name,name)==0) return;
  deletename(name);
  strcpy(v->name,name);
}

/*
  findvar finds the variable name.
  err=0, keep quiet if not found; 1 produce error message.
  *n returns length of variable (if initially>0, checks length)
*/

vari *findvar(name,err,n)
varname name;
INT err, *n;
{ INT i, status;
  vari *v;

  if (strcmp(name,"_NuLl")==0) return(NULL);

  for (i=0; i<root.n; i++)
  { v = viptr(&root,i);
    status = v->stat;
    if ((strcmp(v->name,name)==0) &&
      ((status!=STHIDDEN)&(status!=STEMPTY)) )
    { if (n==NULL) return(v);
      if (*n==-1) *n = vlength(v);
      if ((*n==0) | (*n==vlength(v))) return(v);
      if (err) ERROR(("Variable %s has wrong length",name));
      return(NULL);
    }
  }
  if (err) ERROR(("Variable %s not found",name));
  return(NULL);
}

void deletevar(v) /* delete variable, or top variable if NULL */
vari *v;
{ if (root.n==0) return;
  if (v!=NULL) v->stat = STEMPTY;
  if ((v==NULL) || (v==viptr(&root,root.n-1))) /* top variable */
  { root.n--;
    lfwptr -= ((vari *)viptr(&root,root.n))->bytes;
  }
}

void deleteifhidden(v)
vari *v;
{ if (v==NULL) return;
  if (v->stat == STHIDDEN) deletevar(v);
}

void deletename(name) /* delete variable name, or top variable if NULL */
varname name;
{ vari *v;
  v = findvar(name,0,NULL);
  if (v!=NULL) deletevar(v);
}

vari *createvar(name,status,n,mode)
varname name;
INT status, n, mode;
{
  INT i, bytes;
  vari *v;

  /*
     compute the length of the variable in bytes. some systems
     mess up is this is not a multiple of 8.
  */
  bytes = vbytes(n,mode);
  while ( bytes&8 > 0 ) bytes++;

  if (lf_error) return(NULL);
  if ((status==STSYSTEM)|(status==STREGULAR)|(status==STPLOTVAR))
    deletename(name);

  if (status!=STREADFI)
  { for (i=0; i<root.n; i++) /* does unused variable have space allocated? */
    { v = viptr(&root,i);
      if ((v->stat == STEMPTY) && (v->bytes >= bytes))
      { strcpy(v->name,name);
        vlength(v) = n;
        v->stat = status;
        v->mode = mode;
        return(v);
      }
    }
  }

  /* must allocate next variable. First, is there space? */
  if (root.n==MAXV) ERROR(("Too many variables"));
  if ((status!=STSYSPEC) && (lfwptr+bytes > lf_work))
    ERROR(("Insufficient space for variable creation"));
  if (lf_error) return(NULL);

  v = viptr(&root,root.n);
  strcpy(v->name,name);
  vlength(v) = n;
  v->stat = status;
  v->bytes = bytes;
  v->mode = mode;
  if (status!=STSYSPEC)
  { v->dpr = (double *)(&db[lfwptr]);
    lfwptr += bytes;
  }
  root.n++;
  return(v);
}
