/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */


/* Functions for reading/writing to LFData directory */

#include <unistd.h>
#include "local.h"

#ifdef CVERSION

FILE *lfd=NULL;
extern char *lfhome;
char filename[100];

void closefile()
{ fclose(lfd);
  lfd = NULL;
}

void openfile(mode)
char *mode;
{ if (lfd!=NULL) closefile();
  lfd = fopen(filename,mode);
}

/*
 *  setfilename() places name[] in the filename[] array.
 *    -- quotes are stripped from name.
 *    -- ext is the default extension; .ext is added as an extension.
 *       (unless fp=1).
 *       If ext = "lfd" or "fit", the LFData (and paths) are used.
 *    -- mode is a standard unix mode ("r", "w" etc)
 *    -- fp indicates the full path is given (used by Windoze GUI).
 *    -- checks for validity of filename and mode.
 *    -- returns 0 for successful, 1 for unsuccessful.
 */
INT setfilename(name,ext,mode,fp)
char *name, *ext, *mode;
INT fp;
{ char subp[20];
  int n, quote, use_lfd;

  n = strlen(name);
  quote = ((name[0]=='"') && (name[n-1]=='"'));
  if (quote)
  { name++;
    name[n-2] = '\0';
  }

  use_lfd = (strcmp(ext,"lfd")==0) | (strcmp(ext,"fit")==0);
  if (fp)
    sprintf(filename,"%s",name);
  else
  { if (use_lfd)
      sprintf(subp,"LFData%c",DIRSEP);
    else
      sprintf(subp,"");
    if (strlen(ext)==0)
      sprintf(filename,"%s%s",subp,name);
    else
      sprintf(filename,"%s%s.%s",subp,name,ext);
  }
  if (quote) name[n-2] = '"';

/*
 * If we are writing, check the file is writeable and that the
 * LFData directory exists.
 */
  if ((mode[0]=='w') | (mode[0]=='a'))
  { if (use_lfd)
    { if (access("LFData",F_OK)==-1)
      { if (access(".",W_OK)==0)
        { printf("Creating LFData Directory...\n");
          system("mkdir LFData");
        }
      }
      if (access("LFData",W_OK)==-1)
      { ERROR(("LFData directory not writeable"));
        return(0);
      }
    }
    return(1);  /* definitive test is whether fopen works. */
  }

/*
 *  If we are reading, check the file exists.
 *  If it doesn't and use_lfd is true, also check a defined lfhome.
 */
  if (mode[0]=='r')
  { if (access(filename,R_OK)==0) return(1);

    if ((use_lfd) && (lfhome!=NULL)) /* search system lfhome */
    { if (quote) name[n-2] = '\0';
      sprintf(filename,"%s/%s%s.%s",lfhome,subp,name,ext);
      if (quote) name[n-2] = '"';
      return(access(filename,R_OK)==0);
    }

    return(0);
  }
  ERROR(("setfilename: invalid mode %s",mode));
  return(0);
}

void readchar(c,n)
char *c;
INT n;
{ fread(c,1,n,lfd);
}

void readstr(z)
char *z;
{ while(1)
  { readchar(z,1);
    if (*z=='\0') return;
    z++;
  }
}

void dumpchar(c,n)
char *c;
INT n;
{ fwrite(c,1,n,lfd);
}

void dumpstr(z)
char *z;
{ dumpchar(z,strlen(z)+1);
}

#define LFDATAID -281

void dosavedata(v,fp)
vari *v;
int fp;
{ void (*fn)(), (*fs)();
  INT i, n;
  char *name;
  vari *v1;
  if (argarg(v,0)==NULL)
  { ERROR(("savedata: no filename"));
    return;
  }
  name = argarg(v,0);

  if (setfilename(name,"lfd","wb",fp)==0)
  { ERROR(("savedata: cannot access file %s",filename));
    return;
  }
  openfile("wb");
  if (lf_error) return;
  fn = dumpchar;
  fs = dumpstr;

  i = LFDATAID;
  (*fn)(&i, sizeof(INT));
  n = 0;
  for (i=1; i<v->n; i++) if (!argused(v,i))
  { v1 = findvar(argval(v,i),0,&n);
    if (v==NULL)
    { WARN(("variable %s not found; skipping",argval(v,i)));
    }
    else
    { (*fs)(v1->name);
      (*fn)(&v1->n,sizeof(INT));
      (*fn)(&v1->mode,sizeof(INT)); /* mode indicator for later */
      (*fn)(v1->dpr,v1->n*sizeof(double));
    }
    setused(v,i);
  }
  (*fs)("__end__");
  closefile();
}

void doreaddata(name,fp)
char *name;
int fp;
{ void (*fn)(), (*fs)();
  INT i, k, md, n, of;
  varname vn;
  vari *v;

  if (setfilename(name,"lfd","rb",fp)==0)
  { ERROR(("readdata: cannot access file %s",filename));
    return;
  }
  openfile("rb");
  if (lf_error) return;
  fn = readchar;
  fs = readstr;

  of = 0;
  (*fn)(&i, sizeof(INT));
  if (i!=LFDATAID) /* wrong or old format */
  { if (i==-367)
    { printf("Old format LFData file\n");
      of = 1;
    }
    else
    { ERROR(("not a Locfit data file: %s",name));
    }
  }
  if (lf_error) { closefile(); return; }

  if (of) /* old format: n nv name (10 char) data */
  { (*fn)(&n,sizeof(INT));
    (*fn)(&k,sizeof(INT));
    for (i=0; i<k; i++)
    { (*fn)(vn,10);
      v = createvar(vn,STREGULAR,n,VDOUBLE);
      (*fn)(v->dpr,n*sizeof(double));
    }
  }
  else /* new format: name (str) n mode data __end__ */
  { k = 999999;
    for (i=0; i<k; i++)
    { (*fs)(vn);
      if (strcmp(vn,"__end__")==0) i=999999;
      else
      { (*fn)(&n,sizeof(INT));
        (*fn)(&md,sizeof(INT));
        v = createvar(vn,STREGULAR,n,md);
        (*fn)(v->dpr,n*sizeof(double));
  } } }
  closefile();
}

#define FITID 4395943.3249934

void dosavefit(lf,fi,mode,fp)
lfit *lf;
char *fi, *mode;
int fp;
{ void (*fn)();
  double z;
  INT d = 0, i, k, lm, ld;

  if (fi==NULL) return;
  if (setfilename(fi,"fit",mode,fp)==0)
  { ERROR(("savefit: cannot access file %s.",fi));
    return;
  }

  if (mode[0]=='r')
    fn = readchar;
  else
  { if (lf->mi[MEV]==ENULL) ERROR(("savefit: No fit to save."));
    if (lf_error) return;
    fn = dumpchar;
    z = FITID;
    lm = LENM; ld = LEND;
    d = lf->mi[MDIM];
  }

  openfile(mode);
  (*fn)(&z,sizeof(double));

  if ((mode[0]=='r') && (z!=FITID))
  { ERROR(("readfit: file %s is not an evaluation structure",filename));
    closefile();
    return;
  }

  /* if reading, ensure lf.mi etc are assigned */
  if (mode[0]=='r') fitdefault(lf,0,1);

  (*fn)(&lm,sizeof(INT));
  (*fn)(lf->mi,lm*sizeof(INT));
  (*fn)(&ld,sizeof(INT));
  (*fn)(lf->dp,LEND*sizeof(double));
  (*fn)(&lf->nv,sizeof(INT));
  (*fn)(&lf->nce,sizeof(INT));
  (*fn)(&lf->vc,sizeof(INT));
  (*fn)(&lf->nnl,sizeof(INT)); /* no longer used -- delete sometime! */

  if (mode[0]=='r')
  { d = lf->mi[MDIM];
    trchck(lf,lf->nv,lf->nce,d,lf->mi[MP],lf->vc);
    pcchk(&lf->pc,d,lf->mi[MP],1);
    if ((mode[0]=='r') && (lm<20)) lf->mi[MPC] = 1-noparcomp(lf);
  }
  (*fn)(vdptr(lf->xxev),d*lf->nv*sizeof(double));
  for (i=0; i<3*lf->mi[MDIM]+8; i++)
    (*fn)(&lf->coef[i*lf->nvm],lf->nv*sizeof(double));

  for (i=0; i<d; i++) (*fn)(lf->xname[i],10);
  (*fn)(lf->yname,10);
  (*fn)(lf->bname,10);
  (*fn)(lf->cname,10);
  (*fn)(lf->wname,10);

  (*fn)(lf->sv,lf->nce*sizeof(double));
  (*fn)(lf->fl,2*d*sizeof(double));
  (*fn)(lf->sca,d*sizeof(double));
  (*fn)(lf->ce,lf->nce*lf->vc*sizeof(INT));
  (*fn)(lf->s,lf->nce*sizeof(INT));
  k = 0;
  if ((lf->mi[MEV]==EPHULL) | (lf->mi[MEV]==ETREE)) k = lf->nv;
  if (lf->mi[MEV]==EKDTR) k = lf->nce;
  (*fn)(lf->lo,k*sizeof(INT));
  (*fn)(lf->hi,k*sizeof(INT));
  (*fn)(lf->sty,d*sizeof(INT));
  if (lf->mi[MEV]==EGRID)
    (*fn)(lf->mg,d*sizeof(INT));
  (*fn)(&lf->nd,sizeof(INT));
  (*fn)(lf->deriv,lf->nd*sizeof(INT));

  (*fn)(vdptr(lf->pc.wk),pc_reqd(d,lf->mi[MP])*sizeof(double));
  lf->pc.xtwx.p = lf->mi[MP];
/* MUST save lf->pc.xtwx.sm here */
  lf->pc.xtwx.sm = lf->pc.xtwx.st = JAC_EIGD;

  closefile();
}

#endif
