/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

vari *growvar(vari* vold, int n);

plots *cpl, pl[MAXWIN];
extern device devps, devwin;
INT curwin;
char *psfn;
extern lfit lf;
extern pplot pp;
extern char *curstr;

plots *get_graphics_window(v)
vari *v;
{ int win_no;
  char *w;

  w = getargval(v,"win",0);
  if (w != NULL)
  { sscanf(w,"%d",&win_no);
    if ((win_no<0) | (win_no>=MAXWIN))
    { WARN(("Invalid window %d",win_no));
    }
    else
      curwin = win_no;
  }
  return( &pl[curwin] );
}

char *settype(xyz,type,def)
plxyz *xyz;
char *type, def;
{ if ((type==NULL) || (strlen(type)==0))
  { xyz->type = def;
    return(NULL);
  }
  xyz->type = type[0];
  return(&type[1]);
}

char *pvarname(xyz,c,vn)
plxyz *xyz;
char c;
varname vn;
{ sprintf(vn,"_plv%d%c",xyz->id,c);
  return(vn);
}

plxyz *nextxyz(win,add,ck)
plots *win;
INT add, ck;
{ plxyz *xyz;
  vari *v;
  varname vn;

  if (!add)
  { sprintf(vn,"_xyz%d",curwin);
    v = win->xyzs = createvar(vn,STSYSTEM,5,VXYZ);
    v->n = 0;
    win->xlab[0] = win->ylab[0] = win->zlab[0] = '\0';
  }
  else
    v = win->xyzs = growvar(win->xyzs,vlength(win->xyzs)+ck);

  xyz = (plxyz *)viptr(v,vlength(v));
  xyz->id = (vlength(v) << 4) + win->id;
  xyz->pch = 1;
  v->n++;
  return(xyz);
}

void plotopt(v,re)
vari *v;
INT re;
{ INT i, j, h, w;
  double z[2];
  char *fmt, *ty;
  plxyz *xyz;

  cpl = get_graphics_window(v);

  if (!re)
  { cpl->main[0] = '\0';
    cpl->xl[0] = cpl->xl[1] = cpl->yl[0]
               = cpl->yl[1] = cpl->zl[0] = cpl->zl[1] = 0.0;
    cpl->nsl = 0;
  }

  arvect(getargval(v,"xlim",1), cpl->xl, 2, 2);
  arvect(getargval(v,"ylim",1), cpl->yl, 2, 2);
  arvect(getargval(v,"zlim",1), cpl->zl, 2, 2);

  i = getarg(v,"main",1);
  if (i>0) { strcpy(cpl->main,argval(v,i)); strip(cpl->main); }
  i = getarg(v,"xlab",1);
  if (i>0) { strcpy(cpl->xlab,argval(v,i)); strip(cpl->xlab); }
  i = getarg(v,"ylab",1);
  if (i>0) { strcpy(cpl->ylab,argval(v,i)); strip(cpl->ylab); }
  i = getarg(v,"zlab",1);
  if (i>0) { strcpy(cpl->zlab,argval(v,i)); strip(cpl->zlab); }

  if ( arvect(getargval(v,"view",1), z, 2, 2) == 2 )
  { cpl->theta=z[0];
    cpl->phi  =z[1];
  }

  fmt = getargval(v,"fmt",1);
  if (fmt==NULL) fmt = "xwin";

  i = getarg(v,"split",1);
  if (i>0)
    cpl->nsl = arvect(argval(v,i),cpl->sl,10,1);
  i = getarg(v,"h",1);
  if (i>0) sscanf(argval(v,i),"%d",&h); else h = 0;
  i = getarg(v,"w",1);
  if (i>0) sscanf(argval(v,i),"%d",&w); else w = 0;

  ty = getargval(v,"type",1);
  if (ty != NULL)
  { for (j=0; j<vlength(cpl->xyzs); j++)
    { xyz = (plxyz *)viptr(cpl->xyzs,j);
      ty = settype(xyz,ty,xyz->type);
    }
  }
  if (stm(fmt,"xwin",1)) { plotxwin(cpl,&devwin,curwin,w,h,0); return; }
  if (stm(fmt,"win",1))  { plotxwin(cpl,&devwin,curwin,w,h,0); return; }
  if (stm(fmt,"post",1))
  { psfn = getargval(v,"file",1);
    plotxwin(cpl,&devps,curwin,w,h,0);
    return;
  }
}

void pvari(cmd,xyz,win,ax)
char *cmd, ax;
plxyz *xyz;
plots *win;
{ vari *vv;
  INT k;
  varname vname;
  vv = varith(cmd,pvarname(xyz,ax,vname),STPLOTVAR);
  if (vv==NULL) return;
  k = xyz->id>>4;
  switch(ax)
  { case 'x':
      xyz->x = vv;
      if (k==0) strcpy(win->xlab,cmd);
      return;
    case 'y':
      xyz->y = vv;
      if (k==0) strcpy(win->ylab,cmd);
      return;
    case 'z':
      xyz->z = vv;
      if (k==0) strcpy(win->zlab,cmd);
      return;
  }
  ERROR(("pvari: unknown axis %c",ax));
}

void plotdata(v)
vari *v;
{ INT add, i, j, k;
  plxyz *xyz = NULL, *xyz2 = NULL;
  char *type;

  cpl = get_graphics_window(v);

  i = getarg(v,"add",0);
  add = (i>0) ? getlogic(v,i) : 0;

  type = getargval(v,"type",0);

  i = getarg(v,"data",0);
  if (i>0) doreaddata(argval(v,i),(INT)0);

  i = getarg(v,"pch",0);
  if (i>0) sscanf(argval(v,i),"%d",&xyz->pch);

  xyz = nextxyz(cpl,add,2);
  if (xyz==NULL) return;
  xyz->x = xyz->y = xyz->z = NULL;
  type = settype(xyz,type,'p');

  i = getarg(v,"x",1);
  j = getarg(v,"y",1);
  k = getarg(v,"z",1);

  if (!add) /* set the default view angle */
  { cpl->theta = 45*( 1 - ((j==0)|(k==0)) - 3*(i==0) );
    cpl->phi   = 45*( 1 + ((i==0)|(j==0)) - (k==0) );
  }

  if (i>0) pvari(argval(v,i),xyz,cpl,'x');
  if (j>0) pvari(argval(v,j),xyz,cpl,'y');
  if (k>0) pvari(argval(v,k),xyz,cpl,'z');

  i = getarg(v,"x2",1);
  j = getarg(v,"y2",1);
  k = getarg(v,"z2",1);
  if (i+j+k>0)
  { xyz2= nextxyz(cpl,1,1);
    if (xyz2==NULL) return;
    xyz2->x = xyz->x;
    xyz2->y = xyz->y;
    xyz2->z = xyz->z;
    type = settype(xyz2,type,'s');
    if (i>0) pvari(argval(v,i),xyz2,cpl,'x');
    if (j>0) pvari(argval(v,j),xyz2,cpl,'y');
    if (k>0) pvari(argval(v,k),xyz2,cpl,'z');
  }
  if (lf_error) return;

  cpl->ty |= PLDATA;

  plotopt(v,add);
}

void plotfit(v)
vari *v;
{ INT add, d, dp, i = 0, j = 0, n, sef;
  INT dt, mg[MXDIM], ix, iy;
  double c, sd, xl[2*MXDIM], xll[2];
  char cb;
  varname vn;
  plxyz *xyz, *xyzl, *xyzu, *xyzd;
  char *type;

  cpl = get_graphics_window(v);

  i = getarg(v,"fit",1);
  if (i>0) dosavefit(&lf,argval(v,i),"rb",(INT)0);
  if (nofit()) ERROR(("plotfit: no fit to plot."));
  if (lf_error) return;
  dp = 0;

  d = lf.mi[MDIM];
  for (i=0; i<d; i++)
  { j = getarg(v,lf.xname[i],0);
    if (j>0)
    { j = arvect(argval(v,j),xll,2,1);
      if (j==1)
        xl[i] = xl[i+d] = xll[0];
      else
      { xl[i] = xll[0];
        xl[i+d] = xll[1];
      }
    }
    else
    { xl[i] = lf.fl[i];
      xl[i+d] = lf.fl[i+d];
      j = 2;
    }
    if (j==2)
    { if (dp==2)
      { xl[i] = xl[i+d] = (xl[i]+xl[i+d])/2;
        WARN(("plotfit: fixing %s=%f",lf.xname[i],xl[i]));
        j = 1;
      }
      if (dp==1) { iy = i; dp++; }
      if (dp==0) { ix = i; dp++; }
    }
    mg[i] = 2-j;
  }
  if (dp<=0)
  { ERROR(("No plot variables"));
    return;
  }
  sef = 0; dt = 0;
  i = getarg(v,"data",1);  if (i>0) dt =getlogic(v,i);
  i = getarg(v,"band",1);  cb = (i>0) ? *argval(v,i) : 'n';

  for (i=0; i<d; i++) mg[i] = 1;
  if (dp==1)
    mg[ix] = 100;
  else
    mg[ix] = mg[iy] = 50;
  i = getarg(v,"m",1);
  if (i>0) readilist(mg,argval(v,i),1,lf.mi[MDIM],1);

  i = getarg(v,"add",1);
  add = (i>0) ? getlogic(v,i) : 0;

  type =  getargval(v,"type",1);

  if ((lf.mi[MEV]==EDATA) | (lf.mi[MEV]==ECROS))
    n = setpppoints(&pp,"fitp",mg,xl);
  else
    n = setpppoints(&pp,"grid",mg,xl);
  pp.fit = createvar("_ppfit",STPLOTVAR,n,VDOUBLE);
  if (cb=='n')
    pp.se = NULL;
  else
    pp.se = createvar("_ppsef",STPLOTVAR,n,VDOUBLE);
  if (lf_error) return;
  cpreplot(&pp,v,cb);
  if (lf_error) return;

  xyz = nextxyz(cpl,add,4);
  if (xyz==NULL) return;
  /* set up first predictor variable */
  xyz->x = pp.data[ix];
  setvarname(xyz->x,pvarname(xyz,'x',vn));
  strcpy(cpl->xlab,lf.xname[ix]);

  /* set up second predictor variable */
  if (dp==2)
  { xyz->y = pp.data[iy];
    setvarname(xyz->y,pvarname(xyz,'y',vn));
    strcpy(cpl->ylab,lf.xname[iy]);
  }
  else
  { xyz->y = NULL;
    cpl->ylab[0] = '\0';
  }

  xyz->z = pp.fit;
  setvarname(xyz->z,pvarname(xyz,'z',vn));
  switch(lf.mi[MTG]&63)
  { case TDEN: strcpy(cpl->zlab,"Density"); break;
    case TRAT: strcpy(cpl->zlab,"Rate"); break;
    case THAZ: strcpy(cpl->zlab,"Hazard"); break;
    default: strcpy(cpl->zlab,lf.yname);
  }
  type = settype(xyz,type,(dp==1) ? 'l' : 'c');

  if (pp.se!=NULL)
  { xyzl = nextxyz(cpl,1,3); xyzu = nextxyz(cpl,1,2);
    if ((xyzl!=NULL) & (xyzu!=NULL))
    { sd = sqrt(lf.dp[DRV]);
      xyzl->x = xyzu->x = xyz->x;
      xyzl->y = xyzu->y = xyz->y;
      xyzl->z = createvar(pvarname(xyzl,'z',vn),STPLOTVAR,n,VDOUBLE);
      xyzu->z = createvar(pvarname(xyzu,'z',vn),STPLOTVAR,n,VDOUBLE);
      if (lf_error) return;
      c = docrit(v);
      for (i=0; i<n; i++)
      { vassn(xyzu->z,i,backtr(vitem(pp.fit,i)+c*vitem(pp.se,i),lf.mi,lf.nd));
        vassn(xyzl->z,i,backtr(vitem(pp.fit,i)-c*vitem(pp.se,i),lf.mi,lf.nd));
      }
      type = settype(xyzl,type,(d==1) ? 'l' : 'c');
      type = settype(xyzu,type,(d==1) ? 'l' : 'c');
    }
    deletevar(pp.se);
  }
  if (pp.wh==PCOEF)
    for (i=0; i<vlength(xyz->z); i++)
      vassn(xyz->z,i,backtr(vitem(pp.fit,i),lf.mi,lf.nd));
  if (dt)
  {
    recondat(0,&n);
    if (lf_error) return;
    xyzd = nextxyz(cpl,1,1);
    if (xyzd!=NULL)
    { xyzd->x = createvar(pvarname(xyzd,'x',vn),STPLOTVAR,n,VDOUBLE);
      for (i=0; i<n; i++) vassn(xyzd->x,i,datum(&lf,ix,i));
      if (d==2)
      { xyzd->y = createvar(pvarname(xyzd,'y',vn),STPLOTVAR,n,VDOUBLE);
        for (i=0; i<n; i++) vassn(xyzd->y,i,datum(&lf,iy,i));
      }
      else xyzd->y = NULL;
      xyzd->z = createvar(pvarname(xyzd,'z',vn),STPLOTVAR,n,VDOUBLE);
      for (i=0; i<n; i++)
        vassn(xyzd->z,i,((lf.mi[MTG]&63)==TGAUS) ? resp(&lf,i) : resp(&lf,i)/prwt(&lf,i));
      type = settype(xyzd,type,'p');
    }
  }

  /* now, set default view angle */
  if (!add)
  { if (dp==1) { cpl->theta = 0; cpl->phi = 90; } /* x-z axis */
    else
    { if (xyz->type=='w')
        cpl->theta = cpl->phi = 45; /* wireframes */
      else
        cpl->theta = cpl->phi = 0;  /* x-y plot; e.g. for contours */
    }
  }
  if (lf_error) return;
  cpl->ty |= PLFIT;
  if (dt) cpl->ty |= PLDATA;

  plotopt(v,add);
}

void plottrack(v)
vari *v;
{ INT i, j;
  plxyz *xyz;
  varname vn;

  cpl = get_graphics_window(v);

  if ((cpl->ty & PLTRK)!=PLTRK) /* initialize */
  { xyz = nextxyz(cpl,0,1);
    xyz->x = createvar(pvarname(xyz,'x',vn),STPLOTVAR,100,VDOUBLE);
    xyz->y = createvar(pvarname(xyz,'y',vn),STPLOTVAR,100,VDOUBLE);
    xyz->z = createvar(pvarname(xyz,'z',vn),STPLOTVAR,100,VDOUBLE);
    if (lf_error) return;
    vlength(xyz->x) = vlength(xyz->y) = vlength(xyz->z) = 0;
    settype(xyz,NULL,'p');
    cpl->theta = cpl->phi = 0;
    cpl->ty = PLTRK;
  }
  else
  { vlength(cpl->xyzs) = 0;
    xyz = nextxyz(cpl,1,1);
  }
  j = vlength(xyz->x);
  i = getarg(v,"x",1);
  if (i>0)
  { vassn(xyz->x,j,darith(argval(v,i)));
    strcpy(cpl->xlab,argval(v,i));
    vlength(xyz->x) = j+1;
  }
  i = getarg(v,"y",1);
  if (i>0)
  { vassn(xyz->y,j,darith(argval(v,i)));
    strcpy(cpl->ylab,argval(v,i));
    vlength(xyz->y) = j+1;
  }
  i = getarg(v,"z",1);
  if (i>0)
  { vassn(xyz->z,j,darith(argval(v,i)));
    strcpy(cpl->zlab,argval(v,i));
    vlength(xyz->z) = j+1;
  }
  plotopt(v,0);
}

void setplot(v)
vari *v;
{ INT i, j, w;
  carg *ct;
  varname tname;
  i = getarg(v,"win",1);
  if (i==0)
  { ERROR(("setplot: no win argument"));
    return;
  }
  sscanf(argval(v,i),"%d",&w);
  if ((w<0) | (w>=MAXWIN))
  { ERROR(("setplot: invalid win %s",argval(v,i)));
    return;
  }
  if (vlength(v)==2)
  { deletevar(pl[w].track);
    pl[w].track = NULL;
    return;
  }
  sprintf(tname,"=tpc%d",w);
  pl[w].track = createvar(tname,STSYSTEM,v->n-2,VARGL);
  j = 0;
  pl[w].ty = PLNONE; /* to ensure previous track is cleared */
  for (i=1; i<vlength(v); i++)
  { if (!argused(v,i))
    { ct = (carg *)viptr(pl[w].track,j);
      ct->arg = argarg(v,i);
      ct->val = argval(v,i);
      setused(v,i);
      j++;
    }
  }
  sprintf(tname,"=tps%d",w);
  setvarname(curstr,tname);
}
