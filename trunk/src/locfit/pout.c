/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 *
 *
 * print a plot structure in the format of various graph packages
 */

#include "local.h"

#ifdef CVERSION

#define XCON(x) (INT)(i0+(i1-i0)*(x-px[0])/(px[1]-px[0]))
#define YCON(y) (INT)(k0+(k1-k0)*(y-py[0])/(py[1]-py[0]))

extern double sin(), cos(), sqrt(), atan2(), ceil(), floor(), log10(), pow();
static INT i0, i1, k0, k1, cw, ch;
#ifdef YAWN
static FILE *plf;
#endif
extern INT curwin, lfcm[10];

#ifdef NOPIPES
FILE *popen(const char *cmd,const char *mode) { return(NULL); }
int pclose(FILE *file) { return(0); }
#endif

extern device devps, devwin;

char f2(fmt)
char *fmt;
{ char *z;
  z = strchr(fmt,',');
  if (z==NULL) return('d');
  z++;
  return(*z);
}

INT pretty(xl,k,z)
double *xl, *z;
INT k;
{ double dlt, m;
  INT i, j, u;
  if (k<=0) return(0);
  dlt = (xl[1]-xl[0])/k;
  m = floor(log10(dlt));
  dlt *= pow(10.0,-m);
  if (dlt<2) u = 2;
    else if (dlt<5) u = 5;
      else { u = 1; m++; } /* increments should be u*10^m; */
  dlt = u*pow(10.0,m);
  i = (INT)ceil(xl[0]/dlt);
  j = 0;
  while ((j<k) && ((i+j)*dlt<=xl[1]))
  { z[j] = (i+j)*dlt;
    j++;
  }
  return(j);
}

void isgrid(xyz)
plxyz *xyz;
{ INT i, n;
  vari *x, *y, *z;
  x = xyz->x; y = xyz->y; z = xyz->z;
  if ((x!=NULL) & (y!=NULL) & (z!=NULL))
  { if ((x->n*y->n)==z->n)
    { xyz->nx = x->n;
      xyz->ny = y->n;
      xyz->n = z->n;
      xyz->t = 1;
      return;
    }
    if ((x->n>1) & (y->n>1))
    { i = 0; n = z->n;
      while ((i<y->n) && (vitem(y,0)==vitem(y,i))) i++;
      if ((i>1) && (i*(n/i)==n))
      { xyz->nx = n/i;
        xyz->ny = i;
        xyz->n = n;
        xyz->t = 1;
        return;
      }
    }
  }
  xyz->t = 0;
  xyz->n = 0;
  if ((x!=NULL) && (x->n>xyz->n)) xyz->n = x->n;
  if ((y!=NULL) && (y->n>xyz->n)) xyz->n = y->n;
  if ((z!=NULL) && (z->n>xyz->n)) xyz->n = z->n;
}

void getxyzitem(x,xyz,i)
plxyz *x;
double xyz[3];
int i;
{ xyz[0] = vitem(x->x,i);
  if (x->t==1)
    xyz[1] = vitem(x->y,i/x->x->n);
  else
    xyz[1] = vitem(x->y,i);
  xyz[2] = vitem(x->z,i);
}

static double xl[2], yl[2], zl[2], px[2], py[2];

void project(z,x,theta,phi)
double *z, *x, theta, phi;
{ double z0, z1, z2;
  z0 = (z[0]-xl[0])/(xl[1]-xl[0]);
  z1 = (z[1]-yl[0])/(yl[1]-yl[0]);
  z2 = (z[2]-zl[0])/(zl[1]-zl[0]);
  x[0] = cos(theta)*z0-sin(theta)*z1;
  x[1] = (sin(theta)*z0+cos(theta)*z1)*cos(phi)+sin(phi)*z2;
}

void iproject(z,i,theta,phi)
double *z, theta, phi;
INT *i;
{ double x[2];
  project(z,x,theta,phi);
  i[0] = XCON(x[0]); i[1] = YCON(x[1]);
}

void line3d(z1,z2,theta,phi,dev,col)
double *z1, *z2, theta, phi;
device *dev;
INT col;
{ INT x1[2], x2[2];
  iproject(z1,x1,theta,phi);
  iproject(z2,x2,theta,phi);
  dev->SetColor(lfcm[col]);
  dev->DrawLine(x1[0],x1[1],x2[0],x2[1]);
}

void xyztext(tx,x,ah,av,theta,phi,dev,col)
double *x, theta, phi;
INT ah, av, col;
char *tx;
device *dev;
{ INT xy[2];
  iproject(x,xy,theta,phi);
  dev->SetColor(lfcm[col]);
  dev->DoText(0,xy[0],xy[1],tx,ah,av);
}

int getgreylevel(z)
double z;
{ int c;
  c = 8+11*(z-zl[0])/(zl[1]-zl[0]);
  if (c<8) return(8);
  if (c>18)return(18);
  return(c);
}

void points3d(x,theta,phi,dev,type)
plxyz *x;
double theta, phi;
char type;
device *dev;
{ INT i, xy[2];
  double xyz[3];
  for (i=0; i<x->n; i++)
  {
    getxyzitem(x,xyz,i);
    iproject(xyz,xy,theta,phi);
    if (type=='q') 
      dev->SetColor(getgreylevel(xyz[2]));
    else
      dev->SetColor(lfcm[CPOI]);
    dev->DrawPoint(xy[0],xy[1],x->pch);
  }
}

void lines3d(xyz,theta,phi,dev)
plxyz *xyz;
double theta, phi;
device *dev;
{ INT i;
  double x0[3], x1[3];
  getxyzitem(xyz,x0,0);
  for (i=1; i<xyz->n; i++)
  { if (i&1)
    { getxyzitem(xyz,x1,i);
      line3d(x0,x1,theta,phi,dev,CLIN);
    }
    else 
    { getxyzitem(xyz,x0,i);
      line3d(x1,x0,theta,phi,dev,CLIN);
    }
  }
}

void segments(xyz0,xyz1,theta,phi,dev)
plxyz *xyz0, *xyz1;
double theta, phi;
device *dev;
{ INT i, n;
  double x0[3], x1[3];
  n = xyz0->n;
  if (xyz1->n>n) n = xyz1->n;
  for (i=0; i<n; i++)
  { getxyzitem(xyz0,x0,i);
    getxyzitem(xyz1,x1,i);
    line3d(x0,x1,theta,phi,dev,CSEG);
  }
}

double spl(z0,z1)
double z0, z1;
{ if (z0-z1==0.0) return(0.5);
  return(z0/(z0-z1));
}

void contour3d(x,theta,phi,dev,sl,nsl)
plxyz *x;
double theta, phi, *sl;
INT nsl;
device *dev;
{ INT i, j, k, nx, ny, s;
  double xyz[4][3], u[4], x0[3], x1[3], x2[3], x3[3], r;
  char lb[20];
  if (x->t==0) ERROR(("Contour: not a grid"));
  if (lf_error) return;
  if (nsl==0) nsl = pretty(zl,10,sl);
  nx = x->nx; ny = x->ny;
  for (k=0; k<nsl; k++)
  { x0[2] = x1[2] = x2[2] = x3[2] = sl[k];
    for (i=0; i<nx-1; i++)
      for (j=0; j<ny-1; j++)
      { sprintf(lb,"%g",sl[k]);
        getxyzitem(x,xyz[0],i+j*nx);
        getxyzitem(x,xyz[1],(i+1)+j*nx);
        getxyzitem(x,xyz[2],i+(j+1)*nx);
        getxyzitem(x,xyz[3],i+1+(j+1)*nx);
        u[0] = xyz[0][2]-sl[k];
        u[1] = xyz[1][2]-sl[k];
        u[2] = xyz[2][2]-sl[k];
        u[3] = xyz[3][2]-sl[k];
        if (u[0]*u[1]<=0) /* bottom of cell */
        { r = spl(u[0],u[1]);
          x0[0] = (1-r)*xyz[0][0]+r*xyz[1][0];
          x0[1] = xyz[0][1];
          if (j==0) xyztext(lb,x0,-1,-1,theta,phi,dev,CCLA);
        }
        if (u[1]*u[3]<=0) /* right of cell */
        { r = spl(u[1],u[3]);
          x1[0] = xyz[1][0];
          x1[1] = (1-r)*xyz[1][1]+r*xyz[3][1];
          if (i==nx-2) xyztext(lb,x1,-1,0,theta,phi,dev,CCLA);
        }
        if (u[2]*u[3]<=0) /* top of cell */
        { r = spl(u[2],u[3]);
          x2[0] = (1-r)*xyz[2][0]+r*xyz[3][0];
          x2[1] = xyz[2][1];
          if (j==ny-2) xyztext(lb,x2,-1,1,theta,phi,dev,CCLA);
        }
        if (u[0]*u[2]<=0) /* left of cell */
        { r = spl(u[0],u[2]);
          x3[0] = xyz[0][0];
          x3[1] = (1-r)*xyz[0][1]+r*xyz[2][1];
          if (i==0) xyztext(lb,x3,-1,0,theta,phi,dev,CCLA);
        }

        s = 8*(u[3]>0)+4*(u[2]>0)+2*(u[1]>0)+(u[0]>0);
        switch(s)
        { case 0:
          case 15: break;
          case 1:
          case 14:
            line3d(x0,x3,theta,phi,dev,CCON);
            break;
          case 2:
          case 13:
            line3d(x0,x1,theta,phi,dev,CCON);
            break;
          case 3:
          case 12:
            line3d(x3,x1,theta,phi,dev,CCON);
            break;
          case 4:
          case 11:
            line3d(x3,x2,theta,phi,dev,CCON);
            break;
          case 5:
          case 10:
            line3d(x0,x2,theta,phi,dev,CCON);
            break;
          case 6:
          case 9:
            line3d(x0,x1,theta,phi,dev,CCON);
            break;
          case 7:
          case 8:
            line3d(x1,x2,theta,phi,dev,CCON);
            break;
          default: ERROR(("severe contouring error..."));
        }
      }
  }
}

double angle(x0,x1,x2) /* rotation angle from (x0,x1) to (x0,x2) */
double *x0, *x1, *x2;
/* If x0=0, ||x1=1|| then express
   x2 = u x1 + v y1       where y1=(-x11,x10) = 90 deg anticlk rot.
      i.e. u = <x1,x2>  v = <y1,x2>
   tan(theta) = v/u
   atan2(v,u) returns positive for anticlkws rot;
                      negative for clockwise rot.
*/
{ double u, v;
  u =  (x1[0]-x0[0])*(x2[0]-x0[0]) + (x1[1]-x0[1])*(x2[1]-x0[1]);
  v = -(x1[1]-x0[1])*(x2[0]-x0[0]) + (x1[0]-x0[0])*(x2[1]-x0[1]);
  return(atan2(v,u));
}

void persp3d(xyz,theta,phi,DP,sl,nsl)
plxyz *xyz;
double theta, phi, *sl;
void (*DP)();
INT nsl;
{ INT i, j, k, m, nx, ny, ii, jj, cb, cp, cx[4], cy[4], xhi, yhi;
  double u[4][3], w[4][2], r;
  if (xyz->t==0) ERROR(("persp3d: not a grid"));
  if (lf_error) return;
  /* theta=135 --- 225
            |       |
            45 --- 315;
     i.e start at top right for theta=45 e.t.c. Reverse if sin(phi)<0.
     x starts hi if cos(theta)>0
     y starts hi if sin(theta)>0
  */
  xhi = (cos(theta)*sin(phi)) > 0;
  yhi = (sin(theta)*sin(phi)) > 0;
  nx = xyz->nx; ny = xyz->ny;
  for (i=0; i<nx-1; i++)
    for (j=0; j<ny-1; j++)
    { for (k=0; k<4; k++)
      {  /*     1 -- 2
                |    |
                0 -- 3   */
        ii = (xhi) ? nx-2-i : i;
        jj = (yhi) ? ny-2-j : j;
        ii += (k>1);
        jj += (k==1)+(k==2);
        m = jj*nx+ii;
        getxyzitem(xyz,u[k],m);
        project(u[k],w[k],theta,phi);
        cx[k] = XCON(w[k][0]); cy[k] = YCON(w[k][1]);
      }

      switch(xyz->type)
      { case 'w':
          /* wireframe:
             from top, use color CPA1 for border, CPA2 for patch
             angles are anti-clock from top; r approx 2pi
                             clock from bot; r approx -2pi
          */
          r = angle(w[0],w[3],w[1])+angle(w[1],w[0],w[2])
             +angle(w[2],w[1],w[3])+angle(w[3],w[2],w[0]);
          if (r>0) { cb = lfcm[CPA1]; cp = lfcm[CPA2]; }
              else { cp = lfcm[CPA1]; cb = lfcm[CPA2]; }
          DP(cx,cy,cp,cb);
          break;
        case 'i': /* image */
          if (nsl==0) nsl = pretty(zl,10,sl);
          r = u[0][2] + u[1][2] + u[2][2] + u[3][2];
          cb = cp = getgreylevel(r/4);
          DP(cx,cy,cp,cb);
          break;
      }
    }
}

void updatelim(v,xl)
vari *v;
double *xl;
{ INT i;
  if ((v==NULL) || (vlength(v)==0)) return;
  if (xl[0]==xl[1])
    xl[0] = xl[1] = vitem(v,0);
  for (i=0; i<vlength(v); i++)
  { if (vitem(v,i)<xl[0]) xl[0] = vitem(v,i);
    if (vitem(v,i)>xl[1]) xl[1] = vitem(v,i);
  }
}

void axis(z1,z2,zl,lab,theta,phi,a,s,dev)
double *z1, *z2, *zl, theta, phi;
char *lab;
INT a, s;
device *dev;
{ double x1[3], x2[3], z[50];
  INT i, u0, u1, v0, v1, n, horiz;
  char lb[20];
  dev->SetColor(lfcm[CAXI]);
  project(z1,x1,theta,phi);
  project(z2,x2,theta,phi);
  u0 = XCON(x1[0]); v0 = XCON(x2[0]);
  u1 = YCON(x1[1]); v1 = YCON(x2[1]);
  horiz = abs(v0-u0)>abs(v1-u1);
  dev->DrawLine(u0,u1,v0,v1);
  n = (INT)sqrt((double)((v0-u0)*(v0-u0)+(v1-u1)*(v1-u1)))/((horiz) ? 5*cw : 2*ch);
  if (n>50) n = 50;
  n = pretty(zl,n,z);
  if (n==0) return;
  x1[0] = z1[0];  x1[1] = z1[1]; x1[2] = z1[2];
  if (abs(v0-u0)>abs(v1-u1)) /* horizontal axis */
  { dev->SetColor(lfcm[CTEX]);
    if (lab!=NULL)
      dev->DoText(0,(u0+v0)/2,(u1+v1)/2+s*(dev->ticklength+ch),lab,0,s);
    for (i=0; i<n; i++)
    { x1[a] = z[i];
      project(x1,x2,theta,phi);
      u0 = XCON(x2[0]); u1 = YCON(x2[1]);
      v0 = u0; v1 = u1+s*dev->ticklength;
      sprintf(lb,"%g",z[i]);
      dev->SetColor(lfcm[CAXI]);
      dev->DrawLine(u0,u1,v0,v1);
      dev->SetColor(lfcm[CTEX]);
      dev->DoText(0,v0,v1,lb,0,s);
  } }
  else /* vertical axis */
  { s = 2*((2*v0)>(i0+i1))-1;
    dev->SetColor(lfcm[CTEX]);
    if (lab!=NULL)
      dev->DoText(0,v0,v1-ch,lab,-s,-1);
    for (i=0; i<n; i++)
    { x1[a] = z[i];
      project(x1,x2,theta,phi);
      u0 = XCON(x2[0]); u1 = YCON(x2[1]);
      v0 = u0+s*dev->ticklength; v1 = u1;
      sprintf(lb,"%g",z[i]);
      dev->SetColor(lfcm[CAXI]);
      dev->DrawLine(u0,u1,v0,v1);
      dev->SetColor(lfcm[CTEX]);
      dev->DoText(0,v0,v1,lb,-s,0);
  } }
}

void plotxwin(pl,dev,wn,w,h,rd)
plots *pl;
device *dev;
INT wn, w, h, rd;
{ INT i, j, k, s;
  double z[3], z2[3], xx[2], vx, vy, vz;
  static double theta, phi;
  plxyz *xyz;
  if (pl->ty==PLNONE) return;
  if (h<=0) h = dev->defth;
  if (w<=0) w = dev->deftw;
  if (!dev->makewin(&w,&h,wn,rd)) return;
  dev->TextDim(0,"0",&cw,&ch);
  i0 = 4*cw+dev->ticklength;   i1 = w-2*cw;
  k0 = h-3*ch-dev->ticklength; k1 = 3*ch;
  dev->ClearScreen(lfcm[CBAK]);
  if (pl->xl[0]<pl->xl[1])
  { xl[0] = pl->xl[0]; xl[1] = pl->xl[1]; }
  else
  { xl[0] = xl[1] = 0.0;
    for (i=0; i<pl->xyzs->n; i++)
    { xyz = (plxyz *)viptr(pl->xyzs,i);
      updatelim(xyz->x,xl);
    }
    if (xl[0]==xl[1]) { xl[0] -= 0.5; xl[1] += 0.5; }
  }
  if (pl->yl[0]<pl->yl[1])
  { yl[0] = pl->yl[0]; yl[1] = pl->yl[1]; }
  else
  { yl[0] = yl[1] = 0.0;
    for (i=0; i<pl->xyzs->n; i++)
    { xyz = (plxyz *)viptr(pl->xyzs,i);
      updatelim(xyz->y,yl);
    }
    if (yl[0]==yl[1]) { yl[0] -= 0.5; yl[1] += 0.5; }
  }
  if (pl->zl[0]<pl->zl[1])
  { zl[0] = pl->zl[0]; zl[1] = pl->zl[1]; }
  else
  { zl[0] = zl[1] = 0.0;
    for (i=0; i<pl->xyzs->n; i++)
    { xyz = (plxyz *)viptr(pl->xyzs,i);
      updatelim(xyz->z,zl);
    }
    if (zl[0]==zl[1]) { zl[0] -= 0.5; zl[1] += 0.5; }
  }
  theta = pl->theta*PI/180; phi = pl->phi*PI/180;
  vx = -sin(theta)*sin(phi);
  vy = -cos(theta)*sin(phi);
  vz = cos(phi);

  for (i=0; i<2; i++)
    for (j=0; j<2; j++)
      for (k=0; k<2; k++)
      { z[0] = xl[i]; z[1] = yl[j]; z[2] = zl[k];
        project(z,xx,theta,phi);
        if ((i+j+k==0) | (xx[0]<px[0])) px[0]=xx[0];
        if ((i+j+k==0) | (xx[0]>px[1])) px[1]=xx[0];
        if ((i+j+k==0) | (xx[1]<py[0])) py[0]=xx[1];
        if ((i+j+k==0) | (xx[1]>py[1])) py[1]=xx[1];
      }
  s = 1-2*((cos(phi)<0)^(sin(phi)<0));
  z[0] = xl[0]; z2[0] = xl[1];
  z[1] = z2[1] = yl[(cos(theta)<0)^(sin(phi)<0)];
  z[2] = z2[2] = zl[cos(phi)<0];
  axis(z,z2,xl,pl->xlab,theta,phi,0,s,dev);
  z[0] = z2[0] = xl[(sin(theta)<0)^(sin(phi)<0)];
  z[1] = yl[0]; z2[1] = yl[1];
  z[2] = z2[2] = zl[cos(phi)<0];
  axis(z,z2,yl,pl->ylab,theta,phi,1,s,dev);
  z[0] = z2[0] = xl[cos(theta)<0];
  z[1] = z2[1] = yl[sin(theta)>0];
  z[2] = zl[0]; z2[2] = zl[1];
  axis(z,z2,zl,pl->zlab,theta,phi,2,s,dev);
  if (strlen(pl->main)>0) dev->DoText(1,(i0+i1)/2,2*ch,pl->main,0,-1);
 
  for (i=0; i<pl->xyzs->n; i++)
  { xyz = viptr(pl->xyzs,i);
    isgrid(xyz);
    switch (xyz->type)
    { case 'c':
        contour3d(xyz,theta,phi,dev,pl->sl,pl->nsl);
        break;
      case 'i':
        persp3d(xyz,theta,phi,dev->DrawPatch,pl->sl,pl->nsl);
        break;
      case 'b':
        points3d(xyz,theta,phi,dev,'p');
      case 'l':
        lines3d(xyz,theta,phi,dev);
        break;
      case 'p':
        points3d(xyz,theta,phi,dev,'p');
        break;
      case 'q':
        points3d(xyz,theta,phi,dev,'q');
        break;
      case 's':
        if (i==0) { ERROR(("invalid segements")); }
        else
        segments(viptr(pl->xyzs,i-1),xyz,theta,phi,dev);
        break;
      case 'w':
        persp3d(xyz,theta,phi,dev->DrawPatch,pl->sl,pl->nsl);
        break;
  } }

  dev->wrapup(rd);
}

#ifdef YAWN
void plotmaple(pl)
plots *pl;
{ INT i, j;
  plf = fopen("lfplot","w");
  switch(pl->d)
  { case 1:
      fprintf(plf,"PLOT(\n");
      for (j=0; j<pl->r; j++)
      { fprintf(plf,"CURVES([\n");
        for (i=0; i<pl->mx[j]; i++)
        { if (i>0) fprintf(plf,",\n");
          fprintf(plf,"[%f,%f]",pl->x[j][i],pl->y[j][i]);
        }
        fprintf(plf,"],\nCOLOUR(RGB,0,0,0)),\n");
      }
      if (pl->type[0]=='p') fprintf(plf,"STYLE(POINT),");
      if (pl->main!=NULL) fprintf(plf,"TITLE(%s),",pl->main);
      fprintf(plf,"AXESLABELS(%s,%s)",pl->xlab,pl->ylab);
      fprintf(plf,");\n");
      
      break;
    case 2:
      fprintf(plf,"PLOT3D(GRID(%f..%f,%f..%f,[[\n",pl->x[0][0],pl->x[0][pl->mx[0]-1],pl->y[0][0],pl->y[0][pl->my[0]-1]);
      for (i=0; i<pl->mx[0]; i++)
      { if (i>0) fprintf(plf,"],\n[");
        for (j=0; j<pl->my[0]; j++)
        { if (j>0) fprintf(plf,",\n");
          fprintf(plf,"%f",pl->z[0][i*pl->my[0]+j]);
        }
      }
      fprintf(plf,"]]),\nAXESLABELS(%s,%s,%s),AXESSTYLE(FRAME)",pl->xlab,pl->ylab,pl->zlab);
      if (pl->type[0]=='c') fprintf(plf,",STYLE(CONTOUR),CONTOURS(DEFAULT),ORIENTATION(-90,0.1),COLOUR(ZSHADING)");
      if (pl->main!=NULL) fprintf(plf,",\nTITLE(%s)\n",pl->main);
      fprintf(plf,");\n");
      break;
  }
  fclose(plf);
  printf("Created lfplot file; Maple format.\n");
}

void plotmathe(pl,fmt)
plots *pl;
char *fmt;
{ INT i, j, aut;
  static FILE *plm=NULL;
  aut = f2(fmt)!='m';
#ifdef NOPIPES
  aut = 0;
#endif
  if (aut)
  { if (plm==NULL) plm = (FILE *)popen("math >/dev/null","w");
    plf = plm;
  }
  else
    plf = fopen("lfplot","w");
  switch(pl->d)
  { case 1:
      fprintf(plf,"ListPlot[{{\n");
      for (i=0; i<pl->mx[0]; i++)
      { if (i>0) fprintf(plf,"},\n{");
        fprintf(plf,"%f,%f",pl->x[0][i],pl->y[0][i]);
      }
      fprintf(plf,"}}");
      fprintf(plf,",AxesLabel->{%s,%s}",pl->xlab,pl->ylab);
      if (pl->type[0]=='l') fprintf(plf,",PlotJoined->True");
      break;
    case 2:
      if (pl->type[0]=='c') fprintf(plf,"ListContourPlot[{{");
                    else fprintf(plf,"ListPlot3D[{{");
      for (j=0; j<pl->my[0]; j++)
      { if (j>0) fprintf(plf,"},\n{");
        for (i=0; i<pl->mx[0]; i++)
        { if (i>0) fprintf(plf,",\n");
          fprintf(plf,"%f",pl->z[0][i*pl->my[0]+j]);
        }
      }
      fprintf(plf,"}},\nMeshRange->{{");
      fprintf(plf,"%f,%f},{%f,%f}}\n",pl->x[0][0],pl->x[0][pl->mx[0]-1],pl->y[0][0],pl->y[0][pl->my[0]-1]);
      if (pl->type[0]=='c') fprintf(plf,",FrameLabel->{%s,%s}",pl->xlab,pl->ylab);
                    else fprintf(plf,",AxesLabel->{%s,%s,%s}",pl->xlab,pl->ylab,pl->zlab);
      break;
  }
  if (pl->main!=NULL) fprintf(plf,",PlotLabel->%s\n",pl->main);
  fprintf(plf,"];\n");
  if (aut)
    fflush(plf);
  else
  { fclose(plf);
    printf("Created lfplot file; Mathematica format.\n");
  }
}

void plotmatlb(pl) /* Matlab */
plots *pl;
{ INT i, j;
  plf = fopen("lfplot.m","w");
  switch(pl->d)
  { case 1:
      fprintf(plf,"plot([");
      for (i=0; i<pl->mx[0]; i++)
      { if (i>0) putc(',',plf);
        fprintf(plf,"%f",pl->x[0][i]);
      } 
      fprintf(plf,"],[");
      for (i=0; i<pl->mx[0]; i++)
      { if (i>0) putc(',',plf);
        fprintf(plf,"%f",pl->y[0][i]);
      }
      fprintf(plf,"])\n");
      break;
    case 2:
      if (pl->type[0]=='c') fprintf(plf,"contour([");
                    else fprintf(plf,"mesh([");
      for (i=0; i<pl->my[0]; i++) fprintf(plf,"%f ",pl->y[0][i]); 
      fprintf(plf,"],[");
      for (i=0; i<pl->mx[0]; i++) fprintf(plf,"%f ",pl->x[0][i]); 
      fprintf(plf,"],[\n");
      for (j=0; j<pl->my[0]; j++)
      { fprintf(plf,"[");
        for (i=0; i<pl->mx[0]; i++)
          fprintf(plf,"%f ",pl->z[0][i*pl->my[0]+j]);
        fprintf(plf,"]\n");
      }
      fprintf(plf,"])\n");
      fprintf(plf,"xlabel('%s')\n",pl->xlab);
      fprintf(plf,"ylabel('%s')\n",pl->ylab);
      break;
    default: ERROR(("plotmatlb: invalid dimension %d",pl->d));
  }
  if (pl->main!=NULL) fprintf(plf,"title('%s')\n",pl->main);
  fclose(plf);
  printf("Created lfplot.m file; matlab format.\n");
}

void plotgnup(pl,fmt)
plots *pl;
char *fmt;
{ INT i, j, z, aut;
  char m;
  static FILE *plc=NULL;

  /* first, the data file */
  plf=fopen("lfplot.dat","w");
  switch(pl->d)
  { case 1:
      z = pl->mx[0];
      for (j=0; j<pl->r; j++) if (pl->mx[j]>z) z = pl->mx[j];
      for (i=0; i<z; i++)
      { for (j=0; j<pl->r; j++)
          if (i<pl->mx[j])
            fprintf(plf,"%f %f ",pl->x[j][i],pl->y[j][i]);
          else
            fprintf(plf,"%f %f ",pl->x[j][pl->mx[j]-1],pl->y[j][pl->mx[j]-1]);
        fprintf(plf,"\n");
      }
      break;
    case 2:
      for (j=0; j<pl->my[0]; j++)
      { for (i=0; i<pl->mx[0]; i++)
          fprintf(plf,"%f %f %f\n",pl->x[0][i],pl->y[0][j],pl->z[0][i*pl->my[0]+j]);
        fprintf(plf,"\n");
      }
  }
  fclose(plf);

  /* second, the command file */
  m = f2(fmt);
  aut = (m!='m');
#ifdef NOPIPES
  aut = 0;
#endif
  if (aut)
  { if ((m=='s') && (plc!=NULL)) pclose(plc);
    if ((m=='s') || (plc==NULL))
      plc = (FILE *)popen("gnuplot","w");
    plf = plc;
  }
  else plf = fopen("lfplot","w");
  switch(pl->d)
  { case 1:
      fprintf(plf,"set nokey\n");
      fprintf(plf,"set xlab \"%s\"\n",pl->xlab);
      fprintf(plf,"set ylab \"%s\"\n",pl->ylab);
      if (pl->main != NULL)
        fprintf(plf,"set title \"%s\"\n",pl->main);
      fprintf(plf,"plot ");
      for (i=0; i<pl->r; i++)
      { if (i>0) fprintf(plf,", ");
        fprintf(plf,"\"lfplot.dat\" using %d:%d ",2*i+1,2*i+2);
        switch(pl->type[i])
        { case 'l': fprintf(plf,"with lines"); break;
          case 'p': fprintf(plf,"with points"); break;
          case 'b': fprintf(plf,"with linespoints"); break;
        }
      }
      fprintf(plf,"\n");
      break;
    case 2:
      fprintf(plf,"set xlab \"%s\"\n",pl->xlab);
      fprintf(plf,"set ylab \"%s\"\n",pl->ylab);
      fprintf(plf,"set zlab \"%s\"\n",pl->zlab);
      if (pl->type[0]=='c')
      { fprintf(plf,"set contour\n");
        fprintf(plf,"set nosurface\n");
        fprintf(plf,"set key\n");
      }
      else
      { fprintf(plf,"set nocontour\n");
        fprintf(plf,"set surface\n");
        fprintf(plf,"set nokey\n");
      }
      fprintf(plf,"set view %g,%g\n",pl->phi,pl->theta);
      fprintf(plf,"set parametric\n");
      if (pl->main != NULL)
        fprintf(plf,"set title \"%s\"\n",pl->main);
      fprintf(plf,"splot \"lfplot.dat\" with lines\n");
      break;
  }
  if (aut)
    fflush(plf);
  else
  { fclose(plf);
    printf("Created lfplot, lfplot.dat files; gnuplot format.\n");
  }
}
#endif

#endif
