/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

#include <ctype.h>
#include "local.h"

#ifdef CVERSION

extern lfit lf;
extern int ar_setfunction();

double vadd(double e1,double e2) { return(e1+e2); }
double vsub(double e1,double e2) { return(e1-e2); }
double vmul(double e1,double e2) { return(e1*e2); }
double vdiv(double e1,double e2) { return(e1/e2); }
double vpow(double e1,double e2)
{ if (e2==2) return(e1*e1);
    if (e1<=0) return(0.0);
    return(exp(e2*log(e1)));
}
double vgt(double e1,double e2) { return((double)(e1>e2)); }
double vlt(double e1,double e2) { return((double)(e1<e2)); }
double vge(double e1,double e2) { return((double)(e1>=e2)); }
double vle(double e1,double e2) { return((double)(e1<=e2)); }
double veq(double e1,double e2) { return((double)(e1==e2)); }
double vne(double e1,double e2) { return((double)(e1!=e2)); }

arstruct art;

double lf_exp(double x) { return (x<700.0) ? exp(x) : exp(700.0); }

double vseq(double a,double b,int i,int n) { return(a+(b-a)*i/(n-1)); }

double rsample(v)
vari *v;
{ int i;
    i = (int)( runif() * vlength(v) );
    return( vitem(v,i) );
}

vari *vrep(v1,v2)
vari *v1, *v2;
{ int i, j, k, m, n;
    vari *v;
    n = 0;
    for (i=0; i<vlength(v1); i++) n += vitem(v2,i);
    v = createvar("_vrep",STHIDDEN,n,VDOUBLE);
    k = 0;
    if (vlength(v2)==1)
    { m = vitem(v2,0);
        for (i=0; i<m; i++)
            for (j=0; j<vlength(v1); j++)
                vassn(v,k++,vitem(v1,j));
    }
    else
    { for (i=0; i<vlength(v1); i++)
    { m = vitem(v2,i);
        for (j=0; j<m; j++) vassn(v,k++,vitem(v1,i));
    }
    }
    deleteifhidden(v1);
    deleteifhidden(v2);
    return(v);
}

/*
 Set the arguments on arith structures.
 Argument number p (0,1,2; max args = 3)
 n is the index of source arstruct.
 k is the index of the argument arstruct.
 */
void setnext(va,n,p,k)
vari *va;
int n, p, k;
{ arstruct *rt;
    rt = viptr(va,n);
    if (p>=3)
    { ERROR(("Too many function arguments"));
        return;
    }
    rt->nx[p] = k;
}

void prars(v,i)
vari *v;
int i;
{ arstruct *ars;
    ars = (arstruct *)viptr(v,i);
    printf("%d %c\n",i,ars->cmd);
}

arstruct *cmdptr(v,i)
vari *v;
int i;
{ return((arstruct *)viptr(v,i));
}

int isstring(z,i1,i2)
char *z;
int i1, i2;
{ int i;
    if ((z[i1] != '"') | (z[i2] != '"')) return(0);
    for (i=i1+1; i<i2; i++) if (z[i]=='"') return(0);
    return(1);
}

int isname(z,i1,i2)
char *z;
int i1, i2;
{ int i;
    if (!isalpha(z[i1])) return(0);
    for (i=i1+1; i<=i2; i++)
        if (!isalnum(z[i])) return(0);
    return(1);
}

int isfunction(z,i1,i2)
char *z;
int i1, i2;
{ int i;
    if (z[i2] != ')') return(0);
    i = matchlf(z,i1,i2,'(',')');
    if (lf_error) return(0);
    return(isname(z,i1,i-1));
}

int issubset(z,i1,i2)
char *z;
int i1, i2;
{ int i;
    if (z[i2] != ']') return(0);
    i = matchlf(z,i1,i2,'[',']');
    if (lf_error) return(0);
    return(isname(z,i1,i-1));
}

int isNumber(z,i1,i2,val)
char *z;
int i1, i2;
double *val;
{ char *end;
    *val = strtod(&z[i1],&end);
    return(end==&z[i2+1]);
}

int isoperator(z,i1,i2)
char *z;
int i1, i2;
{ int i;
    
    i = checkrtol(z,i1,i2,",");
    if (i >= i1) return(i);
    
    i = checkrtol(z,i1,i2,"=<>");
    if (i > i1) return(i);
    
    i = checkrtol(z,i1,i2,"+-");
    while ((i>i1) && (strchr("+-*/:^",z[i-1])!=NULL))
        i = checkrtol(z,i1+1,i-1,"+-");
    if (i>i1) return(i);
    
    i = checkrtol(z,i1,i2,"*/");
    if (i >= i1) return(i);
    
    /* looks like a weird priority for : (sequence) but seems to match S. */
    i = checkrtol(z,i1,i2,":");
    if (i >= i1) return(i);
    
    i = checkrtol(z,i1,i2,"^");
    if (i >= i1) return(i);
    
    return(-1);
}

vari *arbuild(z,i1,i2,va,k,dum) /* k=0/1 search variable names? */
char *z;                        /* dum: are dummies x0, x1 etc allowed? */
INT i1, i2, k, dum;
vari *va;
{ INT al, ar, i, j, n, nargs, rargs, inum;
    vari *v;
    arstruct *rt;
    char tmp;
    double val;
    
    if (va==NULL)
    { va = createvar("_varith",STSYSTEM,10,VARC);
        if (lf_error || va == NULL) 
        {
            return(NULL);   
        }
        else
        {
            vlength(va) = 0;
        }
    }
    
    n = vlength(va);
    if (vbytes(n+1,VARC)>va->bytes) /* need to grow */
    { v = va;
        setvarname(va,"_ovarith");
        va = createvar("_varith",STSYSTEM,n+5,VARC);
        vlength(va) = n;
        memcpy(vdptr(va),vdptr(v),vbytes(n,VARC));
        deletevar(v);
    }
    inum = n;
    vlength(va) = n+1;
    
    while ((z[i1]=='(') && (matchrt(z,i1,i2,'(',')')==i2)) { i1++; i2--; }
    
    if (isNumber(z,i1,i2,&val))
    { rt = cmdptr(va,inum);
        rt->cmd = 'D';
        rt->x = val;
        return(va);
    }
    
    if (isstring(z,i1,i2))
    { rt = cmdptr(va,inum);
        rt->cmd = 's';
        rt->vv = createvar("_string",STHIDDEN,i2-i1,VCHAR);
        for (i=0; i<i2-i1-1; i++) ((char *)vdptr(rt->vv))[i] = z[i1+i+1];
        ((char *)vdptr(rt->vv))[i2-i1-1] = '\0';
        return(va);
    }
    
    if (isname(z,i1,i2))
    { tmp = z[i2+1]; z[i2+1] = '\0';
        if (dum) /* search for dummies */
        { for (j=0; j<lf.mi[MDIM]; j++)
            if (strcmp(&z[i1],lf.xname[j])==0)
            { rt = cmdptr(va,inum);
                rt->cmd = 'x';
                rt->m = j;
                z[i2+1] = tmp;
                return(va);
            }
        }
        n = 0;
        v = findvar(&z[i1],1,&n);
        z[i2+1] = tmp;
        if (v==NULL) return(va);
        rt = cmdptr(va,inum);
        rt->cmd = 'v';
        rt->vv = v;
        return(va);
    }
    
    if (isfunction(z,i1,i2))
    {
        /* build the argument list */
        ar = i2;
        al = matchlf(z,i1,i2,'(',')');
        j = al+1;
        nargs = 0;
        
        if (ar>al+1)
        { i = al;
            while (j<=ar)
            { if (z[j]=='(') j = matchrt(z,j,ar-1,'(',')')+1;
                if (z[j]=='[') j = matchrt(z,j,ar-1,'[',']')+1;
                if (lf_error) return(va);
                if ((z[j]==')') | (z[j]==','))
                { setnext(va,n,nargs,vlength(va));
                    va = arbuild(z,i+1,j-1,va,k,dum);
                    nargs++; i = j;
                }
                j++;
            }
        }
        rt = cmdptr(va,inum);
        rt->m = nargs;
        
        rargs = ar_setfunction(rt,&z[i1],al-i1);
        if (rargs != nargs)
            ERROR(("arbuild: wrong number of arguments, %s",&z[i1]));
        return(va);
    }
    
    rt = cmdptr(va,inum);
    
    if (issubset(z,i1,i2))
    { rt->cmd = 'U';
        al = matchlf(z,i1,i2,'[',']');
        setnext(va,n,0,vlength(va));
        va = arbuild(z,i1,al-1,va,k,dum);
        if (lf_error) return(va);
        setnext(va,n,1,vlength(va));
        va = arbuild(z,al+1,i2-1,va,k,dum);
        return(va);
    }
    
    /* that leaves operators */
    
    i = isoperator(z,i1,i2);
    if (i >= i1)
    { rt->cmd = 'O';
        rt->f = NULL;
        al = i-1; ar = i+1;
        if (z[i]==',') rt->cmd = 'C';
        if (z[i]=='>')
        { rt->f = vgt;
            if (z[i-1]=='<') { rt->f = vne; al--; }
        }
        if (z[i]=='<') rt->f = vlt;
        if (z[i]=='=')
        { rt->f = veq;
            if (z[i-1]=='=') al--;
            if (z[i-1]=='<') { rt->f = vle; al--; }
            if (z[i-1]=='>') { rt->f = vge; al--; }
            if (z[i-1]=='!') { rt->f = vne; al--; }
        }
        if (z[i]=='+') rt->f = vadd;
        if (z[i]=='-') rt->f = vsub;
        if (z[i]=='*') rt->f = vmul;
        if (z[i]=='/') rt->f = vdiv;
        if (z[i]==':') rt->cmd = ':';
        if (z[i]=='^') rt->f = vpow;
        
        setnext(va,n,0,vlength(va));
        va = arbuild(z,i1,al,va,k,dum);
        if (lf_error) return(va);
        setnext(va,n,1,vlength(va));
        va = arbuild(z,ar,i2,va,k,dum);
        return(va);
    }
    
    ERROR(("arbuild: unknown expression %s",z));
    return(va);
}

vari *vevop(l,r,f)
vari *l, *r;
double (*f)();
{ INT i, n;
    vari *v;
    double z;
    if ((l==NULL) | (r==NULL)) return(NULL);
    n = vlength(l);
    if (n<vlength(r)) n = vlength(r);
    v = createvar("_vevop",STHIDDEN,n,VDOUBLE);
    if (lf_error) return(NULL);
    for (i=0; i<n; i++)
    { z = f(vitem(l,i),vitem(r,i));
        vassn(v,i,z);
    }
    deleteifhidden(l);
    deleteifhidden(r);
    return(v);
}

vari *vcat(v1,v2)
vari *v1, *v2;
{ INT i, n;
    vari *v;
    n = vlength(v1) + vlength(v2);
    v = createvar("_vcat",STHIDDEN,n,VDOUBLE);
    if (lf_error) return(NULL);
    for (i=0; i<vlength(v1); i++) vassn(v,i,vitem(v1,i));
    for (i=0; i<vlength(v2); i++) vassn(v,i+vlength(v1),vitem(v2,i));
    deleteifhidden(v1);
    deleteifhidden(v2);
    return(v);
}

vari *vrvec(v1,f)
vari *v1;
double (*f)();
{ vari *v;
    INT i, n;
    n = (INT)vitem(v1,0);
    v = createvar("_vrvec",STHIDDEN,n,VDOUBLE);
    if (lf_error) return(NULL);
    for (i=0; i<n; i++) vassn(v,i,f());
    deleteifhidden(v1);
    return(v);
}

vari *vrsca(v,f)
vari *v;
double (*f)();
{ vari *z;
    if (v==NULL) return(NULL);
    z = createvar("_vrsca",STHIDDEN,1,VDOUBLE);
    if (lf_error) return(NULL);
    vassn(z,(INT)0,f(v));
    deleteifhidden(v);
    return(z);
}

vari *vrve2(v1,a1,f)
vari *v1, *a1;
double (*f)();
{ vari *v;
    int i, n;
    n = vitem(v1,0);
    v = createvar("_vrvec",STHIDDEN,n,VDOUBLE);
    if (lf_error) return(NULL);
    for (i=0; i<n; i++) vassn(v,i,f(vitem(a1,i)));
    //deleteifhidden(v1,a1);
    return(v);
}

vari *vvec1(l,f)
vari *l;
double (*f)();
{ vari *v;
    INT i;
    if (l==NULL) 
    {
        ERROR(("vvec1 recieved NULL variable\n"));  
        return NULL;
    }
    v = createvar("_vvec1",STHIDDEN,l->n,VDOUBLE);
    if (lf_error) return(NULL);
    for (i=0; i<vlength(l); i++) vassn(v,i,f(vitem(l,i)));
    deleteifhidden(l);
    return(v);
}

vari *vrsamp(v1,v2)
vari *v1, *v2;
{ vari *v;
    int i, n;
    n = vitem(v2,0);
    v = createvar("_vrsamp",STHIDDEN,n,VDOUBLE);
    if (lf_error) return(NULL);
    for (i=0; i<n; i++) vassn(v,i,rsample(v1));
    deleteifhidden(v1);
    deleteifhidden(v2);
    return(v);
}

vari *vrseq(v1,v2,v3)
vari *v1, *v2, *v3;
{ vari *v;
    double beg, end;
    int i, n;
    n = vitem(v3,0);
    v = createvar("_vrseq",STHIDDEN,n,VDOUBLE);
    if (lf_error) return(NULL);
    beg = vitem(v1,0);
    end = vitem(v2,0);
    for (i=0; i<n; i++) vassn(v,i,vseq(beg,end,i,n));
    deleteifhidden(v1);
    deleteifhidden(v2);
    deleteifhidden(v3);
    return(v);
}

vari *vrse2(v1,v2)
vari *v1, *v2;
{ vari *v;
    double beg, end;
    int i, n;
    beg = vitem(v1,0);
    end = vitem(v2,0);
    n = (beg<=end) ? end-beg+1 : beg-end+1;
    v = createvar("_vrse2",STHIDDEN,n,VDOUBLE);
    if (lf_error) return(NULL);
    if (beg<=end) for (i=0; i<n; i++) vassn(v,i,beg+i);
    else for (i=0; i<n; i++) vassn(v,i,beg-i);
    deleteifhidden(v1);
    deleteifhidden(v2);
    return(v);
}

vari *vsubset(v1,v2)
vari *v1, *v2;
{ vari *v;
    int i, n;
    n = v2->n;
    v = createvar("_vsubs",STHIDDEN,n,VDOUBLE);
    for (i=0; i<n; i++) vassn(v,i,vitem(v1,vitem(v2,i)-1));
    deleteifhidden(v1);
    deleteifhidden(v2);
    return(v);
}

vari *vareval(v,k)
vari *v;
INT k;
{ INT n;
    arstruct *rt;
    rt = viptr(v,k);
    switch (rt->cmd)
    { case 'e': return(NULL);
        case 'v': return(rt->vv);
        case 'O': return(vevop(vareval(v,rt->nx[0]),vareval(v,rt->nx[1]),rt->f));
        case 'C': return(vcat(vareval(v,rt->nx[0]),vareval(v,rt->nx[1])));
        case 'D':
            n = 1;
            rt->vv = createvar("_vevcon",STHIDDEN,n,VDOUBLE);
            if (lf_error) return(NULL);
            vassn(rt->vv,0,rt->x);
            return(rt->vv);
        case 'G': return(vrvec(vareval(v,rt->nx[0]),rt->f));
        case 'H': return(vrve2(vareval(v,rt->nx[0]),vareval(v,rt->nx[1]),rt->f));
        case 'f': return(vvec1(vareval(v,rt->nx[0]),rt->f));
        case 'M': return(vrsamp(vareval(v,rt->nx[0]),vareval(v,rt->nx[1])));
        case 'Q': return(vrseq(vareval(v,rt->nx[0]),vareval(v,rt->nx[1]),
                               vareval(v,rt->nx[2])));
        case ':': return(vrse2(vareval(v,rt->nx[0]),vareval(v,rt->nx[1])));
        case 'S': return(vrsca(vareval(v,rt->nx[0]),rt->f));
        case 'U': return(vsubset(vareval(v,rt->nx[0]),vareval(v,rt->nx[1])));
        case 'R': return(vrep(vareval(v,rt->nx[0]),vareval(v,rt->nx[1])));
        case 'Z': return(vfitted(RMEAN));
        case 'Y': return(vfitted(RDEV));
        case 's': return(rt->vv);
        case 'x':
            ERROR(("Dummy in vareval"));
            return(NULL);
        default : ERROR(("vareval: unknown command %c",rt->cmd));
    }
    return(NULL);
}

vari *saveresult(v,name,status)
vari *v;
int status;
char *name;
{ vari *vr;
    if (v==NULL) return(NULL);
    
    vr = v;
    if (v->stat != STHIDDEN)
    { vr = createvar("_result",STHIDDEN,vlength(v),vmode(v));
        memcpy(vdptr(vr),vdptr(v),vbytes(vlength(v),vmode(v)));
    }
    
    if (name!=NULL)
    { setvarname(vr,name);
        vr->stat = status;
    }
    return(vr);
}

vari *varith(z,name,status)
char *z, *name;
int status;
{ vari *v, *va;
    va = arbuild(z,0,strlen(z)-1,NULL,1,0);
    if (lf_error) return(NULL);
    v = vareval(va,0);
    deletevar(va);
    
    v = saveresult(v,name,status);
    return(v);
}

double dareval(v,k,x)
vari *v;
INT k;
double *x;
{ arstruct *rt;
    rt = viptr(v,k);
    switch (rt->cmd)
    { case 'e': return(0.0);
        case 'v': return(vitem(rt->vv,0));
        case 'O': return(rt->f(dareval(v,rt->nx[0],x),dareval(v,rt->nx[1],x)));
        case 'P': return(rt->f(0.0,dareval(v,rt->nx[1],x)));
        case 'D': return(rt->x);
        case 'G': return(rt->f());
        case 'H': return(rt->f(dareval(v,rt->nx[1],x)));
        case 'f': return(rt->f(dareval(v,rt->nx[0],x)));
        case 'M': return(rsample(vareval(v,rt->nx[0])));
        case 'x': return(x[rt->m]);
        case 'U': return(vitem(vareval(v,rt->nx[0]),(int)dareval(v,rt->nx[1],x)-1));
        case 'Q': ERROR(("sequence in dareval"));
            return(0.0);
        default : ERROR(("dareval: unknown command %c",rt->cmd));
    }
    return(0.0);
}

double darith(z)
char *z;
{ vari *va;
    double y;
    va = arbuild(z,0,strlen(z)-1,NULL,1,0);
    y = dareval(va,0,NULL);
    deletevar(va);
    return(y);
}

INT arvect(z,res,c,a) /* c = no of items to read */
char *z;
INT c, a;
double *res;
{ INT i;
    vari *v;
    
    if (z==NULL) return(0);
    
    v = varith(z,"arvect",STPLOTVAR);
    if (v==NULL || lf_error) 
    {
        return(0);   
    }
    deletevar(v);
    
    for (i=0; (i<c) & (i<vlength(v)); i++) res[i] = vitem(v,i);
    if (i<c) /* insufficient items */
    { switch(a)
        { case 0: /* pad */
                for (; i<c; i++) res[i] = res[0];
                return(c);
            case 1: return(i); /* ignore */
            case 2:
                ERROR(("Insufficient items: %s",z));
                return(i);
        }
    }
    return(i);
}

#endif
