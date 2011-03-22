/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 *
 *  Functions for handling locfit variables in the C version.
 */

extern "C" 
{
#include "local.h"
}

#include <errno.h>

#include "vari.hpp"
//#define MAXV 1000
//#define LF_WORK 102400

//static char *db = NULL;
//static INT lfwptr, lf_work;
//vari root;



#include <map>
#include <string>

using namespace std;


typedef map<string, vari*> VarTable;
VarTable var_table;

void cleardb()
{
    for (VarTable::iterator i = var_table.begin(); i != var_table.end(); ++i)
    {
        if (i->second && i->second->stat != STSYSPEC)
        {
            free(i->second->dpr);
        }
        free(i->second);
    }
    var_table.clear();
}

void initdb() /* initialize locfit's work space */
{ 
    cleardb();
    
//    char *z = NULL;
//    z = getenv("LFWORK");
//    if (z==NULL) lf_work = LF_WORK;
//    else sscanf(z,"%d",&lf_work);
//    lf_work <<= 10;
//    if (db != NULL)
//    {
//        free(db);
//        lfwptr = 0;
//    }
//    db = (char *)calloc(lf_work, 1);
//    if (db == NULL)
//    {
//        fprintf(stderr, "Error: Locfit working space could not be allocated!\n");
//        fprintf(stderr, "Error code %d\n", errno);
//    }
//    
//    root.stat = STSYSTEM;
//    root.mode = VVARI;
//    root.dpr = (double *)db;
//    lfwptr = root.bytes = MAXV*sizeof(vari);
//    root.n = 0;
}


INT vbytes(int n, int mode)
{ 
    switch(mode)
    { 
        case VDOUBLE: return(n*sizeof(double));
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

/* vdptr with NULL check */
double *vdptr(vari* v)

{
    if (v==NULL) 
        return(NULL);
    return(v->dpr);
}

/* return the i'th data item. Cyclic. */
double vitem(vari* v, int i)
{ 
    int index;
    if ((v==NULL) || (vlength(v)==0)) 
        return(0.0);
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

void vassn(vari* v, int i, double x)
{ 
    vdptr(v)[i] = x;
}

vari *growvar(vari* vold, int n)
{ 
    fprintf(stderr, "Error: attempting to grow variable not supported\n");
    return NULL;
//    vari *vnew;
//    int reqd_bytes;
//    
//    if (vold==NULL)
//    { 
//        ERROR(("growvar: NULL old"));
//        return(NULL);
//    }
//    
//    reqd_bytes = vbytes(n, vmode(vold));
//    if (reqd_bytes <= vold->bytes) 
//        return(vold);
//    
//    vnew = createvar("_grow",vold->stat,n,vmode(vold));
//    memcpy(vdptr(vnew),vdptr(vold),vbytes(vlength(vold),vmode(vold)));
//    setvarname(vnew,vold->name);
//    vlength(vnew) = vlength(vold);
//    deletevar(vold);
//    return(vnew);
}

void *viptr(vari* v, int i) /* return pointer to ith data item, take account of mode */
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

void setvarname(vari* v, varname name)
{ 
    if (strcmp(v->name,name)==0) 
        return;
    deletename(name);
    strcpy(v->name,name);
}

/*
 findvar finds the variable name.
 err=0, keep quiet if not found; 1 produce error message.
 *n returns length of variable (if initially>0, checks length)
 */

vari *findvar(varname name, int err, int* n)
{ 
    INT status;
    vari *v;
    
    if (strcmp(name,"_NuLl")==0) return(NULL);
    
    VarTable::iterator i = var_table.find(name);
    
    if (i != var_table.end())
    {
        v = i->second;
        if (v == NULL)
        {
            fprintf(stderr, "Found variable named %s, but data is NULL\n", name);
            return NULL;
        }
        status = v->stat;
        if (status != STHIDDEN && status != STEMPTY)
        {
            if (n == NULL)
                return v;
            if (*n==-1) 
                *n = vlength(v);
            if (*n==0 || *n==vlength(v)) 
                return(v);
            if (err) 
                ERROR(("Variable %s has wrong length",name));
        }
    }

    if (err) 
        ERROR(("Variable %s not found",name));
    return NULL;
}

void deletevar(vari* v) /* delete variable, or top variable if NULL */
{ 
    if (v == NULL)
    {
        fprintf(stderr, "Error: attempting to clear entire table through NULL delete\n");
        return;
    }
    
    VarTable::iterator i = var_table.find(v->name);
    if (i != var_table.end())
    {   
        if (i->second && i->second->stat != STSYSPEC)
        {
            free(i->second->dpr);
        }
        free(i->second);
        var_table.erase(i);
    }
}

void deleteifhidden(vari* v)
{ 
    if (v==NULL) 
        return;
    if (v->stat == STHIDDEN) deletevar(v);
}

void deletename(varname name) /* delete variable name, or top variable if NULL */
{ 
    vari *v;
    v = findvar(name,0,NULL);
    if (v!=NULL) 
        deletevar(v);
}

vari *createvar(varname name, int status, int n, int mode)
{
    int bytes;
    vari *v;
    
    /*
     compute the length of the variable in bytes. some systems
     mess up is this is not a multiple of 8.
     */
    bytes = vbytes(n,mode);
    while ( (bytes & 8) > 0 ) bytes++;
    
    if (lf_error) 
        return(NULL);
    
    // Don't delete the hidden vars
    if (status==STSYSTEM || status==STREGULAR || status==STPLOTVAR)
        deletename(name);
    
    v = findvar(name,0,NULL);
    if (v != NULL)
    {
        fprintf(stderr, "Error: attempting to re-initialize still-live variable %s\n", name);
    }
    
    pair<map<string, vari*>::iterator, bool> inserted;
    string str_name = name;
    pair<string, vari*> p;
    p.first = str_name;
    p.second = (vari*)calloc(1, sizeof(vari));
    inserted = var_table.insert(p);
    
    v = inserted.first->second;
    
    strcpy(v->name,name);
    vlength(v) = n;
    v->stat = status;
    v->bytes = bytes;
    v->mode = mode;
    if (status!=STSYSPEC)
    { 
        v->dpr = (double*)calloc(bytes, 1);
    }

    return(v);
}
