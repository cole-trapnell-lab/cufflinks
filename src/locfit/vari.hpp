//
//  vari.hpp
//  cufflinks
//
//  Created by Cole Trapnell on 3/22/11.
//  Copyright 2011 Cole Trapnell. All rights reserved.
//

#ifdef __cplusplus
extern "C"
{
#endif
    void cleardb();
    void initdb();
    INT vbytes(int n, int mode);
    void setvarname(vari* v, varname name);
    double *vdptr(vari* v);
    double vitem(vari* v, int i);
    void vassn(vari* v, int i, double x);
    vari *findvar(varname name, int err, int* n);
    vari *growvar(vari* vold, int n);
    void *viptr(vari* v, int i);
    void setvarname(vari* v, varname name);
    vari *findvar(varname name, int err, int* n);
    
    void deletevar(vari* v); /* delete variable, or top variable if NULL */
    void deleteifhidden(vari* v);
    vari *createvar(varname name, int status, int n, int mode);
    void deletename(varname name);
    
#ifdef __cplusplus
}
#endif
