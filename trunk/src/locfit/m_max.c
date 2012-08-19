/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *  Routines for maximization of a one dimensional function f()
 *    over an interval [xlo,xhi]. In all cases. the flag argument
 *    controls the return:
 *      flag='x', the maximizer xmax is returned.
 *                otherwise, maximum f(xmax) is returned.
 *
 *  max_grid(f,xlo,xhi,n,flag)
 *    grid maximization of f() over [xlo,xhi] with n intervals.
 *
 *  max_golden(f,xlo,xhi,n,tol,err,flag)
 *    golden section maximization.
 *    If n>2, an initial grid search is performed with n intervals
 *      (this helps deal with local maxima).
 *    convergence criterion is |x-xmax| < tol.
 *    err is an error flag.
 *    if flag='x', return value is xmax.
 *       otherwise, return value is f(xmax).
 *
 *  max_quad(f,xlo,xhi,n,tol,err,flag)
 *    quadratic maximization.
 *
 *  max_nr()
 *    newton-raphson, handles multivariate case.
 *
 *  TODO: additional error checking, non-convergence stop.
 */

#include <math.h>
#include <string.h>
#include "mutil.h"
extern double innerprod();

#define gold_rat 0.6180339887498948482045870
#define max_val(x,y) ((flag=='x') ? x : y)

double max_grid(f,xlo,xhi,n,flag)
double (*f)(), xlo, xhi;
int n;
char flag;
{ int i, mi;
    mi = 0;
    double x, y, mx, my;
    my = 0.0;
    for (i=0; i<=n; i++)
    { x = xlo + (xhi-xlo)*i/n;
        y = f(x);
        if ((i==0) || (y>my))
        { mx = x;
            my = y;
            mi = i;
        }
    }
    if (mi==0) return(max_val(xlo,my));
    if (mi==n) return(max_val(xhi,my));
    return(max_val(mx,my));
}

double max_golden(f,xlo,xhi,n,tol,err,flag)
double (*f)(), xhi, xlo, tol;
int n, *err;
char flag;
{ double x0, x1, x2, x3, y0, y1, y2, y3;
    *err = 0;
    
    if (n>2)
    { x0 = max_grid(f,xlo,xhi,n,'x');
        if (xlo<x0) xlo = x0-1.0/n;
        if (xhi>x0) xhi = x0+1.0/n;
    }
    
    x0 = xlo; y0 = f(xlo);
    x3 = xhi; y3 = f(xhi);
    x1 = gold_rat*x0 + (1-gold_rat)*x3; y1 = f(x1);
    x2 = gold_rat*x3 + (1-gold_rat)*x1; y2 = f(x2);
    
    while (fabs(x3-x0)>tol)
    { if ((y1>=y0) && (y1>=y2))
    { x3 = x2; y3 = y2;
        x2 = x1; y2 = y1;
        x1 = gold_rat*x0 + (1-gold_rat)*x3; y1 = f(x1);
    }
    else if ((y2>=y3) && (y2>=y1))
    { x0 = x1; y0 = y1;
        x1 = x2; y1 = y2;
        x2 = gold_rat*x3 + (1-gold_rat)*x1; y2 = f(x2);
    }
    else
    { if (y3>y0) { x0 = x2; y0 = y2; }
    else { x3 = x1; y3 = y1; }
        x1 = gold_rat*x0 + (1-gold_rat)*x3; y1 = f(x1);
        x2 = gold_rat*x3 + (1-gold_rat)*x1; y2 = f(x2);
    }
    }
    if (y0>=y1) return(max_val(x0,y0));
    if (y3>=y2) return(max_val(x3,y3));
    return((y1>y2) ? max_val(x1,y1) : max_val(x2,y2));
}

double max_quad(f,xlo,xhi,n,tol,err,flag)
double (*f)(), xhi, xlo, tol;
int n, *err;
char flag;
{ double x0, x1, x2, xnew, y0, y1, y2, ynew, a, b;
    *err = 0;
    
    if (n>2)
    { x0 = max_grid(f,xlo,xhi,n,'x');
        if (xlo<x0) xlo = x0-1.0/n;
        if (xhi>x0) xhi = x0+1.0/n;
    }
    
    x0 = xlo; y0 = f(x0);
    x2 = xhi; y2 = f(x2);
    x1 = (x0+x2)/2; y1 = f(x1);
    
    while (x2-x0>tol)
    {
        /* first, check (y0,y1,y2) is a peak. If not,
         * next interval is the halve with larger of (y0,y2).
         */
        if ((y0>y1) | (y2>y1))
        { 
            if (y0>y2) { x2 = x1; y2 = y1; }
            else { x0 = x1; y0 = y1; }
            x1 = (x0+x2)/2;
            y1 = f(x1);
        }
        else /* peak */
        { a = (y1-y0)*(x2-x1) + (y1-y2)*(x1-x0);
            b = ((y1-y0)*(x2-x1)*(x2+x1) + (y1-y2)*(x1-x0)*(x1+x0))/2;
            /* quadratic maximizer is b/a. But first check if a's too
             * small, since we may be close to constant.
             */
            if ((a<=0) | (b<x0*a) | (b>x2*a))
            { /* split the larger halve */
                xnew = ((x2-x1) > (x1-x0)) ? (x1+x2)/2 : (x0+x1)/2;
            }
            else
            { xnew = b/a;
                if (10*xnew < (9*x0+x1)) xnew = (9*x0+x1)/10;
                if (10*xnew > (9*x2+x1)) xnew = (9*x2+x1)/10;
                if (fabs(xnew-x1) < 0.001*(x2-x0))
                {
                    if ((x2-x1) > (x1-x0))
                        xnew = (99*x1+x2)/100;
                    else
                        xnew = (99*x1+x0)/100;
                }
            }
            ynew = f(xnew);
            if (xnew>x1)
            { if (ynew >= y1) { x0 = x1; y0 = y1; x1 = xnew; y1 = ynew; }
            else { x2 = xnew; y2 = ynew; }
            }
            else
            { if (ynew >= y1) { x2 = x1; y2 = y1; x1 = xnew; y1 = ynew; }
            else { x0 = xnew; y0 = ynew; }
            }
        }
    }
    return(max_val(x1,y1));
}

double max_nr(F, coef, old_coef, f1, delta, J, p, maxit, tol, err)
double *coef, *old_coef, *f1, *delta, tol;
int (*F)(), p, maxit, *err;
jacobian *J;
{ double old_f, f, lambda;
    int i, j, fr;
    double nc, nd, cut;
    int rank;
    
    *err = NR_OK;
    J->p = p;
    fr = F(coef, &f, f1, J->Z); J->st = JAC_RAW;
    
    for (i=0; i<maxit; i++)
    { memcpy(old_coef,coef,p*sizeof(double));
        old_f = f;
        rank = jacob_solve(J,f1);
        memcpy(delta,f1,p*sizeof(double));
        
        
        if (rank==0) /* NR won't move! */
            delta[0] = -f/f1[0];
        
        lambda = 1.0;
        
        nc = innerprod(old_coef,old_coef,p);
        nd = innerprod(delta, delta, p);
        cut = sqrt(nc/nd);
        if (cut>1.0) cut = 1.0;
        cut *= 0.0001;
        do
        { for (j=0; j<p; j++) coef[j] = old_coef[j] + lambda*delta[j];
            fr = F(coef, &f, f1, J->Z); J->st = JAC_RAW;
            if (fr==NR_BREAK) return(f);
            
            lambda = (fr==NR_REDUCE) ? lambda/2 : lambda/10.0;
        } while ((lambda>cut) & (f <= old_f - 1.0e-3));
        
        if (f < old_f - 1.0e-3) { *err = NR_NDIV; return(f); }
        if (fr==NR_REDUCE) return(f);
        
        if (fabs(f-old_f) < tol) return(f);
        
    }
    *err = NR_NCON;
    return(f);
}
