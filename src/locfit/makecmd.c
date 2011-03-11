/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   The makecmd() function converts a command line string
 *   into a locfit command variable. If the line has no
 *   commands (for example, a blank line or a comment)
 *   it returns NULL; otherwise, it returns the pointer to
 *   the command variable.
 *
 *   The command line is split into arguments, with arguments
 *   separated by spaces. Exception: A space in a quoted
 *   "str ing" is not split into separate fields.
 *
 *   getcmd() returns a pointer to the next line for processing.
 *     If no lines are ready for processing, it returns NULL.
 *
 *   del_lines()  frees the work space used by processed command lines.
 *
 *   set_forvar(), inc_forvar(), dec_forvar() are used in the
 *     control of for loops.
 */

#include "local.h"

#define isterminator(c) (((c)=='\0') | ((c)=='\n') | ((c)=='#'))
static int clcount = 0;
static int proc_to = 0;
static int del_to = 0;
extern vari *curstr;
extern char filename[];

typedef struct {
  vari *v;
  char *name;
  int line_no;
  int index; } forvar;
static forvar fv[10];
static int for_level = 0, proc_level = 0;

int countfields(z)
char *z;
{ int n, instr;

  n = 0;
  instr = 0;

  while (1)
  { while (*z==' ') z++;
    if (isterminator(*z)) return(n);

    n++;

    while ((instr) || (*z!=' '))
    { if (isterminator(*z))
      {  if (instr) ERROR(("Unterminated String"));
         return(n);
      }
      if (*z=='"') instr = !instr;
      z++;
    }
  }
}

void makefields(z, va, n)
char *z;
vari *va;
int n;
{ int i, instr;
  char *st, *eq;
  carg *curr_arg;

  instr = 0;
  for (i=0; i<n; i++)
  { while (*z==' ') z++;

    curr_arg = (carg *)viptr(va,i);
    curr_arg->val = st = z;
    curr_arg->arg = NULL;
    curr_arg->result = NULL;

    eq = NULL;
    do
    { if (*z=='"') instr = !instr;
      if ((eq==NULL) && (!instr) && (*z=='=')) eq = z;
      z++;
    } while ((instr) || (*z !=' ') && (!isterminator(*z)));
    *z = '\0';

    if (eq != NULL)
    { if (eq==st)
      { ERROR(("command line argument begins with ="));
        return;
      }
      if ((eq[1]!='=') & (strchr("<>!",eq[-1])==NULL))
      { curr_arg->arg = st;
        curr_arg->val = &eq[1];
        *eq = '\0';
      }
    } /* eq != */
    z++;

  }   /* for i */
}

/*
 *  set_forvar and inc_forvar are used to control for loops.
 *  set_forvar is called when the for cmd is built, making the
 *    variable to loop through.
 *  inc_forvar is called when the for is processed, to update
 *    the value of the variable.
 *  dec_forvar is called when the endfor is processed. This
 *    resets the proc_to line count to the beginning of the for loop.
 */
void set_forvar(v,ct)
vari *v;
int ct;
{ 
  varname vn;

  if (vlength(v)<2)
  { ERROR(("for: missing variable"));
    return;
  }

  sprintf(vn,"=forv%d",for_level);
  fv[for_level].v = varith(argval(v,1),vn,STHIDDEN);
  fv[for_level].line_no = ct;
  fv[for_level].index = 0;
  fv[for_level].name = argarg(v,1);
  for_level++;
}

void inc_forvar()
{ varname vn;
  vari *v, *va;
  double x;
  
  if (fv[proc_level].name == NULL)
  { sprintf(vn,"=fv%d",proc_level);
    v = createvar(vn,STHIDDEN,1,VDOUBLE);
  }
  else
    v = createvar(fv[proc_level].name,STREGULAR,1,VDOUBLE);

  va = fv[proc_level].v;
  x = vitem(va, fv[proc_level].index);
  vassn(v,0,x);
  fv[proc_level].index++;
  proc_level++;
}

void dec_forvar()
{ proc_level--;
  if (fv[proc_level].index < vlength(fv[proc_level].v))
    proc_to = fv[proc_level].line_no - 1;
  else
    fv[proc_level].index = 0;
}

void run(va)
vari *va;
{ FILE *runf;
  char cline[256], *z;

  if (vlength(va)<2)
  { ERROR(("run: no command file"));
    return;
  }

  if (!setfilename(argval(va,1),"","r",0))
  { ERROR(("run: cannot read file %s",argval(va,1)));
    return;
  }

  runf = fopen(filename,"r");
  while (1)
  { z = fgets(cline,256,runf);
    if (z==NULL)
    { fclose(runf);
      return;
    }
    makecmd(cline);
  }
}

void makecmd(cmdline)
char *cmdline;
{ 
  varname vn;
  vari *va, *vs;
  int n;

  n = countfields(cmdline);
  if (lf_error) return;
  if (n==0) return;
  clcount++;

  /* vs is used to store the command line string. */
  sprintf(vn,"=clstr%d",clcount);
  vs = createvar(vn,STSYSTEM,1+strlen(cmdline),VCHAR);
  sprintf((char *)vdptr(vs),cmdline);

  /* va is used to store pointers to the command line fields. */
  sprintf(vn,"=cline%d",clcount);
  va = createvar(vn,STSYSTEM,(INT)n,VARGL);
  makefields((char *)vdptr(vs), va, n);

  if (argvalis(va,0,"for")) set_forvar(va,clcount);
  if (argvalis(va,0,"endfor")) for_level--;

/* we want to read in run files here, not when commands are executed.
 * otherwise, we would have problems with run commands in a for loop.
 */
  if (argvalis(va,0,"run")) run(va);

  return;
}

void del_lines()
{ int i;
  varname vn;
  for (i=proc_to; i>del_to; i--)
  { sprintf(vn,"=cline%d",i);
    deletename(vn);
    sprintf(vn,"=clstr%d",i);
    deletename(vn);
  }
  del_to = proc_to;
}

vari *getcmd()
{
  varname vn;
  vari *v;

  if (for_level > 0) return(NULL);

  if (proc_to < clcount)
  { 
    sprintf(vn,"=cline%d",++proc_to);
    v = findvar(vn,1,NULL);
    if (v==NULL) return(v);

/* this nonsense is req'd by setplot and setdef.
 * get rid of it, I hope.
 */
sprintf(vn,"=clstr%d",proc_to);
curstr = findvar(vn,1,NULL);

    return(v);
  }

  return(NULL);
}
