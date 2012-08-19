/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *   The example() function is used to tnterpreting examples
 *   in the locfit.ceg file.
 */

#include "local.h"

#ifdef CVERSION

static FILE *help;
extern char *lfhome;

void example(v)
vari *v;
{ int i, j, run, ti;
  char *z, elin[10], helpfile[100], line[200];

  run = 0;
  i = getarg(v,"ex",1);
  ti = (i==0);
  if (!ti)
  { z = strchr(argval(v,i),'.');
    if (z==NULL)
    { ERROR(("example: invalid number %s",argval(v,i)));
      return;
    }
    j = getarg(v,"run",1);
    if (j>0) run = getlogic(v,j);
  }

  if (lfhome!=NULL)
    sprintf(helpfile,"%s/locfit.ceg",lfhome);
  else
    sprintf(helpfile,"locfit.ceg");

  help = fopen(helpfile,"r");
  if (help==NULL)
  { ERROR(("Can't find locfit.ceg file -- is LFHOME set?"));
    return;
  }

  do
  { z = fgets(line,190,help);
    if (z==NULL) /* reached end of file */
    { if (!ti) ERROR(("No example %s in help file",argval(v,i)));
      fclose(help);
      return;
    }
    if (line[0]=='e')
    { sscanf(&line[2],"%s",elin);
      if (ti) printf("Example %s. ",elin);
    }
    else elin[0] = 0;
    if ((ti) && (line[0]=='t'))
      printf("%s",&line[2]);
  }  while ((ti) || (strcmp(elin,argval(v,i))!=0));

  while(1)
  { z = fgets(line,190,help);
    switch(z[0])
    { case 'f': /* end example */
        fclose(help);
        printf("\n");
        return;
      case 't': /* title */
        printf("\nExample %s. %s\n",argval(v,i),&line[2]);
        break;
      case 'c': /* code */
        printf("  %s",&line[2]);
      case 'h': /* hidden code, usually sleep */
        if (run) makecmd(&line[2]);
        break;
      case 'd': /* discussion */
        printf("%s",&line[2]);
        break;
      case 'n': /* no code */
        printf("There is no code for this example.\n");
        break;
    }
  }
}

#endif
