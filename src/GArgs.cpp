#include <stdlib.h>
#include <string.h>
#include "GArgs.h"
#include <ctype.h>

#define TRACE 1

#include "GBase.h"

//GArgs::is_opt="1"; //just to have a non-NULL value for switch testing


GArgs::GArgs(int argc, char* const argv[], const char* format, bool nodigitopts) {
   /* format is:
       <letter>[:]    for e.g. p:hT    <-  -p testing -ptesting -h -T
       <string>=      for e.g. PID=S=  <-  PID=50 S=3.5
   This means that the = options, if present, must NEVER be given after 
   dashed switches (non-value) directly
   */
  
//parse format string first:
const char* fstr=format;
fmtcount=0;
count=0;
nonOptCount=0;
nonOptPos=0;
optPos=0;
errarg=0;
args=NULL;
fmt=NULL;
int fmtlen=strlen(format);
while (fstr-format < fmtlen ) {
  int l=strcspn(fstr, ":=");
  if (fstr[l]=='\0') { //end of string reached
      //all previous chars are just switches:
       GREALLOC(fmt, (fmtcount+l)*sizeof(fmtdef));
       //store each switches
       for (int i=0; i<l;i++) { 
         GCALLOC(fmt[fmtcount+i].opt, 2); //one char length
         fmt[fmtcount+i].opt[0]=fstr[i];
         fmt[fmtcount+i].type = 0;
         }
       fmtcount+=l;
       break;
     }
   else {
     if (fstr[l]==':') {
         //fstr[l-1] is an argument, but all the previous are just switches
         GREALLOC(fmt, (fmtcount+l)*sizeof(fmtdef));
         //store each switches AND the option
         for (int i=0; i<l;i++) { 
           GCALLOC(fmt[fmtcount+i].opt, 2); //one char length
           fmt[fmtcount+i].opt[0]=fstr[i];
           fmt[fmtcount+i].type = (i==l-1)?1:0;
           }
         fmtcount+=l;
         }
      else { // fstr[l]=='=' case!
         //all these chars are one = style argument
         GREALLOC(fmt, (fmtcount+1)*sizeof(fmtdef));
         GMALLOC(fmt[fmtcount].opt, l+1);
         strncpy(fmt[fmtcount].opt, fstr, l);
         fmt[fmtcount].opt[l]='\0';
         fmt[fmtcount].type=2;
         fmtcount++;
         }
     fstr+=l+1;
     }
  }
//---- that was the parsing of the format string
//now parse the arguments based on given format specification
int p=1; //skip program name
int f=0;
//GMessage("argc=%d\n", argc);
while (p<argc) {
 if (argv[p][0]=='-') { //dashed argument ?
   int cpos=1;
   char c=argv[p][cpos];
   if (c==0 || (nodigitopts && isdigit(c))) { 
      //special case: plain argument '-' or negative number
      GREALLOC(args, (count+1)*sizeof(argdata));
      args[count].opt=NULL;
      if (c==0) {
        GCALLOC(args[count].value, 2);
        args[count].value[0]='-';
        }
       else {
        args[count].value=Gstrdup(argv[p]);
        }
      count++;
      nonOptCount++;
      }
    else { //dashed argument or switch
     COLLAPSED:
      if ((f=validOpt(c))>=0) {
        if (fmt[f].type==0) {//switch type
          GREALLOC(args, (count+1)*sizeof(argdata));
          GCALLOC(args[count].opt, 2);
          args[count].opt[0]=c;
          GCALLOC(args[count].value, 1);
          count++;
          // only switches can be grouped with some other switches or options
          if (argv[p][cpos+1]!='\0') {
             cpos++;
             c=argv[p][cpos];
             goto COLLAPSED;
             }
          }
         else 
           if (fmt[f].type==1) { //dash argument
            GREALLOC(args, (count+1)*sizeof(argdata));
            GCALLOC(args[count].opt, 2);
            args[count].opt[0]=c;
            if (argv[p][cpos+1]=='\0') {
              if (p+1<argc) { //value is the whole next argument
                 p++;
                 GMALLOC(args[count].value, strlen(argv[p])+1);
                 strcpy(args[count].value, argv[p]);
                 }
               else {
                 errarg=p;
                 return;
                 }
              }
             else { //value immediately follows the dash-option
                GMALLOC(args[count].value, strlen(argv[p])-cpos);
                strcpy(args[count].value, (argv[p]+cpos+1));
                //GMessage("args[%d].value = '%s'",count, args[count].value);
              }
            count++;
            }
           else {//inconsistent type
             errarg=p;
             return;
             } 
        } //was validOpt
       else { //option not found in format definition!
         errarg=p;
         return;
         }
      }
   }
 else {//not a dashed argument
   char* e=strchr(argv[p],'=');
   if (e!=NULL && strchr(format,'=')!=NULL && e!=argv[p] && *(e-1)!='\\') { 
     //this must be an '=' option
     //yet the '=' char can be preceded by a '\' in order to not be parsed
     //as a = option
     char part[128];
     strncpy(part, argv[p], e-argv[p]);
     part[e-argv[p]]='\0';
     if ((f=validOpt(part))>=0 && fmt[f].type==2) {
          GREALLOC(args, (count+1)*sizeof(argdata));
          args[count].opt=Gstrdup(part);
          if (strlen(argv[p])-strlen(part)>0) {
            GMALLOC(args[count].value, strlen(argv[p])-strlen(part)+1);
            strcpy(args[count].value, e+1);
            }
           else {
            args[count].value=NULL;
            } 
          count++;
          }
        else { //error - format does not match this '=' argument
         errarg=p;
         return;        
         }
      }
    else { //it seems it's just a plain argument, like a filename, etc.
     GREALLOC(args, (count+1)*sizeof(argdata));
     args[count].opt=NULL; //it's not an option
     args[count].value=Gstrdup(argv[p]);
     count++;
     nonOptCount++;
     }
   }
 p++;//check next arg string
 }
}

GArgs::~GArgs() {
 int i;
 for (i=0; i<fmtcount; i++)
    GFREE(fmt[i].opt);
 GFREE(fmt);
 for (i=0; i<count; i++) {
  GFREE(args[i].opt);
  GFREE(args[i].value);
  }
 GFREE(args);  
}

int GArgs::validOpt(char o) {
 for (int i=0; i<fmtcount; i++) 
  if (fmt[i].opt[0]==o && fmt[i].opt[1]=='\0') return i;
 return -1; 
}

int GArgs::validOpt(char* o) {
 for (int i=0; i<fmtcount; i++) 
  if (strcmp(fmt[i].opt, o)==0) return i;
 return -1;
}

int GArgs::isError() { // returns the offending argv position or 0 if no error
 return errarg;
 }

char* GArgs::getOpt(const char* o) { /* retrieve the value for option o
                   returns 
                       NULL    if option not given at all
                     !=NULL    if boolean option was given
                     opt.value if value option was given
                     */
 for (int i=0; i<count; i++) 
  if (args[i].opt!=NULL && strcmp(args[i].opt, o)==0) 
           return args[i].value;
 return NULL;                    
}

char* GArgs::getOpt(const char o) {
 for (int i=0; i<count; i++) 
  if (args[i].opt!=NULL && args[i].opt[0]==o && args[i].opt[1]=='\0') 
      return args[i].value;
 return NULL;
}

int GArgs::startNonOpt(){ //reset iteration through non-dashed arguments
   //returns the number of non-dashed arguments
nonOptPos=0;
return nonOptCount;   
}
   
   
char* GArgs::nextNonOpt() { //get the next non-dashed argument
               //or NULL if no more 
for (int i=nonOptPos;i<count;i++)
 if (args[i].opt==NULL) {
      nonOptPos=i+1;
      return args[i].value;
      }
return NULL;
}

int GArgs::startOpt(){ //reset iteration through non-dashed arguments
   //returns the number of non-dashed arguments
optPos=0;
return count-nonOptCount;
}
   
   
char* GArgs::nextOpt() { //get the next non-dashed argument
               //or NULL if no more 
for (int i=optPos;i<count;i++)
 if (args[i].opt!=NULL) {
      optPos=i+1;
      return args[i].opt;
      }
return NULL;
}
