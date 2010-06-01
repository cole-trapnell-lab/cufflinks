/*
**
GArgs is a quick'n'dirty object oriented replacement for the standard 
   getopts library call avialable on many unix platforms. 
   it accepts dash style options and = style options
   -<letter>[ ][<value>] 
   <string>=<value>
*/
   
#ifndef G_ARGS_DEFINED
#define G_ARGS_DEFINED



class GArgs {
   //structure for parsing arguments format definition
   struct fmtdef {
    int type; // 0=dashed switch, 1=dashed value, 2='=' value
    char* opt; //switch/opt char/string
    };
   int fmtcount;
   fmtdef* fmt; //this will store format definition after parsing it
   
   struct argdata {
     char*  opt; // this is NULL for non-dashed arguments
                    //one character for dashed style arguments
                    //one string for     = style arguments
     char* value; // is NULL for switches (dashed flags)
     };
   argdata* args; //arguments table after parsing it
   int count; //total count of elements in 'args' array
   int nonOptCount; //count of non-dashed, non= arguments
   int nonOptPos; //current position for nonOpt arguments iterator
   int optPos; //current position for options iterator
   int errarg; //argv error position after parsing
   static const char* is_opt; // = non NULL just for getOpt easy testing
   int validOpt(char o);  //parsing helper function
   int validOpt(char* o);
 public:
 
   GArgs(int argc, char* const argv[], const char* format, bool nodigitopts=false);
   /* format is:
       <letter>[:]    for e.g. p:hT    <=  -p testing -ptesting -h -T
       <string>=      for e.g. PID=S=  <=  PID=50 S=3.5
   This means that the = options, if present, must NEVER be given after 
   dashed (non-value) switches directly
   */
   ~GArgs();
   int isError(); // returns the offending argv position or 0 if no error
   int getCount() { return count; }
   int getFmtCount() { return fmtcount; }
   int getNonOptCount() { return nonOptCount; }
   char* getOpt(const char* o); /* retrieve the value for option o
                   returns 
                       NULL    if option not given at all
                     !=NULL    if boolean option was given
                     opt's value if value option was given
                     */
   char* getOpt(const char o);
   int startOpt(); //init iteration through option arguments
   char* nextOpt(); //get next option argument
   int startNonOpt(); //init iteration through non-option arguments
             //returns the number of non-option arguments
   char* nextNonOpt(); //get the next non-option argument
   
};

#endif
