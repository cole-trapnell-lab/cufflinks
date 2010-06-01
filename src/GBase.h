#ifndef G_BASE_DEFINED
#define G_BASE_DEFINED

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#if defined __WIN32__ || defined _WIN32
  #include <windows.h>
#endif

#ifdef DEBUG
#undef NDEBUG
#endif

typedef unsigned int uint32;
typedef int int32;
typedef unsigned char uchar;
typedef unsigned char byte;

// If long is natively 64 bit, use the regular fseek and ftell
#ifdef _NATIVE_64
 #define ftello ftell
 #define fseeko fseek
#endif

#ifndef MAXUINT
#define MAXUINT ((unsigned int)-1)
#endif

#if defined(_NATIVE_64) || defined(_LP64) || defined(__LP64__)
 typedef long int64;
 typedef unsigned long uint64;
#else
 //assume 32bit environment with long long for int64 stuff
 typedef long long int64;
 typedef unsigned long long uint64;
#endif

/****************************************************************************/

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

/****************************************************************************/
#define ERR_ALLOC "Error allocating memory.\n"
#if defined (__WIN32__) || defined (WIN32)
  #define CHPATHSEP '\\'
  #include <io.h>
  #define ftello ftell
  #define fseeko fseek
 #else
  #define CHPATHSEP '/'
  #include <unistd.h>
#endif

//-------------------

// Debug helpers
#ifndef NDEBUG
 #define GASSERT(exp) ((exp)?((void)0):(void)GAssert(#exp,__FILE__,__LINE__))
 #ifdef TRACE
  #define GTRACE(exp)  (GMessage exp)
 #else
  #define GTRACE(exp)  ((void)0)
 #endif
#else
 #define GASSERT(exp) ((void)0)
 #define GTRACE(exp)  ((void)0)
#endif

#define GERROR(exp) (GError exp)
/**********************************  Macros  ***********************************/
// Abolute value
#define GABS(val) (((val)>=0)?(val):-(val))

// Min and Max
#define GMAX(a,b) (((a)>(b))?(a):(b))
#define GMIN(a,b) (((a)>(b))?(b):(a))

// Min of three
#define GMIN3(x,y,z) ((x)<(y)?GMIN(x,z):GMIN(y,z))

// Max of three
#define GMAX3(x,y,z) ((x)>(y)?GMAX(x,z):GMAX(y,z))

// Return minimum and maximum of a, b
#define GMINMAX(lo,hi,a,b) ((a)<(b)?((lo)=(a),(hi)=(b)):((lo)=(b),(hi)=(a)))

// Clamp value x to range [lo..hi]
#define GCLAMP(lo,x,hi) ((x)<(lo)?(lo):((x)>(hi)?(hi):(x)))

typedef void* pointer;
typedef unsigned int uint;

typedef int GCompareProc(const pointer item1, const pointer item2);
typedef void GFreeProc(pointer item); //usually just delete,
      //but may also support structures with embedded dynamic members

#define GMALLOC(ptr,size)  if (!GMalloc((pointer*)(&ptr),size)) \
                                     GError(ERR_ALLOC)
#define GCALLOC(ptr,size)  if (!GCalloc((pointer*)(&ptr),size)) \
                                     GError(ERR_ALLOC)
#define GREALLOC(ptr,size) if (!GRealloc((pointer*)(&ptr),size)) \
                                     GError(ERR_ALLOC)
#define GFREE(ptr)       GFree((pointer*)(&ptr))

inline char* min(char *arg1, char *arg2) {
    return (strcmp(arg1, arg2) < 0)? arg1 : arg2;
}

inline int iround(double x) {
   return (int)floor(x + 0.5);
}


/****************************************************************************/

inline char* max(char *arg1, char *arg2) {
    return (strcmp(arg2, arg1) < 0)? arg1 : arg2;
}

inline int Gintcmp(int a, int b) {
 //return (a>b)? 1 : ((a==b)?0:-1);
  return a-b;
}

int Gstrcmp(char* a, char* b);
//same as strcmp but doesn't crash on NULL pointers

int Gstricmp(const char* a, const char* b);

inline void swap(int &arg1, int &arg2){
 arg1 ^= arg2 ^= arg1 ^= arg2;
 }

inline void swap(char* &arg1, char* &arg2){
 register char* swp=arg1;
 arg1=arg2; arg2=swp;
 }

inline void swap(unsigned int &arg1, unsigned int &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(short &arg1, short &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(unsigned short &arg1, unsigned short &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(long &arg1, long &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(unsigned long &arg1, unsigned long &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(char &arg1, char &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(unsigned char &arg1, unsigned char &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(bool &arg1, bool &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }


/**************** Memory management ***************************/

bool GMalloc(pointer* ptr, unsigned long size); // Allocate memory
bool GCalloc(pointer* ptr, unsigned long size); // Allocate and initialize memory
bool GRealloc(pointer* ptr,unsigned long size); // Resize memory
void GFree(pointer* ptr); // Free memory, resets ptr to NULL

/********************* debug functions *********************/

void GError(const char* format,...); // Error routine (aborts program)
void GMessage(const char* format,...);// Log message to stderr
// Assert failed routine:- usually not called directly but through GASSERT
void GAssert(const char* expression, const char* filename, unsigned int lineno);


// ****************** string manipulation *************************
char *Gstrdup(const char* str);
//duplicate a string by allocating a copy for it and returning it
char* Gstrdup(const char* sfrom, const char* sto);
//same as GStrdup, but with an early termination (e.g. on delimiter)

char* Gsubstr(const char* str, char* from, char* to=NULL);
//extracts a substring, allocating it, including boundaries (from/to)

int strsplit(char* str, char** fields, int maxfields, const char* delim);
int strsplit(char* str, char** fields, int maxfields, const char delim);
int strsplit(char* str, char** fields, int maxfields); //splits by tab or space

char* replaceStr(char* &str, char* newvalue);

//conversion: to Lower/Upper case
// creating a new string:
char* upCase(const char* str);
char* loCase(const char* str);
// changing string in place:
char* strlower(char * str);
char* strupper(char * str);

//strstr but for memory zones: scans a memory region
//for a substring:
void* Gmemscan(void *mem, unsigned int len,
                  void *part, unsigned int partlen);

// test if a char is in a string:
bool chrInStr(char c, char* str);

char* rstrchr(char* str, char ch);
/* returns a pointer to the rightmost
  occurence of ch in str - like rindex for platforms missing it*/

char* strchrs(char* s, const char* chrs);
//strchr but with a set of chars instead of only one

char* rstrfind(char* str, const char *substr); /* like rindex() but for strings
or like the right side version of strstr()
*/
//reverse character string or
char* reverseChars(char* str, int slen=0);

char* rstrstr(char* rstart, char *lend, char* substr);
/*the reversed, rightside equivalent of strstr: starts searching
 from right end (rstart), going back to left end (lend) and returns
 a pointer to the last (right) matching character in str */

char* strifind(char* str,  const char* substr);
// the case insensitive version of strstr -- finding a string within a strin


//Determines if a string begins with a given prefix
//(returns false when any of the params is NULL,
// but true when prefix is '' (empty string)!)
bool startsWith(char* s, const char* prefix);

// ELF hash function for strings
int strhash(const char* str);



//---- generic base GSeg : genomic segment (interval) --
// coordinates are considered 1-based (so 0 is invalid)
class GSeg {
 public:
  uint start; //start<end always!
  uint end;
  GSeg(uint s=0,uint e=0) {
    if (s>e) { start=e;end=s; }
        else { start=s;end=e; }
    }
  //check for overlap with other segment
  uint len() { return end-start+1; }
  bool overlap(GSeg* d) {
     return start<d->start ? (d->start<=end) : (start<=d->end);
     }

  bool overlap(GSeg& d) {
     return start<d.start ? (d.start<=end) : (start<=d.end);
     }

  bool overlap(GSeg& d, int fuzz) {
     return start<d.start ? (d.start<=end+fuzz) : (start<=d.end+fuzz);
     }

  bool overlap(uint s, uint e) {
    if (s>e) { swap(s,e); }
     return start<s ? (s<=end) : (start<=e);
     }

  //return the length of overlap between two segments
  int overlapLen(GSeg* r) {
     if (start<r->start) {
        if (r->start>end) return 0;
        return (r->end>end) ? end-r->start+1 : r->end-r->start+1;
        }
       else { //r->start<=start
        if (start>r->end) return 0;
        return (r->end<end)? r->end-start+1 : end-start+1;
        }
     }
  int overlapLen(uint rstart, uint rend) {
     if (rstart>rend) { swap(rstart,rend); }
     if (start<rstart) {
        if (rstart>end) return 0;
        return (rend>end) ? end-rstart+1 : rend-rstart+1;
        }
       else { //rstart<=start
        if (start>rend) return 0;
        return (rend<end)? rend-start+1 : end-start+1;
        }
     }

  //fuzzy coordinate matching:
  bool coordMatch(GSeg* s, uint fuzz=0) {
    if (fuzz==0) return (start==s->start && end==s->end);
    uint sd = (start>s->start) ? start-s->start : s->start-start;
    uint ed = (end>s->end) ? end-s->end : s->end-end;
    return (sd<=fuzz && ed<=fuzz);
    }
  //comparison operators required for sorting
  bool operator==(GSeg& d){
      return (start==d.start && end==d.end);
      }
  bool operator>(GSeg& d){
     return (start==d.start)?(end>d.end):(start>d.start);
     }
  bool operator<(GSeg& d){
     return (start==d.start)?(end<d.end):(start<d.start);
     }
};



//--------------------------------------------------------
// ************** simple line reading class for text files

//GLineReader -- text line reading/buffering class
class GLineReader {
   int len;
   int allocated;
   char* buf;
   bool isEOF;
   FILE* file;
   off_t filepos; //current position
   bool pushed; //pushed back
   int lcount; //line counter (read lines)
 public:
   char* chars() { return buf; }
   char* line() { return buf; }
   int readcount() { return lcount; } //number of lines read
   int length() { return len; }
   int size() { return len; } //same as size();
   bool isEof() {return isEOF; }
   bool eof() { return isEOF; }
   off_t getfpos() { return filepos; }
   off_t getFpos() { return filepos; }
   char* nextLine() { return getLine(); }
   char* getLine() { if (pushed) { pushed=false; return buf; }
                            else return getLine(file);  }
   char* getLine(FILE* stream) {
                 if (pushed) { pushed=false; return buf; }
                          else return getLine(stream, filepos); }
   char* getLine(FILE* stream, off_t& f_pos); //read a line from a stream and update
                           // the given file position
   void pushBack() { if (lcount>0) pushed=true; } // "undo" the last getLine request
            // so the next call will in fact return the same line
   GLineReader(FILE* stream=NULL, off_t fpos=0) {
     len=0;
     isEOF=false;
     allocated=1024;
     GMALLOC(buf,allocated);
     lcount=0;
     buf[0]=0;
     file=stream;
     filepos=fpos;
     pushed=false;
     }
   ~GLineReader() {
     GFREE(buf);
     }
};


/* extended fgets() -  to read one full line from a file and
  update the file position correctly !
  buf will be reallocated as necessary, to fit the whole line
  */
char* fgetline(char* & buf, int& buflen, FILE* stream, off_t* f_pos=NULL, int* linelen=NULL);

/*********************** File management functions *********************/

// removes the directory part from a full-path file name
// this is a destructive operation for the given string!
void delFileName(char* filepath);

// returns a pointer to the file name part in a full-path filename
char* getFileName(char* filepath);

int fileExists(const char* fname);
//returns 0 if file entry doesn't exist
//        1 if it's a directory
//        2 if it's a regular file
//        3 otherwise (?)

off_t fileSize(const char* fpath);

//parses the next number found in a string at the current position
//until a non-digit (and not a '.', 'e','E','-','+') is encountered;
//updates the char* pointer to be after the last digit parsed
bool parseNumber(char* &p, double& v);
bool parseDouble(char* &p, double& v); //just an alias for parseNumber

bool parseInt(char* &p, int& i);
bool parseUInt(char* &p, uint& i);
bool parseHex(char* &p,  uint& i);

#endif /* G_BASE_DEFINED */
