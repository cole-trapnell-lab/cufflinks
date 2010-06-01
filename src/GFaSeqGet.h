#ifndef GFASEQGET_H
#define GFASEQGET_H

#include "GBase.h"
#include "GList.hh"


#define MAX_FASUBSEQ 0x10000000
//max 256MB sequence data held in memory at a time

class GSubSeq {
 public:
  uint sqstart; //1-based coord of subseq start on sequence
  uint sqlen;   //length of subseq loaded
  char* sq; //actual subsequence data will be stored here
                // (with end-of-line characters removed)

  /*char* xseq; //the exposed pointer to the last requested subsequence start
  off_t xstart; //the coordinate start for the last requested subseq
  off_t xlen; //the last requested subseq len*/
  GSubSeq() {
     sqstart=0;
     sqlen=0;
     sq=NULL;
     /* xseq=NULL;
     xstart=0;
     xlen=0;*/
     }
  ~GSubSeq() {
     GFREE(sq);
     }
  // genomic, 1-based coordinates:
  void setup(uint sstart, int slen, int sovl=0, int qfrom=0, int qto=0);
    //check for overlap with previous window and realloc/extend appropriately
    //returns offset from seq that corresponds to sstart
    // the window will keep extending until MAX_FASUBSEQ is reached
};

class GFaSeqGet {
  char* fname;
  FILE* fh;
  //raw offset in the file where the sequence actually starts:
  off_t fseqstart;
  int linelen; //length of each sequence line (assumed fixed)
  char lendlen; //length of end-of-line characters between lines
                         //(assumed fixed)
  char lendch; //end-of-line signal character (can only be '\n' or '\r')
  GSubSeq* lastsub;
  void initialParse(off_t fofs=0, bool checkall=true);
  const char* loadsubseq(uint cstart, int& clen);
  void finit(const char* fn, off_t fofs, bool validate);
 public:
  GFaSeqGet() {
    fseqstart=0;
    linelen=0;
    lendch='\0';
    fname=NULL;
    lastsub=NULL;
    }
  GFaSeqGet(const char* fn, off_t fofs, bool validate=false) { 
     finit(fn,fofs,validate); 
     }
  GFaSeqGet(const char* fn, bool validate=false) {
     finit(fn,0,validate);
     }
  /*
  GFaSeqGet(bool readAll, const char* fn, off_t fofs=0);
  GFaSeqGet(bool readAll, FILE* f, off_t fofs=0);
  */
  GFaSeqGet(FILE* f, off_t fofs=0, bool validate=false);
  ~GFaSeqGet() {
    if (fname!=NULL) {
       GFREE(fname);
       fclose(fh);
       }
    delete lastsub;
    }
  const char* subseq(uint cstart, int& clen);
  const char* getRange(uint cstart, uint cend) {
      if (cstart>cend) { swap(cstart, cend); }
      int clen=cend-cstart+1;
      //int rdlen=clen;
      return subseq(cstart, clen);
      }
  //caller is responsible for deallocating copyRange() return string
  char* copyRange(uint cstart, uint cend, bool revCmpl=false, bool upCase=false);

  void loadall() {
    int clen=MAX_FASUBSEQ;
    subseq(1, clen);
    }
  void load(uint cstart, uint cend) {
     //cache as much as possible
      int clen=cend-cstart+1;
      subseq(cstart, clen);
     }
  int getsublen() { return lastsub!=NULL ? lastsub->sqlen : 0 ; }
  off_t getseqofs() { return fseqstart; }
  int getlinelen() { return linelen; }
  int getlendlen() { return lendlen; }
  //reads a subsequence starting at genomic coordinate cstart (1-based)
 };


#endif
