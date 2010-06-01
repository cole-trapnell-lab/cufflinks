#ifndef GFASTAFILE_H
#define GFASTAFILE_H
#include <stdio.h>
#include "gdna.h"

#define CAPINC 64
#define SEQCAPINC 256
#define DEF_FASTA_DELIM ">"

class FastaSeq {  /* fasta record storage */
   public:
     int  id_cap; /* allocated size of the sequence name string*/
     char    *id; /* id only, up to first space */
     int namelen;   // real length of seq name
     char *descr; /* any comment on the defline, after the first space */
     int   d_cap; /* allocated size of the description */
     int descrlen; /* real length of the description */
     //-------actual sequence :
     int s_cap; /* allocated length of the sequence string */
     int   len; /* the actual string length of seq */
     char* seq; /* the sequence buffer itself */
     //----
     FastaSeq(char* cname, char* cdescr=NULL, char* cseq=NULL) {
       //ntCompTableInit();
       if (cname==NULL) {
          FastaSeq();
          return;
          }
       int l=strlen(cname);
       GMALLOC(id, l+1);strcpy(id,cname);
       id_cap=l+1;
       namelen=l;
       if (cdescr==NULL) {
           GMALLOC(descr, CAPINC);
           descr[0]='\0';
           d_cap=CAPINC;
           descrlen=0;
           }
        else {//copy given description
           l=strlen(cdescr);
           GMALLOC(descr, l+1);
           strcpy(descr,cdescr);
           d_cap=l+1;
           descrlen=l;
           }
       if (cseq==NULL) {
           GMALLOC(seq, SEQCAPINC);
           seq[0]='\0';
           len=0;
           s_cap=SEQCAPINC;
           }
         else {
           l=strlen(cseq);
           GMALLOC(seq, l+1);
           strcpy(seq,cseq);
           len=l;
           s_cap=l+1;
           }
       }
     void init(int seqalloc=0) {
       //ntCompTableInit();
       GMALLOC(id, CAPINC);
       id_cap=CAPINC;
       namelen=0;
       id[0]='\0';
       GMALLOC(descr, CAPINC);
       descr[0]='\0';
       d_cap=CAPINC;
       descrlen=0;
       if (seqalloc<=0) {
          s_cap=SEQCAPINC;
          GMALLOC(seq, SEQCAPINC);
          }
        else {
         s_cap=seqalloc;
         GMALLOC(seq, seqalloc);
          }
       seq[0]='\0';
       len=0;
       }
     FastaSeq(int seqalloc=0) {
       init(seqalloc);
       }
     void clear() {
       GFREE(id);id_cap=0;namelen=0;id=NULL;
       GFREE(descr);d_cap=0;descrlen=0;descr=NULL;
       GFREE(seq);s_cap=0;len=0;seq=NULL;
       }
     ~FastaSeq() {
       clear();
       }
     int getNameLen() { return namelen; }
     const char* getName()  { return (const char*) id; }
     const char* name()  { return (const char*) id; }
     const char* getSeqName()  { return (const char*) id; }
     const char* getId()  { return (const char*) id; }
     const char* getDescr()  { return (const char*) descr; }
     int getDescrLen() { return descrlen; }
     const char* getSeq() { return (const char*) seq; }
     int getSeqLen()  { return len; }
     void extendId(char c) {
         if (namelen+1 >= id_cap) {
            id_cap += CAPINC;
            GREALLOC(id, id_cap);
            }
          id[namelen]= c;
          namelen++;
          }
     void extendSeqName(char c) { extendId(c); }
     void extendName(char c) { extendId(c); }
     void extendDescr(char c) {
         if (descrlen+1 >= d_cap) {
            d_cap += CAPINC;
            GREALLOC(descr, d_cap);
            }
          descr[descrlen]= c;
          descrlen++;
          }
     void endId() {   id[namelen]=0;  }
     void endName() {  id[namelen]=0;  }
     void endSeqName() {  id[namelen]=0;  }
     void endDescr() { descr[descrlen]=0; }
     void endSeq() { seq[len]=0; }
     void extendSeq(char c) {
         if (len+1 >= s_cap) {
            s_cap += SEQCAPINC;
            GREALLOC(seq, s_cap);
            }
          seq[len]= c;
          len++;
          }
     void compactIdMem() { if (namelen>0) {
            GREALLOC(id, namelen+1); id_cap=namelen+1;
            }  }
     void compactDescrMem() { if (descrlen>0) {
         GREALLOC(descr, descrlen+1); d_cap=descrlen+1; } }
     void compactSeqMem() { if (len>0) {
         GREALLOC(seq, len+1); s_cap=len+1; }  }
     void compactMem() {
        compactIdMem();
        compactDescrMem();
        compactSeqMem();
        }
     char* detachSeqPtr() { //such that the sequence allocated memory is no longer
           // freed when the FastaSeq object is destroyed
           // the returned pointer MUST be deallocated by the the user, later!
       char* p=seq;
       GMALLOC(seq, SEQCAPINC);
       s_cap=SEQCAPINC;
       len=0;
       return p;
       }
     char* setSeqPtr(char* newseq, int newlen=0, int newcap=0) {
       if (newlen==0) newlen=strlen(newseq);
       if (newcap<=newlen) newcap=newlen+1;
       GFREE(seq);
       seq=newseq;
       len=newlen;
       s_cap=newcap;
       return seq;
       }
     void reset() {// allocated space remains the same!
       namelen=0;id[0]=0;
       descrlen=0;descr[0]=0;
       len=0;seq[0]=0;
      }
     //reverse-complement a nucleotide sequence:
     void reverseComplement() {
      if (len==0) return;
      //ntCompTableInit();
      reverseChars(seq,len);
      for (int i=0;i<len;i++) seq[i]=ntComplement(seq[i]);
     }
  //printing fasta formatted sequence to a file stream
  void fprint(FILE* fout, int line_len=60, bool defline=false) {
       if (defline) {
         if (descrlen>0) fprintf(fout, "%s %s\n", id, descr);
                    else fprintf(fout, ">%s\n", id);
         }
       int l=len;
       char* p=seq;
       while (l>0) {
         int to_write=GMIN(line_len, l);
         fwrite(p,1,to_write,fout);
         fprintf(fout,"\n");
         p+=line_len;
         l-=line_len;
         }
       }
  //
  static void write(FILE *fh, char* seqid, char* descr, char* seq,
                                  const int linelen=60, const int seqlen=0) {
      int i, idlen, ilen;
      //char *s;
      //s = (seqid == NULL) ? (char*)"ANONYMOUS" : seqid;
      // write header line only if given!
      idlen=(seqid==NULL)? 0 : strlen(seqid);
      if (idlen>0) {
        if (*seqid != '>') putc('>', fh);
        fwrite(seqid, 1, idlen, fh);
        i=(descr==NULL)? 0 : strlen(descr);
        if (i>0) {
          putc(' ',fh);
          fwrite(descr, 1, i, fh);
          }
        putc('\n', fh);
        }
      if (linelen>0) { //fixed line length
        int len = (seqlen>0) ? seqlen : strlen(seq);
        ilen=0;
        for (i=0; i < len; i++, ilen++) {
                if (ilen == linelen) {
                     putc('\n', fh);
                     ilen = 0;
                     }
                putc(seq[i], fh);
                } //for
        putc('\n', fh);
        }
       else { //no line length limit
         fprintf(fh, "\n%s\n", seq);
         }
      fflush(fh);
      }


};

typedef int charFunc(char c, int pos, FastaSeq* fseq); //char processing function
/* passes:
   c = current sequence character (generally aminoacid or nucleotide)
   pos = 0-based coordinate of the given character within the sequence
   fseq = FastaSeq pointer (useful for retrieving sequence defline info)
  the return value is not used yet
*/


//(for reading/writing variable length records, etc.)
enum fileMode {
 fmRead,
 fmWrite
 };

class GFastaFile {
  char* fname;
  FILE* fh;
  fileMode fmode;

  long int rec_fpos; //the input stream offset of the current record to be read
  long int cur_fpos;        //the input stream offset of the current byte to be read
  uint seqcoord; //1-based coordinate of the current record's sequence reading position
                 //(updated by getSeqRange() mostly)
 protected:
  void bad_fastafmt() {
      GError("Error parsing file '%s'. Not a Fasta file?\n", fname);
      }
  void check_eof(int c) {
      if (c == EOF) bad_fastafmt();
      }
 public:
  GFastaFile(const char* filename, fileMode filemode=fmRead) {
      fh=NULL;
      cur_fpos=0;
      rec_fpos=0;
      fmode=filemode;
      seqcoord=0;
      const char *mode=(filemode==fmRead) ? "rb" : "wb";
      if (filename == NULL || filename[0]=='\0') {
           fh = (filemode == fmRead) ? stdin : stdout;
           fname=NULL;
           }
         else {
           if ((fh = fopen(filename, mode)) == NULL)
               GError("Cannot open file '%s'!", filename);
           fname=Gstrdup(filename);
           }
       /*
       GCALLOC(curseqid, CAPINC);
       curseqidlen=CAPINC;
       GCALLOC(curdescr, CAPINC);
       curdescrlen=CAPINC;*/
    }

   //attach a GFastaFile object to an already open handle
   GFastaFile(FILE* fhandle, fileMode filemode=fmRead, const char* filename=NULL) {
     fh=fhandle;
     cur_fpos=ftell(fh);
     fmode=filemode;
     rec_fpos=cur_fpos;
     seqcoord=0;
     if (filename == NULL || filename[0]=='\0') {
           fname=NULL;
           }
         else
           fname=Gstrdup(filename);
     }


   void reset() {
    if (fh!=NULL && fh!=stdout && fh!=stdin) {
       fseek(fh,0L, SEEK_SET);
       cur_fpos=0;
       rec_fpos=0;
       }
     else GError("Cannot use GFastaFile::reset() on stdin, stdout or NULL handles.\n");
    }

   void seek(int pos) {
    if (fh!=NULL && fh!=stdout && fh!=stdin) {
      fseek(fh, pos, SEEK_SET);
      cur_fpos=pos;
      seqcoord=0; //seqcoord agnostic after a seek
      }
     else GError("Cannot use GFastaFile::seek() on stdin, stdout or NULL handles.\n");
    }

  ~GFastaFile() {
    if (fh!=NULL && fh!=stdout && fh!=stdin) fclose(fh);
    fh=NULL;
    GFREE(fname);
    /*GFREE(curseqid);
    GFREE(curdescr);*/
    }

   int getReadPos() { return cur_fpos; } /* returns current read position in the
              input stream (can be used within callback) */
   int ReadSeqPos() {return rec_fpos; } /* returns the input stream offset of the last fasta
                                                record processed by getFastaSeq*/
   bool readHeader(FastaSeq& seq) { return (readHeader(&seq)!=NULL); }
   FastaSeq* readSeq(int seqalloc=0) {
     //allocate a new FastaSeq, reads the next record and returns it
     //caller is responsible for deallocating returned FastaSeq memory!
     FastaSeq* r=readHeader(NULL, seqalloc);
     int len=0;
     char before=1; //newline before indicator
     int c=-1;
     //load the whole sequence in FastaSeq
      while ((c = getc(fh)) != EOF && c != '>') {
           cur_fpos++;
           //if (isspace(c) || c<31)
           if (c<=32) {
                  before = (c=='\n' || c=='\r')?1:0;
                  continue; /* skip spaces */
                  }
           if (len >= r->s_cap-1) {
                 GREALLOC(r->seq, r->s_cap + SEQCAPINC);
                 r->s_cap+=SEQCAPINC;
                 }
           r->seq[len] = c;
           before=0;
           len++;
           }
     r->seq[len] = '\0';
     r->len=len;
     return r;
    }
   FastaSeq* readHeader(FastaSeq* seq=NULL, int seqalloc=0) {
    /* reads the Fasta sequence header
    the first character must be '>' for this call, after any spaces,
    if seq is NULL a new FastaSeq object is allocated and returned,
    otherwise id and descr are updated */
     seqcoord=0;
     int* buflen;
     int* buflenstr;
     char** buf;
     int before;
     int c = getc(fh); cur_fpos++;
     while (c!=EOF && c<=32) { c=getc(fh); cur_fpos++; }//skip spaces etc.
     if (c == EOF) return NULL;
     if (c != '>')
            bad_fastafmt();
     if (seq==NULL) seq=new FastaSeq(seqalloc);
         else {  seq->clear(); seq->init(seqalloc); }

     int len = 0; //chars accumulated so far
     buflen=&(seq->id_cap);
     buf=&(seq->id);
     buflenstr=&(seq->namelen);
     before=1;
     while ((c = getc(fh)) != EOF) {
              cur_fpos++;
              if (c=='\n' || c=='\r') break;
              if (len >= *buflen-1) {
                      GREALLOC(*buf, *buflen + CAPINC);
                      *buflen+=CAPINC;
                      }
              if (before && (c<=32)) {
                 // space encountered => seq_name finished
                 before=0;
                 (*buf)[len]='\0';
                 *buflenstr=len;
                 buf=&seq->descr;
                 buflen=&seq->d_cap;
                 buflenstr=&seq->descrlen;
                 len=0;
                 if (c!=1)  // special case, nrdb concatenation
                   continue; // skip this space
                 }
              (*buf)[len]=c;
              len++;
              }
     (*buf)[len]='\0'; /* terminate the comment string */
     *buflenstr = len;
     check_eof(c); /* it's wrong to have eof here */
     seqcoord=1;
     return (seq->namelen==0) ? NULL : seq;
  }

   FastaSeq *getFastaSeq(bool& is_last, FastaSeq* seq, charFunc* callbackFn = NULL ) {
   /* seq must be a pointer to a initialized FastaSeq structure
   if seq is NULL, the sequence is not actually read,
     but just skipped and the file pointer set accordingly, while
     the returned "pointer" will not a valid one, but just NULL (or not NULL if
     end of file was encountered)
   if callbackFn is NULL, the sequence is read entirely in memory in a FastaSeq.seq field
   otherwise only the defline is parsed into FastaSeq::id and FastaSeq::descr but actual
       sequence letters are passed one by one to the callback function
      and the actual sequence is never stored in memory (unless the callback does it)
   */
      int c, len;
      int before;
      rec_fpos=cur_fpos;
      len = 0; //chars accumulated so far
      // -------- read the defline first
      if (seq==NULL) { // navigate only! don't read/parse anything but the record delimiter
          before=1;
          while ((c = getc(fh)) != EOF && c != '\n' && c !='\r') cur_fpos++; // skip defline
          check_eof(c); /* it's wrong to have eof here! */
          cur_fpos++; //to account for the '\n' read
          /*----- read the sequence now: */
          before=1; /* "newline before" flag */
          while ((c = getc(fh)) != EOF && c != '>') {
                cur_fpos++;
                before = (c=='\n' || c=='\r') ? 1 : 0;
                }
         //we should end up at a '>' character here, or EOF
          } /* fasta fmt navigation to next sequence, no seq storage */
       else { // sequence storage:
          readHeader(seq);
          /*----- read the actual sequence now: */
          len=0;
          before=1; //newline before indicator
          if (callbackFn==NULL) { //load the whole sequence in FastaSeq
             while ((c = getc(fh)) != EOF && c != '>') {
                cur_fpos++;
                //if (isspace(c) || c<31)
                if (c<=32) {
                       before = (c=='\n' || c=='\r')?1:0;
                       continue; /* skip spaces */
                       }
                if (len >= seq->s_cap-1) {
                      GREALLOC(seq->seq, seq->s_cap + CAPINC);
                      seq->s_cap+=CAPINC;
                      }
                seq->seq[len] = c;
                before=0;
                len++;
                }
             seq->seq[len] = '\0';
             seq->len=len;
             } /* sequence storage */
          else { //use the callback for each letter, do not store the whole sequence in FastaSeq
             while ((c = getc(fh)) != EOF && c != '>') {
                cur_fpos++;
                if (c<=32) {
                       before = (c=='\n' || c=='\r')?1:0;
                       continue; /* skip spaces within sequence*/
                       }
                (*callbackFn)(c, len, seq); //call the user function for each letter
                before=0;
                len++;
                }
             seq->len=len;
             } /* callback sequence reading (no storage)*/
        } /* sequence parsing */
      if (c=='>') {
         if (!before) bad_fastafmt(); /* '>' must only be at start of line,
                                       never within the sequence ! */
         is_last=false; /* FALSE - not the last one */
         ungetc(c, fh);
         }
        else is_last=true; /* TRUE - eof() here */
      return ((seq==NULL) ? (FastaSeq*)fh : seq); //alwayws return non NULL here!
  } //getFastaSeq

   //simplified call to ignore the is_last flag
  FastaSeq *getFastaSeq(FastaSeq* seq, charFunc* callbackFn = NULL) {
     bool b;
     return getFastaSeq(b, seq, callbackFn);
     }


  uint seqSkip(uint slen, int& c){
   //assumes the header was read !
   //skip exactly slen characters in the actual aa or nt sequence
                      //(spaces are not counted)
   uint skipacc=0;
   while (skipacc<slen && ((c=getc(fh))!= EOF && c != '>')) {
     cur_fpos++;
     if (c<=32) continue; //skip spaces and other non-ASCII characters
     seqcoord++;
     skipacc++;
     }
   return skipacc; //may terminate prematurely
   }

  /* read a sequence range from the current FASTA record
   this is much faster when rcoord>=seqcoord (i.e. when sequence
   ranges are read sequentially)
   if rcoord>=seqcoord assumes the header has been read already!
   Returns the actual length of the sequence returned (0 if rcoord>seq_length)
   and updates seqcoord, cur_fpos accordingly (rec_fpos is unchanged)
 */
  uint getSeqRange(FastaSeq& seq, uint rcoord, uint rlen=0) {
      int c;
      uint len;
      int before;
      rec_fpos=cur_fpos;
      if (!seqcoord || seqcoord>rcoord) {
          // slow -- go back to the beginning of the record
          seek(rec_fpos);
          readHeader(&seq); //this will also reset seqcoord to 1
          }
      if (rcoord!=seqcoord) {
         seqSkip(rcoord-seqcoord, c);
         check_eof(c);
         if (c=='>')
             GError("Error: '>' character found while skipping through sequence!\n");
         }
      len = 0; //chars accumulated so far
      seq.seq[0]='\0';
      seq.len=0;
      //----- read the actual subsequence now:
      len=0;
      before=1; //"newline before" flag
      while ((c = getc(fh)) != EOF && c != '>') {
                cur_fpos++;
                if (c<=32) continue; // skip spaces
                if (len >= (uint) (seq.s_cap-1)) {
                      GREALLOC(seq.seq, seq.s_cap + CAPINC);
                      seq.s_cap+=CAPINC;
                      }
                seq.seq[len] = c;
                len++;
                seqcoord++;
                if (rlen>0 && len==rlen) break;
                }
       seq.seq[len] = '\0';
       seq.len=len;
       if (c=='>') bad_fastafmt(); /* '>' must only be at start of line,
                                       never within the sequence ! */
      return len;
   } //getSeqRange

   //only for writing
   void putFastaSeq(FastaSeq *fa, const int linelen=60) {
      writeFasta(fh, fa->id, fa->descr, fa->seq, linelen);
      }

static void writeFasta(FILE *fh, char* seqid, char* descr, char* seq, const int linelen=60, const int seqlen=0) {
  FastaSeq::write(fh, seqid, descr, seq, linelen, seqlen);
  }

};

// ------------- FASTA parser/handler ----
// REQUIRES the first character processed after init()
// to be the first character of the record delimiter
// (default: ">")


class GFastaCharHandler {
 protected:
  char* recdelim;
  charFunc* seqCallBack;
  bool in_delim;
  int delim_pos;
  bool in_seqname;
  bool in_descr;
  bool in_seq;
  FastaSeq* rec;
  unsigned int seq_pos;
  void reset() {
     in_delim=true;
     delim_pos=0;
     in_seqname=false;
     in_descr=false;
     in_seq=false;
     seq_pos=0;
    }
 public:
  GFastaCharHandler(char* recdel=DEF_FASTA_DELIM) {
    reset();
    rec=NULL;
    recdelim=recdel;
    seqCallBack=NULL;
    }
  GFastaCharHandler(charFunc* chrCallBack, FastaSeq* r=NULL, char* recdel=DEF_FASTA_DELIM) {
    reset();
    rec=r;
    recdelim=recdel;
    seqCallBack=chrCallBack;
    if (rec!=NULL) rec->reset();
    }
   void init() {
     init(rec, seqCallBack);
     }
   void init(charFunc* chrCallBack) {
     init(rec,chrCallBack);
     }
   void init(FastaSeq* r) {
     init(r,seqCallBack);
     }
   void init(FastaSeq* r, charFunc* chrCallBack) {
     rec=r;
     seqCallBack=chrCallBack;
     if (rec==NULL)
        GError("GFastaCharHandler::init() Error: cannot use NULL FastaSeq!\n");
     rec->reset();
     reset();
    }
  void done() {
    if (rec==NULL)
        GError("GFastaCharHandler::done() Error: cannot use NULL FastaSeq!\n");
    rec->endId();
    rec->endDescr();
    rec->endSeq();
    }

  //~GFastaCharHandler();

  void processChar(char c) {
    if (in_delim) { //skip record delimiter -- but it must be there!
      if (recdelim[delim_pos]!=c) {//the only way to detect an Id starting
         in_seqname=true;
         in_delim=false;
         }
      delim_pos++;
      }
    if (in_seqname) {
      if (rec->namelen>0 && c<=32) {
         //breaking out of seq_name
         rec->endId();
         if (c=='\n' || c=='\r') { //end defline
             in_seqname=false;
             in_seq=true;
             }
           else { //seqname break, not defline end
             in_seqname=false;
             in_descr=true;
             }
       } // seqname termination
      else { //seqname continues
       if (c>32) rec->extendId(c);
       }
      return;
      } // in_seqname
    if (in_descr) {
         if (c=='\n' || c=='\r') { //end defline
             rec->endDescr();
             in_descr=false;
             in_seq=true;
             }
          else rec->extendDescr(c);
      return;
      } // in_descr
    if (in_seq && c>32) {
       seq_pos++;  // 1-based sequence position !
       if (seqCallBack==NULL) rec->extendSeq(c);
                   else (*seqCallBack)(c,seq_pos,rec);
       }
   }

 };


#endif
