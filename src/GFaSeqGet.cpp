#include "GFaSeqGet.h"
#include "gdna.h"
#include <ctype.h>

void GSubSeq::setup(uint sstart, int slen, int sovl, int qfrom, int qto) {
     if (sovl==0) {
       GFREE(sq);
       sqstart=sstart;
       sqlen = (slen==0 ? MAX_FASUBSEQ : slen);
       GMALLOC(sq, sqlen);
       return;
       }
  //overlap -- copy the overlapping region
  char* newsq=NULL;
  GMALLOC(newsq, slen);
  memcpy((void*)&newsq[qto], (void*)&sq[qfrom], sovl);
  GFREE(sq);
  sq=newsq;
  sqstart=sstart;
  sqlen=slen;
}

void GFaSeqGet::finit(const char* fn, off_t fofs, bool validate) {
 fh=fopen(fn,"rb");
 if (fh==NULL) {
   GError("Error (GFaSeqGet) opening file '%s'\n",fn);
   }
 fname=Gstrdup(fn);
 initialParse(fofs, validate);
 lastsub=new GSubSeq();
}

GFaSeqGet::GFaSeqGet(FILE* f, off_t fofs, bool validate) {
 if (f==NULL) GError("Error (GFaSeqGet) : null file handle!\n");
 fh=f;
 initialParse(fofs, validate);
 lastsub=new GSubSeq();
}

void GFaSeqGet::initialParse(off_t fofs, bool checkall) {
 static const char gfa_ERRPARSE[]="Error (GFaSeqGet): invalid FASTA file format.\n";
 if (fofs!=0) { fseek(fh,fofs,SEEK_SET); } //e.g. for offsets provided by cdbyank
 //read the first two lines to determine fasta parameters
 fseqstart=fofs;
 int c=getc(fh);
 fseqstart++;
 if (c!='>') GError("Error (GFaSeqGet): not a fasta header?\n");
 while ((c=getc(fh))!=EOF) {
   fseqstart++;
   if (c=='\n' || c=='\r') { break; } //end of defline
   }

 if (c==EOF) GError(gfa_ERRPARSE);
 linelen=0;
 lendlen=0;
 lendch=0;
 while ((c=getc(fh))!=EOF) {
  if (c=='\n' || c=='\r') { //end of line encountered
     if (linelen>0) { //end of the first "sequence" line
        lendch=c;
        lendlen++;
        break;
        }
      else {// another eol char at the end of defline
        fseqstart++;
        continue;
        }
     }// end-of-line characters
  linelen++;
  }
 //we are at the end of first sequence line
 while ((c=getc(fh))!=EOF) {
   if (c=='\n' || c=='\r') lendlen++;
      else {
       ungetc(c,fh);
       break;
       }
   }
 if (c==EOF) return;
 // -- you don't need to check it all if you're sure it's safe
 if (checkall) { //validate the rest of the FASTA record
   int llen=0; //last line length
   int elen=0; //length of last line ending
   bool waseol=true;
   while ((c=getc(fh))!=EOF) {
     if (c=='>' && waseol) { ungetc(c,fh); break; }
     if (c=='\n' ||  c=='\r') {
        // eol char
        elen++;
        if (waseol) continue; //2nd eol char
        waseol=true;
        elen=1;
        continue;
        }
     if (c<=32) GError(gfa_ERRPARSE); //invalid character encountered
     //--- on a seq char here:
     if (waseol) {//beginning of a seq line
       if (elen && (llen!=linelen || elen!=lendlen))
           //GError(gfa_ERRPARSE);
         GError("Error: invalid FASTA format for GFaSeqGet; make sure that\n\
  the sequence lines have the same length (except for the last line)");
       waseol=false;
       llen=0;
       elen=0;
       }
     llen++;
     } //while reading chars
   }// FASTA checking was requested
 fseek(fh,fseqstart,SEEK_SET);
}

const char* GFaSeqGet::subseq(uint cstart, int& clen) {
  //cstart is 1-based genomic coordinate within current fasta sequence
  if (clen>MAX_FASUBSEQ) {
    unsigned int maxlen=MAX_FASUBSEQ;
    GMessage("Error (GFaSeqGet): subsequence cannot be larger than %d\n", maxlen);
    return NULL;
    }
  if (lastsub->sq==NULL || lastsub->sqlen==0) {
    lastsub->setup(cstart, clen);
    loadsubseq(cstart, clen);
    lastsub->sqlen=clen;
    return (const char*)lastsub->sq;
    }
  //allow extension up to MAX_FASUBSEQ
  uint bstart=lastsub->sqstart;
  uint bend=lastsub->sqstart+lastsub->sqlen-1;
  uint cend=cstart+clen-1;
  int qlen=0; //only the extra len to be allocated/appended/prepended
  uint qstart=cstart; //start coordinate of the new seq block of length qlen to be read from file
  int newlen=0; //the new total length of the buffered sequence lastsub->sq
  int kovl=0;
  int czfrom=0;//0-based offsets for copying a previously read sequence chunk
  int czto=0;
  uint newstart=cstart;
  if (cstart>=bstart && cend<=bend) { //new reg contained within existing buffer
     return (const char*) &(lastsub->sq[cstart-bstart]) ;
    }
  //extend downward
  uint newend=GMAX(cend, bend);
  if (cstart<bstart) { //requested start < old buffer start
    newstart=cstart;
    newlen=(newend-newstart+1);
    if (newlen>MAX_FASUBSEQ) {
       newlen=MAX_FASUBSEQ;
       newend=cstart+newlen-1; //keep newstart, set newend
       }
    qlen=bstart-cstart;
    if (newend>bstart) { //overlap
       if (newend>bend) {// new region is larger & around the old one - so we have two regions to update
         kovl=bend-bstart+1;
         czfrom=0;
         czto=bstart-cstart;
         lastsub->setup(newstart, newlen, kovl, czfrom, czto); //this should realloc and copy the kovl subseq
         qlen=bstart-cstart;
         loadsubseq(newstart, qlen);
         qlen=newend-bend;
         int toread=qlen;
         loadsubseq(bend+1, qlen);
         clen-=(toread-qlen);
         lastsub->sqlen=clen;
         return (const char*)lastsub->sq;
         }
        //newend<=bend
       kovl=newend-bstart+1;
       }
     else { //no overlap with previous buffer
       if (newend>bend) kovl=bend-bstart+1;
                   else kovl=newend-bstart+1;
       }
     qlen=bstart-cstart;
     czfrom=0;
     czto=qlen;
    } //cstart<bstart
   else { //cstart>=bstart, possibly extend upwards
    newstart=bstart;
    newlen=(newend-newstart+1);
    if (newlen>MAX_FASUBSEQ) {
       newstart=bstart+(newlen-MAX_FASUBSEQ);//keep newend, assign newstart
       newlen=MAX_FASUBSEQ;
       if (newstart<=bend) { //overlap with old buffer
          kovl=bend-newstart+1;
          czfrom=newstart-bstart;
          czto=0;
          }
       else { //not overlapping old buffer
         kovl=0;
         }
       } //newstart reassigned
    else { //we can extend the buffer to include the old one
      qlen=newend-bend; //how much to read from file
      qstart=bend+1;
      kovl=bend-bstart+1;
      czfrom=0;
      czto=0;
      }
    }
  lastsub->setup(newstart, newlen, kovl, czfrom, czto); //this should realloc but copy any overlapping region
  lastsub->sqlen-=qlen; //appending may result in a premature eof
  int toread=qlen;
  loadsubseq(qstart, qlen); //read the missing chunk, if any
  clen-=(toread-qlen);
  lastsub->sqlen+=qlen;
  return (const char*)(lastsub->sq+(cstart-newstart));
}

char* GFaSeqGet::copyRange(uint cstart, uint cend, bool revCmpl, bool upCase) {
  if (cstart>cend) { swap(cstart, cend); }
  int clen=cend-cstart+1;
  const char* gs=subseq(cstart, clen);
  if (gs==NULL) return NULL;
  char* r=NULL;
  GMALLOC(r,clen+1);
  r[clen]=0;
  memcpy((void*)r,(void*)gs, clen);
  if (revCmpl) reverseComplement(r,clen);
  if (upCase) {
       for (int i=0;i<clen;i++)
            r[i]=toupper(r[i]);
       }
  return r;
 }

const char* GFaSeqGet::loadsubseq(uint cstart, int& clen) {
  //assumes enough lastsub->sq space allocated previously
  //only loads the requested clen chars from file, at offset &lastsub->sq[cstart-lastsub->sqstart]
  int sofs=cstart-lastsub->sqstart;
  char* seqp=lastsub->sq+sofs;
  //find the proper file offset and read the appropriate lines
  uint seqofs=cstart-1;
  uint startlno = seqofs/linelen;
  int lineofs = seqofs % linelen;
  off_t fstart=fseqstart+startlno * (linelen+lendlen);
  fstart+=lineofs;
  fseek(fh, fstart, SEEK_SET);
  int toread=clen;
  if (toread==0) toread=MAX_FASUBSEQ; //read max allowed, or to the end of file
  int actualrlen=0;
  int sublen=0;
  if (lineofs>0) { //read the partial first line
    int reqrlen=linelen-lineofs;
    if (reqrlen>toread) reqrlen=toread; //in case we need to read just a few chars
    actualrlen=fread((void*)seqp, 1, reqrlen, fh);
    if (actualrlen<reqrlen) { //eof reached prematurely
      while (seqp[actualrlen-1]=='\n' || seqp[actualrlen-1]=='\r') actualrlen--;
      //check for new sequences in between
      clen=actualrlen;
      sublen+=actualrlen;
      return (const char*)seqp;
      }
    toread-=reqrlen;
    sublen+=reqrlen;
    fseek(fh, lendlen, SEEK_CUR);
    }
  //read the rest of the lines
  while (toread>=linelen) {
    char* rseqp=&(seqp[sublen]);
    actualrlen=fread((void*)rseqp, 1, linelen, fh);
    if (actualrlen<linelen) {
      while (rseqp[actualrlen-1]=='\n' || rseqp[actualrlen-1]=='\r') actualrlen--;
      sublen+=actualrlen;
      clen=sublen;
      return (const char*)seqp;
      }
    toread-=actualrlen;
    sublen+=actualrlen;
    fseek(fh, lendlen, SEEK_CUR);
    }
  // read the last partial line, if any
  if (toread>0) {
    char* rseqp=&(seqp[sublen]);
    actualrlen=fread((void*)rseqp, 1, toread, fh);
    if (actualrlen<toread) {
      while (rseqp[actualrlen-1]=='\n' || rseqp[actualrlen-1]=='\r')
          actualrlen--;
      }
    //toread-=actualrlen;
    sublen+=actualrlen;
    }
  //lastsub->sqlen+=sublen;
  clen=sublen;
  return (const char*)seqp;
  }


