#include "gff.h"
#include "GStr.h"
#include "GArgs.h"
#include "GHash.hh"
#include "GList.hh"
#include "GFastaIndex.h"
#include "GFaSeqGet.h"
#include <ctype.h>
// don't care about cdb compression
//#ifdef ENABLE_COMPRESSION
//#undef ENABLE_COMPRESSION
//#endif
//#include "GCdbYank.h"

#define USAGE "Usage:\n\
gffread <input_gff> [-g <genomic_seqs_fasta> | <dir>][-s <seq_info.fsize>] \n\
 [-o <outfile.gff>] [-t <tname>] [-r [[<strand>]<chr>:]<start>..<end>] \n\
 [-CTVNMAFGRUVBHSZWTOE] [-w <spl_exons.fa>] [-x <spl_cds.fa>] [-y <tr_cds.fa>]\n\
 [-i <maxintron>] \n\
 Filters and/or converts GFF3/GTF2 records.\n\
 <input_gff> is a GFF file, use '-' if the GFF records will be given at stdin\n\
 \n\
 Options:\n\
  -g  full path to a multi-fasta file with the genomic sequences\n\
      for all input mappings, OR a directory with single-fasta files\n\
      (one per genomic sequence, with file names matching sequence names)\n\
  -s  <seq_info.fsize> is a tab-delimited file providing this info\n\
      for each of the mapped sequences:\n\
      <seq-name> <seq-length> <seq-description>\n\
      (useful for mRNA/EST/protein mappings with -A option)\n\
  -i  discard transcripts having an intron larger than <maxintron>\n\
  -r  only show transcripts crossing coordinate range <start>..<end>\n\
      (on chromosome/contig <chr>, strand <strand> if provided)\n\
  -R  for -r option, discard all transcripts that are not fully \n\
      contained within given range\n\
  -U  discard single-exon transcripts\n\
  -C  discard mRNAs that have no CDS feature\n\
  -F  keep all attributes from last column of GFF/GTF\n\
  -G  only parse additional exon attributes from the first exon\n\
      and move them to the mRNA level (useful for GTF input)\n\
  -A  use the description field from <seq_info.fsize> and add it\n\
      as the value for a 'descr' attribute to the GFF record\n\
  \n\
  -O  process non-transcript GFF records as well (by default non-transcript\
      records are ignored).\n\
  -V  discard any mRNAs with CDS having in-frame stop codons\n\
  -H  for -V option, check and adjust the starting CDS phase\n\
      if the original phase leads to a translation with an \n\
      in-frame stop codon\n\
  -B  for -V option, single-exon transcripts are also checked on the\n\
      opposite strand\n\
  -N  only show multi-exon mRNAs if all their introns have the \n\
      typical splice site consensus ( GT-AG, GC-AG or AT-AC )\n\
  -M  discard any mRNAs that either lack initial START codon\n\
      or the terminal STOP codon, or have an in-frame stop codon\n\
      (only print mRNAs with a fulll, valid CDS)\n\
 \n\
  -E  expose (warn about) duplicate transcript IDs and other potential \n\
      problems with the input GFF/GTF records\n\
  -S  sort output GFF records by genomic sequence and start coordinate\n\
      (this option is automatically enabled by -g option)\n\
  -Z  merge close exons into a single exon (for intron size<4)\n\
  -w  write a fasta file with spliced exons for each GFF transcript\n\
  -x  write a fasta file with spliced CDS for each GFF transcript\n\
  -W  for -w and -x options, also write for each fasta record the exon\n\
      coordinates projected onto the spliced sequence\n\
  -y  write a protein fasta file with the translation of CDS for each record\n\
  -o  the \"filtered\" GFF records will be written to <outfile.gff>\n\
      (use -o- for printing to stdout)\n\
  -t  use <trackname> in the second column of each GFF output line\n\
  -T  -o option will output GTF format instead of GFF3\n\
 "

FILE* ffasta=NULL;
FILE* f_in=NULL;
FILE* f_out=NULL;
FILE* f_w=NULL; //fasta with spliced exons (transcripts)
FILE* f_x=NULL; //fasta with spliced CDS
FILE* f_y=NULL; //fasta with translated CDS
bool wCDSonly=false;

bool validCDSonly=false; // translation with no in-frame STOP
bool bothStrands=false; //for single-exon mRNA validation, check the other strand too
bool altPhases=false; //if original phase fails translation validation,
                     //try the other 2 phases until one makes it
bool mRNAOnly=true; 
bool spliceCheck=false; //only known splice-sites

bool fullCDSonly=false; // starts with START, ends with STOP codon
bool fullattr=false;
bool sortByLoc=false; // if the GFF output should be sorted by location
//GStr gseqpath;
//GStr gcdbfa;
//bool multiGSeq=false; //if a directory or a .cidx file was given to -g option
//GFaSeqGet* faseq=NULL;
//GCdbYank* gcdb=NULL;
//int gseq_id=-1; //current genome sequence ID -- the current GffObj::gseq_id
bool fmtGTF=false;
bool addDescr=false;
//bool protmap=false;
bool multiExon=false;
bool writeExonSegs=false;
char* tracklabel=NULL;
int maxintron=999000000;
bool mergeCloseExons=false;
//range filter:
char* rfltGSeq=NULL;
char rfltStrand=0;
uint rfltStart=0;
uint rfltEnd=MAX_UINT;
bool rfltWithin=false; //check for full containment within given range
bool noExonAttr=false;
class SeqInfo {
 public:
  int len;
  char* descr;
  SeqInfo( int l, char* s) {
   len=l;
   if (s==NULL) {
     descr=NULL;
     }   else {
     descr=Gstrdup(s);
     }
   }
  ~SeqInfo() {
   GFREE(descr);
   }
};

char* getGSeqName(int gseq_id) {
 return GffObj::names->gseqs.getName(gseq_id);
}

//genomic fasta sequence handling
class GFastaHandler {
 public:
  char* fastaPath;
  GFastaIndex* faIdx; //could be a cdb .cidx file
  int last_fetchid;
  GFaSeqGet* faseq;
  //GCdbYank* gcdb;
  char* getFastaFile(int gseq_id) {
     if (fastaPath==NULL) return NULL;
     GStr s(fastaPath);
     s.trimR('/');
     s.appendfmt("/%s",getGSeqName(gseq_id));
     GStr sbase(s);
     if (!fileExists(s.chars())) s.append(".fa");
     if (!fileExists(s.chars())) s.append("sta");
     if (fileExists(s.chars())) return Gstrdup(s.chars());
         else {
             GMessage("Warning: cannot find genomic sequence file %s{.fa,.fasta}\n",sbase.chars());
             return NULL;
             }
     }

   GFastaHandler(const char* fpath=NULL) {
     //gcdb=NULL;
     fastaPath=NULL;
     faseq=NULL;
     faIdx=NULL;
     init(fpath);
     }

   void init(const char* fpath) {
     if (fpath==NULL || fpath[0]==0) return;
     last_fetchid=-1;
     if (!fileExists(fpath))
       GError("Error: file/directory %s does not exist!\n",fpath);
     fastaPath=Gstrdup(fpath);
     GStr gseqpath(fpath);
     /*
     if (gseqpath.rindex(".cidx")==gseqpath.length()-5) {
        //cdbyank index given directly
        gcdb=new GCdbYank(gseqpath.chars());
        if (fileExists(gcdb->getDbName())) {
            gseqpath=gcdb->getDbName();
            } else {
            gseqpath.chomp(".cidx");
            if (!fileExists(gseqpath.chars()))
                GError("Error: cannot locate the fasta file for index %s.cidx !\n",gseqpath.chars());
            }
        GFREE(fastaPath);
        fastaPath=Gstrdup(gseqpath.chars());
        return;
        }
        */
     if (fileExists(fastaPath)>1) { //exists and it's not a directory
            GStr fainame(fastaPath);
            if (fainame.rindex(".fai")==fainame.length()-4) {
               //.fai index file given directly
               fastaPath[fainame.length()-4]=0;
               if (!fileExists(fastaPath))
                  GError("Error: cannot find fasta file for index %s !\n", fastaPath);
               }
              else fainame.append(".fai");
            //GMessage("creating GFastaIndex with fastaPath=%s, fainame=%s\n", fastaPath, fainame.chars());
            faIdx=new GFastaIndex(fastaPath,fainame.chars());
            GStr fainamecwd(fainame);
            int ip=-1;
            if ((ip=fainamecwd.rindex(CHPATHSEP))>=0)
               fainamecwd.cut(0,ip+1);
            if (!faIdx->hasIndex()) { //could not load index
               //try current directory
                  if (fainame!=fainamecwd) {
                    if (fileExists(fainamecwd.chars())>1) {
                       faIdx->loadIndex(fainamecwd.chars());
                       }
                    }
                  } //tried to load index
            if (!faIdx->hasIndex()) {
                 GMessage("No fasta index found for %s. Rebuilding, please wait..\n",fastaPath);
                 faIdx->buildIndex();
                 if (faIdx->getCount()==0) GError("Error: no fasta records found!\n");
                 GMessage("Fasta index rebuilt.\n");
                 FILE* fcreate=fopen(fainame.chars(), "w");
                 if (fcreate==NULL) {
                   GMessage("Warning: cannot create fasta index %s! (permissions?)\n", fainame.chars());
                   if (fainame!=fainamecwd) fcreate=fopen(fainamecwd.chars(), "w");
                   if (fcreate==NULL)
                      GError("Error: cannot create fasta index %s!\n", fainamecwd.chars());
                   }
                 if (faIdx->storeIndex(fcreate)<faIdx->getCount())
                     GMessage("Warning: error writing the index file!\n");
                 } //index created and attempted to store it
            } //multi-fasta
     }
   GFaSeqGet* fetch(int gseq_id, bool checkFasta=false) {
     if (fastaPath==NULL) return NULL;
     if (gseq_id==last_fetchid && faseq!=NULL) return faseq;
     delete faseq;
     faseq=NULL;
     last_fetchid=-1;
     char* gseqname=getGSeqName(gseq_id);
     // DEBUG:
     //GMessage("..processing transcripts on: %s\n",gseqname);
     //genomic sequence given
     /*
     if (gcdb!=NULL) {
       uint32 reclen=0;
       off_t rpos=gcdb->getRecordPos(gseqname, &reclen);
       if (rpos<0) // genomic sequence not found
          GError("Error: cannot find genomic sequence '%s' in %s\n",gseqname, fastaPath);
       // WARNING: does not validate FASTA line-len uniformity!
       faseq=new GFaSeqGet(fastaPath,rpos, false);
       faseq->loadall(reclen); //load the whole sequence, it's faster
       last_fetchid=gseq_id;
       return faseq;
       }
       */
     if (faIdx!=NULL) { //fastaPath was the multi-fasta file name
        GFastaRec* farec=faIdx->getRecord(gseqname);
        if (farec!=NULL) {
             faseq=new GFaSeqGet(fastaPath,farec->seqlen, farec->fpos,
                               farec->line_len, farec->line_blen);
             faseq->loadall(); //just cache the whole sequence, it's faster
             last_fetchid=gseq_id;
             }
        else {
          GMessage("Warning: couldn't find fasta record for '%s'!\n",gseqname);
          return NULL;
          }
        }
     else {
         char* sfile=getFastaFile(gseq_id);
         if (sfile!=NULL) {
            faseq=new GFaSeqGet(sfile,checkFasta);
            faseq->loadall();
            last_fetchid=gseq_id;
            GFREE(sfile);
            }
         } //one fasta file per contig
       return faseq;
     }

   ~GFastaHandler() {
     GFREE(fastaPath);
     //delete gcdb;
     delete faIdx;
     delete faseq;
     }
};


class GSpliceSite {
 public:
  char nt[3];
  GSpliceSite(const char* c, bool revc=false) {
    nt[2]=0;
    if (c==NULL) {
      nt[0]=0;
      nt[1]=0;
      return;
      }
    if (revc) {
      nt[0]=toupper(ntComplement(c[1]));
      nt[1]=toupper(ntComplement(c[0]));
      }
    else {
      nt[0]=toupper(c[0]);
      nt[1]=toupper(c[1]);
      }
    }

  GSpliceSite(const char* intron, int intronlen, bool getAcceptor, bool revc=false) {
    nt[2]=0;
    if (intron==NULL || intronlen==0)
       GError("Error: invalid intron or intron len for GSpliceSite()!\n");
    const char* c=intron;
    if (revc) {
      if (!getAcceptor) c+=intronlen-2;
      nt[0]=toupper(ntComplement(c[1]));
      nt[1]=toupper(ntComplement(c[0]));
      }
    else { //on forward strand
      if (getAcceptor) c+=intronlen-2;
      nt[0]=toupper(c[0]);
      nt[1]=toupper(c[1]);
      }//forward strand
    }

  GSpliceSite(const char n1, const char n2) {
    nt[2]=0;
    nt[0]=toupper(n1);
    nt[1]=toupper(n2);
    }
  bool canonicalDonor() {
    return (nt[0]=='G' && (nt[1]=='C' || nt[1]=='T'));
    }
  bool operator==(GSpliceSite& c) {
    return (c.nt[0]==nt[0] && c.nt[1]==nt[1]);
    }
  bool operator==(GSpliceSite* c) {
    return (c->nt[0]==nt[0] && c->nt[1]==nt[1]);
    }
  bool operator==(const char* c) {
    //return (nt[0]==toupper(c[0]) && nt[1]==toupper(c[1]));
    //assumes given const nucleotides are uppercase already!
    return (nt[0]==c[0] && nt[1]==c[1]);
    }
  bool operator!=(const char* c) {
    //assumes given const nucleotides are uppercase already!
    return (nt[0]!=c[0] || nt[1]!=c[1]);
    }
};

//hash with sequence info
GHash<SeqInfo> seqinfo;
GHash<int> isoCounter; //counts the valid isoforms

//bool debugMode=false;
bool verbose=false;

char* getSeqDescr(char* seqid) {
 static char charbuf[128];
 if (seqinfo.Count()==0) return NULL;
 char* suf=rstrchr(seqid, '.');
 if (suf!=NULL) *suf=0;
 SeqInfo* seqd=seqinfo.Find(seqid);
 if (suf!=NULL) *suf='.';
 if (seqd!=NULL) {
  GStr s(seqd->descr);
  //cleanup some Uniref gunk
  if (s[0]=='[') {
    int r=s.index(']');
    if (r>=0 && r<8 && isdigit(s[1]))
       s.remove(0,r+1);
    }
  if (s.length()>80) {
    int r=s.index(';');
    if (r>5) s.cut(r);
    }
  if (s.length()>127) {
   s.cut(127);
   int r=s.rindex(' ');
   if (r>0) s.cut(r);
   }
  strcpy(charbuf, s.chars());
  return charbuf;
  }
 else return NULL;
}

char* getSeqName(char* seqid) {
  static char charbuf[128];
  char* suf=rstrchr(seqid, '.');
  if (suf!=NULL) *suf=0;
  strcpy(charbuf, seqid);
  if (suf!=NULL) *suf='.';
  return charbuf;
}

void printFasta(FILE* f, GStr& defline, char* seq, int seqlen=-1) {
 if (seq==NULL) return;
 int len=(seqlen>0)?seqlen:strlen(seq);
 if (len<=0) return;
 if (!defline.is_empty())
     fprintf(f, ">%s\n",defline.chars());
 int ilen=0;
 for (int i=0; i < len; i++, ilen++) {
   if (ilen == 70) {
     fputc('\n', f);
     ilen = 0;
     }
   putc(seq[i], f);
   } //for
 fputc('\n', f);
}

void loadSeqInfo(FILE* f, GHash<SeqInfo> &si) {
  GLineReader fr(f);
  while (!fr.isEof()) {
      char* line=fr.getLine();
      if (line==NULL) break;
      char* id=line;
      char* lenstr=NULL;
      char* text=NULL;
      char* p=line;
      while (*p!=0 && !isspace(*p)) p++;
      if (*p==0) continue;
      *p=0;p++;
      while (*p==' ' || *p=='\t') p++;
      if (*p==0) continue;
      lenstr=p;
      while (*p!=0 && !isspace(*p)) p++;
      if (*p!=0) { *p=0;p++; }
      while (*p==' ' || *p=='\t') p++;
      if (*p!=0) text=p; //else text remains NULL
      int len=0;
      if (!parseInt(lenstr,len)) {
         GMessage("Warning: could not parse sequence length: %s %s\n",
                  id, lenstr);
         continue;
         }
      // --- here we have finished parsing the line
      si.Add(id, new SeqInfo(len,text));
      } //while lines
}

GFaSeqGet* fastaSeqGet(GFastaHandler& gfasta, GffObj& mrna) {
  if (gfasta.fastaPath==NULL) return NULL;
  return gfasta.fetch(mrna.gseq_id);
}

int adjust_stopcodon(GffObj& mrna, int adj, GList<GSeg>* seglst=NULL) {
 //adj>0 => extedn CDS,  adj<0 => shrink CDS
 //when CDS is expanded, exons have to be checked too and 
 // expanded accordingly if they had the same boundary
  int realadj=0;
  if (mrna.strand=='-') {
       if ((int)mrna.CDstart>adj) {

           mrna.CDstart-=adj;
           realadj=adj;
           if (adj<0) { //restore
              if (mrna.exons.First()->start==mrna.CDstart+adj) {
                 mrna.exons.First()->start-=adj;
                 mrna.start=mrna.exons.First()->start;
                 mrna.covlen+=adj;
                 }
              }
           else if (mrna.exons.First()->start>=mrna.CDstart) {
                 mrna.exons.First()->start-=adj;
                 mrna.start=mrna.exons.First()->start;
                 mrna.covlen+=adj;
                 }
             }
          }
        else {
         realadj=adj;
         mrna.CDend+=adj;
         if (adj<0) {//restore
           if (mrna.exons.Last()->end==mrna.CDend-adj) {
                        mrna.exons.Last()->end+=adj;
                        mrna.end=mrna.exons.Last()->end;
                        mrna.covlen+=adj;
                        }
          }
         else if (mrna.exons.Last()->end<=mrna.CDend) {
             mrna.exons.Last()->end+=adj;
             mrna.end=mrna.exons.Last()->end;
             mrna.covlen+=adj;
             }
         }
  if (seglst!=NULL) seglst->Last()->end+=adj;
  return realadj;
 }

void process_mRNA(GFastaHandler& gfasta, GffObj& mrna) {
  if (f_out!=NULL && !mRNAOnly && mrna.exons.Count()==0) {
     //a gene or other generic feature without exons
     if (fmtGTF) mrna.printGtf(f_out, tracklabel);
            else mrna.printGff(f_out, tracklabel);
     return;
     }
 GStr gene(mrna.getGene());
 GStr defline(mrna.getID());
 if (!gene.is_empty() && strcmp(gene.chars(),mrna.getID())!=0) {
   int* isonum=isoCounter.Find(gene.chars());
   if (isonum==NULL) {
      isonum=new int(1);
      isoCounter.Add(mrna.getGene(),isonum);
      }
     else (*isonum)++;
  defline.appendfmt(" gene=%s", gene.chars());
  }
  int seqlen=0;

  const char* tlabel=tracklabel;
  if (tlabel==NULL) tlabel=mrna.getTrackName();
  //defline.appendfmt(" track:%s",tlabel);
  char* cdsnt = NULL;
  char* cdsaa = NULL;
  int aalen=0;
  for (int i=1;i<mrna.exons.Count();i++) {
     int ilen=mrna.exons[i]->start-mrna.exons[i-1]->end-1;
     if (ilen>2000000) 
            GMessage("Warning: very large intron (%d) for transcript %s\n",
                           ilen, mrna.getID());
     if (ilen>maxintron) {
         return;
         }
     }
  GList<GSeg> seglst(false,true);
  GFaSeqGet* faseq=fastaSeqGet(gfasta, mrna);
  if (spliceCheck && mrna.exons.Count()>1) {
    //check introns for splice site consensi ( GT-AG, GC-AG or AT-AC )
    if (faseq==NULL) GError("Error: no genomic sequence available!\n");
    int glen=mrna.end-mrna.start+1;
    const char* gseq=faseq->subseq(mrna.start, glen);
    bool revcompl=(mrna.strand=='-');
    bool ssValid=true;
    for (int e=1;e<mrna.exons.Count();e++) {
      const char* intron=gseq+mrna.exons[e-1]->end+1-mrna.start;
      int intronlen=mrna.exons[e]->start-mrna.exons[e-1]->end-1;
      GSpliceSite acceptorSite(intron,intronlen,true, revcompl);
      GSpliceSite    donorSite(intron,intronlen, false, revcompl);
      //GMessage("%c intron %d-%d : %s .. %s\n",
      //           mrna.strand, istart, iend, donorSite.nt, acceptorSite.nt);
      if (acceptorSite=="AG") { // GT-AG or GC-AG
         if (!donorSite.canonicalDonor()) {
            ssValid=false;break;
            }
         }
      else if (acceptorSite=="AC") { //
         if (donorSite!="AT") { ssValid=false; break; }
         }
      else { ssValid=false; break; }
      }
    //GFREE(gseq);
    if (!ssValid) {
      if (verbose)
         GMessage("Invalid splice sites found for '%s'\n",mrna.getID());
      return; //don't print this one!
      }
    }

  bool trprint=true;
  int stopCodonAdjust=0;
  int mCDphase=0;
  bool hasStop=false;
  if (mrna.CDphase=='1' || mrna.CDphase=='2')
      mCDphase = mrna.CDphase-'0';
  if (f_y!=NULL || f_x!=NULL || validCDSonly) {
    if (faseq==NULL) GError("Error: no genomic sequence provided!\n");
    //if (protmap && fullCDSonly) {
    //if (protmap && (fullCDSonly ||  (mrna.qlen>0 && mrna.qend==mrna.qlen))) {
    
    if (validCDSonly) { //make sure the stop codon is always included 
      //adjust_stopcodon(mrna,3);
      stopCodonAdjust=adjust_stopcodon(mrna,3);
      }
    int strandNum=0;
    int phaseNum=0;
  CDS_CHECK:
    cdsnt=mrna.getSpliced(faseq, true, &seqlen,NULL,NULL,&seglst);
    if (cdsnt==NULL) trprint=false;
    if (validCDSonly) {
       cdsaa=translateDNA(cdsnt, aalen, seqlen);
       char* p=strchr(cdsaa,'.');
       hasStop=false;
       if (p!=NULL) {
            if (p-cdsaa>=aalen-2) { //stop found as the last codon
                    *p='0';//remove it
                    hasStop=true;
                    if (aalen-2==p-cdsaa) {
                      //previous to last codon is the stop codon
                      //so correct the CDS stop accordingly
                      adjust_stopcodon(mrna,-3, &seglst);
                      stopCodonAdjust=0; //clear artificial stop adjustment
                      seqlen-=3;
                      cdsnt[seqlen]=0;
                      }
                    aalen=p-cdsaa;
                    }
                 else {//stop found before the last codon
                    trprint=false;
                    }
            }//stop codon found
       if (trprint==false) { //failed CDS validity check
         //in-frame stop codon found
         if (altPhases && phaseNum<3) {
            phaseNum++;
            mrna.CDphase = '0'+((mCDphase+phaseNum)%3);
            GFREE(cdsaa);
            goto CDS_CHECK;
            }
         if (mrna.exons.Count()==1 && bothStrands) {
            strandNum++;
            phaseNum=0;
            if (strandNum<2) {
               GFREE(cdsaa);
               mrna.strand = (mrna.strand=='-') ? '+':'-';
               goto CDS_CHECK; //repeat the CDS check for a different frame
               }
            }
         if (verbose) GMessage("In-frame STOP found for '%s'\n",mrna.getID());
         } //has in-frame STOP
       if (fullCDSonly) {
           if (!hasStop || cdsaa[0]!='M') trprint=false;
           }
       } // CDS check requested
    } //translation or codon check/output was requested
  if (!trprint) {
    GFREE(cdsnt);
    GFREE(cdsaa);
    return;
    }
  if (stopCodonAdjust>0 && !hasStop) {
          //restore stop codon location
          adjust_stopcodon(mrna, -stopCodonAdjust, &seglst);
          if (cdsnt!=NULL && seqlen>0) {
             seqlen-=stopCodonAdjust;
             cdsnt[seqlen]=0;
             }
          if (cdsaa!=NULL) aalen--;
          }

  if (f_out!=NULL) {
     if (fmtGTF) mrna.printGtf(f_out, tracklabel);
            else mrna.printGff(f_out, tracklabel);
     }
  if (f_y!=NULL) { //CDS translation fasta output requested
         //char* 
         if (cdsaa==NULL) { //translate now if not done before
           cdsaa=translateDNA(cdsnt, aalen, seqlen);
           }
         if (fullattr && mrna.attrs!=NULL) {
             //append all attributes found for each transcripts
              for (int i=0;i<mrna.attrs->Count();i++) {
                defline.append(" ");
                defline.append(mrna.getAttrName(i));
                defline.append("=");
                defline.append(mrna.getAttrValue(i));
                }
              }
         printFasta(f_y, defline, cdsaa, aalen);
         }
   if (f_x!=NULL) { //CDS only
         if (writeExonSegs) {
              defline.append(" loc:");
              defline.append(mrna.getGSeqName());
              defline.appendfmt("(%c)",mrna.strand);
              //warning: not CDS coordinates are written here, but the exon ones
              defline+=(int)mrna.start;
              defline+=(char)'-';
              defline+=(int)mrna.end;
              // -- here these are CDS substring coordinates on the spliced sequence:
              defline.append(" segs:");
              for (int i=0;i<seglst.Count();i++) {
                  if (i>0) defline.append(",");
                  defline+=(int)seglst[i]->start;
                  defline.append("-");
                  defline+=(int)seglst[i]->end;
                  }
              }
         if (fullattr && mrna.attrs!=NULL) {
             //append all attributes found for each transcript
              for (int i=0;i<mrna.attrs->Count();i++) {
                defline.append(" ");
                defline.append(mrna.getAttrName(i));
                defline.append("=");
                defline.append(mrna.getAttrValue(i));
                }
              }
         printFasta(f_x, defline, cdsnt, seqlen);
         }
 GFREE(cdsnt);
 GFREE(cdsaa);
 if (f_w!=NULL) { //write spliced exons
    uint cds_start=0;
    uint cds_end=0;
    seglst.Clear();
    char* exont=mrna.getSpliced(faseq, false, &seqlen, &cds_start, &cds_end, &seglst);
    if (exont!=NULL) {
    if (mrna.CDstart>0) {
        defline.appendfmt(" CDS=%d-%d", cds_start, cds_end);
        }
      if (writeExonSegs) {
        defline.append(" loc:");
        defline.append(mrna.getGSeqName());
        defline+=(char)'|';
        defline+=(int)mrna.start;
        defline+=(char)'-';
        defline+=(int)mrna.end;
        defline+=(char)'|';
        defline+=(char)mrna.strand;
        defline.append(" exons:");
        for (int i=0;i<mrna.exons.Count();i++) {
                if (i>0) defline.append(",");
                defline+=(int)mrna.exons[i]->start;
                defline.append("-");
                defline+=(int)mrna.exons[i]->end;
                }
        defline.append(" segs:");
        for (int i=0;i<seglst.Count();i++) {
            if (i>0) defline.append(",");
            defline+=(int)seglst[i]->start;
            defline.append("-");
            defline+=(int)seglst[i]->end;
            }
        }
      if (fullattr && mrna.attrs!=NULL) {
       //append all attributes found for each transcripts
        for (int i=0;i<mrna.attrs->Count();i++) {
          defline.append(" ");
          defline.append(mrna.getAttrName(i));
          defline.append("=");
          defline.append(mrna.getAttrValue(i));
          }
        }
      printFasta(f_w, defline, exont, seqlen);
      GFREE(exont);
      }
    } //writing f_w (spliced exons)
 return;
}

void openfw(FILE* &f, GArgs& args, char opt) {
  GStr s=args.getOpt(opt);
  if (!s.is_empty()) {
      if (s=='-')
       f=stdout;
      else {
       f=fopen(s,"w");
       if (f==NULL) GError("Error creating file: %s\n", s.chars());
       }
     }
}

#define FWCLOSE(fh) if (fh!=NULL && fh!=stdout) fclose(fh)
#define FRCLOSE(fh) if (fh!=NULL && fh!=stdin) fclose(fh)


int main(int argc, char * const argv[]) {
 GArgs args(argc, argv, "hvOUNHWCVMNSXTDAPRZFGEg:i:r:s:t:a:b:o:w:x:y:MINCOV=MINPID=");
 int e;
 if ((e=args.isError())>0)
    GError("%s\nInvalid argument: %s\n", USAGE, argv[e]);
 if (args.getOpt('h')!=NULL) GError("%s\n", USAGE);
 //debugMode=(args.getOpt('D')!=NULL);
 mRNAOnly=(args.getOpt('O')==NULL);
 sortByLoc=(args.getOpt('S')!=NULL);
 addDescr=(args.getOpt('A')!=NULL);
 verbose=(args.getOpt('v')!=NULL);
 wCDSonly=(args.getOpt('C')!=NULL);
 validCDSonly=(args.getOpt('V')!=NULL);
 altPhases=(args.getOpt('H')!=NULL);
 fmtGTF=(args.getOpt('T')!=NULL); //switch output format to GTF
 bothStrands=(args.getOpt('B')!=NULL);
 fullCDSonly=(args.getOpt('M')!=NULL);
 spliceCheck=(args.getOpt('N')!=NULL);
 //protmap=(args.getOpt('P')!=NULL);
 if (fullCDSonly) validCDSonly=true;
 fullattr=(args.getOpt('F')!=NULL);
 if (args.getOpt('G')==NULL) 
    noExonAttr=!fullattr;
   else {
     noExonAttr=true;
     fullattr=true;
     }
 mergeCloseExons=(args.getOpt('Z')!=NULL);
 multiExon=(args.getOpt('U')!=NULL);
 writeExonSegs=(args.getOpt('W')!=NULL);
 tracklabel=args.getOpt('t');
 GFastaHandler gfasta(args.getOpt('g'));
 if (gfasta.fastaPath!=NULL) sortByLoc=true; //enforce sorting by chromosome/contig
 GStr s=args.getOpt('i');
 if (!s.is_empty()) maxintron=s.asInt();
 rfltWithin=(args.getOpt('R')!=NULL);
 s=args.getOpt('r');
 if (!s.is_empty()) {
   s.trim();
   if (s[0]=='+' || s[0]=='-') {
     rfltStrand=s[0];
     s.cut(0,1);
     }
   int isep=s.index(':');
   if (isep>0) { //gseq name given
      if (rfltStrand==0 && (s[isep-1]=='+' || s[isep-1]=='-')) {
        isep--;
        rfltStrand=s[isep];
        s.cut(isep,1);
        }
      if (isep>0) 
          rfltGSeq=Gstrdup((s.substr(0,isep)).chars());
      s.cut(0,isep+1);
      }
   GStr gsend;
   char slast=s[s.length()-1];
   if (rfltStrand==0 && (slast=='+' || slast=='-')) {
      s.chomp(slast);
      rfltStrand=slast;
      }
   if (s.index("..")>=0) gsend=s.split("..");
                    else gsend=s.split('-');
   if (!s.is_empty()) rfltStart=(uint)s.asInt();
   if (!gsend.is_empty()) {
      rfltEnd=(uint)gsend.asInt();
      if (rfltEnd==0) rfltEnd=MAX_UINT;
      }
   
   } //gseq/range filtering
 else {
   if (rfltWithin)
     GError("Error: option -R doesn't make sense without -r!\n");
   }
 s=args.getOpt('s');
 if (!s.is_empty()) {
  FILE* fsize=fopen(s,"r");
  if (fsize==NULL) GError("Error opening info file: %s\n",s.chars());
  loadSeqInfo(fsize, seqinfo);
 }
 /*
 openfw(fgtfok, args, 'a');
 openfw(fgtfbad, args, 'b');
 */
 openfw(f_out, args, 'o');
 //if (f_out==NULL) f_out=stdout;
 if (gfasta.fastaPath==NULL && (validCDSonly || spliceCheck || args.getOpt('w')!=NULL || args.getOpt('x')!=NULL || args.getOpt('y')!=NULL))
  GError("Error: -g option is required for options -w, -x, -y, -V, -N, -M !\n");

 openfw(f_w, args, 'w');
 openfw(f_x, args, 'x');
 openfw(f_y, args, 'y');
 if (f_y!=NULL || f_x!=NULL) wCDSonly=true;
 //useBadCDS=useBadCDS || (fgtfok==NULL && fgtfbad==NULL && f_y==NULL && f_x==NULL);
 GStr infile;
 if (args.startNonOpt()) {
        infile=args.nextNonOpt();
        //GMessage("Given file: %s\n",infile.chars());
        }
 if (!infile.is_empty()) {
    if (infile=="-") f_in=stdin;
      else 
        if ((f_in=fopen(infile, "r"))==NULL)
            GError("Cannot open input file %s!\n",infile.chars());
    }
  else
    f_in=stdin;

 //GffReader* gfreader=new GffReader(f_in,true);
 GffReader* gfreader=new GffReader(f_in, mRNAOnly, sortByLoc);
 if (args.getOpt('E')!=NULL) gfreader->showWarnings(true);
 gfreader->readAll(fullattr, mergeCloseExons, noExonAttr);

 //if (debugMode) GMessage("[D] gff data loaded, now processing each entry\n");
 //GList<GffObj> mrnas(true,true,false);
 for (int m=0;m<gfreader->gflst.Count();m++) {
   GffObj* mrna=gfreader->gflst[m];
   if (mrna->hasErrors() || (mrna->len()+500>GFF_MAX_LOCUS)) { //should probably report these in a file too..
			GMessage("Warning: transcript %s discarded (structural errors found, length=%d).\n", 
			                                        mrna->getID(), mrna->len());
			//gfreader->gflst.freeItem(m);
			continue;
			}
	//GStr feature(mrna->getFeatureName());
	//feature.lower();
	//bool gene_or_locus=(feature.endsWith("gene") ||feature.index("loc")>=0);
	
	//if (mRNAOnly && mrna->monoFeature() && !mrna->isTranscript() &&
	//     (mrna->exons.Count()==0 || gene_or_locus)) {
    if (mRNAOnly && mrna->isDiscarded()) {
	   //discard generic "gene" or "locus" features with no other detailed subfeatures
	   //GMessage("Warning: discarding %s GFF generic gene/locus container %s\n",m->getID());
	   continue;
	   }
   //if (mrna->exons.Count()==0 && (mrna->isTranscript() ||
    //    (!mRNAOnly && !gene_or_locus))) {
    if (mrna->exons.Count()==0) {
      //a non-mRNA feature with no subfeatures
      //just so we get some sequence functions working, add a dummy "exon"-like subfeature here
      mrna->addExon(mrna->start,mrna->end);
      }
   if (rfltGSeq!=NULL) { //filter by gseqName
      if (strcmp(mrna->getGSeqName(),rfltGSeq)!=0) {
        continue;
        }
      }
   if (rfltStrand>0 && mrna->strand !=rfltStrand) {
      continue;
      }
   //check coordinates
   if (rfltStart!=0 || rfltEnd!=MAX_UINT) {
     if (rfltWithin) {
       if (mrna->start<rfltStart || mrna->end>rfltEnd) {
          continue;
          }
       }
     else {
       if (mrna->start>rfltEnd || mrna->end<rfltStart) {
         continue;
         }
       }
     }
   if (multiExon && mrna->exons.Count()<=1) {
       continue;
       }
   if (wCDSonly && mrna->CDstart==0) {
       continue;
       }
   /* if (validCDSonly && mrna->hasErrors) {
       delete mrna;
       gfreader->gflst.Forget(m);
       continue;
       }
   */
   process_mRNA(gfasta, *mrna);
   }
 // M_END:
 delete gfreader;
 seqinfo.Clear();
 //if (faseq!=NULL) delete faseq;
 //if (gcdb!=NULL) delete gcdb;
 GFREE(rfltGSeq);
 FRCLOSE(f_in);
 FWCLOSE(f_out);
 FWCLOSE(f_w);
 FWCLOSE(f_x);
 FWCLOSE(f_y);
 }


