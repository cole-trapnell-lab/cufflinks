#ifndef GFF_H
#define GFF_H

#include "GBase.h"
#include "gdna.h"
#include "codons.h"
#include "GFaSeqGet.h"
#include "GList.hh"
#include "GHash.hh"

/*
const byte exMskMajSpliceL = 0x01;
const byte exMskMajSpliceR = 0x02;
const byte exMskMinSpliceL = 0x04;
const byte exMskMinSpliceR = 0x08;
const byte exMskTag = 0x80;
*/

//reserved Gffnames::feats entries -- basic feature types
extern const int gff_fid_mRNA;
extern const int gff_fid_exon;
extern const int gff_fid_CDS; //never really used, except for display only
                              //use gff_fid_exon instead
//extern bool gff_warns; //show parser warnings - now GffReader member

extern const uint GFF_MAX_LOCUS;
extern const uint GFF_MAX_EXON;
extern const uint GFF_MAX_INTRON;

#define GFF_LINELEN 2048
#define ERR_NULL_GFNAMES "Error: GffObj::%s requires a non-null GffNames* names!\n"


enum GffExonType {
  exgffNone=0,
  exgffStart, //from "start_codon" feature (within CDS)
  exgffStop, //from "stop_codon" feature (may be outside CDS)
  exgffCDS,  //from "CDS" feature
  exgffUTR,  //from "UTR" feature
  exgffExon, //from "exon" feature
};

class GffReader;

class GffLine {
 public:
    char* line;
    int llen;
    char* gseqname;
    char* track;
    char* ftype; //feature name: mRNA/gene/exon/CDS
    char* info; //the last, attributes' field, unparsed
    uint fstart;
    uint fend;
    uint qstart; //overlap coords on query, if available
    uint qend;
    uint qlen; //query len, if given
    double score;
    char strand;
    bool skip;
    bool is_cds; //"cds" and "stop_codon" features
    bool is_exon; //"exon" and "utr" features
    char exontype; // gffExonType
    bool is_mrna;
    char phase;  // '.' , '0', '1' or '2'
    // -- allocated strings:
    char* gname; //gene_id or Name= value (or an otherwise parsed gene denominator
                 //(for grouping isoforms)
    //
    char* Parent; // if a Parent=.. attribute was found
    char* ID;     // if a ID=.. attribute was parsed, or a GTF with 'transcript' line (transcript_id)
    GffLine(GffReader* reader, const char* l); //parse the line accordingly
    void discardParent() { GFREE(Parent); }
    GffLine(GffLine* l) { //a copy constructor
      memcpy((void*)this, (void*)l, sizeof(GffLine));
      line=NULL;
      GMALLOC(line, llen+1);
      memcpy((void*)line, (void*)(l->line), llen+1);
      //--offsets withing line[]
      gseqname=line+(l->gseqname-l->line);
      track=line+(l->track-l->line);
      ftype=line+(l->ftype-l->line);
      info=line+(l->info-l->line);
      //-- allocated string copies:
      Parent=Gstrdup(l->Parent);
      ID=Gstrdup(ID);
      gname=Gstrdup(gname);
      }
    GffLine() {
      line=NULL;
      gseqname=NULL;
      track=NULL;
      ftype=NULL;
      fstart=0;
      fend=0;
      info=NULL;
      Parent=NULL;
      ID=NULL;
      gname=NULL;
      skip=true;
      qstart=0;
      qend=0;
      qlen=0;
      exontype=0;
      is_cds=false;
      is_mrna=false;
      is_exon=false;
      }
    ~GffLine() {
      GFREE(line);
      GFREE(Parent);
      GFREE(ID);
      GFREE(gname);
     }
};

class GffAttr {
 public:
   int attr_id;
   char* attr_val;
   GffAttr(int an_id, const char* av=NULL) {
     attr_id=an_id;
     attr_val=NULL;
     if (av!=NULL) {
       char* lastch = (char*)(av+strlen(av)-1);
       //remove spaces at the end:
       while (*lastch==' ' && lastch!=av) lastch--;
       lastch[1]=0;
       //practical usage: if it doesn't have any spaces just strip those useless double quotes
       if (av[0]=='"' && strpbrk(av+1," ;")==NULL) {
               if (*lastch=='"') *lastch=0;
               attr_val=Gstrdup(av+1);
               }
          else attr_val=Gstrdup(av);
       }

     }
  ~GffAttr() {
     GFREE(attr_val);
     }
  bool operator==(GffAttr& d){
      return (this==&d);
      }
  bool operator>(GffAttr& d){
     return (this>&d);
     }
  bool operator<(GffAttr& d){
    return (this<&d);
    }

 };

class GffNameList;
class GffNames;

class GffNameInfo {
  friend class GffNameList;
protected:
   //unsigned char shared;
   int idx;
 public:
   char* name;
   GffNameInfo() { name=NULL; idx=-1; }
   GffNameInfo(const char* n) {
     name=Gstrdup(n);
     }
   ~GffNameInfo() {
      //if (shared==0)
          GFREE(name);
    }

   bool operator==(GffNameInfo& d){
       return (strcmp(this->name, d.name)==0);
       }
   bool operator>(GffNameInfo& d){
      return (strcmp(this->name, d.name)>0);
      }
   bool operator<(GffNameInfo& d){
     return (strcmp(this->name, d.name)<0);
     }
};

class GffNameList:public GList<GffNameInfo> {
  friend class GffNameInfo;
  friend class GffNames;
protected:
  GHash<GffNameInfo> byName;//hash with shared keys
  int idlast; //fList index of last added/reused name
  void addStatic(const char* tname) {// fast add
     GffNameInfo* f=new GffNameInfo(tname);
     idlast=this->Add(f);
     f->idx=idlast;
     byName.shkAdd(f->name,f);
     }
public:
 GffNameList():GList<GffNameInfo>(false,true,true), byName(false) {
    idlast=-1;
    }
 char* lastNameUsed() { return idlast<0 ? NULL : Get(idlast)->name; }
 int lastNameId() { return idlast; }
 char* getName(int nid) { //retrieve name by its ID
   if (nid<0 || nid>=fCount)
         GError("GffNameList Error: invalid index (%d)\n",nid);
   return fList[nid]->name;
   }

 int addName(const char* tname) {//returns or create an id for the given name
   //check idlast first, chances are it's the same feature name checked
   if (idlast>=0 && strcmp(fList[idlast]->name,tname)==0)
       return idlast;
   GffNameInfo* f=byName.Find(tname);
   int fidx=-1;
   if (f!=NULL) fidx=f->idx;
     else {//add new entry
      f=new GffNameInfo(tname);
      fidx=this->Add(f);
      f->idx=fidx;
      byName.shkAdd(f->name,f);
      }
   idlast=fidx;
   return fidx;
   }
 int getId(const char* tname) { //only returns a name id# if found
    GffNameInfo* f=byName.Find(tname);
    if (f==NULL) return -1;
    return f->idx;
    }
 int removeName() {
   GError("Error: removing names from GffNameList not allowed!\n");
   return -1;
   }
};

class GffNames {
 public:
   int numrefs;
   GffNameList tracks;
   GffNameList gseqs;
   GffNameList attrs;
   GffNameList feats; //feature names - anything except 'mRNA', 'exon', 'CDS' gets stored here
   GffNames():tracks(),gseqs(),attrs(), feats() {
    numrefs=0;
    //the order below is critical!
    //has to match: gff_fid_mRNA, gff_fid_exon, gff_fid_CDS
    feats.addStatic("mRNA");//index 0=gff_fid_mRNA
    feats.addStatic("exon");//index 1=gff_fid_exon
    feats.addStatic("CDS"); //index 2=gff_fid_CDS
    }
};

void gffnames_ref(GffNames* &n);
void gffnames_unref(GffNames* &n);

enum GffPrintMode {
  pgtfAny, //print record as read
  pgtfExon,
  pgtfCDS,
  pgffAny, //print record as read
  pgffExon,
  pgffCDS,
  pgffBoth,
};


class GffAttrs:public GList<GffAttr> {
  public:
    GffAttrs():GList<GffAttr>(false,true,false) { }
    char* getAttr(GffNames* names, const char* attrname) {
      int aid=names->attrs.getId(attrname);
      if (aid>=0)
        for (int i=0;i<Count();i++)
          if (aid==Get(i)->attr_id) return Get(i)->attr_val;
      return NULL;
      }
};


class GffExon : public GSeg {
 public:
  void* uptr; //for later extensions
  GffAttrs* attrs; //other attributes kept for this exon
  double score; // gff score column
  char phase; //GFF phase column - for CDS segments only
             // '.' = undefined (UTR), '0','1','2' for CDS exons
  char exontype; // 1="exon" 2="cds" 3="utr" 4="stop_codon"
  int qstart; // for mRNA/protein exon mappings: coordinates on query
  int qend;
  GffExon(int s=0, int e=0, double sc=0, char fr=0, int qs=0, int qe=0, char et=0) {
    uptr=NULL;
    attrs=NULL;
    if (s<e) {
      start=s;
      end=e;
      }
   else {
     start=e;
     end=s;
    }
   if (qs<qe) {
     qstart=qs;
     qend=qe;
     } else {
     qstart=qe;
     qend=qs;
     }
   score=sc;
   phase=fr;
   exontype=et;
   } //constructor

 char* getAttr(GffNames* names, const char* atrname) {
   if (attrs==NULL || names==NULL || atrname==NULL) return NULL;
   return attrs->getAttr(names, atrname);
   }

 ~GffExon() { //destructor
   if (attrs!=NULL) delete attrs;
   }
};


class GffCDSeg:public GSeg {
 public:
  char phase;
  int exonidx;
};
//one GFF mRNA object -- e.g. a mRNA with its exons and/or CDS segments
class GffObj:public GSeg {
  //utility segment-merging function for addExon() 
  void expandExon(int xovl, uint segstart, uint segend,
       char exontype, double sc, char fr, int qs, int qe);
 protected:
   //coordinate transformation data:
   uint xstart; //absolute genomic coordinates of reference region
   uint xend;
   char xstatus; //coordinate transform status:
            //0 : (start,end) coordinates are absolute
            //'+' : (start,end) coords are relative to xstart..xend region
            //'-' : (start,end) are relative to the reverse complement of xstart..xend region
   //--
   char* gffID; // ID name for mRNA (parent) feature
   char* gname; // value of Name or "gene_id" attribute (if given)
   //-- friends:
   friend class GffReader;
   friend class GffExon;
public:
  bool hasErrors; //overlapping exons, or too short introns, etc.
  static GffNames* names; // common string storage that holds the various attribute names etc.
  int track_id; // index of track name in names->tracks
  int gseq_id; // index of genomic sequence name in names->gseqs
  int ftype_id; // index of this record's feature name in names->feats, or the special gff_fid_mRNA value
  int subftype_id; //index of child subfeature name in names->feats (that subfeature stored in "exons")
                   //if ftype_id==gff_fid_mRNA then this value is ignored
  GList<GffExon> exons; //for non-mRNA entries, these can be any subfeature of type subftype_id
  int udata; //user data, flags etc.
  void* uptr; //user pointer (to a parent object, cluster, locus etc.)
  GffObj* ulink; //link to another GffObj (user controlled field)
  // mRNA specific fields:
  bool isCDS; //just a CDS, no UTRs
  bool partial; //partial CDS
  uint CDstart; //CDS start coord
  uint CDend;   //CDS end coord
  char CDphase; //initial phase for CDS start

  bool ismRNA() { return (ftype_id==gff_fid_mRNA); }

  int addExon(uint segstart, uint segend, double sc=0, char fr='.',
             int qs=0, int qe=0, bool iscds=false, char exontype=0);

  int addExon(GffReader* reader, GffLine* gl, bool keepAttr=false, bool noExonAttr=true);
  void removeExon(int idx);
  char  strand; //true if features are on the reverse complement strand
  double gscore;
  double uscore; //custom, user-computed score, if needed
  int covlen; //total coverage of  reference genomic sequence (sum of maxcf segment lengths)

   //--------- optional data:
  int qlen; //query length, start, end - if available
  int qstart;
  int qend;
  int qcov; //query coverage - percent
  GffAttrs* attrs; //other gff3 attributes found for the main mRNA feature
   //constructor by gff line parsing:
  GffObj(GffReader* gfrd, GffLine* gffline, bool keepAttrs=false, bool noExonAttr=true);
   //if gfline->Parent!=NULL then this will also add the first sub-feature
   // otherwise, only the main feature is created
  void clearAttrs() {
    if (attrs!=NULL) {
      bool sharedattrs=(exons.Count()>0 && exons[0]->attrs==attrs);
      delete attrs; attrs=NULL;
      if (sharedattrs) exons[0]->attrs=NULL;
      }
    }
  GffObj(char* anid=NULL):GSeg(0,0), exons(true,true,false) { //exons: sorted, free, non-unique
       gffID=NULL;
       uptr=NULL;
       ulink=NULL;
       udata=0;
       hasErrors=false;
       ftype_id=-1;
       subftype_id=-1;
       if (anid!=NULL) gffID=Gstrdup(anid);
       gffnames_ref(names);
       qlen=0;
       qstart=0;
       qend=0;
       qcov=0;
       partial=true;
       isCDS=false;
       CDstart=0; // hasCDS <=> CDstart>0
       CDend=0;
       CDphase=0;
       gseq_id=-1;
       track_id=-1;
       xstart=0;
       xend=0;
       xstatus=0;
       strand=0;
       gscore=0;
       uscore=0;
       attrs=NULL;
       covlen=0;
       gname=NULL;
       }
   ~GffObj() {
       GFREE(gffID);
       GFREE(gname);
       clearAttrs();
       gffnames_unref(names);
       }
   //--------------
   GffObj* finalize(GffReader* gfr, bool mergeCloseExons=false); //finalize parsing: must be called in order to merge adjacent/close proximity subfeatures
   void parseAttrs(GffAttrs*& atrlist, char* info, bool noExonAttr=false);
   const char* getSubfName() { //returns the generic feature type of the entries in exons array
     int sid=subftype_id;
     if (sid==gff_fid_exon && isCDS) sid=gff_fid_CDS;
     return names->feats.getName(sid);
     }
   bool monoFeature() {
     return (exons.Count()==0 || 
          (exons.Count()==1 && subftype_id==ftype_id));
     }
   const char* getFeatureName() {
     return names->feats.getName(ftype_id);
     }
   void addAttr(const char* attrname, char* attrvalue);
   int removeAttr(const char* attrname, const char* attrval=NULL);
   const char* getAttrName(int i) {
     if (attrs==NULL) return NULL;
     return names->attrs.getName(attrs->Get(i)->attr_id);
     }
   char* getAttr(const char* attrname, bool checkFirstExon=false) {
     if (names==NULL || attrname==NULL) return NULL;
     char* r=NULL;
     if (attrs==NULL) {
         if (!checkFirstExon) return NULL;
         }
       else r=attrs->getAttr(names, attrname);
     if (r!=NULL) return r;
     if (checkFirstExon && exons.Count()>0) {
        r=exons[0]->getAttr(names, attrname);
        }
     return r; 
     }

   char* getExonAttr(GffExon* exon, const char* attrname) {
      if (exon==NULL || attrname==NULL) return NULL;
      return exon->getAttr(names, attrname);
      }

   char* getExonAttr(int exonidx, const char* attrname) {
      if (exonidx<0 || exonidx>=exons.Count() || attrname==NULL) return NULL;
      return exons[exonidx]->getAttr(names, attrname);
      }

   char* getAttrValue(int i) {
     if (attrs==NULL) return NULL;
     return attrs->Get(i)->attr_val;
     }
   const char* getGSeqName() {
     return names->gseqs.getName(gseq_id);
     }
   const char* getTrackName() {
     return names->tracks.getName(track_id);
     }
   bool exonOverlap(uint s, uint e) {//check if ANY exon overlaps given segment
      //ignores strand!
      if (s>e) swap(s,e);
      for (int i=0;i<exons.Count();i++) {
         if (exons[i]->overlap(s,e)) return true;
         }
      return false;
      }
    bool exonOverlap(GffObj& m) {//check if ANY exon overlaps given segment
      //if (gseq_id!=m.gseq_id) return false;
      // ignores strand and gseq_id, must check in advance
      for (int i=0;i<exons.Count();i++) {
         for (int j=0;j<m.exons.Count();j++) {
            if (exons[i]->start>m.exons[j]->end) continue;
            if (m.exons[j]->start>exons[i]->end) break;
            //-- overlap if we are here:
            return true;
            }
         }
      return false;
      }
    
    int exonOverlapIdx(uint s, uint e, int* ovlen=NULL) {
      //return the exons' index for the overlapping exon
      //ovlen, if given, will return the overlap length
      if (s>e) swap(s,e);
      for (int i=0;i<exons.Count();i++) {
            if (exons[i]->start>e) break;
            if (s>exons[i]->end) continue;
            //-- overlap if we are here:
            if (ovlen!=NULL) {
              int ovlend= (exons[i]->end>e) ? e : exons[i]->end;
              *ovlen= ovlend - ((s>exons[i]->start)? s : exons[i]->start)+1;
              }
            return i;
            } //for each exon
      *ovlen=0;
      return -1;
      }
     
    int exonOverlapLen(GffObj& m) {
      if (start>m.end || m.start>end) return 0;
      int i=0;
      int j=0;
      int ovlen=0;
      while (i<exons.Count() && j<m.exons.Count()) {
        uint istart=exons[i]->start;
        uint iend=exons[i]->end;
        uint jstart=m.exons[j]->start;
        uint jend=m.exons[j]->end;
        if (istart>jend) { j++; continue; }
        if (jstart>iend) { i++; continue; }
        //exon overlap
        uint ovstart=GMAX(istart,jstart);
        if (iend<jend) {
           ovlen+=iend-ovstart+1;
           i++;
           }
        else {
           ovlen+=jend-ovstart+1;
           j++;
           }
        }//while comparing exons
      return ovlen;
      }

    bool exonOverlap(GffObj* m) {
      return exonOverlap(*m);
      }
   //---------- coordinate transformation
   void xcoord(uint grstart, uint grend, char xstrand='+') {
     //relative coordinate transform, and reverse-complement transform if xstrand is '-'
     //does nothing if xstatus is the same already
     if (xstatus) {
          if (xstatus==xstrand && grstart==xstart && grend==xend) return;
          unxcoord();//restore original coordinates
          }
     xstatus=xstrand;
     xstart=grstart;
     xend=grend;
     if (CDstart>0) xcoordseg(CDstart, CDend);
     for (int i=0;i<exons.Count();i++) {
         xcoordseg(exons[i]->start, exons[i]->end);
         }
     if (xstatus=='-') {
       exons.Reverse();
       int flen=end-start;
       start=xend-end+1;
       end=start+flen;
       }
      else {
       start=start-xstart+1;
       end=end-xstart+1;
       }
     }

   //transform an arbitrary segment based on current xstatus/xstart-xend
   void xcoordseg(uint& segstart, uint &segend) {
     if (xstatus==0) return;
     if (xstatus=='-') {
       int flen=segend-segstart;
       segstart=xend-segend+1;
       segend=segstart+flen;
       return;
       }
     else {
       segstart=segstart-xstart+1;
       segend=segend-xstart+1;
       }
     }

   void unxcoord() { //revert back to absolute genomic/gff coordinates if xstatus==true
     if (xstatus==0) return; //nothing to do, no transformation appplied
     if (CDstart>0) unxcoordseg(CDstart, CDend);
     //restore all GffExon intervals too
     for (int i=0;i<exons.Count();i++) {
         unxcoordseg(exons[i]->start, exons[i]->end);
         }
     if (xstatus=='-') {
        exons.Reverse();
        int flen=end-start;
        start=xend-end+1;
        end=start+flen;
        }
      else {
        start=start+xstart-1;
        end=end+xstart-1;
        }
     xstatus=0;
     }
   void unxcoordseg(uint& astart, uint &aend) {
     //restore an arbitrary interval -- does NOT change the transform state!
     if (xstatus==0) return;
     if (xstatus=='-') {
        int flen=aend-astart;
        astart=xend-aend+1;
        aend=astart+flen;
        }
      else {
        astart=astart+xstart-1;
        aend=aend+xstart-1;
        }
     }
   //---------------------
   bool operator==(GffObj& d){
       return (gseq_id==d.gseq_id && start==d.start && end==d.end && strcmp(gffID, d.gffID)==0);
       }
   bool operator>(GffObj& d){
      if (gseq_id!=d.gseq_id) return (gseq_id>d.gseq_id);
      if (start==d.start) {
         if (end==d.end) return strcmp(gffID, d.gffID)>0;
                      else return end>d.end;
         } else return (start>d.start);
      }
   bool operator<(GffObj& d){
     if (gseq_id!=d.gseq_id) return (gseq_id<d.gseq_id);
     if (start==d.start) {
        if (end==d.end) return strcmp(gffID, d.gffID)<0;
                     else return end<d.end;
        } else return (start<d.start);
     }
   char* getID() { return gffID; }
   char* getGene() { return gname; }
   int addSeg(GffLine* gfline);
   int addSeg(int fnid, GffLine* gfline);
   void getCDSegs(GArray<GffCDSeg>& cds);

   void updateExonPhase(); //for CDS-only features, updates GExon::phase

   void printGxfLine(FILE* fout, char* tlabel, char* gseqname, bool iscds,
                                uint segstart, uint segend, int exidx, char phase, bool gff3);
   void printGxf(FILE* fout, GffPrintMode gffp=pgffExon, char* tlabel=NULL);
   void printGtf(FILE* fout, char* tlabel=NULL) {
      printGxf(fout, pgtfAny, tlabel);
      }
   void printGff(FILE* fout, char* tlabel=NULL) {
      printGxf(fout, pgffAny, tlabel);
      }
   void print_mRNAGff(FILE* fout, char* tlabel=NULL, bool showCDS=false) {
      if (ismRNA())
         printGxf(fout, showCDS ? pgffBoth : pgffExon, tlabel);
      }
   void printSummary(FILE* fout=NULL);
   void getCDS_ends(uint& cds_start, uint& cds_end);
   void mRNA_CDS_coords(uint& cds_start, uint& cds_end);
   char* getSpliced(GFaSeqGet* faseq, bool CDSonly=false, int* rlen=NULL,
           uint* cds_start=NULL, uint* cds_end=NULL, GList<GSeg>* seglst=NULL);
    char* getUnspliced(GFaSeqGet* faseq, int* rlen, GList<GSeg>* seglst);
   char* getSplicedTr(GFaSeqGet* faseq, bool CDSonly=true, int* rlen=NULL);
   //bool validCDS(GFaSeqGet* faseq); //has In-Frame Stop Codon ?
   bool empty() { return (start==0); }
};

typedef bool GffRecFunc(GffObj* gobj, void* usrptr1, void* usrptr2);
//user callback after parsing a mapping object:
// Returns: "done with it" status:
//   TRUE if gobj is no longer needed so it's FREEd upon return
//   FALSE if the user needs the gobj pointer and is responsible for
//             collecting and freeing all GffObj objects


//GSeqStat: collect basic stats about a common underlying genomic sequence
//          for multiple GffObj
class GSeqStat {
 public:
   int gseqid; //gseq id in the global static pool of gseqs
   char* gseqname; //just a pointer to the name of gseq
   //int fcount;//number of features on this gseq
   uint mincoord;
   uint maxcoord;
   uint maxfeat_len; //maximum feature length on this genomic sequence
   GffObj* maxfeat; 
   //GList<GffObj> gflst;
   GSeqStat(int id=-1, char* name=NULL) {
        //:gflst(true,false,false) {
     gseqid=id;
     gseqname=name;
     mincoord=MAXUINT;
     maxcoord=0;
     maxfeat_len=0;
     maxfeat=NULL;
     }
   bool operator>(GSeqStat& g) {
    return (gseqid>g.gseqid);
    }
   bool operator<(GSeqStat& g) {
    return (gseqid<g.gseqid);
    }
   bool operator==(GSeqStat& g) {
    return (gseqid==g.gseqid);
    }
};


int gfo_cmpByLoc(const pointer p1, const pointer p2);

class GfList: public GList<GffObj> {
 //just adding the option to sort by genomic sequence and coordinate
   bool mustSort;
 public:
   GfList(bool sortbyloc=false):GList<GffObj>(false,false,false) {
    mustSort=sortbyloc;
    }
   void finalize(GffReader* gfr, bool mergeCloseExons) { //if set, enforce sort by locus
     if (mustSort) { //force (re-)sorting
        this->setSorted(false);
        this->setSorted((GCompareProc*)gfo_cmpByLoc);
        }
     for (int i=0;i<Count();i++) {
       fList[i]->finalize(gfr, mergeCloseExons); 
       //finalize the parsing for this GffObj -- this might also merge close exons/features if requested
       }
     }
};

class GfoHolder {
 public:
   int idx; //position in GffReader::gflst array
   GffObj* gffobj; 
   GfoHolder(GffObj* gfo=NULL, int i=0) {
     idx=i;
     gffobj=gfo;
     }
};

class CNonExon { //utility class used in subfeature promotion
 public:
   int idx;
   GffObj* parent;
   GffExon* exon;
   GffLine* gffline;
   CNonExon(int i, GffObj* p, GffExon* e, GffLine* gl) {
     parent=p;
     exon=e;
     idx=i;
     gffline=new GffLine(gl);
     }
  ~CNonExon() {
     delete gffline;
     }
 };


class GffReader {
  friend class GffObj;
  friend class GffLine;
  char* linebuf;
  off_t fpos;
  int buflen;
 protected:
  bool gff_warns; //warn about duplicate IDs, even when they are on different chromosomes
  FILE* fh;
  char* fname;  //optional fasta file with the underlying genomic sequence to be attached to this reader
  GffNames* names; //just a pointer to the global static Gff names repository in GffObj
  GffLine* gffline;
  bool mrnaOnly; //read only mRNAs ? (exon/CDS features only)
  //bool sortbyloc; //sort by location: genomic sequence and start coordinate
  GHash<GfoHolder> phash; //transcript_id+contig (Parent~Contig) => [gflst index, GffObj]
  GHash<int> tids; //transcript_id uniqueness
  char* gfoBuildId(const char* id, const char* ctg);
  void gfoRemove(const char* id, const char* ctg);
  GfoHolder* gfoAdd(const char* id, const char* ctg, GffObj* gfo, int idx);
  GfoHolder* gfoFind(const char* id, const char* ctg);
 public:
  GfList gflst; //accumulate GffObj being read
  GfoHolder* newGffRec(GffLine* gffline, bool keepAttr, bool noExonAttr, int replaceidx=-1);
  bool addSubFeature(GfoHolder* prevgfo, GffLine* gffline, GHash<CNonExon>& pex, bool noExonAttr);
  GList<GSeqStat> gseqstats; //list of all genomic sequences seen by this reader, accumulates stats
  GffReader(FILE* f, bool justmrna=false, bool sort=false):phash(true), 
                             tids(true), gflst(sort), gseqstats(true,true,true) {
      gff_warns=false;
      names=NULL;
      gffline=NULL;
      mrnaOnly=justmrna;
      fpos=0;
      fname=NULL;
      fh=f;
      GMALLOC(linebuf, GFF_LINELEN);
      buflen=GFF_LINELEN-1;
      }
  GffReader(char* fn, bool justmrna=false, bool sort=false):phash(true),
                             tids(true),gflst(sort),gseqstats(true,true,true) {
      gff_warns=false;
      names=NULL;
      fname=Gstrdup(fn);
      mrnaOnly=justmrna;
      //sortbyloc=sort;
      fh=fopen(fname, "rb");
      fpos=0;
      gffline=NULL;
      GMALLOC(linebuf, GFF_LINELEN);
      buflen=GFF_LINELEN-1;
      }

  virtual ~GffReader() {
      gffline=NULL;
      fpos=0;
      delete gffline;
      gflst.Clear();
      phash.Clear();
      gseqstats.Clear();
      GFREE(fname);
      GFREE(linebuf);
      }

  void showWarnings(bool v=true) {
      gff_warns=v;
      }
  GffLine* nextGffLine();

  // load all subfeatures, re-group them:
  void readAll(bool keepAttr=false, bool mergeCloseExons=false, bool noExonAttr=true);

}; // end of GffReader

#endif
